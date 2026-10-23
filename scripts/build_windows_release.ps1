param(
    [Parameter(Mandatory=$true)][ValidatePattern('^[0-9A-Fa-f]{40}$')][string]$CertificateThumbprint,
    [ValidateSet('Full','Protein')][string]$Edition='Full',
    [string]$SignTool='signtool.exe',
    [string]$InnoSetup='ISCC.exe'
)
$ErrorActionPreference='Stop'
if ([Environment]::OSVersion.Platform -ne 'Win32NT') { throw 'Build Windows releases on Windows.' }
$releaseRoot=Split-Path -Parent $PSScriptRoot
$releasePython=Join-Path $releaseRoot '.venv\Scripts\python.exe'
$releaseCert=Get-Item -LiteralPath ('Cert:\CurrentUser\My\'+$CertificateThumbprint) -ErrorAction SilentlyContinue
if (-not $releaseCert -or -not $releaseCert.HasPrivateKey -or $releaseCert.NotAfter -lt (Get-Date)) {
    throw 'A valid code-signing certificate with a private key is required. No unsigned public installer will be produced.'
}
$signCommand=(Get-Command $SignTool -ErrorAction Stop).Source
$innoCommand=(Get-Command $InnoSetup -ErrorAction Stop).Source
if (-not (Test-Path -LiteralPath $releasePython)) { throw 'Install the build environment into .venv first.' }
$oldEdition=$env:WAVEDOCK_BUILD_EDITION
$oldSigned=$env:WAVEDOCK_SIGNED_RELEASE
Push-Location $releaseRoot
try {
    $env:WAVEDOCK_BUILD_EDITION=$Edition
    $env:WAVEDOCK_SIGNED_RELEASE='1'
    if ($Edition -eq 'Full') {
        & $releasePython -c "from mds.core.capabilities import probe_in_process; import sys; status=probe_in_process(); print(status.reason or status.version); sys.exit(0 if status.available else 1)"
        if ($LASTEXITCODE -ne 0) { throw 'Full release requires working chemistry on the approved build host; skipped chemistry tests cannot qualify it.' }
        & $releasePython -m pytest tests/test_science.py tests/test_workflow_gui.py -q
        if ($LASTEXITCODE -ne 0) { throw 'Full scientific/preparation/docking validation failed.' }
    }
    & $releasePython -m PyInstaller --noconfirm --clean wavedock.spec
    if ($LASTEXITCODE -ne 0) { throw 'PyInstaller failed.' }
    $packageDirectory=Join-Path $releaseRoot 'dist\WaweDock'
    $nativeFiles=Get-ChildItem -LiteralPath $packageDirectory -Recurse -File | Where-Object { $_.Extension -in '.exe','.dll','.pyd' }
    if (-not $nativeFiles) { throw 'No native application files were built.' }
    foreach ($native in $nativeFiles) {
        $signature=Get-AuthenticodeSignature -LiteralPath $native.FullName
        if ($signature.Status -ne 'Valid') {
            & $signCommand sign /sha1 $CertificateThumbprint /fd SHA256 /tr http://timestamp.digicert.com /td SHA256 $native.FullName
            if ($LASTEXITCODE -ne 0) { throw ('Code signing failed: '+$native.FullName) }
        }
        if ((Get-AuthenticodeSignature -LiteralPath $native.FullName).Status -ne 'Valid') {
            throw ('Native file has no trusted valid signature: '+$native.FullName)
        }
    }
    # Test the signed executable from outside the checkout. Full releases must
    # load the actual chemistry backend; Protein releases explicitly omit it.
    $startupFlag=if ($Edition -eq 'Full') { '--startup-check-full' } else { '--startup-check' }
    $releaseExecutable=Join-Path $packageDirectory 'WaweDock.exe'
    $probe=Start-Process -FilePath $releaseExecutable -ArgumentList $startupFlag -WorkingDirectory ([IO.Path]::GetTempPath()) -WindowStyle Hidden -PassThru
    if (-not $probe.WaitForExit(60000)) { Stop-Process -Id $probe.Id; throw 'Frozen startup check timed out.' }
    $probe.Refresh()
    if ($probe.ExitCode -ne 0) { throw ('Frozen '+$Edition+' startup check failed. Installer was not built.') }
    # Inno signs both its generated uninstaller and the final setup executable.
    $innoSign='/Swavedocksign="'+$signCommand+'" sign /sha1 '+$CertificateThumbprint+' /fd SHA256 /tr http://timestamp.digicert.com /td SHA256 $f'
    & $innoCommand $innoSign (Join-Path $releaseRoot 'installer\wavedock.iss')
    if ($LASTEXITCODE -ne 0) { throw 'Installer compilation or signing failed.' }
    $setups=Get-ChildItem -LiteralPath (Join-Path $releaseRoot 'dist\installer') -Filter ('WaweDock-Setup-*-'+$Edition+'.exe') -File
    if (-not $setups) { throw 'Signed installer was not produced.' }
    foreach ($setup in $setups) {
        if ((Get-AuthenticodeSignature -LiteralPath $setup.FullName).Status -ne 'Valid') { throw ('Installer signature is invalid: '+$setup.FullName) }
        Write-Output ('Signed and startup-checked: '+$setup.FullName)
    }
}
finally {
    $env:WAVEDOCK_BUILD_EDITION=$oldEdition
    $env:WAVEDOCK_SIGNED_RELEASE=$oldSigned
    Pop-Location
}
