; Build the frozen app first, then: ISCC.exe installer\wavedock.iss
#define AppVersion "0.1.0"
#define AppEdition GetEnv("WAVEDOCK_BUILD_EDITION")
#if AppEdition == ""
  #define AppEdition "Full"
#endif
[Setup]
AppId=WaweDock.Desktop
AppName=WaweDock
AppVersion={#AppVersion}
DefaultDirName={localappdata}\Programs\WaweDock
DefaultGroupName=WaweDock
PrivilegesRequired=lowest
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64
MinVersion=10.0.17763
OutputDir=..\dist\installer
OutputBaseFilename=WaweDock-Setup-{#AppVersion}-{#AppEdition}
SetupIconFile=..\mds\assets\binding-pocket.ico
UninstallDisplayIcon={app}\WaweDock.exe
Compression=lzma2
SolidCompression=yes
ChangesAssociations=yes
WizardStyle=modern
#if GetEnv("WAVEDOCK_SIGNED_RELEASE") != "1"
  #error "Use scripts/build_windows_release.ps1 to build a signed, startup-checked installer."
#endif
SignTool=wavedocksign
SignedUninstaller=yes

[Files]
Source: "..\dist\WaweDock\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\WaweDock"; Filename: "{app}\WaweDock.exe"; AppUserModelID: "WaweDock.Desktop"
Name: "{autodesktop}\WaweDock"; Filename: "{app}\WaweDock.exe"; AppUserModelID: "WaweDock.Desktop"

[Run]
Filename: "{app}\WaweDock.exe"; Parameters: "--register-file-associations"; Flags: runhidden waituntilterminated
Filename: "{app}\WaweDock.exe"; Description: "Open WaweDock"; Flags: postinstall nowait skipifsilent

[UninstallRun]
Filename: "{app}\WaweDock.exe"; Parameters: "--unregister-file-associations"; RunOnceId: "RemoveWaweDockAssociations"; Flags: runhidden waituntilterminated skipifdoesntexist
