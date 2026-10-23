# Code signing policy

WaweDock's official Windows binaries are code-signed so that Windows SmartScreen
and Smart App Control can verify their origin and integrity.

## Who signs the releases

Free code-signing for this open-source project is provided by the
[SignPath Foundation](https://signpath.org/). The certificate is issued to the
project and the signing operation is performed by the SignPath.io service — the
project maintainer does **not** hold the private key and cannot sign arbitrary
files with it.

- **Certificate subject:** WaweDock (open-source project)
- **Certificate issuer / signing service:** SignPath Foundation (via SignPath.io)

## How a signed release is produced

1. A release is built from tagged source in this repository by the automated
   GitHub Actions workflow (`.github/workflows/release.yml`) — not on a
   maintainer's personal machine. The exact source commit is recorded in the
   build.
2. The workflow produces the application binaries and installer as build
   artifacts.
3. Those artifacts are submitted to SignPath.io, which signs them according to
   this project's signing policy and returns the signed files.
4. Only the signed artifacts are published as a release.

Because builds run from public source through an integrated pipeline, anyone can
reproduce the unsigned binary from the same commit and confirm it matches what
was signed.

## What is signed

- The application executable (`WaweDock.exe`)
- The bundled native libraries (`.dll`, `.pyd`) that ship inside the application
  folder, including third-party scientific libraries (for example RDKit and
  NumPy) that are redistributed under their own permissive licenses
- The Windows installer and its uninstaller

Third-party libraries are signed only to attest that **this project** built and
distributed these exact copies; it does not alter or relicense them.

## How to verify a download

On Windows, right-click the installer or `WaweDock.exe`, choose **Properties ▸
Digital Signatures**, and confirm the signature is present and valid. From
PowerShell:

```powershell
Get-AuthenticodeSignature .\WaweDock-Setup-*.exe | Format-List Status, SignerCertificate
```

`Status` should be `Valid`.

## Reporting a problem

If you find an unofficial or tampered build distributed as WaweDock, please open
an issue in this repository. Report security-sensitive matters privately to the
maintainer rather than in a public issue.

---

*Status: this project is in the process of onboarding to SignPath Foundation.
Until onboarding completes, release artifacts may be unsigned; this document
describes the signing process that applies once it is active.*
