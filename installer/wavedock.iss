; Build the frozen app first, then: ISCC.exe installer\wavedock.iss
#define AppVersion "0.1.0"
[Setup]
AppId=WaweDock.Desktop
AppName=WaweDock
AppVersion={#AppVersion}
DefaultDirName={localappdata}\Programs\WaweDock
DefaultGroupName=WaweDock
PrivilegesRequired=lowest
OutputDir=..\dist\installer
OutputBaseFilename=WaweDock-Setup-{#AppVersion}
SetupIconFile=..\mds\assets\binding-pocket.ico
UninstallDisplayIcon={app}\WaweDock.exe
Compression=lzma2
SolidCompression=yes
ChangesAssociations=yes
WizardStyle=modern

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
