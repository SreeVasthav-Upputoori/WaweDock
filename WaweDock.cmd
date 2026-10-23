@echo off
setlocal
set "WAVEDOCK_ROOT=%~dp0"
if not exist "%WAVEDOCK_ROOT%.venv\Scripts\python.exe" (
    echo WaweDock's Python environment was not found beside this launcher.
    echo Restore or install the .venv folder as described in README.md.
    pause
    exit /b 1
)
"%WAVEDOCK_ROOT%.venv\Scripts\python.exe" "%WAVEDOCK_ROOT%mds\windows_launcher.py" %*
if errorlevel 1 (
    echo WaweDock could not start. See the error dialog and startup-error.txt in:
    echo %LOCALAPPDATA%\WaweDock\WaweDock\logs
    pause
    exit /b 1
)
