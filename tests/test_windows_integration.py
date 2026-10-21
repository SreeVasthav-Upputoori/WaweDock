"""Exercise install, repair and uninstall without changing the user's registry."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock
import pytest
from mds import windows_integration as integration


class Key(str):
    def __enter__(self): return self
    def __exit__(self, *args): pass


class Registry:
    HKEY_CURRENT_USER = 1
    REG_SZ = 1
    KEY_SET_VALUE = 2

    def __init__(self):
        self.values = {}
        self.writes = 0

    def OpenKey(self, root, key, *args):
        if key not in self.values: raise FileNotFoundError(key)
        return Key(key)

    def CreateKey(self, root, key):
        self.values.setdefault(key, {})
        return Key(key)

    def QueryValueEx(self, key, name):
        if name not in self.values[key]: raise FileNotFoundError(name)
        return self.values[key][name], self.REG_SZ

    def SetValueEx(self, key, name, reserved, kind, value):
        self.values[key][name] = value
        self.writes += 1

    def DeleteValue(self, key, name):
        del self.values[key][name]

    def DeleteKey(self, root, key):
        if key not in self.values: raise FileNotFoundError(key)
        if self.values[key] or any(k.startswith(key + "\\") for k in self.values):
            raise OSError("not empty")
        del self.values[key]


@pytest.fixture
def registry(monkeypatch):
    import sys
    registry = Registry()
    monkeypatch.setitem(sys.modules, "winreg", registry)
    monkeypatch.setattr(integration, "sys", SimpleNamespace(platform="win32", executable=r"C:\Program Files\WaweDock\WaweDock.exe", frozen=True))
    monkeypatch.setattr(integration, "assets_dir", lambda: Path("C:/Program Files/WaweDock/_internal/mds/assets"))
    monkeypatch.setattr(integration, "ctypes", SimpleNamespace(windll=SimpleNamespace(shell32=SimpleNamespace(SHChangeNotify=Mock()))))
    return registry


def test_frozen_registration_is_portable_idempotent_and_preserves_defaults(registry):
    choice = r"Software\Microsoft\Windows\CurrentVersion\Explorer\FileExts\.pdb\UserChoice"
    registry.values[choice] = {"ProgId": "OtherApp", "Hash": "existing hash"}
    registry.values[integration.CLASSES + r"\.pdb"] = {"": "OtherApp"}
    assert integration.register()
    writes = registry.writes
    assert not integration.register()
    assert registry.writes == writes
    command = integration.launch_command()
    assert command == '\"C:\\Program Files\\WaweDock\\WaweDock.exe\" \"%1\"'
    assert registry.values[choice]["ProgId"] == "OtherApp"
    assert registry.values[integration.CLASSES + r"\.pdb"][""] == "OtherApp"
    assert not any("tools" in value or "pythonw" in value for _,_,value in integration.registration_entries())
    assert integration.unregister()
    assert registry.values[choice] == {"ProgId": "OtherApp", "Hash": "existing hash"}
    assert registry.values[integration.CLASSES + r"\.pdb"][""] == "OtherApp"
    assert not integration.unregister()


def test_old_install_cannot_unregister_new_install(registry, monkeypatch):
    integration.register()
    monkeypatch.setattr(integration.sys, "executable", r"D:\New location\WaweDock.exe")
    assert integration.register()
    new_command = integration.launch_command()
    monkeypatch.setattr(integration.sys, "executable", r"C:\Program Files\WaweDock\WaweDock.exe")
    assert not integration.unregister()
    key = integration.CLASSES + "\\" + integration.PROGID + r"\shell\open\command"
    assert registry.values[key][""] == new_command


def test_source_registration_uses_packaged_launcher(registry, monkeypatch):
    monkeypatch.setattr(integration.sys, "frozen", False)
    monkeypatch.setattr(integration.sys, "executable", "C:/Python Env/python.exe")
    command = integration.launch_command()
    assert "pythonw.exe" in command and "windows_launcher.py" in command
    assert "tools" not in command
    assert command.endswith(' "%1"')
