"""Startup failures must be visible and logged before native GUI/chemistry imports."""
from pathlib import Path
from unittest.mock import Mock
import os
import subprocess
import sys
import pytest
from mds import app
from mds.core import startup


def test_guard_reports_dependency_import_failure_before_window(tmp_path,monkeypatch):
    error=ImportError('DLL load failed while importing rdinchi: An Application Control policy has blocked this file.')
    monkeypatch.setattr(app,'prepare_startup',lambda:tmp_path)
    def fail(*args,**kwargs): raise error
    monkeypatch.setattr(app,'_run_gui',fail)
    dialog=Mock();monkeypatch.setattr(startup,'show_startup_error',dialog)
    assert app.main(['WaweDock'])==1
    report=tmp_path/'startup-error.txt'
    assert report.exists() and 'rdinchi' in report.read_text()
    assert 'Traceback' in report.read_text() and sys.executable in report.read_text()
    dialog.assert_called_once()
    assert 'properly signed' in dialog.call_args.args[0]
    assert str(report) in dialog.call_args.args[0]


def test_check_mode_returns_failure_without_a_blocking_dialog(tmp_path,monkeypatch,capsys):
    monkeypatch.setattr(app,'prepare_startup',lambda:tmp_path)
    def fail(*args,**kwargs): raise ImportError('Missing native dependency')
    monkeypatch.setattr(app,'_run_gui',fail)
    dialog=Mock();monkeypatch.setattr(startup,'show_startup_error',dialog)
    assert app.main(['WaweDock','--startup-check'])==1
    dialog.assert_not_called()
    assert 'Missing native dependency' in capsys.readouterr().err
    assert (tmp_path/'startup-error.txt').exists()


def test_registration_does_not_import_chemistry_or_start_gui(tmp_path,monkeypatch):
    from mds import windows_integration
    monkeypatch.setattr(app,'prepare_startup',lambda:tmp_path)
    register=Mock();monkeypatch.setattr(windows_integration,'register',register)
    gui=Mock();monkeypatch.setattr(app,'_run_gui',gui)
    assert app.main(['WaweDock','--register-file-associations'])==0
    register.assert_called_once();gui.assert_not_called()


def test_success_exit_code_and_file_arguments_are_preserved(tmp_path,monkeypatch):
    monkeypatch.setattr(app,'prepare_startup',lambda:tmp_path)
    gui=Mock(return_value=0);monkeypatch.setattr(app,'_run_gui',gui)
    args=['WaweDock','C:/molecules/example.pdb']
    assert app.main(args)==0
    gui.assert_called_once_with(args,args[1:],tmp_path,check_startup=False)


def test_logging_failure_still_gets_a_diagnostic_report(tmp_path,monkeypatch):
    monkeypatch.setattr(startup.tempfile,'gettempdir',lambda:str(tmp_path))
    def fail(): raise PermissionError('Application data folder is not writable')
    monkeypatch.setattr(app,'prepare_startup',fail)
    dialog=Mock();monkeypatch.setattr(startup,'show_startup_error',dialog)
    assert app.main(['WaweDock'])==1
    assert (tmp_path/'WaweDock'/'logs'/'startup-error.txt').exists()
    dialog.assert_called_once()


def test_error_dialog_failure_does_not_hide_original_report(tmp_path,monkeypatch):
    def fail(message): raise RuntimeError('Desktop unavailable')
    monkeypatch.setattr(startup,'show_startup_error',fail)
    path=startup.report_startup_failure(ValueError('Original problem'),tmp_path)
    assert path==tmp_path/'startup-error.txt'
    assert 'Original problem' in path.read_text()


def test_import_and_spawn_launcher_do_not_load_native_dependencies(tmp_path):
    root=Path(__file__).resolve().parents[1]
    launcher=root/'mds'/'windows_launcher.py'
    code=(f'import sys,runpy;sys.path.insert(0,{str(root)!r});import mds.app;'
          f'runpy.run_path({str(launcher)!r},run_name="__mp_main__");'
          'assert "rdkit" not in sys.modules;assert "PySide6" not in sys.modules')
    result=subprocess.run([sys.executable,'-c',code],cwd=tmp_path,capture_output=True,text=True,timeout=15)
    assert result.returncode==0,result.stderr


@pytest.mark.parametrize("windowed",[False,True])
def test_fresh_launcher_records_pre_window_import_failure(tmp_path,windowed):
    root=Path(__file__).resolve().parents[1]
    launcher=root/'mds'/'windows_launcher.py'
    code=f'''
import builtins,runpy,sys
if {windowed!r}: sys.stdout=sys.stderr=None
real_import=builtins.__import__
def fail_main_window(name,*args,**kwargs):
    if name=='mds.gui.main_window':
        raise ImportError('Native dependency unavailable during startup')
    return real_import(name,*args,**kwargs)
builtins.__import__=fail_main_window
sys.argv=[{str(launcher)!r},'--startup-check']
runpy.run_path({str(launcher)!r},run_name='__main__')
'''
    env=os.environ.copy();env['LOCALAPPDATA']=str(tmp_path);env['XDG_DATA_HOME']=str(tmp_path)
    result=subprocess.run([sys.executable,'-c',code],cwd=tmp_path,env=env,capture_output=True,text=True,timeout=20)
    assert result.returncode==1,result.stderr
    reports=list(tmp_path.rglob('startup-error.txt'))
    assert reports and 'Native dependency unavailable' in reports[0].read_text()
    if windowed:
        consoles=list(tmp_path.rglob('startup-console.log'))
        assert consoles and 'Diagnostic report' in consoles[0].read_text()
    else:
        assert 'Diagnostic report' in result.stderr
    assert list(tmp_path.rglob('mds.log'))


def test_unwritable_default_log_directory_uses_temporary_fallback(tmp_path):
    root=Path(__file__).resolve().parents[1]
    unavailable=tmp_path/'not-a-folder';unavailable.write_text('occupied')
    fallback=tmp_path/'fallback'
    code=(f'import sys;sys.path.insert(0,{str(root)!r});from pathlib import Path;'
          'from mds.core import startup;'
          f'startup.default_log_dir=lambda:Path({str(unavailable)!r});'
          f'startup.tempfile.gettempdir=lambda:{str(fallback)!r};'
          f'assert startup.prepare_startup()==Path({str(fallback)!r})/"WaweDock"/"logs"')
    result=subprocess.run([sys.executable,'-c',code],cwd=tmp_path,capture_output=True,text=True,timeout=15)
    assert result.returncode==0,result.stderr
    assert (fallback/'WaweDock'/'logs'/'mds.log').exists()
