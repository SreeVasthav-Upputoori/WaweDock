"""Installed feature declarations and fail-closed release prerequisites."""
import builtins
import json
from pathlib import Path
import subprocess
import sys
import pytest
from mds.core import paths,capabilities
from mds import windows_integration as integration
from test_windows_integration import registry


def test_build_edition_is_explicit_and_invalid_metadata_is_rejected(tmp_path,monkeypatch):
    monkeypatch.setattr(paths,'package_root',lambda:tmp_path)
    assert paths.build_edition()=='Full'
    (tmp_path/'build-info.json').write_text(json.dumps({'edition':'Protein'}))
    assert paths.build_edition()=='Protein'
    (tmp_path/'build-info.json').write_text(json.dumps({'edition':'invalid'}))
    with pytest.raises(ValueError,match='edition'): paths.build_edition()


def test_protein_edition_never_attempts_to_import_rdkit(monkeypatch):
    monkeypatch.setattr(paths,'build_edition',lambda:'Protein')
    actual_import=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name.startswith('rdkit'): raise AssertionError('Protein edition must not attempt RDKit import')
        return actual_import(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    status=capabilities.probe_in_process()
    assert not status.available and 'protein tools only' in status.reason


def test_protein_registration_removes_only_owned_unsupported_formats(registry,monkeypatch):
    monkeypatch.setattr(integration,'build_edition',lambda:'Full')
    integration.register()
    sdf=integration.CLASSES+r'\.sdf\OpenWithProgids'
    registry.values[sdf]['OtherApp']=''
    monkeypatch.setattr(integration,'build_edition',lambda:'Protein')
    assert integration.register()
    assert integration.PROGID not in registry.values[sdf]
    assert registry.values[sdf]['OtherApp']==''
    assert integration.PROGID in registry.values[integration.CLASSES+r'\.pdb\OpenWithProgids']
    assert not integration.register()


@pytest.mark.skipif(sys.platform!='win32',reason='Windows release preflight')
def test_release_refuses_to_build_without_a_signing_certificate():
    root=Path(__file__).resolve().parents[1]
    result=subprocess.run(['powershell','-NoProfile','-NonInteractive','-File',str(root/'scripts/build_windows_release.ps1'),
        '-CertificateThumbprint','0'*40],cwd=root,capture_output=True,text=True,timeout=20)
    assert result.returncode!=0
    if 'running scripts is disabled' in result.stderr:
        pytest.skip('Machine PowerShell policy blocks release-script execution; no bypass applied. Syntax is checked separately.')
    assert 'No unsigned public installer will be produced' in result.stderr
