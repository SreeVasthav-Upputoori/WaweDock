"""Optional chemistry availability, probed outside the GUI process."""
from __future__ import annotations
from dataclasses import asdict,dataclass
import threading
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time


@dataclass(frozen=True)
class ChemistryStatus:
    available: bool
    reason: str = ''
    version: str = ''


def probe_in_process():
    """Import real backends; blocked modules are never patched, skipped or emulated."""
    from mds.core.paths import build_edition
    if build_edition()=='Protein':
        return ChemistryStatus(False,'This installation contains protein tools only. The Full edition is required for ligand preparation and docking.')
    try:
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import (Crippen,Descriptors,Lipinski,rdChemReactions,rdDistGeom,
            rdForceFieldHelpers,rdMolDescriptors,rdFingerprintGenerator)
        from rdkit.ML.Cluster import Butina
        from meeko import MoleculePreparation,PDBQTWriterLegacy
        return ChemistryStatus(True,version=rdkit.__version__)
    except Exception as exc:
        return ChemistryStatus(False,f'{type(exc).__name__}: {exc}')


def write_probe(path):
    status=probe_in_process()
    Path(path).write_text(json.dumps(asdict(status)),encoding='utf-8')
    return 0 if status.available else 1


_status=None
_lock=threading.Lock()


def check_chemistry(check_cancelled=lambda:None):
    global _status
    while not _lock.acquire(timeout=.05): check_cancelled()
    try:
        check_cancelled()
        if _status is None:
            try: _status=_probe(check_cancelled)
            except (OSError,ValueError,subprocess.SubprocessError) as exc:
                _status=ChemistryStatus(False,f'Chemistry dependency check failed: {exc}')
        return _status
    finally: _lock.release()


def _probe(check_cancelled):
    """The worker uses an output file so pythonw/frozen builds need no console."""
    with tempfile.TemporaryDirectory(prefix='wavedock-chemistry-') as folder:
        output=Path(folder)/'result.json'
        if getattr(sys,'frozen',False): command=[sys.executable]
        else:
            python=Path(sys.executable)
            if python.name.lower()=='pythonw.exe': python=python.with_name('python.exe')
            command=[str(python),str(Path(__file__).resolve().parents[1]/'windows_launcher.py')]
        command+=['--probe-chemistry',str(output)]
        try:
            process=subprocess.Popen(command,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,
                creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
            try:
                deadline=time.monotonic()+20
                while process.poll() is None:
                    check_cancelled()
                    if time.monotonic()>=deadline:
                        process.terminate();process.wait(timeout=5)
                        return ChemistryStatus(False,'Chemistry dependency check timed out. Restart after checking the installation.')
                    time.sleep(.05)
                if output.exists():
                    data=json.loads(output.read_text(encoding='utf-8'))
                    if isinstance(data.get('available'),bool):
                        return ChemistryStatus(data['available'],str(data.get('reason',''))[:2000],str(data.get('version','')))
                return ChemistryStatus(False,f'Chemistry dependency process exited without a result (code {process.returncode}).')
            finally:
                if process.poll() is None: process.kill();process.wait(timeout=5)
        except (OSError,ValueError,subprocess.SubprocessError) as exc:
            return ChemistryStatus(False,f'Chemistry dependency check failed: {exc}')


class ProbeChemistryJob:
    name='Check optional ligand and docking tools'
    def run(self,ctx):
        ctx.check_cancelled()
        result=check_chemistry(ctx.check_cancelled)
        ctx.check_cancelled()
        return result
