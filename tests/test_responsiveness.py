"""Numerical equivalence, cache invalidation, and nonblocking latest-request jobs."""
from pathlib import Path
import threading
import time
import numpy as np
import pytest
from PySide6.QtCore import QTimer
from conftest import wait_until
from mds.chem import box as boxlib
from mds.chem import contacts as contactlib
from mds.jobs import Job,JobRunner
from mds.jobs.latest import LatestJobQueue

DATA=Path(__file__).parent/'data'


def test_box_cache_reads_once_and_invalidates_after_file_rewrite(tmp_path,monkeypatch):
    path=tmp_path/'protein.pdbqt'
    def write(x):
        path.write_text(f'ATOM      1  CA  ALA A   1    {x:8.3f}{0:8.3f}{0:8.3f}\n')
    write(0)
    original=boxlib.receptor_atom_coords
    reads=[]
    def read(p): reads.append(p); return original(p)
    monkeypatch.setattr(boxlib,'receptor_atom_coords',read)
    cache=boxlib.BoxGeometryCache()
    a=cache.file_coords(path)
    for center in range(20):
        report=boxlib.validate_box(boxlib.DockingBox((center,0,0),(10,10,10)),receptor_xyz=cache.file_coords(path))
        assert bool(report.errors)==(center>5)
    assert len(reads)==1
    with pytest.raises(ValueError): a[0,0]=9
    write(25)
    assert cache.file_coords(path)[0,0]==25
    assert len(reads)==2
    cache.clear()
    cache.file_coords(path)
    assert len(reads)==3


def test_cached_extent_tracks_coordinate_mutation_and_is_bounded(monkeypatch):
    cache=boxlib.BoxGeometryCache(max_extents=2)
    original=boxlib.ligand_extent
    calls=[]
    def compute(coords): calls.append(1); return original(coords)
    monkeypatch.setattr(boxlib,'ligand_extent',compute)
    coords=np.array([[0.,0,0],[3,4,0]])
    assert cache.extent(coords)==5
    assert cache.extent(coords.copy())==5 and len(calls)==1
    coords[1]=[0,0,10]
    assert cache.extent(coords)==10 and len(calls)==2
    cache.extent(np.array([[0.,0,0],[0,0,12]]))
    assert len(cache.extents)==2


def test_chunked_extent_and_contacts_equal_dense_geometry():
    random=np.random.default_rng(17)
    ligand=random.normal(size=(530,3))
    dense=np.linalg.norm(ligand[:,None,:]-ligand[None,:,:],axis=-1)
    assert boxlib.ligand_extent(ligand)==pytest.approx(dense.max(),rel=1e-14)
    protein=random.normal(size=(1200,3))
    blocks=list(contactlib.contact_distance_blocks(protein,ligand,block_size=73))
    assert max(len(distances) for _,distances in blocks)<=73
    assert [start for start,_ in blocks]==list(range(0,len(protein),73))
    np.testing.assert_array_equal(np.concatenate([d for _,d in blocks]),
                                  np.linalg.norm(protein[:,None,:]-ligand[None,:,:],axis=-1))


def test_protein_spatial_buckets_match_dense_polar_carbon_flags():
    protein=contactlib.load_protein_atoms(DATA/'3PTB.pdb')
    hetero=np.array([e in ('N','O','S') for e in protein.element])
    expected=np.zeros(len(protein.xyz),dtype=bool)
    for i,element in enumerate(protein.element):
        if element=='C':
            expected[i]=np.any(np.linalg.norm(protein.xyz[hetero]-protein.xyz[i],axis=1)<1.7)
    np.testing.assert_array_equal(protein.polar_carbon,expected)


def test_protein_cache_invalidates_changed_files(tmp_path,monkeypatch):
    path=tmp_path/'protein.pdb'
    path.write_bytes((DATA/'3PTB.pdb').read_bytes())
    contactlib._cached_protein_atoms.cache_clear()
    original=contactlib.load_protein_atoms
    reads=[]
    def read(p): reads.append(p); return original(p)
    monkeypatch.setattr(contactlib,'load_protein_atoms',read)
    first=contactlib.cached_protein_atoms(path)
    assert contactlib.cached_protein_atoms(path) is first
    path.write_bytes(path.read_bytes()+b'REMARK changed\n')
    assert contactlib.cached_protein_atoms(path) is not first
    assert len(reads)==2


class ControlledJob(Job):
    def __init__(self,name,started,gate,runs,fail=False):
        self.name,self.started,self.gate,self.runs,self.fail=name,started,gate,runs,fail
    def run(self,ctx):
        self.runs.append((self.name,threading.get_ident()))
        self.started.set()
        while not self.gate.wait(0.01): ctx.check_cancelled()
        ctx.check_cancelled()
        if self.fail: raise RuntimeError('deliberate failure')
        return self.name


def test_rapid_requests_keep_latest_and_ui_heartbeat_running(qapp):
    runner=JobRunner(); queue=LatestJobQueue(runner)
    started=threading.Event(); gate=threading.Event(); immediate=threading.Event(); immediate.set()
    runs=[]; results=[]; errors=[]; ticks=[]
    timer=QTimer(); timer.setInterval(1); timer.timeout.connect(lambda:ticks.append(1)); timer.start()
    try:
        queue.submit(ControlledJob('first',started,gate,runs),results.append,errors.append)
        assert wait_until(started.is_set)
        queue.submit(ControlledJob('discard',threading.Event(),immediate,runs),results.append,errors.append)
        queue.submit(ControlledJob('latest',threading.Event(),immediate,runs),results.append,errors.append)
        assert wait_until(lambda:len(ticks)>=3)
        assert not results
        gate.set()
        assert wait_until(lambda:results==['latest'])
        assert [name for name,_ in runs]==['first','latest']
        assert all(thread!=threading.get_ident() for _,thread in runs)
        assert not errors
    finally:
        timer.stop(); gate.set(); runner.wait_for_done(10000)


def test_queue_failure_and_clear_do_not_publish_stale_results(qapp):
    runner=JobRunner(); queue=LatestJobQueue(runner)
    gate=threading.Event(); started=threading.Event(); runs=[]; results=[]; errors=[]
    try:
        queue.submit(ControlledJob('failure',started,gate,runs,True),results.append,errors.append)
        assert wait_until(started.is_set)
        gate.set()
        assert wait_until(lambda:bool(errors))
        assert 'deliberate failure' in errors[0]
        gate.clear(); started.clear()
        queue.submit(ControlledJob('cancelled',started,gate,runs),results.append,errors.append)
        assert wait_until(started.is_set)
        queue.clear(); gate.set()
        assert wait_until(lambda:queue.active is None)
        assert not results and len(errors)==1
    finally:
        gate.set(); runner.wait_for_done(10000)


def test_results_panel_integration_discards_stale_contact_results(qapp,tmp_path,monkeypatch):
    # Requires RDKit; the independent queue and numerical tests run without it.
    from types import SimpleNamespace
    from mds.gui.controller import Controller
    from mds.gui.panels.results_panel import ResultsPanel
    from mds.project import Project
    from test_results_workspace import sample_outcomes
    controller=Controller(JobRunner())
    controller.set_project(Project.create(tmp_path/'project'))
    controller.receptor=SimpleNamespace(prepared_pdb=tmp_path/'receptor.pdb')
    requests=[]
    monkeypatch.setattr(controller,'request_contacts',lambda o,i:requests.append((o,i)))
    def forbidden(*args): raise AssertionError('Synchronous contacts must not run on pose selection')
    monkeypatch.setattr(controller,'contacts_for',forbidden)
    panel=ResultsPanel(controller)
    try:
        controller.outcomes=sample_outcomes()
        panel.show_outcomes(controller.outcomes)
        current=panel._current_pose
        assert requests and requests[-1]==current
        wrong=next(o for o in controller.outcomes if o is not current[0])
        contacts=[contactlib.Contact('ALA1:A',3.2,['hydrophobic contact candidate'])]
        controller.contacts_ready.emit(wrong,0,contacts)
        assert panel.contacts.rowCount()==0
        controller.contacts_ready.emit(current[0],current[1],contacts)
        assert panel.contacts.rowCount()==1
        assert panel.contacts.item(0,0).text()=='ALA1:A'
        panel.search.setText('ligand A')
        assert panel._current_pose[0] is not current[0]
        controller.contacts_ready.emit(current[0],current[1],contacts)
        assert panel.contacts.rowCount()==0
        panel.search.setText('no matching ligand')
        assert panel._current_pose is None
    finally:
        panel.close()
