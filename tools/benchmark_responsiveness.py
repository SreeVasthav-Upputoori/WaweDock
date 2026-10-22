from pathlib import Path
import sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from mds.chem import box as b
from mds.chem import contacts as c
root=Path(__file__).resolve().parents[1]
path=root/'tests/data/1CA2.pdb'
coords=[np.random.default_rng(i).normal(size=(100,3)) for i in range(20)]
box=b.DockingBox((0,0,0),(20,20,20))
def dense_extent(x): return float(np.linalg.norm(x[:,None,:]-x[None,:,:],axis=-1).max())
start=time.perf_counter()
for _ in range(80):
    extent=max(dense_extent(x) for x in coords)
    b.validate_box(box,path,extent)
baseline=time.perf_counter()-start
cache=b.BoxGeometryCache()
start=time.perf_counter()
for _ in range(80):
    extent=max(cache.extent(x) for x in coords)
    b.validate_box(box,ligand_extent=extent,receptor_xyz=cache.file_coords(path))
updated=time.perf_counter()-start
print(f'80 box validations / 20 ligand geometries: baseline={baseline:.3f}s cached={updated:.3f}s speedup={baseline/updated:.1f}x')
# Use translated copies so atom density stays realistic while atom count grows.
base=c.load_protein_atoms(root/'tests/data/3PTB.pdb')
xyz=np.concatenate([base.xyz+[100*i,0,0] for i in range(8)])
el=base.element*8
hetero=np.array([e in ('N','O','S') for e in el])
start=time.perf_counter()
flags=[]
for i,e in enumerate(el):
    if e=='C': flags.append(bool((np.linalg.norm(xyz[hetero]-xyz[i],axis=1)<1.7).any()))
print(f'Dense polar-carbon scan / {len(xyz)} protein atoms: {time.perf_counter()-start:.3f}s')

start=time.perf_counter()
fast=c._polar_carbon_flags(xyz,el)
print(f'Spatial polar-carbon scan / {len(xyz)} protein atoms: {time.perf_counter()-start:.3f}s')
np.testing.assert_array_equal(fast[np.array(el)=='C'],flags)
