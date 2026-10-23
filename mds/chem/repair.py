"""Conservative, additive protein repairs with immutable proposals and provenance.

No loop construction, nonstandard-residue replacement, or movement of existing
coordinates is permitted. PDBFixer is restricted to explicitly listed atoms.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from hashlib import sha256
from importlib.metadata import version, PackageNotFoundError
import io
import json
import math
import random
import threading
from pathlib import Path
from uuid import uuid4

import gemmi

from mds.chem.edit import atom_records
from mds.chem.health import HealthReport, SIDECHAINS, VDW, check_protein
from mds.chem.receptor import read_structure

_hydrogen_lock = threading.Lock()


def fingerprint(path):
    return sha256(Path(path).read_bytes()).hexdigest()


@dataclass
class RepairChange:
    id: str
    kind: str
    chain_index: int
    residue_index: int
    label: str
    names: list[str] = field(default_factory=list)
    conformer: str = ''
    recommended: bool = True


@dataclass
class RepairPlan:
    source: str
    source_sha256: str
    changes: list[RepairChange]
    notes: list[str]
    ph: float
    before: HealthReport
    backend: dict
    segment_starts: list = field(default_factory=list)


@dataclass
class RepairResult:
    plan: RepairPlan
    selected: list[str]
    candidate: str
    after: HealthReport
    added: list[dict]
    removed: list[dict]
    blockers: list[str]
    warnings: list[str]
    candidate_sha256: str
    inherited: list[dict] = field(default_factory=list)
    normalised: list[dict] = field(default_factory=list)

    @property
    def can_accept(self):
        return not self.blockers and bool(self.added or self.removed or self.normalised)

    def export_report(self, path):
        data = asdict(self)
        data['status'] = 'Review required' if self.can_accept else 'Blocked'
        data['parameters'] = {
            'heavy_atom_seed':42,'hydrogen_seed':42,'cpu_threads':2,'existing_coordinates_fixed':True,
            'hydrogen_forcefield':'amber14/protein.ff14SB.xml + amber14/tip3p.xml',
            'general_refinement':False,'loop_reconstruction':False,'implicit_terminal_atoms':False,
            'broad_bond_length_A':{'heavy':[1.0,1.85],'sulfur_heavy':[1.5,2.25],'hydrogen':[0.8,1.55]},
            'broad_bond_angle_degrees':[60,175],'severe_vdw_overlap_A':0.8,
            'donor_H_acceptor_min_angle_degrees':105,'donor_H_acceptor_min_distance_A':1.4,
        }
        data['change_status'] = {c.id: ('Built; review required' if self.can_accept else 'Candidate blocked') if c.id in self.selected else 'Not selected' for c in self.plan.changes}
        data['coordinate_provenance'] = 'Original surviving coordinates unchanged. Added atoms are computationally reconstructed; input coordinate provenance is unknown.'
        Path(path).write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')


def propose_repairs(source, *, selected_residues=None, ph=7.4,
                    progress=lambda *_: None, check_cancelled=lambda: None):
    """Produce a selectable preview. Manual scope uses exact model/chain/residue IDs."""
    if not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError('pH must be finite and between 0 and 14.')
    source = Path(source)
    before = check_protein(source,progress=progress,check_cancelled=check_cancelled)
    st = read_structure(source)
    backend = {}
    for name in ('pdbfixer','openmm','biopython'):
        try: backend[name] = version(name)
        except PackageNotFoundError: backend[name] = 'Not installed'
    available = all(v != 'Not installed' for v in backend.values())
    notes = [
        'Only model 1 is repaired. Other models and all nonstandard residues, cofactors, ligands, metals and waters are preserved.',
        'No loops, missing residues, backbone N/CA/C rebuilding or general geometry refinement in this increment.',
        'Existing coordinates remain fixed. Every reconstructed atom requires manual review; validation does not establish docking readiness.',
    ]
    if not available: notes.append('Install wavedock[repair] for template atom rebuilding and hydrogen addition. Conformer selection remains available.')
    segment_starts=[]
    if st.input_format==gemmi.CoorFormat.Pdb:
        base=atom_records(st)
        split=atom_records(gemmi.read_pdb(str(source),split_chain_on_ter=True))
        if len(base)!=len(split) or any((a['atom'],a['chain'],a['resi'])!=(b['atom'],b['chain'],b['resi']) for a,b in zip(base,split)):
            raise ValueError('Cannot map TER segments safely for repair.')
        previous=None
        for a,b in zip(base,split):
            segment=(b['model'],b['chain_index'])
            if segment!=previous and a['model']==0:
                segment_starts.append((a['chain_index'],a['residue_index']))
            previous=segment
    changes = []
    hydrogens = 0
    protein_keys = set()
    h_possible = available
    # Explicit connectivity to altlocs is not rewritten automatically.
    conect = {serial for pair in st.conect_map.items() for serial in (pair[0], *pair[1])}
    for ci,ch in enumerate(st[0]):
        for ri,res in enumerate(ch):
            check_cancelled()
            info = gemmi.find_tabulated_residue(res.name)
            if not (info and info.is_amino_acid()): continue
            protein_keys.add((0,ci,ri))
            hydrogens += sum(a.is_hydrogen() for a in res)
            if res.name not in SIDECHAINS:
                h_possible = False
                continue
            alts = sorted({a.altloc for a in res if a.altloc not in ('\0',' ')})
            label = f'{ch.name or "(blank)"}:{res.seqid} {res.name}'
            in_scope = selected_residues is None or (0,ci,ri) in selected_residues
            linked = any(a.serial in conect for a in res) or any(
                p.chain_name==ch.name and p.res_id.seqid==res.seqid
                for conn in st.connections for p in (conn.partner1,conn.partner2))
            if alts:
                if linked:
                    h_possible = False
                    notes.append(f'{label}: explicit connectivity present; automatic conformer resolution is unavailable.')
                elif in_scope:
                    occupancy = {alt:sum(a.occ for a in res if a.altloc==alt)/sum(1 for a in res if a.altloc==alt) for alt in alts}
                    best = max(alts,key=lambda alt:(occupancy[alt],alt=='A'))
                    for alt in alts:
                        changes.append(RepairChange(f'alt:{ci}:{ri}:{alt}','altloc',ci,ri,
                            f'{label}: keep conformer {alt} (mean occupancy {occupancy[alt]:.2f}); discard other alternate atoms',
                            conformer=alt,recommended=alt==best))
            shared = {a.name for a in res if a.altloc in ('\0',' ')}
            expected = set(('N CA C O '+SIDECHAINS[res.name]).split())
            # Missing-atom names are resolved again after the chosen conformer is applied.
            missing = sorted(set.union(*(expected-(shared|{a.name for a in res if a.altloc==alt}) for alt in (alts or ['']))))
            anchors = all(name in shared or (bool(alts) and all(any(a.name==name and a.altloc==alt for a in res) for alt in alts))
                          for name in ('N','CA','C'))
            if not anchors:
                h_possible = False
                notes.append(f'{label}: incomplete N/CA/C anchors; requires advanced modelling.')
            elif missing and available and in_scope:
                changes.append(RepairChange(f'heavy:{ci}:{ri}','heavy_atoms',ci,ri,
                    f'{label}: rebuild missing template atoms ({", ".join(missing)})',missing))
            elif missing:
                h_possible = False
    if any(i.code in ('backbone_connectivity','invalid_coordinates') or
           (i.code=='unmodelled_sequence' and i.severity=='Warning') for i in before.issues):
        h_possible = False
        notes.append('Hydrogen addition is unavailable until backbone gaps/invalid coordinates are resolved.')
    if selected_residues is not None:
        h_possible = h_possible and protein_keys <= set(selected_residues)
        notes.append('Manual scope: only the selected canonical residues are proposed; whole-protein hydrogens require all protein residues to be selected.')
    for ci,ch in enumerate(st[0]):
        canonical=[(ri,r) for ri,r in enumerate(ch) if r.name in SIDECHAINS]
        for i,(ri,res) in enumerate(canonical):
            terminal=(i==len(canonical)-1 or (ci,canonical[i+1][0]) in segment_starts or
                      (res.subchain and canonical[i+1][1].subchain and res.subchain!=canonical[i+1][1].subchain))
            if terminal and res.find_atom('OXT','*') is None:
                h_possible=False
                notes.append(f'{ch.name}:{res.seqid}: terminal OXT absent; automatic force-field hydrogen addition is unavailable. No terminal oxygen is added implicitly.')
    if hydrogens:
        h_possible = False
        notes.append('Protein already contains hydrogens; automatic reprotonation/removal of existing hydrogens is not supported.')
    if h_possible and protein_keys:
        changes.append(RepairChange('hydrogens','hydrogens',-1,-1,
            f'Add protein hydrogens at pH {ph:g} using default residue pKas; preserve heavy atoms and non-protein components',recommended=False))
    if fingerprint(source) != before.source_sha256:
        raise ValueError('Source changed during proposal; generate a new preview.')
    return RepairPlan(str(source),before.source_sha256,changes,notes,ph,before,backend,segment_starts)


def _mmcif_with_conect(st):
    """Represent explicit PDB bonds as standard mmCIF struct_conn records."""
    from collections import Counter
    work=st.clone()
    atoms={a.serial:(c,r,a) for c in work[0] for r in c for a in r}
    existing={frozenset((str(c.partner1),str(c.partner2))):c.name for c in work.connections}
    pairs={}
    for serial,partners in work.conect_map.items():
        for other,count in Counter(partners).items():
            key=tuple(sorted((serial,other)))
            pairs[key]=max(pairs.get(key,0),count)
    orders={}
    for i,((a,b),order) in enumerate(sorted(pairs.items())):
        if a not in atoms or b not in atoms or a==b: continue
        aa=gemmi.make_address(*atoms[a]);bb=gemmi.make_address(*atoms[b])
        pair=frozenset((str(aa),str(bb)))
        if pair in existing: name=existing[pair]
        else:
            connection=gemmi.Connection()
            connection.name=f'wavedock_bond_{i}'
            connection.type=gemmi.ConnectionType.Covale
            connection.partner1=aa;connection.partner2=bb
            connection.reported_distance=atoms[a][2].pos.dist(atoms[b][2].pos)
            work.connections.append(connection)
            name=connection.name
        orders[name]=('sing','doub','trip')[min(order,3)-1]
    doc=work.make_mmcif_document()
    block=doc.sole_block()
    category=block.get_mmcif_category('_struct_conn.')
    if orders and category:
        category['pdbx_value_order']=[orders.get(name,'?') for name in category['id']]
        block.set_mmcif_category('_struct_conn.',category)
    return doc


def _pdb_bonds_from_mmcif(st,source):
    """Retain component/explicit bonds in PDB, whose atom serials must be numeric."""
    from itertools import product
    for i,a in enumerate((a for m in st for c in m for r in c for a in r),1): a.serial=i
    block=gemmi.cif.read_file(str(source)).sole_block()
    pair_orders={}
    def bond(a,b,order):
        if a is not None and b is not None and a.serial!=b.serial:
            key=tuple(sorted((a.serial,b.serial)))
            pair_orders[key]=max(pair_orders.get(key,0),order)
    def order(value): return {'DOUB':2,'TRIP':3}.get(str(value).upper(),1)
    table=block.find(['_chem_comp_bond.comp_id','_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2','_chem_comp_bond.value_order'])
    for row in table:
        comp,left,right,value=[gemmi.cif.as_string(v) for v in row]
        for ch in st[0]:
            for res in ch:
                if res.name!=comp or comp in SIDECHAINS: continue
                for a,b in product([a for a in res if a.name==left],[a for a in res if a.name==right]):
                    if a.altloc in ('\0',' ') or b.altloc in ('\0',' ') or a.altloc==b.altloc: bond(a,b,order(value))
    category=block.get_mmcif_category('_struct_conn.')
    explicit_orders=dict(zip(category.get('id',[]),category.get('pdbx_value_order',[])))
    for connection in st.connections:
        endpoints=[]
        for partner in (connection.partner1,connection.partner2):
            cra=st[0].find_cra(partner)
            endpoints.append(cra.atom)
        bond(*endpoints,order(explicit_orders.get(connection.name,'sing')))
    st.clear_conect()
    for (a,b),value in pair_orders.items(): st.add_conect(a,b,value)


def _write_structure(st, path, added, source_hash, source_path=None):
    """Retain original PDB serials/CONECT, or add provenance to the mmCIF document."""
    if path.suffix.lower()=='.pdb':
        if max((a.serial for m in st for c in m for r in c for a in r),default=0)>99999:
            raise ValueError('PDB atom serial limit exceeded; import mmCIF to repair this structure.')
        options = gemmi.PdbWriteOptions()
        options.preserve_serial = True
        options.conect_records = True
        text = st.make_pdb_string(options)
        remarks = ['REMARK 950 WAVEDOCK REPAIR: COMPUTATIONAL ATOMS; REVIEW SIDECAR REPORT',
                   'REMARK 950 SOURCE SHA256 (1/2) '+source_hash[:32],
                   'REMARK 950 SOURCE SHA256 (2/2) '+source_hash[32:]]
        for a in added:
            remarks.append(f'REMARK 950 ADDED {a["chain"] or '(blank)'} {a["resi"]}{a["icode"]} {a["resn"]} {a["atom"]}')
        path.write_text('\n'.join(remarks)+'\n'+text,encoding='utf-8')
    else:
        doc = _mmcif_with_conect(st)
        block = doc.sole_block()
        if source_path is not None and Path(source_path).suffix.lower() in ('.cif','.mmcif'):
            original=gemmi.cif.read_file(str(source_path))
            original_block=original.sole_block()
            for category in ('_atom_site.','_atom_site_anisotrop.'):
                if category in block.get_mmcif_category_names():
                    original_block.set_mmcif_category(category,block.get_mmcif_category(category))
                elif category in original_block.get_mmcif_category_names():
                    original_block.find_mmcif_category(category).erase()
            doc,block=original,original_block
        if '_wavedock_repair_atom.' in block.get_mmcif_category_names():
            block.find_mmcif_category('_wavedock_repair_atom.').erase()
        block.set_pair('_wavedock_repair.source_sha256',source_hash)
        if added:
            loop = block.init_loop('_wavedock_repair_atom.',['chain','residue','insertion_code','name','origin'])
            for a in added: loop.add_row([gemmi.cif.quote(a['chain']),str(a['resi']),gemmi.cif.quote(a['icode'] or '.'),a['atom'],'computational'])
        doc.write_file(str(path))


def reconstructed_atom_indices(path, records=None):
    """Read explicit atom provenance from our output; never infer it from B factors."""
    import re
    path=Path(path)
    records=records if records is not None else atom_records(read_structure(path))
    indices=set()
    sidecar=path.with_suffix(path.suffix+'.repair.json')
    if sidecar.exists() and sidecar.stat().st_size<50_000_000:
        try:
            data=json.loads(sidecar.read_text(encoding='utf-8'))
            if data.get('exported_sha256',data.get('candidate_sha256'))==fingerprint(path):
                for a in data.get('added',[])+data.get('inherited',[]):
                    i=a['atom_index']
                    if 0<=i<len(records) and all(a.get(k)==records[i].get(k) for k in ('model','chain_index','residue_index','atom')): indices.add(i)
        except (OSError,ValueError,KeyError,TypeError): pass
    labels=[]
    if path.suffix.lower() in ('.cif','.mmcif'):
        block=gemmi.cif.read_file(str(path)).sole_block()
        table=block.find(['_wavedock_repair_atom.chain','_wavedock_repair_atom.residue',
                          '_wavedock_repair_atom.insertion_code','_wavedock_repair_atom.name','_wavedock_repair_atom.origin'])
        for row in table:
            values=[gemmi.cif.as_string(v) for v in row]
            if values[4]=='computational':
                try: labels.append((values[0],int(values[1]),values[2] if values[2] not in ('.','?') else '',values[3]))
                except ValueError: pass
    else:
        for line in path.read_text(encoding='utf-8',errors='replace').splitlines():
            if not line.startswith('REMARK 950 ADDED '): continue
            parts=line.split()
            if len(parts)!=7: continue
            match=re.fullmatch(r'(-?\d+)(.*)',parts[4])
            if match: labels.append(('' if parts[3]=='(blank)' else parts[3],int(match[1]),match[2],parts[6]))
    for chain,resi,icode,name in labels:
        hits=[a['atom_index'] for a in records if a['model']==0 and (a['chain'],a['resi'],a['icode'],a['atom'])==(chain,resi,icode,name)]
        if len(hits)==1: indices.add(hits[0])
    return sorted(indices)


def _identity(a):
    return (a['model'],a['chain_index'],a['residue_index'],a['atom'],a['altloc'])


def _work_structure(st, *, protein_only=False, segment_starts=()):
    """Use synthetic IDs only inside the backend; merge new atoms back by exact keys."""
    work = gemmi.Structure()
    model = gemmi.Model('1')
    mapping = {}
    for ci,chain in enumerate(st[0]):
        segment=0
        wc = gemmi.Chain(f'W{ci}_{segment}')
        previous=None
        for ri,res in enumerate(chain):
            if protein_only and res.name not in SIDECHAINS: continue
            if len(wc) and ((ci,ri) in segment_starts or
                    (previous is None or (previous.subchain and res.subchain and previous.subchain!=res.subchain) or previous.find_atom('C','*') is None or res.find_atom('N','*') is None or
                     previous.find_atom('C','*').pos.dist(res.find_atom('N','*').pos)>1.8)):
                model.add_chain(wc)
                segment+=1
                wc=gemmi.Chain(f'W{ci}_{segment}')
            clone = res.clone()
            clone.seqid = gemmi.SeqId(ri+1,' ')
            clone.subchain = wc.name
            wc.add_residue(clone)
            mapping[(wc.name,str(ri+1))]=(ci,ri)
            previous=res
        if len(wc): model.add_chain(wc)
    work.add_model(model)
    work.setup_entities()
    return work,mapping


def _retain_disulfides(fixer,mapping,st):
    """Carry explicit input S-S bonds into the isolated backend topology."""
    lookup={(*mapping[(a.residue.chain.id,a.residue.id)],a.name):a for a in fixer.topology.atoms()}
    present={frozenset((a.index,b.index)) for a,b in fixer.topology.bonds()}
    for connection in st.connections:
        if connection.type!=gemmi.ConnectionType.Disulf: continue
        endpoints=[]
        for partner in (connection.partner1,connection.partner2):
            matches=[(ci,ri,partner.atom_name) for ci,ch in enumerate(st[0]) for ri,r in enumerate(ch)
                     if ch.name==partner.chain_name and r.seqid==partner.res_id.seqid and r.find_atom(partner.atom_name,'*') is not None]
            if len(matches)!=1: raise ValueError('Ambiguous disulfide address; manual review required before repair.')
            endpoints.append(matches[0])
        a,b=endpoints
        aa=st[0][a[0]][a[1]].find_atom(a[2],'*')
        bb=st[0][b[0]][b[1]].find_atom(b[2],'*')
        if not 1.8<=aa.pos.dist(bb.pos)<=2.3:
            raise ValueError('Abnormal explicit disulfide geometry; manual review required.')
        if a in lookup and b in lookup:
            pair=frozenset((lookup[a].index,lookup[b].index))
            if pair not in present:
                fixer.topology.addBond(lookup[a],lookup[b]);present.add(pair)


def _backend_additions(st, changes, ph, progress, check_cancelled, segment_starts=()):
    progress(10,'Loading local protein repair templates')
    check_cancelled()
    import openmm
    from openmm import unit
    from pdbfixer import PDBFixer
    from Bio.PDB.Polypeptide import is_aa
    platform = openmm.Platform.getPlatformByName('CPU')
    platform.setPropertyDefaultValue('Threads','2')
    heavy = {(c.chain_index,c.residue_index):c for c in changes if c.kind=='heavy_atoms'}
    want_h = any(c.kind=='hydrogens' for c in changes)
    work,mapping = _work_structure(st,segment_starts=segment_starts)
    fixer = PDBFixer(pdbxfile=io.StringIO(work.make_mmcif_document().as_string()),platform=platform)
    _retain_disulfides(fixer,mapping,st)
    fixer.missingResidues = {}  # Never call findMissingResidues or replaceNonstandardResidues.
    fixer.findMissingAtoms()
    fixer.missingTerminals = {}  # No implicit terminal oxygens.
    fixer.missingAtoms = {r:[a for a in atoms if a.name in heavy[mapping[(r.chain.id,r.id)]].names]
        for r,atoms in fixer.missingAtoms.items() if mapping[(r.chain.id,r.id)] in heavy and is_aa(r.name,standard=True)}
    progress(30,'Rebuilding selected template atoms; original atoms fixed')
    if any(fixer.missingAtoms.values()): fixer.addMissingAtoms(seed=42)
    check_cancelled()
    additions = []
    for atom,pos in zip(fixer.topology.atoms(),fixer.positions.value_in_unit(unit.angstrom)):
        ci,ri = mapping[(atom.residue.chain.id,atom.residue.id)]
        res = st[0][ci][ri]
        if res.find_atom(atom.name,'*') is None:
            if (ci,ri) not in heavy or atom.name not in heavy[(ci,ri)].names:
                raise ValueError('Backend attempted an unapproved atom addition.')
            additions.append((ci,ri,atom.name,atom.element.symbol,tuple(pos)))
    _merge(st,additions)
    if want_h:
        # Never pass cofactors/nonstandard components to protonation routines.
        for ch in st[0]:
            for res in ch:
                if res.name in SIDECHAINS:
                    if any(a.altloc not in ('\0',' ') for a in res):
                        raise ValueError('Choose conformers for all protein residues before adding hydrogens.')
                    expected = set(('N CA C O '+SIDECHAINS[res.name]).split())
                    if not expected <= {a.name for a in res}:
                        raise ValueError('Complete all canonical heavy atoms before hydrogen addition.')
        work,mapping = _work_structure(st,protein_only=True,segment_starts=segment_starts)
        hfixer = PDBFixer(pdbxfile=io.StringIO(work.make_mmcif_document().as_string()),platform=platform)
        progress(55,'Adding protein hydrogens; non-protein components excluded from protonation')
        _retain_disulfides(hfixer,mapping,st)
        # Force-field hydrogen placement, fixed heavy atoms, reproducible initialisation.
        # The GUI runs this inside its isolated native worker process.
        from openmm.app import ForceField
        forcefield=ForceField('amber14/protein.ff14SB.xml','amber14/tip3p.xml')
        with _hydrogen_lock:
            state=random.getstate()
            try:
                random.seed(42)
                hfixer.addMissingHydrogens(ph,forcefield=forcefield)
            finally: random.setstate(state)
        check_cancelled()
        hs=[]
        for atom,pos in zip(hfixer.topology.atoms(),hfixer.positions.value_in_unit(unit.angstrom)):
            if atom.element.symbol!='H': continue
            ci,ri=mapping[(atom.residue.chain.id,atom.residue.id)]
            if st[0][ci][ri].find_atom(atom.name,'*') is None:
                hs.append((ci,ri,atom.name,'H',tuple(pos)))
        _merge(st,hs)
        additions.extend(hs)
        fixer,mapping = hfixer,mapping
    bonds=[]
    for a,b in fixer.topology.bonds():
        left=mapping.get((a.residue.chain.id,a.residue.id))
        right=mapping.get((b.residue.chain.id,b.residue.id))
        if left is not None and right is not None:
            bonds.append(((*left,a.name),(*right,b.name)))
    return additions,bonds


def _merge(st, additions):
    serial = max([a.serial for m in st for c in m for r in c for a in r]+[n for k,values in st.conect_map.items() for n in (k,*values)]+[0])
    for ci,ri,name,element,xyz in additions:
        serial+=1
        atom = gemmi.Atom()
        atom.name,atom.element,atom.pos = name,gemmi.Element(element),gemmi.Position(*xyz)
        atom.serial,atom.occ,atom.b_iso = serial,1.0,0.0
        st[0][ci][ri].add_atom(atom)


def _validate_additions(st, added_keys, bonds):
    """Check real merged geometry, including contacts with untouched cofactors."""
    atoms = {}
    for ci,ch in enumerate(st[0]):
        for ri,res in enumerate(ch):
            for a in res:
                # Repair targets have no unresolved conformers; retain all other conformers for contact checks below.
                atoms[(ci,ri,a.name,a.altloc)]=a
    graph=defaultdict(set)
    canonical={(ci,ri,name):a for (ci,ri,name,alt),a in atoms.items() if alt in ('\0',' ')}
    def label(key):
        ci,ri,name=key[:3]
        chain=st[0][ci];res=chain[ri]
        return f'{res.name} {chain.name or "(blank)"}:{res.seqid} {name}'
    blockers=[]
    for a,b in bonds:
        graph[a].add(b);graph[b].add(a)
        if not ({a,b}&added_keys) or a not in canonical or b not in canonical: continue
        aa,bb=canonical[a],canonical[b]
        d=aa.pos.dist(bb.pos)
        elements={aa.element.name,bb.element.name}
        low,high = ((0.8,1.55) if 'H' in elements else (1.5,2.25) if 'S' in elements else (1.0,1.85))
        if not math.isfinite(d) or not low<=d<=high:
            blockers.append(f'Abnormal reconstructed bond {label(a)} / {label(b)}: {d:.3f} A.')
    from itertools import combinations
    from Bio.PDB.vectors import Vector, calc_angle
    for center,neighbors in graph.items():
        if center not in canonical: continue
        for left,right in combinations(sorted(neighbors),2):
            if not {left,center,right}&added_keys or left not in canonical or right not in canonical: continue
            positions=[Vector(*tuple(canonical[k].pos)) for k in (left,center,right)]
            angle=math.degrees(calc_angle(*positions))
            if not math.isfinite(angle) or not 60<=angle<=175:
                blockers.append(f'Abnormal reconstructed bond angle {label(left)} / {label(center)} / {label(right)}: {angle:.2f} degrees.')
    for key in added_keys:
        if not graph[key]: blockers.append(f'Reconstructed atom has no bonded neighbours: {label(key)}.')
    radii={**VDW,'H':1.1}
    grid=defaultdict(list)
    for key,atom in atoms.items():
        if atom.occ<=0 or atom.element.name not in radii: continue
        if not all(math.isfinite(v) for v in (atom.pos.x,atom.pos.y,atom.pos.z)):
            if key[:3] in added_keys: blockers.append(f'Non-finite reconstructed coordinates: {label(key)}.')
            continue
        cell=tuple(math.floor(v/4) for v in (atom.pos.x,atom.pos.y,atom.pos.z))
        grid[cell].append((key,atom))
    for key in added_keys:
        atom=canonical[key]
        excluded=set(graph[key])|{key}
        for neighbor in graph[key]: excluded.update(graph[neighbor])
        p=tuple(math.floor(v/4) for v in (atom.pos.x,atom.pos.y,atom.pos.z))
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    for other,b in grid.get((p[0]+dx,p[1]+dy,p[2]+dz),()):
                        if other[:3] in excluded: continue
                        distance=atom.pos.dist(b.pos)
                        if atom.element.name=='S' and b.element.name=='S' and 1.8<=distance<=2.3: continue
                        # A normal donor-H...acceptor contact is not a clash.
                        if atom.element.name=='H' and b.element.name in ('N','O','S') and distance>=1.4:
                            parents=[canonical[p] for p in graph[key] if p in canonical and canonical[p].element.name in ('N','O','S')]
                            if any(math.degrees(gemmi.calculate_angle(parent.pos,atom.pos,b.pos))>=105 for parent in parents): continue
                        if radii[atom.element.name]+radii[b.element.name]-distance>0.8:
                            blockers.append(f'Severe reconstructed-atom overlap: {label(key)} / {label(other)}, {distance:.3f} A.')
                        if len(blockers)>=40: return blockers
    return list(dict.fromkeys(blockers))


def build_repairs(plan: RepairPlan, selected, output_dir, *, progress=lambda *_: None,
                  check_cancelled=lambda: None):
    """Build an isolated candidate. Accepting it is a separate controller action."""
    source=Path(plan.source)
    if fingerprint(source)!=plan.source_sha256: raise ValueError('Repair proposal is stale; generate a new preview.')
    choices={c.id:c for c in plan.changes}
    if not selected or set(selected)-choices.keys(): raise ValueError('Select valid proposed corrections first.')
    changes=[choices[key] for key in dict.fromkeys(selected)]
    st=read_structure(source)
    original=atom_records(st)
    origin_indices=reconstructed_atom_indices(source,original)
    def physical_key(a): return (a['model'],a['chain_index'],a['residue_index'],a['atom'],a['x'],a['y'],a['z'])
    inherited_keys={physical_key(original[i]) for i in origin_indices}
    removed=[]
    normalised=[]
    resolved=set()
    for c in changes:
        check_cancelled()
        if c.kind!='altloc': continue
        key=(c.chain_index,c.residue_index)
        if key in resolved: raise ValueError('Choose only one conformer for each residue.')
        resolved.add(key)
        res=st[0][c.chain_index][c.residue_index]
        if c.conformer not in {a.altloc for a in res}: raise ValueError('Conformer is no longer present.')
        normalised.extend(r for r in original if r['model']==0 and (r['chain_index'],r['residue_index'])==key and r['altloc']==c.conformer)
        removed.extend(r for r in original if r['model']==0 and (r['chain_index'],r['residue_index'])==key and r['altloc'] not in ('',c.conformer))
        for ai in reversed(range(len(res))):
            if res[ai].altloc not in ('\0',' ',c.conformer): del res[ai]
        for atom in res: atom.altloc='\0'
        if len({a.name for a in res})!=len(res): raise ValueError('Chosen conformer leaves duplicate atom names; manual review required.')
    for c in changes:
        if c.kind=='heavy_atoms' and any(a.altloc not in ('\0',' ') for a in st[0][c.chain_index][c.residue_index]):
            raise ValueError('Choose a conformer for this residue before rebuilding atoms.')
    expected_before=atom_records(st)
    additions,bonds=([],[])
    if any(c.kind in ('heavy_atoms','hydrogens') for c in changes):
        additions,bonds=_backend_additions(st,changes,plan.ph,progress,check_cancelled,plan.segment_starts)
    check_cancelled()
    added_keys={(ci,ri,name) for ci,ri,name,_,_ in additions}
    added=[a for a in atom_records(st) if a['model']==0 and (a['chain_index'],a['residue_index'],a['atom']) in added_keys]
    if not added and not removed and not normalised:
        raise ValueError('Selected changes produced no atom additions or alternate-atom removals.')
    directory=Path(output_dir)
    directory.mkdir(parents=True,exist_ok=True)
    suffix='.cif' if source.suffix.lower() in ('.cif','.mmcif') else '.pdb'
    candidate=directory/f'protein_{uuid4().hex[:12]}_repaired{suffix}'
    inherited=[a for a in atom_records(st) if physical_key(a) in inherited_keys]
    _write_structure(st,candidate,added+inherited,plan.source_sha256,source)
    saved=read_structure(candidate)
    surviving=[a for a in atom_records(saved) if not (a['model']==0 and (a['chain_index'],a['residue_index'],a['atom']) in added_keys)]
    if len(surviving)!=len(expected_before) or any(
        _identity(a)!=_identity(b) or any(abs(a[k]-b[k])>1e-5 for k in ('x','y','z'))
        for a,b in zip(expected_before,surviving)):
        raise ValueError('Candidate changed existing atom identities or coordinates; repair rejected.')
    blockers=_validate_additions(saved,added_keys,bonds)
    after=check_protein(candidate,progress=lambda p,m:progress(70+p*.3,m),check_cancelled=check_cancelled)
    for issue in after.issues:
        if issue.code=='missing_atoms' and any((a['chain_index'],a['residue_index']) in
            {(c.chain_index,c.residue_index) for c in changes if c.kind=='heavy_atoms'} for a in issue.atoms):
            blockers.append('Selected residue remains incomplete: '+issue.location)
    def issue_key(issue):
        # Connectivity depends on C/N endpoints, not newly completed side chains.
        atoms = [a for a in issue.atoms if a['atom'] in ('C','N')] if issue.code=='backbone_connectivity' else issue.atoms
        return (issue.code,tuple(sorted({(a['model'],a['chain_index'],a['residue_index'],a['atom'],a['x'],a['y'],a['z']) for a in atoms})))
    previous_issues={issue_key(i) for i in plan.before.issues}
    modified=resolved|{(ci,ri) for ci,ri,_,_,_ in additions}
    for issue in after.issues:
        if issue.code in ('steric_overlap','bond_geometry','angle_geometry','backbone_connectivity','invalid_coordinates') and issue_key(issue) not in previous_issues and any((a['chain_index'],a['residue_index']) in modified for a in issue.atoms):
            blockers.append('New geometry finding in repaired region: '+issue.location+' / '+issue.code)
    if fingerprint(source)!=plan.source_sha256: raise ValueError('Source changed while building; candidate was not activated.')
    result=RepairResult(plan,list(selected),str(candidate),after,added,removed,blockers,[
        'Review all computationally reconstructed atoms. Existing warnings may remain.',
        'No loops, Ramachandran/rotamer validation, charge assignment or general force-field validation were performed.',
        'Hydrogen placement uses Amber14/ff14SB, seed 42 and default residue pKas, not PROPKA. Fixed heavy atoms, metal coordination and binding-site protonation need review.'
    ],fingerprint(candidate),inherited,normalised)
    result.export_report(candidate.with_suffix(candidate.suffix+'.repair.json'))
    progress(100,'Candidate ready for review' if result.can_accept else 'Candidate blocked by geometry checks')
    return result


def export_repaired_structure(result: RepairResult, target):
    """Export a validated candidate/revision with its full provenance sidecar."""
    import shutil
    target=Path(target)
    candidate=Path(result.candidate)
    if not result.can_accept: raise ValueError('Geometry-blocked candidates cannot be exported as repaired structures.')
    if fingerprint(candidate)!=result.candidate_sha256: raise ValueError('Candidate changed; rebuild before exporting.')
    if target.resolve() in (candidate.resolve(),Path(result.plan.source).resolve()):
        raise ValueError('Choose a separate export path; original/revision files are immutable.')
    if target.suffix.lower() not in ('.pdb','.cif','.mmcif'): raise ValueError('Choose PDB or mmCIF format.')
    if target.suffix.lower()==candidate.suffix.lower(): shutil.copyfile(candidate,target)
    else:
        structure=read_structure(candidate)
        if target.suffix.lower()=='.pdb' and candidate.suffix.lower() in ('.cif','.mmcif'):
            _pdb_bonds_from_mmcif(structure,candidate)
        _write_structure(structure,target,result.added+result.inherited,result.plan.source_sha256,candidate)
    report_path=target.with_suffix(target.suffix+'.repair.json')
    result.export_report(report_path)
    data=json.loads(report_path.read_text(encoding='utf-8'))
    data['exported_structure']=str(target)
    data['exported_sha256']=fingerprint(target)
    report_path.write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')
