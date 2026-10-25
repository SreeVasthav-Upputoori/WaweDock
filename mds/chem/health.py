"""Read-only protein health screening; no coordinate changes or inferred repairs.

Conservative geometric screens are not a substitute for MolProbity validation.
Atom identities use the same Gemmi traversal as the editable molecular viewer.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass, field
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
from typing import Callable

import gemmi

from mds.chem.disulfides import find_disulfides
from mds.chem.edit import atom_records
from mds.chem.ramachandran import compute_ramachandran, summarize as summarize_ramachandran
from mds.chem.receptor import read_structure
from mds.chem.rotamers import compute_rotamers, summarize as summarize_rotamers

# Canonical heavy atoms only. OXT is deliberately optional at chain termini.
SIDECHAINS = {
    'ALA': 'CB', 'ARG': 'CB CG CD NE CZ NH1 NH2', 'ASN': 'CB CG OD1 ND2',
    'ASP': 'CB CG OD1 OD2', 'CYS': 'CB SG', 'GLN': 'CB CG CD OE1 NE2',
    'GLU': 'CB CG CD OE1 OE2', 'GLY': '', 'HIS': 'CB CG ND1 CD2 CE1 NE2',
    'ILE': 'CB CG1 CG2 CD1', 'LEU': 'CB CG CD1 CD2', 'LYS': 'CB CG CD CE NZ',
    'MET': 'CB CG SD CE', 'PHE': 'CB CG CD1 CD2 CE1 CE2 CZ', 'PRO': 'CB CG CD',
    'SER': 'CB OG', 'THR': 'CB OG1 CG2', 'TRP': 'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR': 'CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL': 'CB CG1 CG2',
}
# Explicit standard topology avoids mistaking bonded and 1-3 pairs for clashes.
SIDE_BONDS = {
    'ALA': '', 'ARG': 'CB-CG CG-CD CD-NE NE-CZ CZ-NH1 CZ-NH2',
    'ASN': 'CB-CG CG-OD1 CG-ND2', 'ASP': 'CB-CG CG-OD1 CG-OD2',
    'CYS': 'CB-SG', 'GLN': 'CB-CG CG-CD CD-OE1 CD-NE2',
    'GLU': 'CB-CG CG-CD CD-OE1 CD-OE2', 'GLY': '',
    'HIS': 'CB-CG CG-ND1 CG-CD2 ND1-CE1 CE1-NE2 NE2-CD2',
    'ILE': 'CB-CG1 CB-CG2 CG1-CD1', 'LEU': 'CB-CG CG-CD1 CG-CD2',
    'LYS': 'CB-CG CG-CD CD-CE CE-NZ', 'MET': 'CB-CG CG-SD SD-CE',
    'PHE': 'CB-CG CG-CD1 CG-CD2 CD1-CE1 CD2-CE2 CE1-CZ CE2-CZ',
    'PRO': 'CB-CG CG-CD CD-N', 'SER': 'CB-OG', 'THR': 'CB-OG1 CB-CG2',
    'TRP': 'CB-CG CG-CD1 CG-CD2 CD1-NE1 NE1-CE2 CE2-CD2 CD2-CE3 CE3-CZ3 CZ3-CH2 CH2-CZ2 CZ2-CE2',
    'TYR': 'CB-CG CG-CD1 CG-CD2 CD1-CE1 CD2-CE2 CE1-CZ CE2-CZ CZ-OH',
    'VAL': 'CB-CG1 CB-CG2',
}
VDW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80}


@dataclass
class HealthIssue:
    code: str
    severity: str
    location: str
    message: str
    recommendation: str
    atoms: list[dict] = field(default_factory=list)
    status: str = 'Needs review'


@dataclass
class HealthReport:
    source: str
    source_sha256: str
    issues: list[HealthIssue] = field(default_factory=list)
    coverage: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    backend: str = field(default_factory=lambda: f'Gemmi {gemmi.__version__}; WaweDock geometric screening v1')
    ramachandran: list = field(default_factory=list)
    ramachandran_summary: dict = field(default_factory=dict)
    disulfides: list = field(default_factory=list)
    free_cysteines: list = field(default_factory=list)
    rotamers: list = field(default_factory=list)
    rotamers_summary: dict = field(default_factory=dict)

    def export(self, path: Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2, allow_nan=False), encoding='utf-8')

    @property
    def summary(self) -> str:
        warnings = sum(i.severity == 'Warning' for i in self.issues)
        errors = sum(i.severity == 'Error' for i in self.issues)
        return f'{errors} error(s), {warnings} warning(s), {len(self.issues) - errors - warnings} note(s). Review screening limits.'


def _compatible(a, b):
    return not a.altloc.strip('\0 ') or not b.altloc.strip('\0 ') or a.altloc == b.altloc


def check_protein(source: Path, *, progress: Callable = lambda *_: None,
                  check_cancelled: Callable = lambda: None) -> HealthReport:
    """Screen the first model, preserving input bytes and all structural components.

    Sequence gaps require a complete entity sequence and consistent label_seq
    alignment. Numbering jumps alone never establish missing residues.
    """
    source = Path(source)
    digest = sha256(source.read_bytes()).hexdigest()
    st = read_structure(source)
    records = atom_records(st)
    for record in records:
        for axis in ("x", "y", "z"):
            if not math.isfinite(record[axis]): record[axis] = None
    report = HealthReport(str(source), digest)
    report.coverage = [
        'All chains in the first model: canonical amino-acid heavy-atom completeness (N, CA, C, O and side chains; OXT optional).',
        'Within-segment C-N connectivity and broad backbone bond-length / N-CA-C angle screens.',
        'Severe heavy-atom overlap screen: van der Waals overlap > 0.8 A; excludes topology neighbours through two bonds, disulfides, metals and incompatible altlocs.',
        'Alternate conformers, protein hydrogen presence and nonstandard/component inventory.',
        'Unmodelled sequence regions only when entity sequence and label_seq alignment agree.',
        'Backbone phi/psi Ramachandran region screen (favored/allowed/outlier) for general, glycine and proline residues.',
        'Disulfide bonds from cysteine SG-SG distance (1.8-2.5 A) and an inventory of free cysteines.',
        'Side-chain rotamer screen: sp3 chi torsions versus ideal staggered positions (favored/allowed/outlier).',
    ]
    report.limitations = [
        'Geometric screening is not a full stereochemical validation or a docking-readiness certificate.',
        'The Ramachandran screen uses fixed, generous phi/psi region boxes, not MolProbity percentile contours; borderline conformations are not outliers and pre-proline is treated as general.',
        'Disulfides are assigned by SG-SG distance only; oxidation-state labels, cross-link chemistry and non-Cys covalent links are not verified.',
        'The rotamer screen compares only sp3 chi torsions to ideal staggered angles (+/-60/180); it is not a rotamer-library or MolProbity percentile analysis, and planar/aromatic terminal chi angles and proline are not assessed.',
        'Hydrogen-bond, charge, pKa and force-field validation were not performed.',
        'No crystallographic symmetry contacts, complete ligand topology or hydrogen clash analysis. Metal coordination is not validated.',
        'Sequence absence describes unmodelled coordinates; experimental disorder cannot be inferred from coordinates alone.',
        'No atoms or residues were rebuilt, minimised, removed or moved. All reported coordinates come from the input; their experimental/modelled provenance is unknown.',
    ]
    if len(st) > 1:
        report.limitations.append(f'Only model 1 of {len(st)} was screened; other models were not analysed.')
    supplied_labels = {(ci,ri): res.label_seq for ci,ch in enumerate(st[0]) for ri,res in enumerate(ch)}
    # Align only against supplied entity sequences, never invent a sequence.
    try:
        st.assign_label_seq_id()
    except (RuntimeError, ValueError):
        report.limitations.append('Sequence alignment unavailable; missing sequence regions were not classified.')
    segment_starts = set()
    if st.input_format == gemmi.CoorFormat.Pdb:
        # The default parser merges consecutive same-ID chains across TER.
        # Read a second segmentation only for boundaries; never replace the
        # canonical atom records used by the viewer and edit history.
        split = gemmi.read_pdb(str(source), split_chain_on_ter=True)
        split_records = atom_records(split)
        identity_keys = ("model", "chain", "resi", "icode", "resn", "atom", "altloc")
        if len(split_records) == len(records) and all(
                all(a[k] == b[k] for k in identity_keys) for a,b in zip(records,split_records)):
            previous_segment = None
            for atom in split_records:
                segment = (atom["model"],atom["chain_index"])
                if segment != previous_segment: segment_starts.add(atom["atom_index"])
                previous_segment = segment
        else:
            report.limitations.append('TER segmentation could not be mapped to canonical atoms; backbone connectivity screening was omitted.')
            segment_starts = None
    index, residues, heavy, graph = 0, [], [], defaultdict(set)
    protein_h = protein_heavy = 0
    issues_by_residue = {}

    def add(code, severity, location, message, recommendation, atoms=()):
        # Keep UI/report memory bounded for pathological input, retaining all atom selectors of each issue.
        if len(report.issues) < 5000:
            report.issues.append(HealthIssue(code, severity, location, message, recommendation, list(atoms)))
        elif not any('5000' in item for item in report.limitations):
            report.limitations.append('Issue list truncated at 5000 entries; review the original structure with an external validator.')

    progress(5, 'Inspecting residues and atom completeness')
    for ci, chain in enumerate(st[0]):
        for ri, res in enumerate(chain):
            check_cancelled()
            info = gemmi.find_tabulated_residue(res.name)
            protein = bool(info and info.is_amino_acid())
            location = f'{chain.name or "(blank)"}:{res.seqid} {res.name}'
            entries = [(index + ai, atom) for ai, atom in enumerate(res)]
            selectors = records[index:index+len(res)]
            index += len(res)
            residues.append((chain, res, entries, selectors, protein))
            valid = [(i, a) for i, a in entries if all(math.isfinite(v) for v in (a.pos.x, a.pos.y, a.pos.z))]
            if len(valid) != len(entries):
                add('invalid_coordinates', 'Error', location, 'Non-finite atomic coordinates.', 'Correct or replace the source model.', selectors)
            altlocs = sorted({a.altloc for _, a in entries if a.altloc not in ('\0', ' ')})
            if altlocs:
                add('alternate_conformations', 'Warning', location, 'Alternate conformations: ' + ', '.join(altlocs),
                    'Review occupancy and choose a consistent conformer before preparation.', selectors)
            if protein:
                protein_h += sum(a.is_hydrogen() for _, a in entries)
                protein_heavy += sum(not a.is_hydrogen() for _, a in entries)
                if res.name in SIDECHAINS:
                    # Check every conformer independently, including shared blank atoms.
                    shared = {a.name for _, a in entries if a.altloc in ('\0', ' ')}
                    expected = set(('N CA C O ' + SIDECHAINS[res.name]).split())
                    for alt in altlocs or ['']:
                        names = shared | {a.name for _, a in entries if a.altloc == alt}
                        missing = sorted(expected - names)
                        if missing:
                            add('missing_atoms', 'Warning', location,
                                ('Conformer ' + alt + ': ' if alt else '') + 'Missing heavy atoms: ' + ', '.join(missing),
                                'Rebuild using a validated residue template; inspect the result.', selectors)
                    bonds = 'N-CA CA-C C-O C-OXT ' + ('CA-CB ' if res.name != 'GLY' else '') + SIDE_BONDS[res.name]
                    for pair in bonds.split():
                        left, right = pair.split('-')
                        for i, a in valid:
                            if a.name != left: continue
                            for j, b in valid:
                                if b.name == right and _compatible(a, b):
                                    graph[i].add(j); graph[j].add(i)
                    for left, right, low, high in [('N','CA',1.25,1.65), ('CA','C',1.30,1.75), ('C','O',1.05,1.45)]:
                        for i, a in valid:
                            if a.name != left: continue
                            for j, b in valid:
                                if b.name == right and _compatible(a,b):
                                    distance = a.pos.dist(b.pos)
                                    if not low <= distance <= high:
                                        add('bond_geometry', 'Warning', location,
                                            f'{left}-{right} distance {distance:.2f} A is outside the broad {low}-{high} A screen.',
                                            'Inspect coordinates and connectivity; consider restrained refinement.', [records[i],records[j]])
                    for _, a in valid:
                        if a.name != 'CA': continue
                        for _, n in valid:
                            if n.name != 'N' or not _compatible(a,n): continue
                            for _, c in valid:
                                if c.name != 'C' or not _compatible(a,c) or not _compatible(n,c): continue
                                angle = math.degrees(gemmi.calculate_angle(n.pos,a.pos,c.pos))
                                if math.isfinite(angle) and not 95 <= angle <= 125:
                                    add('angle_geometry','Warning',location,f'N-CA-C angle {angle:.1f} degrees outside broad 95-125 degree screen.',
                                        'Inspect backbone geometry; consider restrained refinement.',selectors)
                else:
                    add('nonstandard_residue','Note',location,'Nonstandard amino acid; canonical atom templates were not applied.',
                        'Preserve unless explicitly authorised; review chemistry and downstream parameters.',selectors)
            elif not res.is_water():
                add('structural_component','Note',location,'Non-protein component retained; its chemistry and geometry were not validated.',
                    'Review essential cofactors, ligands and ions before docking preparation.',selectors)
            for i, a in valid:
                if not a.is_hydrogen() and a.occ > 0 and a.element.name in VDW:
                    heavy.append((i,a,ci,ri,res.name,protein))
            issues_by_residue[(ci,ri)] = (location, selectors)
    if not protein_heavy:
        add('no_protein','Error','Structure','No recognised amino-acid heavy atoms found.','Load a protein structure.')
    elif protein_h == 0:
        add('hydrogens','Note','Protein','No protein hydrogen atoms are present; common for experimental structures.',
            'Assign protonation and add hydrogens during preparation at the intended pH.')
    else:
        report.limitations.append(f'{protein_h} protein hydrogen atom(s) found; hydrogen completeness and protonation states were not validated.')
    progress(35, 'Checking backbone segments and sequence evidence')
    # Use chain indices: Gemmi can create new Python wrappers for the same chain.
    for ci, chain in enumerate(st[0]):
        poly = [(res,entries,selectors) for ch,res,entries,selectors,protein in residues
                if protein and selectors and selectors[0]['chain_index'] == ci]
        for (a,aa,asel),(b,bb,bsel) in zip(poly,poly[1:]):
            check_cancelled()
            # TER/subchain boundaries and separate chains are legitimate termini.
            if a.subchain and b.subchain and a.subchain != b.subchain: continue
            if segment_starts is None or (bsel and bsel[0]["atom_index"] in segment_starts): continue
            pairs = [(i,c,j,n) for i,c in aa for j,n in bb if c.name=='C' and n.name=='N' and _compatible(c,n)]
            if not pairs: continue  # missing backbone already reported
            distances = [c.pos.dist(n.pos) for _,c,_,n in pairs]
            for i,c,j,n in pairs:
                if c.pos.dist(n.pos) <= 1.8:
                    graph[i].add(j); graph[j].add(i)
            if all(math.isfinite(d) for d in distances) and (min(distances)>1.8 or max(distances)<1.1):
                add('backbone_connectivity','Warning',f'{chain.name or "(blank)"}:{a.seqid}-{b.seqid}',
                    f'Peptide C-N distances {min(distances):.2f}-{max(distances):.2f} A. A break, unresolved region or unusual geometry needs review; numbering alone does not establish missing residues.',
                    'Review sequence/termination metadata. Use sequence-guided modelling for genuine gaps; never simply join distant atoms.',asel+bsel)
        # Only entity-backed, fully consistent sequence alignment supports missing regions.
        if not poly: continue
        entity = st.get_entity_of(chain.get_polymer())
        sequence = list(entity.full_sequence) if entity else []
        mapped = [(r.label_seq,r.name) for r,_,_ in poly]
        consistent = bool(sequence and mapped) and all(
            n is not None and 1 <= n <= len(sequence) and name in sequence[n-1].split(',') for n,name in mapped)
        # Sequence-only alignments must be unique. Repeated motifs can otherwise
        # place an unresolved segment on the wrong side of an observed fragment.
        authoritative = all(supplied_labels[(ci,sel[0]["residue_index"])] is not None
                            for _,_,sel in poly if sel)
        if consistent and not authoritative:
            forward, backward = [], []
            cursor = 0
            for _,name in mapped:
                while cursor < len(sequence) and name not in sequence[cursor].split(','): cursor += 1
                forward.append(cursor); cursor += 1
            cursor = len(sequence)-1
            for _,name in reversed(mapped):
                while cursor >= 0 and name not in sequence[cursor].split(','): cursor -= 1
                backward.append(cursor); cursor -= 1
            consistent = forward == list(reversed(backward))
        if not consistent or len({n for n,_ in mapped}) != len(mapped):
            report.limitations.append(f'Chain {chain.name or "(blank)"}: no unambiguous full-sequence alignment; missing sequence regions were not classified.')
            continue
        observed = {n for n,_ in mapped}
        absent = sorted(set(range(1,len(sequence)+1))-observed)
        runs = []
        for n in absent:
            if runs and n == runs[-1][-1]+1: runs[-1].append(n)
            else: runs.append([n])
        for run in runs:
            neighbours = [sel for r,_,sel in poly if r.label_seq in (run[0]-1,run[-1]+1)]
            terminal = run[0]==1 or run[-1]==len(sequence)
            add('unmodelled_sequence','Note' if terminal else 'Warning',f'{chain.name or "(blank)"}:sequence {run[0]}-{run[-1]}',
                f'{len(run)} entity-sequence residue(s) have no modelled coordinates' + (' at a sequence terminus.' if terminal else '.'),
                'Review experimental metadata and target sequence; reconstruction is computational and requires independent validation.',
                [item for group in neighbours for item in group])
    progress(55, 'Screening severe nonbonded heavy-atom overlaps')
    cells = defaultdict(list)
    clashes = defaultdict(list)
    offsets = list(product((-1,0,1),repeat=3))
    for step,(i,a,ci,ri,resn,protein) in enumerate(heavy):
        if step % 256 == 0:
            check_cancelled(); progress(55 + 40*step/max(len(heavy),1),'Screening severe nonbonded heavy-atom overlaps')
        cell = tuple(math.floor(v/4.0) for v in (a.pos.x,a.pos.y,a.pos.z))
        excluded = set(graph[i])
        for neighbour in graph[i]: excluded.update(graph[neighbour])
        for dx,dy,dz in offsets:
            for j,b,cj,rj,bresn,bprotein in cells.get((cell[0]+dx,cell[1]+dy,cell[2]+dz), ()):
                if j in excluded or not (protein or bprotein) or not _compatible(a,b): continue
                # Ligand intramolecular topology is not available; skip those pairs.
                if ci==cj and ri==rj and resn not in SIDECHAINS: continue
                distance = a.pos.dist(b.pos)
                if a.element.name=='S' and b.element.name=='S' and 1.8 <= distance <= 2.3: continue
                overlap = VDW[a.element.name]+VDW[b.element.name]-distance
                if overlap > 0.8:
                    key = tuple(sorted(((ci,ri),(cj,rj))))
                    if len(clashes[key]) < 3: clashes[key].append((i,j,distance,overlap))
        cells[cell].append((i,a,ci,ri,resn,protein))
    for keys,pairs in clashes.items():
        left,right = (issues_by_residue[k][0] for k in keys)
        i,j,distance,overlap = max(pairs,key=lambda p:p[3])
        add('steric_overlap','Warning',f'{left} / {right}',
            f'Severe overlap candidate: {records[i]["atom"]}-{records[j]["atom"]}, {distance:.2f} A, van der Waals overlap {overlap:.2f} A.',
            'Inspect alternate conformers, covalent links and ligand chemistry before any local refinement.',[records[k] for pair in pairs for k in pair[:2]])
    check_cancelled()
    progress(96, 'Screening backbone torsions (Ramachandran region screen)')
    report.ramachandran = compute_ramachandran(st, check_cancelled=check_cancelled)
    report.ramachandran_summary = summarize_ramachandran(report.ramachandran)
    residue_selectors = {(ch.name or '', res.seqid.num, (res.seqid.icode or '').strip()): sel
                         for ch, res, _, sel, protein in residues if protein}
    for point in report.ramachandran:
        if point.classification != 'outlier': continue
        location = f'{point.chain or "(blank)"}:{point.seqid}{point.icode} {point.resname}'
        add('ramachandran_outlier', 'Warning', location,
            f'Backbone torsion outlier: phi {point.phi:.0f}, psi {point.psi:.0f} degrees, outside the allowed-region screen.',
            'Review backbone geometry; consider restrained refinement. The region screen is coarse and is not a substitute for MolProbity.',
            residue_selectors.get((point.chain, point.seqid, point.icode), []))
    progress(98, 'Detecting disulfide bonds')
    report.disulfides, report.free_cysteines = find_disulfides(st, check_cancelled=check_cancelled)
    for ss in report.disulfides:
        location = (f'{ss.chain_a or "(blank)"}:{ss.seqid_a}{ss.icode_a} - '
                    f'{ss.chain_b or "(blank)"}:{ss.seqid_b}{ss.icode_b}')
        selectors = (residue_selectors.get((ss.chain_a, ss.seqid_a, ss.icode_a), [])
                     + residue_selectors.get((ss.chain_b, ss.seqid_b, ss.icode_b), []))
        add('disulfide_bond', 'Note', location,
            f'Disulfide bond: cysteine SG-SG distance {ss.distance:.2f} A.',
            'Preserve this bridge during preparation; treat both cysteines as oxidised (no added thiol hydrogen).',
            selectors)
    if report.free_cysteines:
        free_selectors = [item for fc in report.free_cysteines
                          for item in residue_selectors.get((fc['chain'], fc['seqid'], fc['icode']), [])]
        add('free_cysteine', 'Note', 'Protein',
            f'{len(report.free_cysteines)} free cysteine(s) not in a disulfide bond within the 1.8-2.5 A window.',
            'Review reduced-thiol protonation and possible metal or ligand coordination before preparation.',
            free_selectors)
    progress(99, 'Screening side-chain rotamers')
    report.rotamers = compute_rotamers(st, check_cancelled=check_cancelled)
    report.rotamers_summary = summarize_rotamers(report.rotamers)
    for rot in report.rotamers:
        if rot.classification != 'outlier': continue
        location = f'{rot.chain or "(blank)"}:{rot.seqid}{rot.icode} {rot.resname}'
        chi_text = ', '.join(f'{value:.0f}' for value in rot.chis)
        add('rotamer_outlier', 'Warning', location,
            f'Side-chain rotamer outlier: chi angles [{chi_text}] degrees are far from ideal staggered positions.',
            'Review the side-chain conformation; consider a rotamer correction. The screen checks staggered sp3 chi torsions only, not a rotamer library.',
            residue_selectors.get((rot.chain, rot.seqid, rot.icode), []))
    if sha256(source.read_bytes()).hexdigest() != digest:
        raise RuntimeError('Structure changed while screening; run Check Protein again.')
    progress(100,'Protein screening complete; review findings and limits')
    return report
