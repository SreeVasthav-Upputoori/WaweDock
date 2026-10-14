"""Protein-ligand interaction fingerprinting (ProLIF) and a 2D interaction diagram.

Classifies contacts into named interaction types (hydrogen bond, hydrophobic, pi-stacking,
salt bridge, etc.) using ProLIF's geometric+chemical rules (distance *and* angle, not distance
alone). This is still a **calculated** classification from a docked pose, not an experimentally
verified interaction — see :data:`CAVEAT`.

Building the protein as an RDKit molecule
------------------------------------------
ProLIF needs correct bonds/aromaticity on the protein. RDKit's own PDB reader assigns bonds by
distance only and frequently mis-kekulizes rings; MDAnalysis's RDKit converter fails outright on
some residues; Gemmi's own monomer-library bond templates (used elsewhere in this app) only cover
the standard residue names and not PDB2PQR's protonation-variant names (HID/HIE/ASH/CYX/...), so
they cannot be reused here without a much larger bundled template set. Open Babel's residue-aware
bond perception handles every variant reliably, so it builds the bonds/aromaticity (as a MOL2
block), while residue names/numbers/chains are re-attached from Gemmi's parse of the *same*
file — the two readers visit atoms in the same, file order (verified on real structures), so
atoms are matched by index.

Open Babel is GPL-2.0 licensed. It is therefore an *optional* dependency
(``pip install wavedock[interactions]``) and is never bundled into the distributed executable
(see ``wavedock.spec``'s ``excludes``) — this feature is simply unavailable in that build unless
Open Babel is installed separately into the same environment.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
from rdkit import Chem

from mds.chem.report import PreparationError

CAVEAT = ("Interaction types are calculated from a docked pose using ProLIF's geometric and "
          "chemical rules (distance and angle, not distance alone). They are classifications of "
          "a predicted pose, not experimentally verified interactions.")

_INTERACTION_LABELS = {
    "HBDonor": "H-bond (ligand donor)", "HBAcceptor": "H-bond (ligand acceptor)",
    "Hydrophobic": "Hydrophobic contact", "PiStacking": "Pi-stacking", "PiCation": "Cation-pi",
    "CationPi": "Pi-cation", "Cationic": "Salt bridge (ligand cation)",
    "Anionic": "Salt bridge (ligand anion)", "XBDonor": "Halogen bond (ligand donor)",
    "XBAcceptor": "Halogen bond (ligand acceptor)", "MetalDonor": "Metal coordination",
    "MetalAcceptor": "Metal coordination", "VdWContact": "Van der Waals contact",
    "EdgeToFace": "Pi-stacking (edge-to-face)", "FaceToFace": "Pi-stacking (face-to-face)",
}


@dataclass
class ResidueInteractions:
    residue: str                         # e.g. "ASP25:A"
    interactions: list[str] = field(default_factory=list)   # human-readable labels


@dataclass
class InteractionResult:
    residues: list[ResidueInteractions]
    network_html: str | None             # standalone HTML page for the 2D diagram, or None
    n_interactions: int


def protein_to_rdkit(pdb_path: Path) -> Chem.Mol:
    """Build an RDKit molecule for a prepared receptor PDB with correct bonds and residue labels."""
    try:
        from openbabel import pybel
    except ImportError as exc:
        raise PreparationError(
            "Interaction analysis needs Open Babel, which is not installed. It is GPL-licensed "
            "and kept optional/separate from this program; install it with "
            "'pip install openbabel-wheel' (or 'pip install wavedock[interactions]') into this "
            "same Python environment, then try again. It is not available in the packaged .exe."
        ) from exc

    pdb_path = Path(pdb_path)
    st = gemmi.read_structure(str(pdb_path))
    labels = [(ch.name, res.name, res.seqid.num) for ch in st[0] for res in ch for _ in res]

    ob_mol = next(pybel.readfile("pdb", str(pdb_path)), None)
    if ob_mol is None:
        raise PreparationError(f"Open Babel could not read {pdb_path.name}")
    mol = Chem.MolFromMol2Block(ob_mol.write("mol2"), removeHs=False, sanitize=False)
    if mol is None or mol.GetNumAtoms() != len(labels):
        raise PreparationError(f"Atom count mismatch while building the interaction model for "
                               f"{pdb_path.name} ({'parse failed' if mol is None else f'{mol.GetNumAtoms()} vs {len(labels)}'})")
    # Skip kekulization: Open Babel already assigned correct aromaticity; RDKit's own
    # kekulization occasionally fails on unusual ring systems even when the structure is valid.
    flags = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
    problems = Chem.SanitizeMol(mol, sanitizeOps=flags, catchErrors=True)
    if problems:
        raise PreparationError(f"Could not sanitise the protein model for interaction analysis "
                               f"({problems!r})")

    rw = Chem.RWMol(mol)
    for i, (chain, resname, resnum) in enumerate(labels):
        info = Chem.AtomPDBResidueInfo()
        info.SetChainId(chain)
        info.SetResidueName(f"{resname:<3}")
        info.SetResidueNumber(int(resnum))
        info.SetName(f" {mol.GetAtomWithIdx(i).GetSymbol():<3}")
        rw.GetAtomWithIdx(i).SetMonomerInfo(info)
    return rw.GetMol()


def compute_interactions(protein_rdkit: Chem.Mol, pose: Chem.Mol, *, diagram: bool = True) -> InteractionResult:
    """Run a ProLIF fingerprint for one docked pose against the prepared receptor."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import prolif
        from prolif.molecule import Molecule

    prot_mol = Molecule(protein_rdkit)
    lig_mol = Molecule.from_rdkit(Chem.Mol(pose))

    fp = prolif.Fingerprint()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # n_jobs=1: a single pose is fast enough sequentially, and this method may run inside a
        # worker thread (QThreadPool), where spawning a multiprocessing pool is unsafe on Windows.
        fp.run_from_iterable([lig_mol], prot_mol, residues="all", progress=False, n_jobs=1)

    by_residue: dict[str, list[str]] = {}
    n = 0
    if fp.ifp:
        for data in fp.ifp[0].interactions():
            key = _format_residue(data.protein)
            human = _INTERACTION_LABELS.get(data.interaction, data.interaction)
            by_residue.setdefault(key, [])
            if human not in by_residue[key]:
                by_residue[key].append(human)
            n += 1

    residues = [ResidueInteractions(res, sorted(types)) for res, types in sorted(by_residue.items())]

    html = None
    if diagram and fp.ifp:
        try:
            html = _build_network_html(fp, lig_mol, prot_mol)
        except Exception:  # noqa: BLE001 — diagram is a bonus; never fail the analysis for it
            html = None
    return InteractionResult(residues, html, n)


def _format_residue(prot_res) -> str:
    # ProLIF ResidueId string form is like "ASP25.A" -> normalise to "ASP25:A"
    text = str(prot_res)
    return text.replace(".", ":") if ":" not in text else text


def _build_network_html(fp, lig_mol, _prot_mol) -> str:
    from prolif.plotting.network import LigNetwork

    # prolif.molecule.Molecule IS an rdkit.Chem.Mol subclass; pass it directly.
    net = LigNetwork.from_fingerprint(fp, lig_mol, kind="frame", frame=0, display_all=False)
    html = net._get_html()
    # ProLIF's template loads vis-network's JS and CSS from unpkg.com; replace both with the
    # vendored copies (mds/assets/viewer/vis-network.min.{js,css}) so the diagram renders with
    # no network access. The caller sets the QWebEngineView base URL to that same folder.
    html = re.sub(r'src="https://unpkg\.com/vis-network[^"]*\.js"', 'src="vis-network.min.js"', html)
    html = re.sub(r'href="https://unpkg\.com/vis-network[^"]*\.css"', 'href="vis-network.min.css"', html)
    return html
