"""Ligand reading, validation, 3D generation and PDBQT preparation (RDKit + Meeko).

Scientific policy
-----------------
* The protonation/tautomeric state in the input is **respected**.  Common ionisable
  groups are detected and reported so the user can judge whether the input state is
  appropriate for pH 7.4.  An optional ``ionize="simple"`` mode applies a small set of
  textbook rules (carboxylic acids, aliphatic amines, amidines/guanidines); it is not
  a pKa predictor and says so in the report.
* Undefined stereocentres are never silently resolved: a warning is issued because the
  3D embedding then picks an arbitrary configuration.
* Failure to sanitise a structure (e.g. invalid valence) is an error, not a warning.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from rdkit import Chem

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

# Chemistry imports are local to operations, so protein tools do not require RDKit.

from mds.chem.report import PreparationError, Report

log = logging.getLogger("mds.chem.ligand")


def _configure_rdkit_logging():
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")


# Elements for which AutoDock Vina atom types exist.
VINA_ELEMENTS = frozenset({"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "B", "Si", "Se"})
SUPPORTED_SUFFIXES = (".sdf", ".mol", ".mol2", ".smi", ".smiles", ".txt", ".csv")

_IONISABLE = {
    "carboxylic acid (neutral)": "[CX3](=O)[OX2H1]",
    "carboxylate (anionic)": "[CX3](=O)[OX1-]",
    "aliphatic amine (neutral)": "[NX3;H2,H1,H0;!$(NC=[O,S,N]);!$(N-a);!$(N-[SX4]);!$(N#*);!$(N=*);!$(N-N);!$(N-O)]([#6])",
    "amidine/guanidine (neutral)": "[NX3;!$(NC=O)]C(=[NX2;!$(N-a)])",
    "ammonium/amidinium (cationic)": "[N+;!$([N+][O-]);!$([N+]=O)]",
    "phosphate/sulfonate (acidic)": "[P,S](=O)[OX2H1]",
}


@dataclass
class LigandEntry:
    """A ligand as read from a file, before preparation."""
    name: str
    mol: Chem.Mol | None
    error: str = ""
    source: str = ""


@dataclass
class PreparedLigand:
    name: str
    pdbqt: str                          # PDBQT text for the primary conformer
    mol: Chem.Mol                       # hydrogenated 3D molecule
    smiles: str
    properties: dict[str, float | int | str]
    report: Report
    conformer_pdbqts: list[str] = field(default_factory=list)  # one per conformer docked
    n_rotatable: int = 0

    def write(self, directory: Path) -> Path:
        directory.mkdir(parents=True, exist_ok=True)
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in self.name) or "ligand"
        path = directory / f"{safe}.pdbqt"
        path.write_text(self.pdbqt)
        return path


# ---------------------------------------------------------------- reading
def read_ligands(path: Path) -> list[LigandEntry]:
    """Read every molecule of an SDF/MOL/MOL2/SMILES/CSV file."""
    from rdkit import Chem
    _configure_rdkit_logging()
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix not in SUPPORTED_SUFFIXES:
        raise PreparationError(f"Unsupported ligand format '{suffix}'. "
                               f"Supported: {', '.join(SUPPORTED_SUFFIXES)}")
    entries: list[LigandEntry] = []
    if suffix == ".sdf":
        for i, mol in enumerate(Chem.SDMolSupplier(str(path), removeHs=False)):
            name = mol.GetProp("_Name") if mol is not None and mol.HasProp("_Name") else ""
            entries.append(LigandEntry(name or f"{path.stem}_{i + 1}", mol,
                                       "" if mol is not None else "RDKit could not parse record",
                                       path.name))
    elif suffix == ".mol":
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
        entries.append(LigandEntry(path.stem, mol, "" if mol else "RDKit could not parse file",
                                   path.name))
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(path), removeHs=False)
        entries.append(LigandEntry(path.stem, mol, "" if mol else "RDKit could not parse file",
                                   path.name))
    elif suffix == ".csv":
        entries = _read_csv(path)
    else:  # .smi / .smiles / .txt : "SMILES [name]"
        for i, line in enumerate(path.read_text().splitlines()):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            entries.append(_from_smiles(parts[0], parts[1] if len(parts) > 1 else f"{path.stem}_{i + 1}",
                                        path.name))
    if not entries:
        raise PreparationError(f"No molecules found in '{path.name}'")
    return entries


def _from_smiles(smiles: str, name: str, source: str = "") -> LigandEntry:
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles) if smiles.strip() else None
    return LigandEntry(name, mol, "" if mol else f"Invalid SMILES (valence or syntax error): {smiles}",
                       source)


def _read_csv(path: Path) -> list[LigandEntry]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower().strip(): c for c in (reader.fieldnames or [])}
        smi_col = next((cols[k] for k in cols if k in ("smiles", "canonical_smiles", "smile")), None)
        if smi_col is None:
            raise PreparationError("CSV needs a 'smiles' column (header row required).")
        name_col = next((cols[k] for k in cols if k in ("name", "id", "title", "compound", "molecule")),
                        None)
        return [_from_smiles((row.get(smi_col) or "").strip(),
                             ((row.get(name_col) or "").strip() if name_col else "") or f"{path.stem}_{i + 1}",
                             path.name)
                for i, row in enumerate(reader)]


def entry_from_smiles(smiles: str, name: str = "ligand") -> LigandEntry:
    _configure_rdkit_logging()
    return _from_smiles(smiles.strip(), name, "SMILES input")


# ------------------------------------------------------------- properties
def compute_properties(mol: Chem.Mol) -> dict[str, float | int | str]:
    """Descriptors, Lipinski Rule-of-Five and Veber checks (on the heavy-atom graph)."""
    from rdkit import Chem
    from rdkit.Chem import Crippen
    from rdkit.Chem import Descriptors
    from rdkit.Chem import Lipinski
    from rdkit.Chem import rdMolDescriptors
    m = Chem.RemoveHs(mol)
    mw = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    hbd, hba = Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    rot = rdMolDescriptors.CalcNumRotatableBonds(m)
    ro5 = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return {
        "MW": round(mw, 2), "LogP": round(logp, 2), "TPSA": round(tpsa, 1),
        "HBD": hbd, "HBA": hba, "RotB": rot, "HeavyAtoms": m.GetNumHeavyAtoms(),
        "Rings": rdMolDescriptors.CalcNumRings(m), "FormalCharge": Chem.GetFormalCharge(m),
        "Ro5_violations": ro5, "Lipinski_pass": "yes" if ro5 <= 1 else "no",
        "Veber_pass": "yes" if (rot <= 10 and tpsa <= 140) else "no",
    }


# ------------------------------------------------------------ preparation
def _ionisation_notes(mol: Chem.Mol, report: Report, *, input_kept: bool = True) -> None:
    from rdkit import Chem
    found = []
    for label, smarts in _IONISABLE.items():
        patt = Chem.MolFromSmarts(smarts)
        n = len({m[0] for m in mol.GetSubstructMatches(patt)}) if patt is not None else 0
        if n:
            found.append(f"{n}x {label}")
    if found:
        state = "The input protonation state was kept. " if input_kept else "After heuristic ionisation. "
        report.warn("Ionisable groups detected: " + "; ".join(found) +
                    ". " + state + "Confirm the state is appropriate for the target pH.")


def _apply_simple_ionisation(mol: Chem.Mol, report: Report) -> Chem.Mol:
    from rdkit import Chem
    # Edit atoms in place on a copy: reaction products can reorder atoms and
    # mishandle explicit hydrogen counts. Neutral-only matches make this idempotent.
    rules = [
        ("carboxylic acid -> carboxylate", "[CX3](=O)[OX2H1+0]", 2, -1),
        ("aliphatic amine -> ammonium",
         "[NX3+0;H2,H1,H0;!$(NC=[O,S,N]);!$(N-a);!$(N-[SX4]);!$(N#*);!$(N=*);!$(N-N);!$(N-O);!$(N-C=C)]([#6])", 0, 1),
        ("amidine/guanidine -> amidinium",
         "[NX3+0;!$(NC=O)][CX3]=[NX2+0;!$(N-a);H1,H0]", 2, 1),
    ]
    out = Chem.RemoveHs(mol)
    for label, smarts, target, charge in rules:
        indices = sorted({match[target] for match in out.GetSubstructMatches(Chem.MolFromSmarts(smarts))})
        for index in indices:
            atom = out.GetAtomWithIdx(index)
            hydrogens = atom.GetTotalNumHs()
            atom.SetFormalCharge(charge)
            atom.SetNumExplicitHs(hydrogens + charge)
            atom.SetNoImplicit(True)
        if indices:
            try:
                Chem.SanitizeMol(out)
            except Exception as exc:
                report.error(f"Simple ionisation failed ({label}): {exc}")
                raise PreparationError("Simple ionisation produced an invalid structure", report) from exc
            report.info(f"Simple ionisation rule applied ({label}) x{len(indices)}.")
    report.warn("Heuristic ionisation uses fixed textbook rules, not predicted pKa values; "
                "verify charged states for multi-ionisable or unusual molecules.")
    return out


def prepare_ligand(
    entry: LigandEntry,
    *,
    seed: int = 42,
    n_conformers: int = 1,
    ionize: str = "none",
    optimize: bool = True,
    keep_input_3d: bool = False,
) -> PreparedLigand:
    """Validate and prepare a ligand. Keeping input 3D skips optimisation.

    If no input 3D conformer exists, generate and optionally optimise conformers.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDistGeom

    _configure_rdkit_logging()
    rep = Report()
    if ionize not in {"none", "simple"}:
        rep.error(f"Unknown ionisation mode: {ionize}")
        raise PreparationError(rep.errors[-1], rep)
    if isinstance(n_conformers, bool) or not isinstance(n_conformers, int) or n_conformers < 1:
        rep.error("Conformer count must be a positive integer.")
        raise PreparationError(rep.errors[-1], rep)
    if entry.mol is None:
        rep.error(entry.error or "Molecule could not be read.")
        raise PreparationError(f"{entry.name}: {entry.error}", rep)

    mol = Chem.Mol(entry.mol)
    if mol.GetNumHeavyAtoms() == 0:
        rep.error("Molecule contains no heavy atoms.")
        raise PreparationError(f"{entry.name}: empty ligand", rep)
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # noqa: BLE001
        rep.error(f"Sanitisation failed (invalid valence or aromaticity): {exc}")
        raise PreparationError(f"{entry.name}: invalid structure ({exc})", rep) from exc

    # keep the largest fragment (counter-ions / solvent removal), reported.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        biggest = max(frags, key=lambda f: f.GetNumHeavyAtoms())
        rep.warn(f"{len(frags)} fragments found (salt/solvent?); kept the largest "
                 f"({biggest.GetNumHeavyAtoms()} heavy atoms).")
        mol = Chem.Mol(biggest)

    bad = sorted({a.GetSymbol() for a in mol.GetAtoms()} - VINA_ELEMENTS)
    if bad:
        msg = f"Elements without AutoDock Vina atom types: {', '.join(bad)}"
        rep.error(msg)
        raise PreparationError(f"{entry.name}: {msg}", rep)

    # Stereochemistry
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    undefined = [i for i, tag in centres if tag == "?"]
    if undefined:
        rep.warn(f"{len(undefined)} undefined stereocentre(s) (atoms {undefined}). The 3D embedding "
                 "picks an arbitrary configuration; specify stereochemistry if it matters.")
    elif centres:
        rep.info(f"{len(centres)} defined stereocentre(s) preserved.")

    if ionize == "simple":
        mol = _apply_simple_ionisation(mol, rep)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    _ionisation_notes(mol, rep, input_kept=ionize == "none")

    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    props = compute_properties(mol)

    # 3D coordinates
    has_3d = mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()
    mol_h = Chem.AddHs(mol, addCoords=has_3d and keep_input_3d)
    if not (keep_input_3d and has_3d):
        mol_h.RemoveAllConformers()
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = seed
        params.useSmallRingTorsions = True
        n_req = n_conformers
        cids = list(rdDistGeom.EmbedMultipleConfs(mol_h, numConfs=n_req, params=params))
        if not cids:
            params.useRandomCoords = True
            cids = list(rdDistGeom.EmbedMultipleConfs(mol_h, numConfs=n_req, params=params))
        if not cids:
            rep.error("3D embedding failed (strained or unusual structure).")
            raise PreparationError(f"{entry.name}: 3D embedding failed", rep)
        rep.info(f"Generated {len(cids)} 3D conformer(s) with ETKDGv3 (seed {seed}).")
    else:
        rep.info("Input 3D coordinates were kept.")
    if optimize and not (keep_input_3d and has_3d):
        _optimise(mol_h, rep)
    elif optimize:
        rep.info("Optimisation skipped to preserve input 3D coordinates.")

    # Meeko
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        prep = MoleculePreparation()
    except Exception as exc:
        rep.error(f"Meeko could not initialise: {exc}")
        raise PreparationError(f"{entry.name}: Meeko unavailable", rep) from exc
    pdbqts: list[str] = []
    n_rot = 0
    for conf in mol_h.GetConformers():
        single = _single_conformer(mol_h, conf.GetId())
        try:
            setups = prep.prepare(single)
            if not setups:
                raise ValueError("Meeko returned no setup")
            text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        except Exception as exc:
            rep.error(f"Meeko preparation failed: {exc}")
            raise PreparationError(f"{entry.name}: Meeko preparation failed", rep) from exc
        if not ok:
            rep.error(f"Meeko PDBQT writing failed: {err}")
            raise PreparationError(f"{entry.name}: {err}", rep)
        pdbqts.append(text)
        n_rot = sum(1 for ln in text.splitlines() if ln.startswith("BRANCH"))
    props["Torsions(Vina)"] = n_rot
    if n_rot > 20:
        rep.warn(f"{n_rot} torsions: docking very flexible ligands is unreliable; "
                 "increase exhaustiveness substantially.")
    rep.info(f"PDBQT written ({n_rot} active torsions). Ligand ready for docking.")

    return PreparedLigand(entry.name, pdbqts[0], mol_h, smiles, props, rep, pdbqts, n_rot)


def _single_conformer(mol: Chem.Mol, conf_id: int) -> Chem.Mol:
    from rdkit import Chem
    single = Chem.Mol(mol)
    single.RemoveAllConformers()
    single.AddConformer(Chem.Conformer(mol.GetConformer(conf_id)), assignId=True)
    return single


def _optimise(mol_h: Chem.Mol, rep: Report) -> None:
    from rdkit.Chem import rdForceFieldHelpers
    if rdForceFieldHelpers.MMFFHasAllMoleculeParams(mol_h):
        method = "MMFF94"
        optimise = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs
    elif rdForceFieldHelpers.UFFHasAllMoleculeParams(mol_h):
        method = "UFF"
        optimise = rdForceFieldHelpers.UFFOptimizeMoleculeConfs
        rep.warn("MMFF94 parameters unavailable for this molecule; UFF was used instead.")
    else:
        rep.error("Neither MMFF94 nor UFF has all parameters for this molecule.")
        raise PreparationError("Cannot optimise ligand: missing force-field parameters", rep)
    try:
        res = optimise(mol_h, maxIters=2000, numThreads=0)
    except Exception as exc:
        rep.error(f"{method} optimisation failed: {exc}")
        raise PreparationError("Ligand optimisation failed", rep) from exc
    if len(res) != mol_h.GetNumConformers() or any(c < 0 for c, _ in res):
        rep.error(f"{method} could not optimise all conformers.")
        raise PreparationError("Ligand optimisation failed", rep)
    converged = sum(1 for c, _ in res if c == 0)
    rep.info(f"{method} optimisation: {converged}/{len(res)} conformer(s) converged.")
    if converged < len(res):
        rep.warn(f"Some conformers did not converge within 2000 {method} iterations.")
