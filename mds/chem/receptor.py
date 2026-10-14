"""Receptor (protein) preparation for AutoDock Vina.

Pipeline
--------
1. **Gemmi** reads PDB/mmCIF and applies explicit, *reported* decisions: first
   model only, chain selection, alternate conformations, waters, ions,
   co-crystallised ligands (extracted to separate files, never silently lost).
2. **PDB2PQR + PROPKA** rebuild missing heavy atoms, predict side-chain pKa at the
   requested pH, choose histidine tautomers / protonation variants (HID/HIE/HIP,
   ASH, GLH, LYN, CYX), add hydrogens and optimise the hydrogen-bond network.
3. Retained metal ions are re-attached to the prepared structure.
4. **Meeko** converts the hydrogenated PDB into a receptor PDBQT (AutoDock atom
   types + Gasteiger charges from residue templates).

Every operation and warning is written to a :class:`Report`.  Limitations are
stated, not hidden: cofactors such as heme/NAD are removed (they need their own
atom parameters, which this version does not provide), and protonation states
are *predictions* that should be reviewed for important binding-site residues.
"""

from __future__ import annotations

import contextlib
import io
import logging
import math
import sys
import threading
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import gemmi

from mds.chem.report import PreparationError, Report

log = logging.getLogger("mds.chem.receptor")

WATERS = frozenset({"HOH", "DOD", "WAT", "H2O"})
# Metal ions that Meeko has atom-type templates for (alias -> template name).
SUPPORTED_METALS = {"ZN": "ZN", "MG": "MG", "CA": "CA", "MN": "MN", "FE": "FE",
                    "FE2": "FE", "CU": "CU", "CU1": "CU"}

_meeko_lock = threading.Lock()  # Meeko's CLI reads sys.argv, so serialise calls.


@dataclass
class ReceptorOptions:
    chains: list[str] | None = None        # None = all polymer chains
    keep_metals: bool = True               # retain supported catalytic/structural metals
    ph: float = 7.4
    drop_waters: bool = True


@dataclass
class HeteroResidue:
    name: str
    chain: str
    number: int
    atoms: int
    kind: str                              # "ligand", "metal", "ion", "water", "modified"
    extracted_to: Path | None = None

    @property
    def label(self) -> str:
        return f"{self.name} {self.chain}:{self.number}"


@dataclass
class ReceptorResult:
    pdbqt: Path
    prepared_pdb: Path
    report: Report
    hetero: list[HeteroResidue] = field(default_factory=list)
    n_atoms: int = 0
    n_residues: int = 0


# --- inspection
def read_structure(path: Path) -> gemmi.Structure:
    try:
        st = gemmi.read_structure(str(path))
    except Exception as exc:  # gemmi raises RuntimeError / ValueError subclasses
        raise PreparationError(f"Cannot read structure '{path.name}': {exc}") from exc
    st.setup_entities()
    if len(st) == 0 or sum(len(ch) for ch in st[0]) == 0:
        raise PreparationError(f"'{path.name}' contains no atoms")
    return st


def _is_polymer_residue(res: gemmi.Residue) -> bool:
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and info.is_amino_acid())


def classify_hetero(st: gemmi.Structure) -> list[HeteroResidue]:
    """List every non-standard residue in the first model without modifying it."""
    out: list[HeteroResidue] = []
    for ch in st[0]:
        for res in ch:
            if res.name in WATERS:
                kind = "water"
            elif res.het_flag != "H" and _is_polymer_residue(res):
                continue
            elif res.name == "MSE":
                kind = "modified"
            elif res.name in SUPPORTED_METALS and len(res) == 1:
                kind = "metal"
            elif len(res) == 1:
                kind = "ion"
            elif _is_polymer_residue(res):
                continue
            else:
                kind = "ligand"
            out.append(HeteroResidue(res.name, ch.name, res.seqid.num, len(res), kind))
    return out


def find_chain_breaks(st: gemmi.Structure) -> list[str]:
    """Describe gaps where consecutive residues are not covalently connected."""
    breaks = []
    for ch in st[0]:
        poly = [r for r in ch if _is_polymer_residue(r)]
        for a, b in zip(poly, poly[1:]):
            c, n = a.find_atom("C", "*"), b.find_atom("N", "*")
            if c is None or n is None:
                continue
            if c.pos.dist(n.pos) > 2.0:
                missing = b.seqid.num - a.seqid.num - 1
                extra = f" ({missing} residue(s) absent from the model)" if missing > 0 else ""
                breaks.append(f"chain {ch.name}: break between {a.name}{a.seqid.num} and "
                              f"{b.name}{b.seqid.num}{extra}")
    return breaks


# ------------------------------------------------------------------ helpers
def _single_residue_structure(res: gemmi.Residue, chain_name: str) -> gemmi.Structure:
    ns = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain(chain_name)
    chain.add_residue(res)
    model.add_chain(chain)
    ns.add_model(model)
    return ns


def _heavy_atom_counts(st: gemmi.Structure) -> Counter:
    counts: Counter = Counter()
    for ch in st[0]:
        for res in ch:
            counts[(ch.name, res.seqid.num, res.seqid.icode.strip())] = sum(
                1 for a in res if not a.is_hydrogen())
    return counts


def _run_meeko(argv: list[str]) -> tuple[int, str]:
    """Call Meeko's ``mk_prepare_receptor`` in-process (works inside a frozen .exe)."""
    from meeko.cli import mk_prepare_receptor

    buf = io.StringIO()
    with _meeko_lock:
        old_argv = sys.argv
        sys.argv = ["mk_prepare_receptor", *argv]
        code = 0
        try:
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                mk_prepare_receptor.main()
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
        except Exception as exc:  # noqa: BLE001 — surface Meeko's message to the user
            buf.write(f"\n{type(exc).__name__}: {exc}")
            code = 1
        finally:
            sys.argv = old_argv
    return code, buf.getvalue()


class _CollectWarnings(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


def _run_pdb2pqr(pdb_in: Path, pqr_out: Path, pdb_out: Path, ph: float) -> list[str]:
    from pdb2pqr.main import run_pdb2pqr

    collector = _CollectWarnings()
    root = logging.getLogger("pdb2pqr")
    root.addHandler(collector)
    try:
        run_pdb2pqr([
            "--ff=AMBER", "--ffout=AMBER", f"--with-ph={ph}",
            "--titration-state-method=propka", f"--pdb-output={pdb_out}",
            "--keep-chain", str(pdb_in), str(pqr_out),
        ])
    except SystemExit as exc:
        raise PreparationError(f"PDB2PQR terminated (exit {exc.code})") from exc
    except Exception as exc:  # noqa: BLE001
        raise PreparationError(f"PDB2PQR failed: {exc}") from exc
    finally:
        root.removeHandler(collector)
    return collector.messages


# ------------------------------------------------------------------ main API
def prepare_receptor(
    source: Path,
    out_dir: Path,
    options: ReceptorOptions | None = None,
    report: Report | None = None,
) -> ReceptorResult:
    """Prepare ``source`` (PDB/mmCIF) and write ``receptor.pdbqt`` into ``out_dir``."""
    opts = options or ReceptorOptions()
    rep = report or Report()
    out_dir = Path(out_dir)
    work = out_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    lig_dir = out_dir / "extracted_ligands"

    st = read_structure(Path(source))
    rep.info(f"Read '{Path(source).name}': {len(st)} model(s), "
             f"{len(st[0])} chain(s), {sum(len(c) for c in st[0])} residues.")

    # -- model selection
    if len(st) > 1:
        rep.warn(f"Structure has {len(st)} models (e.g. NMR ensemble); only model 1 is used.")
        while len(st) > 1:
            del st[len(st) - 1]

    # -- alternate conformations
    n_alt = sum(1 for ch in st[0] for r in ch for a in r if a.altloc != "\0")
    if n_alt:
        st.remove_alternative_conformations()
        rep.warn(f"{n_alt} atoms had alternate conformations; the first conformer of each was "
                 "kept (occupancy-weighted choice is not applied).")

    # -- chain selection
    if opts.chains:
        available = [c.name for c in st[0]]
        missing = [c for c in opts.chains if c not in available]
        if missing:
            raise PreparationError(f"Chain(s) {missing} not found; available: {available}")
        for name in available:
            if name not in opts.chains:
                st[0].remove_chain(name)
                rep.info(f"Chain {name} removed (not selected).")

    # -- selenomethionine -> methionine
    for ch in st[0]:
        for res in ch:
            if res.name == "MSE":
                res.name = "MET"
                res.het_flag = "A"
                for a in res:
                    if a.name.strip() == "SE":
                        a.name, a.element = "SD", gemmi.Element("S")
                rep.info(f"Selenomethionine {ch.name}:{res.seqid.num} converted to methionine.")

    # -- hetero residues
    hetero = classify_hetero(st)
    kept_ions = gemmi.Structure()
    kept_model = gemmi.Model("1")
    water_count = 0
    for ch in st[0]:
        drop: list[int] = []
        for i, res in enumerate(ch):
            h = next((x for x in hetero if (x.chain, x.number, x.name) ==
                      (ch.name, res.seqid.num, res.name)), None)
            if h is None:
                if not _is_polymer_residue(res):
                    rep.warn(f"Non-protein residue {res.name} {ch.name}:{res.seqid.num} "
                             "removed (nucleic acids/unknown polymers are not supported).")
                    drop.append(i)
                continue
            if h.kind == "water":
                water_count += 1
                drop.append(i)
            elif h.kind == "metal" and opts.keep_metals:
                template = SUPPORTED_METALS[res.name]
                res.name = template
                res[0].name = template
                tmp_chain = gemmi.Chain(ch.name)
                tmp_chain.add_residue(res)
                kept_model.add_chain(tmp_chain)
                rep.info(f"Metal ion {h.label} retained (parameterised by Meeko).")
                drop.append(i)  # re-attached after PDB2PQR
            elif h.kind == "ligand":
                lig_dir.mkdir(parents=True, exist_ok=True)
                path = lig_dir / f"{h.name}_{h.chain}{h.number}.pdb"
                _single_residue_structure(res, ch.name).write_pdb(str(path))
                h.extracted_to = path
                rep.warn(f"Hetero group {h.label} ({h.atoms} atoms) removed from the receptor "
                         f"and saved to '{path.name}'. If it is a cofactor essential for "
                         "binding (heme, NAD, FAD, etc.), it is NOT in the docking receptor — "
                         "cofactor parameterisation is not supported in this version.")
                drop.append(i)
            else:
                rep.warn(f"Ion {h.label} removed (no AutoDock parameters bundled for it).")
                drop.append(i)
        for i in reversed(drop):
            del ch[i]
    if water_count:
        rep.info(f"{water_count} water molecule(s) removed "
                 "(explicit waters are not used in standard Vina docking).")
    kept_ions.add_model(kept_model)

    # -- drop chains left empty
    for name in [c.name for c in st[0] if len(c) == 0]:
        st[0].remove_chain(name)
    n_res_in = sum(len(c) for c in st[0])
    if n_res_in == 0:
        raise PreparationError("No protein residues remain after cleaning.", rep)

    for line in find_chain_breaks(st):
        rep.warn(f"Chain break — {line}. Missing segments are NOT rebuilt; residues next to "
                 "the gap get artificial termini.")

    before = _heavy_atom_counts(st)
    cleaned = work / "cleaned.pdb"
    st.write_minimal_pdb(str(cleaned))

    # -- hydrogens, pKa, repair
    rep.info(f"Running PDB2PQR/PROPKA at pH {opts.ph} (AMBER naming).")
    pqr, pdb_h = work / "pdb2pqr.pqr", work / "pdb2pqr.pdb"
    for msg in dict.fromkeys(_run_pdb2pqr(cleaned, pqr, pdb_h, opts.ph)):  # de-duplicated
        if msg.startswith("Missing atom OXT"):
            continue  # terminal carboxylate oxygen: rebuilt routinely, not a data problem
        rep.warn(f"PDB2PQR: {msg}")

    after_st = gemmi.read_structure(str(pdb_h))
    after = _heavy_atom_counts(after_st)
    rebuilt = {k: after[k] - before[k] for k in before if k in after and after[k] > before[k]}
    for (chain, num, _), n in sorted(rebuilt.items()):
        res_after = next((r for c in after_st[0] if c.name == chain for r in c
                          if r.seqid.num == num), None)
        name = res_after.name if res_after is not None else "?"
        if n == 1 and res_after is not None and res_after.find_atom("OXT", "*") is not None                 and before[(chain, num, "")] == len([a for a in res_after if not a.is_hydrogen()]) - 1:
            rep.info(f"Residue {name} {chain}:{num}: C-terminal OXT oxygen added.")
            continue
        rep.warn(f"Residue {name} {chain}:{num}: {n} missing heavy atom(s) rebuilt by PDB2PQR "
                 "(side-chain placement is approximate).")
    states = Counter(r.name for c in after_st[0] for r in c
                     if r.name in {"HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX"} or
                     (len(r.name) == 4 and r.name[1:] in {"HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX"}))
    if states:
        rep.info("Predicted non-default protonation/bonding states at pH "
                 f"{opts.ph}: " + ", ".join(f"{n}x{s}" for s, n in sorted(states.items())) +
                 ". Review histidines and acids in the binding site before trusting results.")

    # -- re-attach metals and write final PDB
    prepared = out_dir / "receptor_prepared.pdb"
    lines = [l for l in pdb_h.read_text().splitlines() if not l.startswith(("END", "TER"))]
    if len(kept_ions[0]):
        ion_text = _minimal_pdb_text(kept_ions, work / "ions.pdb")
        lines += [l for l in ion_text.splitlines() if l.startswith("HETATM")]
    prepared.write_text("\n".join(lines) + "\nEND\n")

    # -- Meeko: PDB -> PDBQT
    pdbqt = out_dir / "receptor.pdbqt"
    code, output = _run_meeko(["--read_pdb", str(prepared), "-p", str(pdbqt)])
    if code != 0 or not pdbqt.exists():
        tail = "\n".join(output.strip().splitlines()[-8:])
        rep.error(f"Meeko could not build the receptor PDBQT:\n{tail}")
        raise PreparationError("Meeko receptor preparation failed", rep)
    for line in output.splitlines():
        if "warning" in line.lower():
            rep.warn(f"Meeko: {line.strip()}")

    n_atoms = validate_pdbqt(pdbqt, rep)
    rep.info(f"Receptor PDBQT written: {pdbqt.name} ({n_atoms} atoms).")
    return ReceptorResult(pdbqt, prepared, rep, hetero, n_atoms, n_res_in)


def _minimal_pdb_text(st: gemmi.Structure, scratch: Path) -> str:
    st.write_minimal_pdb(str(scratch))
    return scratch.read_text()


def validate_pdbqt(path: Path, report: Report) -> int:
    """Sanity-check a PDBQT file; raises :class:`PreparationError` if unusable."""
    n = 0
    types: Counter = Counter()
    total_charge = 0.0
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            n += 1
            try:
                x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                q = float(line[70:76])
            except ValueError as exc:
                raise PreparationError(f"Malformed PDBQT atom record: {line!r}") from exc
            if not all(map(math.isfinite, (x, y, z, q))):
                raise PreparationError(f"Non-finite value in PDBQT record: {line!r}")
            total_charge += q
            types[line[77:79].strip()] += 1
    if n == 0:
        raise PreparationError(f"'{Path(path).name}' contains no atoms")
    report.info("Receptor atom types: " + ", ".join(f"{t}:{c}" for t, c in sorted(types.items())) +
                f"; sum of partial charges = {total_charge:+.2f} e.")
    return n
