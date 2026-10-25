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


# Non-default states AMBER naming assigns, with the standard residue each comes
# from and what it means. "alters_charge" marks the ones that change the formal
# charge rather than only a tautomer or bonding choice, because those are what
# shift binding-site electrostatics and most need review.
PROTONATION_STATES = {
    "HID": ("HIS", "neutral histidine, proton on ND1 (delta tautomer)", False),
    "HIE": ("HIS", "neutral histidine, proton on NE2 (epsilon tautomer)", False),
    "HIP": ("HIS", "doubly protonated histidine, carries +1", True),
    "ASH": ("ASP", "neutral (protonated) aspartate; normally -1", True),
    "GLH": ("GLU", "neutral (protonated) glutamate; normally -1", True),
    "LYN": ("LYS", "neutral (deprotonated) lysine; normally +1", True),
    "CYX": ("CYS", "cysteine in a disulfide bond, no thiol hydrogen", False),
}


@dataclass
class ProtonationState:
    chain: str
    seqid: int
    residue: str      # standard name, e.g. HIS
    state: str        # assigned state, e.g. HIP
    meaning: str
    alters_charge: bool


def collect_protonation(st: gemmi.Structure) -> list["ProtonationState"]:
    """Per-residue non-default protonation/bonding states from AMBER naming."""
    found: list[ProtonationState] = []
    for chain in st[0]:
        for res in chain:
            # AMBER prefixes terminal residues, e.g. NHID / CHID.
            key = res.name if res.name in PROTONATION_STATES else (
                res.name[1:] if len(res.name) == 4 and res.name[1:] in PROTONATION_STATES else None)
            if key is None:
                continue
            standard, meaning, alters = PROTONATION_STATES[key]
            found.append(ProtonationState(chain.name or "", res.seqid.num, standard, key, meaning, alters))
    return found


@dataclass
class ReceptorResult:
    pdbqt: Path
    prepared_pdb: Path
    report: Report
    hetero: list[HeteroResidue] = field(default_factory=list)
    n_atoms: int = 0
    n_residues: int = 0
    protonation: list[ProtonationState] = field(default_factory=list)


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
                extra = f" (numbering gap of {missing}; missing sequence is unconfirmed)" if missing > 0 else ""
                breaks.append(f"chain {ch.name}: break between {a.name}{a.seqid.num} and "
                              f"{b.name}{b.seqid.num}{extra}")
    return breaks


def _check_options(opts: ReceptorOptions, report: Report) -> None:
    if not math.isfinite(opts.ph) or not 0 <= opts.ph <= 14:
        report.error("Protonation pH must be finite and between 0 and 14.")
    if not opts.drop_waters:
        report.error("Retaining explicit waters is not supported by this receptor pipeline. "
                     "Use drop_waters=True; waters will not be silently discarded.")


def inspect_receptor(source: Path, options: ReceptorOptions | None = None) -> Report:
    """Read-only review of the exact chains/options selected for preparation."""
    opts, rep = options or ReceptorOptions(), Report()
    _check_options(opts, rep)
    st = read_structure(Path(source))
    if len(st) > 1:
        rep.warn(f"{len(st)} models present; preparation uses only the first model.")
    available = {ch.name for ch in st[0]}
    missing = set(opts.chains or []) - available
    if missing:
        rep.error(f"Chain(s) {sorted(missing)} not found; available: {sorted(available)}")
        return rep
    if opts.chains:
        for i in reversed(range(len(st[0]))):
            if st[0][i].name not in opts.chains:
                del st[0][i]
    polymer = [(ch.name, r) for ch in st[0] for r in ch if _is_polymer_residue(r)]
    if not polymer:
        rep.error("No protein residues in the selected chains.")
    rep.info("Selected chains: " + ", ".join(ch.name or "(blank)" for ch in st[0]) +
             f"; {len(polymer)} protein residue(s).")
    for chain, residue in polymer:
        label = f"{residue.name} {chain}:{residue.seqid}"
        absent = [name for name in ("N", "CA", "C", "O") if residue.find_atom(name, "*") is None]
        if absent:
            rep.warn(f"Incomplete backbone at {label}: missing {', '.join(absent)}. "
                     "Preparation may rebuild atoms or remove an unreconstructable residue.")
        alternatives = sorted({a.altloc for a in residue if a.altloc != "\0"})
        if alternatives:
            rep.warn(f"Alternate conformations at {label}: {', '.join(alternatives)}. "
                     "Preparation keeps the first conformer, not an occupancy-weighted choice.")
    for line in find_chain_breaks(st):
        rep.warn(f"Chain break: {line}. Missing segments are not rebuilt.")
    hetero = classify_hetero(st)
    waters = sum(h.kind == "water" for h in hetero)
    rep.info(f"Waters: {waters}; all will be removed for standard Vina preparation.")
    for h in hetero:
        if h.kind == "metal":
            rep.info(f"Metal {h.label}: " + ("retained with a Meeko template." if opts.keep_metals
                                          else "removed by the selected Keep metal ions setting."))
        elif h.kind == "ion":
            rep.warn(f"Ion {h.label}: removed; no supported metal template identified.")
        elif h.kind == "ligand":
            rep.warn(f"Hetero group {h.label}: extracted and excluded from the receptor. "
                     "Review whether it is an essential cofactor before docking.")
    rep.info(f"Protonation: PDB2PQR/PROPKA will predict states at pH {opts.ph:g}. "
             "Review binding-site histidines, acids and metal-coordinating residues after preparation. "
             "This check does not predict pKa or rebuild missing structure.")
    return rep


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
    _check_options(opts, rep)
    if rep.errors:
        raise PreparationError("Invalid receptor preparation options", rep)
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
            elif h.kind == "metal":
                rep.info(f"Metal ion {h.label} removed (Keep metal ions is disabled).")
                drop.append(i)
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

    # -- drop residues PDB2PQR cannot handle (incomplete backbone at chain gaps/ends).
    # Without this, a single stub residue causes "Too few atoms present to reconstruct..."
    # and the whole workflow fails even though the stub is usually far from the binding site.
    _drop_backbone_stubs(st, hetero, rep)
    n_res_in = sum(len(c) for c in st[0])
    if n_res_in == 0:
        raise PreparationError("No protein residues remain after cleaning.", rep)

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
    # Record which residue got which state, not just how many: a reviewer needs
    # to know that His57 specifically is +1, and where it sits.
    protonation = collect_protonation(after_st)
    states = Counter(p.state for p in protonation)
    if states:
        rep.info("Predicted non-default protonation/bonding states at pH "
                 f"{opts.ph}: " + ", ".join(f"{n}x{s}" for s, n in sorted(states.items())) +
                 ". Review histidines and acids in the binding site before trusting results.")
        charge_changing = [p for p in protonation if p.alters_charge]
        for p in charge_changing[:25]:
            rep.info(f"{p.residue} {p.chain}:{p.seqid} assigned {p.state}: {p.meaning}.")
        if len(charge_changing) > 25:
            rep.info(f"{len(charge_changing) - 25} further charge-changing assignment(s) "
                     "are listed in the preparation result.")

    # -- re-attach metals and write final PDB
    prepared = out_dir / "receptor_prepared.pdb"
    lines = [l for l in pdb_h.read_text().splitlines() if not l.startswith(("END", "TER"))]
    if len(kept_ions[0]):
        ion_text = _minimal_pdb_text(kept_ions, work / "ions.pdb")
        lines += [l for l in ion_text.splitlines() if l.startswith("HETATM")]
    prepared.write_text("\n".join(lines) + "\nEND\n")

    # -- Meeko: PDB -> PDBQT.
    # First attempt is strict. If Meeko refuses because a handful of residues do not match
    # its templates (common on real PDBs: non-standard C-terminal residues, odd protonation
    # states, modified residues at chain ends), retry with --allow_bad_res so those residues
    # are dropped and the rest of the receptor still gets built. Each dropped residue is
    # named in the report with its distance to the nearest extracted ligand, so the user
    # can judge whether it affects their binding site.
    pdbqt = out_dir / "receptor.pdbqt"
    code, output = _run_meeko(["--read_pdb", str(prepared), "-p", str(pdbqt)])
    if code != 0 or not pdbqt.exists():
        bad_res = _parse_meeko_bad_residues(output)
        reason = ("template mismatch" if bad_res else
                  _meeko_crash_residue(output))
        # Even without a clear residue list we retry once with --allow_bad_res: Meeko sometimes
        # raises an exception on a single oddball residue instead of listing it as a mismatch.
        if bad_res or reason:
            rep.info("Meeko refused the receptor on the strict first pass; retrying with "
                     "--allow_bad_res to drop problematic residues and keep the rest.")
            for res_id, dist in _annotate_distances(bad_res, prepared, hetero):
                near = " (NEAR the extracted ligand — may affect your binding site!)" if dist is not None and dist < 10 else ""
                dist_text = f", {dist:.1f} A from nearest extracted ligand" if dist is not None else ""
                rep.warn(f"Residue {res_id} dropped by Meeko (template mismatch{dist_text}){near}")
            if reason and not bad_res:
                rep.warn(f"Meeko exception on first pass: {reason}. Retrying with --allow_bad_res.")
            code, output = _run_meeko(["--read_pdb", str(prepared), "-p", str(pdbqt), "--allow_bad_res"])
            if code == 0 and pdbqt.exists():
                # The retry output lists exactly which residues it dropped; surface them.
                dropped = _parse_meeko_bad_residues(output)
                for res_id, dist in _annotate_distances(dropped, prepared, hetero):
                    if res_id in bad_res:
                        continue  # already reported above
                    near = " (NEAR the extracted ligand — may affect your binding site!)" if dist is not None and dist < 10 else ""
                    dist_text = f", {dist:.1f} A from nearest extracted ligand" if dist is not None else ""
                    rep.warn(f"Residue {res_id} dropped by Meeko --allow_bad_res{dist_text}{near}")
        # If Meeko still crashes on a specific residue (physically broken geometry — alt-conf
        # debris, bad bond lengths — that PDB2PQR trusted but Meeko rejects), delete that
        # residue by hand and try one more time.
        if (code != 0 or not pdbqt.exists()) and (crash_res := _meeko_crashed_on(output)):
            dist = next((d for rid, d in _annotate_distances([crash_res], prepared, hetero)), None)
            near = " (NEAR the extracted ligand — may affect your binding site!)" if dist is not None and dist < 10 else ""
            dist_text = f", {dist:.1f} A from nearest extracted ligand" if dist is not None else ""
            rep.warn(f"Residue {crash_res} removed: physically broken geometry that Meeko cannot "
                     f"process (likely a disordered side chain from the X-ray data){dist_text}{near}.")
            fixed = work / "receptor_prepared_fixed.pdb"
            _write_pdb_without_residue(prepared, fixed, crash_res)
            code, output = _run_meeko(["--read_pdb", str(fixed), "-p", str(pdbqt), "--allow_bad_res"])
        if code != 0 or not pdbqt.exists():
            tail = "\n".join(output.strip().splitlines()[-8:])
            rep.error(f"Meeko could not build the receptor PDBQT:\n{tail}")
            raise PreparationError("Meeko receptor preparation failed", rep)
    for line in output.splitlines():
        if "warning" in line.lower():
            rep.warn(f"Meeko: {line.strip()}")

    n_atoms = validate_pdbqt(pdbqt, rep)
    rep.info(f"Receptor PDBQT written: {pdbqt.name} ({n_atoms} atoms).")
    return ReceptorResult(pdbqt, prepared, rep, hetero, n_atoms, n_res_in, protonation)


def _minimal_pdb_text(st: gemmi.Structure, scratch: Path) -> str:
    st.write_minimal_pdb(str(scratch))
    return scratch.read_text()


# Backbone atom names required by PDB2PQR to keep (or reconstruct a cap for) a residue.
_BACKBONE = ("N", "CA", "C", "O")


def _drop_backbone_stubs(st: gemmi.Structure, hetero: list[HeteroResidue], report: Report) -> None:
    """Drop polymer residues whose backbone is too incomplete for PDB2PQR to handle.

    These typically sit at the ragged edges of crystallographic chain gaps or the chain's own
    terminus. Keeping them makes PDB2PQR fail with "Too few atoms present to reconstruct";
    dropping just them lets the rest of the receptor through.
    """
    hetero_centroids = _hetero_centroids(hetero)
    for ch in list(st[0]):
        to_remove: list[int] = []
        for i, res in enumerate(ch):
            if not _is_polymer_residue(res):
                continue
            missing = [n for n in _BACKBONE if res.find_atom(n, "*") is None]
            # A C-terminal residue with just O missing is a routine case PDB2PQR handles via
            # OXT; only treat it as a stub if N or CA or C are missing.
            if not any(a in missing for a in ("N", "CA", "C")):
                continue
            to_remove.append(i)
            dist = _nearest_hetero_distance(res, hetero_centroids)
            near = " (NEAR an extracted ligand — may affect your binding site!)" \
                if dist is not None and dist < 10 else ""
            dist_text = f", {dist:.1f} A from nearest extracted ligand" if dist is not None else ""
            report.warn(f"Residue {res.name} {ch.name}:{res.seqid.num} dropped: incomplete "
                        f"backbone (missing {', '.join(missing)}); PDB2PQR cannot repair it"
                        f"{dist_text}{near}.")
        for i in reversed(to_remove):
            del ch[i]
    # chains left empty after removing stubs
    for name in [c.name for c in st[0] if len(c) == 0]:
        st[0].remove_chain(name)


def _parse_meeko_bad_residues(output: str) -> list[str]:
    """Extract residue IDs Meeko named as template-match failures, e.g. "A:639"."""
    import re
    out: list[str] = []
    for pattern in (
        r"Template matching failed for:\s*\[([^\]]+)\]",   # strict run
        r"Template matching failed for:\s*\['([^']+)'\]",  # retry output format
    ):
        for m in re.finditer(pattern, output):
            out.extend(token.strip().strip("'\"") for token in m.group(1).split(",")
                       if token.strip())
    # also Meeko --allow_bad_res summary lines: "'A:639' Ignored due to allow_bad_res"
    for m in re.finditer(r"\['([^']+)'\]\s+Ignored due to allow_bad_res", output):
        out.append(m.group(1))
    # de-duplicate while preserving order
    seen: set[str] = set()
    deduped = []
    for r in out:
        if r not in seen:
            deduped.append(r)
            seen.add(r)
    return deduped


def _meeko_crash_residue(output: str) -> str:
    """Return a short description if Meeko raised a Python exception on a specific residue."""
    import re
    m = re.search(r"RuntimeError:\s*unable to build rdkit mol for residue \S+ corresponding to "
                  r"key (\S+)", output)
    if m:
        return f"could not build molecule for {m.group(1)}"
    m = re.search(r"RuntimeError:\s*(.+)", output)
    return m.group(1).strip()[:160] if m else ""


def _meeko_crashed_on(output: str) -> str | None:
    """Return the chain:number key of the residue that triggered a Meeko RuntimeError, if any."""
    import re
    m = re.search(r"unable to build rdkit mol for residue \S+ corresponding to key (\S+)", output)
    return m.group(1) if m else None


def _write_pdb_without_residue(src: Path, dst: Path, res_id: str) -> None:
    """Copy src -> dst, dropping every ATOM/HETATM line of the chain:number residue in res_id."""
    try:
        chain, num = res_id.split(":", 1)
        num_i = int(num)
    except ValueError:
        dst.write_text(src.read_text())
        return
    kept: list[str] = []
    for line in src.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                line_chain = line[21].strip()
                line_num = int(line[22:26])
            except (IndexError, ValueError):
                kept.append(line)
                continue
            if line_chain == chain and line_num == num_i:
                continue
        kept.append(line)
    dst.write_text("\n".join(kept) + "\n")


def _annotate_distances(residue_ids: list[str], prepared_pdb: Path,
                        hetero: list[HeteroResidue]) -> list[tuple[str, float | None]]:
    """For each 'A:639'-style residue ID, return (res_id, distance_to_nearest_hetero_A)."""
    centroids = _hetero_centroids(hetero)
    if not centroids:
        return [(rid, None) for rid in residue_ids]
    try:
        st = gemmi.read_structure(str(prepared_pdb))
    except Exception:  # noqa: BLE001
        return [(rid, None) for rid in residue_ids]
    out = []
    for rid in residue_ids:
        try:
            chain, num = rid.split(":", 1)
            num = int(num)
        except ValueError:
            out.append((rid, None))
            continue
        res = next((r for ch in st[0] if ch.name == chain for r in ch if r.seqid.num == num), None)
        out.append((rid, _nearest_hetero_distance(res, centroids) if res else None))
    return out


def _hetero_centroids(hetero: list[HeteroResidue]) -> dict[str, tuple[float, float, float]]:
    """Compute centroids of each extracted ligand/metal from the files on disk."""
    out: dict[str, tuple[float, float, float]] = {}
    for h in hetero:
        if h.extracted_to is None or not h.extracted_to.exists():
            continue
        xs, ys, zs = [], [], []
        for line in h.extracted_to.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    continue
        if xs:
            out[h.label] = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    return out


def _nearest_hetero_distance(res: gemmi.Residue | None,
                             centroids: dict[str, tuple[float, float, float]]) -> float | None:
    if res is None or not centroids:
        return None
    atoms = [(a.pos.x, a.pos.y, a.pos.z) for a in res]
    if not atoms:
        return None
    best = math.inf
    for cx, cy, cz in centroids.values():
        for ax, ay, az in atoms:
            d = math.sqrt((ax - cx) ** 2 + (ay - cy) ** 2 + (az - cz) ** 2)
            if d < best:
                best = d
    return best if math.isfinite(best) else None


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
