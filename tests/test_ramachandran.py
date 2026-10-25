"""Tests for Ramachandran backbone-torsion screening and plotting."""
from pathlib import Path

from mds.chem.health import check_protein
from mds.chem.ramachandran import classify, compute_ramachandran, render_svg, summarize

DATA = Path(__file__).parent / 'data' / '3PTB.pdb'


def test_classify_canonical_basins():
    assert classify(-63, -43, 'general') == 'favored'   # right-handed helix
    assert classify(-120, 130, 'general') == 'favored'  # beta sheet
    assert classify(60, 45, 'general') == 'favored'     # left-handed helix
    assert classify(0, 0, 'general') == 'outlier'       # forbidden centre


def test_glycine_is_more_permissive_than_general():
    # A point forbidden for a general residue but reachable by glycine symmetry.
    assert classify(150, 150, 'general') == 'outlier'
    assert classify(150, 150, 'glycine') in ('allowed', 'favored')


def test_proline_phi_is_restricted():
    assert classify(-63, -30, 'proline') == 'favored'
    assert classify(-63, 150, 'proline') == 'favored'
    # Proline cannot adopt a positive-phi (left-handed) conformation.
    assert classify(60, 60, 'proline') == 'outlier'


def test_real_structure_distribution_is_realistic():
    report = check_protein(DATA)
    s = report.ramachandran_summary
    assert s['total'] > 150
    # A well-refined 1.7 A trypsin should be mostly favored with few outliers.
    assert s['percent']['favored'] > 85
    assert s['percent']['outlier'] < 2.5
    assert s['counts']['favored'] + s['counts']['allowed'] + s['counts']['outlier'] == s['total']


def test_report_carries_ramachandran_and_coverage():
    report = check_protein(DATA)
    assert report.ramachandran, 'expected per-residue phi/psi points'
    assert any('Ramachandran' in line for line in report.coverage)
    assert any('MolProbity' in line for line in report.limitations)
    # Outliers are surfaced as reviewable issues that can highlight in the viewer.
    outlier_issues = [i for i in report.issues if i.code == 'ramachandran_outlier']
    assert len(outlier_issues) == report.ramachandran_summary['counts']['outlier']


def test_angles_only_computed_across_real_peptide_bonds():
    # Every computed residue must carry finite, in-range phi/psi.
    points = compute_ramachandran(check_protein_structure(DATA))
    assert points
    for p in points:
        assert -180 <= p.phi <= 180 and -180 <= p.psi <= 180
        assert p.classification in ('favored', 'allowed', 'outlier')
        assert p.residue_class in ('general', 'glycine', 'proline')


def test_render_svg_is_wellformed_and_plots_every_point():
    report = check_protein(DATA)
    svg = render_svg(report.ramachandran)
    assert svg.startswith('<svg') and svg.rstrip().endswith('</svg>')
    assert svg.count('<circle') >= len(report.ramachandran)  # points + legend dots
    for label in ('favored', 'allowed', 'outlier', 'phi', 'psi'):
        assert label in svg


def test_summarize_handles_empty():
    s = summarize([])
    assert s['total'] == 0 and s['percent'] == {'favored': 0.0, 'allowed': 0.0, 'outlier': 0.0}


def check_protein_structure(path):
    """Read the structure the same way check_protein does, for direct compute."""
    from mds.chem.receptor import read_structure
    return read_structure(path)
