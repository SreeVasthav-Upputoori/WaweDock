"""Tests for the side-chain rotamer (staggered chi-torsion) screen."""
from pathlib import Path

from mds.chem.health import check_protein
from mds.chem.rotamers import classify, compute_rotamers, summarize

DATA = Path(__file__).parent / 'data' / '3PTB.pdb'


def test_classify_from_deviations():
    assert classify([5, 12]) == 'favored'      # both near a staggered minimum
    assert classify([5, 35]) == 'allowed'      # one slightly off
    assert classify([55]) == 'outlier'         # eclipsed, far from +/-60/180


def test_real_structure_distribution_is_realistic():
    report = check_protein(DATA)
    s = report.rotamers_summary
    assert s['total'] > 120
    # A well-refined trypsin should be mostly favored with few outliers.
    assert s['percent']['favored'] > 85
    assert s['percent']['outlier'] < 6
    assert s['counts']['favored'] + s['counts']['allowed'] + s['counts']['outlier'] == s['total']


def test_non_rotameric_residues_are_skipped():
    points = compute_rotamers(_read(DATA))
    classified = {(p.chain, p.seqid) for p in points}
    # ALA, GLY and PRO have no screened chi angles and must not appear.
    for p in points:
        assert p.resname not in ('ALA', 'GLY', 'PRO')
        assert p.chis  # at least chi1 computed


def test_report_carries_rotamers_and_issues():
    report = check_protein(DATA)
    assert report.rotamers
    assert any('rotamer' in line.lower() for line in report.coverage)
    assert any('rotamer-library' in line or 'rotamer screen' in line.lower() for line in report.limitations)
    outlier_issues = [i for i in report.issues if i.code == 'rotamer_outlier']
    assert len(outlier_issues) == report.rotamers_summary['counts']['outlier']
    assert all(issue.atoms for issue in outlier_issues)  # highlightable in viewer


def test_summarize_handles_empty():
    s = summarize([])
    assert s['total'] == 0 and s['percent'] == {'favored': 0.0, 'allowed': 0.0, 'outlier': 0.0}


def _read(path):
    from mds.chem.receptor import read_structure
    return read_structure(path)
