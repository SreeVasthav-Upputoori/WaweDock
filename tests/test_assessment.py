"""Tests for the consolidated protein-report assessment and score."""
from pathlib import Path

from mds.chem.health import build_assessment, check_protein

DATA = Path(__file__).parent / 'data' / '3PTB.pdb'


def test_real_structure_scores_well():
    report = check_protein(DATA)
    a = report.assessment
    assert 0 <= a['score'] <= 100
    assert a['grade'] in ('Excellent', 'Good', 'Fair', 'Poor')
    # A clean 1.7 A trypsin should score high.
    assert a['score'] >= 85
    m = a['metrics']
    assert m['residues'] > 150
    assert m['disulfides'] == 6
    # Metrics must match the underlying screens, not be recomputed differently.
    assert m['ramachandran_outliers'] == report.ramachandran_summary['counts']['outlier']
    assert m['rotamer_outliers'] == report.rotamers_summary['counts']['outlier']


def test_score_penalises_defects_and_clamps():
    report = check_protein(DATA)
    clean = report.assessment['score']

    class FakeIssue:
        def __init__(self, code): self.code = code
    # Inject many clashes/missing atoms: score must drop and never go below 0.
    report.issues = [FakeIssue('steric_overlap') for _ in range(500)] + \
                    [FakeIssue('missing_atoms') for _ in range(500)]
    degraded = build_assessment(report, report.assessment['metrics']['residues'])
    assert degraded['score'] < clean
    assert degraded['score'] >= 0
    assert degraded['grade'] == 'Poor'


def test_assessment_documents_itself_as_heuristic():
    report = check_protein(DATA)
    a = report.assessment
    assert 'formula' in a['score_formula'].lower() or '100' in a['score_formula']
    assert 'MolProbity' in a['note']


def test_zero_residues_does_not_divide_by_zero():
    report = check_protein(DATA)
    a = build_assessment(report, 0)
    assert 0 <= a['score'] <= 100
