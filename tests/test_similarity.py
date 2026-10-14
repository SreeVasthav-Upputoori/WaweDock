"""Similarity and clustering: validated against known drug chemistry."""

from __future__ import annotations

import pytest
from rdkit import Chem

from mds.chem.similarity import cluster, compare_to_reference, most_similar_pairs, tanimoto

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
SALICYLIC_ACID = "Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
NAPROXEN = "COc1ccc2cc(ccc2c1)C(C)C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"

NAMES = ["aspirin", "salicylic_acid", "ibuprofen", "naproxen", "caffeine", "paracetamol"]
SMILES = [ASPIRIN, SALICYLIC_ACID, IBUPROFEN, NAPROXEN, CAFFEINE, PARACETAMOL]


def _mols():
    return [Chem.MolFromSmiles(s) for s in SMILES]


def test_tanimoto_identity_and_ordering():
    mols = _mols()
    assert tanimoto(mols[0], mols[0]) == pytest.approx(1.0)
    close = tanimoto(mols[0], mols[1])          # aspirin vs salicylic acid: close analogues
    unrelated = tanimoto(mols[0], mols[4])      # aspirin vs caffeine: unrelated scaffolds
    assert close > unrelated
    assert 0.0 <= unrelated <= 1.0 and 0.0 <= close <= 1.0


def test_clustering_groups_structural_analogues():
    mols = _mols()
    # Loose cutoff (join if similarity >= 0.4): analogue pairs should share a cluster.
    result = cluster(mols, NAMES, cutoff=0.6)
    by_name = dict(zip(result.names, result.clusters))
    assert by_name["aspirin"] == by_name["salicylic_acid"]
    assert by_name["ibuprofen"] == by_name["naproxen"]
    assert by_name["caffeine"] != by_name["aspirin"]
    assert len(result.cluster_centroids) == len(set(result.clusters))


def test_clustering_tight_cutoff_separates_everything():
    mols = _mols()
    result = cluster(mols, NAMES, cutoff=0.1)
    assert len(set(result.clusters)) == len(NAMES)


def test_compare_to_reference_ranks_by_similarity():
    mols = _mols()
    reference = Chem.MolFromSmiles(ASPIRIN)
    result = compare_to_reference(mols, NAMES, reference, "aspirin")
    assert result.reference_similarity[NAMES.index("aspirin")] == pytest.approx(1.0)
    ranked = sorted(zip(NAMES, result.reference_similarity), key=lambda t: -t[1])
    assert ranked[0][0] == "aspirin"
    assert ranked[1][0] == "salicylic_acid"   # closest analogue, should rank second


def test_most_similar_pairs_are_sorted_descending():
    mols = _mols()
    result = cluster(mols, NAMES, cutoff=0.4)
    pairs = most_similar_pairs(result, top_n=5)
    scores = [s for _, _, s in pairs]
    assert scores == sorted(scores, reverse=True)
    assert {"aspirin", "salicylic_acid"} == set(pairs[0][:2])


def test_mismatched_lengths_raise():
    with pytest.raises(ValueError):
        cluster(_mols(), NAMES[:-1])


def test_single_and_empty_inputs():
    assert cluster([], []).names == []
    one = cluster(_mols()[:1], NAMES[:1])
    assert one.clusters == [0]
