"""Molecular fingerprint similarity and clustering (RDKit only — no new dependencies).

Uses Morgan (ECFP-like) fingerprints and Tanimoto similarity, the standard combination for
whole-molecule similarity in cheminformatics. Clustering is Butina's algorithm (RDKit
``ML.Cluster.Butina``): a distance-threshold clustering that does not require choosing a
number of clusters up front, appropriate for the small-to-medium ligand sets this tool handles.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from rdkit import Chem

from dataclasses import dataclass

import numpy as np

DEFAULT_RADIUS = 2       # ECFP4-equivalent
DEFAULT_N_BITS = 2048
DEFAULT_CLUSTER_CUTOFF = 0.4   # Tanimoto DISTANCE threshold (1 - similarity); smaller = tighter clusters

from functools import lru_cache

@lru_cache(maxsize=1)
def _fingerprint_generator():
    from rdkit.Chem import rdFingerprintGenerator
    return rdFingerprintGenerator.GetMorganGenerator(radius=DEFAULT_RADIUS, fpSize=DEFAULT_N_BITS)


@dataclass
class SimilarityResult:
    names: list[str]
    matrix: np.ndarray                    # (N, N) symmetric Tanimoto similarity, diagonal = 1.0
    clusters: list[int]                   # cluster id per ligand (0-indexed, by position in `names`)
    cluster_centroids: dict[int, str]      # cluster id -> name of its centroid (most central member)
    reference_name: str | None = None
    reference_similarity: list[float] | None = None   # Tanimoto similarity to the reference, per ligand


def fingerprint(mol: Chem.Mol):
    """Morgan fingerprint of the heavy-atom graph (stereochemistry-aware via default settings)."""
    from rdkit import Chem
    return _fingerprint_generator().GetFingerprint(Chem.RemoveHs(mol) if mol.GetNumAtoms() else mol)


def tanimoto(mol_a: Chem.Mol, mol_b: Chem.Mol) -> float:
    from rdkit import DataStructs
    return DataStructs.TanimotoSimilarity(fingerprint(mol_a), fingerprint(mol_b))


def similarity_matrix(mols: list[Chem.Mol]) -> np.ndarray:
    from rdkit import DataStructs
    fps = [fingerprint(m) for m in mols]
    n = len(fps)
    mat = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            s = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            mat[i, j] = mat[j, i] = s
    return mat


def cluster(mols: list[Chem.Mol], names: list[str],
           cutoff: float = DEFAULT_CLUSTER_CUTOFF) -> SimilarityResult:
    """Butina-cluster ``mols`` by Morgan/Tanimoto distance."""
    from rdkit.ML.Cluster import Butina
    from rdkit import DataStructs
    if len(mols) != len(names):
        raise ValueError("mols and names must be the same length")
    if not mols:
        return SimilarityResult([], np.zeros((0, 0)), [], {})

    fps = [fingerprint(m) for m in mols]
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1 - s for s in sims)

    raw_clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True) if n > 1 else ((0,),)

    cluster_of = [0] * n
    centroids: dict[int, str] = {}
    for cid, members in enumerate(raw_clusters):
        for idx in members:
            cluster_of[idx] = cid
        centroids[cid] = names[members[0]]   # Butina's first member is the cluster centroid

    mat = similarity_matrix(mols)
    return SimilarityResult(list(names), mat, cluster_of, centroids)


def compare_to_reference(mols: list[Chem.Mol], names: list[str], reference: Chem.Mol,
                         reference_name: str = "reference", *,
                         cutoff: float = DEFAULT_CLUSTER_CUTOFF) -> SimilarityResult:
    """Cluster ``mols`` and additionally rank them by similarity to ``reference``."""
    from rdkit import DataStructs
    result = cluster(mols, names, cutoff)
    ref_fp = fingerprint(reference)
    result.reference_name = reference_name
    result.reference_similarity = [DataStructs.TanimotoSimilarity(ref_fp, fingerprint(m)) for m in mols]
    return result


def most_similar_pairs(result: SimilarityResult, top_n: int = 10) -> list[tuple[str, str, float]]:
    """The ``top_n`` most similar distinct pairs, most similar first."""
    n = len(result.names)
    pairs = [(result.names[i], result.names[j], float(result.matrix[i, j]))
             for i in range(n) for j in range(i + 1, n)]
    pairs.sort(key=lambda t: t[2], reverse=True)
    return pairs[:top_n]
