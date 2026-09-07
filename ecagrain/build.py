"""Build: coordinates per coarse label pooled across samples; grouping inside sample × label blocks.
Grouping is SuperCell's algorithm: kNN graph → walktrap → cut into round(n/γ) communities."""

from __future__ import annotations

import igraph as ig
import numpy as np
import scanpy as sc
from anndata import AnnData
from sklearn.neighbors import NearestNeighbors

from .hvg import vst_hvg


def label_embedding(counts, samples, lateral, n_hvg, n_pcs, seed):
    """counts: csr cells×genes of one coarse label (all samples). Returns (coords, tested_gene_mask).
    HVG per sample when >1 sample (batch-only genes drop out); blocked genes (lateral + heat shock + mt/ribo/Malat1) are
    removed before ranking, so n_hvg genes are still selected."""
    a = AnnData(X=counts.copy())
    hv = vst_hvg(counts, n_hvg, batch=np.asarray(samples, dtype=str), exclude=lateral)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sub = a[:, hv].copy()
    sc.pp.scale(sub, max_value=10)
    n_comp = int(max(2, min(n_pcs, sub.n_obs // 5, sub.n_vars - 1)))
    # full LAPACK SVD: BLAS-3, ~6x faster than arpack on these dense n×2000 matrices, identical components
    sc.tl.pca(sub, n_comps=n_comp, random_state=seed, svd_solver="full")
    return np.asarray(sub.obsm["X_pca"]), hv


def knn_graph(coords, k):
    n = len(coords)
    k = min(k, n - 1)
    idx = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(coords, return_distance=False)
    src, dst = np.repeat(np.arange(n), k), idx[:, 1:].ravel()
    edges = np.unique(np.stack([np.minimum(src, dst), np.maximum(src, dst)], 1), axis=0)
    return ig.Graph(n=n, edges=edges.tolist())


def cut(g, n_clusters):
    """walktrap dendrogram cut; never fewer clusters than connected components."""
    n_comp = len(g.connected_components())
    return np.asarray(g.community_walktrap(steps=4).as_clustering(max(n_clusters, n_comp)).membership)


def partition_block(coords, gamma, k, cap=None):
    """Membership (0..m-1) for one block and its kNN graph (None when n < 3).
    γ is a target mean: n cells → max(1, round(n/γ)) groups, sizes decided by the tree;
    walktrap communities are unbalanced (1..2.5γ seen), so communities above `cap` are split in place."""
    n = len(coords)
    n_mc = max(1, int(n / gamma + 0.5))
    if n < 3:
        return np.zeros(n, int), None
    g = knn_graph(coords, k)
    memb = np.zeros(n, int) if n_mc == 1 else cut(g, n_mc)
    if cap:
        memb = enforce_cap(g, memb, gamma, cap)
    return memb, g


def enforce_cap(g, memb, gamma, cap):
    """Split every community larger than `cap` on its own subgraph until none exceeds it (each cut strictly shrinks)."""
    memb = memb.copy()
    while True:
        sizes = np.bincount(memb)
        big = np.flatnonzero(sizes > cap)
        if not len(big):
            return memb
        for c in big:
            members = np.flatnonzero(memb == c)
            parts = cut(
                g.induced_subgraph(members.tolist()),
                max(2, int(len(members) / gamma + 0.5)),
            )
            memb[members[parts > 0]] = memb.max() + parts[parts > 0]


def split_in_place(g, members):
    """Recheck: walktrap on the induced subgraph of one metacell, cut in two (more only if disconnected)."""
    return cut(g.induced_subgraph([int(m) for m in members]), 2)
