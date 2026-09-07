"""Build: coordinates per coarse label pooled across samples; grouping inside sample × label blocks.
Grouping is SuperCell's algorithm: kNN graph → walktrap → cut into round(n/γ) communities.
Heavy products go through scipy's BLAS (`scipy.linalg.blas`): numpy in some environments carries no BLAS at all."""

from __future__ import annotations

import igraph as ig
import numpy as np
import scipy.linalg
from scipy import sparse
from scipy.linalg import blas
from sklearn.neighbors import NearestNeighbors

from .hvg import vst_hvg


def lognorm(counts):
    """CP10k + log1p on a csr counts matrix (= scanpy normalize_total(target_sum=1e4) + log1p)."""
    counts = sparse.csr_matrix(counts, dtype=np.float32)
    tot = np.asarray(counts.sum(1)).ravel()
    tot[tot == 0] = 1.0
    return (sparse.diags((1e4 / tot).astype(np.float32)) @ counts).log1p().tocsr()


def pca_scores(L, n_comp, max_value=10.0, chunk=20000):
    """PCA scores of the scaled matrix without ever densifying it.
    z = clip((x − mean) / std, ±max_value) per gene (scanpy `scale(max_value=10)`, std with ddof 1, constant genes
    untouched), then scores = (z − mean_z) V with V the top eigenvectors of the covariance of z. Two passes over
    row chunks; memory is chunk × p plus p × p. Equal to sc.pp.scale + sc.tl.pca(svd_solver="full") up to sign."""
    n, p = L.shape
    mean = np.asarray(L.mean(0)).ravel()
    sq = np.asarray(L.multiply(L).mean(0)).ravel()
    var = (sq - mean**2) * (n / (n - 1) if n > 1 else 1.0)
    std = np.sqrt(np.maximum(var, 0.0))
    std[std == 0] = 1.0

    def z(i):
        x = L[i : i + chunk].toarray().astype(np.float64)
        x -= mean
        x /= std
        return np.clip(x, -max_value, max_value, out=x)

    S, s = np.zeros((p, p), order="F"), np.zeros(p)
    for i in range(0, n, chunk):
        x = z(i)
        S = blas.dsyrk(1.0, x, beta=1.0, c=S, trans=1, overwrite_c=1)  # upper triangle of Σ xᵀx
        s += x.sum(0)
    S = np.triu(S) + np.triu(S, 1).T
    mu = s / n
    C = (S - n * np.outer(mu, mu)) / (n - 1)
    _, V = scipy.linalg.eigh(C, subset_by_index=[p - n_comp, p - 1])
    V = np.ascontiguousarray(V[:, ::-1])  # descending variance
    out = np.empty((n, n_comp), np.float32)
    for i in range(0, n, chunk):
        out[i : i + chunk] = blas.dgemm(1.0, z(i) - mu, V)
    return out


def label_embedding(counts, samples, lateral, n_hvg, n_pcs):
    """counts: csr cells×genes of one coarse label (all samples). Returns (coords, tested_gene_mask).
    HVG per sample when >1 sample (batch-only genes drop out); blocked genes (lateral + heat shock + mt/ribo/Malat1) are
    removed before ranking, so n_hvg genes are still selected."""
    hv = vst_hvg(counts, n_hvg, batch=np.asarray(samples, dtype=str), exclude=lateral)
    L = lognorm(counts)[:, hv]
    n_comp = int(max(2, min(n_pcs, L.shape[0] // 5, L.shape[1] - 1)))
    return pca_scores(L, n_comp), hv


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
            parts = cut(g.induced_subgraph(members.tolist()), max(2, int(len(members) / gamma + 0.5)))
            memb[members[parts > 0]] = memb.max() + parts[parts > 0]


def split_in_place(g, members):
    """Recheck: walktrap on the induced subgraph of one metacell, cut in two (more only if disconnected)."""
    return cut(g.induced_subgraph([int(m) for m in members]), 2)
