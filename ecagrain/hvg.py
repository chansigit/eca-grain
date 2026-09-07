"""Seurat-vst style highly variable genes on raw counts, without scikit-misc (its loess binary is broken in the venv).
Trend log10(var) ~ log10(mean) fitted with statsmodels lowess (local linear, span 0.3) instead of loess (local quadratic);
standardized variance with values clipped at sqrt(N). With `batch`, genes are ranked per batch and combined the
scanpy way (number of batches where the gene is in the top n, then median rank)."""

from __future__ import annotations

import numpy as np
from scipy import sparse
from statsmodels.nonparametric.smoothers_lowess import lowess


def standardized_variance(X, span=0.3):
    """X: csr cells × genes raw counts → norm_gene_var per gene (0 for constant genes)."""
    X = sparse.csr_matrix(X, dtype=np.float64)
    n = X.shape[0]
    mean = np.asarray(X.mean(0)).ravel()
    sq = np.asarray(X.multiply(X).mean(0)).ravel()
    var = (sq - mean**2) * n / (n - 1)
    ok = var > 0
    out = np.zeros(X.shape[1])
    if ok.sum() < 10:
        return out
    x, y = np.log10(mean[ok]), np.log10(var[ok])
    fit = lowess(y, x, frac=span, it=0, delta=0.01 * (x.max() - x.min()), return_sorted=True)
    reg_std = np.sqrt(10 ** np.interp(x, fit[:, 0], fit[:, 1]))
    clip = mean[ok] + reg_std * np.sqrt(n)
    Xc = X[:, ok].copy()
    Xc.data = np.minimum(Xc.data, clip[Xc.indices])
    s1 = np.asarray(Xc.sum(0)).ravel()
    s2 = np.asarray(Xc.multiply(Xc).sum(0)).ravel()
    m = mean[ok]
    out[ok] = (n * m**2 + s2 - 2 * s1 * m) / ((n - 1) * reg_std**2)
    return out


def vst_hvg(X, n_top=2000, batch=None, min_cells_per_batch=10):
    """Boolean mask of the n_top most variable genes."""
    X = sparse.csr_matrix(X)
    if batch is None or len(np.unique(batch)) < 2:
        v = standardized_variance(X)
        mask = np.zeros(X.shape[1], bool)
        mask[np.argsort(-v, kind="stable")[:n_top]] = True
        return mask
    batch = np.asarray(batch)
    ranks, in_top = [], []
    for b in np.unique(batch):
        idx = np.flatnonzero(batch == b)
        if len(idx) < min_cells_per_batch:
            continue
        v = standardized_variance(X[idx])
        r = np.empty(len(v))
        r[np.argsort(-v, kind="stable")] = np.arange(len(v))
        ranks.append(r)
        in_top.append(r < n_top)
    if not ranks:
        return vst_hvg(X, n_top)
    nb, med = np.sum(in_top, 0), np.median(ranks, 0)
    order = np.lexsort((med, -nb))  # most batches first, then best median rank
    mask = np.zeros(X.shape[1], bool)
    mask[order[:n_top]] = True
    return mask
