"""mcRigor DETECT (Liu & Li, Nat Commun 2025) re-implemented from the R package source
(R/mcRigor_function.R + src/mc_test_stats.cpp). No R involved.

Per grain: z-score its cells × HVG matrix per gene, T = ||corr − I||_F / sqrt(p (p − 0.5)).
TT_div = T(data) / T(each gene column shuffled across cells)  — the grain's heterogeneity score.
Null    = T(each cell row shuffled across genes) / T(that, column-shuffled), Nrep times.
Threshold per size = (1 − cutoff) quantile of the null, lowess-smoothed over size; TT_div above it = dubious.
Deviation from R kept on purpose: R uses ALL genes instead of the HVGs whenever n_genes <= n_cells; we always use HVGs."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.linalg import blas
from statsmodels.nonparametric.smoothers_lowess import lowess

from .hvg import vst_hvg


def lognorm_hvg(counts, n_hvg=2000):
    """Seurat NormalizeData (CP10k, log1p) + FindVariableFeatures vst on the whole unit → (logdata csr, hvg idx)."""
    a = AnnData(X=counts.copy())
    hv = np.flatnonzero(vst_hvg(counts, n_hvg))
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a.X.tocsr(), hv


def t_stat(dat):
    """T = ||corr − I||_F / sqrt(p (p − 0.5)) with corr = cᵀc/(n−1), evaluated through the n×n Gram matrix G = c cᵀ:
    ||cᵀc||_F = ||c cᵀ||_F and tr(cᵀc) = tr(G), so no p×p matrix is formed (n²p instead of np² work; equal to 1e-15)."""
    n, p = dat.shape
    c = dat - dat.mean(0)
    G = blas.dgemm(1.0, c, c, trans_b=1)  # scipy BLAS: numpy may have none
    fro2 = (G * G).sum() / (n - 1) ** 2 - 2.0 * np.trace(G) / (n - 1) + p
    return float(np.sqrt(max(fro2, 0.0)) / np.sqrt(p * (p - 0.5)))


def scale_cols(m):
    c = m - m.mean(0)
    sd = np.sqrt((c**2).sum(0) / (m.shape[0] - 1))
    nz = sd > 0
    c[:, nz] /= sd[nz]
    return c


def grain_stats(logdata_hvg, members, rng, gene_filter=0.1, nrep=1):
    """logdata_hvg: csr cells × hvg. None when size < 2 or fewer than 2 genes pass the expressed-in-10% filter."""
    n = len(members)
    if n < 2:
        return None
    sub = logdata_hvg[members].toarray()
    keep = (sub > 0).sum(0) > n * gene_filter
    if keep.sum() < 2:
        return None
    dat = scale_cols(sub[:, keep].astype(np.float64))
    res = {
        "size": n,
        "n_genes": int(keep.sum()),
        "T_org": t_stat(dat),
        "T_colperm": t_stat(rng.permuted(dat, axis=0)),
    }
    for r in range(1, nrep + 1):
        x = rng.permuted(dat, axis=1)
        res[f"T_rowperm{r}"] = t_stat(x)
        res[f"T_bothperm{r}"] = t_stat(rng.permuted(x, axis=0))
    return res


def tt_div(df):
    return (df["T_org"] / df["T_colperm"]).fillna(1.0)


def threshold_table(df, test_cutoff=0.05, frac=1 / 6):
    """(size, thre): per size the (1 − cutoff) quantile of T_rowperm/T_bothperm, then R-style lowess (f=1/6, iter=3)."""
    rows = sorted(c for c in df if c.startswith("T_rowperm"))
    boths = sorted(c for c in df if c.startswith("T_bothperm"))
    null = df[rows].to_numpy(float) / df[boths].to_numpy(float)
    null = np.where(np.isfinite(null), null, 1.0)
    sizes = df["size"].to_numpy()
    tab = [(s, np.quantile(null[sizes == s].ravel(), 1 - test_cutoff)) for s in np.unique(sizes) if s >= 2]
    thre = pd.DataFrame(tab, columns=["size", "thre"])
    if len(thre) >= 3:
        x, y = thre["size"].to_numpy(float), thre["thre"].to_numpy(float)
        sm = lowess(y, x, frac=frac, it=3, delta=0.01 * (x.max() - x.min()), return_sorted=True)
        thre = pd.DataFrame({"size": sm[:, 0], "thre": sm[:, 1]})
    return thre


def flag_dubious(st, test_cutoff=0.05):
    """(dubious bool Series, threshold table) for a stats table. Grains without T_org (untested) are never dubious;
    a table without any tested grain gets an empty threshold table."""
    tested = st["T_org"].notna() if "T_org" in st else pd.Series(False, index=st.index)
    if not tested.any():
        return pd.Series(False, index=st.index), pd.DataFrame(columns=["size", "thre"])
    thre = threshold_table(st[tested], test_cutoff)
    return pd.Series(classify(st["size"], st["TT_div"], thre), index=st.index) & tested, thre


def classify(size, ttd, thre):
    """Dubious when TT_div exceeds the threshold at that size (interpolated, edge-clamped); size 1 never dubious."""
    size = np.asarray(size, float)
    t = np.interp(size, thre["size"].to_numpy(float), thre["thre"].to_numpy(float))
    return (np.asarray(ttd, float) > t) & (size > 1)
