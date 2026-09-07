"""Cell-level outliers, adapted from MetaCells 2 find_deviant_cells (gaps policy), fixed thresholds, one pass.
Tested genes = the label's HVGs. A cell is deviant on a gene when, inside its grain, a ≥ `fold` (log2) gap
separates it (with at most `allowed` other cells) from the rest. Genes that flag > `noisy_frac` of the block's
cells are bursty, not informative, and are ignored. Outliers are capped at `max_cell_frac` per block."""

from __future__ import annotations

import numpy as np


def gap_flags(V, fold, allowed):
    """V: cells×genes log2 values of one grain → bool cells×genes (deviant on gene via top or bottom tail)."""
    n, p = V.shape
    order = np.argsort(V, axis=0, kind="stable")
    Vs = np.take_along_axis(V, order, axis=0)
    big = (Vs[1:] - Vs[:-1]) >= fold  # gap i lies between sorted ranks i and i+1
    rank = np.arange(n - 1)
    top_ok = big & ((n - 1 - rank) <= allowed)[:, None]  # cells above gap i
    bot_ok = big & ((rank + 1) <= allowed)[:, None]  # cells at or below gap i
    flag_sorted = np.zeros((n, p), bool)
    for j in np.flatnonzero(top_ok.any(0)):
        flag_sorted[np.flatnonzero(top_ok[:, j])[0] + 1 :, j] = True
    for j in np.flatnonzero(bot_ok.any(0)):
        flag_sorted[: np.flatnonzero(bot_ok[:, j])[-1] + 1, j] = True
    flags = np.zeros((n, p), bool)
    np.put_along_axis(flags, order, flag_sorted, axis=0)
    return flags


def find_outliers(
    counts,
    tested,
    members_list,
    min_size=5,
    fold=3.0,
    max_gap_cells=3,
    max_gap_frac=0.1,
    noisy_frac=0.03,
    max_cell_frac=0.25,
    depth_cap=1e4,
    min_flag_genes=3,
    gap_quantile=0.999,
):
    """counts: csr block cells × all genes; tested: gene indices; members_list: local index arrays per grain.
    Returns (outlier bool per block cell, number of flagging genes per cell, fold threshold used).
    Two guards on top of MC2's rule, both needed for small (γ≈20) grains and non-UMI data:
    the fold threshold is max(fold, block-wide `gap_quantile` of all within-grain gaps), so technologies with
    noisier counts (Smart-seq2 amplification, dropouts) do not flag half the cells; and a cell must be flagged
    on ≥ `min_flag_genes` genes, since a single-gene burst is noise while a foreign cell deviates on many."""
    n_cells = counts.shape[0]
    totals = np.maximum(np.asarray(counts.sum(1)).ravel(), 1.0)
    depth = min(float(np.median(totals)), depth_cap)  # common depth before the +1 regularization
    X = counts[:, tested].toarray().astype(np.float32)
    V = np.log2(X * (depth / totals)[:, None].astype(np.float32) + 1.0)
    tested_mc = [m for m in members_list if len(m) >= min_size]
    if not tested_mc:
        return np.zeros(n_cells, bool), np.zeros(n_cells, int), fold
    gaps = np.concatenate([np.diff(np.sort(V[m], axis=0), axis=0).ravel() for m in tested_mc])
    fold_used = float(max(fold, np.quantile(gaps, gap_quantile)))
    flags = np.zeros(V.shape, bool)
    for m in tested_mc:
        n = len(m)
        allowed = max(max_gap_cells, int(np.floor(max_gap_frac * n)))
        flags[m] = gap_flags(V[m], fold_used, allowed)
    flags[:, flags.mean(0) > noisy_frac] = False
    n_flag = flags.sum(1)
    out = n_flag >= min_flag_genes
    cap = int(np.floor(max_cell_frac * n_cells))
    if out.sum() > cap:
        keep = np.argsort(-n_flag, kind="stable")[:cap]
        out = np.zeros(n_cells, bool)
        out[keep] = True
    return out, n_flag, fold_used
