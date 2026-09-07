"""Check the mcRigor port against the R package's outputs from the June Tabula Sapiens run
(supercell2.0/outputs/<tissue>/: input/counts.mtx, supercell/membership.tsv, mcrigor/mcRigor_TabMC.tsv, threshold.tsv).
T_org is deterministic and must agree; TT_div and the threshold carry permutation noise; verdicts should mostly agree.
The remaining systematic difference is HVG selection (scanpy seurat_v3 vs Seurat vst)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

from . import rigor


def validate(tissue_dir, nrep=1, seed=0, test_cutoff=0.05):
    d = Path(tissue_dir)
    counts = mmread(d / "input" / "counts.mtx").T.tocsr().astype(np.float32)  # genes×cells → cells×genes
    cells = pd.read_csv(d / "input" / "cells.tsv", sep="\t")["cell"].astype(str)
    memb = pd.read_csv(d / "supercell" / "membership.tsv", sep="\t").set_index("cell")["metacell"].reindex(cells)
    R = pd.read_csv(d / "mcrigor" / "mcRigor_TabMC.tsv", sep="\t").set_index("metacell_original")
    R_thre = pd.read_csv(d / "mcrigor" / "mcRigor_threshold.tsv", sep="\t")

    logdata, hv = rigor.lognorm_hvg(counts)
    L = logdata[:, hv].tocsr()
    rng = np.random.default_rng(seed)
    rows = {
        m: r
        for m, idx in memb.groupby(memb).indices.items()
        if (r := rigor.grain_stats(L, np.asarray(idx), rng, 0.1, nrep))
    }
    st = pd.DataFrame.from_dict(rows, orient="index")
    st["TT_div"] = rigor.tt_div(st)
    thre = rigor.threshold_table(st, test_cutoff)
    st["dubious"] = rigor.classify(st["size"], st["TT_div"], thre)

    c = st.join(
        R[["size", "T_org", "T_colperm", "TT_div", "mcRigor"]],
        rsuffix="_R",
        how="inner",
    )
    agree = (c["dubious"] == (c["mcRigor"] == "dubious")).mean()
    print(
        f"metacells compared: {len(c)} (python {len(st)}, R {len(R)}); sizes equal: {(c['size'] == c['size_R']).all()}"
    )
    print(
        f"T_org   : max|Δ| {np.abs(c['T_org'] - c['T_org_R']).max():.4g}  pearson {np.corrcoef(c['T_org'], c['T_org_R'])[0, 1]:.4f}"
    )
    print(
        f"TT_div  : max|Δ| {np.abs(c['TT_div'] - c['TT_div_R']).max():.4g}  pearson {np.corrcoef(c['TT_div'], c['TT_div_R'])[0, 1]:.4f}"
    )
    print(
        f"verdict : agreement {agree:.1%}  python dubious {int(c['dubious'].sum())}  R dubious {int((c['mcRigor'] == 'dubious').sum())}"
    )
    t = thre.merge(R_thre, on="size", how="inner", suffixes=("_py", "_R"))
    print(f"threshold at shared sizes ({len(t)}): max|Δ| {np.abs(t['thre_py'] - t['thre_R']).max():.4g}")
    return c, t
