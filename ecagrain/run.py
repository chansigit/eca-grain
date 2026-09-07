"""One-through pipeline for one released unit. Fails only on missing input or a failed conservation check.
Blocks = sample × label (default label = zmip_ann_coarse); coordinates per label pooled across samples.
Stages: load → ① build → ② outliers → ③ diagnose → ④ recheck → ⑤ deliver (conservation, files, report page)."""

from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from . import __version__, build, outlier, rigor
from .genes import blocked_mask

DEFAULTS = dict(
    gamma=20,
    k=15,
    max_size_factor=1.5,
    n_hvg=2000,
    n_pcs=30,
    min_pool=100,
    min_test_size=5,
    nrep=20,
    test_cutoff=0.05,
    gene_filter=0.1,
    seed=0,
)
COLS = dict(
    sample_col="eca_sample_id",
    label_col="zmip_ann_coarse",  # block boundary (user decision 2026-09-07: coarse lineage)
    audit_col="zmip_ann_fine",  # recorded per metacell (majority + purity), never used for grouping
    counts_layer="counts",
    embed_key="X_pca_harmony",
)


def load(h5ad, sample_col, label_col, audit_col, counts_layer, embed_key):
    a = ad.read_h5ad(h5ad)
    for c in (sample_col, label_col, audit_col):
        if c not in a.obs:
            raise SystemExit(f"missing obs column {c!r} in {h5ad}")
    if counts_layer != "X" and counts_layer not in a.layers:
        raise SystemExit(f"missing layer {counts_layer!r} in {h5ad}")
    counts = sparse.csr_matrix(a.X if counts_layer == "X" else a.layers[counts_layer], dtype=np.float32)
    obs = pd.DataFrame(
        {
            "cell": a.obs_names.astype(str),
            "sample": a.obs[sample_col].astype(str).to_numpy(),
            "label": a.obs[label_col].astype(str).to_numpy(),
            "audit": a.obs[audit_col].astype(str).to_numpy(),
        }
    )
    if obs["sample"].str.contains("|", regex=False).any():
        raise SystemExit(f"sample ids must not contain '|' (block separator): {sample_col!r} in {h5ad}")
    obs["block"] = obs["sample"] + "|" + obs["label"]
    gpca = np.asarray(a.obsm[embed_key]) if embed_key in a.obsm else None
    ghv = a.var["highly_variable"].to_numpy() if "highly_variable" in a.var else None
    umap = np.asarray(a.obsm["X_umap"]) if "X_umap" in a.obsm else None
    return counts, obs, a.var.copy(), gpca, ghv, umap


def build_stage(counts, obs, blocked, gpca, ghv, p, embed_key):
    """① Coordinates per label (all samples pooled), groups per sample × label block.
    Returns build id per cell, block of each build id, tested genes per label, (kNN graph, cells) per block, block table."""
    mc_of = np.full(len(obs), -1)
    mc_block, tested_of, graphs, blocks = [], {}, {}, []
    cap = int(p["max_size_factor"] * p["gamma"])
    for label, pool in obs.groupby("label").indices.items():
        pool = np.asarray(pool)
        if len(pool) >= p["min_pool"]:
            coords, hv = build.label_embedding(
                counts[pool], obs["sample"].to_numpy()[pool], blocked, p["n_hvg"], p["n_pcs"]
            )
            src = "label"
        else:
            if gpca is None or ghv is None:
                raise SystemExit(f"label {label!r} has {len(pool)} cells < min_pool and no {embed_key}")
            coords, hv, src = gpca[pool], ghv & ~blocked, "global"
        tested_of[label] = np.flatnonzero(hv)
        pos = {c: i for i, c in enumerate(pool)}
        for block, idx in obs.iloc[pool].groupby("block").indices.items():
            cells = pool[np.asarray(idx)]
            memb, g = build.partition_block(coords[[pos[c] for c in cells]], p["gamma"], p["k"], cap=cap)
            graphs[block] = (g, cells)
            mc_of[cells] = memb + len(mc_block)
            n_mc = int(memb.max()) + 1
            mc_block += [block] * n_mc
            blocks.append(
                dict(
                    block=block,
                    sample=block.split("|", 1)[0],
                    label=label,
                    embedding=src,
                    n_cells=len(cells),
                    n_metacells_build=n_mc,
                )
            )
    return mc_of, mc_block, tested_of, graphs, pd.DataFrame(blocks).set_index("block")


def outlier_stage(counts, obs, mc_of, graphs, tested_of, blocks):
    """② MC2 gaps rule with guards, per block; tested genes = the label's HVGs. Adds fold and counts to `blocks`."""
    is_out, n_flag, fold_used = np.zeros(len(obs), bool), np.zeros(len(obs), int), {}
    for block, (_g, cells) in graphs.items():
        mcs = mc_of[cells]
        members = [np.flatnonzero(mcs == m) for m in np.unique(mcs)]
        o, nf, fold_used[block] = outlier.find_outliers(counts[cells], tested_of[block.split("|", 1)[1]], members)
        is_out[cells], n_flag[cells] = o, nf
    blocks["outlier_fold"] = pd.Series(fold_used)
    blocks["n_outliers"] = pd.Series(is_out).groupby(obs["block"].to_numpy()).sum()
    return is_out, n_flag


def grain_stats(L, groups, ids, rng, p):
    """③ mcRigor statistics per grain id (size only when below min_test_size)."""
    rows = {}
    for m in ids:
        n = len(groups[m])
        r = rigor.metacell_stats(L, groups[m], rng, p["gene_filter"], p["nrep"]) if n >= p["min_test_size"] else None
        rows[m] = r or {"size": n}
    st = pd.DataFrame.from_dict(rows, orient="index")
    st["TT_div"] = rigor.tt_div(st) if "T_org" in st else 1.0
    return st


def recheck_stage(st, groups, graphs, mc_block, mc_of, floor, min_test_size):
    """④ Split dubious grains in place (walktrap on their own subgraph), floor γ/2, one level.
    Children get new ids appended to `mc_block`; returns the status of every build id."""
    status = {}
    for m in st.index:
        if not st.loc[m, "dubious"]:
            status[m] = "trustworthy" if len(groups[m]) >= min_test_size else "untested"
            continue
        g, cells = graphs[mc_block[m]]
        members = groups[m]
        if g is None or len(members) < 2 * floor:
            status[m] = "residual_dubious"
            continue
        pos = {c: i for i, c in enumerate(cells)}
        parts = build.split_in_place(g, [pos[c] for c in members])
        sizes = np.bincount(parts)
        if sizes.min() < floor:
            status[m] = "residual_dubious"
            continue
        for h in range(len(sizes)):
            mc_of[members[parts == h]] = len(mc_block)
            mc_block.append(mc_block[m])
        status[m] = "split"
    return status


def run(h5ad, outdir, **kw):
    p = {**DEFAULTS, **{k: v for k, v in kw.items() if k in DEFAULTS}}
    cols = {**COLS, **{k: v for k, v in kw.items() if k in COLS}}
    t0 = time.time()

    def log(msg):
        print(f"[ecagrain +{time.time() - t0:6.1f}s] {msg}", flush=True)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(p["seed"])
    counts, obs, var, gpca, ghv, umap = load(h5ad, **cols)
    n_cells = len(obs)
    blocked = blocked_mask(var)  # never drive grouping: cell cycle, stress, heat shock, mt, ribo, Malat1
    floor = p["gamma"] // 2
    log(
        f"loaded {n_cells} cells × {counts.shape[1]} genes, {obs['sample'].nunique()} samples, {obs['label'].nunique()} labels"
    )

    mc_of, mc_block, tested_of, graphs, blocks = build_stage(counts, obs, blocked, gpca, ghv, p, cols["embed_key"])
    log(f"① build: {len(mc_block)} grains in {len(blocks)} blocks")

    is_out, n_flag = outlier_stage(counts, obs, mc_of, graphs, tested_of, blocks)
    build_id = mc_of.copy()
    mc_of[is_out] = -1
    log(f"② outliers: {int(is_out.sum())} ({is_out.mean():.1%})")

    logdata, hv_unit = rigor.lognorm_hvg(counts, p["n_hvg"])
    L = logdata[:, hv_unit].tocsr()
    groups = pd.Series(np.arange(n_cells)).groupby(mc_of).indices  # build-id -> cell idx
    ids0 = [m for m in sorted(groups) if m >= 0]
    st = grain_stats(L, groups, ids0, rng, p)
    st["dubious"], thre = rigor.flag_dubious(st, p["test_cutoff"])
    st["level"], st["parent"] = 0, -1
    n_dub0 = int(st["dubious"].sum())
    log(f"③ diagnose: {n_dub0} dubious of {len(ids0)} ({n_dub0 / max(len(ids0), 1):.1%})")

    status = recheck_stage(st, groups, graphs, mc_block, mc_of, floor, p["min_test_size"])
    groups = pd.Series(np.arange(n_cells)).groupby(mc_of).indices
    new_ids = [m for m in sorted(groups) if m >= 0 and m not in status]
    if new_ids:
        st1 = grain_stats(L, groups, new_ids, rng, p)
        st1["dubious"] = rigor.classify(st1["size"], st1["TT_div"], thre)
        st1["level"] = 1
        st1["parent"] = [build_id[groups[m][0]] for m in new_ids]
        for m in new_ids:  # children are ≥ floor ≥ min_test_size, so always tested
            status[m] = "residual_dubious" if st1.loc[m, "dubious"] else "trustworthy"
        st = pd.concat([st, st1])
    log(f"④ recheck: split {sum(s == 'split' for s in status.values())}, {len(new_ids)} children")

    # ---- ⑤ aggregate + conservation + write
    final = [m for m in sorted(groups) if m >= 0 and status[m] != "split"]
    row_of = {m: i for i, m in enumerate(final)}
    row = np.array([row_of.get(m, -1) for m in mc_of])
    assigned = row >= 0
    if not np.array_equal(assigned, ~is_out):
        raise SystemExit("conservation failed: assigned cells != non-outlier cells")
    M = sparse.csr_matrix(
        (np.ones(assigned.sum()), (row[assigned], np.flatnonzero(assigned))), shape=(len(final), n_cells)
    )
    size = np.asarray(M.sum(1)).ravel().astype(int)
    if size.sum() + is_out.sum() != n_cells:
        raise SystemExit("conservation failed: sum(size) + outliers != n_cells")
    per_block = pd.Series(size).groupby(pd.Series([mc_block[m] for m in final])).sum()
    check = per_block.add(blocks["n_outliers"], fill_value=0).reindex(blocks.index).fillna(0)
    if not np.array_equal(check.to_numpy().astype(int), blocks["n_cells"].to_numpy()):
        raise SystemExit("conservation failed inside a block")

    mc_id = np.array([f"mc{i:05d}" for i in range(len(final))])
    mem = obs.assign(row=row)
    aud = mem[assigned].groupby("row")["audit"].agg(lambda s: s.value_counts().idxmax())
    purity = mem[assigned].groupby("row")["audit"].agg(lambda s: s.value_counts().iloc[0] / len(s))
    mobs = pd.DataFrame(
        {
            "size": size,
            "block": [mc_block[m] for m in final],
            "sample": [mc_block[m].split("|", 1)[0] for m in final],
            "label": [mc_block[m].split("|", 1)[1] for m in final],
            "audit_majority": aud.reindex(range(len(final))).to_numpy(),
            "audit_purity": purity.reindex(range(len(final))).to_numpy(),
            "level": st.loc[final, "level"].to_numpy(),
            "build_id": final,
            "parent_build_id": st.loc[final, "parent"].to_numpy(),
            "mcRigor": [status[m] for m in final],
            "TT_div": st.loc[final, "TT_div"].to_numpy(),
            "n_test_genes": st.loc[final, "n_genes"].to_numpy() if "n_genes" in st else np.nan,
            "gamma": p["gamma"],
        },
        index=mc_id,
    )
    mca = ad.AnnData(X=(M @ counts).tocsr(), obs=mobs, var=var)
    Mn = sparse.diags(1.0 / size) @ M
    if umap is not None:
        mca.obsm["X_umap_mean"] = np.asarray(Mn @ umap)
    if gpca is not None:
        mca.obsm[f"{cols['embed_key']}_mean"] = np.asarray(Mn @ gpca)
    mca.uns["ecagrain"] = {"version": __version__, "params": p, "columns": cols, "input": str(h5ad)}
    mca.write_h5ad(outdir / "metacells.h5ad")

    membership = obs.assign(
        metacell_id=np.where(assigned, mc_id[np.maximum(row, 0)], None),
        status=np.where(is_out, "outlier", "member"),
        n_flag_genes=n_flag,
        build_id=build_id,
        level=np.where(assigned, st["level"].reindex(mc_of).fillna(-1).to_numpy(), -1).astype(int),
        mcRigor=[status.get(m) for m in mc_of],
        umap_1=umap[:, 0] if umap is not None else np.nan,  # cell UMAP kept here so the report page needs no input
        umap_2=umap[:, 1] if umap is not None else np.nan,
    )
    membership.to_parquet(outdir / "membership.parquet", index=False)
    thre.to_csv(outdir / "threshold.tsv", sep="\t", index=False)

    mc_by_block = mobs.groupby("block")
    blocks["n_metacells_final"] = mc_by_block.size()
    blocks["n_residual_dubious"] = mc_by_block["mcRigor"].apply(lambda s: int((s == "residual_dubious").sum()))
    blocks = blocks.fillna({"n_metacells_final": 0, "n_residual_dubious": 0})
    summary = {
        "version": __version__,
        "input": str(h5ad),
        "params": p,
        "columns": cols,
        "n_cells": int(n_cells),
        "n_blocks": int(len(blocks)),
        "n_labels": int(obs["label"].nunique()),
        "n_samples": int(obs["sample"].nunique()),
        "n_metacells_build": int(len(ids0)),
        "n_outliers": int(is_out.sum()),
        "outlier_rate": float(is_out.mean()),
        "n_dubious_build": n_dub0,
        "dubious_rate_build": n_dub0 / max(len(ids0), 1),
        "n_split": int(sum(s == "split" for s in status.values())),
        "n_metacells_final": int(len(final)),
        "n_residual_dubious": int(sum(status[m] == "residual_dubious" for m in final)),
        "n_untested": int(sum(status[m] == "untested" for m in final)),
        "cells_by_status": {
            s: int(size[[status[m] == s for m in final]].sum()) for s in ("trustworthy", "residual_dubious", "untested")
        },
        "size_quantiles": {q: float(np.quantile(size, q)) for q in (0, 0.1, 0.5, 0.9, 1)},
        "audit_purity_median": float(np.nanmedian(mobs["audit_purity"])),
        "conservation": "ok",
        "elapsed_s": round(time.time() - t0, 1),
        "blocks": blocks.reset_index().to_dict(orient="records"),
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    log(f"⑤ deliver: {len(final)} grains, conservation ok, files written")
    from .figures import write_html  # per-unit page = the agreed template (two figures)

    write_html([outdir], outdir / "report.html")
    summary["elapsed_pipeline_s"], summary["elapsed_s"] = summary["elapsed_s"], round(time.time() - t0, 1)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (outdir / "report.md").write_text(report(summary, blocks))
    log("report page written")
    return summary


def report(s, blocks):
    L = [
        f"# ecagrain {s['version']} · {s['input']}",
        "",
        f"- cells {s['n_cells']} · samples {s['n_samples']} · labels {s['n_labels']} · blocks {s['n_blocks']}",
        f"- final metacells {s['n_metacells_final']} (size 0/10/50/90/100%: "
        + " / ".join(f"{v:.0f}" for v in s["size_quantiles"].values())
        + ")",
        f"- outliers {s['n_outliers']} ({s['outlier_rate']:.1%}); audit-label purity median {s['audit_purity_median']:.3f}; "
        f"conservation {s['conservation']}; {s['elapsed_s']} s",
        f"- internals: built {s['n_metacells_build']}, dubious at build {s['n_dubious_build']} ({s['dubious_rate_build']:.1%}), "
        f"split {s['n_split']}, residual dubious {s['n_residual_dubious']}, untested {s['n_untested']}; cells by status "
        + ", ".join(f"{k} {v} ({v / s['n_cells']:.1%})" for k, v in s["cells_by_status"].items()),
        "",
        "| block | embedding | cells | outliers | metacells |",
        "|---|---|---:|---:|---:|",
    ]
    for b, r in blocks.iterrows():
        L.append(
            f"| {b} | {r['embedding']} | {r['n_cells']} | {int(r['n_outliers'])} | {int(r['n_metacells_final'])} |"
        )
    return "\n".join(L) + "\n"
