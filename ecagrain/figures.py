"""Self-contained results page showing the final grains only, two figures per unit:
1. UMAP panel — left: grains (opaque, dot size ∝ cells) on the raw cells (translucent), both coloured by coarse
   lineage; right: UMAP computed on the grains, Leiden clusters on top, coarse-lineage islands (concave hulls,
   translucent) underneath; one shared lineage legend.
2. Split marker heatmap — rows = de-duplicated top-5 Wilcoxon markers of the grain clusters; columns split by
   cluster, hierarchically ordered inside each split, splits ordered by clustering the split means; strip = lineage.
python -m ecagrain figures <run_dir>... --out results.html"""

from __future__ import annotations

import base64
import html as html_mod
import io
import json
from pathlib import Path

import anndata as ad
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage

from .hvg import vst_hvg

MAX_LEGEND = 20  # more lineages than this → numbered centroids + a table
LEIDEN_RES = 1.0
TOP_GENES = 5
HULL_RATIO = 0.2
F_TITLE, F_LEGEND, F_GENE, F_STRIP, F_NUM = 12, 10, 10, 9, 10  # font sizes


def _png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _cell_umap(mem, input_h5ad):
    """Cell UMAP from membership.parquet (written since 0.3.0); older run dirs fall back to the input h5ad."""
    if "umap_1" in mem and mem["umap_1"].notna().any():
        return mem[["umap_1", "umap_2"]].to_numpy(float)
    try:
        with h5py.File(input_h5ad) as f:
            return np.asarray(f["obsm/X_umap"]) if "obsm/X_umap" in f else None
    except OSError:
        return None


def _short(s, n=40):
    return s if len(s) <= n else s[: n - 1] + "…"


def _number(ax, labels, xy, lab_of_point):
    """Number each lineage at its centroid (used only when there are too many lineages for a legend)."""
    for j, c in enumerate(labels, 1):
        i = lab_of_point == c
        if i.any():
            ax.text(
                np.median(xy[i, 0]),
                np.median(xy[i, 1]),
                str(j),
                fontsize=F_NUM,
                fontweight="bold",
                ha="center",
                va="center",
                bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="none", alpha=0.75),
                zorder=7,
            )


def _hier_order(Z):
    """Leaf order from average-linkage hierarchical clustering; identity for < 3 items."""
    if len(Z) < 3:
        return np.arange(len(Z))
    return leaves_list(linkage(Z, method="average", metric="euclidean"))


def _islands(ax, pts, color, ratio=HULL_RATIO):
    """Split `pts` into islands (connected components at 6× the median nearest-neighbour distance) and draw each
    island with ≥ 4 points as a concave hull: translucent fill (padded 2.2% of the axis width) underneath,
    outline (padded 1.2%) above it."""
    import shapely
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    if len(pts) < 4:
        return
    tree = cKDTree(pts)
    r = 6.0 * np.median(tree.query(pts, k=2)[0][:, 1])
    pairs = np.array(list(tree.query_pairs(r)))
    if len(pairs):
        adj = coo_matrix(
            (np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
            shape=(len(pts), len(pts)),
        )
        comp = connected_components(adj, directed=False)[1]
    else:
        comp = np.arange(len(pts))
    pad = 0.012 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    for c in np.unique(comp):
        idx = np.flatnonzero(comp == c)
        if len(idx) < 4:
            continue
        core = shapely.concave_hull(shapely.MultiPoint(pts[idx]), ratio=ratio, allow_holes=False)
        for geom, kind in (
            (core.buffer(pad * 1.8), "fill"),
            (core.buffer(pad), "line"),
        ):
            polys = [geom] if geom.geom_type == "Polygon" else list(getattr(geom, "geoms", []))
            for poly in polys:
                if poly.geom_type != "Polygon":
                    continue
                x, y = poly.exterior.xy
                if kind == "fill":
                    ax.fill(x, y, color=color, alpha=0.22, lw=0, zorder=0)
                else:
                    ax.plot(
                        x,
                        y,
                        color=color,
                        lw=0.9,
                        alpha=0.9,
                        solid_capstyle="round",
                        zorder=1,
                    )


def figures(run_dir):
    d = Path(run_dir)
    s = json.loads((d / "summary.json").read_text())
    mc = ad.read_h5ad(d / "grains.h5ad")
    mo = mc.obs
    mem = pd.read_parquet(d / "membership.parquet")
    lin_cells = mem["label"].to_numpy()  # block label = coarse lineage
    lin_mc = mo["label"].astype(str).to_numpy()
    labels = list(pd.Series(lin_cells).value_counts().index)  # biggest first
    cm = plt.get_cmap("tab20")
    cmap = {c: cm(i % 20) for i, c in enumerate(labels)}
    numbered = len(labels) > MAX_LEGEND
    out = {}
    sz = 6 + mo["size"].to_numpy() * 0.8

    def lab(c, j):
        return f"{j}. {_short(c)}" if numbered else _short(c)

    # grain re-analysis: normalize → HVG → PCA → kNN → UMAP + Leiden
    a = mc.copy()
    hv = vst_hvg(a.X, 2000)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    lognorm = a.copy()
    a = a[:, hv].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=int(min(30, a.n_obs - 1, a.n_vars - 1)), random_state=0)
    sc.pp.neighbors(a, n_neighbors=min(15, a.n_obs - 1), random_state=0)
    sc.tl.umap(a, random_state=0)
    sc.tl.leiden(
        a, resolution=LEIDEN_RES, random_state=0, key_added="cl", flavor="igraph", n_iterations=2, directed=False
    )
    cl = a.obs["cl"].astype(int).to_numpy()
    pd.DataFrame({"grain_id": mo.index, "viz_cluster": cl}).to_csv(d / "viz_clusters.tsv", sep="\t", index=False)
    cl_ids = sorted(np.unique(cl))
    cl_cmap = {c: plt.get_cmap("tab20b" if i % 2 else "tab20c")(i % 20) for i, c in enumerate(cl_ids)}
    xy = np.asarray(a.obsm["X_umap"])

    # ---- figure 1: two square UMAPs sharing the lineage legend
    umap = _cell_umap(mem, s["input"])
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.8), gridspec_kw={"wspace": 0.04})
    if "X_umap_mean" in mc.obsm:
        mxy = np.asarray(mc.obsm["X_umap_mean"])
        if umap is not None:  # raw cells: small translucent dots in their lineage colour
            for c in labels:
                i = lin_cells == c
                axes[0].scatter(
                    umap[i, 0],
                    umap[i, 1],
                    s=1.5,
                    color=cmap[c],
                    alpha=0.25,
                    linewidths=0,
                    zorder=1,
                )
        for c in labels:  # grains: opaque, bigger, on top
            i = lin_mc == c
            axes[0].scatter(
                mxy[i, 0],
                mxy[i, 1],
                s=sz[i],
                color=cmap[c],
                alpha=1.0,
                edgecolors="none",
                linewidths=0,
                zorder=3,
            )
        if numbered:
            _number(axes[0], labels, mxy, lin_mc)
    axes[0].set_title(
        f"Grains (opaque, dot size ∝ cell#) on raw cells (translucent)\n{len(mo)} grains",
        fontsize=F_TITLE,
    )
    for c in cl_ids:
        i = cl == c
        axes[1].scatter(
            xy[i, 0],
            xy[i, 1],
            s=sz[i],
            color=cl_cmap[c],
            alpha=1.0,
            edgecolors="none",
            linewidths=0,
            zorder=3,
        )
        axes[1].text(
            np.median(xy[i, 0]),
            np.median(xy[i, 1]),
            str(c),
            fontsize=F_NUM,
            fontweight="bold",
            ha="center",
            va="center",
            bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="none", alpha=0.8),
            zorder=6,
        )
    for c in labels:
        i = lin_mc == c
        if i.sum() >= 4:
            _islands(axes[1], xy[i], cmap[c])
    axes[1].set_title(
        f"Grain UMAP · Leiden clusters (res {LEIDEN_RES}, n = {len(cl_ids)})\nislands = coarse lineages",
        fontsize=F_TITLE,
    )
    axes[1].legend(
        handles=[
            Patch(
                facecolor=cmap[c],
                alpha=0.6,
                edgecolor=cmap[c],
                label=f"{lab(c, j)} ({(lin_mc == c).sum()})",
            )
            for j, c in enumerate(labels, 1)
            if (lin_mc == c).sum() > 0
        ],
        title="coarse lineage (grains)",
        title_fontsize=F_LEGEND,
        fontsize=F_LEGEND,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        handlelength=1.2,
        borderaxespad=0.2,
    )
    for ax in axes:
        ax.set_xticks([]), ax.set_yticks([])
        ax.set_box_aspect(1)
    out["umap"] = _png(fig)

    # ---- figure 2: split heatmap
    lognorm.obs["cl"] = pd.Categorical([str(c) for c in cl], categories=[str(c) for c in cl_ids])
    sc.tl.rank_genes_groups(lognorm, "cl", method="wilcoxon", n_genes=TOP_GENES)
    genes, gene_cl = [], []
    for c in cl_ids:
        for g in lognorm.uns["rank_genes_groups"]["names"][str(c)]:
            if g not in genes:
                genes.append(g), gene_cl.append(c)
    gene_cl = np.asarray(gene_cl)
    X = lognorm[:, genes].X.toarray()
    Z = np.clip((X - X.mean(0)) / (X.std(0) + 1e-9), -3, 3)
    means = np.vstack([Z[cl == c].mean(0) for c in cl_ids])
    split_order = [cl_ids[i] for i in _hier_order(means)]
    col_blocks = []
    for c in split_order:
        idx = np.flatnonzero(cl == c)
        col_blocks.append(idx[_hier_order(Z[idx])])
    row_idx = np.concatenate([np.flatnonzero(gene_cl == c) for c in split_order])
    n_rows = len(row_idx)
    widths = [max(len(b), 6) for b in col_blocks]
    row_h = 0.24
    fig = plt.figure(figsize=(15, row_h * n_rows + 1.6))
    gs = fig.add_gridspec(
        3,
        len(col_blocks),
        width_ratios=widths,
        height_ratios=[0.4, 0.3, row_h * n_rows],
        wspace=0.02,
        hspace=0.015,
        left=0.1,
        right=0.9,
        top=0.985,
        bottom=0.015,
    )
    lin_codes = np.array([labels.index(c) for c in lin_mc])
    strip_cmap = matplotlib.colors.ListedColormap([cmap[c] for c in labels])
    im = None
    for j, (c, idx) in enumerate(zip(split_order, col_blocks, strict=True)):
        ax_t = fig.add_subplot(gs[0, j])  # cluster id and size, sitting right on the strip
        ax_t.axis("off")
        ax_t.text(
            0.5,
            0.0 if j % 2 == 0 else 0.5,
            f"C{c}\n{len(idx)}",
            ha="center",
            va="bottom",
            fontsize=F_STRIP,
            color=cl_cmap[c],
            fontweight="bold",
            linespacing=0.9,
        )
        ax_s = fig.add_subplot(gs[1, j])
        ax_s.imshow(
            lin_codes[idx][None, :],
            aspect="auto",
            cmap=strip_cmap,
            vmin=-0.5,
            vmax=len(labels) - 0.5,
            interpolation="nearest",
        )
        ax_s.set_xticks([]), ax_s.set_yticks([])
        ax_h = fig.add_subplot(gs[2, j])
        im = ax_h.imshow(
            Z[idx][:, row_idx].T,
            aspect="auto",
            cmap="RdBu_r",
            vmin=-3,
            vmax=3,
            interpolation="nearest",
        )
        ax_h.set_xticks([])
        if j == 0:
            ax_h.set_yticks(np.arange(n_rows))
            ax_h.set_yticklabels([f"{genes[r]}  (C{gene_cl[r]})" for r in row_idx], fontsize=F_GENE)
        else:
            ax_h.set_yticks([])
        for spine in ax_h.spines.values():
            spine.set_linewidth(0.4)
    fig.legend(
        handles=[Patch(color=cmap[c], label=lab(c, j)) for j, c in enumerate(labels, 1) if (lin_mc == c).sum() > 0],
        loc="upper left",
        bbox_to_anchor=(0.903, 0.93),
        fontsize=F_LEGEND - 1,
        frameon=False,
        title="coarse lineage",
        title_fontsize=F_LEGEND - 1,
        handlelength=1.2,
    )
    cax = fig.add_axes([0.905, 0.03, 0.012, 0.22])
    fig.colorbar(im, cax=cax).ax.tick_params(labelsize=F_STRIP)
    fig.suptitle(
        f"Grain marker heatmap · rows = top-{TOP_GENES} Wilcoxon markers per grain cluster (de-duplicated, z-score ±3) · "
        f"columns split by cluster (C#, n), hierarchical order inside splits, splits ordered by clustering split means · strip = coarse lineage",
        fontsize=F_LEGEND,
        y=0.998,
    )
    out["heat"] = _png(fig)

    table = None
    if numbered:
        cells = pd.Series(lin_cells).value_counts().reindex(labels).fillna(0)
        mcn = pd.Series(lin_mc).value_counts().reindex(labels).fillna(0)
        outl = pd.Series(lin_cells[(mem["status"] == "outlier").to_numpy()]).value_counts().reindex(labels).fillna(0)
        rows = "".join(
            f"<tr><td>{j}</td><td><span style='display:inline-block;width:10px;height:10px;background:{matplotlib.colors.to_hex(cmap[c])};margin-right:6px'></span>{html_mod.escape(c)}</td>"
            f"<td>{int(cells[c])}</td><td>{int(mcn[c])}</td><td>{int(outl[c])}</td></tr>"
            for j, c in enumerate(labels, 1)
        )
        table = f"<table><tr><th>#</th><th>coarse lineage</th><th>cells</th><th>grains</th><th>outliers</th></tr>{rows}</table>"
    return s, out, table


def write_html(run_dirs, out_path):
    parts = []
    for rd in run_dirs:
        s, figs, table = figures(rd)
        q = " / ".join(f"{v:.0f}" for v in s["size_quantiles"].values())
        lines = (
            f"<li>{s['n_cells']} cells, {s['n_samples']} samples, {s['n_labels']} lineages, {s['n_blocks']} blocks</li>"
            f"<li>{s['n_grains']} grains; size quantiles 0/10/50/90/100%: {q} (γ = {s['params']['gamma']})</li>"
            f"<li>{s['n_outliers']} outliers ({s['outlier_rate']:.1%}); cell conservation {s['conservation']}</li>"
        )
        imgs = "".join(f'<figure><img src="data:image/png;base64,{b}"></figure>' for b in figs.values())
        tbl = f"<details open><summary>Lineage numbers</summary>{table}</details>" if table else ""
        parts.append(
            f"<section><h2>{Path(rd).name} <small>{s['input']}</small></h2><ul>{lines}</ul><div class='grid'>{imgs}</div>{tbl}</section>"
        )
    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>eca-grain results</title>
<style>body{{font:15px/1.6 -apple-system,"PingFang SC","Noto Sans CJK SC",sans-serif;color:#1f2933;max-width:1500px;margin:0 auto;padding:20px}}
h2{{border-top:1px solid #d9dee3;padding-top:10px;margin:24px 0 2px}} h2 small{{font-size:12px;color:#52606d;font-weight:400;display:block}}
.grid{{display:grid;grid-template-columns:1fr;gap:4px}} figure{{margin:0}}
img{{width:100%;height:auto;border:1px solid #e5e7eb;border-radius:4px;display:block}} ul{{font-size:14px;color:#52606d;margin:2px 0 6px}}
table{{border-collapse:collapse;font-size:12px;margin-top:8px}} td,th{{border:1px solid #d9dee3;padding:2px 8px;text-align:left}} th{{background:#f1f4f7}}
details summary{{cursor:pointer;color:#52606d;font-size:14px;margin-top:8px}}</style></head><body>
<h1>eca-grain results (final grains)</h1>
<p>Grains are built inside sample × coarse-lineage blocks. For display, the grain matrix is re-analysed (HVG → PCA → Leiden 1.0) to give grain clusters. Two figures per unit:
top left = grains (opaque, dot size ∝ cells) over the single cells (translucent) on the unit UMAP, both coloured by coarse lineage; top right = UMAP of the grains themselves (Leiden cluster numbers on top, translucent islands underneath = lineages), sharing the legend on the right.
Bottom = marker heatmap: rows = top-5 markers per grain cluster (de-duplicated); columns split by cluster, hierarchically ordered inside each split, splits ordered by clustering their mean profiles; strip = coarse lineage.</p>
{"".join(parts)}</body></html>"""
    Path(out_path).write_text(html)
    return out_path
