"""End-to-end run on a small synthetic unit: every cell is accounted for, grains never cross a block,
summed counts match the members, and the report page carries the two template figures."""

import base64
import re

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from ecagrain.run import run


@pytest.fixture(scope="module")
def unit(tmp_path_factory):
    rng = np.random.default_rng(0)
    n_genes, per_block = 500, 160
    samples, labels = ["s1", "s2"], ["T cell", "Fibroblast"]
    obs, rows = [], []
    for s in samples:
        for j, lab in enumerate(labels):
            mu = rng.gamma(0.5, 2.0, n_genes)
            mu[j * 50 : (j + 1) * 50] *= 8  # 50 label-specific genes
            for _ in range(per_block):
                rows.append(rng.poisson(mu * rng.gamma(4, 0.25)))
                obs.append((s, lab, f"{lab} sub{rng.integers(2)}"))
    counts = sparse.csr_matrix(np.array(rows, dtype=np.float32))
    a = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(obs, columns=["eca_sample_id", "zmip_ann_coarse", "zmip_ann_fine"]),
        var=pd.DataFrame(index=[f"G{i:04d}" for i in range(n_genes)]),
    )
    a.obs_names = [f"cell{i}" for i in range(a.n_obs)]
    a.layers["counts"] = counts
    a.obsm["X_umap"] = rng.normal(size=(a.n_obs, 2))
    a.obsm["X_pca_harmony"] = rng.normal(size=(a.n_obs, 10))
    p = tmp_path_factory.mktemp("unit") / "final.h5ad"
    a.write_h5ad(p)
    return p, a


def test_run_conserves_cells_and_writes_report(unit, tmp_path):
    h5ad, a = unit
    out = tmp_path / "grain"
    s = run(h5ad, out, nrep=3)
    mc = ad.read_h5ad(out / "metacells.h5ad")
    mem = pd.read_parquet(out / "membership.parquet")

    assert s["conservation"] == "ok"
    assert mc.n_vars == a.n_vars and list(mc.var_names) == list(a.var_names)  # whole genome kept
    assert len(mem) == a.n_obs
    assert mc.obs["size"].sum() + s["n_outliers"] == a.n_obs
    assert set(mem.loc[mem.status == "member", "metacell_id"]) == set(mc.obs_names)
    assert mem.loc[mem.status == "outlier", "metacell_id"].isna().all()

    # a grain lives inside one sample × label block
    m = mem[mem.status == "member"]
    assert (m.groupby("metacell_id")["block"].nunique() == 1).all()
    assert (mc.obs.loc[m["metacell_id"], "block"].to_numpy() == m["block"].to_numpy()).all()

    # summed counts equal the members' counts
    g = mc.obs_names[0]
    cells = m.loc[m.metacell_id == g, "cell"]
    expect = np.asarray(a[cells].layers["counts"].sum(0)).ravel()
    assert np.allclose(mc[g].X.toarray().ravel(), expect)
    assert mc.obs.loc[g, "size"] == len(cells)

    assert set(mc.obs["mcRigor"]) <= {"trustworthy", "residual_dubious", "untested"}
    assert "X_umap_mean" in mc.obsm and "X_pca_harmony_mean" in mc.obsm
    assert (out / "threshold.tsv").exists() and (out / "summary.json").exists()

    html = (out / "report.html").read_text()
    pngs = re.findall(r'src="data:image/png;base64,([A-Za-z0-9+/=]+)"', html)
    assert len(pngs) == 2
    for p in pngs:
        assert base64.b64decode(p)[:8] == b"\x89PNG\r\n\x1a\n"


def test_missing_column_fails_loudly(unit, tmp_path):
    h5ad, _ = unit
    with pytest.raises(SystemExit, match="missing obs column"):
        run(h5ad, tmp_path / "x", label_col="nope")
