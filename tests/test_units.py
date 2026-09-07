"""Unit checks: HVG selection fills n_top after exclusion, the Gram-form T statistic equals the p×p form,
and the ribosomal regex hits ribosomal proteins only."""

import re

import numpy as np
from scipy import sparse

from ecagrain.genes import FAMILY_REGEX
from ecagrain.hvg import vst_hvg
from ecagrain.rigor import scale_cols, t_stat


def _counts(n=300, p=400, seed=0):
    rng = np.random.default_rng(seed)
    mu = rng.gamma(0.5, 2.0, p)
    return sparse.csr_matrix(rng.poisson(mu, size=(n, p)).astype(np.float32)), rng


def test_hvg_fills_n_top_after_exclusion():
    X, rng = _counts()
    exclude = np.zeros(X.shape[1], bool)
    exclude[:100] = True
    for batch in (None, rng.integers(3, size=X.shape[0])):
        hv = vst_hvg(X, 200, batch=batch, exclude=exclude)
        assert hv.sum() == 200
        assert not (hv & exclude).any()
    assert vst_hvg(X, 350, exclude=exclude).sum() == 300  # only 300 candidates left


def test_t_stat_gram_equals_pxp():
    rng = np.random.default_rng(1)
    for n, p in ((5, 40), (20, 800), (37, 1500)):
        x = rng.normal(size=(n, p))
        x[:, :3] = 0  # constant genes stay unscaled, as in grain_stats
        dat = scale_cols(x)
        c = dat - dat.mean(0)
        cov = c.T @ c / (n - 1)
        cov[np.diag_indices(p)] -= 1.0
        ref = np.linalg.norm(cov) / np.sqrt(p * (p - 0.5))
        assert abs(t_stat(dat) - ref) < 1e-12


def test_ribo_regex():
    rx = re.compile(FAMILY_REGEX["ribo"], re.I)
    for g in ("RPS3A", "RPL7", "Rpl36al", "RPLP0", "RPSA", "RPS4Y1", "RPL22L1", "Rps27l"):
        assert rx.match(g), g
    for g in ("RPS6KA1", "RPS6KB2", "Rps23rg1", "Rpl15-ps6", "RPL34-DT", "RPS10-NUDT3", "RPN1"):
        assert not rx.match(g), g


def test_gap_flags_top_and_bottom_tails():
    from ecagrain.outlier import gap_flags

    # gene 0: one cell far above the rest; gene 1: one cell far below; gene 2: a 2-cell top tail
    V = np.array([[0, 4, 0], [0, 4, 0], [0, 4, 0], [0, 4, 0], [0, 4, 5], [5, 0, 5]], float)
    flags = gap_flags(V, fold=3.0, allowed=1)
    expect = np.zeros_like(flags)
    expect[5, 0] = expect[5, 1] = True  # the 2-cell tail of gene 2 exceeds allowed=1
    assert (flags == expect).all()
    flags2 = gap_flags(V, fold=3.0, allowed=2)
    assert flags2[:, 2].tolist() == [False, False, False, False, True, True]
    assert (flags2[:, :2] == expect[:, :2]).all()
    assert not gap_flags(V, fold=6.0, allowed=2).any()  # gaps of 4 and 5 are below a fold of 6


def test_enforce_cap_bounds_every_group():
    from ecagrain.build import partition_block

    rng = np.random.default_rng(0)
    coords = np.vstack([rng.normal(c, 1.0, size=(200, 5)) for c in (0, 8, 16)])
    memb, g = partition_block(coords, gamma=20, k=15, cap=30)
    sizes = np.bincount(memb)
    assert sizes.max() <= 30 and sizes.min() >= 1
    assert memb.min() == 0 and memb.max() == len(sizes) - 1  # contiguous ids
    assert g.vcount() == len(coords)


def test_flag_dubious_never_flags_untested():
    import pandas as pd

    from ecagrain.rigor import flag_dubious

    tested = pd.DataFrame(
        {
            "size": [5, 10, 20],
            "T_org": [1.0, 1.0, 1.0],
            "T_colperm": [0.5, 0.5, 0.5],
            "T_rowperm1": [0.5] * 3,
            "T_bothperm1": [1.0] * 3,
        },
        index=[0, 1, 2],
    )
    tested["TT_div"] = tested["T_org"] / tested["T_colperm"]  # 2.0, far above the null ratio 0.5
    untested = pd.DataFrame({"size": [3], "TT_div": [1.0]}, index=[3])
    st = pd.concat([tested, untested])
    dub, thre = flag_dubious(st, 0.05)
    assert dub.loc[[0, 1, 2]].all() and not dub.loc[3]
    assert (thre["thre"] < 1.0).all()  # the guard, not the threshold, protects the untested grain
    dub0, thre0 = flag_dubious(untested, 0.05)
    assert not dub0.any() and thre0.empty


def test_pca_scores_matches_scanpy_path():
    import anndata as ad
    import scanpy as sc

    from ecagrain.build import lognorm, pca_scores

    X, _ = _counts(n=400, p=300, seed=2)
    a = ad.AnnData(X=X.copy())
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    L = lognorm(X)
    assert np.allclose(L.toarray(), a.X.toarray(), atol=1e-5)
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=10, svd_solver="full")
    ref = np.asarray(a.obsm["X_pca"])
    for chunk in (400, 128):  # single chunk and several chunks
        new = pca_scores(L, 10, chunk=chunk)
        r = [abs(np.corrcoef(ref[:, i], new[:, i])[0, 1]) for i in range(10)]
        assert min(r) > 0.9999
        assert np.allclose(new.var(0), ref.var(0), rtol=1e-4)
