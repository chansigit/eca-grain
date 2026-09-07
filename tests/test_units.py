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
        x[:, :3] = 0  # constant genes stay unscaled, as in metacell_stats
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
