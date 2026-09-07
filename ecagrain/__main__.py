"""python -m ecagrain run <final.h5ad> <outdir> [--gamma 20 ...]
python -m ecagrain validate-rigor <supercell2.0/outputs/<tissue>> [--nrep 1]"""

from __future__ import annotations

import argparse
import os

# umap/pynndescent JIT-compile on import (15 s cold); a persistent numba cache cuts the report page to ~half.
os.environ.setdefault(
    "NUMBA_CACHE_DIR",
    os.path.join(os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache"), "ecagrain", "numba"),
)

from .run import COLS, DEFAULTS, run


def main():
    ap = argparse.ArgumentParser(prog="ecagrain")
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run", help="one-through grains for one released unit")
    r.add_argument("h5ad")
    r.add_argument("outdir")
    for k, v in {**DEFAULTS, **COLS}.items():
        r.add_argument(f"--{k.replace('_', '-')}", type=type(v), default=v)
    v = sp.add_parser(
        "validate-rigor",
        help="compare the mcRigor port with the R outputs of one June tissue",
    )
    v.add_argument("tissue_dir")
    v.add_argument("--nrep", type=int, default=1)
    f = sp.add_parser("figures", help="self-contained results page for one or more run dirs")
    f.add_argument("run_dirs", nargs="+")
    f.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "figures":
        from .figures import write_html

        print(write_html(a.run_dirs, a.out))
    elif a.cmd == "run":
        s = run(a.h5ad, a.outdir, **{k: getattr(a, k) for k in {**DEFAULTS, **COLS}})
        print(
            f"{s['n_cells']} cells → {s['n_grains']} grains, {s['n_outliers']} outliers, "
            f"{s['n_residual_dubious']} residual dubious, {s['elapsed_s']} s → {a.outdir}"
        )
    else:
        from .validate import validate

        validate(a.tissue_dir, a.nrep)


if __name__ == "__main__":
    main()
