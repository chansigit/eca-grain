<p align="center">
  <img src="assets/logo.svg" alt="eca-grain logo: cells grouped into grains inside a lineage island" width="176" height="176">
</p>

<h1 align="center">eca-grain: GRAIN</h1>

<p align="center">
  <strong>Guarded Rigor-Audited In-lineage Neighborhoods.</strong><br>
  One-pass, fully audited grains (metacell-style aggregates) from a released ECA-RSI unit.
</p>

<p align="center">
  <a href="pyproject.toml"><img src="https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat&amp;logo=python&amp;logoColor=white" alt="Python 3.10 or newer"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-173B49?style=flat" alt="MIT license"></a>
  <a href="https://github.com/chansigit/eca-rsi"><img src="https://img.shields.io/badge/Ecosystem-ECA--RSI-258B81?style=flat" alt="Part of the ECA-RSI ecosystem"></a>
</p>

<p align="center">
  <a href="#why-grains">Why grains</a> ·
  <a href="#how-it-works">How it works</a> ·
  <a href="#get-started">Get started</a> ·
  <a href="#read-your-results">Results</a> ·
  <a href="#documentation">Documentation</a>
</p>

eca-grain takes the `release/final.h5ad` of an [ECA-RSI](https://github.com/chansigit/eca-rsi)
unit and groups its cells into **grains**: about γ = 20 expression-similar cells
each, drawn from one sample and one coarse lineage, with counts summed. It runs
once with fixed parameters, drops outlier cells the way MetaCells 2 does,
audits every grain with a Python port of mcRigor, splits the doubtful ones in
place, and keeps every cell on the ledger. The output is meant for downstream
work such as gene-regulatory-network inference, where grain counts should stay
roughly proportional to cell counts rather than being balanced across types.

## Why grains

Single cells are sparse and noisy; aggregating neighbours gives denser
profiles without losing population structure. Existing tools either need
several rounds of human tuning (MetaCells 2) or never drop a cell (SuperCell,
SEACells). eca-grain fixes the workflow instead: one pass, fixed parameters,
outliers removed but recorded, and every grain checked for hidden mixing.

Three guarantees hold on every run:

- **Guarded grouping.** A grain never crosses a sample or a coarse lineage.
  Coordinates are computed per lineage (samples pooled, sample-aware HVGs),
  grouping happens inside each sample × lineage block.
- **Rigor audit.** Each grain of at least 5 cells gets an mcRigor
  heterogeneity score against a unit-level threshold curve. Dubious grains are
  split once at a γ/2 floor; whatever still fails is delivered with the
  `residual_dubious` flag, never silently.
- **Cell conservation.** Σ grain sizes + outliers equals the input cell count,
  globally and per block, or the run fails.

## How it works

| Stage | What happens | Module |
| --- | --- | --- |
| Build | Per coarse lineage: own vst HVG selection (per-sample ranking, blocked genes removed before ranking) → PCA of the scaled matrix computed from a chunked covariance (no dense matrix, memory independent of cell count). Per sample × lineage block: kNN graph → walktrap → cut into round(n/γ) communities; communities above 1.5γ are re-cut on their own subgraph. Blocks below γ form one grain. Lineages under 100 cells fall back to the unit's `X_pca_harmony`. | `build.py`, `hvg.py`, `genes.py` |
| Outlier | MetaCells 2 "gaps" deviant rule with three guards: depth rescaling, an adaptive fold from the block's gap distribution, and at least three deviant genes. At most 25% of a block can be dropped. Outliers keep their row in the membership table. | `outlier.py` |
| Diagnose | mcRigor DETECT in numpy (matches the R package to 1e-9 on the same HVGs): T = ‖corr − I‖_F / √(p(p − 0.5)), TT_div = T / T(column-permuted), null from row permutation, per-size quantile + lowess threshold, Nrep = 20. Grains under 5 cells are `untested`. | `rigor.py` |
| Recheck | Dubious grains are split in two on their own subgraph (floor γ/2, one level) and re-tested against the saved threshold. Failures stay as `residual_dubious`. | `run.py` |
| Deliver | Conservation check, summed counts, mean embeddings, membership ledger, summary and the report page. | `run.py`, `figures.py` |

The output keeps every gene of the input; blocked genes affect HVG selection
only and are never removed from the output matrix: Tirosh cell-cycle genes, the van den Brink dissociation-stress
panel (taken from [OSP](https://github.com/chansigit/osp)), heat-shock genes
(`HSPA/B/D/E/H`, `HSP90`, `DNAJA/B`), mitochondrial and ribosomal genes, and
`MALAT1`. Haemoglobin genes are deliberately kept so erythrocyte datasets
group normally.

## Get started

### 1. Install

```bash
python -m pip install "eca-grain @ git+https://github.com/chansigit/eca-grain.git"
```

Python 3.10 or newer. The dependency on `osp-sc` provides the dissociation
gene panel; no agent runtime is used.

### 2. Run one unit

```bash
eca-grain run <unit>/release/final.h5ad <outdir>
# equivalently: python -m ecagrain run ...
```

The input is an ECA-RSI release: a `counts` layer, `obs["eca_sample_id"]`,
`obs["zmip_ann_coarse"]` (block boundary) and `obs["zmip_ann_fine"]` (recorded
per grain as majority label and purity, never used for grouping), optionally
`obsm["X_pca_harmony"]` and `obsm["X_umap"]`. Older releases without
`eca_sample_id` can pass `--sample-col <column>`. Every entry of the fixed
parameter set can be overridden on the command line (`--gamma`, `--k`,
`--n-hvg`, `--nrep`, …) but the pipeline is designed to run unchanged across
datasets.

### 3. Combine pages

```bash
eca-grain figures <outdir1> <outdir2> ... --out results.html
```

Each run already writes its own `report.html`; this stitches several units
into one page. The report step imports umap-learn, whose numba compilation
costs about 15 s cold; the command line sets `NUMBA_CACHE_DIR` (under
`XDG_CACHE_HOME` or `~/.cache/ecagrain`) when it is unset so this happens once
per machine.

### 4. Run many units

Released units live at `<dataset>/eca-pp/<Tissue>/rsi/units/<unit>/release/final.h5ad`.
The convention in this project is to write grains next to `rsi/`, in
`<dataset>/eca-pp/<Tissue>/grain/`, one directory per unit:

```bash
IN=.../eca-pp/Liver/rsi/units/liver/release/final.h5ad
eca-grain run "$IN" "${IN%/rsi/units/*}/grain"
```

Units are independent, so a scheduler array is the whole batch driver. On Slurm,
one task per line of a manifest:

```bash
ls */eca-pp/*/rsi/units/*/release/final.h5ad > units.txt
sbatch --array=1-$(wc -l < units.txt) --cpus-per-task=4 --mem=24G --time=02:00:00 \
       --wrap 'IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" units.txt); eca-grain run "$IN" "${IN%/rsi/units/*}/grain"'
```

Budget about 4 CPUs and 24 GB per task. In a 102-unit, 1.16-million-cell run the
largest unit (69k cells) took 90 s and peaked near 12 GB; the median unit took
under a minute. Set `NUMBA_CACHE_DIR` and `MPLCONFIGDIR` to a shared writable
path so the tasks do not each recompile umap or rebuild the font cache.

## Read your results

| File | Content |
| --- | --- |
| `grains.h5ad` | One row per grain. `X` = summed counts. `obs`: `size`, `sample`, `label` (lineage), `block`, `audit_majority` / `audit_purity` (fine label), `level`, `build_id` / `parent_build_id`, `mcRigor` status (`trustworthy`, `residual_dubious`, `untested`), `TT_div`, `n_test_genes`, `gamma`. `obsm`: `X_umap_mean`, `X_pca_harmony_mean`. `uns["ecagrain"]`: version, parameters, columns, input path. |
| `membership.parquet` | One row per input cell: `cell`, `sample`, `label`, `audit`, `block`, `grain_id` (empty for outliers), `status` (`member` / `outlier`), `n_flag_genes`, `build_id`, `level`, `mcRigor`, `umap_1` / `umap_2` (the unit's cell UMAP, so the report page needs no input file). |
| `threshold.tsv` | The unit-level mcRigor threshold curve by grain size. |
| `summary.json`, `report.md` | Counts by stage and status, size quantiles, audit purity, timing (`elapsed_s` includes the report page, `elapsed_pipeline_s` does not), and a per-block table with the outlier fold used. |
| `report.html` | Self-contained page: two square UMAPs sharing a lineage legend (grains over translucent raw cells; a grain-level UMAP with Leiden clusters over concave-hull lineage islands) and a marker heatmap split by grain cluster. |
| `viz_clusters.tsv` | The Leiden cluster of each grain on the report page. |

Grain counts follow cell counts: a lineage with ten times more cells gets
about ten times more grains. Intermediate states (dubious, split, threshold)
are recorded in the files but do not appear on the page.

### Using grains downstream

`grains.h5ad` holds summed raw counts over the input's complete gene set, so it
drops into any workflow that expects a counts matrix:

```python
import anndata as ad, scanpy as sc

g = ad.read_h5ad(".../grain/grains.h5ad")     # grains × all genes, summed counts
g = g[g.obs["mcRigor"] == "trustworthy"]      # optional; see the note below
sc.pp.normalize_total(g, target_sum=1e4)      # grains differ in size, so normalise
sc.pp.log1p(g)
```

Three things worth knowing before filtering or weighting:

- **`size` is the number of cells**, so it is the natural weight for a
  regression and the reason to normalise before comparing grains.
- **`mcRigor` is advice, not a verdict.** `trustworthy` passed the heterogeneity
  test, `residual_dubious` failed it and could not be split at the γ/2 floor,
  `untested` had fewer than 5 cells. Dropping `residual_dubious` costs roughly a
  fifth of the cells and biases against small blocks; keeping everything is the
  default for network inference.
- **`membership.parquet` traces every cell**, including the outliers that no
  grain contains, so any grain-level result can be pushed back to cells.

## Fixed parameters

γ 20 · size cap 1.5γ · split floor γ/2 · minimum tested size 5 · 2000 HVGs ·
up to 30 PCs · k 15 · lineage pool below 100 cells uses the global embedding ·
mcRigor Nrep 20, cutoff 0.05, gene filter 0.1 · seed 0.

## Known limitations

- Smart-seq2 units (Tabula Muris FACS) lose about 12% of cells as outliers;
  UMI units lose 0.1% to 2.4%.
- Roughly a fifth of the cells of a unit sit in `residual_dubious` grains,
  mostly blocks too small to split at the γ/2 floor. They are delivered with
  the flag.
- mcRigor's dubious calls track grain size more than fine-label mixing.
- Figure text is English only (no CJK font is assumed on compute nodes).
- The report page has a fixed cost of about 20 s per unit (scanpy and umap-learn
  imports plus UMAP itself); the pipeline proper takes about 1 s per 1000 cells.

## Command reference

| Command | What it does |
| --- | --- |
| `eca-grain run <final.h5ad> <outdir>` | The whole pipeline for one unit. Flags: `--sample-col` / `--label-col` / `--audit-col` / `--counts-layer` / `--embed-key` for input columns, and one flag per fixed parameter (`--gamma`, `--k`, `--max-size-factor`, `--n-hvg`, `--n-pcs`, `--min-pool`, `--min-test-size`, `--nrep`, `--test-cutoff`, `--gene-filter`, `--seed`). |
| `eca-grain figures <run_dir>... --out page.html` | Rebuild the report page for one run, or stitch several runs into one page. Reads only the run directories. |
| `eca-grain validate-rigor <tissue_dir>` | Compare the mcRigor port against R outputs of the same tissue (`input/counts.mtx`, `supercell/membership.tsv`, `mcrigor/*.tsv`). |

`python -m ecagrain ...` is equivalent to `eca-grain ...` and needs no
installation when the repository is on `PYTHONPATH`. Each stage prints a line
with its elapsed time, so a stalled batch task shows where it stopped.

## Troubleshooting

| Symptom | Cause and fix |
| --- | --- |
| `missing obs column 'eca_sample_id'` | An older release that predates the column. Pass another per-sample column, for example `--sample-col channel` or `--sample-col project`. |
| `sample ids must not contain '\|'` | Block names are `sample\|lineage`. Choose a different sample column. |
| `label 'X' has N cells < min_pool and no X_pca_harmony` | A tiny lineage needs the unit's global embedding as a fallback. Re-release with `obsm["X_pca_harmony"]`, or raise `--min-pool` so the lineage computes its own coordinates. |
| `conservation failed …` | A genuine bug: the run aborts instead of writing a partial result. Report the unit and the message. |
| The report page has no cell UMAP | The input had no `obsm["X_umap"]`, or the run predates 0.3.0 and its input has moved. Rerun `eca-grain run` to store the coordinates in `membership.parquet`. |
| Everything is slow, especially PCA | A numpy built without BLAS makes matrix products two orders of magnitude slower. eca-grain routes its own heavy products through `scipy.linalg.blas`, but scanpy and scikit-learn calls in the same environment stay slow. Check with `numpy.show_config()`. |

## Methods and credits

- Grouping follows the SuperCell algorithm (Bilous et al., BMC Bioinformatics 2022): kNN graph, walktrap, cut.
- The outlier rule adapts `find_deviant_cells` (gaps policy) from MetaCells 2 (Ben-Kiki et al., Genome Biology 2022).
- `rigor.py` re-implements mcRigor DETECT (Liu & Li, Nature Communications 2025; R package by Pan Liu, MIT licence) from its published source; it is validated against the R outputs with `eca-grain validate-rigor`.
- The dissociation-stress panel is van den Brink et al., Nature Methods 2017, as curated in OSP.

## Documentation

- [`docs/design.md`](docs/design.md): the dated design record (Chinese), including the decisions on coarse-lineage blocks, gene blocking and the report-page template.
- [`docs/design-2026-09-06.html`](docs/design-2026-09-06.html): the earlier plain-language walkthrough of the workflow.
- [`CHANGELOG.md`](CHANGELOG.md).

## Position in the ecosystem

[ECA-PP](https://github.com/chansigit/eca-pp) prepares inputs; [ECA-RSI](https://github.com/chansigit/eca-rsi)
drives [OSP](https://github.com/chansigit/osp), [MSP](https://github.com/chansigit/msp) and
[ZMIP](https://github.com/chansigit/zmip) to a released unit; eca-grain runs after release and does not modify it.
