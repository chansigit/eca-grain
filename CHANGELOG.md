# Changelog

All notable changes to eca-grain. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## 0.2.0 - 2026-09-07

### Changed
- HVG selection removes blocked genes before ranking, so every lineage gets the
  full 2000 HVGs (previously blocked genes took slots and were dropped
  afterwards, leaving about 1850). Grain partitions shift slightly.
- Per-lineage PCA uses the full LAPACK SVD instead of arpack: identical
  components, about six times faster on dense matrices.
- mcRigor T statistic is evaluated through the n×n Gram matrix instead of the
  p×p correlation matrix: same value to 1e-15, tens of times faster.
- `elapsed_s` in `summary.json` now covers the report page too
  (`elapsed_pipeline_s` keeps the old measurement).

### Fixed
- Grains too small to test can no longer be labelled `residual_dubious`; a
  unit without any testable grain no longer crashes.
- Sample ids containing `|` are rejected with a clear message instead of
  corrupting block names.
- Ribosomal regex (used only when the input lacks OSP's `ribo` column) no
  longer matches RPS6K kinases or pseudogenes and now includes RPLP0/1/2 and
  RPSA.
- `membership.parquet` `level` is an integer column; the heatmap legend lists
  only lineages that have grains.

## 0.1.0 - 2026-09-07

### Added
- One-pass grain construction from a released ECA-RSI unit: SuperCell-style
  graph partition inside sample × coarse-lineage blocks, MetaCells 2 gaps
  outlier rule with depth, adaptive-fold and multi-gene guards, a numpy port of
  mcRigor DETECT (matches the R package to 1e-9 on the same HVGs), and a
  split-only recheck of dubious grains at a γ/2 floor.
- Cell conservation is enforced globally and per block: Σ size + outliers equals
  the input cell count, and every cell has a row in `membership.parquet`.
- Per-unit `report.html` with two figures: grains over raw cells on the unit
  UMAP, a grain-level UMAP with Leiden clusters over concave-hull lineage
  islands, and a marker heatmap split by grain cluster.
- Gene blocking for HVG selection only (cell cycle, dissociation stress, heat
  shock, mitochondrial, ribosomal, Malat1); haemoglobin genes are kept.
- `eca-grain run | figures | validate-rigor` command line.
