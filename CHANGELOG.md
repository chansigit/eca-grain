# Changelog

All notable changes to eca-grain. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

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
