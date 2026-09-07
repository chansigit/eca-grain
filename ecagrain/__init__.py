"""ecagrain — GRAIN: Guarded Rigor-Audited In-lineage Neighborhoods.

Aggregate a released ECA-RSI unit into grains (about γ expression-similar cells each, counts summed) in one pass:
build (SuperCell algorithm) → outlier (MetaCells 2 gaps rule) → diagnose (mcRigor port) → recheck (split in place).
Design record: docs/design.md."""

__version__ = "0.4.0"
