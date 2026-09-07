# How well does mcRigor detect heterogeneous grains at γ = 20?

Measured 2026-09-07, because the `residual_dubious` flag showed almost no correlation with fine-label purity on real
grains, which raised the question of whether the test carries any signal at our granularity.

## Design

Grains cannot cross a coarse lineage by construction, so the realistic failure mode is a grain that mixes two fine
states inside one sample × coarse lineage. Synthetic grains of known composition were drawn from two released liver
units, and separated by TT_div alone (AUC, threshold-free, 200 replicates per condition):

- `pure` — all cells from one fine label,
- `mix_fine_XX` — XX% of the cells replaced by a second fine label of the same coarse lineage and sample,
- `mix_lineage` — 50/50 from two different coarse lineages of one sample (a positive control that the pipeline
  cannot actually produce).

Script: `eca-grain-jobs/rigor_power/power.py` (outside the repository; it reads released units).

## Result

AUC of TT_div separating each condition from `pure`:

| grain size | mix_fine_10 | mix_fine_25 | mix_fine_50 | mix_lineage |
| ---: | ---: | ---: | ---: | ---: |
| 20 (our γ) | 0.62 / 0.65 | 0.73 / 0.71 | 0.77 / 0.77 | 0.97 / 0.98 |
| 50 | 0.63 / 0.67 | 0.73 / 0.74 | 0.77 / 0.77 | 0.94 / 1.00 |
| 100 | 0.66 / 0.63 | 0.76 / 0.71 | 0.79 / 0.79 | 0.94 / 1.00 |
| 200 | 0.66 / 0.86 | 0.74 / 0.97 | 0.76 / 1.00 | 0.97 / — |

Two values per cell: mca3.0 Liver (34,708 cells) / mca1.1 Liver (4,725 cells). The mca1.1 column at size 200 rests on
few independent pools in a small unit and should not be read as a size effect.

## Reading

- The statistic is not noise at γ = 20. It separates grossly mixed grains almost perfectly and half-and-half fine-state
  mixtures with AUC ≈ 0.77.
- Power is flat in grain size for the failure mode that matters. Testing at a coarser granularity would not help, so
  the diagnosis stays where the grains are built.
- Power falls off with the amount of contamination: at 10% foreign cells, AUC ≈ 0.65. Real grains sit at that end of
  the range, which explains why the per-grain flag correlates so weakly with fine-label purity. The flag is a real but
  weak per-grain signal, and a meaningful one in aggregate.
- Consequence for users: do not treat `residual_dubious` as a verdict on a single grain, and do not filter on it by
  default. Use it to compare units or blocks, or as a covariate.
