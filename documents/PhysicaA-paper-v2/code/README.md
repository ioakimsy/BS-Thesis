# Learning dynamics as a heterogeneous cellular automaton: code and data

This repository contains the simulation and analysis code for

> C. I. T. Sy, R. M. L. Roxas-Villanueva, J. Y. Bantang, "Learning dynamics as a heterogeneous cellular automaton: comparing broadcast and nearest-neighbor interaction regimes in a classroom model", *Physica A* (under review).

It also contains the parameter table, including the random seed of every run, and the summary data behind every figure and table. The full simulation output is regenerated exactly from these seeds in a few minutes.

## Model

The model is a probabilistic cellular automaton on an L×L lattice with open boundaries:

- **States:** each agent is unlearned or learned, and learning is irreversible.
- **Learning rates:** each agent has λ = λ₀ ± δλ, with λ₀ = 0.5. The two classes are assigned independently with probability ½ and fixed for the run.
- **Update:** synchronous.
- **Broadcast regime (traditional instruction):** an unlearned agent learns with probability ρ₀λ.
- **Nearest-neighbor regime (peer instruction):** with k learned agents in its Moore neighborhood, it learns with probability 1 − (1 − ρ₀λ)ᵏ.
- **Run length:** every run continues until all agents have learned (absorption). The cap of 10⁶ steps is a safety limit, never reached.
- **Twin runs:** each broadcast run shares its random seed with a nearest-neighbor run, so the two start from the same initially learned agents and the same learning-rate assignment.

The full definition is in the header of `scripts/2DBPCA-IH-common-v2.jl` and in Section 2 of the paper.

## Requirements

- **Julia 1.11** (developed with 1.11.0). Newer Julia versions cannot load the pinned `Manifest.toml`. With [juliaup](https://github.com/JuliaLang/juliaup), run `juliaup add 1.11`, then `juliaup override set 1.11` in the repository folder, so that `julia` there means 1.11. Alternatively, start each command with `julia +1.11`.
- **Packages:** pinned in `Project.toml` and `Manifest.toml`. Install them once from the repository root:

  ```
  julia --project=. -e "using Pkg; Pkg.instantiate()"
  ```

- **OpenGL:** the 3D surface figures use GLMakie, which needs a graphics context. All other figures use CairoMakie.

## Running

Run every script **from the repository root**: paths are relative to it. Use several threads for the simulation.

| Step | Command | Time | Output |
|---|---|---|---|
| 1. Parameters and seeds | `julia --project=. scripts/2DBPCA-IH-parameter-generation-v2.jl` | seconds | `output/.../parameter_combinations.csv` (already included; rerunning keeps the existing seeds) |
| 2. Simulation | `julia -t auto --project=. scripts/2DBPCA-IH-simulation-v2.jl` | ~2.5 min on 16 threads | `output/.../raw/*.parquet` (~720 MB), `analysis/run_summary.csv` |
| 3. Validation (optional) | `julia -t auto --project=. scripts/2DBPCA-IH-validation-v2.jl` | a few min | prints PASS/FAIL for the correctness checks |
| 4. Return maps | `julia -t auto --project=. scripts/2DBPCA-IH-return-map-v2.jl` | ~1 min | `analysis/return-map-*.csv`, return-map figures |
| 5. Broadcast gain | `julia --project=. scripts/2DBPCA-IH-ti-gain-v2.jl` | ~1 min | `analysis/ti-gain-drift.csv`, Fig. 8 |
| 6. Figures | `julia -t auto --project=. scripts/2DBPCA-IH-plotting-v2.jl` | a few min | all remaining figures |
| 7. Quoted numbers | `julia --project=. scripts/2DBPCA-IH-paper-numbers-v2.jl` | seconds | prints Table 5 and the Section 3.3 values |

Here `output/...` stands for `output/2D-Binary-PCA-IH-v2`. Steps 4–6 read the raw simulation output, so step 2 must be run first. Step 7 needs only `run_summary.csv`, which is included.

`run_summary.csv` has one row per run: parameters, seed, and the times t_max, t₅₀, t₉₀, t₉₅ at which all, 50%, 90%, 95% of agents had learned. Each raw Parquet file stores, for every agent, the step at which it learned, so the full trajectory of any run can be reconstructed.

**Reproducibility check.** Rerunning step 2 from the included seeds reproduces the included `run_summary.csv` exactly for all 100,800 runs, except for the `runtime_s` column. Steps 4 and 5 reproduce the included `return-map-*.csv` and `ti-gain-drift.csv` files byte for byte.

## Paper figures and tables

Figure files are written to `output/2D-Binary-PCA-IH-v2/analysis/plots/`. Script file names follow the original submission's figure numbering, so they differ from the revised numbering below. Each figure is also saved without legends in `plots/clean/`; the published figures were labeled from these.

| Paper | File (under `analysis/plots/`) | Script |
|---|---|---|
| Fig. 1 (flowchart), Fig. 2 (seating arrangements) | drawn by hand | – |
| Fig. 3 | `learned-t/figure3-learned-t-v2.png`, `clean/figure3-quad-clean.png` | plotting |
| Fig. 4 | `class-split-learning-curves.png` | plotting |
| Fig. 5 | `figure4-snapshots-v2.png` | plotting |
| Fig. 6 | `figure5-tmax-t95-surface-v2.png` | plotting |
| Fig. 7 | `return-map/figure6-option1-student-level.png`, `clean/figure6-<L>-<ρ₀>-0.5-<δλ>-clean.png` | return map |
| Fig. 8 | `return-map/ti-gain-vs-exact.png` | broadcast gain |
| Fig. 9 | `return-map/figure6-option2-class-level.png` (requires the empirical data) | return map |
| Fig. 10 | `t95-vs-tmax.png` | plotting |
| Fig. 11 | `return-map/coefficients-vs-L.png` | return map |
| Fig. 12 | `seed-density-comparison.png` | plotting |
| Table 2 | printed by `t95_ratio_table()` | plotting |
| Table 3 (simulation rows) | `analysis/return-map-summary.csv` | return map |
| Table 4 | `analysis/ti-gain-drift.csv` | broadcast gain |
| Table 5, Section 3.3 values | printed | paper numbers |

## Empirical data

The classroom pre- and post-test scores used in Section 3.4 (Figs. 7 and 9, and the empirical rows of Table 3) are not included. They are available from the corresponding author on reasonable request. To include them, place the files `in_8_<k>.txt` and `out_8_<k>.txt` (k = 1, …, 14) in `data/empirical/`. Without them, the return-map script plots the simulations only and skips Fig. 9.

## Contents

```
Project.toml, Manifest.toml        Julia environment
scripts/
  2DBPCA-IH-common-v2.jl           model, parameters, file layout (included by every script)
  2DBPCA-IH-parameter-generation-v2.jl
  2DBPCA-IH-simulation-v2.jl
  2DBPCA-IH-validation-v2.jl
  2DBPCA-IH-return-map-v2.jl
  2DBPCA-IH-ti-gain-v2.jl
  2DBPCA-IH-plotting-v2.jl
  2DBPCA-IH-paper-numbers-v2.jl
output/2D-Binary-PCA-IH-v2/
  parameter_combinations.csv       5,040 parameter sets × 20 seeds
  analysis/*.csv                   run summary and analysis tables
data/empirical/                    (not included; see above)
```

## License and citation

[License to be chosen.] If you use this code, please cite the paper above [and the Zenodo DOI, to be added].
