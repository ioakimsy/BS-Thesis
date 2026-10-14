# PHYSA-263545 — Synthesis of Reviewer Comments

Source: `PHYSA-263545-reviews.pdf` (three reviewers, Sept 2026).

| Reviewer | Date | Recommendation | Scope of comments |
|---|---|---|---|
| R1 | 08/09/2026 | **Major revision** | Longest; methodological (controls, seed density, stopping rule, statistics, reproducibility) plus interpretation |
| R2 | 22/09/2026 | Revision (not stated as major/minor; tone is "clarify and present more cautiously") | Mostly framing and interpretation; overlaps heavily with R1 |
| R3 | 23/09/2026 | Revision before further consideration | Technical consistency of the model; includes an analytical benchmark for the broadcast regime |

All three agree the topic fits Physica A and that the manuscript is well organized. R2 explicitly calls the multi-stage relaxation in the heterogeneous broadcast regime "an interesting result." The criticism is about **whether the conclusions are controlled and robust**, and whether the **claims go beyond what was simulated**.

---

## Consensus at a glance

| # | Theme | R1 | R2 | R3 | Kind of work needed |
|---|---|:-:|:-:|:-:|---|
| 1 | Stopping criterion: sensitivity and how terminated runs are handled | ✔ | ✔ | ✔ | **New simulations** + text |
| 2 | "Optimal" mixed TI+PI strategy is not demonstrated | ✔ | ✔ | ✔ | Text |
| 3 | Empirical / Hake-gain agreement is overstated | ✔ | ✔ | ✔ | Text (maybe a small analysis) |
| 4 | Return-map Fig. 6 / Table 3 conflates system size with other parameters | ✔ | | ✔ | **New simulations** or reframing |
| 5 | Missing implementation details (update scheme, λ-class assignment, seats, code) | ✔ | | ✔ | Text + code deposit |
| 6 | Flowchart (Fig. 1) is inconsistent with the text | ✔ | | ✔ | Figure fix |
| 7 | "Percolation-driven" lacks quantitative support | ✔ | ✔ | | Terminology change or **new analysis** |
| 8 | Framing: model dynamics vs. real classrooms (sizes, "small class", "safe choice") | ✔ | ✔ | | Text |
| 9 | TI vs. PI comparison is not controlled (initial state + update rule differ) | ✔ | | | Discussion, possibly **new simulations** |
| 10 | Size scaling under a fixed 4 seeds means seed density falls as L grows | ✔ | | | **New simulations** |
| 11 | Broadcast return map: a constant gain (Eq. 6) is inconsistent with heterogeneous Eq. 3 | | | ✔ | Analysis + text |
| 12 | Quadratic PI return map is a fit, not an evolution equation | | | ✔ | Text |
| 13 | Uncertainty reporting (min–max bands; correlated return-map points) | ✔ | | | Re-analysis of existing data |
| 14 | Abstract claims a mean-learning-rate dependence that was not studied | ✔ | | | Text (or new sweep) |
| 15 | Research gap in the Introduction is not stated sharply | | ✔ | | Text |

---

## 1. Stopping criterion and relaxation times — *all three reviewers*

**The concern:** A run stops after a window of consecutive steps with no transition, whether or not it is fully absorbed. The reported t_max values are therefore lower bounds. This matters most in exactly the regime that yields the headline multi-stage relaxation result: the strongly heterogeneous broadcast case. R1 cites **17.4 %** of runs terminated by the stopping rule rather than absorption, and R2 names **δλ = 0.49**.

**What they ask for:**
- A **sensitivity analysis over the stopping-window length**, at least for representative cases in the affected range (R1, R2, R3).
- A clear statement of **how terminated runs enter the averaged learning curves and t_max** (R3), and a separation of absorbed from prematurely terminated runs in the relaxation-time analysis (R1).
- R3 also asks that the **four-stage behavior be checked against the analytical broadcast curve** (see §11) and against the **log-time axis**. The aim is to separate genuine two-timescale structure from artifacts of plotting or of early termination.

**Related terminology point (R1):** The Fig. 1 flowchart says "steady state achieved," but the Methods say the stopping rule does not guarantee a steady state. These need to agree.

> **Note for response planning (internal, not from reviewers):** As recorded in the repo's `CLAUDE.md`, the manuscript's stated "100 consecutive time steps" does not match the data. The δλ ≤ 0.4 sweep used a fixed **20-step** window. The δλ ∈ {0.45, 0.49} runs used **20 % of L²**. This has to be resolved before any sensitivity analysis is written up, because the reviewers' most-cited numbers (the 17.4 % at δλ = 0.49) come from the subset that used the scaled window.

## 2. "Optimal" combined strategy is not demonstrated — *all three*

None of the simulations combines broadcast and nearest-neighbor dynamics, so the conclusion that mixing them is "optimal" is unsupported. Every reviewer asks for it to be recast as a **hypothesis or future direction**. R1 additionally objects to calling nearest-neighbor interaction a generally **"safe choice"**, given that the model is binary-state with irreversible learning.

## 3. Empirical comparison and Hake gain — *all three*

- **Moderate "qualitative agreement" or "validation" language** (R1, R2, R3). The points they cite:
  - The empirical data have large scatter (**R² ≈ 0.63**).
  - The **sign of the fitted quadratic coefficient differs** between the empirical data and the nearest-neighbor simulations. The manuscript already acknowledges this.
  - The **seating-arrangement result differs from previous empirical findings** (R2).
- **Say explicitly which features are reproduced and which are not** (R1, R2).
- **Justify the mapping from simulation step to assessment interval** (R1). Hake gain is defined between two assessments (pre/post), whereas the return map compares consecutive time steps. What one step represents in classroom time is not established.

## 4. Return-map size comparison is confounded — *R1, R3*

The "small" and "large" system examples in Fig. 6 / Table 3 differ in more than L. R3 lists the differences:
- **ρ₀ goes from 0.9 to 0.1.**
- The learning-rate distribution goes from **homogeneous to heterogeneous.**

The coefficient differences therefore cannot be attributed to system size. The options are:
- (a) **redo the comparison at fixed ρ₀, δλ, and initialization** (both reviewers prefer this), or
- (b) state that Fig. 6 only illustrates representative parameter combinations, and drop the size attribution (R3 accepts this as an alternative).

## 5. Implementation details and reproducibility — *R1, R3*

The reviewers ask for the following to be stated explicitly:
- Whether the update is **synchronous or asynchronous** (R1, R3).
- How the **two-level λ distribution** is built (R1, R3):
  - Is the split exactly half/half or probabilistic?
  - Is the assignment spatially random?
  - Does it stay fixed within a run?
  - Is it regenerated for each trial?
- The **precise seat coordinates or rules** for each seating arrangement, since boundary effects may matter (R1).
- A **public, permanent code and data repository** in place of "available upon request" (R1).

## 6. Flowchart (Fig. 1) inconsistencies — *R1, R3*

- It shows **zero initially learned agents** for PI, but the text says four. Under Eq. 4 an all-unlearned PI lattice would never change (R3). The flowchart should show each regime's actual initial condition.
- It says **"steady state achieved"** (R1; see §1).

## 7. "Percolation-driven" terminology — *R1, R2*

The manuscript presents no cluster-size distributions, spanning probabilities, critical thresholds, or finite-size scaling. What the results show is local stochastic spreading, not a percolation transition. The authors should either **add quantitative percolation analysis** or switch to **"percolation-like spreading" / "local propagation dynamics."** Both reviewers offer the same alternative wording.

## 8. Framing relative to real classrooms — *R1, R2*

- The simulated populations (smallest > 1000 agents) are far larger than real classes. **"Small class"** language is therefore misleading. Distinguish **lattice size vs. population size** (R1).
- TI and PI are represented by **simplified topologies**, not full descriptions of either pedagogy. Conclusions should be framed in terms of **model dynamics rather than educational effectiveness** (R2).

## 9. Fairness of the TI vs. PI comparison — *R1 only, but substantial*

The two regimes differ in more than topology:
- **Initial state:** TI starts with 0 learned agents; PI starts with 4 seeds.
- **Exposure mechanism:** in TI, a global source acts independently on each agent. In PI, the learning probability compounds over several learned neighbors.

R1 asks:
- What does "controlled comparison" mean here?
- Either run simulations with **comparable initial learned fractions**, or **discuss explicitly** how the different initializations and update rules shape the conclusions.

## 10. Seed density in PI size scaling — *R1 only, but substantial*

With a fixed 4 seeds, a larger L means a lower seed fraction and a longer mean distance to the nearest seed. That alone could produce the observed growth of t_max with N. R1 asks for **at least part of the size-scaling analysis to be repeated at a fixed initial learned fraction**. The goal is to see whether the scaling is intrinsic to nearest-neighbor dynamics or comes from seed dilution.

## 11. Broadcast return map vs. heterogeneous dynamics — *R3*

Eq. 6 implies a **constant normalized gain**. With two λ-classes, however, the fast learners are depleted first, so the remaining unlearned pool shifts toward slow learners and the aggregate gain should change over time. R3 supplies the benchmark:

$$F(t) = 1 - w(1-p_+)^t - (1-w)(1-p_-)^t,\qquad p_\pm = \rho_0(\lambda_0 \pm \delta\lambda)$$

**Asked for:**
- Reconcile Eq. 6 with this expression.
- State whether Eq. 6 is an approximation valid over a particular parameter or fitting range.
- Show the normalized increments or the fit residuals.

This expression is also the "analytical broadcast curve" that R3 wants the four-stage behavior checked against (§1).

## 12. Status of the quadratic PI return map — *R3*

The quadratic return map for PI is **fitted** to simulated trajectories. The manuscript should not present it as a general evolution equation for the learned fraction.

## 13. Statistical uncertainty — *R1*

- The uncertainty bands are **min–max over 20 runs**. Replace or supplement them with **SD, SE, or confidence intervals**. This matters especially at high heterogeneity, where t_max distributions may be broad or skewed.
- **Return-map fits:** consecutive points from one trajectory are **correlated**, so treating them as independent understates the fit uncertainty. The authors should address this.

## 14. Mean learning rate in the abstract — *R1*

The abstract says relaxation is studied as a function of the **mean learning rate**, but λ₀ is fixed in the Methods and only δλ varies. Either **add results that vary λ₀** or **correct the abstract** and related statements.

## 15. Research gap — *R2*

The Introduction reviews CA, heterogeneous-agent, and classroom-learning work but does not state **what remains unknown about the combined effects of interaction topology and agent heterogeneity**. It also does not say how this study goes beyond prior models. The gap needs to be stated sharply.

---

## Triage: what needs new computation vs. text only

**Requires new simulations or analysis (the bulk of the effort):**
1. Stopping-window sensitivity sweep for high-δλ broadcast, with absorbed and terminated runs separated (§1). *Resolve the 20-step / 20 %-of-L² / "100" discrepancy first.*
2. Return-map comparison at fixed ρ₀, δλ, and initialization across sizes (§4). Alternatively, reframe Fig. 6 as illustrative only.
3. PI size scaling at a fixed initial learned fraction (§10).
4. Comparable-initial-fraction TI vs. PI runs (§9). Alternatively, a careful discussion in place of new runs.
5. Analytical broadcast benchmark F(t) overlaid on simulations, with increments or residuals (§11, §1).
6. Re-computed uncertainty (SE/CI) from existing trial data, plus treatment of correlation in the return-map fits (§13). *Existing data likely suffice; the `2DBPCA-IH-return-map-stats.jl` and `2DBPCA-IH-absorption-check.jl` scripts are natural starting points.*
7. *(Optional)* Quantitative percolation measures (§7), unless the terminology is softened instead.
8. *(Optional)* A λ₀ sweep (§14), unless the abstract is corrected instead.

**Text, figure, or administrative changes only:**
- Soften "optimal," "safe choice," "validation," and "percolation-driven" (§2, §3, §7).
- Explain the Hake-gain mapping and its limits, and list which features are and are not reproduced (§3).
- Fix Fig. 1: per-regime initial conditions and the steady-state wording (§6).
- Add implementation details: update scheme, λ-class assignment, seat coordinates (§5).
- Deposit code and data publicly, e.g. Zenodo/GitHub with a DOI (§5).
- Use lattice vs. population size consistently and remove "small class" (§8).
- Correct the abstract's λ₀ claim (§14).
- Sharpen the research gap in the Introduction (§15).
- Qualify the PI quadratic as a fitted description (§12).

---

## Response plan (2026-10-14)

### Precondition: fix the PI update rule before any rerun

**Discovered while checking the code; no reviewer raised it.** In `scripts/2DBPCA-IH-simulation.jl`, `generate_next_generation` uses each **neighbor's** λ (`λ_neighborhood`) in the PI learning probability. The paper's Eq. 4 and the thesis (`3BPCAIH.tex`) both specify the **learner's** λᵢⱼ. TI already uses the learner's λ.

**Decision:** fix the code to match the paper, so P = 1 − ∏[1 − λᵢⱼ ρ₀ sₙ].

**What this changes:**
- PI results at δλ = 0 are unaffected, since every λ is identical there. Those runs double as a regression test for the new code.
- The effect grows with δλ. The paper's "PI is robust to heterogeneity" claim must be re-checked, because the neighbor-λ rule may have been averaging heterogeneity out.
- Once fixed, a PI agent with one learned neighbor has exactly the TI per-step probability λᵢⱼρ₀. That gives a clean normalization argument for §9.

### Full rerun

Everything is rerun with the corrected rule, **run to full absorption** (no stopping window), and with **fixed per-run seeds**. All figures, Tables 1–3, and the numbers quoted in the text are regenerated. The two existing data subsets (δλ ≤ 0.4 with a 20-step window; δλ ∈ {0.45, 0.49} with a window of 20 % of L²) are replaced by one consistent sweep, which also resolves the "100 steps" discrepancy. Code changes are tracked separately.

### Per-theme plan

| § | Theme | Plan | Work |
|---|---|---|---|
| 1 | Stopping criterion | Remove the stopping window and run every trial to absorption (f = 1), which every run reaches with probability 1 since λ₋ > 0. Keep a safety step cap; any capped run is flagged and reported. Tables 1–2 become a single statement that all runs absorbed. Overlay R3's analytical F(t) on the TI curves (log-t axis) to show the multi-stage structure is the two-rate mixture, not a stopping or plotting artifact. Consider reporting t₉₅ alongside t_max, since at high δλ t_max is set by the slowest agent (roughly ln N / ρ₀λ₋ for TI). | Rerun + analysis + text |
| 2 | "Optimal" mixed strategy | Reframe as a hypothesis or future work in the Results, Conclusion, and abstract. | Text |
| 3 | Empirical / Hake gain | Soften "validation" and "qualitative agreement." List explicitly what is and isn't reproduced: the sign of a₂, the scatter (R² ≈ 0.63), and the seating-arrangement result that disagrees with Roxas (inner vs. outer corner). Answer the step-vs-assessment-interval point: if one interval spans k steps, the TI map composed k times is still linear with constant gain 1 − (1 − a₀)ᵏ, and the PI map stays concave. The shape class is therefore independent of the unknown mapping; only the coefficients depend on it. | Text |
| 4 | Return-map size confound | Fig. 6 actually compares the fastest corner (L=32, ρ₀=0.9, δλ=0) against the slowest (L=128, ρ₀=0.1, δλ=0.4) of the δλ ≤ 0.4 sweep. Relabel it **best vs. worst case** rather than small vs. large, in the caption, the text, and Table 3. Existing plots at fixed (ρ₀, δλ) across L (`output/2D-Binary-PCA-IH/analysis/plots/return-map/`, 45 combinations) show nearly identical shapes across sizes. Add a sentence or supplementary table saying return-map shape is governed by ρ₀ and δλ, not L. Also fix the panel label "δλ = 0.5 ± 0.4" (should be λ) and the panel order. | Text + regenerate figure |
| 5 | Implementation details | State: **synchronous** update; λ-class assignment is **independent Bernoulli(½) per agent**, fixed within a run and redrawn for each trial; **open boundaries** (the grid is padded with permanently unlearned cells); exact seed placement rule for each seating arrangement (see "Implementation (v2)" below). Fixed per-run seeds via StableRNGs.jl. Public code and data repository (Zenodo DOI) deferred until after curation, but the rerun code is written with this in mind (shared config, seeds). | Text + code |
| 6 | Flowchart | Show the per-regime initial condition (TI: 4 learned to match PI, see §9; PI: 4 seeds). Replace "steady state achieved" with "all agents learned (absorbed)." | Figure |
| 7 | Percolation | Adopt "percolation-like spreading" / "local propagation dynamics." No new percolation analysis. | Text |
| 8 | Classroom framing | Frame conclusions as model dynamics. Distinguish lattice side L from population N = L². Remove "small class" and "safe choice." | Text |
| 9 | Fairness of comparison | Concede the word "controlled" and redefine it: same λ distribution, ρ₀, lattice, and initial learned set. **TI now starts with the same 4 learned agents as PI.** Every TI run is a **twin** of a PI run: it shares the same RNG seed, so seed positions and λ-class grid are identical and only the update rule differs. The difference in exposure mechanism is intentional and is what distinguishes the two pedagogies. Argue this explicitly: in the code the two regimes differ by a single line (`probability_table`), and a PI agent with one learned neighbor has exactly the TI probability ρ₀λᵢⱼ. | Rerun (folded in) + text |
| 10 | Seed density | A **second full run set** at a **fixed initial learned fraction** 4/32², so n = 4, 9, 16, 36, 64 for L = 32 … 128. It covers **every seating arrangement** and both regimes as TI/PI twins. Seeds fill outward from each arrangement's anchor points, nearest cell first, with random tie-breaking; this is the same scheme as `grid_order` in the Classroom Dynamics Research repo. L = 32 (n = 4) is shared with the main set. An early test hint (PI, random seating, L=128, ρ₀=0.1, δλ=0, with an earlier seed draw) showed t_max dropping from ~374 to ~129 with 64 seeds, so the effect on t_max is **not** negligible. What R1 asks is whether the *scaling with L* changes, which this run set answers. | New run set |
| 11 | Broadcast return map (R3) | Use R3's F(t). For δλ = 0 it reduces to 1 − (1 − p)ᵗ, whose return map is **exactly** h(f) = f + p(1 − f) with p = ρ₀λ₀. The existing fit a₀ = 0.4454 vs. 0.9 × 0.5 = 0.45 confirms this. For δλ > 0 the effective gain decays from ρ₀λ₀ toward ρ₀λ₋ as fast learners are depleted. State that Eq. 6 is exact at δλ = 0 and an effective average otherwise, and show F(t) overlays and residuals or normalized increments. | Analysis + text |
| 12 | PI quadratic | Call it a fitted description, not an evolution equation. Also soften the highlight "Return map yields analytic Hake gain expressions." | Text |
| 13 | Statistics | Replace min–max bands with **mean ± SD** for learning curves. For t_max use **median + IQR**, since high-δλ distributions will be skewed. Return-map coefficient uncertainty: **fit each of the 20 trials separately** and report mean ± SE across trials, because trials are independent and points within a trajectory are not. | Analysis |
| 14 | λ₀ in abstract | Remove the λ₀ dependence claim from the abstract and text, and point to future work. | Text |
| 15 | Research gap | Frame the paper as a mechanistic, bottom-up model built from interaction rules, as opposed to predictive or phenomenological gain models (Hake, Nitta, Pritchard). Gap: prior classroom PCA models are homogeneous, and heterogeneous-agent studies don't compare instruction topologies. Contribution: how topology and heterogeneity jointly shape relaxation. Acknowledge that the reality-to-model mapping is an interpretive choice. | Text |

### Also to update (easy to miss)
- Highlights and abstract: "qualitative agreement," "analytic Hake gain expressions," "percolation-driven," and the λ₀ claim.
- Methods: rewrite the stopping-rule paragraph (line ~250) for run-to-absorption, and the "min–max" statements (lines ~289, 306, 322).
- Methods: specify the δλ set used in the unified sweep.
- Results: the bands are described as "min-max range across runs" but were actually mean ± SD in v1 (and are ±SD in v2). Fix lines ~289, 306, 322.

---

## Implementation (v2): status 2026-10-14

Branch **`physa-revision`**. All new code is in `-v2` files next to the v1 scripts; v1 is untouched and kept for comparison. Output goes to `output/2D-Binary-PCA-IH-v2/`. **Engine validated; both full run sets completed 2026-10-14 (commit `39ffeb5`).**

### Files
| File | Role |
|---|---|
| `scripts/2DBPCA-IH-common-v2.jl` | Shared definitions, `include`d by every v2 script: model, engine, seed placement, parameter keys and paths, Parquet I/O. |
| `scripts/2DBPCA-IH-parameter-generation-v2.jl` | Writes `parameter_combinations.csv` (the single source of the parameter grid) with per-trial seeds. |
| `scripts/2DBPCA-IH-simulation-v2.jl` | Threaded sweep over the parameter table. Resumes (skips rows whose raw file exists). |
| `scripts/2DBPCA-IH-validation-v2.jl` | Validation suite (below). Re-runnable; prints PASS/FAIL. |
| `scripts/2DBPCA-IH-plotting-v2.jl` | Figs. 3 and 5 (annotated and clean), v1-vs-v2 comparison, seed-density and fast/slow figures, classroom-evolution animations. Fig. 5 uses GLMakie (needed for correct overlap of the 3D surfaces); everything else uses CairoMakie. |

Environment: added `StableRNGs` 1.1.0 and `DuckDB` 1.5.2. DuckDB is pinned because the version the resolver picked, 1.1.0, failed to load. The pin updated many *transitive* dependencies; among directly used packages only `GR` moved (0.72.8 → 0.73.27). The full environment (Plots/Makie stack included) precompiles cleanly on Julia 1.11.0.

### Model as implemented (text for Methods, §5)
- **Lattice:** L×L, binary states, irreversible learning. **Open boundaries**: off-lattice cells count as unlearned neighbors.
- **λ classes:** each agent is independently "fast" (λ₀+δλ) or "slow" (λ₀−δλ) with probability ½. Class sizes are therefore Binomial, not exactly half. Classes are fixed within a run and redrawn for each trial.
- **Update:** **synchronous**. Every step-t transition is decided from the step-(t−1) state, then all are applied.
  - TI: P = ρ₀λᵢⱼ.
  - PI: P = 1 − (1 − ρ₀λᵢⱼ)ᵏ, where k = number of learned Moore neighbors. This uses the **learner's** λ, matching Eq. 4. v1 used each neighbor's λ.
- **Termination:** runs continue until **every agent has learned** (absorption). There is no stopping window. A safety cap of 10⁶ steps exists; any run that hits it is flagged (`absorbed = false`, T = −1 for unlearned agents). None did in validation.
- **Seed placement:** the n cells nearest the arrangement's anchor points, with random tie-breaking.
  - Anchor points: center → lattice center (L/2+½, L/2+½); outer corner → the 4 corner cells; inner corner → the 4 cells at L÷4 from each edge (rows/cols q and L−q+1); random → uniform without replacement.
  - For n = 4 this is exactly v1's placement (verified).
  - For n = 9 the split across the 4 anchors is unavoidably uneven (4 anchor cells + 5 ties).

### Parameter grid
- **Main set** (n = 4): regime {TI, PI} × SA {outer_corner, inner_corner, center, random} × L {32, 48, 64, 96, 128} × ρ₀ {0.1 … 1.0} × δλ {0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49}, with λ₀ = 0.5 → 2,800 rows.
- **Seed-density set** (n = round(4L²/32²) = 9, 16, 36, 64): same grid without L = 32 → 2,240 rows.
- 5,040 rows × 20 trials = **100,800 runs**.

### Seeds and reproducibility
- `seed_1 … seed_20` columns in `parameter_combinations.csv`, drawn from `RandomDevice`. This follows the convention of the Classroom Dynamics Research repo, including seed protection: rerunning the generator keeps existing seeds.
- **TI/PI twins share seeds** (keyed on everything except regime).
- Each run uses one `StableRNG(seed)` stream, consumed in a fixed order: seed placement → λ classes → dynamics.
- StableRNGs (rather than `Xoshiro`) guarantees the same stream across Julia versions, which matters for a public deposit.

### Output
- `raw/<regime>/<SA>-<L>-<n>/<ρ₀>-<λ₀>-<δλ>.parquet` (ZSTD, via DuckDB): one file per parameter row, one row per agent (column-major index). Columns:
  - `T_<trial>`: **learning-time grid**. 0 = seed, t = learned at step t, −1 = never learned. The full trajectory, every snapshot, and f(t) are recoverable from it.
  - `class_<trial>`: 1 = fast, 2 = slow.
  - Gitignored; regenerable from the seeds.
- `analysis/run_summary.csv`: one row per run with parameters, trial, seed, `absorbed`, `t_max`, `t_50`, `t_90`, `t_95`. f(t) is normalized by N = L². v1 normalized by the final learned count.
- Dropped from v1: per-generation state CSVs (the source of the ~95 GB), and the power-law fit (unused in the paper).
- Expected raw size: roughly 0.5 MB per parameter row at L=128, about 2 GB in total.

### Validation results (`2DBPCA-IH-validation-v2.jl`, all PASS)
1. **Seed placement:** n = 4 reproduces v1's cells for every SA and L. For n > 4, every seed is at least as close to an anchor as every non-seed.
2. **Twins:** TI and PI with the same seed have identical seed positions and λ grids. The parameter table pairs every row.
3. **Reproducibility:**
   - The same seed gives a bit-identical run, repeated and across threads.
   - The Parquet round-trip is exact.
   - The run fingerprint `f15c851058812532` is identical in separate processes with 1, 8, and 16 threads.
4. **Engine vs. a literal transcription of Eqs. 3–4:** a naive full-grid implementation (product over the neighborhood, no lookup tables) gives **bit-identical** learning-time grids in 54/54 runs (TI and PI; random, outer-corner, and center seating; δλ ∈ {0, 0.4, 0.49}).
5. **TI vs. exact analytics, 200 runs per case:**
   - The mean learned fraction is within 0.0005–0.0014 of R3's F(t) at every step.
   - Mean t_max matches the exact expected maximum of independent geometric times (|z| ≤ 1.9). Example: L=96, ρ₀=0.1, δλ=0.49 gives 9,018 simulated vs. 9,008 exact.
6. **PI vs. v1 data at δλ = 0** (where the λ fix is a no-op), 5 cases, v1 with 20 trials and v2 with 100 runs: t_max and t₅₀ agree (|z| ≤ 1.1 for t_max; ≤ 2.2 for t₅₀, the largest being outer_corner L=128 ρ₀=0.1: 209.8 vs 206.5).
7. **Absorption:** the slowest corner (L=128, ρ₀=0.1, δλ=0.49) absorbs in every run. TI t_max is 7,400–12,700 and PI t_max is 1,600–2,300, far below the 10⁶ cap. The worst case runs in about 0.02 s.

### Early indications (small test subsets only; not for the paper)
- **Effect of the λ fix (PI, δλ = 0.4, vs. v1):** t₅₀ is unchanged. t_max is 10–25% longer (e.g. center L=64 ρ₀=0.3: 102.8 → 115.0), because slow learners no longer borrow neighbors' fast λ. Curve shapes should hold; the high-δλ tail is what changes.
- **Run-to-absorption at L=128, ρ₀=0.1, δλ=0.49:** TI t_max ≈ 9,500 vs. PI ≈ 1,750 (t₉₅ ≈ 2,300 vs. 650–780). PI's advantage at high heterogeneity grows once runs are not truncated.

### Quick look at the v2 data vs. the paper's claims (2026-10-14; means of 20 runs, not final analysis)

| Claim in the paper | v2 verdict | Evidence |
|---|---|---|
| Fig. 5: PI wins at high δλ (gap grows with δλ); TI can win at large L, low δλ, high ρ₀ | **Holds** | TI/PI winner flips in only 5 of 350 (L, ρ₀, δλ) cells (inner corner), all near-ties, e.g. L=32 ρ₀=0.9 δλ=0: 12.5 vs 12.6. TI wins in 2/70 cells at L=32 and 25/70 at L=128. |
| Asymmetry: "when TI wins the difference is small, when PI wins it can be very large" → PI is a "safe choice" | **Only true in absolute steps** | Max excess when TI wins: 170 steps. When PI wins: 8,068 steps. **As ratios it's nearly symmetric**: PI up to 5.7× slower (all SAs; 2.5× for inner corner, e.g. L=128 ρ₀=0.9 δλ=0: 17.6 vs 44.0) vs. TI up to 6.1× slower. Reviewers already asked to drop "safe choice"; restate the asymmetry as absolute, or drop it. |
| TI "largely unaffected by L" | **Mostly holds; refine wording** | TI t₅₀ = 3.0 at every L, but t_max grows slowly, 42 → 60 for L=32 → 128 (ρ₀=0.5, δλ=0.2). That's the logarithmic growth of the slowest agent's learning time. |
| PI: larger L delays onset, shape preserved | **Holds** | Inner corner t₅₀ 9.8 → 33.9 and t_max 20 → 63 for L=32 → 128, matching v1 (20.7 → 63.0). |
| Inner corner best; center ≈ outer corner; random widest spread | **Holds strongly** | Inner corner is best in 332/350 cells. At L=64 ρ₀=0.5 δλ=0.2: inner 35.1, center 60.9, outer 58.2, random 61.5 ± 8.8 (largest SD). |
| ρ₀ effect nonlinear; TI earlier onset; PI catches up and overtakes | **Holds** | L=64 δλ=0.2: PI t_max 130 / 35 / 24 for ρ₀ = 0.1 / 0.5 / 0.9. TI t₅₀ < PI t₅₀ at every ρ₀, while PI t_max < TI t_max at every ρ₀. |
| **PI: heterogeneity gives "only a gradual uniform delay… no qualitative change in curve shape"; percolation "averages over" λ** | **Does not hold — changed by the bug fix** | PI now develops a slow-learner tail. L=64 ρ₀=0.5 inner corner, δλ 0 → 0.49: t₅₀ 17 → 31, t₉₅ 26 → 107, t_max 32 → 311 (v1 had 195). t_max/t₅₀ goes from 1.9 to 10. The "averaging" was largely the neighbor-λ bug. PI is **still far more robust than TI** (TI t_max/t₅₀ ≈ 276, t_max 1,603), but for a different, genuine reason: a slow learner with k learned neighbors learns with probability 1 − (1 − ρ₀λ₋)ᵏ, up to ~8× its solo rate. This mechanism needs rewriting in §3.2 and the abstract. |
| PI t_max at high δλ (all Fig. 5 values) | **Numbers change** | The fix raises PI t_max at δλ ≥ 0.4, e.g. L=32 ρ₀=0.1 δλ=0.49: 437 → 1,317. Running to absorption raises TI where v1 truncated, e.g. L=32 ρ₀=0.1 δλ=0.4: 377 → 659. All figures and quoted numbers need regenerating, but no qualitative Fig. 5 conclusion flips. |
| TI "four-stage" relaxation at δλ=0.49 "not fully understood"; stopping-rule Table 1 | **Explained / moot** | TI matches R3's two-rate F(t) exactly (validation check 5). On a log-time axis a two-rate mixture rises, plateaus at ~½, rises again, then saturates, which is the fast–slow–fast–slow pattern. Table 1 is replaced by "all 100,800 runs absorbed". |
| Size dependence of PI (R1 seed-density concern) | **Depends on seating** | Log-log slope of PI t_max vs L, n=4 vs fixed fraction: inner corner 0.81 vs 0.78; center 0.91 vs 0.89; outer corner 0.93 vs 0.89 (ρ₀=0.5 δλ=0.2; similar at ρ₀=0.1 δλ=0.4). **Random: 0.91 vs 0.18** (0.65 vs 0.09 at ρ₀=0.1 δλ=0.4). For clustered seating the size dependence is genuine: the front must cross a distance ~L. For random seating it is mostly seed dilution, as R1 suspected. |
| Uncertainty / skew (R1 §13) | **Mild skew** | L=128 ρ₀=0.1 δλ=0.49: TI t_max mean 9,818, median 9,662, SD 1,281, range 7,754–13,171. PI: mean 1,750, median 1,688, range 1,386–2,777. Right-skewed; median + IQR is appropriate. |
| Return map (Fig. 6, Table 3) | **Not checked yet** | PI at δλ=0 is unchanged by the fix (Fig. 6b case). The Fig. 6a case (L=128 ρ₀=0.1 δλ=0.4) shifts modestly (PI t_max 272 → 332). Needs the per-trial refit planned for §13. |

### Regenerated figures and mechanism exploration (2026-10-14, v2 data)

All figures live under `output/2D-Binary-PCA-IH-v2/analysis/plots/`. They are generated by `scripts/2DBPCA-IH-plotting-v2.jl`; animations come from `animate_evolution(...)` in the same script.

| Figure | File | Purpose |
|---|---|---|
| Fig. 3 (learning curves), annotated like the paper | `learned-t/figure3-learned-t-v2.png` | Check against the published version |
| Fig. 3, **clean** (no series labels; 8000×6000, same canvas as v1 `size-SA_ρ₀-δλ_quad.png`) | `clean/figure3-quad-clean.png` | Drop-in replacement in PowerPoint |
| Fig. 5 (t_max surfaces), annotated | `figure5-tmax-surface-v2.png` | Check against the published version |
| Fig. 5, **clean** panels (2500×2000, same as v1 `{32,128}_truncated.png`) | `clean/figure5-{32,128}-clean.png` | Drop-in replacements in PowerPoint |
| v1 vs v2 learning curves (L=64, ρ₀=0.5) | `learned-t/learned-t-v1-vs-v2.png` | Shows the effect of the fix; TI overlays R3's F(t) |
| Seed density, 4 seeds vs fixed fraction (R1) | `seed-density-comparison.png` | Candidate supplementary figure |
| Fast vs slow learners, TI vs PI | `class-split-learning-curves.png` | Mechanism figure; candidate paper panel |
| Classroom evolution, fastest vs slowest run | `evolution/*.mp4` and `*-strip.png` | Exploration only |

**Figures 3 and 5: qualitatively unchanged.**
- Panels (a)–(c) of Fig. 3 are practically identical to the published version.
- Fig. 3(d): PI at δλ = 0.45 and 0.49 now has a slow-learner tail, finishing around 2⁸ instead of about 2⁶·⁵.
- Fig. 5 keeps the same shape and the same L=128 crossover. The TI peaks are taller because runs are no longer truncated: about 650 vs 377 at L=32, and about 950 vs 650 at L=128.
- **Text error found:** the paper says the bands are "min-max range across runs", but the v1 plotting script drew **mean ± SD**. v2 also uses ±SD, so the text needs correcting (lines ~289, 306, 322).

**Mechanism findings.** These feed the text rewrites and the responses.

1. **Seed coverage sets t_max** at low to moderate δλ.
   - Learning fronts in the 3×3 (Moore) neighborhood advance at a roughly constant speed. t_max therefore tracks the *covering radius*: the largest distance (counted in 3×3 neighborhood steps, so diagonals count as one) from any agent to its nearest seed.
   - For random seating across 20 runs, the correlation between covering radius and t_max is 0.82 (L=64) and 0.88 (L=128) at ρ₀=0.5, δλ=0.2. It is 0.79 at ρ₀=0.1, δλ=0.4. Spreading takes about 1.5 steps per cell at ρ₀=0.5 and about 6 steps per cell at ρ₀=0.1, δλ=0.4.
   - **This explains the seating-arrangement ranking (§3.1).**
     - Inner corner places its 4 seeds at the quadrant centers. That is the optimal 4-point covering of a square, with radius L/4.
     - Center and outer corner both leave a farthest point about L/2 away, hence their near-identical t_max.
     - This should replace the paper's "consequence of the isotropic coupling assumption".
   - **It also frames the Roxas discrepancy (R2, §3):** in the model only coverage geometry matters, and real classrooms have factors the model omits.
   - Random seating's large run-to-run spread is luck of coverage. Clumped seeds waste a seed and leave a far corner, as the `PI-random-64-4-0.5-0.2` animation shows.

2. **Seed density (R1, §10) depends on seating arrangement.**
   - Log–log slope of PI t_max vs L, 4 seeds vs fixed fraction:
     - Clustered arrangements are essentially unchanged across all 70 (ρ₀, δλ) cells. For example, inner corner 0.81 → 0.78 and outer corner 0.93 → 0.89.
     - **Random collapses**: 0.91 → 0.18, and fixed-fraction slopes are about 0.1–0.2 everywhere. That is close to TI's slopes of 0.2–0.3, which come from TI's t_max growing roughly with ln N.
   - Interpretation:
     - Clustered seating has a genuine size dependence, because the front must travel a distance proportional to L from the anchors.
     - For random seating at fixed density the covering radius stops growing with L, so R1's concern is valid for that arrangement.
   - **Consequence:** for random seating, the TI-vs-PI winner at large L depends on the seeding scheme. At L=128, ρ₀=0.5, δλ=0.2, with 4 seeds TI wins (about 58 vs 112). With a fixed fraction PI wins (about 39 vs 58). Fig. 5 uses inner corner, where seeding doesn't matter, so its conclusion stands, but the R1 response should state this.
   - Candidate supplementary figure: `seed-density-comparison.png`, panels (a) and (d).

3. **At high δλ, PI relaxes in two stages, and t_max becomes one straggler's waiting time.** At ρ₀=0.1, δλ=0.49:
   - **Stage 1, to about t=150–250:** fronts branch along connected fast learners and bypass slow learners.
   - **Stage 2, to t_max:** slow learners are mopped up. Even with 8 learned neighbors, a slow learner learns with only about 0.8% probability per step, and clusters of slow learners take longest.
   - For random seating, coverage still predicts **t₉₅** (corr 0.84 at L=128) but **not t_max** (corr −0.02 at L=128 and −0.12 at L=64).
   - At this δλ, t_max measures the unluckiest slow learner, not the dynamics. This is strong support for **reporting t₉₅ alongside t_max** (§1, §13).

4. **Fast vs slow learners: the actual robustness mechanism.** From `class-split-learning-curves.png`, inner corner, L=64:
   - In TI, fast and slow learners move on separate timescales (ratio ≈ λ₊/λ₋, about 100× at δλ=0.49). The whole-class curve is just their average, which is the plateau in Fig. 3(d).
   - In PI, the slow-learner curve trails the fast-learner curve by only a factor of about 2–3 at δλ=0.49, and is indistinguishable at δλ=0.2. Fast learners propagate the front, and slow learners are then surrounded by learned neighbors, so they learn with probability 1 − (1 − ρ₀λ₋)ᵏ.
   - **TI's fast learners always finish before PI's fast learners. At high δλ, PI's slow learners finish well before TI's slow learners.** PI's advantage comes entirely from the slow learners.
   - Use: this is the replacement for the "percolation averages over λ" sentence (§3.2, abstract) that the bug fix invalidated. The fast/slow lag (e.g. ratio of the two groups' t₅₀) is a candidate single-number equity/robustness metric. It is also the concrete motivation for proposing a mixed TI/PI strategy as **future work** (§2), without claiming it is optimal.

5. **Possible quantitative percolation link (§7): unverified hypothesis.**
   - Fast learners occupy half the sites. Site percolation on the square lattice with Moore neighbors has threshold p_c ≈ 0.41, so the fast subpopulation should form a spanning cluster, which is what stage 1 runs along.
   - Test: vary the fast fraction across p_c and check whether stage 1 stops sweeping the room below it. That would be a model extension, since the model currently fixes the fraction at ½.
   - Optional. The agreed §7 response (softer "percolation-like" wording) does not depend on it.

### Next
- [x] Review the v2 code and validation, then commit on `physa-revision` (`39ffeb5`).
- [x] Run the main set and the seed-density set. 100,800 runs in **2 min 23 s** wall time on 16 threads (Ryzen 7 7700), 1,254 s of summed simulation time, median **~6 ms per run** in both regimes. **All runs absorbed**, none hit the cap. Raw output 719 MB. Validation was re-run just before the sweep (all PASS), and 40 randomly sampled stored runs re-simulated from their seeds match exactly. Timings are in `analysis/run_log.txt` and the `runtime_s` column of `run_summary.csv`.
- [ ] Rewire analysis and plotting to read `run_summary.csv` and the Parquet grids. This includes the per-trial return-map fits (§13), F(t) overlays (§1, §11), and mean ± SD / median + IQR bands.
- [x] Regenerate Figs. 3 and 5 from v2 data, plus clean copies at v1 dimensions for PowerPoint re-labeling.
- [ ] Re-label the clean Figs. 3 and 5 in PowerPoint (manual). PI δλ=0.49 in Fig. 3(d) now ends further right; the Fig. 5 z-ranges are taller.
- [ ] Return map (Fig. 6, Table 3): per-trial refits on v2 data, best/worst-case relabeling (§4), R3's F(t) comparison (§11). Not checked yet.
- [ ] Add t₉₅ alongside t_max in the summary analysis and figures (mechanism finding 3).
- [ ] Decide which exploration figures go into the paper or supplement: seed density (a)+(d) for R1 §10; fast/slow class split for §3.2.
- [ ] Text rewrites driven by v2:
  - §3.2 heterogeneity mechanism and the abstract (finding 4).
  - Seating-arrangement explanation via coverage (finding 1).
  - "TI unaffected by L" → t₅₀ unchanged, t_max grows ~ln N.
  - Asymmetry / "safe choice" (absolute vs ratio).
  - Four-stage → two-rate F(t).
  - Bands min–max → ±SD.
- [ ] Optional: fast-fraction percolation check (finding 5).
