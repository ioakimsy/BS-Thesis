# Response to Reviewers

**Manuscript:** PHYSA-263545, "Learning dynamics as a heterogeneous cellular automaton: comparing broadcast and nearest-neighbor interaction regimes in a classroom model"

**Authors:** C. I. T. Sy, R. M. L. Roxas-Villanueva, J. Y. Bantang

We thank the three reviewers for their careful reading and constructive comments. Their concerns overlapped substantially, and addressing them led us to rerun every simulation and to revise the analysis and presentation throughout. Below we first summarize the main changes, then reply to each comment in turn. Reviewer comments are quoted as indented blocks, abridged with [...] where long. Where a comment was raised by more than one reviewer, each reply is self-contained, with a cross-reference to the other replies.

All equation, figure, table, and section numbers below refer to the revised manuscript unless stated otherwise.

## Summary of main changes

1. **No stopping rule.** Every simulation now runs until all agents have learned (the absorbing state). All 100,800 runs reached it; none came near the safety cap of 10⁶ steps. The former Tables 1 and 2 (runs terminated by the stopping rule) are removed.
2. **Matched comparison.** Broadcast and nearest-neighbor runs now start from the same four learned agents, with the same assignment of learning rates and the same random seed, so that paired runs differ only in the update rule.
3. **New run set at a fixed initial learned fraction** (44,800 runs), to separate system size from seed density (Appendix C).
4. **Corrected implementation.** While rerunning, we found that the code computed the nearest-neighbor probability with the neighbors' learning rates rather than the learner's own, as Eq. (4) specifies. The equation in the manuscript was correct; the code now matches it. Results at δλ = 0 are unaffected, and the comparison between the two regimes is unchanged.
5. **New analyses:**
   - the time t₉₅ for 95% of agents to learn, alongside t_max (Eq. 5, Table 2, Fig. 6, Appendix A);
   - the exact broadcast gain for a heterogeneous population (Eq. 8, Fig. 8, Table 4);
   - the learning of fast and slow agents tracked separately (Fig. 4);
   - return-map coefficients at fixed parameters across lattice sizes (Fig. 7, Table 3, Appendix B);
   - a class-level comparison with the empirical data (Fig. 9).
6. **Revised claims.** The "optimal" mixed strategy is now a hypothesis; "percolation-driven" is replaced by "local propagation"; "safe choice" is replaced by measured differences; the empirical comparison is stated as consistency rather than validation; the claimed dependence on the mean learning rate is removed.
7. **Reproducibility.** The update scheme, learning-rate assignment, boundary conditions, and seating rules (Table 1) are specified, every run's random seed is recorded, and code and data are deposited in a public repository [repository URL and Zenodo DOI to be added].

---

## Reviewer 1

> This manuscript presents a probabilistic cellular automaton model to compare two classroom interaction regimes: a broadcast regime representing traditional instruction and a nearest-neighbor regime representing peer instruction. [...] I therefore recommend major revision.

We thank the reviewer for the detailed assessment. We have addressed each point below; in particular, both regimes now start from the same initial conditions, the size scaling is repeated at a fixed initial learned fraction, the stopping rule is removed, and the percolation and Hake-gain claims are moderated.

### Comment 1.1: Controlled comparison

> The comparison between the two interaction regimes is not fully controlled. The broadcast simulations start with all agents in the unlearned state, whereas the nearest-neighbor simulations start with four agents already in the learned state. In addition, the transition rule is different in the two cases. [...] The authors should clarify what is meant by a controlled comparison between the two interaction regimes. It would be useful to include simulations with comparable initial learned fractions, or otherwise discuss explicitly how the different initialization and update rules influence the conclusions.

We agree that the original initialization differed between the two regimes. In the revised simulations:

- **Both regimes start with the same four learned agents**, placed by the same seating rule.
- **Each broadcast run is paired with a nearest-neighbor run** that uses the same random seed, so the two share the initial learned agents and the assignment of learning rates. Paired runs therefore differ only in the update rule (Section 2).

The difference in update rule is intended: it is the contrast between the two modes of instruction. We now state this explicitly (Section 2). The rules differ in which agents are exposed (all agents in the broadcast regime, only neighbors of learned agents in the nearest-neighbor regime) and in how exposure accumulates. For an agent with exactly one learned neighbor, the two rules give the same probability, ρ₀λᵢⱼ, so each learned neighbor acts as a source of the same strength as the broadcast source.

Four initially learned agents are the minimum that seeds every seating arrangement symmetrically. They make up at most 0.4% of the population and have a negligible effect in the broadcast regime, where they do not influence other agents.

Following the reviewer's suggestion, we also repeated the sweep with a comparable initial learned fraction across lattice sizes (see Comment 1.2). In the broadcast regime, this leaves the size dependence of t_max unchanged (change in scaling exponent 0.00 on average; Table 5).

### Comment 1.2: Seed density and system size

> The system-size dependence of the nearest-neighbor regime may be strongly affected by the fixed number of initial learned agents. [...] I suggest repeating at least part of the size-scaling analysis using a fixed initial fraction of learned agents rather than a fixed number of four seeds. This would help determine whether the reported size dependence is a genuine property of nearest-neighbor dynamics or mainly a consequence of decreasing seed density.

We repeated the full sweep for L ≥ 48 with the initial learned fraction fixed at that of the L = 32 lattice, i.e. with 4L²/32² initially learned agents (9, 16, 36, and 64), placed by the same seating rules (Section 2; 44,800 additional runs).

We quantify the size dependence by the scaling exponent γ of t_max ∝ L^γ between L = 32 and L = 128 (Section 3.1). The answer depends on the seating arrangement (Appendix C, Table 5, Fig. 12):

| Arrangement | γ, four seeds | γ, fixed fraction | Mean change over all ρ₀, δλ |
|---|---|---|---|
| Inner corner | 0.81 | 0.78 | −0.03 |
| Center | 0.91 | 0.89 | −0.03 |
| Outer corner | 0.93 | 0.89 | −0.04 |
| Random | 0.91 | 0.18 | −0.62 |

(γ values at ρ₀ = 0.5, δλ = 0.2.)

- **Clustered arrangements:** the size dependence is essentially unchanged. The learned front must cross a distance proportional to L, so the size dependence is intrinsic to nearest-neighbor propagation.
- **Random arrangement:** a fixed fraction removes most of the size dependence. There, the size dependence comes mostly from the dilution of the initially learned agents, as the reviewer anticipated.

The return-map coefficients behave the same way: their dependence on L changes by less than 0.01 in the exponent when the fraction is held fixed (Section 3.4).

### Comment 1.3: Mean learning rate

> The treatment of the mean learning rate should be clarified. The abstract states that the relaxation time is studied as a function of population size, mean learning rate, learning-rate heterogeneity, and coupling strength. However, in the Methods section the mean learning rate appears to be fixed [...]

The reviewer is correct: λ₀ = 0.5 throughout. We have revised the abstract and the Introduction to state that the relaxation time is studied as a function of lattice side length, heterogeneity about a fixed mean learning rate, and coupling strength. Varying λ₀ is now listed as a future extension (Section 4).

### Comment 1.4: Stopping criterion

> The stopping criterion may affect the reported relaxation times. [...] The authors should provide a sensitivity analysis using different stopping-window lengths and should distinguish clearly between fully absorbed runs and prematurely terminated runs in the relaxation-time analysis.

We have removed the stopping criterion. Every run now continues until all agents have learned, so t_max is the absorption time rather than a lower bound. Since every agent has a nonzero transition probability, absorption occurs with probability one. All 100,800 runs were absorbed, none near the safety cap of 10⁶ steps (Section 2). The former Tables 1 and 2 are removed, and no sensitivity analysis is needed.

The multi-stage relaxation at δλ = 0.49 is unaffected. It is the two-rate relaxation of the high- and low-λ subpopulations: the mean of the simulated runs agrees with the exact broadcast result F(t) (Reviewer 3, comment 1) to within 1.5 × 10⁻³ at every time step (Section 3.2).

All figures and quoted values were regenerated from the new runs. The relative timing across parameter values and between the two regimes is unchanged: across all 350 combinations of lattice size, coupling strength, and heterogeneity, the faster regime changed in only five, all near-ties.

During the rerun we also corrected the implementation error described in item 4 of the summary. Its main visible effect is that, at the strongest heterogeneity (δλ ≥ 0.45), the nearest-neighbor relaxation curve keeps its shape but approaches full learning more slowly, which we now state (abstract, Section 3.2, Section 4).

### Comment 1.5: "Percolation-driven"

> The use of the term "percolation-driven" requires stronger quantitative support. [...] I suggest either providing additional quantitative percolation analysis or using more cautious terminology such as "percolation-like spreading" or "local propagation dynamics."

We agree that our results show local stochastic spreading, not a percolation transition. We have replaced "percolation" throughout with "local propagation", defined in Section 3.1 as the spread of the learned state outward from the initially learned agents through contact between neighbors. The keyword is now "local propagation dynamics".

### Comment 1.6: Return map, Hake gain, and empirical comparison

> The interpretation of the return-map analysis and Hake gain requires further justification. [...] the physical or educational meaning of one simulation step relative to the interval between assessments is not sufficiently established. [...] the statement that the model shows qualitative agreement with empirical classroom data should be moderated. [...] the small-system and large-system examples used in the return-map comparison appear to involve different parameter values. [...] A direct comparison at fixed coupling strength, heterogeneity, and initialization would be more informative.

We address the three parts in turn.

**Meaning of one time step.** One time step is one learning step: the interval over which an unlearned agent learns with probability ρ₀λᵢⱼ from the broadcast source, or from each learned neighbor. A step therefore has no fixed duration in real time, and the model does not fix how many steps lie between two assessments (Section 2). We show that this does not affect the comparison we make. If the interval between assessments spans k steps, the k-step gain (Eq. 10) is increasing in f whenever the one-step gain is, and a constant one-step gain gives a constant k-step gain. The shape of the gain is therefore independent of k, while its magnitude is not, so we compare shapes only (Section 3.4).

**Empirical comparison.** We agree that the original wording overstated the agreement, and we have reanalyzed the data:

- The original fit treated each student as a point, whereas the learned fraction f of the model is a class-level quantity. Within every class, the score change decreases with the pre-test score. This is expected from regression to the mean and from the ceiling on scores, although the data cannot separate these from a genuine dependence on prior score. It is this within-class trend that makes the student-level fit convex, and the resulting coefficient is not significantly different from zero (a₂ = 0.142 ± 0.086).
- At the class level (new Fig. 9), the quadratic coefficient is negative, a₂ = −0.67 ± 0.34, matching the sign of the nearest-neighbor regime. The class Hake gain increases with the class pre-test score (r = 0.80). This excludes class 9, whose mean score decreased, which an absorbing learned state cannot produce. With class 9 included, the coefficient remains negative but is not constrained (a₂ = −0.20 ± 0.67, r = 0.25).
- We also note that this rests on fourteen classes falling into two groups of pre-test scores. We therefore describe the class-level data as consistent with, rather than confirming, the nearest-neighbor gain shape.
- We now list what the model does not reproduce: the score decreases in the data, the large student-level scatter (R² = 0.63), and the empirical seating-arrangement ranking.
- "Qualitative agreement" and "validate" are replaced by these specific statements in the abstract, Introduction, Section 3.4, and Section 4.

**Size comparison.** The reviewer is right that the two original examples differed in coupling strength and heterogeneity as well as in size, so the size attribution was not supported. Fig. 7 now shows L = 32 and L = 128 at the same ρ₀ and δλ in each row, so the panels in a row differ only in size. Table 3 gives the fitted coefficients. With this comparison:

- The broadcast gain does not depend on L, as expected for agents that transition independently. The earlier statement that it decreases with L came from the change in ρ₀ and δλ and has been removed.
- The nearest-neighbor coefficients do decrease with L at fixed parameters, with |a₂| ∝ L^−β and β between 0.5 and 0.9 (Appendix B, Fig. 11), also at a fixed initial learned fraction.

### Comment 1.7: Statistical uncertainty

> The statistical uncertainty should be presented more carefully. [...] I suggest reporting standard deviations, standard errors, or confidence intervals. [...] For the return-map fitting, the authors should also consider whether the fitted data points can be regarded as statistically independent, since consecutive points from the same simulation trajectory are correlated.

**Bands.** On checking, the bands in the original figures already showed ±1 standard deviation across runs; the text incorrectly described them as min–max ranges. We have corrected the text. Section 2 now states that bands show one standard deviation across runs, while uncertainties of fitted coefficients are standard errors of the mean, and each caption states which is shown.

**Skewness of t_max.** In the highly heterogeneous regime, the run-to-run distributions of t_max are only mildly skewed (mean/median ≤ 1.12), so we report means. Fig. 10 shows t_max and t₉₅ with their standard deviations up to δλ = 0.49.

**Return-map fits.** We agree that consecutive points of one trajectory are correlated. We now fit each of the 20 runs separately and report the mean of the coefficients ± the standard error across runs, treating the independent runs rather than the points as samples (Section 3.4, Table 3, Fig. 11).

### Comment 1.8: Implementation details and code

> Several implementation details should be clarified to improve reproducibility. The manuscript should explicitly state whether the cellular automaton is updated synchronously or asynchronously. [...] Finally, since the study is primarily computational, depositing the simulation and analysis code, together with representative datasets, in a permanent public repository would improve reproducibility.

We have added the following to Section 2:

- **Update:** synchronous. All transition probabilities at step t are computed from the states at step t − 1, and all transitions are applied together.
- **Learning rates:** each agent is assigned to the high or low group independently with probability ½. The two groups are therefore spatially random and equal in size only on average. Assignments are fixed during a run and redrawn for each run.
- **Seating arrangements:** Table 1 gives the rule for each arrangement. Initially learned agents are the sites nearest to anchor points (the lattice corners, the points (L/4, L/4) and their symmetric counterparts, or the lattice center), with random tie-breaking; the random arrangement draws them uniformly.
- **Boundaries:** open. Sites outside the lattice count as unlearned neighbors.
- **Reproducibility:** every run uses its own random number stream from a recorded seed and can be reproduced exactly.

The simulation and analysis code and the simulation output, including all seeds, are deposited at [repository URL and Zenodo DOI to be added] (Code and Data availability statements).

### Comment 1.9: Flowchart terminology; lattice vs population size

> The flowchart refers to a "steady state achieved" condition, while the Methods section explicitly states that the no-transition stopping rule does not guarantee a true steady state. The terminology should be made consistent. The manuscript should also distinguish more carefully between lattice size and population size. Descriptions such as "small class" may be misleading [...]

- **Flowchart:** Fig. 1 now ends when all agents have learned, consistent with running to absorption. It also shows four initially learned agents. The caption states that updates are synchronous and which equation gives the probability in each regime.
- **Size terminology:** we now distinguish the lattice side length L from the population N = L² and use "lattice side length" throughout. "Small class" and "large class" are removed.

### Comment 1.10: "Optimal" combination and "safe choice"

> The conclusion that combining broadcast and nearest-neighbor instruction would be "optimal" is not directly demonstrated in the present study. [...] Similarly, statements suggesting that nearest-neighbor interaction is a generally "safe choice" should be softened [...]

**Combined strategy.** We agree. The combined strategy is now presented as a hypothesis motivated by the complementary strengths of the two regimes, not tested in this study, in Sections 3.2, 3.3, and 4 and in the future work.

**"Safe choice."** This is removed. In its place we report the measured asymmetry (Section 3.3). Within the range of Fig. 6, the broadcast regime is faster by at most 27 iterations, whereas the nearest-neighbor regime can be faster by more than 600. As ratios, however, the two are comparable: at most 2.8 times slower for the nearest-neighbor regime and 3.9 times for the broadcast regime.

We also report the time t₉₅ for 95% of agents to learn, which shows that the answer depends on the measure. For the larger lattice, the broadcast regime reaches 95% sooner over most of the parameter range, while the nearest-neighbor regime leaves fewer agents far behind (Eq. 5, Table 2, Fig. 6, Appendix A). All such conclusions are framed as properties of the model dynamics, consistent with the reviewer's point about the simplified binary model (Section 2).

> Overall, the manuscript contains potentially interesting results, but several important methodological issues need to be addressed [...]

We thank the reviewer again. As summarized above, we have matched the initial conditions and documented the update rules, tested seed density, removed the stopping criterion, and moderated the percolation and Hake-gain interpretations.

---

## Reviewer 2

> This manuscript presents a probabilistic cellular automaton model for classroom learning dynamics [...] However, several aspects would benefit from further clarification and more cautious presentation.

We thank the reviewer for the positive assessment and the helpful comments, which we address below.

### Comment 2.1: Research gap

> The research motivation is generally clear, but the specific research gap should be stated more explicitly. [...] The authors should clarify what is still unknown about the combined effects of interaction topology and agent heterogeneity, and how the present study addresses this gap beyond previous models.

We added a paragraph to the Introduction that states the gap explicitly. It reviews what existing models cover and what they leave open:

- **Hake-gain models** (Hake; Pritchard et al.; Nitta) describe class-level score changes without specifying which students interact, so they cannot distinguish interaction topologies.
- **Agent-based classroom models** (Bordogna and Albano; Ormazábal et al.; Velásquez-Rojas and Laguna) include students who differ in prior achievement, attention, or aptitude, and some compare lectures alone with lectures supplemented by group work. In these models, however, peers interact in small groups or clusters without a spatial arrangement, and group work supplements the lecture rather than replacing it, so broadcast and local interaction are not compared as alternatives.
- **Seating studies** address spatial placement, and our earlier proceeding compared broadcast and nearest-neighbor coupling on a lattice, but with identical agents. Comparisons of local and mean-field spreading outside the classroom (e.g. Fukś and Boccara) likewise assume identical agents.

The open question is therefore how learning-rate heterogeneity interacts with interaction topology on a spatial lattice:

- whether heterogeneity affects broadcast and local coupling differently (Section 3.2, Fig. 4);
- how the resulting relaxation times depend on system size (Section 3.1, Appendix C);
- how long the slowest agents are left behind (t₉₅ versus t_max, Section 3.3).

We also state that the model is bottom-up, with class-level dynamics following from the interaction rules, and that the correspondence between model quantities and classroom observables is an interpretive choice. New references: Ormazábal et al. (2021), Velásquez-Rojas and Laguna (2021), Koponen et al. (2017), Singh et al. (2009), Cirillo et al. (2024), and Fukś and Boccara (1997).

### Comment 2.2: "Percolation-driven"

> The interpretation of the nearest-neighbor dynamics as "percolation-driven" needs stronger quantitative support. [...] The authors should either provide additional evidence or use more cautious terminology such as "percolation-like spreading."

We agree; our results show local stochastic spreading, not a percolation transition. "Percolation" is replaced throughout by "local propagation", defined in Section 3.1 as the spread of the learned state outward from the initially learned agents through contact between neighbors (see also Reviewer 1, comment 1.5).

### Comment 2.3: Stopping criterion

> The robustness of the stopping criterion should be examined. At high heterogeneity, especially for delta lambda = 0.49, a non-negligible fraction of broadcast simulations terminate because of the stopping rule rather than full absorption. [...]

We have removed the stopping criterion. Every run now continues until all agents have learned, and all 100,800 runs reached this absorbing state (Section 2). The multi-stage relaxation at δλ = 0.49 is therefore not affected by premature termination. It matches the exact two-rate broadcast result to within 1.5 × 10⁻³ at every time step (Section 3.2). All figures and values were regenerated, and the main conclusions are unchanged (see also Reviewer 1, comment 1.4).

### Comment 2.4: Connection to real classrooms

> The connection between the idealized model and real classroom instruction should be presented more cautiously. [...] The conclusions should therefore be framed primarily in terms of model dynamics rather than direct educational effectiveness.

We agree. Section 2 now states that the two interaction topologies are deliberately simplified representations of traditional and peer instruction, and that the comparisons concern the relaxation dynamics of the model, not the educational effectiveness of the two modes.

We use "lattice side length" rather than "class size" throughout and have removed "small class" and "safe choice". Conclusions about which regime is faster are now stated as measured, model-level differences, including t₉₅ alongside t_max (Section 3.3). The proposed combination of regimes is presented as a hypothesis (see Comment 2.6).

### Comment 2.5: Empirical comparison

> The empirical comparison should be described as qualitative rather than as strong validation. The empirical data show substantial scatter, the fitted quadratic coefficient differs in sign from the nearest-neighbor simulations, and the seating-arrangement result also differs from previous empirical findings. [...]

We agree and have revised the comparison (Section 3.4):

- **Sign of the quadratic coefficient.** The original fit treated individual students as points, while the learned fraction of the model is class-level. Within each class the score change decreases with pre-test score, consistent with regression to the mean and the ceiling on scores, and this within-class trend makes the student-level fit convex (a₂ = 0.142 ± 0.086, not significantly different from zero). Fitted to class means (new Fig. 9), the coefficient is negative, a₂ = −0.67 ± 0.34, matching the nearest-neighbor regime. Given the fourteen classes and their two groups of pre-test scores, we describe this as consistent with, rather than confirming, the nearest-neighbor gain shape.
- **Discrepancies, stated explicitly.** The model does not reproduce the score decreases present in the data, the student-level scatter (R² = 0.63), or the empirical seating-arrangement ranking (inner versus outer corner), which we attribute to the isotropic coupling and the absence of aptitude-similarity effects.
- **Wording.** "Validate" and "qualitative agreement" are replaced by these specific statements in the abstract, Introduction, and Section 4 (see also Reviewer 1, comment 1.6).

### Comment 2.6: "Optimal" mixed strategy

> The claim that a mixed broadcast and nearest-neighbor strategy would be "optimal" is not directly supported by the current simulations. [...]

We agree. The combined strategy is now presented as a hypothesis motivated by the complementary strengths of the two regimes and is not tested here (Sections 3.2, 3.3, and 4). The new t₉₅ analysis makes these strengths concrete. The broadcast regime brings most agents to the learned state sooner on larger lattices, while the nearest-neighbor regime leaves fewer agents far behind (Section 3.3). A mixed-regime model that tests the hypothesis is listed as future work.

---

## Reviewer 3

> The manuscript examines the effects of learning-rate heterogeneity on classroom learning dynamics using a probabilistic cellular automaton model. [...] My specific comments are given below.

We thank the reviewer for the careful technical reading. The analytical benchmark suggested in Comment 1 proved especially useful and is now part of the paper.

### Comment 3.1: Broadcast return map and heterogeneous dynamics

> The treatment of the broadcast return map in Section 3.4 requires further justification. Equation (6) implies a constant normalized gain, whereas the heterogeneous update rule in Eq. (3) does not generally lead to this result. [...] The authors should reconcile the constant-gain interpretation with the heterogeneous dynamics and clarify whether Eq. (6) is intended as an approximation over a particular parameter or fitting range. Examining the normalized increments or fitting residuals would help establish this point.

We thank the reviewer for this point and the benchmark. The original Eq. (6) is now Eq. (7). Following the reviewer's suggestion, we reconcile the two in Section 3.4:

- **Homogeneous case.** For δλ = 0 the constant gain is exact: every unlearned agent transitions with p = ρ₀λ₀, so 1 − f(t+1) = (1 − p)[1 − f(t)], giving a₀ = ρ₀λ₀. Our fits recover this value.
- **Heterogeneous case.** For δλ > 0, the reviewer's F(t) gives the exact gain at step t (new Eq. 8). It is the mean of p₊ and p₋ over the agents still unlearned, so it decreases from ρ₀λ₀ to ρ₀(λ₀ − δλ) as the fast agents are depleted. Its initial and final values differ by the factor λ₀/(λ₀ − δλ), independent of ρ₀.
- **Normalized increments and residuals.** New Fig. 8 compares the per-step gain measured in the simulations, normalized by ρ₀λ₀, with Eq. (8), and shows the residuals. The two agree to within 0.022 ρ₀λ₀ at every f.
- **Range of validity.** New Table 4 lists how far the fitted constant departs from the exact gain. The constant gain is exact at δλ = 0, a fair approximation for δλ ≤ 0.2 (largest relative difference about 0.36), and only an effective average for stronger heterogeneity. We state this in Section 3.4 and in the Table 3 caption.
- **Mechanism.** New Fig. 4 shows the learning of the fast and slow agents separately. In the broadcast regime they learn on separate timescales, which is the depletion described by Eq. (8). In the nearest-neighbor regime the two groups learn nearly together, because a slow agent is usually surrounded by learned neighbors by the time the front reaches it (Section 3.2).

We also note that the heterogeneous broadcast gain decreases with f, opposite to the trend in the class-level empirical data.

### Comment 3.2: System-size attribution in Table 3 and Fig. 6

> The interpretation of the fitted coefficients in Table 3 as evidence of a system-size effect is not supported by the comparison shown in Fig. 6. [...] This attribution should be revised, or supported by a comparison at fixed coupling strength and heterogeneity. [...]

We agree that the original comparison varied ρ₀ and δλ together with L. The return-map figure (Fig. 6 in the original, now Fig. 7) now shows L = 32 and L = 128 at the same ρ₀ and δλ in each row, for an unfavorable (ρ₀ = 0.1, δλ = 0.4) and a favorable (ρ₀ = 0.9, δλ = 0) parameter set, so panels in a row differ only in size. Table 3 is rebuilt accordingly. With this comparison:

- The broadcast gain does not depend on L. The earlier statement that it does came from the change in ρ₀ and δλ and has been removed.
- The nearest-neighbor coefficients do decrease with L at fixed parameters, with |a₂| ∝ L^−β and β between 0.5 and 0.9 for every combination of ρ₀ ∈ {0.1, 0.5, 0.9} and δλ ∈ {0, 0.2, 0.4} (new Appendix B, Fig. 11). The same holds at a fixed initial learned fraction (exponents change by less than 0.01).

The coefficients are now fitted separately for each run (see Reviewer 1, comment 1.7).

### Comment 3.3: Stopping criterion and four-stage behavior

> The influence of the stopping criterion on the reported relaxation times needs clarification. [...] Please specify this treatment and check the sensitivity of the relevant results to the stopping window for representative cases in the affected parameter range. The reported four-stage behavior should also be examined in relation to the analytical broadcast curve and the logarithmic time axis. [...]

**Stopping criterion.** We have removed it. Every run now continues until all agents have learned, and all 100,800 runs reached this absorbing state. There are therefore no terminated runs to treat separately, and t_max is the absorption time (Section 2). The former Tables 1 and 2 are removed.

**Four-stage behavior.** As the reviewer suggested, this is the two-rate relaxation described by F(t), not a separate mechanism or an artifact of termination. The mean of the simulated broadcast runs agrees with F(t) to within 1.5 × 10⁻³ at every time step. On the logarithmic time axis, the curve rises as the high-λ subpopulation learns, plateaus near f = 1/2, rises again as the low-λ subpopulation learns, and then saturates (Section 3.2). We now describe it in these terms rather than as a behavior that is not fully understood. The effect of the logarithmic axis on how the stages appear is noted in Section 3.1.

All figures and values were regenerated from the new runs (see also Reviewer 1, comment 1.4).

### Comment 3.4: Initialization in Fig. 1 and model specification

> There is an inconsistency between the initialization described in the text and that shown in Fig. 1. [...] The model description should also specify the proportions and spatial assignment of the two learning-rate classes, whether learning rates remain fixed during each run, and whether all states are updated synchronously.

**Fig. 1.** The flowchart now shows four initially learned agents, which is the initialization in both regimes in the revised simulations (see Reviewer 1, comment 1.1). It also ends when all agents have learned. The caption states which equation gives the transition probability in each regime and that updates are synchronous.

**Model specification.** Section 2 now states the following:

- Each agent is assigned to the high or low learning-rate class independently with probability ½, so the assignment is spatially random and the classes are equal in size only on average.
- The learning rates are fixed during a run and redrawn for each run.
- All states are updated synchronously.
- Boundaries are open, and the seating rules are given in Table 1.

### Comment 3.5: Scope of conclusions

> Some conclusions extend beyond what is established by the present analysis. The quadratic PI return map is fitted to the simulated trajectories [...] Similarly, the separate simulations of broadcast and nearest-neighbor interactions do not establish that a combined protocol is optimal. [...] The wording concerning empirical validation should also reflect the discrepancies already acknowledged in Section 3.4.

- **Quadratic return map.** We now describe the nearest-neighbor return map as well described by a fitted quadratic. Section 3.4 states that, unlike the broadcast form at δλ = 0, it is not derived from the update rule: the next learned fraction depends on the positions of the learned agents, not on f alone. It is therefore a fitted description of the simulated return maps over the range of Fig. 7 (δλ ≤ 0.4), used to summarize the shape of the gain. "Analytic" is removed from the highlights and the Conclusion.
- **Combined protocol.** This is now presented as a hypothesis for future investigation, not a conclusion (Sections 3.2, 3.3, and 4).
- **Empirical wording.** "Validate" and "qualitative agreement" are replaced by specific statements of what the model does and does not reproduce (see Reviewer 1, comment 1.6, and Reviewer 2, comment 2.5).

We thank the reviewer again; addressing these points has made the analysis considerably more consistent.
