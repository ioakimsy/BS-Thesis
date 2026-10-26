Concern 1: Stopping rule makes t_max a lower bound (all)

Response: We have removed the stopping rule and now run all simulations to full absorption. Dropped table 1 and 2 that detailed the percentage of runs terminated by the stopping rule. Simulations were all reran and all figures and values were updated accordingly. The relative timing between parameter values and between TI and PI, which was our main concern, remains the same regardless of the t_max values: the faster regime changed in only 5 of 350 parameter combinations, all near-ties. While rerunning, we also found and corrected an implementation error in which the nearest-neighbor probability used the neighbors' rather than the learner's λ (Eq. 4 in the manuscript was correct). Results at δλ = 0 are unaffected and the TI–PI comparison is unchanged. At the strongest heterogeneity (δλ ≥ 0.45), the corrected PI relaxation keeps its shape but has a slower final approach to full learning, which we now state (abstract, §3.2, conclusion).

Concern 2: "Optimal mixed strategy" not shown (all)

Response: now presented as a hypothesis for future work instead of a natural conclusion or assertion.

Concern 3: Empirical / Hake agreement overstated; what does one step correspond to? (all)

Response: added a class-average version of the empirical plot (new Fig. 9), since the model's f is per class, not per student. per student, a2 is positive, but at the class level it's negative (same as PI), and class gain increases with pre-test score. we call this consistent with PI, not a confirmation. class 9 was dropped as an outlier.

for the time step question: one time step is one learning step: the interval over which a student learns with probability ρ₀λ (per learned neighbor in PI; Eqs. 3 and 4). a step thus has no fixed duration in real time, and the number of steps between two assessments is not fixed by the model. we added Eq. 10 (k-step gain), showing that the shape of the gain doesn't depend on the number of time steps between assessments, only its size. since we only compare shapes, this choice does not affect the comparison. (stated in Methods, §2, and §3.4)

also replaced "qualitative agreement" with specific claims, and listed what the model does and doesn't reproduce.

Concern 4: Return map mixes size with other parameters (R1, R3)

response: We agree the original comparison mixed L with ρ₀ and δλ (the two cases differed in all three), so the claim that the TI gain decreases with L was wrong. it came from the change in ρ₀ and δλ, and it has been removed. Fig. 7 now compares L = 32 and 128 at the same ρ₀ and δλ (each row differs only in size). Together with Table 3 (return-map fits) and Appendix B (coefficients vs L), this shows the TI gain doesn't depend on L, while the PI coefficients decrease with L at fixed parameters, and the same holds with a fixed seed fraction.

Concern 5: Missing implementation details and code (R1, R3)

Response: Added into the paper. Synchronous update, two \lambda levels are assigned randomly with probability 1/2, spatially random, fixed within a run, and redrawn for each run. Added seat coordinates/rules in table 1. Stated open boundary. Code availability edit to follow. Seeds are recorded for reproducibility.

Concern 6: Fig. 1 flowchart inconsistent (R1, R3)

Response: figure updated to n_0 = 4 to match the simulation. Changed end condition to reflect removal of stopping condition. Caption now states the synchronous update and which equation gives P for each regime.

Concern 7: "Percolation-driven" unsupported (R1, R2)

Response: changed all mentions of percolation to local propagation (keyword: local propagation dynamics), defined once in §3.1 as learning spreading outward from the initially learned agents through contact between neighbors.

Concern 8: Real-classroom framing: "small class", "safe choice" (R1, R2)

Response: wording is now lattice side length consistently, dropped the use of small/large N in favor of smaller or larger lattice side length. PI being a safe choice has been removed, we now show the absolute and relative differences in PI vs TI and emphasize its asymmetry and symmetry respectively. Added in the methodology that these are simplified models and that the results don't necessarily reflect teaching effectiveness but model dynamics. Additionally we show time to 95% learned (t95), which shows that which regime is faster depends on the measure: TI gets 95% of agents there sooner on large lattices, while PI leaves fewer agents behind. This replaces the "safe choice" claim (Fig. 6(c,d), Eq. 5, Table 2, Appendix A).

Concern 9: TI vs PI not a controlled comparison (R1)

Response: TI vs PI implementations are now twinned such that the TI now begins with 4 learned students like PI, that \lambda assignments are the same, and the same seeds are used. The only difference should be the update rules.

Concern 10: Seed density falls as L grows (R1)

Response: we resimulated for L ≥ 48 with the initial learned fraction fixed at that of L = 32 (9, 16, 36, 64 learned students). for the fixed-position clusters and TI, this did not change the size dependence of tmax on L (scaling exponent γ changed by −0.04 to 0.00 on average); for the clusters, the learned front still has to cross a distance that grows with L. however, the random seating arrangement for PI showed a big difference (γ fell from 0.91 to 0.18 at ρ₀ = 0.5, δλ = 0.2), so most of its size dependence came from the falling seed density. (Methods, §3.1, Appendix C with Table 5)

Concern 11: TI constant gain vs heterogeneous model (R3)

Response: we agree, the constant gain of Eq. 7 is exact only for a homogeneous class (δλ = 0), where a₀ = ρ₀λ₀; our fits recover this. for δλ > 0, we now derive the exact gain from F(t) (new Eq. 8, following reviewer suggestion): it is the mean transition probability of the agents still unlearned, so it starts at ρ₀λ₀ and decreases to ρ₀(λ₀ − δλ) as the fast learners are used up. the simulations agree with it to within 0.022 ρ₀λ₀ (new Fig. 8). the text now states that the fitted a₀ is exact at δλ = 0, a fair approximation up to δλ ≈ 0.2, and only an effective average for stronger heterogeneity (new Table 4). we also show the mechanism by tracking fast and slow learners separately (new Fig. 4): in TI they learn on separate timescales, while in PI slow learners are carried along by the front. a heterogeneous TI gain also falls with f, opposite to the class-level empirical trend. (§3.2, §3.4, Table 3 caption)

Concern 12: PI quadratic is a fit, not an evolution law (R3)

Response: changed the wording of the claim that PI doesn't yield a quadratic return map, but a quadratic fit describes it well

Concern 13: Uncertainty: min–max bands; correlated return-map points (R1)

Response: min-max bands changed to SD or SE bands (check fig or table description). To address the correlated return map points, we've instead fitted each run separately, and took the mean and standard error of the parameters. This applied to three places (table 3, figure 7, and appendix B; reasoning in 3.4)

Concern 14: Abstract claims a mean-λ dependence that wasn't studied (R1)

Response: removed from abstract and other descriptions that claim it, moved to the possible future studies direction

Concern 15: Research gap not sharp (R2)

Response: added a gap paragraph to the introduction. Gain models (Hake, Pritchard, Nitta) are class-level and don't specify who interacts. Classroom agent-based models (Bordogna & Albano 2001/2003; Ormazábal et al. 2021; Velásquez-Rojas & Laguna 2021) include heterogeneous students and compare lectures with lectures plus group work, but their groups have no spatial arrangement and group work supplements the lecture, so broadcast and local interaction are never compared as alternatives. Seating studies (Roxas 2010; Monterola 2009) relate placement to one pre/post session through a neural network, not a dynamical model. Our SPP proceeding compared broadcast and nearest-neighbor coupling, but with identical agents, and local vs mean-field spreading outside education (e.g. Fukś & Boccara 1997) also assumes identical agents. The gap: how learning-rate heterogeneity and interaction topology act together on a spatial lattice, including the effect on relaxation times, their size dependence, and how long the slowest agents are left behind. Also framed the model as bottom-up (dynamics follow from interaction rules), with the model-to-classroom mapping stated as an interpretive choice. New references: Ormazábal 2021, Velásquez-Rojas 2021, Koponen 2017, Singh 2009, Cirillo 2024, Fukś & Boccara 1997.