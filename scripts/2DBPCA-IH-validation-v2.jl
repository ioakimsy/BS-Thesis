# ! USE BS THESIS ENVIRONMENT (file paths will break if not)
#=
Validation of the v2 engine (2DBPCA-IH-common-v2.jl) before the full sweep. Read-only except
for a scratch folder under output/2D-Binary-PCA-IH-v2/validation/, which it deletes afterwards.
Prints PASS/FAIL per check; run with several threads (julia -t auto) for check 3.

  1. Seed placement: n_learned = 4 reproduces v1's seed cells for every SA and L; larger
     n_learned fills outward from the SA's points of interest.
  2. TI/PI twins: same seed -> identical seed positions and λ-class grid.
  3. Reproducibility: same seed -> identical run, repeated and across threads; the Parquet
     write/read round-trip is exact; and a fingerprint of fixed runs is printed so separate
     Julia processes (or machines) can be compared.
  4. Engine vs. a literal transcription of Eq. 4: a naive full-grid implementation (product
     over the Moore neighborhood, no lookup tables or neighbor counts) consuming the RNG in
     the same order must give identical learning-time grids.
  5. TI vs. analytics: mean learned fraction vs. F(t) (reviewer R3) and mean t_max vs. the
     exact expected maximum of independent geometric learning times.
  6. PI vs. v1 data at δλ = 0, where the learner-λ fix changes nothing: t_50 and t_max must
     agree within sampling error.
  7. Worst cases absorb well below MAX_STEPS.
=#

begin
    using Pkg
    Pkg.activate(".")
    include("2DBPCA-IH-common-v2.jl")
    using Statistics, Printf
end

report(name, ok, detail = "") = (println(ok ? "PASS" : "FAIL", "  ", name, isempty(detail) ? "" : "  — " * detail); ok)

#* --- 1. Seed placement -------------------------------------------------------
function check_placement()
    lin(L, r, c) = LinearIndices((L, L))[r, c]
    v1(SA, L) = let h = L ÷ 2, q = L ÷ 4
        SA == "center"       ? [lin(L, h, h), lin(L, h + 1, h), lin(L, h, h + 1), lin(L, h + 1, h + 1)] :
        SA == "outer_corner" ? [lin(L, 1, 1), lin(L, 1, L), lin(L, L, 1), lin(L, L, L)] :
                               [lin(L, q, q), lin(L, q, L - q + 1), lin(L, L - q + 1, q), lin(L, L - q + 1, L - q + 1)]
    end
    ok = all(sort(initial_seeds(SA, L, 4, StableRNG(s))) == sort(v1(SA, L))
             for SA in ["center", "outer_corner", "inner_corner"], L in [32, 48, 64, 96, 128], s in 1:5)
    report("4 seeds reproduce v1 placement (all SA, all L)", ok)

    #* Larger n: every seed is at least as close to the points of interest as every non-seed.
    ok2 = true
    for SA in ["center", "outer_corner", "inner_corner"], (L, n) in [(48, 9), (64, 16), (96, 36), (128, 64)]
        pts = SA_points(SA, L); cart = CartesianIndices((L, L))
        d(i) = minimum(hypot(cart[i][1] - p[1], cart[i][2] - p[2]) for p in pts)
        seeds = initial_seeds(SA, L, n, StableRNG(1))
        rest = setdiff(1:L^2, seeds)
        ok2 &= length(unique(seeds)) == n && maximum(d.(seeds)) <= minimum(d.(rest)) + 1e-12
    end
    report("n > 4 seeds fill outward from SA points (nearest-first)", ok2)

    println("      seed layouts, 9×9 window around one point of interest (■ = seed):")
    for SA in ["center", "outer_corner", "inner_corner"], (L, n) in [(48, 9), (64, 16)]
        S = falses(L, L); S[initial_seeds(SA, L, n, StableRNG(2))] .= true
        p = round.(Int, SA_points(SA, L)[1]); r0 = clamp(p[1] - 4, 1, L - 8); c0 = clamp(p[2] - 4, 1, L - 8)
        rows = [join(S[r, c] ? "■" : "·" for c in c0:c0 + 8) for r in r0:r0 + 8]
        println("      $(rpad(SA, 13)) L=$L n=$n (rows $r0-$(r0 + 8), cols $c0-$(c0 + 8)), $(count(S)) seeds total")
        foreach(r -> println("        ", r), rows)
    end
    return ok && ok2
end

#* --- 2. Twins ---------------------------------------------------------------
function check_twins()
    ok = true
    for SA in ["random", "center"], seed in 1:5
        rng1 = StableRNG(seed); s1 = initial_seeds(SA, 64, 16, rng1); c1 = initial_λ_class(64, rng1)
        _, cTI, _ = simulate_run("TI", SA, 64, 16, 0.3, 0.5, 0.4, seed)
        _, cPI, _ = simulate_run("PI", SA, 64, 16, 0.3, 0.5, 0.4, seed)
        ok &= c1 == cTI == cPI
    end
    #* Seeds are recoverable from T == 0; check TI and PI share them.
    T1, _, _ = simulate_run("TI", "random", 64, 16, 0.3, 0.5, 0.4, 7)
    T2, _, _ = simulate_run("PI", "random", 64, 16, 0.3, 0.5, 0.4, 7)
    ok &= findall(==(0), vec(T1)) == findall(==(0), vec(T2))
    #* And the parameter table actually pairs them.
    params = read_parameter_combinations()
    twin_groups = groupby(params, [:SA, :L, :n_learned, :ρ₀, :λ₀, :δλ])
    ok_table = all(nrow(sg) == 2 && sg[1, :seed_1] == sg[2, :seed_1] && sg[1, :seed_20] == sg[2, :seed_20] for sg in twin_groups)
    report("TI/PI twins share seed positions and λ grid", ok) & report("parameter table pairs every TI row with its PI twin", ok_table)
end

#* --- 3. Reproducibility -------------------------------------------------------
function check_reproducibility()
    run1 = simulate_run("PI", "random", 64, 16, 0.3, 0.5, 0.4, 12345)
    run2 = simulate_run("PI", "random", 64, 16, 0.3, 0.5, 0.4, 12345)
    ok1 = run1 == run2
    runs = Vector{Any}(undef, 4 * Threads.nthreads())
    Threads.@threads for j in eachindex(runs)
        runs[j] = simulate_run(isodd(j) ? "TI" : "PI", "center", 64, 4, 0.3, 0.5, 0.4, 999)
    end
    ok2 = all(runs[j] == runs[isodd(j) ? 1 : 2] for j in eachindex(runs))

    dir = "$(OUTPUT_DIR)/validation"
    row = (; regime = "PI", SA = "random", L = 48, n_learned = 9, ρ₀ = 0.3, λ₀ = 0.5, δλ = 0.4)
    T, c, _ = simulate_run(row..., 2024)
    path = "$(dir)/roundtrip.parquet"
    write_parquet(DataFrame(T_1 = vec(T), class_1 = vec(c)), path)
    raw = read_parquet(path)
    ok3 = reshape(Vector{Int32}(raw.T_1), 48, 48) == T && reshape(Vector{Int8}(raw.class_1), 48, 48) == c
    rm(dir; recursive = true)

    #* Fingerprint: compare this line between separate runs of the script.
    fp = hash([simulate_run(r, SA, 32, 4, 0.5, 0.5, 0.45, s)[1] for r in ("TI", "PI"), SA in ("random", "inner_corner"), s in 1:3])
    println("      fingerprint (must be identical across processes/machines): ", string(fp, base = 16))
    report("same seed -> identical run", ok1) & report("identical across $(Threads.nthreads()) threads", ok2) &
        report("Parquet round-trip exact", ok3)
end

#* --- 4. Engine vs. literal Eq. 4 --------------------------------------------
#* Same RNG consumption order as simulate_run: seeds, λ classes, then per step one draw per
#* unlearned agent with P > 0, in increasing linear index.
#* Lattice sites around (r, c), excluding (r, c) itself and sites off the lattice
moore_neighbors(r, c, L) = [(r2, c2) for r2 in max(1, r - 1):min(L, r + 1), c2 in max(1, c - 1):min(L, c + 1) if (r2, c2) != (r, c)]

function naive_run(regime, SA, L, n, ρ₀, λ₀, δλ, seed)
    rng = StableRNG(seed)
    seeds = initial_seeds(SA, L, n, rng)
    λ = [c == 1 ? λ₀ + δλ : λ₀ - δλ for c in initial_λ_class(L, rng)]
    s = zeros(Int, L, L); s[seeds] .= 1
    T = fill(Int32(-1), L, L); T[seeds] .= 0
    t = 0
    while any(==(0), s)
        t += 1
        next = copy(s)
        for i in 1:L^2
            s[i] == 1 && continue
            r, c = Tuple(CartesianIndices((L, L))[i])
            P = if regime == "TI"
                ρ₀ * λ[i]
            else
                1 - prod((1 - λ[i] * ρ₀ * s[r2, c2]) for (r2, c2) in moore_neighbors(r, c, L); init = 1.0)
            end
            if P > 0 && rand(rng) < P
                next[i] = 1; T[i] = t
            end
        end
        s = next
    end
    return T
end

function check_vs_naive()
    n_same = 0; n_total = 0
    for regime in ("TI", "PI"), SA in ("random", "outer_corner", "center"), δλ in (0.0, 0.4, 0.49), seed in 1:3
        n_total += 1
        n_same += naive_run(regime, SA, 32, 4, 0.3, 0.5, δλ, seed) == simulate_run(regime, SA, 32, 4, 0.3, 0.5, δλ, seed)[1]
    end
    report("engine == literal Eq. 3/4 implementation (bit-identical T)", n_same == n_total, "$n_same / $n_total runs identical")
end

#* --- 5. TI vs. analytics ------------------------------------------------------
#* E[t_max] for independent geometric times with per-step probabilities p_i:
#* E[max] = Σ_{t≥0} (1 - Π_i (1 - (1-p_i)^t)).
function expected_tmax(n_fast, n_slow, p₊, p₋)
    E = 0.0; t = 0
    while true
        cdf = (1 - (1 - p₊)^t)^n_fast * (1 - (1 - p₋)^t)^n_slow
        E += 1 - cdf
        (1 - cdf < 1e-12 && t > 0) && break
        t += 1
    end
    return E
end

function check_TI_analytics(; n_runs = 200)
    ok = true
    for (L, ρ₀, δλ) in [(128, 0.1, 0.4), (64, 0.5, 0.0), (32, 0.9, 0.49), (96, 0.1, 0.49)]
        p₊, p₋ = ρ₀ * (0.5 + δλ), ρ₀ * (0.5 - δλ)
        runs = [simulate_run("TI", "random", L, 4, ρ₀, 0.5, δλ, s) for s in 1:n_runs]
        tmax = maximum(maximum(r[1]) for r in runs)
        meanf = mean(learned_fraction(r[1], tmax) for r in runs)
        s0 = 4 / L^2
        F(t) = s0 + (1 - s0) * (1 - 0.5 * (1 - p₊)^t - 0.5 * (1 - p₋)^t)
        err_f = maximum(abs(meanf[t + 1] - F(t)) for t in 0:tmax)
        #* exact E[t_max] given each run's own fast/slow counts among non-seed agents
        sim_tmax = [maximum(r[1]) for r in runs]
        exp_tmax = [begin
            T, c, _ = r; free = T .> 0
            expected_tmax(count(free .& (c .== 1)), count(free .& (c .== 2)), p₊, p₋)
        end for r in runs]
        se = std(sim_tmax) / sqrt(n_runs)
        z = (mean(sim_tmax) - mean(exp_tmax)) / se
        okc = err_f < 0.005 && abs(z) < 3.5
        ok &= okc
        @printf("      L=%-3d ρ₀=%.1f δλ=%-4s  max|⟨f⟩-F| = %.4f   ⟨t_max⟩ sim %.1f vs exact %.1f (z = %+.2f)\n",
            L, ρ₀, δλ, err_f, mean(sim_tmax), mean(exp_tmax), z)
    end
    report("TI matches F(t) and exact E[t_max]", ok, "$n_runs runs per case")
end

#* --- 6. PI vs. v1 data at δλ = 0 ---------------------------------------------
v1_path(SA, L, ρ, tr) = "./output/2D-Binary-PCA-IH/$SA-$L/$ρ-0.5-0.0/trial_$tr/data/2DBPCAIH-$SA-$L-$ρ-0.5-0.0-trial_$tr-fit_params.csv"

function check_vs_v1()
    ok = true
    if !isfile(v1_path("center", 64, 0.3, 1)) #* v1 output is not part of this repository
        println("SKIP  PI at δλ = 0 agrees with v1 data  — v1 output not available")
        return true
    end
    for (SA, L, ρ) in [("center", 64, 0.3), ("outer_corner", 128, 0.1), ("inner_corner", 32, 0.9), ("center", 128, 0.5), ("inner_corner", 96, 0.2)]
        isfile(v1_path(SA, L, ρ, 1)) || (println("      (v1 data missing for $SA-$L ρ₀=$ρ, skipped)"); continue)
        old = [CSV.read(v1_path(SA, L, ρ, tr), DataFrame).learned_per_gen for tr in 1:20]
        new = [learned_fraction(simulate_run("PI", SA, L, 4, ρ, 0.5, 0.0, 10_000 + s)[1]) for s in 1:100]
        t_old = [length(c) - 1 for c in old]; t_new = [length(c) - 1 for c in new]
        t50_old = [findfirst(>=(0.5), c) - 1 for c in old]; t50_new = [findfirst(>=(0.5), c) - 1 for c in new]
        se(a, b) = sqrt(var(a) / length(a) + var(b) / length(b))
        z_max = (mean(t_new) - mean(t_old)) / se(t_new, t_old)
        z_50 = (mean(t50_new) - mean(t50_old)) / max(se(t50_new, t50_old), 0.5)
        okc = abs(z_max) < 3.5 && abs(z_50) < 3.5
        ok &= okc
        @printf("      %-13s L=%-3d ρ₀=%.1f  t_max v1 %6.1f±%4.1f  v2 %6.1f±%4.1f (z %+.2f) | t_50 v1 %5.1f v2 %5.1f (z %+.2f)\n",
            SA, L, ρ, mean(t_old), std(t_old), mean(t_new), std(t_new), z_max, mean(t50_old), mean(t50_new), z_50)
    end
    report("PI at δλ = 0 agrees with v1 data", ok, "v1: 20 trials, v2: 100 runs")
end

#* --- 7. Worst cases absorb -----------------------------------------------------
function check_absorption()
    worst = Tuple{String, Int}[]
    ok = true
    for regime in ("TI", "PI"), SA in ("outer_corner", "random"), s in 1:10
        T, _, absorbed = simulate_run(regime, SA, 128, 4, 0.1, 0.5, 0.49, s)
        ok &= absorbed
        push!(worst, (regime, Int(maximum(T))))
    end
    for regime in ("TI", "PI")
        ts = [t for (r, t) in worst if r == regime]
        @printf("      %s  L=128 ρ₀=0.1 δλ=0.49: t_max %d–%d (cap %d)\n", regime, minimum(ts), maximum(ts), MAX_STEPS)
    end
    report("slowest corner of the grid absorbs every run", ok)
end

#! Main
begin
    println("Threads: ", Threads.nthreads())
    results = [
        check_placement(),
        check_twins(),
        check_reproducibility(),
        check_vs_naive(),
        check_TI_analytics(),
        check_vs_v1(),
        check_absorption(),
    ]
    println(all(results) ? "\nALL CHECKS PASSED" : "\nSOME CHECKS FAILED")
end
