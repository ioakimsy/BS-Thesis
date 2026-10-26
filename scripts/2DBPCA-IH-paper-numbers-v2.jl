# ! USE BS THESIS ENVIRONMENT (file paths will break if not)
#=
2DBPCA-IH v2: summary numbers quoted in the paper that are not printed by the figure scripts.
Reads only analysis/run_summary.csv.
  * Table 5 (Appendix C): scaling exponent γ of t_max ∝ L^γ (log–log slope, L = 32 … 128) with four
    initially learned agents (γ₄) and with a fixed initial learned fraction (γ_f), and Δγ = γ_f − γ₄
    over all 70 (ρ₀, δλ) combinations.
  * Section 3.3: largest lead of each regime in t_max, in steps and as a ratio (inner corner,
    L ∈ {32, 128}, δλ ≤ 0.4, i.e. the range of Fig. 6).
Table 2 (t_max/t₉₅ ratios) is printed by t95_ratio_table() in the plotting script.
=#

begin
    using Pkg
    Pkg.activate(".")
    include("2DBPCA-IH-common-v2.jl")
    using Printf, Statistics
end

const Ls = [32, 48, 64, 96, 128]
n_frac(L) = round(Int, 4 * L^2 / 32^2)
loglog_slope(t) = (x = log.(Ls); y = log.(t); cov(x, y) / var(x))

function scaling_exponent(g, regime, SA, ρ₀, δλ; fixed_fraction)
    n_learned(L) = fixed_fraction ? n_frac(L) : 4
    t = [only(rows_where(g; regime, SA, L, n_learned = n_learned(L), ρ₀, δλ).t_max) for L in Ls]
    return loglog_slope(t)
end

function table5(summary)
    mean_tmax = combine(groupby(summary, [:regime, :SA, :L, :n_learned, :ρ₀, :δλ]), :t_max => mean => :t_max)
    ρs = round.(0.1:0.1:1.0; digits = 1); δλs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49]
    println("Table 5: γ₄ / γ_f at (ρ₀, δλ) = (0.5, 0.2) and (0.1, 0.4); Δγ mean (range) over all ρ₀, δλ")
    for (regime, SA) in (("PI", "inner_corner"), ("PI", "center"), ("PI", "outer_corner"), ("PI", "random"), ("TI", "inner_corner"))
        γ(ρ₀, δλ, f) = scaling_exponent(mean_tmax, regime, SA, ρ₀, δλ; fixed_fraction = f)
        Δ = [γ(ρ₀, δλ, true) - γ(ρ₀, δλ, false) for ρ₀ in ρs, δλ in δλs]
        @printf("  %-2s %-13s %.2f / %.2f   %.2f / %.2f   Δγ %+.2f (%+.2f to %+.2f)\n", regime, SA,
            γ(0.5, 0.2, false), γ(0.5, 0.2, true), γ(0.1, 0.4, false), γ(0.1, 0.4, true), mean(Δ), minimum(Δ), maximum(Δ))
    end
end

function asymmetry(summary)
    runs = summary[(summary.n_learned .== 4) .& (summary.SA .== "inner_corner") .& in.(summary.L, Ref((32, 128))) .& (summary.δλ .<= 0.4), :]
    mean_tmax = combine(groupby(runs, [:regime, :L, :ρ₀, :δλ]), :t_max => mean => :t_max)
    by_regime = unstack(mean_tmax, :regime, :t_max)
    ti_faster = by_regime[by_regime.TI .< by_regime.PI, :]; pi_faster = by_regime[by_regime.PI .< by_regime.TI, :]
    println("Section 3.3 (inner corner, L ∈ {32, 128}, δλ ≤ 0.4):")
    @printf("  broadcast faster in %d of %d cells, by at most %.1f steps; nearest-neighbor faster by at most %.1f steps\n",
        nrow(ti_faster), nrow(by_regime), maximum(ti_faster.PI .- ti_faster.TI), maximum(pi_faster.TI .- pi_faster.PI))
    @printf("  as ratios: nearest-neighbor at most %.2f× slower, broadcast at most %.2f× slower\n",
        maximum(ti_faster.PI ./ ti_faster.TI), maximum(pi_faster.TI ./ pi_faster.PI))
end

#! Main
begin
    summary = CSV.read(SUMMARY_PATH, DataFrame)
    table5(summary)
    asymmetry(summary)
end
