# ! USE BS THESIS ENVIRONMENT (file paths will break if not)
#=
Physica A revision, reviewer comment 7: report fitted equations, standard errors,
R^2, and N (number of return-map points) for Figure 6 (return map), and make
explicit how f=1 points are handled.

This script does NOT modify scripts/2DBPCA-IH-trad-PI-comparison.jl. It re-reads
the same raw simulation output and the same Roxas experiment-data files that
script uses, and refits with two changes:

1. For the TI (traditional/broadcast) and PI (inner_corner/nearest-neighbor)
   SIMULATION fits, h(1)=1 is imposed as a hard constraint rather than left
   free. This is not just a fitting convenience: f=1 (all agents learned) is
   an exact absorbing state of the model (see generate_next_generation in
   2DBPCA-IH-simulation.jl -- a learned agent can never become unlearned), so
   h(1)=1 must hold identically for simulation data. Enforcing it removes the
   0/0 ambiguity in the analytic gain g(f) = (h(f)-f)/(1-f) at f=1 (the limit
   is finite and equals the derivative-based expressions already used in the
   paper: g=a0 for TI, g=a0-a2*f for PI).
   - TI: h(f) = f + a0*(1-f)                          (1 free parameter: a0=m)
   - PI: h(f) = f + (1-f)*(a0 - a2*f)                  (2 free parameters: a0, a2)
2. For the EMPIRICAL (Roxas) fit, h(1)=1 is NOT imposed, because real
   pre/post assessment pairs can and do show S1=1 with S2<1 (score decreases
   are present in the dataset, per the paper's own discussion) -- there is no
   structural absorbing-state guarantee for real classrooms. The empirical
   fit therefore stays a free 3-parameter quadratic, exactly as in
   2DBPCA-IH-trad-PI-comparison.jl, and its derived gain g(f) is not
   evaluated/plotted at f=1 for this reason.

Output: a single CSV of fit coefficients + standard errors + R^2 + N, one row
per (dataset, regime) pair, written to
./output/2D-Binary-PCA-IH/analysis/return-map-stats.csv
=#

begin
    using Pkg
    Pkg.activate(".")

    using CSV
    using DataFrames
    using LsqFit
    using DelimitedFiles
    using Statistics
end

function param_in_fail_df(fail_df, params)
    seat_config, ρ₀, class_size, trial, δλ = params
    any(row -> row.seat_config == seat_config &&
                row.ρ₀ == ρ₀ &&
                row.class_size == class_size &&
                row.trial == trial &&
                row.δλ == δλ, eachrow(fail_df))
end

function read_time_series_data_raw(class_size, SA, ρ₀, λ₀, δλ; n_learned=4)
    learned_raw = []
    fail_df = CSV.read("./output/2D-Binary-PCA-IH/analysis/failures.csv", DataFrame)

    for trial in 1:20
        if param_in_fail_df(fail_df, [SA, ρ₀, class_size, trial, δλ])
            continue
        end

        if SA == "random"
            path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
        else
            path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"
        end

        params_df = CSV.read(path * "fit_params.csv", DataFrame)
        learned = params_df[!, "learned_per_gen"]
        push!(learned_raw, learned)
    end

    return learned_raw
end

function return_map_pairs(class_size, SA, ρ₀, λ₀, δλ)
    learned_data = read_time_series_data_raw(class_size, SA, ρ₀, λ₀, δλ)

    f_t1 = Float64[]
    f_t = Float64[]
    for series in learned_data
        append!(f_t1, series[1:end-1])
        append!(f_t, series[2:end])
    end
    return f_t1, f_t
end

#* R^2 relative to the h(f)=f identity is not meaningful here; use the usual
#* definition relative to the mean of the fitted (dependent) variable.
function r_squared(y_data, y_fit)
    ȳ = mean(y_data)
    SS_tot = sum((y_data .- ȳ) .^ 2)
    SS_res = sum((y_data .- y_fit) .^ 2)
    return 1 - SS_res / SS_tot
end

#* TI: h(f) = f + a0*(1-f)  -- forces h(1)=1 by construction
ti_model(f, p) = f .+ p[1] .* (1 .- f)

#* PI: h(f) = f + (1-f)*(a0 - a2*f) -- forces h(1)=1 by construction
pi_model(f, p) = f .+ (1 .- f) .* (p[1] .- p[2] .* f)

#* Empirical: free quadratic, no h(1)=1 constraint (see header note)
empirical_model(f, p) = p[1] .* f .^ 2 .+ p[2] .* f .+ p[3]

function fit_row(label, model, f_t1, f_t, p0)
    fit = curve_fit(model, f_t1, f_t, p0)
    se = try
        stderror(fit)
    catch
        fill(NaN, length(p0))
    end
    r2 = r_squared(f_t, model(f_t1, fit.param))
    return (
        dataset = label,
        n_points = length(f_t1),
        params = fit.param,
        stderrors = se,
        r_squared = r2,
    )
end

begin
    configs = [
        ("large (L=128, ρ0=0.1, δλ=0.4)", 128, 0.1, 0.5, 0.4),
        ("small (L=32, ρ0=0.9, δλ=0.0)", 32, 0.9, 0.5, 0.0),
    ]

    rows = NamedTuple[]

    for (label, L, ρ₀, λ₀, δλ) in configs
        ti_f_t1, ti_f_t = return_map_pairs(L, "traditional", ρ₀, λ₀, δλ)
        pi_f_t1, pi_f_t = return_map_pairs(L, "inner_corner", ρ₀, λ₀, δλ)

        push!(rows, fit_row("TI, " * label, ti_model, ti_f_t1, ti_f_t, [0.5]))
        push!(rows, fit_row("PI, " * label, pi_model, pi_f_t1, pi_f_t, [0.0, 1.0]))
    end

    #* Empirical (Roxas) data -- unconstrained quadratic, see header note
    input_data = Float64[]
    output_data = Float64[]
    for i in 1:14
        append!(input_data, readdlm("./output/2D-Binary-PCA-IH/analysis/PI experiment data/in_8_$(i).txt", '\n')[:, 1])
        append!(output_data, readdlm("./output/2D-Binary-PCA-IH/analysis/PI experiment data/out_8_$(i).txt", '\n')[:, 1])
    end
    push!(rows, fit_row("Empirical (Roxas, 14 classrooms)", empirical_model, input_data, output_data, [0.0, 1.0, 0.0]))

    #* Flatten into a CSV: one row per fit, with padded coefficient/SE columns
    max_p = maximum(length(r.params) for r in rows)

    out_df = DataFrame(
        dataset = String[],
        n_points = Int[],
        r_squared = Float64[],
    )
    for i in 1:max_p
        out_df[!, "coef_$(i)"] = Union{Missing,Float64}[]
        out_df[!, "stderr_$(i)"] = Union{Missing,Float64}[]
    end

    for r in rows
        row = Dict{String,Any}(
            "dataset" => r.dataset,
            "n_points" => r.n_points,
            "r_squared" => r.r_squared,
        )
        for i in 1:max_p
            row["coef_$(i)"] = i <= length(r.params) ? r.params[i] : missing
            row["stderr_$(i)"] = i <= length(r.stderrors) ? r.stderrors[i] : missing
        end
        push!(out_df, row)
    end

    mkpath("./output/2D-Binary-PCA-IH/analysis")
    CSV.write("./output/2D-Binary-PCA-IH/analysis/return-map-stats.csv", out_df)

    println(out_df)
end
