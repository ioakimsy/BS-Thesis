# ! USE BS THESIS ENVIRONMENT (file paths will break if not)
#=
Builds output/2D-Binary-PCA-IH-v2/parameter_combinations.csv: one row per parameter set,
with columns seed_1..seed_N_TRIALS holding each run's RNG seed.

Seeds are drawn from an unseeded source (RandomDevice), one set per twin_key: a TI row and
its PI twin share seeds, so they start from identical seed positions and λ-class grids. Seeds
are PROTECTED: rows already in the file keep their seeds (matched on twin_key), and only new
rows get fresh ones. Set REPLACE_SEEDS = true to deliberately draw a fresh set for every row.

Every downstream v2 script reads the parameter grid from this file instead of
redefining the parameter lists.
=#

begin
    using Pkg
    Pkg.activate(".")
    include("2DBPCA-IH-common-v2.jl")
end

const REPLACE_SEEDS = false

function build_parameter_grid()
    regimes = ["TI", "PI"]
    SAs = ["outer_corner", "inner_corner", "center", "random"]
    Ls = [32, 48, 64, 96, 128]
    Ρs = round.(collect(0.1:0.1:1.0); digits = 4)
    δλs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49]
    λ₀ = 0.5

    rows = NamedTuple[]
    #* Main sweep: 4 initially learned agents for every SA and L (both regimes start identically)
    for regime in regimes, SA in SAs, L in Ls, ρ₀ in Ρs, δλ in δλs
        push!(rows, (; regime, SA, L, n_learned = 4, ρ₀, λ₀, δλ))
    end
    #* Seed-density set (reviewer R1): initial learned fraction fixed at 4/32², every SA.
    #* Seeds fill outward from each SA's points of interest (see initial_seeds). L = 32 has
    #* n_learned = 4, so it is shared with the main sweep rather than duplicated.
    for regime in regimes, SA in SAs, L in Ls, ρ₀ in Ρs, δλ in δλs
        n_learned = round(Int, 4 * L^2 / 32^2)
        n_learned == 4 && continue
        push!(rows, (; regime, SA, L, n_learned, ρ₀, λ₀, δλ))
    end
    return DataFrame(rows)
end

function attach_seeds!(param_sets, path; n_trials = N_TRIALS, replace_seeds = REPLACE_SEEDS)
    seed_cols = [seed_col(t) for t in 1:n_trials]

    existing = Dict{String, Vector{Int64}}()
    if isfile(path) && !replace_seeds
        previous = read_parameter_combinations(path)
        if all(c -> String(c) in names(previous), seed_cols)
            for r in eachrow(previous)
                existing[twin_key(r)] = [r[c] for c in seed_cols]
            end
        end
    end

    rng = RandomDevice()
    seeds = Matrix{Int64}(undef, nrow(param_sets), n_trials)
    n_kept = 0
    for (i, r) in enumerate(eachrow(param_sets))
        key = twin_key(r)
        if haskey(existing, key)
            n_kept += 1
        else
            existing[key] = rand(rng, 1:typemax(Int64), n_trials) #* the twin row reuses these
        end
        seeds[i, :] = existing[key]
    end
    for (t, c) in enumerate(seed_cols)
        param_sets[!, c] = seeds[:, t]
    end
    println("Seeds: $n_kept rows reused existing or twin seeds, $(nrow(param_sets) - n_kept) rows drew fresh ones",
        replace_seeds ? " (REPLACE_SEEDS = true)" : "")
    return param_sets
end

#! Main
begin
    param_sets = build_parameter_grid()
    @assert allunique(param_key.(eachrow(param_sets)))
    attach_seeds!(param_sets, PARAMETER_COMBINATIONS_PATH)
    mkpath(dirname(PARAMETER_COMBINATIONS_PATH))
    CSV.write(PARAMETER_COMBINATIONS_PATH, param_sets)
    println("Wrote $(nrow(param_sets)) parameter rows × $(N_TRIALS) trials to $(PARAMETER_COMBINATIONS_PATH)")
end
