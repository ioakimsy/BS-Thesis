# Run from the repository root (paths are relative to it); see README.md.
#=
2DBPCA-IH v2: shared definitions for the Physica A revision rerun.

Included (not run) by every v2 script, so the model, parameter schema, and file
layout are defined once:
    include("scripts/2DBPCA-IH-common-v2.jl")

Model (matches PhysicaA-paper.tex Eqs. 3-4):
  * L×L lattice, binary states, learning is irreversible.
  * Each agent's λ is λ₀+δλ ("fast", class 1) or λ₀-δλ ("slow", class 2), drawn
    independently with probability ½ per agent, fixed for the run.
  * Synchronous update: every agent's step-t probability uses the step-(t-1) state.
  * TI (broadcast):        P = ρ₀ λᵢⱼ
    PI (nearest-neighbor): P = 1 - (1 - ρ₀ λᵢⱼ)^k, k = learned agents in the Moore
                           neighborhood. Uses the LEARNER's λ; v1 used each neighbor's.
  * Open boundaries: cells off the lattice count as unlearned neighbors.
  * Runs continue until every agent has learned (absorption), or MAX_STEPS.

Per-run output is a learning-time grid T: T = 0 for seeds, T = t if the agent
learned at step t, T = -1 if it never learned (only when the run hit MAX_STEPS).
The state at step t is T .>= 0 .&& T .<= t, so the full trajectory is recoverable.
=#

begin
    using Random, StableRNGs
    using CSV, DataFrames, DuckDB
end

#* ---------------------------------------------------------------------------
#* Parameters and paths
#* ---------------------------------------------------------------------------
const MODEL_NAME = "2D-Binary-PCA-IH-v2"
const OUTPUT_DIR = "./output/$(MODEL_NAME)"
const PARAMETER_COMBINATIONS_PATH = "$(OUTPUT_DIR)/parameter_combinations.csv"
const SUMMARY_PATH = "$(OUTPUT_DIR)/analysis/run_summary.csv"

const N_TRIALS = 20
const MAX_STEPS = 1_000_000 #* safety cap only; every run should absorb well before this

#* Stringified exactly as in the output paths; keep these rounded so paths are stable.
fmt(x::Real) = string(round(Float64(x); digits = 4))

#* Key shared by the parameter table, seed protection, and the raw-output path.
param_key(regime, SA, L, n_learned, ρ₀, λ₀, δλ) =
    "$(regime)/$(SA)-$(L)-$(n_learned)/$(fmt(ρ₀))-$(fmt(λ₀))-$(fmt(δλ))"
param_key(row) = param_key(row.regime, row.SA, row.L, row.n_learned, row.ρ₀, row.λ₀, row.δλ)

#* TI/PI twins: rows that differ only in regime share seeds, so each TI run starts from exactly
#* the same seed positions and λ-class grid as its PI twin (both are drawn from the RNG before
#* any dynamics).
twin_key(row) = "$(row.SA)-$(row.L)-$(row.n_learned)/$(fmt(row.ρ₀))-$(fmt(row.λ₀))-$(fmt(row.δλ))"

#* One Parquet file per parameter row, holding all trials (columns T_<trial>, class_<trial>,
#* one row per agent in column-major linear index order).
raw_path(row) = "$(OUTPUT_DIR)/raw/$(param_key(row)).parquet"

seed_col(trial) = Symbol("seed_$(trial)")

#* A parameter row as used by raw_path / read_trial (the analysis scripts mostly use the inner corner).
param_row(regime, SA, L, ρ₀, δλ; n_learned = 4, λ₀ = 0.5) = (; regime, SA, L, n_learned, ρ₀, λ₀, δλ)

#* Rows of `df` whose columns match every keyword; floating-point values are compared with isapprox.
#* Example: rows_where(summary; regime = "TI", L = 128, ρ₀ = 0.5)
function rows_where(df; conditions...)
    keep = trues(nrow(df))
    for (col, value) in conditions
        keep .&= value isa AbstractFloat ? isapprox.(df[!, col], value) : (df[!, col] .== value)
    end
    return df[keep, :]
end

function read_parameter_combinations(path = PARAMETER_COMBINATIONS_PATH)
    #* read the bytes first: a memory-mapped file stays locked on Windows, so it could not be rewritten
    return CSV.read(read(path), DataFrame;
        types = Dict(
            :regime => String, :SA => String, :L => Int, :n_learned => Int,
            :ρ₀ => Float64, :λ₀ => Float64, :δλ => Float64,
        ))
end

#* ---------------------------------------------------------------------------
#* Initial conditions
#* ---------------------------------------------------------------------------

#* Points of interest (row, col) that seeds cluster around, per seating arrangement.
#* For n_learned = 4 the nearest cells are exactly v1's seed cells: the central 2×2 block,
#* the four corners, and the four inner-corner cells at L÷4 from each edge. (The inner-corner
#* points are lattice sites, not cell-block centers, so 4 seeds land one per corner.)
function SA_points(SA::AbstractString, L::Int)
    if SA == "center"
        return [(L / 2 + 0.5, L / 2 + 0.5)]
    elseif SA == "outer_corner"
        return [(1, 1), (1, L), (L, 1), (L, L)]
    elseif SA == "inner_corner"
        q = L ÷ 4
        return [(q, q), (q, L - q + 1), (L - q + 1, q), (L - q + 1, L - q + 1)]
    end
    error("Not a valid seating arrangement: $SA")
end

#* Linear (column-major) indices of the initially learned agents: the n_learned cells nearest
#* to the SA's points of interest, ties broken at random (same scheme as grid_order in the
#* Classroom Dynamics Research repo). "random" takes n_learned cells uniformly at random.
function initial_seeds(SA::AbstractString, L::Int, n_learned::Int, rng::AbstractRNG)
    iseven(L) || error("L must be even (got $L)")
    order = shuffle(rng, collect(1:L^2))
    if SA != "random"
        points = SA_points(SA, L)
        cart = CartesianIndices((L, L))
        dist(i) = minimum(hypot(cart[i][1] - p[1], cart[i][2] - p[2]) for p in points)
        sort!(order; by = dist) #* stable sort, so the shuffle breaks ties
    end
    return order[1:n_learned]
end

#* λ class per agent: 1 = fast (λ₀+δλ), 2 = slow (λ₀-δλ), each with probability ½.
initial_λ_class(L::Int, rng::AbstractRNG) = Int8[rand(rng) < 0.5 ? 1 : 2 for _ in 1:L, _ in 1:L]

#* ---------------------------------------------------------------------------
#* Dynamics
#* ---------------------------------------------------------------------------

#* Per-step learning probability, tabulated as ptable[λ class, k+1] for k = 0..8 learned
#* neighbors. The two regimes differ only here.
function probability_table(regime::AbstractString, ρ₀, λ₀, δλ)
    λs = (λ₀ + δλ, λ₀ - δλ)
    if regime == "TI"
        return [ρ₀ * λs[c] for c in 1:2, k in 0:8]               #* global source; k ignored
    elseif regime == "PI"
        return [1 - (1 - ρ₀ * λs[c])^k for c in 1:2, k in 0:8]   #* learner's λ, compounded over k
    end
    error("Not a valid regime: $regime")
end

#* Simulates one run to absorption. Returns (T, λ_class, absorbed).
#* RNG draws happen in a fixed order (seeds, then λ classes, then dynamics) from one stream,
#* so a run is fully determined by its seed.
function simulate_run(regime, SA, L, n_learned, ρ₀, λ₀, δλ, seed; max_steps = MAX_STEPS)
    rng = StableRNG(seed)
    seeds = initial_seeds(SA, L, n_learned, rng)
    λ_class = initial_λ_class(L, rng)
    ptable = probability_table(regime, ρ₀, λ₀, δλ)

    T = fill(Int32(-1), L, L)
    k = zeros(Int8, L, L) #* learned Moore neighbors of each agent
    cart = CartesianIndices((L, L))

    #* Marks agent i as learned at step t and adds 1 to the learned-neighbor count of each neighbor.
    function mark_learned!(i, t)
        T[i] = t
        row, col = Tuple(cart[i])
        for dc in -1:1, dr in -1:1
            (dr == 0 && dc == 0) && continue
            r, c = row + dr, col + dc
            if 1 <= r <= L && 1 <= c <= L
                k[r, c] += Int8(1)
            end
        end
    end

    foreach(i -> mark_learned!(i, 0), seeds)
    unlearned = [i for i in 1:L^2 if T[i] < 0]
    newly = Int[]

    t = 0
    while !isempty(unlearned) && t < max_steps
        t += 1
        #* Synchronous: decide every transition from the step-(t-1) state, then apply them.
        empty!(newly)
        for i in unlearned
            p = ptable[λ_class[i], k[i] + 1]
            #* agents with p = 0 draw no random number (this fixes the RNG stream, so keep the order)
            if p > 0 && rand(rng) < p
                push!(newly, i)
            end
        end
        if !isempty(newly)
            foreach(i -> mark_learned!(i, t), newly)
            filter!(i -> T[i] < 0, unlearned)
        end
    end
    return T, λ_class, isempty(unlearned)
end

#* Learned fraction at steps 0..t_end from a learning-time grid (normalized by N = L², not
#* by the final learned count as in v1).
function learned_fraction(T::AbstractMatrix, t_end = maximum(T))
    counts = zeros(Int, t_end + 1)
    for τ in T
        τ >= 0 && (counts[τ + 1] += 1)
    end
    return cumsum(counts) ./ length(T)
end

#* Learned fraction of one λ class (1 = fast, 2 = slow) at steps 0..t_end, normalized by that
#* class's size (seeds included).
function learned_fraction_by_class(T::AbstractMatrix, λ_class::AbstractMatrix, class, t_end = maximum(T))
    counts = zeros(Int, t_end + 1)
    n = 0
    for (τ, c) in zip(T, λ_class)
        c == class || continue
        n += 1
        τ >= 0 && (counts[τ + 1] += 1)
    end
    return cumsum(counts) ./ max(n, 1)
end

#* First step at which the learned fraction reaches `q` (-1 if never).
function time_to_fraction(T::AbstractMatrix, q)
    f = learned_fraction(T)
    i = findfirst(>=(q - 1e-12), f)
    return i === nothing ? -1 : i - 1
end

#* ---------------------------------------------------------------------------
#* Raw-output I/O (Parquet via DuckDB, same approach as the Classroom Dynamics Research repo)
#* ---------------------------------------------------------------------------

#* A path as a quoted SQL string literal.
sql_string(path) = "'" * replace(path, "'" => "''") * "'"

#* One in-memory DuckDB per write; .tmp then rename, so a crash never leaves a
#* truncated Parquet under the final name.
function write_parquet(df::DataFrame, path::AbstractString)
    mkpath(dirname(path))
    tmp = path * ".tmp"
    db = DuckDB.DB()
    con = DBInterface.connect(db)
    try
        DuckDB.register_data_frame(con, df, "tbl")
        DBInterface.execute(con, "COPY tbl TO $(sql_string(tmp)) (FORMAT PARQUET, COMPRESSION ZSTD)")
    finally
        DBInterface.close!(con)
        DBInterface.close!(db) #* release the file handle now (Windows can't move/delete open files)
    end
    mv(tmp, path; force = true)
    return path
end

function read_parquet(path::AbstractString)
    db = DuckDB.DB()
    con = DBInterface.connect(db)
    df = try
        DataFrame(DBInterface.execute(con, "SELECT * FROM read_parquet($(sql_string(path)))"))
    finally
        DBInterface.close!(con)
        DBInterface.close!(db)
    end
    return disallowmissing!(df; error = false)
end

#* Learning-time grid and λ-class grid of one trial from a parameter row's raw file.
function read_trial(row, trial; raw = read_parquet(raw_path(row)))
    T = reshape(Vector{Int32}(raw[!, "T_$(trial)"]), row.L, row.L)
    λ_class = reshape(Vector{Int8}(raw[!, "class_$(trial)"]), row.L, row.L)
    return T, λ_class
end

#* (T, λ_class) of every trial of a parameter row, reading its raw file once.
function read_trials(row)
    raw = read_parquet(raw_path(row))
    return [read_trial(row, trial; raw) for trial in 1:N_TRIALS]
end
