# Run from the repository root (paths are relative to it); see README.md.
#=
2DBPCA-IH v2 sweep: runs every row of parameter_combinations.csv × N_TRIALS to absorption.

Outputs (under output/2D-Binary-PCA-IH-v2/):
  * raw/<param_key>.parquet — one file per parameter row; one row per agent (column-major
    linear index), columns T_<trial> (learning time; 0 = seed, -1 = never learned) and
    class_<trial> (1 = fast λ₀+δλ, 2 = slow λ₀-δλ).
  * analysis/run_summary.csv — one row per run: parameters, trial, seed, absorbed flag,
    t_max, threshold times t_50/t_90/t_95, and runtime_s (wall time of the simulation itself,
    excluding file I/O; NaN for rows reloaded from an earlier sweep).
  * analysis/run_log.txt — appended per sweep: start/end timestamps, elapsed time, threads,
    Julia version, git commit, rows simulated vs. reloaded.

Rows whose raw file already exists are skipped unless OVERWRITE = true, so an interrupted
sweep resumes where it stopped. Needs a multithreaded session (julia -t auto, or the VS Code
julia.NumThreads setting).
=#

begin
    using Pkg
    Pkg.activate(".")
    include("2DBPCA-IH-common-v2.jl")
    using ProgressMeter, Dates
end

const OVERWRITE = false

#* One row of run_summary.csv.
function summary_row(row, trial, T, absorbed, runtime_s)
    return (;
        regime = row.regime, SA = row.SA, L = row.L, n_learned = row.n_learned,
        ρ₀ = row.ρ₀, λ₀ = row.λ₀, δλ = row.δλ, trial, seed = row[seed_col(trial)], absorbed,
        t_max = absorbed ? Int(maximum(T)) : -1,
        t_50 = time_to_fraction(T, 0.50),
        t_90 = time_to_fraction(T, 0.90),
        t_95 = time_to_fraction(T, 0.95),
        runtime_s,
    )
end

#* Runs all trials of one parameter row and writes its raw file. Returns its summary rows.
function run_parameter_row(row)
    raw = DataFrame()
    summary = NamedTuple[]
    for trial in 1:N_TRIALS
        seed = row[seed_col(trial)]
        runtime_s = @elapsed (T, λ_class, absorbed) = simulate_run(row.regime, row.SA, row.L,
            row.n_learned, row.ρ₀, row.λ₀, row.δλ, seed)
        raw[!, "T_$(trial)"] = vec(T)
        raw[!, "class_$(trial)"] = vec(λ_class)
        push!(summary, summary_row(row, trial, T, absorbed, runtime_s))
    end
    write_parquet(raw, raw_path(row))
    return summary
end

#* Recomputes a row's summary from its existing raw file (used for skipped rows).
function summarize_existing_row(row)
    summary = NamedTuple[]
    for (trial, (T, _)) in enumerate(read_trials(row))
        push!(summary, summary_row(row, trial, T, all(>=(0), T), NaN))
    end
    return summary
end

const RUN_LOG_PATH = "$(OUTPUT_DIR)/analysis/run_log.txt"

#* Git commit of the code, for the run log ("unknown" outside a git repository).
git_commit() = try strip(read(`git rev-parse --short HEAD`, String)) catch; "unknown" end
scripts_modified() = try !isempty(strip(read(`git status --porcelain -- scripts`, String))) catch; false end

function run_sweep(param_sets; overwrite = OVERWRITE)
    t_start = now()
    commit = git_commit()
    dirty = scripts_modified()
    println("Sweep started $(t_start) — $(nrow(param_sets)) rows × $(N_TRIALS) trials, $(Threads.nthreads()) threads")

    results = Vector{Vector{NamedTuple}}(undef, nrow(param_sets))
    simulated = fill(false, nrow(param_sets)) #* Vector{Bool}, not BitVector: written from threads
    prog = Progress(nrow(param_sets); showspeed = true)
    Threads.@threads for i in 1:nrow(param_sets)
        row = param_sets[i, :]
        if !overwrite && isfile(raw_path(row))
            results[i] = summarize_existing_row(row)
        else
            results[i] = run_parameter_row(row)
            simulated[i] = true
        end
        next!(prog)
    end
    summary_df = DataFrame(reduce(vcat, results))
    mkpath(dirname(SUMMARY_PATH))
    CSV.write(SUMMARY_PATH, summary_df)
    t_end = now()

    n_capped = count(!, summary_df.absorbed)
    elapsed = Dates.canonicalize(Dates.CompoundPeriod(t_end - t_start))
    sim_total = sum(filter(!isnan, summary_df.runtime_s); init = 0.0)
    log_text = """
    ---
    sweep start   $(t_start)
    sweep end     $(t_end)
    elapsed       $(elapsed)
    threads       $(Threads.nthreads())
    julia         $(VERSION)
    git commit    $(commit)$(dirty ? " (scripts/ has uncommitted changes)" : "")
    rows          $(nrow(param_sets)) ($(count(simulated)) simulated, $(count(!, simulated)) reloaded from existing raw files)
    runs          $(nrow(summary_df)); hit MAX_STEPS: $(n_capped)
    sim time      $(round(sim_total; digits = 1)) s summed over simulated runs (excludes file I/O)
    """
    open(io -> write(io, log_text), RUN_LOG_PATH, "a")
    print(log_text)
    return summary_df
end

#! Main
begin
    Threads.nthreads() == 1 && @warn "Running single-threaded; start Julia with -t auto for the full sweep"
    param_sets = read_parameter_combinations()
    summary_df = run_sweep(param_sets)
end
