# ! USE BS THESIS ENVIRONMENT (file paths will break if not)
# ! Run multithreaded, e.g.:  julia -t auto scripts/2DBPCA-IH-absorption-check.jl
#=
Diagnostic: for the inhomogeneous (2D-Binary-PCA-IH) simulations, count how many
runs terminated by reaching full absorption (every agent learned, final learned
count == N = L^2) versus how many were terminated by the stopping/censoring rule
(no transitions for the tolerance window, but some agents still unlearned).

Fast approach:
  Each ...-data.csv stores agent states with one ROW per agent and one COLUMN per
  generation; simulate_steady_state() trims the constant tail, so the LAST column
  is the terminal configuration. Because learning is monotonic, an agent is
  learned at the end iff the last field of its row is '1'. The last field of a row
  is simply the byte immediately before that row's newline, so we can count learned
  agents with a single raw byte scan per file --- no CSV parsing. Files are read in
  parallel across threads (this is pure I/O plus a trivial byte loop, so it scales
  well with thread count on an SSD).

This script does not modify any existing script or data. Output:
  ./output/2D-Binary-PCA-IH/analysis/absorption-summary.csv   (per-parameter-set counts)
plus a printed summary broken down by regime and by heterogeneity delta-lambda.
=#

begin
    using Pkg
    Pkg.activate(".")

    using CSV
    using DataFrames
    using Printf
    using Base.Threads
end

#* Count learned agents in the terminal (last) generation by scanning raw bytes.
#* Returns (n_learned_final, n_agents). The header line is skipped; '\r\n' and a
#* possibly-missing trailing newline are both handled.
function final_learned_bytes(path)
    bytes = read(path)
    n = length(bytes)
    #* skip header line
    i = 1
    @inbounds while i <= n && bytes[i] != 0x0a
        i += 1
    end
    learned = 0
    rows = 0
    @inbounds for j in (i + 1):n
        if bytes[j] == 0x0a
            k = j - 1
            if k >= 1 && bytes[k] == 0x0d      # handle CRLF
                k -= 1
            end
            if k >= i + 1
                bytes[k] == 0x31 && (learned += 1)   # '1'
                rows += 1
            end
        end
    end
    #* final row without a trailing newline
    if n >= i + 1 && bytes[n] != 0x0a
        k = n
        bytes[k] == 0x0d && (k -= 1)
        bytes[k] == 0x31 && (learned += 1)
        rows += 1
    end
    return learned, rows
end

function data_path(SA, class_size, ρ₀, λ₀, δλ, trial; n_learned=4)
    if SA == "random"
        return "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)-$(n_learned)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(n_learned)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-data.csv"
    else
        return "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-data.csv"
    end
end

function run_absorption_check()
    #* Parameter grid (mirrors the simulation / analysis scripts)
    sizes = [32, 48, 64, 96, 128]
    seat_configs = ["traditional", "inner_corner", "outer_corner", "center", "random"]
    Ρs = collect(0.1:0.1:1.0)
    δλs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49]
    n_trials = 20
    λ₀ = 0.5
    n_learned = 4

    #* Enumerate every existing run first, so the heavy work is a flat parallel loop.
    meta = NamedTuple[]
    paths = String[]
    for SA in seat_configs, class_size in sizes, ρ₀ in Ρs, δλ in δλs, trial in 1:n_trials
        p = data_path(SA, class_size, ρ₀, λ₀, δλ, trial; n_learned=n_learned)
        if isfile(p)
            push!(meta, (; SA, class_size, ρ₀, δλ))
            push!(paths, p)
        end
    end

    nfiles = length(paths)
    println("Found $(nfiles) run files; scanning with $(nthreads()) thread(s)...")

    learned_final = Vector{Int}(undef, nfiles)
    agents = Vector{Int}(undef, nfiles)
    done = Threads.Atomic{Int}(0)

    t0 = time()
    @threads for idx in 1:nfiles
        l, r = final_learned_bytes(paths[idx])
        learned_final[idx] = l
        agents[idx] = r
        d = Threads.atomic_add!(done, 1) + 1
        if d % 5000 == 0
            @printf("  ...%d / %d files (%.0fs)\n", d, nfiles, time() - t0)
        end
    end
    @printf("Scan complete in %.1fs\n", time() - t0)

    #* Aggregate per (SA, size, ρ, δλ)
    absorbed = learned_final .== agents
    scan = DataFrame(
        seat_config = [m.SA for m in meta],
        class_size = [m.class_size for m in meta],
        ρ₀ = [m.ρ₀ for m in meta],
        δλ = [m.δλ for m in meta],
        absorbed = absorbed,
    )

    summary = combine(groupby(scan, [:seat_config, :class_size, :ρ₀, :δλ]),
        :absorbed => count => :n_absorbed,
        :absorbed => length => :n_found,
    )
    summary.n_stopped = summary.n_found .- summary.n_absorbed
    sort!(summary, [:seat_config, :class_size, :ρ₀, :δλ])

    mkpath("./output/2D-Binary-PCA-IH/analysis")
    CSV.write("./output/2D-Binary-PCA-IH/analysis/absorption-summary.csv", summary)

    #* ---- Printed summary ----
    total_found = nrow(scan)
    total_absorbed = count(absorbed)
    total_stopped = total_found - total_absorbed

    println("\n=== Absorption vs. stopping-rule summary (2D-Binary-PCA-IH) ===")
    @printf("Runs scanned: %d\n", total_found)
    @printf("Reached full absorption: %d (%.1f%%)\n", total_absorbed, 100 * total_absorbed / total_found)
    @printf("Terminated by stopping rule: %d (%.1f%%)\n", total_stopped, 100 * total_stopped / total_found)

    println("\n-- By instruction regime --")
    for (label, mask) in [
        ("Broadcast (traditional)", scan.seat_config .== "traditional"),
        ("Nearest-neighbor (PI, all SAs)", scan.seat_config .!= "traditional"),
    ]
        f = count(mask); a = count(absorbed .& mask); s = f - a
        @printf("%-32s found=%6d  absorbed=%6d (%5.1f%%)  stopped=%6d (%5.1f%%)\n",
            label, f, a, 100a / f, s, 100s / f)
    end

    println("\n-- Stopping-rule rate, regime x delta-lambda (% stopped) --")
    @printf("%-8s%12s%12s\n", "δλ", "Broadcast", "Near-neigh")
    δλs_sorted = sort(unique(scan.δλ))
    for δλ in δλs_sorted
        tb = (scan.δλ .== δλ) .& (scan.seat_config .== "traditional")
        tp = (scan.δλ .== δλ) .& (scan.seat_config .!= "traditional")
        fb = count(tb); sb = fb - count(absorbed .& tb)
        fp = count(tp); sp = fp - count(absorbed .& tp)
        (fb == 0 && fp == 0) && continue
        @printf("%-8s%11s %11s\n", string(δλ),
            fb == 0 ? "-" : @sprintf("%.1f%%", 100sb / fb),
            fp == 0 ? "-" : @sprintf("%.1f%%", 100sp / fp))
    end

    return summary
end

run_absorption_check()
