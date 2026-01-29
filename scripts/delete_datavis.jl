# List of parameters
sizes = [32,48,64,96,128]
seat_configs = ["inner_corner", "outer_corner", "center", "random", "traditional"]
Λs = collect(0.1:0.1:1.0)
δλs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.49]
n_trials = 20
n_learned = 4

for size in sizes, SA in seat_configs, lambda in Λs, delta_lambda in δλs, trial in 1:n_trials
    if SA != "random"
        filepath = normpath("./output/2D-Binary-PCA-IH/$(SA)-$(size)/$(lambda)-0.5-$(delta_lambda)/trial_$(trial)/images/")
    else
        filepath = normpath("./output/2D-Binary-PCA-IH/$(SA)-$(size)-$(n_learned)/$(lambda)-0.5-$(delta_lambda)/trial_$(trial)/images/")
    end

    if isdir(filepath)
        files = readdir(filepath)
        for file in files
            rm(joinpath(filepath, file))
        end
    end
end
