begin
    println("Loading packages")
    using Pkg
    #Pkg.activate(".")

    using CSV
    using DataFrames
    using Plots
    using Statistics
    using Measurements
    using LsqFit
    using LaTeXStrings
    using Interpolations
    using Random
    using DelimitedFiles
    using Alert
    using ProgressMeter

    # using Images
    using ImageMorphology
    #Pkg.status()
    println("Done loading packages")
end

function max_distance_from_centroid(box, centroid)
    max_distance = 0.0
    centroid_x, centroid_y = centroid
    
    for idx in box
        idx = Tuple(idx)
        distance = sqrt((idx[1] + sqrt(2)/2 - centroid_x)^2 + (idx[2] + sqrt(2)/2 - centroid_y)^2)
        if distance > max_distance
            max_distance = distance
        end
    end
    
    return max_distance
end

begin
    class_size = 128
    SA = "inner_corner"
    λ₀ = 0.5
    ρ₀ = 0.1
    δλ = 0.4
    trial = 1

    path = "./output/2D-Binary-PCA-IH/$(SA)-$(class_size)/$(ρ₀)-$(λ₀)-$(δλ)/trial_$(trial)/data/2DBPCAIH-$(SA)-$(class_size)-$(ρ₀)-$(λ₀)-$(δλ)-trial_$(trial)-"

    data_df = CSV.read(path * "data.csv", DataFrame)
    num_generations = size(data_df)[2]
    generations = [data_df[!,"Generation $(i)"] for i in 1:num_generations]
    generations = [Matrix(reshape(generations[i], (class_size,class_size))) for i in 1:num_generations]

    circularities = []
    for generation in 1:length(generations)
        gen = generations[generation]
        label = ImageMorphology.label_components(gen, trues(3,3))
        if maximum(label)==4
            boxes = component_boxes(label)
            centroids = ImageMorphology.component_centroids(label)

            radii = [max_distance_from_centroid(boxes[i], centroids[i]) for i in 1:length(boxes)-1]
            areas  = π .* radii .^ 2
            pixel_count = ImageMorphology.component_lengths(label)[1:end]

            circularity = (pixel_count ./ areas)
            replace!(circularity, Inf => 1)
            push!(circularities, circularity)
        else
            break
        end
    end
end

df_cols = ["generation $(i)" for i in 1:length(circularities)]
circularities_df = DataFrame(circularities, df_cols)


max_values = []
max_indices = []
for row in eachrow(circularities_df)
    min_val = maximum(row)
    min_idx = argmax(row)
    push!(max_values, min_val)
    push!(max_indices, min_idx)
end

println("Minimum values for each row: ", max_values)
println("Column indices of minimum values: ", max_indices)

circularities_df

hcat(["cluster $(i)" for i in 1:4], String.(max_indices), max_values)

labels =  label_components(generations[2], trues(3,3))

component_boxes(labels)
component_centroids(labels)
max_distance_from_centroid(component_boxes(labels)[1], component_centroids(labels)[1])