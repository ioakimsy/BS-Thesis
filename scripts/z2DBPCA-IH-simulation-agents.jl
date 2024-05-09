#* Probabilistic cellular automata model of peer to peer learning in a classroom with inhomogeneous learning rates
begin
    using Pkg

    using Agents

    using Random
    using Statistics
    using DelimitedFiles
    using CSV
    using DataFrames
    using LsqFit

    using Alert
    using ProgressMeter

    using CairoMakie

end


#* Initializing learning factor due to spatial factors
#* Use get(Ρ_dict, (i,j), 0) to get the value of ρ at (i,j) if it exists, else 0
#* By default, homogenous learning rate is assumed. If Ρ₀ is provided, then the learning rate is inhomogeneous
function initialize_spatial_factor(ρ₀::Float64, Ρ₀::Matrix{Float64}=Matrix{Float64}(undef, 0, 0))
    if isempty(Ρ₀)
        Ρ₀ = Matrix{Float64}([
            ρ₀ ρ₀ ρ₀;
            ρ₀ 1 ρ₀;
            ρ₀ ρ₀ ρ₀
            ])
        end
        
        if size(Ρ₀)[1] == size(Ρ₀)[2]
            n = size(Ρ₀)[1]
            relpos = Tuple.(CartesianIndices((-(n-1)÷2:(n-1)÷2,-(n-1)÷2:(n-1)÷2)))
            Ρ_dict =  Dict(relpos .=> Ρ₀)
    else
        error("ERROR: Ρ₀ should be a square matrix")
    end

    return Ρ_dict
end

#* Defining evolution rules
function class_step!(student, model)
    if student.learned
        return
    end
    
    Ρ_dict = model.Ρ
    
    neighbors = nearby_agents(student, model, 1) # Vector of agents
    neighbors_learned = [neighbor.learned for neighbor in neighbors] # Vector of student states
    neighbors_relpos = Tuple.([neighbor.pos .- student.pos for neighbor in neighbors]) # Vector of relative positions of neighbors
    neighbors_ρ = [get(Ρ_dict, neighbor_relpos, 0) for neighbor_relpos in neighbors_relpos] # Vector of ρ values for each neighbor
    
    #* Probability of learning
    #* Pᵢ = 1 - Πⱼ(1-λᵢ⋅ρⱼ⋅sⱼ) where i is the student of interest and j are the neighbors
    learn_prob = 1 - prod(1 .- neighbors_ρ[neighbors_learned .== true] .* student.learning_rate)
    
    if rand() < learn_prob
        student.learned = true
    end

    return
end

begin
    L = 64
    ρ₀ = 0.5
    n_learned = 4
    
    @agent struct Student(GridAgent{2})
        learned::Bool
        learning_rate::Float64
    end

    Ρ_dict = initialize_spatial_factor(ρ₀)
    properties = Dict(:Ρ => Ρ_dict)
    
    classroom = GridSpaceSingle((L,L), periodic = false, metric = :chebyshev)
    
    

    model = StandardABM(Student, classroom; agent_step! = class_step!, model_step! = dummystep, properties)

    for i in 1:L^2
        add_agent_single!(model; pos = Tuple(CartesianIndices((1:L, 1:L))[i]), learned = false, learning_rate = rand()<0.5 ? 0.25 : 0.75)
    end

    for i in 1:L^2
        if model[i].pos in list_learned
            model[i].learned = true
        end
    end

    adf, mdf = run!(model, 0; adata = [(:learned, sum)])
end

studentcolor(a) = a.learned ? :black : :white
studentmarker(a) = a.learning_rate == 0.25 ? :circle : :rect

for i in 1:20
    run!(model, 1)
    fig, _ = abmplot(model; 
        agent_color = studentcolor, 
        agent_marker = studentmarker,
        agent_markersize = 0.1,
    )
    display(fig)
end

list_learned = [(1,1), (1,L), (L,1), (L,L)]

model[64].pos in [Tuple(CartesianIndices((1:L÷2, 1:L÷2))[i]) for i in 1:L]
[model[i].pos in list_learned for i in 1:L^2]