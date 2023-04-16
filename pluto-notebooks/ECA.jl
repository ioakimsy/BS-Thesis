### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5810c94d-e480-4057-8d82-5849df2b5251
begin
	using Pkg
	Pkg.activate(".")
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
end

# ╔═╡ eb5e6216-05f3-4432-ba22-bd268cc1b556
begin
	using Plots
	using PlutoUI
end

# ╔═╡ de9aceb7-7c0b-4af2-a674-4f517dba7bf0
#Conversion of rule number to state set
function rule_to_state(rule)
    answer = zeros(Int64,1,8)
    if rule>=0 && rule<=255
        for i in 7:-1:0
            if rule>=2^i
                answer[i+1]=1
                rule = rule-2^i
            end
        end
        return reverse(answer)
    end
end

# ╔═╡ 593e2c29-9873-4472-8294-64e37624bda7
function generate_next_generation(prev_generation,array_size,state_array)
    states = [
            [1,1,1],
            [1,1,0],
            [1,0,1],
            [1,0,0],
            [0,1,1],
            [0,1,0],
            [0,0,1],
            [0,0,0]
            ]
    next_generation = zeros(Int64,1,array_size)
    for i in range(2,array_size-1)
        neighborhood = prev_generation[i-1:i+1]
        for state_number in 1:1:8
            if neighborhood == states[state_number]
                next_generation[i]=state_array[state_number]
            end
        end
    end
    return next_generation
end


# ╔═╡ 1b5a3d3c-b924-41e3-94b5-fea01aab3a09
function generate_plot(grid,rule)
    max_generations, array_size = size(grid)
    state_array = rule_to_state(rule)


    current_generation_number = 1
    while current_generation_number < max_generations-1
        current_generation = grid[current_generation_number,:]
        next_generation = generate_next_generation(current_generation,array_size,state_array)
        grid[current_generation_number+1,:] = next_generation #error here, check if next generation is being produced correctly  

        current_generation_number = current_generation_number + 1

    end

    return grid, current_generation_number

end

# ╔═╡ 12d87f0c-c549-4a72-90b4-c48045774896
function wolfram_ECA(array_size,max_generations,rule)
    #Setting up first generation
    ECA_grid = zeros(Int64,max_generations,array_size)

    middle = (array_size+1)/2
    middle = convert(Int,middle)

    first_generation = zeros(Int64, 1, array_size)
    first_generation[middle] = 1

    ECA_grid[1,:] = first_generation

    ECA_grid, number_of_generated_generations = generate_plot(ECA_grid,rule)
    return heatmap(ECA_grid,
            plot_title = "Wolfram ECA Rule $(rule) with $(number_of_generated_generations+1) \ngenerations of $(array_size) length arrays",
            #plot_titlefontsize = 24,
            yflip=true,
            aspect_ratio=:equal,
            xlims=(0.5,array_size+0.5),
            showaxis = false, #Removes the axes' black lines
            axis=([], false), #Removes the numbrs on the axis
            ylims=(0.5,number_of_generated_generations+0.5),
            #ylabel = "Number of Generations",
            #dpi=150,
            size=((array_size+50)*5,(number_of_generated_generations+75)*5),
            cbar=false,
            c = palette(:grays,rev=true)
            )

end

# ╔═╡ a10479ff-765b-4296-bbdc-cb2c64f5c563
md"""
Rule number: $(@bind rule NumberField(0:255))
Maximum Generations: $(@bind max_generations NumberField(5:500))
Array Size: $(@bind array_size NumberField(1:201))
"""

# ╔═╡ 0e8ad584-bbd0-4ff8-b304-fe2a0678b64d
wolfram_ECA(array_size,max_generations,rule)

# ╔═╡ Cell order:
# ╠═5810c94d-e480-4057-8d82-5849df2b5251
# ╠═eb5e6216-05f3-4432-ba22-bd268cc1b556
# ╠═de9aceb7-7c0b-4af2-a674-4f517dba7bf0
# ╠═593e2c29-9873-4472-8294-64e37624bda7
# ╠═1b5a3d3c-b924-41e3-94b5-fea01aab3a09
# ╠═12d87f0c-c549-4a72-90b4-c48045774896
# ╠═a10479ff-765b-4296-bbdc-cb2c64f5c563
# ╠═0e8ad584-bbd0-4ff8-b304-fe2a0678b64d
