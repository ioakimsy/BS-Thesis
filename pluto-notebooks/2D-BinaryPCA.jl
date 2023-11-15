### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ bcc8a210-fa01-11ed-1db7-d7cfebd38074
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 20ca43a0-7499-4e6e-b5af-44a56b3f83c7
# ╠═╡ disabled = true
#=╠═╡
begin
	Pkg.add("Plots")
	Pkg.add("Random")
	Pkg.add("BenchmarkTools")
	Pkg.add("CurveFit")
	Pkg.add("DataFrames")
	Pkg.add("CSV")
end
  ╠═╡ =#

# ╔═╡ c54a1e59-e363-4567-b956-5b05ea26f172
begin
	using Random
	println("Random")
	using Plots
	using Plots.PlotMeasures
	println("Plots")
	using BenchmarkTools
	println("Benchmark tools")
	using CurveFit
	println("CurveFit")
	using DataFrames
	println("DataFrames")
	using DelimitedFiles
	println("DelimatedFiles")
	using CSV
	println("CSV")
end

# ╔═╡ a5ae8b1a-2249-4266-8dd3-7e6a45e637f3
md"""
# 2D Binary PCA
"""

# ╔═╡ 5798f6f7-4e75-4276-8a66-5d4518cde118
md"""
## Setting up packages
"""

# ╔═╡ 33b8010b-e92a-4a8e-a738-5441a2371119
# ╠═╡ disabled = true
#=╠═╡
Pkg.update()
  ╠═╡ =#

# ╔═╡ 520615ae-62e5-4482-8473-9027c12eed8f
default(aspect_ratio=:equal,
	cbar = false,
	showaxis = false,
	c = palette(:grays,rev=true),
	size = (512,512),
	yflip=true,
)

# ╔═╡ 8bf82d25-6643-4245-bee3-546c560708ef
function initiate_grid(type::String, L::Int=8)
    # May add line to ensure that L is even.
	#type should be ::Symbol. Symbol is :xyz
    if L%2 == 1
        L=L+1
    end

    grid = zeros(Int,L,L)

    if type == "center"
        grid[L÷2:L÷2+1,L÷2:L÷2+1] .= 1
    elseif type == "outer_corner"
        grid[1,1] = 1
        grid[1,end] = 1
        grid[end,1] = 1
        grid[end,end] = 1 
    elseif type == "inner_corner"
        grid[L÷4,L÷4] = 1
        grid[L÷4,end-L÷4+1] = 1
        grid[end-L÷4+1,L÷4] = 1
        grid[end-L÷4+1,end-L÷4+1] = 1
    else
        error("ERROR: Not a valid initital seating arrangement type")
    end

    return grid
end

# ╔═╡ caea6e08-89c4-429c-b850-3e91b026ba80
function initiate_grid_rand(num_learned::Int=4,L::Int=8)
    grid = zeros(Int,L,L)

    while sum(grid) < num_learned
        grid[rand(1:L),rand(1:L)] = 1
    end

    return grid
end

# ╔═╡ 4b621695-f682-4bb8-b9d0-17284457e040
#Note: adjacency_matrix[x,y] is not necessarily adjacency_matrix[y,x]
#Note: only works for 3x3 λ's
function generate_neighbor_list(grid,λ)
    num_students = prod(size(grid))
    adjacency_matrix = zeros(num_students,num_students)

    #Generation of adjacency matrix
    for row in 1:size(grid)[1], col in 1:size(grid)[2]
        student = grid[row,col]

        if row==1 && col==1 #TL
            adjacency_matrix[student,grid[row+1,col]] = λ[6]
            adjacency_matrix[student,grid[row,col+1]] = λ[8]
            adjacency_matrix[student,grid[row+1,col+1]] = λ[9]

        elseif row==1 && col==size(grid)[2] #TR
            adjacency_matrix[student,grid[row+1,col]] = λ[6]
            adjacency_matrix[student,grid[row,col-1]] = λ[2]
            adjacency_matrix[student,grid[row+1,col-1]] = λ[3]

        elseif row==size(grid)[1] && col==size(grid)[2] #BR
            adjacency_matrix[student,grid[row-1,col]] = λ[4]
            adjacency_matrix[student,grid[row-1,col-1]] = λ[1]
            adjacency_matrix[student,grid[row,col-1]] = λ[2]

        elseif row==size(grid)[1] && col==1 #BL
            adjacency_matrix[student,grid[row-1,col]] = λ[4]
            adjacency_matrix[student,grid[row-1,col+1]] = λ[7]
            adjacency_matrix[student,grid[row,col+1]] = λ[8]

        elseif row==1
            adjacency_matrix[student,grid[row,col-1]] = λ[2]
            adjacency_matrix[student,grid[row+1,col-1]] = λ[3]
            adjacency_matrix[student,grid[row+1,col]] = λ[6]
            adjacency_matrix[student,grid[row,col+1]] = λ[8]
            adjacency_matrix[student,grid[row+1,col+1]] = λ[9]
            
        elseif row==size(grid)[1]
            adjacency_matrix[student,grid[row-1,col-1]] = λ[1]
            adjacency_matrix[student,grid[row,col-1]] = λ[2]
            adjacency_matrix[student,grid[row-1,col]] = λ[4]
            adjacency_matrix[student,grid[row-1,col+1]] = λ[7]
            adjacency_matrix[student,grid[row,col+1]] = λ[8]
            
        elseif col==1
            adjacency_matrix[student,grid[row-1,col]] = λ[4]
            adjacency_matrix[student,grid[row+1,col]] = λ[6]
            adjacency_matrix[student,grid[row-1,col+1]] = λ[7]
            adjacency_matrix[student,grid[row,col+1]] = λ[8]
            adjacency_matrix[student,grid[row+1,col+1]] = λ[9]

        elseif col==size(grid)[2]
            adjacency_matrix[student,grid[row-1,col-1]] = λ[1]
            adjacency_matrix[student,grid[row,col-1]] = λ[2]
            adjacency_matrix[student,grid[row+1,col-1]] = λ[3]
            adjacency_matrix[student,grid[row-1,col]] = λ[4]
            adjacency_matrix[student,grid[row+1,col]] = λ[6]
            
        else
            adjacency_matrix[student,grid[row-1,col-1]] = λ[1]
            adjacency_matrix[student,grid[row,col-1]] = λ[2]
            adjacency_matrix[student,grid[row+1,col-1]] = λ[3]
            adjacency_matrix[student,grid[row-1,col]] = λ[4]
            adjacency_matrix[student,grid[row+1,col]] = λ[6]
            adjacency_matrix[student,grid[row-1,col+1]] = λ[7]
            adjacency_matrix[student,grid[row,col+1]] = λ[8]
            adjacency_matrix[student,grid[row+1,col+1]] = λ[9]
        end
    end

    return adjacency_matrix
end

# ╔═╡ 6680e9cb-40e7-47c3-8c34-cf8a7df1e323
function generate_next_generation(initial_grid::Matrix{Int},λ::Matrix{Float64})
	border_size = (size(λ)[1]-1)÷2
	
	if size(λ)[1] != size(λ)[2]
		error("λ not a square matrix")
	end
	
    grid = zeros(Int,size(initial_grid)[1]+border_size+1,size(initial_grid)[1]+border_size+1)
	
    grid[border_size+1:end-border_size,border_size+1:end-border_size] .= initial_grid

    next_gen = zeros(size(grid))
	
    for col in border_size+1:size(initial_grid)[1]+border_size, row in border_size+1:size(initial_grid)[2]+border_size

        if grid[row,col] == 1
            next_gen[row,col] = 1
            continue
        end

        neighborhood = grid[row-1:row+1, col-1:col+1]
        learn_prob = 1 - prod(1 .- (neighborhood .* λ))

        random_number = rand()

        if random_number <= learn_prob
            next_gen[row,col] = 1
        end
        
    end
    return next_gen[border_size+1:end-border_size,border_size+1:end-border_size]
	
end

# ╔═╡ 0dbd7206-cf15-40d4-b520-36455a048a6d
function generate_next_generation2(initial_grid::Array{Int},adj_matrix::Array{Float64})
	next_gen_grid = zeros(size(initial_grid))
	initial_grid_size = size(initial_grid)
	one_row_size = (1,prod(initial_grid_size))
	one_row_init_grid = collect(reshape(initial_grid,one_row_size))

	adj_dict = Dict( [(student,adj_matrix[student,:]) for student in 1:size(adj_matrix)[1]])


	for student in 1:prod(initial_grid_size)
		student_adj_mat = reshape(adj_dict[student],one_row_size) 
		# ^ reshape converts vector into matrix: use dot()
		if initial_grid[student] == 1
			next_gen_grid[student]=1
		end
		
		learn_prob = 1 - prod(1 .- (student_adj_mat .* one_row_init_grid))

		rand_num = rand()

		if rand_num <= learn_prob
			next_gen_grid[student] = 1
		end

	end

	return next_gen_grid
end

# ╔═╡ 31d1deea-ba04-4b88-92c5-14cbe4733a6d
function generate_next_generation3(initial_grid::Array{Int},Λ::Array{Float64})
	initial_grid_size = size(initial_grid)	   #remove
	one_row_size = (1,prod(initial_grid_size)) #length
	s⃗ = reshape(initial_grid,one_row_size)	#make this a column vector using vec()

	#don't hardcode
	#set vector of ones to be a variable; has to be row column
	# Λ should be a sparse matrix; use sparse(Λ)
	# Λ .* (ones(64) * s⃗')
	# map function: see file from 5-31 RM
	p = 1 .- (Λ .* (ones(64) * s⃗))

	p⃗ = 1 .- prod(p,dims=2)

	u⃗ = rand(64)

	rand_results = transpose(u⃗ .<= p⃗)

	s⃗ᵗ⁺¹ = BitMatrix(s⃗).||rand_results

	output = reshape(s⃗ᵗ⁺¹,initial_grid_size)
	#Make it so that input and output are s⃗ and not a matrix
	return output
end

# ╔═╡ 64e02837-8603-49f1-a18f-0b0391ccb430
function simulate_steady_state(seat_config, class_size, λ, steady_state_tolerance)
	# Simulate until reaching a steady state given some steady state tolerance
	initial_class = initiate_grid(seat_config, class_size)
	generations = [initial_class]
	
	steady_state = false
	num_generations = 1

	while steady_state == false
			next_gen = generate_next_generation(generations[end],λ)	
			push!(generations, next_gen)
			num_generations = num_generations + 1
		
		if generations[end] == generations[max(1,length(generations)-steady_state_tolerance)] 
			#Identifying when steady state is reached should be changed for λ << 1 because there might be generations that don't change just by chance
			steady_state = true
			#print("steady state")
		end
		
	end

	generations = generations[begin:end-steady_state_tolerance]

	num_generations = length(generations)
	
	return generations, num_generations
end

# ╔═╡ 186b2799-ede3-46d1-825c-55e67a56f4d3
function simulate(seat_config, class_size, λ, max_gen::Int)
	#Simulates the classroom until a set amount of generations
	initial_class = initiate_grid(seat_config, class_size)
	generations = [initial_class]
	
	steady_state = false
	num_generations = 1
	
	while num_generations < max_gen
		
		next_gen = generate_next_generation(generations[end],λ)
		push!(generations, next_gen)
		num_generations = num_generations + 1
		print(num_generations)
		
	end

	return generations, num_generations
end

# ╔═╡ da8c68a1-a62a-4f49-a500-e58cced203d9
function class_simulation(sizes::Vector{Int}, seat_configs::Vector{String},Λs::Vector{Float64}, steady_state_tolerance::Int)
	
	# seat config in string (eg. "outer_corners")
	# Λ is the uniform spread probability

	l = @layout [
			a{0.75h}
			b{}
		]

	for seat_config in seat_configs, λ₀ in Λs, class_size in sizes
		λ = Float64.( Matrix(
		[ 	λ₀ 		λ₀ 		λ₀;
			λ₀ 		0 		λ₀;
			λ₀ 		λ₀ 		λ₀]
		))

		generations, num_generations = simulate_steady_state(seat_config, class_size, λ, steady_state_tolerance)

		# Saving raw data
		df_cols = ["Generation $(i)" for i in 1:length(generations)]
		df_data = vec.(generations)
		
		# df_data -> student number -> student state at generation
		# state per student over time
		#df_data = [[df_data[g][s] for g in 1:length(df_data)] for s in 1:class_size^2]
	
		# row = ith student; column = jth generation
		student_states_df = DataFrame(df_data,df_cols)
		CSV.write("./../output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/data/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-data.csv",student_states_df)

		learned = map(x->sum(x), generations)

		#Set up in case need to truncate outliers
		learned_y = learned[1:end]
		generation_domain = 1:length(learned_y)

		#axᵇ
		power_coeffs = power_fit(generation_domain, learned_y)
		power_vals = power_coeffs[1] .* generation_domain .^ power_coeffs[2]

		#aeᵇˣ
		#exp_coeffs = exp_fit(generation_domain, learned_y)
		#exp_vals = exp_coeffs[1] .* exp.(exp_coeffs[2].*generation_domain)

		#polynomial function
		poly_terms = 4 #How many terms to consider in polynomial function
		poly_coeffs = poly_fit(generation_domain, learned_y, poly_terms)
		poly_vals = []
		for i in 1:poly_terms
			push!(poly_vals,poly_coeffs[i] .* generation_domain .^ (i-1))
		end
		poly_vals = sum(poly_vals)

		#Writing parameters to CSV file; will break if length of generation is longer than length of fit coeffs - fixable
		fit_params_df = DataFrame(learned_per_gen=learned,
			power_fit=[power_coeffs...; [missing for _ in 1:length(learned)-length(power_coeffs)]],
			polynomial_fit=[poly_coeffs...; [missing for _ in 1:length(learned)-length(poly_coeffs)]],
		)

		CSV.write("./../output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/data/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-fit_params.csv", fit_params_df)

		# plotting per set of parameters
		class_plots = []
		default()

		# Generating each frame
		for i in 1:num_generations
			class_plot = heatmap(generations[i], 
				title = "Classroom over time \n λ=$(λ₀) \n size: $(class_size)",
				aspect_ratio=:equal,
				cbar = true,
				showaxis = false,
				c = palette(:grays,rev=true),
				size = (512,512),
				yflip=true,
				dpi = 300,
				axis=([], false)
			)
	
			learned_plot = scatter(learned, legend=:topleft,
				xlabel = "Generations",
				ylabel = "Number of learned",
				title = "Learned over time",
				label = "Learned students",
				legend_font_pointsizes = 5,
				yrot = 0,
				ytickfontsize = 4,
				dpi = 300,
				leftmargin = 5mm,
				rightmargin = 5mm,
				markersize = 3,
				scale = :log10,
				xlims = (1,num_generations),
				ylims = (1, maximum(learned))
			);
	
			learned_plot = vline!([i], label = "Current generation: $(i)");
	
			learned_plot = plot!(power_vals, 
				label="Power fit: $(round(power_coeffs[1], digits=2))⋅x^$(round(power_coeffs[2],digits=2))"
			);
	
			#learned_plot = plot!(exp_vals)
	
			learned_plot = plot!(poly_vals,
				label = "Polynomial fit: $(round.(poly_coeffs,digits=2))"
			);
		
			class_plot = plot(class_plot,learned_plot, 
				layout = l,
				size=(512,512/0.75),
				dpi = 300
			);

			savefig(class_plot,"./../output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/images/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-$(i).png");
			push!(class_plots, class_plot)
		end
	
		anim = @animate for i in 1:num_generations
			plot(class_plots[i]);
		end
	
		mp4(anim, "./../output/2D-Binary-PCA/$(seat_config)-$(class_size)/$(λ₀)/2DBPCA-$(seat_config)-$(class_size)-$(λ₀)-animation.mp4", fps = 16)	
	end
end

# ╔═╡ 9def329a-c1fd-4dd2-869d-1d61c6fb425d
# ╠═╡ disabled = true
#=╠═╡
begin #main
	# List of parameters
	sizes = [32,64,128]
	seat_configs = ["inner_corner","center","outer_corner"]
	Λs = [0.25,0.5,0.75]
	steady_state_tolerance = 10

	# Making the directories
	folders = ["images","data"]
	for seat_config in seat_configs, 
		size in sizes, 
		λ in Λs, 
		folder in folders
		
		mkpath("./../output/2D-Binary-PCA/$(seat_config)-$(size)/$(λ)/$(folder)")
	end
	
	class_simulation(sizes,
		seat_configs,
		Λs,
		steady_state_tolerance
	)
	
end
  ╠═╡ =#

# ╔═╡ 1716999b-c986-4bfc-a249-6eeb22182b96
md"
### Curve fit parameters:
- We use the Linear Least Square (method?) to find the coefficients via `CurveFit.jl`
- For powerfit `power_fit(x, y)`: find coefficient a and b
$y_i=ax_i^b$
- For sum of exponential fit `expsum_fit(x, y, 2, withconst = true)`: finds coefficients k, p and λ
$y_i = k + p_1\cdot\exp{(λ_1 \cdot x_i)} + p_2 \cdot \exp{(\lambda_2 \cdot x_i)}$
- polynomial fit `polyfit(x,y,n)`: finds coefficient a[k]
$y_i = a_1 + a_2 \cdot x_i + a_3 \cdot x_i^2 + ... + a_{n+1} \cdot x_i^n$
"

# ╔═╡ df8e5d9f-f8cb-49dd-9156-a016ceee32dc
md"
### Getting the effective radius of the spreading from the simulation
```math
A = n \pi r^2
```
* where n is the number of seeds in the initial set up
```math
r = \sqrt{\frac{A}{n\pi}}
```
* where A is the number of learned per generation
* r is the effective radius of the learned
---
### Derivation
```math
A(r) = \pi r^2
```
```math
\frac{\,dA}{\,dt} = k \cdot 2\pi r
```
```math
\frac{\,d}{\,dt}(\pi r^2) = k \cdot C
```
```math
2\pi r \cdot \dot{r} = k \cdot 2\pi r
```
```math
\dot{r} = k
```
```math
\therefore \text{the rate of expansion of the circle is constant}
```
"

# ╔═╡ 25c26b82-9913-438b-9080-5ddb5cc0ffe0
md"""
# Testing
"""

# ╔═╡ d40d6b87-0dc9-4078-81b7-8b6bb81a43c1
md"""
## Calculations
"""

# ╔═╡ 4c37d313-1a61-4352-9ae6-c91d138add70
#calculations scratch
begin
	class_size = 128
	max_generations = 50
	
	λ₀ = 1

	λ = Float64.( Matrix(
	[ 	λ₀ 		λ₀ 		λ₀;
		λ₀ 		0 		λ₀;
		λ₀ 		λ₀ 		λ₀]
	))

	#Go until steady state
	generations, num_generations = simulate_steady_state("inner_corner",class_size,λ,1)
	learned = map(x->sum(x), generations)
	
	learned_y = learned[1:end]
	generation_domain = 1:length(learned_y)
	
	power_coeffs = power_fit(generation_domain, learned_y)
	power_vals = power_coeffs[1] .* generation_domain .^ power_coeffs[2]

	exp_coeffs = exp_fit(generation_domain, learned_y)
	exp_vals = exp_coeffs[1] .* exp.(exp_coeffs[2].*generation_domain)

	poly_terms = 4
	poly_coeffs = poly_fit(generation_domain, learned_y, poly_terms)
	poly_vals = []
	for i in 1:poly_terms
		push!(poly_vals,poly_coeffs[i] .* generation_domain .^ (i-1))
	end
	poly_vals = sum(poly_vals)

end

# ╔═╡ c5b4fd60-cb59-494c-ae99-800bb6ddb893
begin

	circle_radii = 1:length(learned)
	circle_area = 4 * (π .* circle_radii .^ 2)
	
	default()
	
	hline([class_size^2], label="finite size effect cap", linestyles=:dash)
	
	plot!(learned, 
		label = "Simulation",
		legend =:topleft,
		#scale =:log10,
		ylabel = "Number of learned",
		xlabel = "Generation number",
		#ylabel = "Number of learned (log₁₀)",
		#xlabel = "Generation number (log₁₀)",
		dpi = 300,
		minorticks = true,
		minorgrid = true
	)

	plot!(circle_area,
		label = "y = 4πx²"
	)
end

# ╔═╡ fbb39644-4521-44b7-9977-aee3fad519e8
begin
	circle_multiplier = sum(generations[1])
	effective_r = sqrt.(learned ./ (circle_multiplier*π))

	plot(effective_r
		, legend = false
		, ylabel = "Effective radius"
		, xlabel = "Generation count"
		, title = "Growth of effective radius over each generation"
		, dpi = 300
	)
end

# ╔═╡ cc33628a-ec2f-4095-8e0a-dc289dfbe208
md"""
## Plotting
"""

# ╔═╡ 9fc867f8-c574-4559-97ea-d816df11b4f3
#plotting scratch
begin #plotting part of the main function
	class_plots = []
	default()

	l = @layout [
		a{0.75h}
		b{}
	]
	
	for i in 1:num_generations
		class_plot = heatmap(generations[i], title = "Classroom over time $(size)",
			aspect_ratio=:equal,
			cbar = true,
			showaxis = false,
			c = palette(:grays,rev=true),
			size = (512,512),
			yflip=true,
			dpi = 300,
			axis=([], false),
		)

		learned_plot = scatter(learned, legend=:topleft,
			xlabel = "Generations",
			ylabel = "Number of learned",
			title = "Learned over time",
			label = "Learned students",
			legend_font_pointsizes = 5,
			yrot = 0,
			ytickfontsize = 4,
			dpi = 300,
			leftmargin = 5mm,
			rightmargin = 5mm,
			markersize = 3,
			scale =:log10
		)

		learned_plot = vline!([i], label = "Current generation: $(i)")

		learned_plot = plot!(power_vals, 
			label="Power fit: $(round(power_coeffs[1], digits=2))⋅x^$(round(power_coeffs[2],digits=2))"
		)

		#learned_plot = plot!(exp_vals)

		learned_plot = plot!(poly_vals,
			label = "Polynomial fit: $(round.(poly_coeffs,digits=2))"
		)
	
		class_plot = plot(class_plot,learned_plot, 
			layout = l,
			#size=(512,512/0.75),
			dpi = 300
		)
		
		push!(class_plots, class_plot)
	end

	class_plots[5]

	#=
	anim = @animate for i in 1:num_generations
		plot(class_plots[i])
	end

	gif(anim, "test_animation_scratch.gif", fps = 16)
	=#
	
end

# ╔═╡ 5b33e153-dc48-488f-8008-f5482d8cad54
begin
	learned_plot = scatter(learned, legend=:topleft,
			xlabel = "Generations",
			ylabel = "Number of learned",
			title = "Learned over time",
			label = "Learned students",
			legend_font_pointsizes = 5,
			yrot = 0,
			ytickfontsize = 4,
			dpi = 300,
			leftmargin = 5mm,
			rightmargin = 5mm,
			markersize = 3,
			scale =:log10,
			xlims = (1, num_generations),
			ylims = (1,maximum(learned))
		)

		learned_plot = vline!([4], label = "Current generation: $(4)")

		learned_plot = plot!(power_vals, 
			label="Power fit: $(round(power_coeffs[1], digits=2))⋅x^$(round(power_coeffs[2],digits=2))",
			#scale =:log10
		)

		#learned_plot = plot!(exp_vals)

		learned_plot = plot!(poly_vals,
			label = "Polynomial fit: $(round.(poly_coeffs,digits=2))",
			#scale =:log10
		)
end

# ╔═╡ f4809e61-c955-4c4c-9a08-283339d1f07f
learned_plot;

# ╔═╡ 6f326290-cbb6-446e-98a8-bb2a1079b848
md"""
## Raw data -> CSV
"""

# ╔═╡ 074e7215-83f7-4343-94b3-572a262277f9
begin #scratch for exporting csv
	df_cols = ["Generation $(i)" for i in 1:length(generations)]
	df_data = vec.(generations)
	
	# df_data -> student number -> student state at generation
	# state per student over time
	#df_data = [[df_data[g][s] for g in 1:length(df_data)] for s in 1:class_size^2]

	# row = ith student; column = jth generation
	student_states_df = DataFrame(df_data,df_cols)
	CSV.write("test_data_scratch.csv",student_states_df)

	#fit_params = ["power fit", power_coeffs...]

end

# ╔═╡ 538929f8-973d-4258-a10d-1cc8f8e87ed3
md"""
## Parameters -> CSV
"""

# ╔═╡ c9bd742e-9908-43ee-807c-9d2f0ba7473e
begin
	fit_params_df = DataFrame(learned_per_gen=learned,
		power_fit=[power_coeffs...; [missing for _ in 1:length(learned)-length(power_coeffs)]],
		polynomial_fit=[poly_coeffs...; [missing for _ in 1:length(learned)-length(poly_coeffs)]],
	)

	CSV.write("test_fit_params_scratch.csv", fit_params_df)
end

# ╔═╡ 54dc0a71-2ac7-4db8-b3c6-f513f31a4165
begin
	CSV.read("test_fit_params_scratch.csv",DataFrame)
end

# ╔═╡ b4a2643b-38c6-4b20-8d24-4fc151aaf106
md"
### Change log:
#### May 24, 2023:
* Added adjacency matrix method to generating the next generation. Slower than convolution. Maybe have better implementation of dictionaries. Make sure that you are multiplying matrices when using `.*`
### May 31, 2023:
* Added comments in `generate_next_generation3()` function for optimization. Also, can change the name to `next!()`
"

# ╔═╡ 0e6d9178-f9b7-48ca-b2cd-985f6ee2c21c
md"
#### 5-25-2023 notes
* Remove calculating the probability of learning from a loop. Duplicate the weighted adjacency matrix with M(number of student) copies of S(state of students) so that the operation is `64x64 matrix .* 64x64 matrix`. 
* Then take the product of each row (makes a column vector), and use it to get a vector P where Pᵢ is the probability of student i to learn.
* Make a function `decide((randᵢ,pᵢ))` where in the tuple is a made from the `zip()` function where the inputs are (1) vector of M random numbers, and (2) P. This allows us to use `decide.()` on the inputs to have a element-wise or function instead of `.||` (?)
* Make the input of `generate_next_generation` function the **`initial grid`** and **`adj matrix`**, instead of **`initial grid`** and **`λ`** so that it doesn't compute the adjacency matrix (does not change over time, double check) for every generation.
---
Make slides discussing the following
* Main equation: (for generating the next generation)
```math 
\vec{S}^{t+1} = \vec{S}^t \text{ or } (\vec{U} \leq \vec{P}^t )
```
where U is a vector of random numbers, S is the vector of states, and P the vector of the probabilities of student i to learn
* Agent description/characteristics (student ID and states)
* Neighborhood and the connection of the agents to their neighbors (λ)
* Rules of interaction (see equation above)
* Observables (purpose) - **DO THIS FIRST** (number of learned over time/rate of learning, spatio-temporal distribution of learned, etc)
---
Other things:
* Analytical reason of inner corner efficiency based on geometry where the area of learned is proportional to the circumference/perimeter
* Can use mean field approximation?
* Solution to #learned vs time is sigmoid (dS = λS(1-S))
---
* Check phone pics for whiteboard pictures
"

# ╔═╡ Cell order:
# ╟─a5ae8b1a-2249-4266-8dd3-7e6a45e637f3
# ╟─5798f6f7-4e75-4276-8a66-5d4518cde118
# ╠═bcc8a210-fa01-11ed-1db7-d7cfebd38074
# ╠═20ca43a0-7499-4e6e-b5af-44a56b3f83c7
# ╠═33b8010b-e92a-4a8e-a738-5441a2371119
# ╠═c54a1e59-e363-4567-b956-5b05ea26f172
# ╠═520615ae-62e5-4482-8473-9027c12eed8f
# ╠═8bf82d25-6643-4245-bee3-546c560708ef
# ╠═caea6e08-89c4-429c-b850-3e91b026ba80
# ╟─4b621695-f682-4bb8-b9d0-17284457e040
# ╟─6680e9cb-40e7-47c3-8c34-cf8a7df1e323
# ╟─0dbd7206-cf15-40d4-b520-36455a048a6d
# ╟─31d1deea-ba04-4b88-92c5-14cbe4733a6d
# ╠═64e02837-8603-49f1-a18f-0b0391ccb430
# ╟─186b2799-ede3-46d1-825c-55e67a56f4d3
# ╠═da8c68a1-a62a-4f49-a500-e58cced203d9
# ╠═9def329a-c1fd-4dd2-869d-1d61c6fb425d
# ╟─1716999b-c986-4bfc-a249-6eeb22182b96
# ╟─df8e5d9f-f8cb-49dd-9156-a016ceee32dc
# ╟─c5b4fd60-cb59-494c-ae99-800bb6ddb893
# ╟─fbb39644-4521-44b7-9977-aee3fad519e8
# ╟─25c26b82-9913-438b-9080-5ddb5cc0ffe0
# ╟─d40d6b87-0dc9-4078-81b7-8b6bb81a43c1
# ╠═4c37d313-1a61-4352-9ae6-c91d138add70
# ╟─cc33628a-ec2f-4095-8e0a-dc289dfbe208
# ╠═9fc867f8-c574-4559-97ea-d816df11b4f3
# ╠═f4809e61-c955-4c4c-9a08-283339d1f07f
# ╠═5b33e153-dc48-488f-8008-f5482d8cad54
# ╟─6f326290-cbb6-446e-98a8-bb2a1079b848
# ╠═074e7215-83f7-4343-94b3-572a262277f9
# ╟─538929f8-973d-4258-a10d-1cc8f8e87ed3
# ╠═c9bd742e-9908-43ee-807c-9d2f0ba7473e
# ╠═54dc0a71-2ac7-4db8-b3c6-f513f31a4165
# ╟─b4a2643b-38c6-4b20-8d24-4fc151aaf106
# ╟─0e6d9178-f9b7-48ca-b2cd-985f6ee2c21c
