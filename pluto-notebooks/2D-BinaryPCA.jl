### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ bcc8a210-fa01-11ed-1db7-d7cfebd38074
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 20ca43a0-7499-4e6e-b5af-44a56b3f83c7
begin
	Pkg.add("Plots")
	Pkg.add("Random")
	Pkg.add("BenchmarkTools")
end

# ╔═╡ c54a1e59-e363-4567-b956-5b05ea26f172
begin
	using Random
	using Plots
	using BenchmarkTools
end

# ╔═╡ a5ae8b1a-2249-4266-8dd3-7e6a45e637f3
md"""
# 2D Binary PCA
"""

# ╔═╡ 5798f6f7-4e75-4276-8a66-5d4518cde118
md"""
## Setting up packages
"""

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
	initial_grid_size = size(initial_grid)
	one_row_size = (1,prod(initial_grid_size))
	s⃗ = reshape(initial_grid,one_row_size)

	p = 1 .- (Λ .* (fill(1,64) * s⃗))

	p⃗ = 1 .- prod(p,dims=2)

	u⃗ = rand(64,1)

	#dont know what to name this
	rand_results = transpose(u⃗ .<= p⃗)

	s⃗ᵗ⁺¹ = BitMatrix(s⃗).||rand_results

	output = reshape(s⃗ᵗ⁺¹,initial_grid_size)

	return output
end

# ╔═╡ 4c37d313-1a61-4352-9ae6-c91d138add70
#main
begin
	class_size = 8
	initial_class = initiate_grid("inner_corner", class_size)

	λ = Float64.( Matrix(
	[ 	0.5 		0.5 		0.5;
		0.5 		0.0 		0.5;
		0.5 		0.5 		0.5]
))

	#Currently only works for square classrooms
	Λ = generate_neighbor_list(reshape(1:class_size^2,class_size,class_size),λ)

	output = generate_next_generation3(initial_class,Λ)
end

# ╔═╡ 824dd823-18c2-42cb-a238-df2eade0ed5a
generate_next_generation3(initial_class,Λ)

# ╔═╡ 2b519268-a947-40e0-8818-75de7148fd41
heatmap(output)

# ╔═╡ 7ec1daf4-a1b1-4b9e-912f-123e35ffc866
test_grid = initiate_grid("inner_corner",8)

# ╔═╡ 29f8c82c-9bfc-456a-93ac-a4028887fb80
test_λ = Float64.( Matrix(
	[ 	1 		1 		1;
		1 		0 		1;
		1 		1 		1]
))

# ╔═╡ 3f4efea6-526f-41f4-9416-6a4d19094f1e
@benchmark generate_next_generation(initial_class,λ)

# ╔═╡ 14140ee9-a968-46f5-865c-987a250073ae
@benchmark generate_next_generation2(initial_class,Λ)

# ╔═╡ 18f74419-f32b-4418-9e5a-d5846c9cf279
@benchmark generate_next_generation3(initial_class,Λ)

# ╔═╡ db01cf96-02aa-48ab-8077-831b08ea0dcb
#=╠═╡
test_adjacency_matrix = generate_neighbor_list(reshape(1:64,8,8),test_λ)
  ╠═╡ =#

# ╔═╡ 70bdc2eb-f29b-4953-bfbc-60617d98c79a
#=╠═╡
heatmap(test_adjacency_matrix, yflip=false, cbar=true, c=cgrad(:grays,rev=true), clims=(0,1))
  ╠═╡ =#

# ╔═╡ 7bd55f31-4804-4b83-846f-d6587df0d861
#=╠═╡
test_adjacency_matrix
  ╠═╡ =#

# ╔═╡ 2458c1ae-a1a2-4a96-a324-326ec003ad16
#=╠═╡
test_dict = Dict( [(student,test_adjacency_matrix[student,:]) for student in 1:size(test_adjacency_matrix)[1]])
  ╠═╡ =#

# ╔═╡ b4a2643b-38c6-4b20-8d24-4fc151aaf106
md"
### Change log:
#### May 24, 2023:
* Added adjacency matrix method to generating the next generation. Slower than convolution. Maybe have better implementation of dictionaries. Make sure that you are multiplying matrices when using `.*`
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
# ╠═c54a1e59-e363-4567-b956-5b05ea26f172
# ╠═520615ae-62e5-4482-8473-9027c12eed8f
# ╟─8bf82d25-6643-4245-bee3-546c560708ef
# ╟─caea6e08-89c4-429c-b850-3e91b026ba80
# ╟─4b621695-f682-4bb8-b9d0-17284457e040
# ╠═6680e9cb-40e7-47c3-8c34-cf8a7df1e323
# ╠═0dbd7206-cf15-40d4-b520-36455a048a6d
# ╠═31d1deea-ba04-4b88-92c5-14cbe4733a6d
# ╠═824dd823-18c2-42cb-a238-df2eade0ed5a
# ╠═4c37d313-1a61-4352-9ae6-c91d138add70
# ╠═2b519268-a947-40e0-8818-75de7148fd41
# ╠═7ec1daf4-a1b1-4b9e-912f-123e35ffc866
# ╠═29f8c82c-9bfc-456a-93ac-a4028887fb80
# ╠═3f4efea6-526f-41f4-9416-6a4d19094f1e
# ╠═14140ee9-a968-46f5-865c-987a250073ae
# ╠═18f74419-f32b-4418-9e5a-d5846c9cf279
# ╠═db01cf96-02aa-48ab-8077-831b08ea0dcb
# ╠═70bdc2eb-f29b-4953-bfbc-60617d98c79a
# ╠═7bd55f31-4804-4b83-846f-d6587df0d861
# ╠═2458c1ae-a1a2-4a96-a324-326ec003ad16
# ╟─b4a2643b-38c6-4b20-8d24-4fc151aaf106
# ╟─0e6d9178-f9b7-48ca-b2cd-985f6ee2c21c
