### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 14604187-05fa-40ce-8545-f49dd2883db3
using Pkg

# ╔═╡ 2847efe9-de0f-40f7-a7b9-4f3d7217551f
Pkg.activate(".")

# ╔═╡ 073789a1-40da-4eaa-b714-1120a1c0fc17
# ╠═╡ show_logs = false
begin
	Pkg.add("Plots")
	Pkg.add("FFTW")
	Pkg.add("Distributions")
	Pkg.add("KernelDensity")
	Pkg.add("Random")
end

# ╔═╡ 7a89dcb0-7e4d-4094-a2ce-d03e50f03d70
begin
	using Plots
	using FFTW
	using Distributions
	using KernelDensity
	using Random
end

# ╔═╡ 48cb029f-a45d-48e1-9c08-56af67c78d9f
begin
	activities = 18
	ppa = [3,2,1] .* 10 # points per activity - reduced to GCF
	total_points = 18 * sum(ppa)
	average = total_points/3
	standard_req =  0.7 * total_points
	honor_req =  0.9 * total_points

	runs = 1000
	
	award = []
	finals1 = []
	finals2 = []
	finals3 = []
	
	for i in 1:runs
		points1 = []
		points2 = []
		points3 = []
		
		for steps in 1:activities
			places = shuffle([1,2,3])
			push!(points1,ppa[places[1]])
			push!(points2,ppa[places[2]])
			push!(points3,ppa[places[3]])
		end

		final1 = sum(points1)
		final2 = sum(points2)
		final3 = sum(points3)

		
		if final1>honor_req && final2>honor_req && final3>honor_req
			push!(award,2)
			
		elseif final1>standard_req && final2>standard_req && final2>standard_req
			push!(award,1)
			
		else
			push!(award,0)
		end

		push!(finals1,sum(points1))
		push!(finals2,sum(points2))
		push!(finals3,sum(points3))
		
	end

	
end

# ╔═╡ b3b99668-b19a-4f13-ba22-aefff8b44d35
ppa

# ╔═╡ 9f2017d6-d8be-4b42-8a54-81072ae8688f
total_points, standard_req,honor_req

# ╔═╡ 415ae00b-0b4c-4eaf-8d16-5b02c477cd8a
# Number of items members achieve honor
honor_count = count(==(2),award)/runs

# ╔═╡ e87659dd-a17b-4fe4-a752-29bcb25c7967
# Number of times members achieve standard
standard_count = count(==(1),award)/runs

# ╔═╡ 8eb886da-002b-4850-937e-84b9b16399d7
# Number of times members achieve nothing
no_award_count = count(==(0),award)/runs

# ╔═╡ 805351fb-15c2-45fe-9cf9-597d8a4205d7
mean(finals1)

# ╔═╡ 4e869802-485d-47db-a245-39e2c5e74f73
std(finals1)

# ╔═╡ 2553a2e7-db6f-4e37-bb9f-9bf51222033e
histogram(finals1)

# ╔═╡ 2e8e460e-d90e-42de-bc9f-50c8a18571f5
quantile(finals1,[0.9])

# ╔═╡ Cell order:
# ╠═14604187-05fa-40ce-8545-f49dd2883db3
# ╠═2847efe9-de0f-40f7-a7b9-4f3d7217551f
# ╠═073789a1-40da-4eaa-b714-1120a1c0fc17
# ╠═7a89dcb0-7e4d-4094-a2ce-d03e50f03d70
# ╠═48cb029f-a45d-48e1-9c08-56af67c78d9f
# ╠═b3b99668-b19a-4f13-ba22-aefff8b44d35
# ╠═9f2017d6-d8be-4b42-8a54-81072ae8688f
# ╠═415ae00b-0b4c-4eaf-8d16-5b02c477cd8a
# ╠═e87659dd-a17b-4fe4-a752-29bcb25c7967
# ╠═8eb886da-002b-4850-937e-84b9b16399d7
# ╠═805351fb-15c2-45fe-9cf9-597d8a4205d7
# ╠═4e869802-485d-47db-a245-39e2c5e74f73
# ╠═2553a2e7-db6f-4e37-bb9f-9bf51222033e
# ╠═2e8e460e-d90e-42de-bc9f-50c8a18571f5
