### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 5d60e89c-2806-45c9-a2b9-b5ef76745e13
start_num = 0

# ╔═╡ c999cc2e-3192-11ee-315e-8b1e61925b8c
max_num = 8

# ╔═╡ 264185a2-f69f-4e24-b68e-5c7d30589149
exponent = 2

# ╔═╡ e1f7baf4-0693-4b1c-8159-86415c53cdae
size = 2^exponent

# ╔═╡ 7ad69827-0b5c-4c5f-ad72-692dd8f3b6d4
begin
	out = []
	for x in start_num:max_num-1 
		#julia ranges are inclusive (unlike python, kaya may -1)
		if x%size < size÷2
			push!(out,x)
		end
	end
end

# ╔═╡ b3d309aa-6455-456f-9428-8ef529354cbc
out

# ╔═╡ 79f0a6e2-3142-4989-b779-e7c1685e5922
out .+ size÷2

# ╔═╡ 2978903e-3e95-4d68-8548-4e890faa7de1
s = [x for x in start_num:max_num-1 if x%size < size÷2]

# ╔═╡ 6705811b-2a3a-4bdc-bc81-5497e9ab7635
s .+ 2

# ╔═╡ Cell order:
# ╠═5d60e89c-2806-45c9-a2b9-b5ef76745e13
# ╠═c999cc2e-3192-11ee-315e-8b1e61925b8c
# ╠═264185a2-f69f-4e24-b68e-5c7d30589149
# ╠═e1f7baf4-0693-4b1c-8159-86415c53cdae
# ╠═7ad69827-0b5c-4c5f-ad72-692dd8f3b6d4
# ╠═b3d309aa-6455-456f-9428-8ef529354cbc
# ╠═79f0a6e2-3142-4989-b779-e7c1685e5922
# ╠═2978903e-3e95-4d68-8548-4e890faa7de1
# ╠═6705811b-2a3a-4bdc-bc81-5497e9ab7635
