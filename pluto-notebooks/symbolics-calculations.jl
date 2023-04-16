### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ fa8e8760-6ad0-11ed-0dd7-7da0e425ea82
using Pkg

# ╔═╡ 63958d35-aacd-4554-96b5-80396b739449
begin
	Pkg.activate(".")
	Pkg.add("Symbolics")
	Pkg.add("SymbolicUtils")
	Pkg.add("Latexify")
end

# ╔═╡ 1a021caf-3f5f-4cea-8804-e7fd050cb334
begin
	using Symbolics
	using SymbolicUtils
	using Latexify
end

# ╔═╡ 420e8001-6ecb-4faf-bd41-d7b8d53913c4
@variables x y z R L C S

# ╔═╡ 32b1a086-77b7-4252-bdf2-e0cab378cee9
t = (1/S*L + S*C)^(-1) * ((1/R + S*C)^(-1) + (1/(S*L) + S*C)^(-1))

# ╔═╡ 79808bbf-5779-447a-8de5-e285e66d01fb
simplified_t = simplify(t; expand=true)

# ╔═╡ 57a54bf4-26c8-4ec3-9471-f18157cf4aa5
denom = (1+S^2*L*C) * (L+ S^2 * C) * (1+C*R*S)

# ╔═╡ a3ce1c30-2bc6-4105-b29c-034e9dae0225
simplify(denom; expand=true)

# ╔═╡ 87b907dd-c23b-4df3-a294-0a74d6d16a22
final_t = substitute.(simplified_t, (Dict(
	R => 1000.0, 
	L => 0.1, 
	C => 10*10^(-6)
),))

# ╔═╡ 0e38ad96-5a60-45dd-90ee-89c8e1a937d8


# ╔═╡ 969ad0bd-13f2-49fd-9057-82006f30dc0d
Z_L = S*L

# ╔═╡ 8853d977-9601-47ec-a81d-84a7ab738c27
Z_RC = (1/R + S*C)^(-1)

# ╔═╡ 666335a5-d15d-4c2c-a639-aa1eedcdba0d
transfer_function = (Z_L)/(Z_RC+Z_L)

# ╔═╡ 8f90ed4f-531c-40f2-bfe2-7d8d1ca3dd3a
simplified_transfer_function = simplify(transfer_function; expand=true)

# ╔═╡ bc6784c2-e299-40d9-afbd-81894cfe3c94
final_transfer_function = substitute.(simplified_transfer_function, (Dict(
	R => 1000.0, +
	L => 0.1, 
	C => 100*10^(-7)
),))

# ╔═╡ 1dc092f2-b1d3-4649-8af5-93ff022bdd2a
G = 2x/(S*(S+3))

# ╔═╡ c04aef1e-2b78-4ce2-8d4f-07e99bb0fc36
simplify(G/(1+G); expand = true)

# ╔═╡ Cell order:
# ╠═fa8e8760-6ad0-11ed-0dd7-7da0e425ea82
# ╠═63958d35-aacd-4554-96b5-80396b739449
# ╠═1a021caf-3f5f-4cea-8804-e7fd050cb334
# ╠═420e8001-6ecb-4faf-bd41-d7b8d53913c4
# ╠═32b1a086-77b7-4252-bdf2-e0cab378cee9
# ╠═79808bbf-5779-447a-8de5-e285e66d01fb
# ╠═57a54bf4-26c8-4ec3-9471-f18157cf4aa5
# ╠═a3ce1c30-2bc6-4105-b29c-034e9dae0225
# ╠═87b907dd-c23b-4df3-a294-0a74d6d16a22
# ╠═0e38ad96-5a60-45dd-90ee-89c8e1a937d8
# ╠═969ad0bd-13f2-49fd-9057-82006f30dc0d
# ╠═8853d977-9601-47ec-a81d-84a7ab738c27
# ╠═666335a5-d15d-4c2c-a639-aa1eedcdba0d
# ╠═8f90ed4f-531c-40f2-bfe2-7d8d1ca3dd3a
# ╠═bc6784c2-e299-40d9-afbd-81894cfe3c94
# ╠═1dc092f2-b1d3-4649-8af5-93ff022bdd2a
# ╠═c04aef1e-2b78-4ce2-8d4f-07e99bb0fc36
