### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 14604187-05fa-40ce-8545-f49dd2883db3
using Pkg

# ╔═╡ 073789a1-40da-4eaa-b714-1120a1c0fc17
# ╠═╡ show_logs = false
begin
	Pkg.activate(".")
	Pkg.add("Plots")
	Pkg.add("FFTW")
end

# ╔═╡ 7a89dcb0-7e4d-4094-a2ce-d03e50f03d70
begin
	using Plots
	using FFTW
end

# ╔═╡ d27136f0-e89a-4678-9c53-7a6c8b0b6f7d
begin
	dt = 0.01
	N = 50
	f = 5
	t = range(1, length=N, step=dt)
end

# ╔═╡ d6f53aae-80bb-47f1-b777-99a4da447069
begin
	y = sin.(2*pi*f.*t.+0.3)
	Fourier = fft(y)
end

# ╔═╡ 195e5e38-d767-4fce-b99e-2465a2b0bdd5
begin
	plot(t,y)
	xlabel!("time")
	ylabel!("amplitude")
	title!("signal")
end

# ╔═╡ 61baefb8-b1c4-4906-8dd5-c0f304c2d027
begin
	omega = exp.((-(t.-0.25).^2)./(2 .* 0.05))
	plot(t,omega)
	ynew = omega .* y
	plot!(t,ynew)
end

# ╔═╡ 61ad0489-b13b-4a54-9e35-e6da53e5c3bd
begin
	df = 1 / ((N-1)*dt)
	fmax = 1/(2*dt)
	fnew = range(-fmax, step=df, stop=fmax)
	plot(fnew,fftshift(abs.(Fourier)))
end

# ╔═╡ Cell order:
# ╠═14604187-05fa-40ce-8545-f49dd2883db3
# ╠═073789a1-40da-4eaa-b714-1120a1c0fc17
# ╠═7a89dcb0-7e4d-4094-a2ce-d03e50f03d70
# ╠═d27136f0-e89a-4678-9c53-7a6c8b0b6f7d
# ╠═d6f53aae-80bb-47f1-b777-99a4da447069
# ╠═195e5e38-d767-4fce-b99e-2465a2b0bdd5
# ╠═61baefb8-b1c4-4906-8dd5-c0f304c2d027
# ╠═61ad0489-b13b-4a54-9e35-e6da53e5c3bd
