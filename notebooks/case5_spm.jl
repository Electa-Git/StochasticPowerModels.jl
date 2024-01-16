### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 24944b1e-3420-4cb3-9d32-a5f00108f0c9
begin
	import Pkg
    # url for StochasticPowerModels
	spm_url = "https://github.com/timmyfaraday/StochasticPowerModels.jl.git"
	
	# activate a temporary environment
    Pkg.activate(mktempdir())
	
	Pkg.add(Pkg.PackageSpec(name="Distributions"))
	Pkg.add(Pkg.PackageSpec(name="GR"))
	Pkg.add(Pkg.PackageSpec(name="Images"))
	Pkg.add(Pkg.PackageSpec(name="Ipopt"))
	Pkg.add(Pkg.PackageSpec(name="PowerModels"))
	Pkg.add(Pkg.PackageSpec(name="Plots"))
	
    Pkg.add(Pkg.PackageSpec(url=spm_url))
	
	Pkg.build("GR")
	
    using Distributions, Images, Ipopt, Plots
	using PowerModels, StochasticPowerModels
	
	const PM  = PowerModels
	const SPM = StochasticPowerModels
	
	solver = Ipopt.Optimizer
	
	md"""
	# Notebook on PCE-CC-OPF (ACR vs IVR)
	"""
end

# ╔═╡ 014e2e15-a4a3-4570-95ce-cccfc669bc43
md"""
## Network Specification

In this notebook, a [five-bus system](https://github.com/power-grid-lib/pglib-opf/blob/master/pglib_opf_case5_pjm.m) is considered, with three loads, four generators and six branches. 
"""

# ╔═╡ 95a2af1a-cc26-4bc0-bb8f-7e7e7b275b3a
load("figures/five_bus_system.png")

# ╔═╡ e27c7b66-f9f8-419a-b491-8ccfdbe9de21
md"""
A stochastic germ is considered comprised of two sources of uncertainty: 
- Beta distribution 𝓑(2.0,2.0) for load at bus two
- Normal distribution 𝓝(0.0,1.0) for load at bus three
	
| i | distribution 	| μₚ [MW] | σₚ [MW] 	| Q [Mvar]  |
|---|---------------|-----------|---------------|-----------|
| 2 | 𝓑(2.0,2.0)	 | 300.0 	 | 30.0 		 | 98.61 	 |
| 3 | 𝓝(0.0,1.0)    | 300.0     | 30.0 		 | 98.61     |
| 4 | -  			| 400.0 	| - 			| 131.47 	|
"""

# ╔═╡ 0171d114-086d-4dda-9e7f-094c7621ba28
begin
	grid = "case5_spm.m"
	file = joinpath(SPM.BASE_DIR, "test/data/matpower", grid)
end;

# ╔═╡ 171640e4-b830-4c9d-b9ef-a7d59862982c
md"""
## Problem Specification
1. Problem Formulation `frm`:
    - ACR `acr`, 
    - IVR `ivr`, or 
    - Reduced IVR `red_ivr`
2. Auxiliary Variables `aux`: 
    - `true`, or 
    - `false`
3. PCE Degree `deg`: 
    - `1`, or
    - `2`
"""

# ╔═╡ fcf621e6-fe1b-4e86-b233-a0b245e40cb7
begin
	frm = "ivr"
	aux = true
	deg = 1
end;

# ╔═╡ 2ca1c1e9-2333-4e52-b106-750c1c190126
begin
	case 	   = String[]
	case_study = String[]
	objective  = Real[]
	solve_time = Real[]
	voltage    = Vector{Real}[]
	
md"""
## Run the Appropriate Convenience Functions for Stochastic OPF
"""
end

# ╔═╡ 622523f5-6cb4-4c1b-a0bc-493a92f40465
begin
	if aux
		push!(case, "$(frm) w. aux (deg=$(deg))")
		push!(case_study, "$(frm) \n w. aux \n deg=$(deg)")
	else
		push!(case, "$(frm) w/o aux (deg=$(deg))")
		push!(case_study, "$(frm) \n w/o aux \n deg=$(deg)")
	end
	
	if frm == "acr"
		pmf = PM.ACRPowerModel
		res = run_sopf_acr(file, pmf, solver, aux=aux, deg=deg)
	elseif frm == "ivr"
		pmf = PM.IVRPowerModel
		res = run_sopf_iv(file, pmf, solver, aux=aux, deg=deg, red=false)
	elseif frm == "red_ivr"
		pmf = PM.IVRPowerModel
		res = run_sopf_iv(file, pmf, solver, aux=aux, deg=deg, red=true)
	end		
end;

# ╔═╡ 7ce7de61-9d9c-41ba-86bf-8722b1ab5aea
@assert res["termination_status"] == PM.LOCALLY_SOLVED

# ╔═╡ 9276dba7-4889-4f1f-9fc7-f7df592f2fe0
begin
	push!(objective, res["objective"] / 1e6)
	Plots.scatter(1:length(objective),objective,
					title="Objective Value",
					legend=false,
					xaxis=false,
					xgrid=false,
					xrotation=90,
					xticks=(collect(1:length(objective)),case_study),
					ylabel="Objective Value [M€]",
					ylim=(0.99*minimum(objective),1.01*maximum(objective))
				)
end

# ╔═╡ 959304fc-7c4a-48c8-ab94-87ec968ed422
begin
	push!(solve_time, res["solve_time"])
	Plots.scatter(1:length(solve_time),solve_time,
					title="Solve Time",
					legend=false,
					xaxis=false,
					xgrid=false,
					xrotation=90,
					xticks=(collect(1:length(objective)),case_study),
					ylabel="Solve Time [s]",
					ylim=(0.0,1.25*maximum(solve_time)),
					yminorgrid=true
				)
end

# ╔═╡ 4f33f0e4-2b0c-41c7-8a9f-f6a9b1d835ea
begin
	labels = Array{String}(undef,1,length(case))
	for (ni,nc) in enumerate(case) labels[1,ni] = nc end
	push!(voltage, pdf(SPM.density(res, "bus", 2, "vm", sample_size=25000), 0.96:0.0001:1.02))
	Plots.plot(0.96:0.0001:1.02,voltage,
				title="Voltage Magnitude at Bus Two",
				labels=labels,
				xaxis="Voltage Magnitude [V]"
	)
end

# ╔═╡ Cell order:
# ╟─24944b1e-3420-4cb3-9d32-a5f00108f0c9
# ╟─014e2e15-a4a3-4570-95ce-cccfc669bc43
# ╟─95a2af1a-cc26-4bc0-bb8f-7e7e7b275b3a
# ╟─e27c7b66-f9f8-419a-b491-8ccfdbe9de21
# ╠═0171d114-086d-4dda-9e7f-094c7621ba28
# ╟─171640e4-b830-4c9d-b9ef-a7d59862982c
# ╠═fcf621e6-fe1b-4e86-b233-a0b245e40cb7
# ╟─2ca1c1e9-2333-4e52-b106-750c1c190126
# ╠═622523f5-6cb4-4c1b-a0bc-493a92f40465
# ╠═7ce7de61-9d9c-41ba-86bf-8722b1ab5aea
# ╟─9276dba7-4889-4f1f-9fc7-f7df592f2fe0
# ╟─959304fc-7c4a-48c8-ab94-87ec968ed422
# ╟─4f33f0e4-2b0c-41c7-8a9f-f6a9b1d835ea
