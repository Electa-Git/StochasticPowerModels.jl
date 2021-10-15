################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# input -- subject to change
deg  = 1
aux  = true
red  = true
case = "matpower/case30_spm.m"

# using pkgs
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels

# pkg constants 
const _PMs = PowerModels
const _SPM = StochasticPowerModels

# data
path  = joinpath(_SPM.BASE_DIR,"test/data/$case")
data  = _PMs.parse_file(path)
sdata = build_stochastic_data(data, deg)

# initialize solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)

# solve problem
result = _SPM.run_sopf_iv(sdata, _PMs.IVRPowerModel, solver, aux=aux, deg=deg, red=red)