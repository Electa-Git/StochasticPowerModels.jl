################################################################################
#  Copyright 2021, Tom Van Acker, Arpan Koirala                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# input -- subject to change
deg  = 1
aux  = true
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
data  = _PM.parse_file(path)
sdata = _SPM.build_stochastic_data(data, deg)

# initialize solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)

# solve problem
result = _SPM.run_sopf_acr(sdata, _PM.ACRPowerModel, solver, aux=aux, deg=deg)
