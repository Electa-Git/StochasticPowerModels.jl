################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# using pkgs
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels

# constants 
const _PMs = PowerModels
const _SPM = StochasticPowerModels

# data

# path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case118_spm.m")
# path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case3_spm.m")
path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case30_spm_muhlpfordt.m")

# include(joinpath(_SPM.BASE_DIR,"test/io.jl"))
data = _PMs.parse_file(path)

for (b, branch) in data["branch"]
    branch["b_fr"] = 0
    branch["b_to"] = 0
end
for (s, shunt) in data["shunt"]
    shunt["gs"] = 0
    shunt["bs"] = 0
end


# add_stochastic_data!(data)

# initialize solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# solve problem
result_stc = run_sopf_iv(data, _PMs.IVRPowerModel, solver, deg = 1)
# deg = 1
# sdata = _SPM.build_stochastic_data(data, deg)
# pm = instantiate_model(sdata, _PMs.IVRPowerModel, _SPM.build_sopf_iv)
# result = optimize_model!(pm, optimizer=solver)
    
##
# solve problem iteratively
(result_dtr, result_itr) = run_sopf_iv_itr(data, _PMs.IVRPowerModel, solver, deg = 1);

# assert
@assert isapprox(result_stc["objective"], result_itr["objective"], rtol=1e-6)