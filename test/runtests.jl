################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# load pkgs 
using Test
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels

# constants 
const _PM = PowerModels
const _SPM = StochasticPowerModels

# solvers
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer,  "max_cpu_time"  => 600.0,
                                                                "tol"           => 1e-8,
                                                                "print_level"   => 0)

@testset "StochasticPowerModels.jl" begin

    include("prob.jl")
    include("util.jl")

end
