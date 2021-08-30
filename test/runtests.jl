################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# load pkgs 
using Test
using Ipopt
using PolyChaos
using PowerModels
using StochasticPowerModels

# constants 
const _PMs = PowerModels
const _SPM = StochasticPowerModels

# solvers
ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer,"max_cpu_time"=>300.0,
                                                         "tol"=>1e-9,
                                                         "print_level"=>0)

@testset "StochasticPowerModels.jl" begin

    # include("form.jl")
    include("prob.jl")

end
