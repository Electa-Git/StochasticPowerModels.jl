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

    include("form.jl")

end
