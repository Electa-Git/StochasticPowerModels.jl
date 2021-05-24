################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels(Distribution).jl for Stochastic (Optimal)#
# Power Flow                                                                   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

module StochasticPowerModels

    # import pkgs
    import Distributions
    import InfrastructureModels
    import JuMP
    import LinearAlgebra
    import PolyChaos
    import PowerModels
    import PowerModelsDistribution

    # import types
    import PowerModels: AbstractPowerModel, AbstractACRModel, AbstractIVRModel
    import PowerModels: comp_start_value
    import InfrastructureModels: ids, ref, var, con, sol, nw_ids, nws, sol_component_value

    # pkgs const
    const _DST = Distributions
    const _IMs = InfrastructureModels
    const _PCE = PolyChaos
    const _PMs = PowerModels
    const _PMD = PowerModelsDistribution

    # const 
    const nw_id_default = 1

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # include
    include("core/constraint_template.jl")
    include("core/objective.jl")
    include("core/variable.jl")

    include("form/iv.jl")
    include("form/acr.jl")

    include("prob/sopf_iv.jl")
    include("prob/sopf_iv_itr.jl")
    include("prob/sopf_acr.jl")

    include("util/util.jl")

    # export
    export BASE_DIR

    export run_sopf_iv, run_sopf_iv_itr

end 
