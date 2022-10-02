################################################################################
#  Copyright 2021, Tom Van Acker, Frederik Geth                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# NOTE: dc lines are omitted from the current formulation                      #
################################################################################

""
function solve_sopf_iv(data::Dict, model_constructor, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    @assert _IM.ismultinetwork(data) == false "The data supplied is multinetwork, it should be single-network"
    @assert model_constructor <: _PM.AbstractIVRModel "This problem type only supports the IVRModel"
    
    sdata = build_stochastic_data(data, deg)
    result = _PM.solve_model(sdata, model_constructor, optimizer, build_sopf_iv; multinetwork=true, solution_processors=solution_processors, kwargs...)
    result["mop"] = sdata["mop"]
    
    return result
end

""
function solve_sopf_iv(file::String, model_constructor, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    data = _PM.parse_file(file)
    
    return solve_sopf_iv(data, model_constructor, optimizer; deg=deg, solution_processors=solution_processors, kwargs...)
end

""
function build_sopf_iv(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        variable_bus_voltage(pm, nw=n)

        variable_branch_current(pm, nw=n)

        variable_gen_power(pm, nw=n, bounded=false)                             # enforcing bounds alters the objective 
        variable_gen_current(pm, nw=n, bounded=false)                           # enforcing bounds makes problem infeasible
        variable_load_current(pm, nw=n)
    end

    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
            constraint_gp_bus_voltage_magnitude_squared(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_voltage_drop(pm, b, nw=n)

            constraint_gp_branch_series_current_magnitude_squared(pm, b, nw=n)
        end

        for g in _PM.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    for i in _PM.ids(pm, :bus, nw=1)
        constraint_cc_bus_voltage_magnitude_squared(pm, i, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_cc_branch_series_current_magnitude_squared(pm, b, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_cc_gen_power(pm, g, nw=1)
    end

    objective_min_expected_generation_cost(pm)
end