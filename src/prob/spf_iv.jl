################################################################################
#  Copyright 2024, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# NOTE: dc lines are omitted from the current formulation                      #
################################################################################

""
function solve_spf_iv(data::Dict, model_constructor, optimizer; solution_processors=[sol_data_model!], kwargs...)
    result = _PM.solve_model(data, model_constructor, optimizer, build_spf_iv; multinetwork=true, solution_processors=solution_processors, kwargs...)
    result["mop"] = data["mop"]

    return result
end 

""
function build_spf_iv(pm::_PM.AbstractIVRModel)
    for (n,~) in _PM.nws(pm)
        variable_bus_voltage(pm, nw=n, bounded=false)

        variable_branch_current(pm, nw=n, bounded=false)

        variable_gen_power(pm, nw=n, bounded=false)
        variable_gen_current(pm, nw=n, bounded=false)
        variable_load_current(pm, nw=n, bounded=false)
    end

    objective_load_flow(pm)

    for (n,~) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_voltage_drop(pm, b, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    println(pm.model)
end