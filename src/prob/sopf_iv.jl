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
function run_sopf_iv(sdata, model_constructor, optimizer; aux::Bool=true, deg::Int=1, red::Bool=true, kwargs...)
    if aux && red
        return _PM.run_model(sdata, model_constructor, optimizer, build_sopf_iv_reduced_with_aux; multinetwork=true, kwargs...)
    elseif aux && !red    
        return _PM.run_model(sdata, model_constructor, optimizer, build_sopf_iv_with_aux; multinetwork=true, kwargs...)
    elseif !aux &&  red
        return _PM.run_model(sdata, model_constructor, optimizer, build_sopf_iv_reduced_without_aux; multinetwork=true, kwargs...)
    elseif !aux && !red
        return _PM.run_model(sdata, model_constructor, optimizer, build_sopf_iv_without_aux; multinetwork=true, kwargs...)
    end
end

""
function build_sopf_iv_with_aux(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        variable_bus_voltage(pm, nw=n, aux=true)

        variable_branch_current(pm, nw=n, aux=true)

        variable_gen_power(pm, nw=n, bounded=false)                             # enforcing bounds alters the objective 
        variable_gen_current(pm, nw=n, bounded=false)                           # enforcing bounds makes problem infeasible
        variable_load_current(pm, nw=n)
    end

    for i in _PM.ids(pm, :bus, nw=1)
        constraint_bus_voltage_squared_cc_limit(pm, i, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_gen_power_cc_limit(pm, g, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_branch_series_current_squared_cc_limit(pm, b, nw=1)
    end

    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
            constraint_gp_bus_voltage_squared(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_current_from(pm, b, nw=n)
            _PM.constraint_current_to(pm, b, nw=n)

            _PM.constraint_voltage_drop(pm, b, nw=n)

            constraint_gp_branch_series_current_squared(pm, b, nw=n)
        end

        for g in _PM.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    objective_min_expected_generation_cost(pm)
end

""
function build_sopf_iv_without_aux(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        variable_bus_voltage(pm, nw=n, aux=false)

        variable_branch_current(pm, nw=n, aux=false)

        variable_gen_power(pm, nw=n, bounded=false)                             # enforcing bounds alters the objective 
        variable_gen_current(pm, nw=n, bounded=false)                           # enforcing bounds makes problem infeasible
        variable_load_current(pm, nw=n)
    end

    for i in _PM.ids(pm, :bus, nw=1)
        constraint_bus_voltage_cc_limit(pm, i, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_gen_power_cc_limit(pm, g, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_branch_series_current_cc_limit(pm, b, nw=1)
    end

    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_current_from(pm, b, nw=n)
            _PM.constraint_current_to(pm, b, nw=n)

            _PM.constraint_voltage_drop(pm, b, nw=n)
        end

        for g in _PM.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    objective_min_expected_generation_cost(pm)
end

""
function build_sopf_iv_reduced_with_aux(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        variable_bus_voltage(pm, nw=n, aux=true)

        variable_branch_current_reduced(pm, nw=n, aux=true)

        variable_gen_power(pm, nw=n, bounded=false)                             # enforcing bounds alters the objective 
        variable_gen_current(pm, nw=n, bounded=false)                           # enforcing bounds makes problem infeasible
        variable_load_current(pm, nw=n)
    end

    for i in _PM.ids(pm, :bus, nw=1)
        constraint_bus_voltage_squared_cc_limit(pm, i, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_gen_power_cc_limit(pm, g, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_branch_series_current_squared_cc_limit(pm, b, nw=1)
    end

    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
            constraint_gp_bus_voltage_squared(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_voltage_drop(pm, b, nw=n)

            constraint_gp_branch_series_current_squared(pm, b, nw=n)
        end

        for g in _PM.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    objective_min_expected_generation_cost(pm)
end

""
function build_sopf_iv_reduced_without_aux(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        variable_bus_voltage(pm, nw=n, aux=false)

        variable_branch_current_reduced(pm, nw=n, aux=false)

        variable_gen_power(pm, nw=n, bounded=false)                             # enforcing bounds alters the objective 
        variable_gen_current(pm, nw=n, bounded=false)                           # enforcing bounds makes problem infeasible
        variable_load_current(pm, nw=n)
    end

    for i in _PM.ids(pm, :bus, nw=1)
        constraint_bus_voltage_cc_limit(pm, i, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_gen_power_cc_limit(pm, g, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_branch_series_current_cc_limit(pm, b, nw=1)
    end

    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_voltage_drop(pm, b, nw=n)
        end

        for g in _PM.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    objective_min_expected_generation_cost(pm)
end