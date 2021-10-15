################################################################################
#  Copyright 2021, Arpan Koirala, Tom Van Acker                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# NOTE: dc lines are omitted from the current formulation                      #
################################################################################

""
function run_sopf_acr(sdata, model_constructor::Type, optimizer; aux::Bool=true,  deg::Int=1, kwargs...)
    if aux
        return _PMs.run_model(sdata, model_constructor, optimizer, build_sopf_acr_with_aux; multinetwork=true, kwargs...)
    else
        return _PMs.run_model(sdata, model_constructor, optimizer, build_sopf_acr_without_aux; multinetwork=true, kwargs...)
    end
end

""
function build_sopf_acr_with_aux(pm::AbstractPowerModel)
    for (n, network) in _PMs.nws(pm) 
        variable_bus_voltage(pm, nw=n, aux=true)

        variable_gen_power(pm, nw=n, bounded=false)

        variable_branch_power(pm, nw=n, bounded=false)
        variable_branch_current(pm, nw=n, bounded=false, aux=true) 
    end

    for i in _PMs.ids(pm, :bus,nw=1)
        constraint_bus_voltage_squared_cc_limit(pm, i, nw=1)
    end

    for g in _PMs.ids(pm, :gen, nw=1)
        constraint_gen_power_cc_limit(pm, g, nw=1)
    end

    for b in _PMs.ids(pm, :branch, nw=1)
       constraint_branch_series_current_squared_cc_limit(pm, b, nw=1)
    end

    for (n, network) in _PMs.nws(pm)
        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
            constraint_gp_bus_voltage_squared(pm, i, nw=n)
        end

        for b in _PMs.ids(pm, :branch, nw=n)                                     
            constraint_gp_power_branch_to(pm, b, nw=n)
            constraint_gp_power_branch_from(pm, b, nw=n)

            constraint_branch_voltage(pm, b, nw=n)

            constraint_gp_current_squared(pm, b, nw=n)
        end
    end

    objective_min_expected_generation_cost(pm)
end

""
function build_sopf_acr_without_aux(pm::AbstractPowerModel)
    for (n, network) in _PMs.nws(pm) 
        variable_bus_voltage(pm, nw=n, aux=false)

        variable_gen_power(pm, nw=n, bounded=false)

        variable_branch_power(pm, nw=n, bounded=false)
        variable_branch_current(pm, nw=n, bounded=true, aux=false) 
    end

    for i in _PMs.ids(pm, :bus,nw=1)
        constraint_bus_voltage_cc_limit(pm, i, nw=1)
    end

    for g in _PMs.ids(pm, :gen, nw=1)
        constraint_gen_power_cc_limit(pm, g, nw=1)
    end

    for b in _PMs.ids(pm, :branch, nw=1)
       constraint_branch_series_current_cc_limit(pm, b, nw=1)
    end

    for (n, network) in _PMs.nws(pm)

        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            constraint_gp_power_branch_to(pm, b, nw=n)
            constraint_gp_power_branch_from(pm, b, nw=n)

            constraint_branch_voltage(pm, b, nw=n)
        end
    end

    objective_min_expected_generation_cost(pm)
end