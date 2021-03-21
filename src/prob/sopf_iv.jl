################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels(Distribution).jl for Stochastic (Optimal)#
# Power Flow                                                                   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

""
function run_sopf_iv(file, model_constructor, optimizer; kwargs...)
    return run_model(file, model_constructor, optimizer, build_sopf_iv; kwargs...)
end

""
function build_sopf_iv(pm::AbstractPowerModel)
    for (n, network) in nws(pm) 
        variable_bus_voltage(pm, nw=n, bounded=false)
        variable_branch_current(pm, nw=n, bounded=false)

        variable_load_current(pm, nw=n, bounded=false)
        variable_gen_current(pm, nw=n, bounded=false)
        variable_gen_power(pm, nw=n, bounded=false)

        for i in ids(pm, :ref_buses)
            _PMs.constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus)
            constraint_current_balance(pm, i, nw=n)
            constraint_gp_squared_bus_voltage(pm, i, nw=n)
        end

        for b in ids(pm, :branch)
            _PMs.constraint_current_from(pm, b, nw=n)
            _PMs.constraint_current_to(pm, b, nw=n)

            _PMs.constraint_voltage_drop(pm, b, nw=n)

            constraint_gp_squared_branch_current(pm, b, nw=n)
        end

        for g in ids(pm, :gen)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in ids(pm, :load)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    for i in ids(pm, :bus)
        constraint_bus_voltage_squared_cc_limit(pm, i)
    end

    for g in ids(pm, :gen)
        constraint_gen_power_cc_limit(pm, g)
    end

    for b in ids(pm, :branch)
        constraint_branch_series_current_squared_cc_limit(pm, b)
    end

    objective_min_expected_fuel_cost(pm)
end