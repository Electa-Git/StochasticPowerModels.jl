################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels(Distribution).jl for Stochastic (Optimal)#
# Power Flow                                                                   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
"sOPF with ACPPowerModel"
""
function run_ac_sopf(file, optimizer; kwargs...)
    return run_sopf(file, ACPPowerModel, optimizer; kwargs...)
end

""
function run_sopf_iv(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_sopf_iv; kwargs...)
end

function build_sopf_iv(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        # variables
        _PMs.variable_bus_voltage(pm, nw=n, bounded=false)

        variable_unit_current(pm, nw=n, bounded=false)
        variable_branch_current(pm, nw=n, bounded=false)

        # equality constraints
        for i in ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
        end
        for i in ids(pm, :branch, nw=n)
            constraint_current_from(pm, i, nw=n)
            constraint_current_to(pm, i, nw=n)
            constraint_voltage_drop(pm, i, nw=n)
        end
        
        # bounds 
        for i in ids(pm, :gen, nw=n) if active_bounds(i,n)
            constraint_gen_bounds(pm, i, nw=n)
        end end
        for i in ids(pm, :load, nw=n) if active_bounds(i,n)
            constraint_load_bounds(pm, i, nw=n)
        end end
        for i in ids(pm, :branch, nw=n) if active_bounds(i,n)
            constraint_thermal_limit_from(pm, i, n)
            constraint_thermal_limit_to(pm, i, n)
        end end
    end

    objective_min_fuel_and_flow_cost(pm)

end