################################################################################
#  Copyright 2020, Tom Van Acker, Arpan Koirala                                #
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
function run_sopf(file, model_type::Type, optimizer; kwargs...)
    return run_model(file, model_type, optimizer, build_sopf; kwargs...)
end

function build_sopf(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_gen_power(pm, nw=n, bounded=false)
        variable_bus_voltage(pm, nw=n, bounded=false)
        variable_branch_power(pm, nw=n, bounded=false)
        variable_dcline_power(pm, nw=n, bounded=false)

        constraint_model_voltage(pm, nw=n)

        for i in ids(pm, :ref_buses, nw=n)
            constraint_theta_ref(pm, i, nw=n)
        end

        for i in ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
        end

        for i in ids(pm, :branch, nw=n)
            constraint_ohms_yt_from(pm, i, nw=n)
            constraint_ohms_yt_to(pm, i, nw=n)

            constraint_voltage_angle_difference(pm, i, nw=n)

            constraint_thermal_limit_from(pm, i, nw=n)
            constraint_thermal_limit_to(pm, i, nw=n)
        end

        for i in ids(pm, :dcline, nw=n)
            constraint_dcline_power_losses(pm, i, nw=n)
        end
    end

    objective_min_fuel_and_flow_cost(pm)
end