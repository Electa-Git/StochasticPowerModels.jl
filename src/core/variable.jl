################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# bus
""
function variable_bus_voltage(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_bus_voltage_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"variable: `vms[i]` for `i` in `bus`"
function variable_bus_voltage_magnitude_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vms = _PM.var(pm, nw)[:vms] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_vms",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "vms_start", 1.0)
    )

    if bounded
        for (i, bus) in _PM.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vms[i], -2.0 * bus["vmax"]^2)
            JuMP.set_upper_bound(vms[i],  2.0 * bus["vmax"]^2)
        end
    end

    report && _PM.sol_component_value(pm, nw, :bus, :vms, _PM.ids(pm, nw, :bus), vms)
end

# branch
"variable: `cmss[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_magnitude_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cmss = _PM.var(pm, nw)[:cmss] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_cmss",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "cmss_start", 0.0)
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        branch = _PM.ref(pm, nw, :branch)

        for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
            b = branch[l]
            ub = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"] * b["tap"]
                y_fr = abs(b["g_fr"] + im * b["b_fr"])
                y_to = abs(b["g_to"] + im * b["b_to"])
                shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
                series_current = max(rate / bus[i]["vmin"], rate / bus[j]["vmin"])
                ub = series_current + shunt_current
            end
            if haskey(b, "c_rating_a")
                total_current = b["c_rating_a"]
                y_fr = abs(b["g_fr"] + im * b["b_fr"])
                y_to = abs(b["g_to"] + im * b["b_to"])
                shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
                ub = total_current + shunt_current
            end

            if !isinf(ub)
                JuMP.set_lower_bound(cmss[l], -2.0 * ub^2)
                JuMP.set_upper_bound(cmss[l],  2.0 * ub^2)
            end
        end
    end

    report && _PM.sol_component_value(pm, nw, :branch, :cmss, _PM.ids(pm, nw, :branch), cmss)
end

# generator
""
function variable_gen_power(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_gen_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end