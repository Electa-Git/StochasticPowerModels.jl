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

function variable_active_dcbranch_flow(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=false, report::Bool=true, kwargs...)
    _PMACDC.variable_active_dcbranch_flow(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_dcbranch_current(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=false, report::Bool=true, kwargs...)
    _PMACDC.variable_dcbranch_current(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_dcgrid_voltage_magnitude(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=false, report::Bool=true, kwargs...)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_dc_converter(pm::_PM.AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true,report::Bool=false, kwargs...)  
    _PMACDC.variable_filter_voltage_real(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_filter_voltage_imaginary(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_voltage_real(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_voltage_imaginary(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_transformer_current_real_from(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_transformer_current_real_to(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_transformer_current_imaginary_from(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_transformer_current_imaginary_to(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_reactor_current_real_from(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_reactor_current_real_to(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_reactor_current_imaginary_from(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_reactor_current_imaginary_to(pm, nw=nw, bounded=bounded,kwargs...)
    _PMACDC.variable_converter_current_real(pm, nw=nw, bounded=bounded,  kwargs...)
    _PMACDC.variable_converter_current_imaginary(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_current_dc(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_current_lin(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_active_power(pm, nw=nw, bounded=bounded,  kwargs...)
    _PMACDC.variable_converter_reactive_power(pm, nw=nw, bounded=bounded,  kwargs...)
    _PMACDC.variable_dcside_power(pm, nw=nw, bounded=bounded,  kwargs...)


    ################
    _PMACDC.variable_conv_tranformer_flow(pm, nw=nw, bounded=bounded,  kwargs...)
    _PMACDC.variable_conv_reactor_flow(pm, nw=nw, bounded=bounded, kwargs...)

    _PMACDC.variable_converter_active_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMACDC.variable_converter_reactive_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMACDC.variable_acside_current(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_dcside_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMACDC.variable_converter_firing_angle(pm, nw=nw, bounded=bounded, kwargs...)

    _PMACDC.variable_converter_filter_voltage(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_internal_voltage(pm, nw=nw, bounded=bounded, kwargs...)

    _PMACDC.variable_converter_to_grid_active_power(pm, nw=nw, bounded=bounded, kwargs...)
    _PMACDC.variable_converter_to_grid_reactive_power(pm, nw=nw, bounded=bounded, kwargs...)
end

function variable_dc_converter_squared(pm::_PM.AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true,report::Bool=true, kwargs...)
    variable_filter_voltage_squared(pm, nw=nw, bounded=bounded, kwargs...)
    variable_converter_voltage_squared(pm, nw=nw, bounded=bounded, kwargs...)
    variable_transformer_current_from_squared(pm, nw=nw, bounded=bounded,  kwargs...)
    variable_transformer_current_to_squared(pm, nw=nw, bounded=bounded,  kwargs...)
    variable_reactor_current_from_squared(pm, nw=nw, bounded=bounded,  kwargs...)
    variable_reactor_current_to_squared(pm, nw=nw, bounded=bounded,  kwargs...)
    variable_converter_current_lin_squared(pm, nw=nw, bounded=bounded,  kwargs...)
end


function variable_filter_voltage_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vk_s = _PM.var(pm, nw)[:vk_s] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_vk_s",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "vr_start", 1.0)
    )

    if bounded
        for (i, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(vk_s[i], (1/2)*convdc["Vmmin"]^2)
            # JuMP.set_lower_bound(vk_s[i], -2*convdc["Vmmax"]^2)
            JuMP.set_upper_bound(vk_s[i], 2*convdc["Vmmax"]^2)
        end
    end
    report && _PM.sol_component_value(pm, nw, :convdc, :vk_s, _PM.ids(pm, nw, :convdc), vk_s)
    #report && _IM.sol_component_value(pm, _PM.pm_it_sym, nw, :convdc, :vk_s, _PM.ids(pm, nw, :convdc), vk_s)
end

function variable_converter_voltage_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

    vc_s = _PM.var(pm, nw)[:vc_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_vc_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "v_start", 1.0)
    )

    if bounded
        for (i, convdc) in _PM.ref(pm, nw, :convdc)
            JuMP.set_lower_bound(vc_s[i], (1/2)*convdc["Vmmin"]^2)
            # JuMP.set_lower_bound(vc_s[i], -2*convdc["Vmmax"]^2)
            JuMP.set_upper_bound(vc_s[i], 2*convdc["Vmmax"]^2)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :vc_s, _PM.ids(pm, nw, :convdc), vc_s)

end

function variable_transformer_current_from_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 2;
    vpu = 1;
    iik_s = _PM.var(pm, nw)[:iik_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_iik_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            # JuMP.set_lower_bound(iik_s[c],  -(convdc["Pacrated"]^2)/vpu * bigM)
            JuMP.set_lower_bound(iik_s[c],  0)
            JuMP.set_upper_bound(iik_s[c],  (convdc["Pacrated"]^2)/vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :iik_s, _PM.ids(pm, nw, :convdc), iik_s)
end


function variable_transformer_current_to_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 2;
    vpu = 1;
    iki_s = _PM.var(pm, nw)[:iki_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_iki_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            # JuMP.set_lower_bound(iki_s[c],  -(convdc["Pacrated"]^2)/vpu * bigM)
            JuMP.set_lower_bound(iki_s[c],  0)
            JuMP.set_upper_bound(iki_s[c],  (convdc["Pacrated"]^2)/vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :iki_s, _PM.ids(pm, nw, :convdc), iki_s)
end

function variable_reactor_current_from_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 2;
    vpu = 1;
    ikc_s = _PM.var(pm, nw)[:ikc_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_ikc_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            # JuMP.set_lower_bound(ikc_s[c],  -(convdc["Pacrated"]^2)/vpu * bigM)
            JuMP.set_lower_bound(ikc_s[c],  0)
            JuMP.set_upper_bound(ikc_s[c],  (convdc["Pacrated"]^2)/vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :ikc_s, _PM.ids(pm, nw, :convdc), ikc_s)
end

function variable_reactor_current_to_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 2;
    vpu = 1;
    ick_s = _PM.var(pm, nw)[:ick_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_ick_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            # JuMP.set_lower_bound(ick_s[c],  -(convdc["Pacrated"]^2)/vpu * bigM)
            JuMP.set_lower_bound(ick_s[c],  0)
            JuMP.set_upper_bound(ick_s[c],  (convdc["Pacrated"]^2)/vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :ick_s, _PM.ids(pm, nw, :convdc), ick_s)
end

function variable_converter_current_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 2;
    vpu = 1;
    ic_s = _PM.var(pm, nw)[:ic_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_ic_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            # JuMP.set_lower_bound(ic_s[c],  -(convdc["Imax"]^2) * bigM)
            JuMP.set_lower_bound(ic_s[c],  0)
            JuMP.set_upper_bound(ic_s[c],  (convdc["Imax"]^2) * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :ic_s, _PM.ids(pm, nw, :convdc), ic_s)
end


function variable_converter_current_lin_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool = true, report::Bool=true)
    bigM = 2;
    vpu = 1;
    iconv_lin_s = _PM.var(pm, nw)[:iconv_lin_s] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :convdc)], base_name="$(nw)_iconv_lin_s",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :convdc, i), "P_g", 1.0)
    )
    if bounded
        for (c, convdc) in _PM.ref(pm, nw, :convdc)
            # JuMP.set_lower_bound(iconv_lin_s[c],  -(convdc["Imax"]^2) * vpu * bigM)
            JuMP.set_lower_bound(iconv_lin_s[c],  0)
            JuMP.set_upper_bound(iconv_lin_s[c],  (convdc["Imax"]^2) * vpu * bigM)
        end
    end

    report && _PM.sol_component_value(pm, nw, :convdc, :iconv_lin_s, _PM.ids(pm, nw, :convdc), iconv_lin_s)
end

function variable_RES_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_RES_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_RES_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_RES_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crd_RES = _PM.var(pm, nw)[:crd_RES] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :RES)], base_name="$(nw)_crd_RES",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :RES, i), "crd_RES_start")
    )
    # if bounded
    #     for (i, RES) in _PM.ref(pm, nw, :RES)
    #         JuMP.set_lower_bound(crd_RES[i], 0)
    #     end
    # end
    
    report && _PM.sol_component_value(pm, nw, :RES, :crd_RES, _PM.ids(pm, nw, :RES), crd_RES)
end


function variable_RES_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cid_RES = _PM.var(pm, nw)[:cid_RES] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :RES)], base_name="$(nw)_cid_RES",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :RES, i), "cid_RES_start")
    )
    # if bounded
    #     for (i, RES) in _PM.ref(pm, nw, :RES)
    #         JuMP.set_lower_bound(cid_RES[i], 0)
    #     end
    # end

    report && _PM.sol_component_value(pm, nw, :RES, :cid_RES, _PM.ids(pm, nw, :RES), cid_RES)
end

function variable_redispatch(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    re_pg = _PM.var(pm, nw)[:re_pg] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_re_pg",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "re_pg_start")
    )

    # for (g, gen) in _PM.ref(pm, nw, :gen)
    #     JuMP.set_lower_bound(re_pg[g], 0)
    # end

    report && _PM.sol_component_value(pm, nw, :gen, :re_pg, _PM.ids(pm, nw, :gen), re_pg)
end