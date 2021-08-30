################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

"variable: `vs[i]` for `i` in `bus`es"
function variable_bus_voltage_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false)
    vs = _PMs.var(pm, nw)[:vs] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_vs",
        start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "vs_start", 1.0)
    )

    if bounded
        for (i, bus) in _PMs.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vs[i], bus["vmin"]^2)
            JuMP.set_upper_bound(vs[i], bus["vmax"]^2)
        end
    end
    
    if aux_fix 
        JuMP.fix.(vs, 1.0; force=true)
    end

    report && sol_component_value(pm, nw, :bus, :vs, _PMs.ids(pm, nw, :bus), vs)
    # report && _IMs.sol_component_value(pm, nw, :bus, :vs, _PMs.ids(pm, nw, :bus), vs)
end

"variable: `crd[j]` for `j` in `load`"
function variable_load_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crd = _PMs.var(pm, nw)[:crd] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :load)], base_name="$(nw)_crd",
        start = comp_start_value(_PMs.ref(pm, nw, :load, i), "crd_start")
    )

    report && sol_component_value(pm, nw, :load, :crd, _PMs.ids(pm, nw, :load), crd)
    # report && _IMs.sol_component_value(pm, nw, :load, :crd, _PMs.ids(pm, nw, :load), crd)
end


"variable: `cid[j]` for `j` in `load`"
function variable_load_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cid = _PMs.var(pm, nw)[:cid] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :load)], base_name="$(nw)_cid",
        start = comp_start_value(_PMs.ref(pm, nw, :load, i), "cid_start")
    )

    report && sol_component_value(pm, nw, :load, :cid, _PMs.ids(pm, nw, :load), cid)
    # report && _IMs.sol_component_value(pm, nw, :load, :cid, _PMs.ids(pm, nw, :load), cid)
end

"variable: `css[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, aux_fix::Bool=false, report::Bool=true)
    css = _PMs.var(pm, nw)[:css] = JuMP.@variable(pm.model,
        [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_css",
        start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "css_start", 0.0)
    )

    if bounded
        bus = _PMs.ref(pm, nw, :bus)
        branch = _PMs.ref(pm, nw, :branch)

        for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
            b = branch[l]
            ub = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"]*b["tap"]
                y_fr = abs(b["g_fr"] + im*b["b_fr"])
                y_to = abs(b["g_to"] + im*b["b_to"])
                shunt_current = max(y_fr*bus[i]["vmax"]^2, y_to*bus[j]["vmax"]^2)
                series_current = max(rate/bus[i]["vmin"], rate/bus[j]["vmin"])
                ub = series_current + shunt_current
            end
            if haskey(b, "c_rating_a")
                total_current = b["c_rating_a"]
                y_fr = abs(b["g_fr"] + im*b["b_fr"])
                y_to = abs(b["g_to"] + im*b["b_to"])
                shunt_current = max(y_fr*bus[i]["vmax"]^2, y_to*bus[j]["vmax"]^2)
                ub = total_current + shunt_current
            end

            if !isinf(ub)
                JuMP.set_lower_bound(css[l], -2.0 * ub^2)
                JuMP.set_upper_bound(css[l],  2.0 * ub^2)
            end
        end
    end

    if aux_fix
        JuMP.fix.(css, 0.0; force=true)
    end

    report && sol_component_value(pm, nw, :branch, :css, _PMs.ids(pm, nw, :branch), css)
    # report && _IMs.sol_component_value(pm, nw, :branch, :css, _PMs.ids(pm, nw, :branch), css)
end


"variable: voltage drop real for `(l,i,j)` in `arcs_from`"
function variable_branch_voltage_drop_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vbdr = _PMs.var(pm, nw)[:vbdr] = JuMP.@variable(pm.model,
        [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_vbdr",
        start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "vbdr_start", 0.0)
    )
    
    report && sol_component_value(pm, nw, :branch, :vbdr, _PMs.ids(pm, nw, :branch), vbdr)
    # report && _IMs.sol_component_value(pm, nw, :branch, :vbdr, _PMs.ids(pm, nw, :branch), vbdr)
end
"variable: voltage drop img for `(l,i,j)` in `arcs_from`"
function variable_branch_voltage_drop_img(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

    vbdi = _PMs.var(pm, nw)[:vbdi] = JuMP.@variable(pm.model,
        [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_vbdi",
        start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "vbdi_start", 0.0)
    )

    report && sol_component_value(pm, nw, :branch, :vbdi, _PMs.ids(pm, nw, :branch), vbdi)
    # report && _IMs.sol_component_value(pm, nw, :branch, :vbdi, _PMs.ids(pm, nw, :branch), vbdi)
end