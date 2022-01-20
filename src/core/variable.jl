################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# bus voltage
"variable: `ve[i]` for `i` in `bus`es"
function variable_bus_voltage_expectation(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false)
    ve = _PM.var(pm, nw)[:ve] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_ve",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "ve_start", 1.0)
    )

    if bounded
        for (i, bus) in _PM.ref(pm, nw, :bus)
            JuMP.set_lower_bound(ve[i], 0.0)
        end
    end
    
    if aux_fix 
        JuMP.fix.(ve, 1.0; force=true)
    end

    report && _PM.sol_component_value(pm, nw, :bus, :ve, _PM.ids(pm, nw, :bus), ve)
end

"variable: `vv[i]` for `i` in `bus`es"
function variable_bus_voltage_variance(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false)
    vv = _PM.var(pm, nw)[:vv] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_vv",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "vv_start", 1.0)
    )

    if bounded
        for (i, bus) in _PM.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vv[i],  0.0)
        end
    end
    
    if aux_fix 
        JuMP.fix.(vv, 1.0; force=true)
    end

    report && _PM.sol_component_value(pm, nw, :bus, :vv, _PM.ids(pm, nw, :bus), vv)
end

"variable: `vs[i]` for `i` in `bus`es"
function variable_bus_voltage_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false)
    vs = _PM.var(pm, nw)[:vs] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_vs",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "vs_start", 1.0)
    )

    if bounded
        for (i, bus) in _PM.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vs[i], -2.0 * bus["vmax"]^2)
            JuMP.set_upper_bound(vs[i],  2.0 * bus["vmax"]^2)
        end
    end
    
    if aux_fix 
        JuMP.fix.(vs, 1.0; force=true)
    end

    report && _PM.sol_component_value(pm, nw, :bus, :vs, _PM.ids(pm, nw, :bus), vs)
end

# load current
"variable: `crd[j]` for `j` in `load`"
function variable_load_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crd = _PM.var(pm, nw)[:crd] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_crd",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "crd_start")
    )

    report && _PM.sol_component_value(pm, nw, :load, :crd, _PM.ids(pm, nw, :load), crd)
end


"variable: `cid[j]` for `j` in `load`"
function variable_load_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cid = _PM.var(pm, nw)[:cid] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_cid",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cid_start")
    )

    report && _PM.sol_component_value(pm, nw, :load, :cid, _PM.ids(pm, nw, :load), cid)
end

# load current
"variable: `crd[j]` for `j` in `load`"
function variable_PV_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crd_pv = _PM.var(pm, nw)[:crd_pv] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_crd_pv",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "crd_pv_start")
    )

    report && _PM.sol_component_value(pm, nw, :load, :crd_pv, _PM.ids(pm, nw, :load), crd_pv)
end


"variable: `cid[j]` for `j` in `load`"
function variable_PV_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cid_pv = _PM.var(pm, nw)[:cid_pv] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_cid_pv",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cid_pv_start")
    )

    report && _PM.sol_component_value(pm, nw, :load, :cid_pv, _PM.ids(pm, nw, :load), cid_pv)
end

# branch current
"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cr = _PM.var(pm, nw)[:cr] = Dict()

    bus = _PM.ref(pm, nw, :bus)
    branch = _PM.ref(pm, nw, :branch)

    for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PM.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PM.var(pm, nw, :vr, i)
        vi_fr = _PM.var(pm, nw, :vi, i)
    
        vr_to = _PM.var(pm, nw, :vr, j)
        vi_to = _PM.var(pm, nw, :vi, j)
    
        csr_fr = _PM.var(pm, nw, :csr, l)
        csi_fr = _PM.var(pm, nw, :csi, l)

        cr[(l,i,j)] = (tr * csr_fr - ti * csi_fr + g_sh_fr * vr_fr - b_sh_fr * vi_fr) / tm^2
        cr[(l,j,i)] = -csr_fr + g_sh_to * vr_to - b_sh_to * vi_to

        # ub = Inf
        # if haskey(b, "rate_a")
        #     rate_fr = b["rate_a"]*b["tap"]
        #     rate_to = b["rate_a"]
        #     ub = max(rate_fr/bus[i]["vmin"], rate_to/bus[j]["vmin"])
        # end
        # if haskey(b, "c_rating_a")
        #     ub = b["c_rating_a"]
        # end

        # if !isinf(ub)
        #     JuMP.@constraint(pm.model, cr[(l,i,j)] >= -ub)
        #     JuMP.@constraint(pm.model, cr[(l,i,j)] <= ub)

        #     JuMP.@constraint(pm.model, cr[(l,j,i)] >= -ub)
        #     JuMP.@constraint(pm.model, cr[(l,j,i)] <= ub)
        # end
    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :cr_fr, :cr_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), cr)
end

"variable: `ci[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ci = _PM.var(pm, nw)[:ci] = Dict()

    bus = _PM.ref(pm, nw, :bus)
    branch = _PM.ref(pm, nw, :branch)

    for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PM.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PM.var(pm, nw, :vr, i)
        vi_fr = _PM.var(pm, nw, :vi, i)
    
        vr_to = _PM.var(pm, nw, :vr, j)
        vi_to = _PM.var(pm, nw, :vi, j)
    
        csr_fr = _PM.var(pm, nw, :csr, l)
        csi_fr = _PM.var(pm, nw, :csi, l)

        ci[(l,i,j)] = (tr * csi_fr + ti * csr_fr + g_sh_fr * vi_fr + b_sh_fr * vr_fr) / tm^2
        ci[(l,j,i)] = -csi_fr + g_sh_to * vi_to + b_sh_to * vr_to

        # ub = Inf
        # if haskey(b, "rate_a")
        #     rate_fr = b["rate_a"]*b["tap"]
        #     rate_to = b["rate_a"]
        #     ub = max(rate_fr/bus[i]["vmin"], rate_to/bus[j]["vmin"])
        # end
        # if haskey(b, "c_rating_a")
        #     ub = b["c_rating_a"]
        # end

        # if !isinf(ub)
        #     JuMP.@constraint(pm.model, ci[(l,i,j)] >= -ub)
        #     JuMP.@constraint(pm.model, ci[(l,i,j)] <= ub)

        #     JuMP.@constraint(pm.model, ci[(l,j,i)] >= -ub)
        #     JuMP.@constraint(pm.model, ci[(l,j,i)] <= ub)
        # end
    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :ci_fr, :ci_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), ci)
end

"variable: `cse[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_expectation(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, aux_fix::Bool=false, report::Bool=true)
    cse = _PM.var(pm, nw)[:cse] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_cse",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "cse_start", 0.0)
    )

    if bounded
        for l in _PM.ids(pm, nw, :branch)
            JuMP.set_lower_bound(cse[l], 0.0)
        end
    end

    if aux_fix
        JuMP.fix.(cse, 0.0; force=true)
    end

    report && _PM.sol_component_value(pm, nw, :branch, :cse, _PM.ids(pm, nw, :branch), cse)
end

"variable: `csv[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_variance(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, aux_fix::Bool=false, report::Bool=true)
    csv = _PM.var(pm, nw)[:csv] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_csv",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "csv_start", 0.0)
    )

    if bounded
        for l in _PM.ids(pm, nw, :branch)
                JuMP.set_lower_bound(csv[l], 0.0)
        end
    end

    if aux_fix
        JuMP.fix.(csv, 0.0; force=true)
    end

    report && _PM.sol_component_value(pm, nw, :branch, :csv, _PM.ids(pm, nw, :branch), csv)
end

"variable: `css[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, aux_fix::Bool=false, report::Bool=true)
    css = _PM.var(pm, nw)[:css] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_css",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "css_start", 0.0)
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
                JuMP.set_lower_bound(css[l], -2.0 * ub^2)
                JuMP.set_upper_bound(css[l],  2.0 * ub^2)
            end
        end
    end

    if aux_fix
        JuMP.fix.(css, 0.0; force=true)
    end

    report && _PM.sol_component_value(pm, nw, :branch, :css, _PM.ids(pm, nw, :branch), css)
end

# branch voltage drop
"variable: `vbdr[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_voltage_drop_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vbdr = _PM.var(pm, nw)[:vbdr] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_vbdr",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "vbdr_start", 0.01)
    )
    
    report && _PM.sol_component_value(pm, nw, :branch, :vbdr, _PM.ids(pm, nw, :branch), vbdr)
end

"variable: `vbdi[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_voltage_drop_img(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

    vbdi = _PM.var(pm, nw)[:vbdi] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_vbdi",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "vbdi_start", 0.01)
    )

    report && _PM.sol_component_value(pm, nw, :branch, :vbdi, _PM.ids(pm, nw, :branch), vbdi)
end

"expression: voltage drop real"
function expression_branch_voltage_drop_real(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    if !haskey(_PM.var(pm, nw), :vbdr)
         vbdr=_PM.var(pm, nw)[:vbdr] = Dict()
    end

        branch = _PM.ref(pm, nw, :branch, b)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        
        #vbdr[b] = (vr_fr-vr_to)

        expression_branch_voltage_drop_real(pm, nw, b, f_bus, t_bus)
end


"expression: voltage drop img"
function expression_branch_voltage_drop_img(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    if !haskey(_PM.var(pm, nw), :vbdi)
        vbdi=_PM.var(pm, nw)[:vbdi] = Dict()
    end
        #branch = _PMs.ref(pm, nw, :branch)
    branch = _PM.ref(pm, nw, :branch, b)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    
    expression_branch_voltage_drop_img(pm, nw, b, f_bus, t_bus)
end

