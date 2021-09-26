################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

"variable: `ve[i]` for `i` in `bus`es"
function variable_bus_voltage_expectation(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false)
    ve = _PMs.var(pm, nw)[:ve] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_ve",
        start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "ve_start", 1.0)
    )

    if bounded
        for (i, bus) in _PMs.ref(pm, nw, :bus)
            JuMP.set_lower_bound(ve[i],  0.0)
        end
    end
    
    if aux_fix 
        JuMP.fix.(ve, 1.0; force=true)
    end

    report && sol_component_value(pm, nw, :bus, :ve, _PMs.ids(pm, nw, :bus), ve)
    # report && _IMs.sol_component_value(pm, nw, :bus, :vs, _PMs.ids(pm, nw, :bus), vs)
end

"variable: `vv[i]` for `i` in `bus`es"
function variable_bus_voltage_variance(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false)
    vv = _PMs.var(pm, nw)[:vv] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_vv",
        start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "vv_start", 1.0)
    )

    if bounded
        for (i, bus) in _PMs.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vv[i],  0.0)
        end
    end
    
    if aux_fix 
        JuMP.fix.(vv, 1.0; force=true)
    end

    report && sol_component_value(pm, nw, :bus, :vv, _PMs.ids(pm, nw, :bus), vv)
    # report && _IMs.sol_component_value(pm, nw, :bus, :vs, _PMs.ids(pm, nw, :bus), vs)
end

"variable: `vs[i]` for `i` in `bus`es"
function variable_bus_voltage_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false)
    vs = _PMs.var(pm, nw)[:vs] = JuMP.@variable(pm.model,
        [i in _PMs.ids(pm, nw, :bus)], base_name="$(nw)_vs",
        start = comp_start_value(_PMs.ref(pm, nw, :bus, i), "vs_start", 1.0)
    )

    if bounded
        for (i, bus) in _PMs.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vs[i], -2.0 * bus["vmax"]^2)
            JuMP.set_upper_bound(vs[i],  2.0 * bus["vmax"]^2)
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

"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cr = _PMs.var(pm, nw)[:cr] = Dict()

    bus = _PMs.ref(pm, nw, :bus)
    branch = _PMs.ref(pm, nw, :branch)

    for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PMs.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PMs.var(pm, nw, :vr, i)
        vi_fr = _PMs.var(pm, nw, :vi, i)
    
        vr_to = _PMs.var(pm, nw, :vr, j)
        vi_to = _PMs.var(pm, nw, :vi, j)
    
        csr_fr = _PMs.var(pm, nw, :csr, l)
        csi_fr = _PMs.var(pm, nw, :csi, l)

        cr[(l,i,j)] = (tr*csr_fr - ti*csi_fr + g_sh_fr*vr_fr - b_sh_fr*vi_fr)/tm^2
        cr[(l,j,i)] = -csr_fr + g_sh_to*vr_to - b_sh_to*vi_to

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

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :branch, :cr_fr, :cr_to, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), cr)
end

"variable: `ci[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ci = _PMs.var(pm, nw)[:ci] = Dict()

    bus = _PMs.ref(pm, nw, :bus)
    branch = _PMs.ref(pm, nw, :branch)

    for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PMs.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PMs.var(pm, nw, :vr, i)
        vi_fr = _PMs.var(pm, nw, :vi, i)
    
        vr_to = _PMs.var(pm, nw, :vr, j)
        vi_to = _PMs.var(pm, nw, :vi, j)
    
        csr_fr = _PMs.var(pm, nw, :csr, l)
        csi_fr = _PMs.var(pm, nw, :csi, l)

        ci[(l,i,j)] = (tr*csi_fr + ti*csr_fr + g_sh_fr*vi_fr + b_sh_fr*vr_fr)/tm^2
        ci[(l,j,i)] = -csi_fr + g_sh_to*vi_to + b_sh_to*vr_to

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

    report && _IMs.sol_component_value_edge(pm, _PMs.pm_it_sym, nw, :branch, :ci_fr, :ci_to, _PMs.ref(pm, nw, :arcs_from), _PMs.ref(pm, nw, :arcs_to), ci)
end

"variable: `cse[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_expectation(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, aux_fix::Bool=false, report::Bool=true)
    cse = _PMs.var(pm, nw)[:cse] = JuMP.@variable(pm.model,
        [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_cse",
        start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "cse_start", 0.0)
    )

    if bounded
        for l in _PMs.ids(pm, nw, :branch)
                #JuMP.set_lower_bound(cse[l], 0.0)
        end
    end

    if aux_fix
        JuMP.fix.(cse, 0.0; force=true)
    end

    report && sol_component_value(pm, nw, :branch, :cse, _PMs.ids(pm, nw, :branch), cse)
    # report && _IMs.sol_component_value(pm, nw, :branch, :css, _PMs.ids(pm, nw, :branch), css)
end

"variable: `csv[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_variance(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, aux_fix::Bool=false, report::Bool=true)
    csv = _PMs.var(pm, nw)[:csv] = JuMP.@variable(pm.model,
        [l in _PMs.ids(pm, nw, :branch)], base_name="$(nw)_csv",
        start = comp_start_value(_PMs.ref(pm, nw, :branch, l), "csv_start", 0.0)
    )

    if bounded
        for l in _PMs.ids(pm, nw, :branch)
                JuMP.set_lower_bound(csv[l], 0.0)
        end
    end

    if aux_fix
        JuMP.fix.(csv, 0.0; force=true)
    end

    report && sol_component_value(pm, nw, :branch, :csv, _PMs.ids(pm, nw, :branch), csv)
    # report && _IMs.sol_component_value(pm, nw, :branch, :css, _PMs.ids(pm, nw, :branch), css)
end


"variable:"
function variable_gen_power_real_variance(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, aux_fix::Bool=false, report::Bool=true)
    pgv = _PMs.var(pm, nw)[:pgv] = JuMP.@variable(pm.model,
        [g in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_pgv",
        start = comp_start_value(_PMs.ref(pm, nw, :gen, g), "pgv_start", 0.0)
    )

    if bounded
        for g in _PMs.ids(pm, nw, :gen)
                JuMP.set_lower_bound(pgv[g], 0.0)
                JuMP.set_upper_bound(pgv[g], 10.0)
        end
    end

    report && sol_component_value(pm, nw, :gen, :pgv, _PMs.ids(pm, nw, :gen), pgv)
    # report && _IMs.sol_component_value(pm, nw, :branch, :css, _PMs.ids(pm, nw, :branch), css)
end

"variable:"
function variable_gen_power_img_variance(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, aux_fix::Bool=false, report::Bool=true)
    qgv = _PMs.var(pm, nw)[:qgv] = JuMP.@variable(pm.model,
        [g in _PMs.ids(pm, nw, :gen)], base_name="$(nw)_qgv",
        start = comp_start_value(_PMs.ref(pm, nw, :gen, g), "qgv_start", 0.0)
    )

    if bounded
        for g in _PMs.ids(pm, nw, :gen)
                JuMP.set_lower_bound(qgv[g], 0.0)
                JuMP.set_upper_bound(qgv[g], 10.0)
        end
    end

    report && sol_component_value(pm, nw, :gen, :qgv, _PMs.ids(pm, nw, :gen), qgv)
    # report && _IMs.sol_component_value(pm, nw, :branch, :css, _PMs.ids(pm, nw, :branch), css)
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