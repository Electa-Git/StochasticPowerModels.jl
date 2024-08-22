###################################################################################
#  Copyright 2024, Kaan Yurtseven                                                 #
###################################################################################
# StochasticPowerModels.jl                                                        #
# An extention package of PowerModels.jl for Stochastic Power System Optimization #
# See http://github.com/Electa-Git/StochasticPowerModels.jl                       #
###################################################################################

# variable
## branch 
""
function variable_branch_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    expression_variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_branch_series_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_branch_current_on_off(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    expression_variable_branch_current_real_on_off(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_variable_branch_current_imaginary_on_off(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_branch_series_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
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

function expression_variable_branch_current_real_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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

        z_branch = _PM.var(pm, _FP.first_id(pm,nw,:PCE_coeff), :z_branch, l)

        cr[(l,i,j)] = z_branch * ((tr * csr_fr - ti * csi_fr + g_sh_fr * vr_fr - b_sh_fr * vi_fr)) / tm^2
        cr[(l,j,i)] = z_branch * (-csr_fr + g_sh_to * vr_to - b_sh_to * vi_to)

    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :cr_fr, :cr_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), cr)
end
"variable: `ci[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_imaginary_on_off(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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
        
        z_branch = _PM.var(pm, _FP.first_id(pm,nw,:PCE_coeff), :z_branch, l)

        ci[(l,i,j)] = z_branch * ((tr * csi_fr + ti * csr_fr + g_sh_fr * vi_fr + b_sh_fr * vr_fr)) / tm^2
        ci[(l,j,i)] = z_branch * (-csi_fr + g_sh_to * vi_to + b_sh_to * vr_to)

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
    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :ci_fr, :ci_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), ci)
end

## load
""
function variable_load_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_load_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_load_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
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


## generator
""
function variable_gen_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# general constraints
## bus
""
function constraint_current_balance(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)

    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

# galerkin projection constraint
## branch
""
function constraint_gp_branch_series_current_magnitude_squared(pm::AbstractIVRModel, n::Int, i, T2, T3)
    cmss  = _PM.var(pm, n, :cmss, i)
    csr = Dict(nw => _PM.var(pm, nw, :csr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    csi = Dict(nw => _PM.var(pm, nw, :csi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * cmss
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (csr[n1] * csr[n2] + csi[n1] * csi[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_branch_series_current_magnitude_squared_contingency(pm::AbstractIVRModel, n::Int, i, T2, T3)
    cmss  = _PM.var(pm, n, :cmss, i)
    csr  = _PM.var(pm, n, :csr, i)
    csi  = _PM.var(pm, n, :csi, i)


    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * cmss
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (csr[n1] * csr[n2] + csi[n1] * csi[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * csr
                                ==
                                0
                    )

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * csi
                                ==
                                0
                    )
end

## generator
""
function constraint_gp_gen_power_real(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    pg  = _PM.var(pm, n, :pg, g)

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)
    
    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * pg
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                (vr[n1] * crg[n2] + vi[n1] * cig[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_gen_power_real_redispatch(pm::AbstractIVRModel, n::Int, i, g, T2, T3)

    pg  = _PM.var(pm, n, :pg, g)

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    pg_base = _PM.var(pm, coeff_idx, :pg, g)

    relax_param = 0.01
    
    # JuMP.@constraint(pm.model,  pg == pg_base)

    JuMP.@constraint(pm.model,  pg >= pg_base * (1 - relax_param))
    JuMP.@constraint(pm.model,  pg <= pg_base * (1 + relax_param))



                    
end

function constraint_gen_redispatch(pm::AbstractIVRModel, g::Int; nw::Int=nw_id_default)
    # display(nw)
    pg  = _PM.var(pm, nw, :pg, g)
    re_pg  = _PM.var(pm, nw, :re_pg, g)

    coeff_idx = _FP.coord(pm, nw, :PCE_coeff)

    pg_base = _PM.var(pm, coeff_idx, :pg, g)

    JuMP.@constraint(pm.model, re_pg >= pg - pg_base) 
    
end

function constraint_voltage_square_redispatch(pm::AbstractIVRModel, i::Int; nw::Int=nw_id_default)
    # display(nw)
    vms  = _PM.var(pm, nw, :vms, i)
    
    coeff_idx = _FP.coord(pm, nw, :PCE_coeff)

    vms_base = _PM.var(pm, coeff_idx, :vms, i)


    JuMP.@constraint(pm.model, vms == vms_base) 
    
end
""
function constraint_gp_gen_power_imaginary(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    qg  = _PM.var(pm, n, :qg, g)

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)
    
    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * qg
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vi[n1] * crg[n2] - vr[n1] * cig[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

## load
""
function constraint_gp_load_power_real(pm::AbstractIVRModel, n::Int, i, l, pd, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    crd = Dict(nw => _PM.var(pm, nw, :crd, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cid = Dict(nw => _PM.var(pm, nw, :cid, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * pd
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vr[n1] * crd[n2] + vi[n1] * cid[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_load_curt_power_real(pm::AbstractIVRModel, n::Int, i, l, pd_curt, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    crd_curt = Dict(nw => _PM.var(pm, nw, :crd_curt, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cid_curt = Dict(nw => _PM.var(pm, nw, :cid_curt, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * pd_curt
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vr[n1] * crd_curt[n2] + vi[n1] * cid_curt[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end
""
function constraint_gp_load_power_imaginary(pm::AbstractIVRModel, n::Int, i, l, qd, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    crd = Dict(nw => _PM.var(pm, nw, :crd, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cid = Dict(nw => _PM.var(pm, nw, :cid, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * qd
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vi[n1] * crd[n2] - vr[n1] * cid[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_load_curt_power_imaginary(pm::AbstractIVRModel, n::Int, i, l, qd_curt, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    crd_curt = Dict(nw => _PM.var(pm, nw, :crd_curt, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cid_curt = Dict(nw => _PM.var(pm, nw, :cid_curt, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * qd_curt
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vi[n1] * crd_curt[n2] - vr[n1] * cid_curt[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * qd_curt
                                == 
                                0                                
                    )
end

# solution
""
function sol_data_model!(pm::AbstractIVRModel, solution::Dict)
    _PM.apply_pm!(_sol_data_model_ivr!, solution)
end

""
function _sol_data_model_ivr!(solution::Dict)
    if haskey(solution, "bus")
        for (i, bus) in solution["bus"]
            if haskey(bus, "vr") && haskey(bus, "vi")
                bus["vm"] = hypot(bus["vr"], bus["vi"])
                bus["va"] = atan(bus["vi"], bus["vr"])
            end
        end
    end
end