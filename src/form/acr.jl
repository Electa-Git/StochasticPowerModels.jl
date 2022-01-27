################################################################################
#  Copyright 2021, Arpan Koirala, Tom Van Acker                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# variables
""
function variable_bus_voltage(pm::AbstractACRModel; nw::Int=nw_id_default, aux::Bool=true, aux_fix::Bool=false, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    if aux
        variable_bus_voltage_squared(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    else
        if nw == nw_id_default
            variable_bus_voltage_expectation(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
            variable_bus_voltage_variance(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
        end 
    end
end

""
function variable_gen_power(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_gen_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_branch_power(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_branch_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

"" # this should be renamed to more accurately reflect the variable, i.e., voltage drop
function variable_branch_current(pm::AbstractACRModel; nw::Int=nw_id_default, aux::Bool=true, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false, kwargs...)
    #variable_branch_voltage_drop_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    #variable_branch_voltage_drop_img(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    if aux
        variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    else
        if nw == nw_id_default
            variable_branch_series_current_expectation(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
            variable_branch_series_current_variance(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
        end 
    end
end

# general constraints
""
function constraint_bus_voltage_ref(pm::AbstractACRModel, n::Int, i::Int)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    vn = ifelse(n == 1, 1.0, 0.0)

    JuMP.@constraint(pm.model, vr == vn)
    JuMP.@constraint(pm.model, vi == 0.0)
end

function constraint_power_balance(pm::AbstractACRModel, n::Int, i::Int, bus_arcs, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    p    = _PM.get(_PM.var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = _PM.get(_PM.var(pm, n),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    
    pg   = _PM.get(_PM.var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = _PM.get(_PM.var(pm, n),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")

    JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2)
    )
    JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2)
    )
end

""
function constraint_branch_voltage(pm::AbstractACRModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    
    vbdr  = _PM.var(pm, nw, :vbdr, i)
    vbdi  = _PM.var(pm, nw, :vbdi, i)
    
    vr_fr = _PM.var(pm, nw, :vr, f_bus)
    vr_to = _PM.var(pm, nw, :vr, t_bus)
    vi_fr = _PM.var(pm, nw, :vi, f_bus)
    vi_to = _PM.var(pm, nw, :vi, t_bus)

    JuMP.@constraint(pm.model,  vbdr == (vr_fr-vr_to))
    JuMP.@constraint(pm.model,  vbdi == (vi_fr-vi_to))         
end

# galerkin projection
""
function constraint_gp_bus_voltage_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    vs  = _PM.var(pm, n, :vs, i)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vs 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_power_branch_to(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
    p_to = _PM.var(pm, n, :p, t_idx)
    q_to = _PM.var(pm, n, :q, t_idx)  
    
    vr_fr = Dict(nw => _PM.var(pm, nw, :vr, f_bus) for nw in _PM.nw_ids(pm))
    vr_to = Dict(nw => _PM.var(pm, nw, :vr, t_bus) for nw in _PM.nw_ids(pm))
    vi_fr = Dict(nw => _PM.var(pm, nw, :vi, f_bus) for nw in _PM.nw_ids(pm))
    vi_to = Dict(nw => _PM.var(pm, nw, :vi, t_bus) for nw in _PM.nw_ids(pm))
   
    JuMP.@constraint(pm.model,  p_to * T2.get([n-1,n-1])
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    ((g + g_to) * (vr_to[n1] * vr_to[n2] + vi_to[n1] * vi_to[n2]) + 
                                     (-g * tr - b * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-b * tr + g * ti) / tm^2 * (-(vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2]))
                                    )
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
   
    JuMP.@constraint(pm.model,  q_to * T2.get([n-1,n-1])
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (-(b + b_to) * (vr_to[n1] * vr_to[n2] + vi_to[n1] * vi_to[n2]) - 
                                     (-b * tr + g * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-g * tr - b * ti) / tm^2 * (-(vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2]))
                                    )
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

""
function constraint_gp_power_branch_from(pm::AbstractACRModel, n::Int,f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
    p_fr = _PM.var(pm, n, :p, f_idx)
    q_fr = _PM.var(pm, n, :q, f_idx)

    vr_fr = Dict(nw => _PM.var(pm, nw, :vr, f_bus) for nw in _PM.nw_ids(pm))
    vr_to = Dict(nw => _PM.var(pm, nw, :vr, t_bus) for nw in _PM.nw_ids(pm))
    vi_fr = Dict(nw => _PM.var(pm, nw, :vi, f_bus) for nw in _PM.nw_ids(pm))
    vi_to = Dict(nw => _PM.var(pm, nw, :vi, t_bus) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  p_fr * T2.get([n-1,n-1])
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    ((g + g_fr) / tm^2 * (vr_fr[n1] * vr_fr[n2] + vi_fr[n1] * vi_fr[n2]) + 
                                     (-g * tr + b * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-b * tr - g * ti) / tm^2 * (vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2])
                                    )
                                for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )

    JuMP.@constraint(pm.model,  q_fr * T2.get([n-1,n-1])
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (-(b + b_fr) / tm^2 * (vr_fr[n1] * vr_fr[n2] + vi_fr[n1] * vi_fr[n2]) - 
                                     (-b * tr - g * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-g * tr + b * ti) / tm^2 * (vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2]) 
                                    )
                                for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

""
function constraint_gp_current_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    css = _PM.var(pm, n, :css, i)
    
    branch = _PM.ref(pm, n, :branch, i)
    g, b  = _PM.calc_branch_y(branch)
    
    vbdr = Dict(nw => _PM.var(pm, nw, :vbdr, i) for nw in _PM.nw_ids(pm))
    vbdi = Dict(nw => _PM.var(pm, nw, :vbdi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * css
                                ==
                                (g^2 + b^2) * sum(  T3.get([n1-1,n2-1,n-1]) * 
                                                    (vbdr[n1] * vbdr[n2] + vbdi[n1] * vbdi[n2])
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

# chance constraints
""
function constraint_bus_voltage_cc_limit(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, T4)
    ntws = _PM.nw_ids(pm)

    ve  = _PM.var(pm, nw_id_default, :ve, i)
    vv  = _PM.var(pm, nw_id_default, :vv, i)

    vr  = Dict(n => _PM.var(pm, n, :vr, i) for n in ntws)
    vi  = Dict(n => _PM.var(pm, n, :vi, i) for n in ntws)
    
    T44 = Dict((n1,n2,n3,n4) => T4.get([n1-1,n2-1,n3-1,n4-1]) for n1 in ntws, n2 in ntws, n3 in ntws, n4 in ntws)

    # expectation
    JuMP.@constraint(pm.model, ve == sum((vr[n]^2 + vi[n]^2) * T2.get([n-1,n-1]) for n in ntws))
    # 'variance'
    JuMP.@NLconstraint(pm.model, ve^2 + vv^2 
                                 ==
                                 sum(
                                    (vr[n1] * vr[n2] * vr[n3] * vr[n4] + 
                                     2 * vr[n1] * vr[n2] * vi[n3] * vi[n4] +
                                     vi[n1] * vi[n2] * vi[n3] * vi[n4]
                                    ) *
                                    T44[(n1,n2,n3,n4)]
                                    for n1 in ntws, n2 in ntws, n3 in ntws, n4 in ntws
                                 ) 
                    )  
    # chance constraint bounds
    JuMP.@constraint(pm.model,  vmin^2
                                <=
                                ve - λmin * vv
                    )
    JuMP.@constraint(pm.model,  ve + λmax * vv                                  # not used by Muhlpfordt
                               <=
                               vmax^2
                   )
end

""
function constraint_bus_voltage_squared_cc_limit(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vs  = [_PM.var(pm, n, :vs, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vs, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vs, mop) <= vmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(vs, T2)
                                <=
                               ((_PCE.mean(vs, mop) - vmin^2) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(vs, T2)
                               <=
                                ((vmax^2 - _PCE.mean(vs, mop)) / λmax)^2
                    )
end

""
function constraint_gen_power_real_cc_limit(pm::AbstractACRModel, g, pmin, pmax, λmin, λmax, T2, mop)
    pg  = [_PM.var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)]

     # bounds on the expectation 
     JuMP.@constraint(pm.model,  pmin <= _PCE.mean(pg, mop))
     JuMP.@constraint(pm.model,  _PCE.mean(pg, mop) <= pmax)
     # chance constraint bounds
     JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
                                 <=
                                ((_PCE.mean(pg, mop) - pmin) / λmin)^2
                   )
     JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
                                 <=
                                 ((pmax - _PCE.mean(pg, mop)) / λmax)^2
                   )
end

""
function constraint_gen_power_imaginary_cc_limit(pm::AbstractACRModel, g, qmin, qmax, λmin, λmax, T2, mop)
    qg  = [_PM.var(pm, nw, :qg, g) for nw in sorted_nw_ids(pm)]

    # bounds on the expectation 
    JuMP.@constraint(pm.model,  qmin <= _PCE.mean(qg, mop))
    JuMP.@constraint(pm.model,  _PCE.mean(qg, mop) <= qmax)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((_PCE.mean(qg,mop) - qmin) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                               <=
                                ((qmax - _PCE.mean(qg,mop)) / λmax)^2
                    )
end

""
function constraint_branch_series_current_cc_limit(pm::AbstractACRModel, b, cmax, λmax, T2, T4, gs, bs)
    ntws = _PM.nw_ids(pm)

    cse  = _PM.var(pm, nw_id_default, :cse, b)
    csv  = _PM.var(pm, nw_id_default, :csv, b)

    vbdr = Dict(n => _PM.var(pm, n, :vbdr, b) for n in ntws)
    vbdi = Dict(n => _PM.var(pm, n, :vbdi, b) for n in ntws)
    
    T44 = Dict((n1,n2,n3,n4) => T4.get([n1-1,n2-1,n3-1,n4-1]) for n1 in ntws, n2 in ntws, n3 in ntws, n4 in ntws)

    # expectation
    JuMP.@constraint(pm.model, cse == (gs^2+bs^2) * sum((vbdr[n]^2 + vbdi[n]^2) * T2.get([n-1,n-1]) for n in ntws))
    # 'variance'
    JuMP.@NLconstraint(pm.model, cse^2 + csv^2 
                                 ==
                                 (gs^2+bs^2)^2 * sum(
                                    (vbdr[n1] * vbdr[n2] * vbdr[n3] * vbdr[n4] + 
                                     2 * vbdr[n1] * vbdr[n2] * vbdi[n3] * vbdi[n4] +
                                     vbdi[n1] * vbdi[n2] * vbdi[n3] * vbdi[n4]
                                    ) *
                                    T44[(n1,n2,n3,n4)]
                                    for n1 in ntws, n2 in ntws, n3 in ntws, n4 in ntws
                                 ) 
                    )  
    # chance constraint bounds
    JuMP.@constraint(pm.model,  cse + λmax * csv
                                <=
                                cmax^2
                    )
end

""
function constraint_branch_series_current_squared_cc_limit(pm::AbstractACRModel, b, cmax, λcmax, T2, mop)
    css = [_PM.var(pm, nw, :css, b) for nw in sorted_nw_ids(pm)]
    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(css, mop) <= cmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(css,T2)
                                <=
                                ((cmax^2 - _PCE.mean(css,mop)) / λcmax)^2
                    )
end

""
function sol_data_model!(pm::AbstractACRModel, solution::Dict)
    _PM.apply_pm!(_sol_data_model_acr!, solution)
end


""
function _sol_data_model_acr!(solution::Dict)
    if haskey(solution, "bus")
        for (i, bus) in solution["bus"]
            if haskey(bus, "vr") && haskey(bus, "vi")
                bus["vm"] = hypot(bus["vr"], bus["vi"])
                bus["va"] = atan(bus["vi"], bus["vr"])
            end
        end
    end
end


""
function  expression_branch_voltage_drop_real(pm::AbstractACRModel, n::Int, b::Int, f_bus, t_bus)
    vr_fr = _PM.var(pm, n, :vr, f_bus)
    vr_to = _PM.var(pm, n, :vr, t_bus)
    _PM.var(pm, n, :vbdr)[b]=vr_fr-vr_to
end

""
function  expression_branch_voltage_drop_img(pm::AbstractACRModel, n::Int, b::Int, f_bus, t_bus)
    vi_fr = _PM.var(pm, n, :vi, f_bus)
    vi_to = _PM.var(pm, n, :vi, t_bus)
    _PM.var(pm, n, :vbdi)[b] = (vi_fr-vi_to)
end