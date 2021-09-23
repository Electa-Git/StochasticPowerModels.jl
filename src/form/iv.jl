################################################################################
#  Copyright 2021, Tom Van Acker, Frederik Geth                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# variables
""
function variable_bus_voltage(pm::AbstractIVRModel; nw::Int=nw_id_default, aux::Bool=true, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false, kwargs...)
    _PMs.variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    if aux
        variable_bus_voltage_squared(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    else
        if nw == nw_id_default
            variable_bus_voltage_expectation(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
            variable_bus_voltage_variance(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    end end
end

""
function variable_branch_current(pm::AbstractIVRModel; nw::Int=nw_id_default, aux::Bool=true, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false, kwargs...)   
    _PMs.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    _PMs.variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    if aux
        variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    else
        if nw == nw_id_default
            variable_branch_series_current_expectation(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
            variable_branch_series_current_variance(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    end end
end

""
function variable_branch_current_reduced(pm::AbstractIVRModel; nw::Int=nw_id_default, aux::Bool=true, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false, kwargs...)
    _PMs.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    expression_variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    if aux
        variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    end
end

""
function variable_load_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_load_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_load_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_gen_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_gen_power(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_gen_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_gen_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# constraints
""
function constraint_bus_voltage_ref(pm::AbstractIVRModel, n::Int, i::Int)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    vn = ifelse(n == 1, 1.0, 0.0)

    JuMP.@constraint(pm.model, vr == vn)
    JuMP.@constraint(pm.model, vi == 0.0)
end

""
function constraint_current_balance(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    cr = _PMs.var(pm, n, :cr)
    ci = _PMs.var(pm, n, :ci)

    crd = _PMs.var(pm, n, :crd)
    cid = _PMs.var(pm, n, :cid)
    crg = _PMs.var(pm, n, :crg)
    cig = _PMs.var(pm, n, :cig)

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

""
function constraint_gp_bus_voltage_squared(pm::AbstractIVRModel, n::Int, i, T2, T3)
    vs  = _PMs.var(pm, n, :vs, i)
    vr  = Dict(nw => _PMs.var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm))
    vi  = Dict(nw => _PMs.var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vs 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

""
function constraint_gp_branch_series_current_squared(pm::AbstractIVRModel, n::Int, i, T2, T3)
    css  = _PMs.var(pm, n, :css, i)
    csr = Dict(nw => _PMs.var(pm, nw, :csr, i) for nw in _PMs.nw_ids(pm))
    csi = Dict(nw => _PMs.var(pm, nw, :csi, i) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * css
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (csr[n1] * csr[n2] + csi[n1] * csi[n2]) 
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

""
function constraint_gp_gen_power_real(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PMs.var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm))
    vi  = Dict(nw => _PMs.var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm))
    
    crg = Dict(nw => _PMs.var(pm, nw, :crg, g) for nw in _PMs.nw_ids(pm))
    cig = Dict(nw => _PMs.var(pm, nw, :cig, g) for nw in _PMs.nw_ids(pm))

    pg  = _PMs.var(pm, n, :pg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pg
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * crg[n2] + vi[n1] * cig[n2])
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

""
function constraint_gp_gen_power_imaginary(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PMs.var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm))
    vi  = Dict(nw => _PMs.var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm))
    
    crg = Dict(nw => _PMs.var(pm, nw, :crg, g) for nw in _PMs.nw_ids(pm))
    cig = Dict(nw => _PMs.var(pm, nw, :cig, g) for nw in _PMs.nw_ids(pm))

    qg  = _PMs.var(pm, n, :qg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qg
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crg[n2] - vr[n1] * cig[n2])
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

""
function constraint_gp_load_power_real(pm::AbstractIVRModel, n::Int, i, l, pd, T2, T3)
    vr  = Dict(nw => _PMs.var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm))
    vi  = Dict(nw => _PMs.var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm))

    crd = Dict(nw => _PMs.var(pm, nw, :crd, l) for nw in _PMs.nw_ids(pm))
    cid = Dict(nw => _PMs.var(pm, nw, :cid, l) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vr[n1] * crd[n2] + vi[n1] * cid[n2])
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

""
function constraint_gp_load_power_imaginary(pm::AbstractIVRModel, n::Int, i, l, qd, T2, T3)
    vr  = Dict(n => _PMs.var(pm, n, :vr, i) for n in _PMs.nw_ids(pm))
    vi  = Dict(n => _PMs.var(pm, n, :vi, i) for n in _PMs.nw_ids(pm))

    crd = Dict(n => _PMs.var(pm, n, :crd, l) for n in _PMs.nw_ids(pm))
    cid = Dict(n => _PMs.var(pm, n, :cid, l) for n in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crd[n2] - vr[n1] * cid[n2])
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

# chance constraints
""
function constraint_bus_voltage_cc_limit(pm::AbstractIVRModel, i, vmin, vmax, λmin, λmax, T2, T4)
    ntws = _PMs.nw_ids(pm)

    ve  = _PMs.var(pm, nw_id_default, :ve, i)
    vv  = _PMs.var(pm, nw_id_default, :vv, i)

    vr  = Dict(n => _PMs.var(pm, n, :vr, i) for n in ntws)
    vi  = Dict(n => _PMs.var(pm, n, :vi, i) for n in ntws)
    
    T22 = Dict(n => T2.get([n-1,n-1]) for n in ntws)
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
    JuMP.@constraint(pm.model,  ve + λmax * vv
                                <=
                                vmax^2
                    )
end

""
function constraint_bus_voltage_squared_cc_limit(pm::AbstractIVRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vs = [_PMs.var(pm, nw, :vs, i) for nw in sorted_nw_ids(pm)]

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
function constraint_branch_series_current_cc_limit(pm::AbstractIVRModel, b, cmax, λmax, T2, T4,gs,bs)
    ntws = _PMs.nw_ids(pm)

    cse  = _PMs.var(pm, nw_id_default, :cse, b)
    csv  = _PMs.var(pm, nw_id_default, :csv, b)

    csr  = Dict(n => _PMs.var(pm, n, :csr, b) for n in ntws)
    csi  = Dict(n => _PMs.var(pm, n, :csi, b) for n in ntws)
    
    T22 = Dict(n => T2.get([n-1,n-1]) for n in ntws)
    T44 = Dict((n1,n2,n3,n4) => T4.get([n1-1,n2-1,n3-1,n4-1]) for n1 in ntws, n2 in ntws, n3 in ntws, n4 in ntws)

    # expectation
    JuMP.@constraint(pm.model, cse == sum((csr[n]^2 + csi[n]^2) * T2.get([n-1,n-1]) for n in ntws))
    # 'variance'
    JuMP.@NLconstraint(pm.model, cse^2 + csv^2 
                                 ==
                                 sum(
                                    (csr[n1] * csr[n2] * csr[n3] * csr[n4] + 
                                     2 * csr[n1] * csr[n2] * csi[n3] * csi[n4] +
                                     csi[n1] * csi[n2] * csi[n3] * csi[n4]
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
function constraint_branch_series_current_squared_cc_limit(pm::AbstractIVRModel, b, imax, λmax, T2, mop)
    css = [_PMs.var(pm, nw, :css, b) for nw in sorted_nw_ids(pm)]

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(css, mop) <= imax^2)
    # chance constraint bound
    JuMP.@constraint(pm.model,  _PCE.var(css,T2)
                                <=
                                ((imax^2 - _PCE.mean(css, mop)) / λmax)^2
                    )
end

""
function constraint_gen_power_real_cc_limit(pm::AbstractIVRModel, g, pmin, pmax, λmin, λmax, T2, mop)
    pg  = [_PMs.var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)]

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
function constraint_gen_power_imaginary_cc_limit(pm::AbstractIVRModel, g, qmin, qmax, λmin, λmax, T2, mop)
    qg  = [_PMs.var(pm, nw, :qg, g) for nw in sorted_nw_ids(pm)]

    # bounds on the expectation 
    JuMP.@constraint(pm.model,  qmin <= _PCE.mean(qg, mop))
    JuMP.@constraint(pm.model,  _PCE.mean(qg, mop) <= qmax)
    #
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((_PCE.mean(qg,mop) - qmin) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((qmax - _PCE.mean(qg,mop)) / λmax)^2
                    )
end