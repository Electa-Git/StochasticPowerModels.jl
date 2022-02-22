################################################################################
#  Copyright 2021, Tom Van Acker, Frederik Geth                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# variables
""
function variable_branch_current(pm::AbstractIVRModel; nw::Int=nw_id_default, aux::Bool=true, bounded::Bool=true, report::Bool=true, aux_fix::Bool=false, kwargs...)   
    _PM.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    _PM.variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
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
    _PM.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    expression_variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    if aux
        variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    else
        if nw == nw_id_default
            variable_branch_series_current_expectation(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
            variable_branch_series_current_variance(pm, nw=nw, bounded=bounded, report=report, aux_fix=aux_fix; kwargs...)
    end end
end

""
function variable_load_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_load_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_load_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_PV_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_PV_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_PV_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_gen_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# general constraints
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


# current balance with PV
""
function constraint_current_balance_with_PV(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs, bus_PV)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)

    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

    crd_pv = _PM.var(pm, n, :crd_pv)
    cid_pv = _PM.var(pm, n, :cid_pv)
    p_size = _PM.var(pm, 1, :p_size)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                + sum(crd_pv[p] for p in bus_PV)  #*p_size[p]
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                + sum(cid_pv[p] for p in bus_PV) #*p_size[p]
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end
# galerkin projection
""
function constraint_gp_branch_series_current_squared(pm::AbstractIVRModel, n::Int, i, T2, T3)
    css  = _PM.var(pm, n, :css, i)
    csr = Dict(nw => _PM.var(pm, nw, :csr, i) for nw in _PM.nw_ids(pm))
    csi = Dict(nw => _PM.var(pm, nw, :csi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * css
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (csr[n1] * csr[n2] + csi[n1] * csi[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

""
function constraint_gp_gen_power_real(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _PM.nw_ids(pm))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _PM.nw_ids(pm))

    pg  = _PM.var(pm, n, :pg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pg
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * crg[n2] + vi[n1] * cig[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

""
function constraint_gp_gen_power_imaginary(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _PM.nw_ids(pm))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _PM.nw_ids(pm))

    qg  = _PM.var(pm, n, :qg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qg
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crg[n2] - vr[n1] * cig[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

""
function constraint_gp_load_power_real(pm::AbstractIVRModel, n::Int, i, l, pd, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    crd = Dict(nw => _PM.var(pm, nw, :crd, l) for nw in _PM.nw_ids(pm))
    cid = Dict(nw => _PM.var(pm, nw, :cid, l) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vr[n1] * crd[n2] + vi[n1] * cid[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

""
function constraint_gp_pv_power_real(pm::AbstractIVRModel, n::Int, i, p, pd, T2, T3, p_size)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    crd_pv = Dict(nw => _PM.var(pm, nw, :crd_pv, p) for nw in _PM.nw_ids(pm))
    cid_pv = Dict(nw => _PM.var(pm, nw, :cid_pv, p) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd * p_size
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vr[n1] * crd_pv[n2] + vi[n1] * cid_pv[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

""
function constraint_gp_load_power_imaginary(pm::AbstractIVRModel, n::Int, i, l, qd, T2, T3)
    vr  = Dict(n => _PM.var(pm, n, :vr, i) for n in _PM.nw_ids(pm))
    vi  = Dict(n => _PM.var(pm, n, :vi, i) for n in _PM.nw_ids(pm))

    crd = Dict(n => _PM.var(pm, n, :crd, l) for n in _PM.nw_ids(pm))
    cid = Dict(n => _PM.var(pm, n, :cid, l) for n in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crd[n2] - vr[n1] * cid[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end


""
function constraint_gp_pv_power_imaginary(pm::AbstractIVRModel, n::Int, i, p, qd, T2, T3, p_size)
    vr  = Dict(n => _PM.var(pm, n, :vr, i) for n in _PM.nw_ids(pm))
    vi  = Dict(n => _PM.var(pm, n, :vi, i) for n in _PM.nw_ids(pm))

    crd_pv = Dict(n => _PM.var(pm, n, :crd_pv, p) for n in _PM.nw_ids(pm))
    cid_pv = Dict(n => _PM.var(pm, n, :cid_pv, p) for n in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd * p_size
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crd_pv[n2] - vr[n1] * cid_pv[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end
# chance constraints
""
function constraint_branch_series_current_cc_limit(pm::AbstractIVRModel, b, cmax, λmax, T2, T4, gs, bs)
    ntws = _PM.nw_ids(pm)

    cse  = _PM.var(pm, nw_id_default, :cse, b)
    csv  = _PM.var(pm, nw_id_default, :csv, b)

    csr  = Dict(n => _PM.var(pm, n, :csr, b) for n in ntws)
    csi  = Dict(n => _PM.var(pm, n, :csi, b) for n in ntws)
    
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