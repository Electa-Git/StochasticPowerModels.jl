################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# variables
""
function variable_bus_voltage(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_bus_voltage_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_branch_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    _PMs.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
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

function variable_gen_power(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_gen_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_gen_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# constraints
""
function constraint_current_balance(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_gs, bus_bs)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    cr =  var(pm, n, :cr)
    ci =  var(pm, n, :ci)
    crdc = var(pm, n, :crdc)
    cidc = var(pm, n, :cidc)

    crd = var(pm, n, :crd)
    cid = var(pm, n, :cid)
    crg = var(pm, n, :crg)
    cig = var(pm, n, :cig)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                + sum(crdc[d] for d in bus_arcs_dc)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

""
function constraint_gp_bus_voltage_squared(pm::AbstractIVRModel, n::Int, i, T2, T3)
    vs  = var(pm, n, :vs, i)
    vr  = Dict(nw => var(pm, nw, :vr, i) for nw in nws(pm))
    vi  = Dict(nw => var(pm, nw, :vi, i) for nw in nws(pm))

    JuMP.@constraint(pm.model,  T2.get([n,n]) * vs 
                                ==
                                sum(T3.get([n1,n2,n]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in nws(pm), n2 in nws(pm))
                    )
end

""
function constraint_gp_branch_series_current_squared(pm::AbstractIVRModel, n::Int, i, T2, T3)
    css  = var(pm, n, :css, i)
    csr = Dict(nw => var(pm, nw, :csr, i) for nw in nws(pm))
    csi = Dict(nw => var(pm, nw, :csi, i) for nw in nws(pm))

    JuMP.@constraint(pm.model,  T2.get([n,n]) * cs
                                ==
                                sum(T3.get([n1,n2,n]) * 
                                    (csr[n1] * csr[n2] + csi[n1] * csi[n2]) 
                                    for n1 in nws(pm), n2 in nws(pm))
                    )
end

""
function constraint_gp_gen_power_real(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => var(pm, nw, :vr, i) for nw in nws(pm))
    vi  = Dict(nw => var(pm, nw, :vi, i) for nw in nws(pm))
    
    crg = Dict(nw => var(pm, nw, :crg, g) for nw in nws(pm))
    cig = Dict(nw => var(pm, nw, :cig, g) for nw in nws(pm))

    pg  = var(pm, n, :pg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n,n]) * pg
                                ==
                                sum(T3.get([n1,n2,n]) * 
                                    (vr[n1] * crg[n2] + vi[n1] * cig[n2])
                                    for n1 in nws(pm), n2 in nws(pm))
                    )
end

""
function constraint_gp_gen_power_imaginary(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => var(pm, nw, :vr, i) for nw in nws(pm))
    vi  = Dict(nw => var(pm, nw, :vi, i) for nw in nws(pm))
    
    crg = Dict(nw => var(pm, nw, :crg, g) for nw in nws(pm))
    cig = Dict(nw => var(pm, nw, :cig, g) for nw in nws(pm))

    qg  = var(pm, n, :qg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n,n]) * qg
                                ==
                                sum(T3.get([n1,n2,n]) *
                                    (vi[n1] * crg[n2] - vr[n1] * cig[n2])
                                    for n1 in nws(pm), n2 in nws(pm))
                    )
end

""
function constraint_gp_load_power_real(pm::AbstractIVRModel, n::Int, i, l, pd, T2, T3)
    vr  = Dict(nw => var(pm, nw, :vr, i) for nw in nws(pm))
    vi  = Dict(nw => var(pm, nw, :vi, i) for nw in nws(pm))

    crd = Dict(nw => var(pm, nw, :crd, l) for nw in nws(pm))
    cid = Dict(nw => var(pm, nw, :cid, l) for nw in nws(pm))

    JuMP.@constraint(pm.model,  T2.get([n,n]) * pd
                                ==
                                sum(T3.get([n1,n2,n]) *
                                    (vr[n1] * crd[n2] + vi[n1] * cid[n2])
                                    for n1 in nws(pm), n2 in nws(pm))
                    )
end

""
function constraint_gp_load_power_imaginary(pm::AbstractIVRModel, n::Int, i, l, qd, T2, T3)
    vr  = Dict(n => var(pm, n, :vr, i) for n in nws(pm))
    vi  = Dict(n => var(pm, n, :vi, i) for n in nws(pm))

    crd = Dict(n => var(pm, n, :crd, l) for n in nws(pm))
    cid = Dict(n => var(pm, n, :cid, l) for n in nws(pm))

    JuMP.@constraint(pm.model,  T2.get([n,n]) * qd
                                ==
                                sum(T3.get([n1,n2,n]) *
                                    (vi[n1] * crd[n2] - vr[n1] * cid[n2])
                                    for n1 in nws(pm), n2 in nws(pm))
                    )
end

""
function constraint_bus_voltage_squared_cc_limit(pm::AbstractIVRModel, i, vmin, vmax, λ, T2, mop)
    vs  = [var(pm, n, :vs, i) for n in nws(pm)]

    JuMP.@constraint(pm.model,  _PCE.var(vs,T2)
                                <=
                                ((_PCE.mean(vs,mop) - vmin^2) / λ)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(vs,T2)
                                <=
                                ((vmax^2 - _PCE.mean(vs,mop)) / λ)^2
                    )
end

""
function constraint_branch_current_cc_limit(pm::AbstractIVRModel, b, imax, λ, T2, mop)
    css  = [var(pm, nw, :css, b) for nw in nws(pm)]

    JuMP.@constraint(pm.model,  _PCE.var(css,T2)
                                <=
                                ((imax^2 - _PCE.mean(css,mop)) / λ)^2
                    )
end

""
function constraint_gen_power_real_cc_limit(pm::AbstractIVRModel, g, pgmin, pgmax, λ, T2, mop)
    pg  = [var(pm, nw, :pg, g) for nw in nws(pm)]

    JuMP.@constraint(pm.model,  _PCE.var(pg,T2)
                                <=
                                ((_PCE.mean(pg,mop) - pgmin) / λ)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(pg,T2)
                                <=
                                ((pgmax - _PCE.mean(pg,mop)) / λ)^2
                    )
end

""
function constraint_gen_power_imaginary_cc_limit(pm::AbstractIVRModel, g, qgmin, qgmax, λ, T2, mop)
    qg  = [var(pm, nw, :qg, g) for nw in nws(pm)]

    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((_PCE.mean(qg,mop) - qgmin) / λ)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((qgmax - _PCE.mean(qg,mop)) / λ)^2
                    )
end