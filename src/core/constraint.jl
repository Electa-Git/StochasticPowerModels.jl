################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# general constraints
## bus
""
function constraint_bus_voltage_ref(pm::AbstractACRModel, n::Int, i::Int)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    vn = ifelse(n == 1, 1.0, 0.0)

    JuMP.@constraint(pm.model, vr == vn)
    JuMP.@constraint(pm.model, vi == 0.0)
end

# galerkin projection constraints
## bus
""
function constraint_gp_bus_voltage_magnitude_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    vms = _PM.var(pm, n, :vms, i)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vms 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

# chance constraints
## bus
""
function constraint_cc_bus_voltage_magnitude_squared(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vms  = [_PM.var(pm, n, :vms, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vms, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vms, mop) <= vmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(vms, T2)
                                <=
                               ((_PCE.mean(vms, mop) - vmin^2) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(vms, T2)
                               <=
                                ((vmax^2 - _PCE.mean(vms, mop)) / λmax)^2
                    )
end

## branch
""
function constraint_cc_branch_series_current_magnitude_squared(pm::AbstractACRModel, b, cmax, λcmax, T2, mop)
    cmss = [_PM.var(pm, nw, :cmss, b) for nw in sorted_nw_ids(pm)]

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(cmss, mop) <= cmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(cmss,T2)
                                <=
                                ((cmax^2 - _PCE.mean(cmss,mop)) / λcmax)^2
                    )
end

## generator
""
function constraint_cc_gen_power_real(pm::AbstractACRModel, g, pmin, pmax, λmin, λmax, T2, mop)
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
function constraint_cc_gen_power_imaginary(pm::AbstractACRModel, g, qmin, qmax, λmin, λmax, T2, mop)
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