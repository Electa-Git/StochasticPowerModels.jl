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

    # vn = ifelse(n == 1, 1.0, 0.0)

    if _FP.is_first_id(pm,n,:PCE_coeff)
        vn = 1.0
    else
        vn = 0.0
    end

    JuMP.@constraint(pm.model, vr == vn)
    JuMP.@constraint(pm.model, vi == 0.0)
end


# galerkin projection constraints
## bus
""
function constraint_gp_bus_voltage_magnitude_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    vms = _PM.var(pm, n, :vms, i)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    # prev_n = _FP.prev_id(pm, n, :PCE_coeff)
    coeff_idx = _FP.coord(pm, n, :PCE_coeff)


    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * vms 
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
    
    
#     JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vms 
#                                 ==
#                                 sum(T3.get([n1-1,n2-1,n-1]) * 
#                                     (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
#                                     for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
#                     )

end

# chance constraints
## bus
""
function constraint_cc_bus_voltage_magnitude_squared(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop, nw)
    # vms  = [_PM.var(pm, n, :vms, i) for n in sorted_nw_ids(pm)]

    vms  = [_PM.var(pm, n, :vms, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]
    

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
function constraint_cc_branch_series_current_magnitude_squared(pm::AbstractACRModel, b, cmax, λcmax, T2, mop, nw)
    # cmss = [_PM.var(pm, nw, :cmss, b) for nw in sorted_nw_ids(pm)]

    cmss = [_PM.var(pm, n, :cmss, b) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]
    

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(cmss, mop) <= cmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(cmss,T2)
                                <=
                                ((cmax^2 - _PCE.mean(cmss,mop)) / λcmax)^2
                    )
end

function constraint_cc_branch_series_current_magnitude_squared_on_off(pm::AbstractACRModel, b, cmax, λcmax, T2, mop, nw)
    # cmss = [_PM.var(pm, nw, :cmss, b) for nw in sorted_nw_ids(pm)]

    cmss = [_PM.var(pm, n, :cmss, b) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]

    z_branch = _PM.var(pm, _FP.first_id(pm, nw, :PCE_coeff), :z_branch, b)
    

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(cmss, mop) <= (z_branch*cmax^2))
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(cmss,T2)
                                <=
                                (((z_branch * cmax^2) - _PCE.mean(cmss,mop)) / λcmax)^2
                    )
end


## generator
""
function constraint_cc_gen_power_real(pm::AbstractACRModel, g, pmin, pmax, λmin, λmax, T2, mop, nw)
    # pg  = [_PM.var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)]

    pg  = [_PM.var(pm, n, :pg, g) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]

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

function constraint_cc_load_curt_power_real(pm::AbstractACRModel, l, pmin, pmax, λmin, λmax, T2, mop, nw)

    pd_curt  = [_PM.var(pm, n, :pd_curt, l) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]

     # bounds on the expectation 
     JuMP.@constraint(pm.model,  pmin <= _PCE.mean(pd_curt, mop))
     JuMP.@constraint(pm.model,  _PCE.mean(pd_curt, mop) <= pmax)
     # chance constraint bounds
     JuMP.@constraint(pm.model,  _PCE.var(pd_curt, T2)
                                 <=
                                ((_PCE.mean(pd_curt, mop) - pmin) / λmin)^2
                   )
     JuMP.@constraint(pm.model,  _PCE.var(pd_curt, T2)
                                 <=
                                 ((pmax - _PCE.mean(pd_curt, mop)) / λmax)^2
                   )
end

function constraint_cc_RES_curt_power(pm::AbstractACRModel, p, pmin, pmax, λmin, λmax, T2, mop, nw)

    p_RES_curt  = [_PM.var(pm, n, :p_RES_curt, p) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]

     # bounds on the expectation 
     JuMP.@constraint(pm.model,  pmin <= _PCE.mean(p_RES_curt, mop))
     JuMP.@constraint(pm.model,  _PCE.mean(p_RES_curt, mop) <= pmax)
     # chance constraint bounds
     JuMP.@constraint(pm.model,  _PCE.var(p_RES_curt, T2)
                                 <=
                                ((_PCE.mean(p_RES_curt, mop) - pmin) / λmin)^2
                   )
     JuMP.@constraint(pm.model,  _PCE.var(p_RES_curt, T2)
                                 <=
                                 ((pmax - _PCE.mean(p_RES_curt, mop)) / λmax)^2
                   )
end

""
function constraint_cc_gen_power_imaginary(pm::AbstractACRModel, g, qmin, qmax, λmin, λmax, T2, mop, nw)
    # qg  = [_PM.var(pm, nw, :qg, g) for nw in sorted_nw_ids(pm)]

    qg  = [_PM.var(pm, n, :qg, g) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]

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

function constraint_current_balance_with_RES(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs, bus_RES)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)
    cidc = _PM.var(pm, n, :cidc)

    iik_r = _PM.var(pm, n, :iik_r)
    iik_i = _PM.var(pm, n, :iik_i)

    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

    crd_RES = _PM.var(pm, n, :crd_RES)
    cid_RES = _PM.var(pm, n, :cid_RES)
    #p_size = _PM.var(pm, 1, :p_size)

    
    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs) + sum(iik_r[c] for c in bus_convs_ac)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                + sum(crd_RES[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs) + sum(iik_i[c] for c in bus_convs_ac)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                + sum(cid_RES[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

function constraint_current_balance_with_RES_curt_load(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs, bus_RES)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)
    cidc = _PM.var(pm, n, :cidc)

    iik_r = _PM.var(pm, n, :iik_r)
    iik_i = _PM.var(pm, n, :iik_i)

    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

    crd_RES = _PM.var(pm, n, :crd_RES)
    cid_RES = _PM.var(pm, n, :cid_RES)
    
    crd_curt = _PM.var(pm, n, :crd_curt)
    cid_curt = _PM.var(pm, n, :cid_curt)

    
    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs) + sum(iik_r[c] for c in bus_convs_ac)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads) + sum(crd_curt[d] for d in bus_loads)
                                + sum(crd_RES[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs) + sum(iik_i[c] for c in bus_convs_ac)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads) + sum(cid_curt[d] for d in bus_loads)
                                + sum(cid_RES[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

function constraint_current_balance_with_RES_curt_RES(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs, bus_RES)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)
    cidc = _PM.var(pm, n, :cidc)

    iik_r = _PM.var(pm, n, :iik_r)
    iik_i = _PM.var(pm, n, :iik_i)

    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

    crd_RES = _PM.var(pm, n, :crd_RES)
    cid_RES = _PM.var(pm, n, :cid_RES)
    
    crd_RES_curt = _PM.var(pm, n, :crd_RES_curt)
    cid_RES_curt = _PM.var(pm, n, :cid_RES_curt)
    

    
    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs) + sum(iik_r[c] for c in bus_convs_ac)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                + sum(crd_RES[p] for p in bus_RES) - sum(crd_RES_curt[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs) + sum(iik_i[c] for c in bus_convs_ac)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                + sum(cid_RES[p] for p in bus_RES) - sum(cid_RES_curt[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

function constraint_current_balance_with_RES_curt_all(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs, bus_RES)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)
    cidc = _PM.var(pm, n, :cidc)

    iik_r = _PM.var(pm, n, :iik_r)
    iik_i = _PM.var(pm, n, :iik_i)

    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

    crd_RES = _PM.var(pm, n, :crd_RES)
    cid_RES = _PM.var(pm, n, :cid_RES)
    
    crd_RES_curt = _PM.var(pm, n, :crd_RES_curt)
    cid_RES_curt = _PM.var(pm, n, :cid_RES_curt)
    
    crd_curt = _PM.var(pm, n, :crd_curt)
    cid_curt = _PM.var(pm, n, :cid_curt)

    
    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs) + sum(iik_r[c] for c in bus_convs_ac)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads) + sum(crd_curt[d] for d in bus_loads)
                                + sum(crd_RES[p] for p in bus_RES) - sum(crd_RES_curt[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs) + sum(iik_i[c] for c in bus_convs_ac)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads) + sum(cid_curt[d] for d in bus_loads)
                                + sum(cid_RES[p] for p in bus_RES) - sum(cid_RES_curt[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

function constraint_current_balance_with_RES_ac(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_gs, bus_bs, bus_RES)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)
    cidc = _PM.var(pm, n, :cidc)


    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

    crd_RES = _PM.var(pm, n, :crd_RES)
    cid_RES = _PM.var(pm, n, :cid_RES)
    #p_size = _PM.var(pm, 1, :p_size)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                + sum(crd_RES[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                + sum(cid_RES[p] for p in bus_RES)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

function constraint_current_balance_dc(pm::_PM.AbstractIVRModel, n::Int, bus_arcs_dcgrid, bus_convs_dc, pd)
    
    igrid_dc = _PM.var(pm, n, :igrid_dc)
    iconv_dc = _PM.var(pm, n, :iconv_dc)


    JuMP.@constraint(pm.model, sum(igrid_dc[a] for a in bus_arcs_dcgrid) + sum(iconv_dc[c] for c in bus_convs_dc) == 0) # deal with pd

end

function constraint_gp_ohms_dc_branch(pm::AbstractACRModel, n::Int, i, T2, T3, f_bus, t_bus, f_idx, t_idx, r, p)

    p_fr  = _PM.var(pm, n,  :p_dcgrid, f_idx)
    p_to  = _PM.var(pm, n,  :p_dcgrid, t_idx)

    #Dict(nw => _PM.var(pm, nw, :vk_r, i) for nw in _PM.nw_ids(pm))
    vmdc_fr = Dict(nw => _PM.var(pm, nw, :vdcm, f_bus) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vmdc_to = Dict(nw => _PM.var(pm, nw, :vdcm, t_bus) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))


    # for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))

    i_dc_fr = Dict(nw => _PM.var(pm, nw, :igrid_dc, f_idx) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    i_dc_to = Dict(nw => _PM.var(pm, nw, :igrid_dc, t_idx) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * p_fr 
                                ==  
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                (vmdc_fr[n1] * i_dc_fr[n2]) 
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )


    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * p_to 
                                ==  
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vmdc_to[n1] * i_dc_to[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )


end




function constraint_ohms_dc_branch(pm::AbstractACRModel, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, p)

    vmdc_fr = _PM.var(pm, n,  :vdcm, f_bus)
    vmdc_to = _PM.var(pm, n,  :vdcm, t_bus)
    i_dc_fr = _PM.var(pm, n,  :igrid_dc, f_idx)
    i_dc_to = _PM.var(pm, n,  :igrid_dc, t_idx)


    if r == 0
        JuMP.@constraint(pm.model, i_dc_fr + i_dc_to == 0)
    else
        JuMP.@constraint(pm.model, vmdc_to ==  vmdc_fr - 1/p * r * i_dc_fr)
        JuMP.@constraint(pm.model, vmdc_fr ==  vmdc_to - 1/p * r * i_dc_to)
    end
end


function constraint_gp_filter_voltage_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    vk_s = _PM.var(pm, n, :vk_s, i)
    vk_r  = Dict(nw => _PM.var(pm, nw, :vk_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vk_i  = Dict(nw => _PM.var(pm, nw, :vk_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    
    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * vk_s 
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (vk_r[n1] * vk_r[n2] + vk_i[n1] * vk_i[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_converter_voltage_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    vc_s = _PM.var(pm, n, :vc_s, i)
    vc_r  = Dict(nw => _PM.var(pm, nw, :vc_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vc_i  = Dict(nw => _PM.var(pm, nw, :vc_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1])  * vc_s 
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (vc_r[n1] * vc_r[n2] + vc_i[n1] * vc_i[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_transformer_current_from_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    iik_s = _PM.var(pm, n, :iik_s, i)
    iik_r  = Dict(nw => _PM.var(pm, nw, :iik_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    iik_i  = Dict(nw => _PM.var(pm, nw, :iik_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * iik_s 
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (iik_r[n1] * iik_r[n2] + iik_i[n1] * iik_i[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_transformer_current_to_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    iki_s = _PM.var(pm, n, :iki_s, i)
    iki_r  = Dict(nw => _PM.var(pm, nw, :iki_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    iki_i  = Dict(nw => _PM.var(pm, nw, :iki_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * iki_s 
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (iki_r[n1] * iki_r[n2] + iki_i[n1] * iki_i[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_reactor_current_from_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    ikc_s = _PM.var(pm, n, :ikc_s, i)
    ikc_r  = Dict(nw => _PM.var(pm, nw, :ikc_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    ikc_i  = Dict(nw => _PM.var(pm, nw, :ikc_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * ikc_s 
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (ikc_r[n1] * ikc_r[n2] + ikc_i[n1] * ikc_i[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_reactor_current_to_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    ick_s = _PM.var(pm, n, :ick_s, i)
    ick_r  = Dict(nw => _PM.var(pm, nw, :ick_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    ick_i  = Dict(nw => _PM.var(pm, nw, :ick_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * ick_s 
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (ick_r[n1] * ick_r[n2] + ick_i[n1] * ick_i[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_converter_current_squared(pm::AbstractACRModel, n::Int, i, T2, T3)

    ic_s = _PM.var(pm, n, :ic_s, i)
    ic_r  = Dict(nw => _PM.var(pm, nw, :ic_r, i) for nw in _PM.nw_ids(pm))
    ic_i  = Dict(nw => _PM.var(pm, nw, :ic_i, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * ic_s 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (ic_r[n1] * ic_r[n2] + ic_i[n1] * ic_i[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_iconv_lin_squared_1(pm::AbstractACRModel, n::Int, i, T2, T3)

    iconv_lin_s = _PM.var(pm, n, :iconv_lin_s, i)
    ic_r  = Dict(nw => _PM.var(pm, nw, :ic_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    ic_i  = Dict(nw => _PM.var(pm, nw, :ic_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * iconv_lin_s
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (ic_r[n1] * ic_r[n2] + ic_i[n1] * ic_i[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end



function constraint_gp_converter_dc_power(pm::AbstractACRModel, n::Int, i, T2, T3, b_idx)

    pconv_dc = _PM.var(pm, n, :pconv_dc, i)
    iconv_dc  = Dict(nw => _PM.var(pm, nw, :iconv_dc, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vdcm  = Dict(nw => _PM.var(pm, nw, :vdcm, b_idx) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    #Eq. (43)
    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * pconv_dc 
                                ==  
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                (vdcm[n1] * iconv_dc[n2]) 
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )

end

function constraint_gp_converter_ac_power(pm::_PM.AbstractIVRModel, n::Int, i::Int, T2, T3)
    vc_r = Dict(nw => _PM.var(pm, nw, :vc_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vc_i = Dict(nw => _PM.var(pm, nw, :vc_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    ic_r = Dict(nw => _PM.var(pm, nw, :ic_r, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    ic_i = Dict(nw => _PM.var(pm, nw, :ic_i, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    pconv_ac = _PM.var(pm, n, :pconv_ac, i)
    qconv_ac = _PM.var(pm, n, :qconv_ac, i)

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)


    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * pconv_ac 
                                ==  
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                (vc_r[n1] * ic_r[n2]) 
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                                +
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                (vc_i[n1] * ic_i[n2]) 
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * qconv_ac 
                                ==  
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                (vc_i[n1] * ic_r[n2]) 
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                                -
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                (vc_r[n1] * ic_i[n2]) 
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )

end



function constraint_gp_converter_losses(pm::_PM.AbstractIVRModel, n::Int, i, T2, T3, a, b, c, plmax)
    iconv_lin = _PM.var(pm, n, :iconv_lin, i)
    iconv_lin_s = _PM.var(pm, n, :iconv_lin_s, i)
    pconv_ac = _PM.var(pm, n, :pconv_ac, i)
    pconv_dc = _PM.var(pm, n, :pconv_dc, i)

    JuMP.@constraint(pm.model, pconv_ac + pconv_dc == a + b * iconv_lin + c * iconv_lin_s)

end




function constraint_gp_iconv_lin_squared_2(pm::AbstractACRModel, n::Int, i, T2, T3)

    iconv_lin_s = _PM.var(pm, n, :iconv_lin_s, i)
    iconv_lin  = Dict(nw => _PM.var(pm, nw, :iconv_lin, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * iconv_lin_s 
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) * 
                                    (iconv_lin[n1] * iconv_lin[n2]) 
                                    for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_RES_power_real(pm::AbstractIVRModel, n::Int, i, p, pd, T2, T3, p_size)
        
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    crd_RES = Dict(nw => _PM.var(pm, nw, :crd_RES, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cid_RES = Dict(nw => _PM.var(pm, nw, :cid_RES, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    p_RES  = _PM.var(pm, n, :p_RES, p)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * pd * p_size
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vr[n1] * crd_RES[n2] + vi[n1] * cid_RES[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * p_RES
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vr[n1] * crd_RES[n2] + vi[n1] * cid_RES[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
        
end

function constraint_gp_RES_curt_power_real(pm::AbstractIVRModel, n::Int, i, p, pd, T2, T3, p_size)
        
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    crd_RES_curt = Dict(nw => _PM.var(pm, nw, :crd_RES_curt, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cid_RES_curt = Dict(nw => _PM.var(pm, nw, :cid_RES_curt, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    p_RES_curt  = _PM.var(pm, n, :p_RES_curt, p)
    p_RES  = _PM.var(pm, n, :p_RES, p)

    JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * p_RES_curt
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vr[n1] * crd_RES_curt[n2] + vi[n1] * cid_RES_curt[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )

    # JuMP.@constraint(pm.model,  T2.get([coeff_idx-1,coeff_idx-1]) * p_RES_curt
    #                             >=
    #                             0
    #                 )   
        
end

function constraint_gp_RES_power_imaginary(pm::AbstractIVRModel, n::Int, i, p, qd, T2, T3, q_size)
    
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    crd_RES = Dict(nw => _PM.var(pm, nw, :crd_RES, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cid_RES = Dict(nw => _PM.var(pm, nw, :cid_RES, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    q_RES  = _PM.var(pm, n, :q_RES, p)

    JuMP.@constraint(pm.model, T2.get([coeff_idx-1,coeff_idx-1]) * qd * q_size
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vi[n1] * crd_RES[n2] - vr[n1] * cid_RES[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )

    JuMP.@constraint(pm.model, T2.get([coeff_idx-1,coeff_idx-1]) * q_RES
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vi[n1] * crd_RES[n2] - vr[n1] * cid_RES[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )
end

function constraint_gp_RES_curt_power_imaginary(pm::AbstractIVRModel, n::Int, i, p, qd, T2, T3, q_size)
    
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    crd_RES_curt = Dict(nw => _PM.var(pm, nw, :crd_RES_curt, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
    cid_RES_curt = Dict(nw => _PM.var(pm, nw, :cid_RES_curt, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

    coeff_idx = _FP.coord(pm, n, :PCE_coeff)

    q_RES_curt  = _PM.var(pm, n, :q_RES_curt, p)

    JuMP.@constraint(pm.model, T2.get([coeff_idx-1,coeff_idx-1]) * q_RES_curt
                                ==
                                sum(T3.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n2, :PCE_coeff)-1, coeff_idx-1]) *
                                (vi[n1] * crd_RES_curt[n2] - vr[n1] * cid_RES_curt[n2])
                                for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)), n2 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    )

    JuMP.@constraint(pm.model, T2.get([coeff_idx-1,coeff_idx-1]) * q_RES_curt
                                ==
                                0
                    )
end

function constraint_conv_transformer(pm::_PM.AbstractIVRModel, n::Int, i::Int, rtf, xtf, acbus, tm, transformer)
    vi_r = _PM.var(pm, n, :vr, acbus)
    vi_i = _PM.var(pm, n, :vi, acbus)
    vk_r = _PM.var(pm, n, :vk_r, i)
    vk_i = _PM.var(pm, n, :vk_i, i)

    iik_r = _PM.var(pm, n, :iik_r, i)
    iik_i = _PM.var(pm, n, :iik_i, i)
    iki_r = _PM.var(pm, n, :iki_r, i)
    iki_i = _PM.var(pm, n, :iki_i, i)

    #TODO add transformation ratio.....
    if transformer
        JuMP.@constraint(pm.model, vk_r == vi_r - rtf * iik_r + xtf * iik_i) #(24)
        JuMP.@constraint(pm.model, vk_i == vi_i - rtf * iik_i - xtf * iik_r) #(25)
        JuMP.@constraint(pm.model, vi_r == vk_r - rtf * iki_r + xtf * iki_i) #reverse
        JuMP.@constraint(pm.model, vi_i == vk_i - rtf * iki_i - xtf * iki_r) #reverse
        
    else
        JuMP.@constraint(pm.model, vk_r == vi_r)
        JuMP.@constraint(pm.model, vk_i == vi_i)
        JuMP.@constraint(pm.model, iik_r + iki_r == 0)
        JuMP.@constraint(pm.model, iik_i + iki_i == 0)
    end

end

function constraint_conv_reactor(pm::_PM.AbstractIVRModel, n::Int, i::Int, rc, xc, reactor)
    vk_r = _PM.var(pm, n, :vk_r, i)
    vk_i = _PM.var(pm, n, :vk_i, i)
    vc_r = _PM.var(pm, n, :vc_r, i)
    vc_i = _PM.var(pm, n, :vc_i, i)

    ikc_r = _PM.var(pm, n, :ikc_r, i)
    ikc_i = _PM.var(pm, n, :ikc_i, i)
    ick_r = _PM.var(pm, n, :ick_r, i)
    ick_i = _PM.var(pm, n, :ick_i, i)
    ic_r = _PM.var(pm, n, :ic_r, i)
    ic_i = _PM.var(pm, n, :ic_i, i)

    JuMP.@constraint(pm.model, ick_r + ic_r == 0) #(20)
    JuMP.@constraint(pm.model, ick_i + ic_i == 0) #(21)

    if reactor
        JuMP.@constraint(pm.model, vc_r == vk_r - rc * ikc_r + xc * ikc_i) #(28)
        JuMP.@constraint(pm.model, vc_i == vk_i - rc * ikc_i - xc * ikc_r) #(29)
        JuMP.@constraint(pm.model, vk_r == vc_r - rc * ick_r + xc * ick_i) #reverse
        JuMP.@constraint(pm.model, vk_i == vc_i - rc * ick_i - xc * ick_r) #reverse
    else
        JuMP.@constraint(pm.model, vk_r == vc_r)
        JuMP.@constraint(pm.model, vk_i == vc_i)
        JuMP.@constraint(pm.model, ikc_r + ick_r == 0)
        JuMP.@constraint(pm.model, ikc_i + ick_i == 0)
    end
end

function constraint_conv_filter(pm::_PM.AbstractIVRModel, n::Int, i::Int, bv, filter)
    iki_r = _PM.var(pm, n, :iki_r, i)
    iki_i = _PM.var(pm, n, :iki_i, i)
    ikc_r = _PM.var(pm, n, :ikc_r, i)
    ikc_i = _PM.var(pm, n, :ikc_i, i)

    vk_r = _PM.var(pm, n, :vk_r, i)
    vk_i = _PM.var(pm, n, :vk_i, i)

    JuMP.@constraint(pm.model,   iki_r + ikc_r + bv * filter * vk_i == 0)
    JuMP.@constraint(pm.model,   iki_i + ikc_i - bv * filter * vk_r == 0)
end

function constraint_cc_filter_voltage_squared(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop, nw)
    # vk_s  = [_PM.var(pm, n, :vk_s, i) for n in sorted_nw_ids(pm)]

    vk_s  = [_PM.var(pm, n, :vk_s, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vk_s, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vk_s, mop) <= vmax^2)
    # chance constraint bounds
    vk_s_min_cc = JuMP.@constraint(pm.model,  _PCE.var(vk_s, T2)
                                <=
                               ((_PCE.mean(vk_s, mop) - vmin^2) / λmin)^2
                    )
    vk_s_max_cc = JuMP.@constraint(pm.model,  _PCE.var(vk_s, T2)
                               <=
                                ((vmax^2 - _PCE.mean(vk_s, mop)) / λmax)^2
                    )

    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_filter_voltage_max] = vk_s_max_cc
        _PM.sol(pm, 1, :convdc, i)[:dual_filter_voltage_min] = vk_s_min_cc
    end
end


function constraint_cc_converter_voltage_squared(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop, nw)
    # vc_s  = [_PM.var(pm, n, :vc_s, i) for n in sorted_nw_ids(pm)]

    vc_s  = [_PM.var(pm, n, :vc_s, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vc_s, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vc_s, mop) <= vmax^2)
    # chance constraint bounds
    vc_s_min_cc = JuMP.@constraint(pm.model,  _PCE.var(vc_s, T2)
                                <=
                               ((_PCE.mean(vc_s, mop) - vmin^2) / λmin)^2
                    )
    vc_s_max_cc = JuMP.@constraint(pm.model,  _PCE.var(vc_s, T2)
                               <=
                                ((vmax^2 - _PCE.mean(vc_s, mop)) / λmax)^2
                    )

    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_converter_voltage_max] = vc_s_max_cc
        _PM.sol(pm, 1, :convdc, i)[:dual_converter_voltage_min] = vc_s_min_cc
    end
end

#Eq. (33)
function constraint_cc_transformer_current_from_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop, nw)
    # iik_s  = [_PM.var(pm, n, :iik_s, i) for n in sorted_nw_ids(pm)]

    iik_s  = [_PM.var(pm, n, :iik_s, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iik_s, mop) <= Imax^2)

    # chance constraint bounds
    iik_s_max_cc = JuMP.@constraint(pm.model,  _PCE.var(iik_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(iik_s, mop)) / λmax)^2
                    )

    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_transformer_current_from_max] = iik_s_max_cc
    end
end

#Eq. (34)
function constraint_cc_transformer_current_to_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop, nw)
    # iki_s  = [_PM.var(pm, n, :iki_s, i) for n in sorted_nw_ids(pm)]

    iki_s  = [_PM.var(pm, n, :iki_s, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iki_s, mop) <= Imax^2)

    # chance constraint bounds
    iki_s_max_cc = JuMP.@constraint(pm.model,  _PCE.var(iki_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(iki_s, mop)) / λmax)^2
                    )

    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_transformer_current_to_max] = iki_s_max_cc
    end
end

#Eq. (35)
function constraint_cc_reactor_current_from_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop, nw)
    # ikc_s  = [_PM.var(pm, n, :ikc_s, i) for n in sorted_nw_ids(pm)]

    ikc_s  = [_PM.var(pm, n, :ikc_s, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(ikc_s, mop) <= Imax^2)

    # chance constraint bounds
    ikc_s_max_cc = JuMP.@constraint(pm.model,  _PCE.var(ikc_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(ikc_s, mop)) / λmax)^2
                    )

    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_transformer_current_from_max] = ikc_s_max_cc
    end
end

#Eq. (36)
function constraint_cc_reactor_current_to_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop, nw)
    # ick_s  = [_PM.var(pm, n, :ick_s, i) for n in sorted_nw_ids(pm)]

    ick_s  = [_PM.var(pm, n, :ick_s, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(ick_s, mop) <= Imax^2)

    # chance constraint bounds
    ick_s_max_cc = JuMP.@constraint(pm.model,  _PCE.var(ick_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(ick_s, mop)) / λmax)^2
                    )

    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_transformer_current_to_max] = ick_s_max_cc
    end
end

function constraint_cc_converter_current_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop)
    iconv_lin_s  = [_PM.var(pm, n, :iconv_lin_s, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iconv_lin_s, mop) <= Imax^2)

    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(iconv_lin_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(iconv_lin_s, mop)) / λmax)^2
                    )
end


function constraint_cc_dc_branch_current(pm::AbstractACRModel, i, Imax, Imin, λmax, λmin, f_idx, t_idx, T2, mop, nw)
    # i_dc_fr = [_PM.var(pm, n, :igrid_dc, f_idx) for n in sorted_nw_ids(pm)]
    # i_dc_to = [_PM.var(pm, n, :igrid_dc, t_idx) for n in sorted_nw_ids(pm)]

    i_dc_fr = [_PM.var(pm, n, :igrid_dc, f_idx) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]
    i_dc_to = [_PM.var(pm, n, :igrid_dc, t_idx) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]


    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(i_dc_fr, mop) <= Imax)

    JuMP.@constraint(pm.model, Imin <= _PCE.mean(i_dc_fr, mop))

    # chance constraint bounds
    i_dc_fr_max_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_fr, T2)
                               <=
                                ((Imax - _PCE.mean(i_dc_fr, mop)) / λmax)^2
                    )

    
    i_dc_fr_min_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_fr, T2)
                                <=
                               ((_PCE.mean(i_dc_fr, mop) - Imin) / λmin)^2
    )
    

    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(i_dc_to, mop) <= Imax)
    JuMP.@constraint(pm.model, Imin <= _PCE.mean(i_dc_to, mop))

    # chance constraint bounds
    i_dc_to_max_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_to, T2)
                               <=
                                ((Imax - _PCE.mean(i_dc_to, mop)) / λmax)^2
                    )
    
    i_dc_to_min_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_to, T2)
                                <=
                               ((_PCE.mean(i_dc_to, mop) - Imin) / λmin)^2
    )
    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_from_min] = i_dc_fr_min_cc
        _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_from_max] = i_dc_fr_max_cc
        _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_to_min] = i_dc_to_min_cc
        _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_to_max] = i_dc_to_max_cc
    end
end

function constraint_cc_dc_branch_current_on_off(pm::AbstractACRModel, i, Imax, Imin, λmax, λmin, f_idx, t_idx, T2, mop, nw)
    # i_dc_fr = [_PM.var(pm, n, :igrid_dc, f_idx) for n in sorted_nw_ids(pm)]
    # i_dc_to = [_PM.var(pm, n, :igrid_dc, t_idx) for n in sorted_nw_ids(pm)]

    i_dc_fr = [_PM.var(pm, n, :igrid_dc, f_idx) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]
    i_dc_to = [_PM.var(pm, n, :igrid_dc, t_idx) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]

    z_branch_dc = _PM.var(pm, _FP.first_id(pm, nw, :PCE_coeff), :igrid_dc, t_idx)

    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(i_dc_fr, mop) <= z_branch_dc*Imax)

    JuMP.@constraint(pm.model, z_branch_dc*Imin <= _PCE.mean(i_dc_fr, mop))

    # chance constraint bounds
    i_dc_fr_max_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_fr, T2)
                               <=
                                ((z_branch_dc*Imax - _PCE.mean(i_dc_fr, mop)) / λmax)^2
                    )

    
    i_dc_fr_min_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_fr, T2)
                                <=
                               ((_PCE.mean(i_dc_fr, mop) - z_branch_dc*Imin) / λmin)^2
    )
    

    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(i_dc_to, mop) <= z_branch_dc*Imax)
    JuMP.@constraint(pm.model, z_branch_dc*Imin <= _PCE.mean(i_dc_to, mop))

    # chance constraint bounds
    i_dc_to_max_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_to, T2)
                               <=
                                ((z_branch_dc*Imax - _PCE.mean(i_dc_to, mop)) / λmax)^2
                    )
    
    i_dc_to_min_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_to, T2)
                                <=
                               ((_PCE.mean(i_dc_to, mop) - z_branch_dc*Imin) / λmin)^2
    )
    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_from_min] = i_dc_fr_min_cc
        _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_from_max] = i_dc_fr_max_cc
        _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_to_min] = i_dc_to_min_cc
        _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_to_max] = i_dc_to_max_cc
    end
end

function constraint_cc_iconv_lin_squared(pm::AbstractACRModel, i, Imax, λmax, T2, mop, nw)
    # iconv_lin_s  = [_PM.var(pm, n, :iconv_lin_s, i) for n in sorted_nw_ids(pm)]

    iconv_lin_s  = [_PM.var(pm, n, :iconv_lin_s, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]    

    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iconv_lin_s, mop) <= Imax^2)

    # chance constraint bounds
    iconv_lin_s_max_cc = JuMP.@constraint(pm.model,  _PCE.var(iconv_lin_s, T2)
                               <=
                                ((Imax^2 - _PCE.mean(iconv_lin_s, mop)) / λmax)^2
                    )
    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_iconv_lin_max] = iconv_lin_s_max_cc
    end
end

function constraint_cc_iconv_lin(pm::AbstractACRModel, i, Imax, Imin, λmax, λmin, T2, mop, nw)
    # iconv_lin  = [_PM.var(pm, n, :iconv_lin, i) for n in sorted_nw_ids(pm)]

    iconv_lin  = [_PM.var(pm, n, :iconv_lin, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iconv_lin, mop) <= Imax)
    JuMP.@constraint(pm.model, Imin <= _PCE.mean(iconv_lin, mop))
    # chance constraint bounds
    iconv_lin_max_cc = JuMP.@constraint(pm.model,  _PCE.var(iconv_lin, T2)
                               <=
                                ((Imax - _PCE.mean(iconv_lin, mop)) / λmax)^2
                    )
    iconv_lin_min_cc = JuMP.@constraint(pm.model,  _PCE.var(iconv_lin, T2)
                                <=
                                ((_PCE.mean(iconv_lin, mop) - Imin) / λmin)^2
                    )
    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_iconv_lin_max] = iconv_lin_max_cc
        _PM.sol(pm, 1, :convdc, i)[:dual_iconv_lin_min] = iconv_lin_min_cc
    end
end


function constraint_cc_conv_ac_power_real(pm, i, Pacmin, Pacmax, λmin, λmax, T2, mop, nw)
    # pconv_ac  = [_PM.var(pm, n, :pconv_ac, i) for n in sorted_nw_ids(pm)]

    pconv_ac  = [_PM.var(pm, n, :pconv_ac, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 


    # bounds on the expectation
    JuMP.@constraint(pm.model, Pacmin <= _PCE.mean(pconv_ac, mop))
    JuMP.@constraint(pm.model, _PCE.mean(pconv_ac, mop) <= Pacmax)

    # chance constraint bounds
    pconv_ac_min_cc = JuMP.@constraint(pm.model,  _PCE.var(pconv_ac, T2)
                            <=
                            ((_PCE.mean(pconv_ac, mop) - Pacmin) / λmin)^2
                    )

    pconv_ac_max_cc = JuMP.@constraint(pm.model,  _PCE.var(pconv_ac, T2)
                            <=
                            ((Pacmax - _PCE.mean(pconv_ac, mop)) / λmax)^2
                    )
    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_conv_ac_power_real_max] = pconv_ac_max_cc
        _PM.sol(pm, 1, :convdc, i)[:dual_conv_ac_power_real_min] = pconv_ac_min_cc
    end
end

function constraint_cc_conv_ac_power_imaginary(pm, i, Qacmin, Qacmax, λmin, λmax, T2, mop, nw)
    # qconv_ac  = [_PM.var(pm, n, :qconv_ac, i) for n in sorted_nw_ids(pm)]

    qconv_ac  = [_PM.var(pm, n, :qconv_ac, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 

    # bounds on the expectation
    JuMP.@constraint(pm.model, Qacmin <= _PCE.mean(qconv_ac, mop))
    JuMP.@constraint(pm.model, _PCE.mean(qconv_ac, mop) <= Qacmax)

    # chance constraint bounds
    qconv_ac_min_cc = JuMP.@constraint(pm.model,  _PCE.var(qconv_ac, T2)
                            <=
                            ((_PCE.mean(qconv_ac, mop) - Qacmin) / λmin)^2
                    )

    qconv_ac_max_cc = JuMP.@constraint(pm.model,  _PCE.var(qconv_ac, T2)
                            <=
                            ((Qacmax - _PCE.mean(qconv_ac, mop)) / λmax)^2
                    )

    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_conv_ac_power_imaginary_max] = qconv_ac_max_cc
        _PM.sol(pm, 1, :convdc, i)[:dual_conv_ac_power_imaginary_min] = qconv_ac_min_cc
    end

end

function constraint_cc_conv_dc_power(pm, i, Pdcmin, Pdcmax, λmin, λmax, T2, mop, nw)
    # pconv_dc  = [_PM.var(pm, n, :pconv_dc, i) for n in sorted_nw_ids(pm)]

    pconv_dc  = [_PM.var(pm, n, :pconv_dc, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 


    # bounds on the expectation
    JuMP.@constraint(pm.model, Pdcmin <= _PCE.mean(pconv_dc, mop))
    JuMP.@constraint(pm.model, _PCE.mean(pconv_dc, mop) <= Pdcmax)

    # chance constraint bounds
    pconv_dc_min_cc = JuMP.@constraint(pm.model,  _PCE.var(pconv_dc, T2)
                            <=
                            ((_PCE.mean(pconv_dc, mop) - Pdcmin) / λmin)^2
                    )

    pconv_dc_max_cc = JuMP.@constraint(pm.model,  _PCE.var(pconv_dc, T2)
                            <=
                            ((Pdcmax - _PCE.mean(pconv_dc, mop)) / λmax)^2
                    )
    
    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_conv_dc_power_max] = pconv_dc_max_cc
        _PM.sol(pm, 1, :convdc, i)[:dual_conv_dc_power_min] = pconv_dc_min_cc
    end
end




function constraint_cc_converter_dc_current(pm::AbstractACRModel, i, Imax, Imin, λmax, λmin, T2, mop, nw)
    # iconv_dc = [_PM.var(pm, n, :iconv_dc, i) for n in sorted_nw_ids(pm)]

    iconv_dc = [_PM.var(pm, n, :iconv_dc, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))] 

    # bounds on the expectation
    JuMP.@constraint(pm.model, _PCE.mean(iconv_dc, mop) <= Imax)

    JuMP.@constraint(pm.model, Imin <= _PCE.mean(iconv_dc, mop))

    # chance constraint bounds
    iconv_dc_max_cc = JuMP.@constraint(pm.model,  _PCE.var(iconv_dc, T2)
                               <=
                                ((Imax - _PCE.mean(iconv_dc, mop)) / λmax)^2
                    )

    
    iconv_dc_min_cc = JuMP.@constraint(pm.model,  _PCE.var(iconv_dc, T2)
                                <=
                               ((_PCE.mean(iconv_dc, mop) - Imin) / λmin)^2
    )
    
    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :convdc, i)[:dual_conv_dc_current_max] = iconv_dc_max_cc
        _PM.sol(pm, 1, :convdc, i)[:dual_conv_dc_current_min] = iconv_dc_min_cc
    end

end


function constraint_cc_conv_voltage_magnitude(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop, nw)
    
    # vdcm  = [_PM.var(pm, n, :vdcm, i) for n in sorted_nw_ids(pm)]

    vdcm  = [_PM.var(pm, n, :vdcm, i) for n in _FP.similar_ids(pm, nw; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin <= _PCE.mean(vdcm, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vdcm, mop) <= vmax)
    # chance constraint bounds
    
    vdcm_min_cc = JuMP.@constraint(pm.model,  _PCE.var(vdcm, T2)
                                <=
                               ((_PCE.mean(vdcm, mop) - vmin) / λmin)^2
                    )
    vdcm_max_cc = JuMP.@constraint(pm.model,  _PCE.var(vdcm, T2)
                               <=
                                ((vmax - _PCE.mean(vdcm, mop)) / λmax)^2
                    )

    # JuMP.@constraint(pm.model,  _PCE.var(vdcm, T2)
    #                             <=
    #                             (0.00001)
    #                 )
    
    if _IM.report_duals(pm)
        _PM.sol(pm, 1, :busdc, i)[:dual_conv_voltage_min] = vdcm_min_cc
        _PM.sol(pm, 1, :busdc, i)[:dual_conv_voltage_max] = vdcm_max_cc
    end

end



# function constraint_cc_dc_branch_current_on_off(pm::AbstractACRModel, i, Imax, Imin, λmax, λmin, f_idx, t_idx, T2, mop)
#     i_dc_fr_on_off = [_PM.var(pm, n, :igrid_dc_on_off, f_idx) for n in sorted_nw_ids(pm)]
#     i_dc_to_on_off = [_PM.var(pm, n, :igrid_dc_on_off, t_idx) for n in sorted_nw_ids(pm)]

#     # bounds on the expectation
#     JuMP.@constraint(pm.model, _PCE.mean(i_dc_fr_on_off, mop) <= Imax)

#     JuMP.@constraint(pm.model, Imin <= _PCE.mean(i_dc_fr_on_off, mop))

#     # chance constraint bounds
#     i_dc_fr_max_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_fr_on_off, T2)
#                                <=
#                                 ((Imax - _PCE.mean(i_dc_fr_on_off, mop)) / λmax)^2
#                     )

    
#     i_dc_fr_min_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_fr_on_off, T2)
#                                 <=
#                                ((_PCE.mean(i_dc_fr_on_off, mop) - Imin) / λmin)^2
#     )
    

#     # bounds on the expectation
#     JuMP.@constraint(pm.model, _PCE.mean(i_dc_to_on_off, mop) <= Imax)
#     JuMP.@constraint(pm.model, Imin <= _PCE.mean(i_dc_to_on_off, mop))

#     # chance constraint bounds
#     i_dc_to_max_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_to_on_off, T2)
#                                <=
#                                 ((Imax - _PCE.mean(i_dc_to_on_off, mop)) / λmax)^2
#                     )
    
#     i_dc_to_min_cc = JuMP.@constraint(pm.model,  _PCE.var(i_dc_to_on_off, T2)
#                                 <=
#                                ((_PCE.mean(i_dc_to_on_off, mop) - Imin) / λmin)^2
#     )
#     if _IM.report_duals(pm)
#         _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_from_min] = i_dc_fr_min_cc
#         _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_from_max] = i_dc_fr_max_cc
#         _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_to_min] = i_dc_to_min_cc
#         _PM.sol(pm, 1, :branchdc, i)[:dual_dc_brach_current_to_max] = i_dc_to_max_cc
#     end
# end


function constraint_gp_dc_branch_indicator(pm::AbstractIVRModel, i::Int; nw::Int=nw_id_default)
    err = 1e-6

    z_branch_dc = _PM.var(pm, _FP.first_id(pm, nw, :PCE_coeff), :z_branch_dc, i)

    JuMP.@constraint(pm.model, z_branch_dc * (1 - z_branch_dc) <= err)

    # JuMP.fix(z_branch_dc[1], 1; force=true)

end

function constraint_gp_ac_branch_indicator(pm::AbstractIVRModel, i::Int; nw::Int=nw_id_default)
    err = 1e-6

    z_branch = _PM.var(pm, _FP.first_id(pm, nw, :PCE_coeff), :z_branch, i)

    JuMP.@constraint(pm.model, z_branch * (1 - z_branch) <= err)

    # JuMP.fix(z_branch_dc[1], 1; force=true)

end