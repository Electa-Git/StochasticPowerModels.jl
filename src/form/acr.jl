################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

#sorted_nw_ids(pm) = sort(collect(_PMs.nw_ids(pm)))

# variables
""
function variable_bus_voltage(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    #variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_bus_voltage_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_branch_power(pm::AbstractACRModel;nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_branch_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_branch_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    #variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_branch_current(pm::AbstractACRModel;nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    #_PMs.variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    #_PMs.variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_voltage_drop_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_voltage_drop_img(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end



function variable_gen_power(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_gen_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_gen_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end


""
function variable_bus_voltage_squared(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
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

report && _PMs.sol_component_value(pm, nw, :bus, :vs, _PMs.ids(pm, nw, :bus), vs)
end
# constraints

""
function constraint_voltage_setpoint(pm::AbstractACRModel, i::Int , nw::Int, vm)
    vr  = _PMs.var(pm, nw, :vr, i) 
    vi  = _PMs.var(pm, nw, :vi, i)
   
    if nw ==1
        JuMP.@constraint(pm.model, (vr^2 + vi^2) == vm^2)
    end
end


"""
function constraint_gp_current_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    #old one
    branch = _PMs.ref(pm, n, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    
    css  = _PMs.var(pm, n, :css, i)
    
    g, b = _PMs.calc_branch_y(branch)
    
    vr_fr = Dict(nw => _PMs.var(pm, nw, :vr, f_bus) for nw in _PMs.nw_ids(pm))
    vr_to = Dict(nw => _PMs.var(pm, nw, :vr, t_bus) for nw in _PMs.nw_ids(pm))
    vi_fr = Dict(nw => _PMs.var(pm, nw, :vi, f_bus) for nw in _PMs.nw_ids(pm))
    vi_to = Dict(nw => _PMs.var(pm, nw, :vi, t_bus) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * css
                                ==
                                (g^2+b^2) *
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                ((vr_fr[n1]-vr_to[n1])*(vr_fr[n2]-vr_to[n2])+
                                #2*(vr_fr[n1]-vr_to[n1])*(vi_fr[n2]-vi_to[n2])+
                                (vi_fr[n1]-vi_to[n1])*(vi_fr[n2]-vi_to[n2]))
                                for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end
"""


""
function constraint_branch_voltage(pm::AbstractACRModel,i; nw::Int=nw_id_default)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    
    vbdr  = _PMs.var(pm, nw, :vbdr, i)
    vbdi  = _PMs.var(pm, nw, :vbdi, i)
    #g, b = _PMs.calc_branch_y(branch)
    
    vr_fr = _PMs.var(pm, nw, :vr, f_bus)
    vr_to = _PMs.var(pm, nw, :vr, t_bus)
    vi_fr = _PMs.var(pm, nw, :vi, f_bus)
    vi_to = _PMs.var(pm, nw, :vi, t_bus)

    JuMP.@constraint(pm.model,  vbdr==vr_fr-vr_to)
    JuMP.@constraint(pm.model,  vbdi==vi_fr-vi_to)
                
end

"modified with seperate voltage drop variable"
function constraint_gp_current_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    branch = _PMs.ref(pm, n, :branch, i)
    css  = _PMs.var(pm, n, :css, i)
    
    g, b = _PMs.calc_branch_y(branch)
    
    vbdr = Dict(nw => _PMs.var(pm, nw, :vbdr, i) for nw in _PMs.nw_ids(pm))
    vbdi = Dict(nw => _PMs.var(pm, nw, :vbdi, i) for nw in _PMs.nw_ids(pm))
    

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * css
                                ==
                                (g^2+b^2) *
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                ((vbdr[n1])*(vbdr[n2])+
                                #2*(vr_fr[n1]-vr_to[n1])*(vi_fr[n2]-vi_to[n2])+
                                (vbdi[n1]*vbdi[n2]))
                                for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

function constraint_theta_ref(pm::AbstractACRModel, n::Int, i::Int, vm)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)

    vn = ifelse(n == 1, vm, 0.0)

    JuMP.@constraint(pm.model, vr == vn)
    JuMP.@constraint(pm.model, vi == 0.0)
end

"`vmin <= vm[i] <= vmax` relax the voltage limit so that CC can be applied
not sure if this is required or not"
function constraint_voltage_magnitude_bounds(pm::AbstractACRModel, n::Int, i, vmin, vmax)
    #@assert vmin <= vmax
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)
    vs = _PMs.var(pm, n, :vs, i)
if n==1
    JuMP.@constraint(pm.model, 0.7 <= vr <= 1.3)
    JuMP.@constraint(pm.model, -0.3 <= vi <= 0.3)
    #JuMP.@constraint(pm.model, 0.8 <= vi <= 1.3)
end

end
""
function constraint_gp_bus_voltage_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    vs  = _PMs.var(pm, n, :vs, i)
    vr  = Dict(nw => _PMs.var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm))
    vi  = Dict(nw => _PMs.var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vs 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * vr[n2] + 
                                    #2*vr[n1]*vi[n2]+
                                    vi[n1] * vi[n2]) 
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end


"Expression for branch power from"

function expression_branch_power_ohms_yt_from(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)

    _PMs.var(pm, n, :p)[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
    _PMs.var(pm, n, :q)[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
end


"Expression for branch power to"
function expression_branch_power_ohms_yt_to(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)

    _PMs.var(pm, n, :p)[t_idx] =  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
    _PMs.var(pm, n, :q)[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
end



function constraint_gp_power_branch_from(pm::AbstractACRModel, n::Int,f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
    #complete with shunt and taps
    vr_fr = Dict(nw => _PMs.var(pm, nw, :vr, f_bus) for nw in _PMs.nw_ids(pm))
    vr_to = Dict(nw => _PMs.var(pm, nw, :vr, t_bus) for nw in _PMs.nw_ids(pm))
    vi_fr = Dict(nw => _PMs.var(pm, nw, :vi, f_bus) for nw in _PMs.nw_ids(pm))
    vi_to = Dict(nw => _PMs.var(pm, nw, :vi, t_bus) for nw in _PMs.nw_ids(pm))

    cstr_p = JuMP.@constraint(pm.model,
        _PMs.var(pm, n, :p)[f_idx]* T2.get([n-1,n-1])
        ==
        sum(T3.get([n1-1,n2-1,n-1]) *
        (g+g_fr)/tm^2*(vr_fr[n1]*vr_fr[n2] + vi_fr[n1]*vi_fr[n2]) + (-g*tr+b*ti)/tm^2*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + (-b*tr-g*ti)/tm^2*(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2])
        for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
    )

    cstr_q = JuMP.@constraint(pm.model,
    _PMs.var(pm, n, :q)[f_idx]*T2.get([n-1,n-1])
                    ==
                    sum(T3.get([n1-1,n2-1,n-1])*
                    -(b+b_fr)/tm^2*(vr_fr[n1]*vr_fr[n2] + vi_fr[n1]*vi_fr[n2]) - (-b*tr-g*ti)/tm^2*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + (-g*tr+b*ti)/tm^2*(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2])
                        for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                            )
end

function constraint_gp_power_branch_to(pm::AbstractACRModel, n::Int,f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
    #complete with shunt and taps   
       vr_fr = Dict(nw => _PMs.var(pm, nw, :vr, f_bus) for nw in _PMs.nw_ids(pm))
       vr_to = Dict(nw => _PMs.var(pm, nw, :vr, t_bus) for nw in _PMs.nw_ids(pm))
       vi_fr = Dict(nw => _PMs.var(pm, nw, :vi, f_bus) for nw in _PMs.nw_ids(pm))
       vi_to = Dict(nw => _PMs.var(pm, nw, :vi, t_bus) for nw in _PMs.nw_ids(pm))
   
       cstr_p = JuMP.@constraint(pm.model,
           _PMs.var(pm, n, :p)[t_idx]* T2.get([n-1,n-1])
           ==
           sum(T3.get([n1-1,n2-1,n-1]) *
           (g+g_to)*(vr_to[n1]*vr_to[n2] + vi_to[n1]*vi_to[n2]) + (-g*tr-b*ti)/tm^2*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + (-b*tr+g*ti)/tm^2*(-(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2]))
               for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
       )
   
       cstr_q = JuMP.@constraint(pm.model,
       _PMs.var(pm, n, :q)[t_idx]*T2.get([n-1,n-1])
                       ==
                       sum(T3.get([n1-1,n2-1,n-1])*
                       -(b+b_to)*(vr_to[n1]*vr_to[n2] + vi_to[n1]*vi_to[n2]) - (-b*tr+g*ti)/tm^2*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + (-g*tr-b*ti)/tm^2*(-(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2]))
                           for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                               )
   end
   

function constraint_gp_power_branch_from_simplified(pm::AbstractACRModel, n::Int,f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
## simplified  without shunt and taps
    vr_fr = Dict(nw => _PMs.var(pm, nw, :vr, f_bus) for nw in _PMs.nw_ids(pm))
    vr_to = Dict(nw => _PMs.var(pm, nw, :vr, t_bus) for nw in _PMs.nw_ids(pm))
    vi_fr = Dict(nw => _PMs.var(pm, nw, :vi, f_bus) for nw in _PMs.nw_ids(pm))
    vi_to = Dict(nw => _PMs.var(pm, nw, :vi, t_bus) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,
        _PMs.var(pm, n, :p)[f_idx]* T2.get([n-1,n-1])
        ==
        sum(T3.get([n1-1,n2-1,n-1]) *
        g*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + b*(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2])
        for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
    )

    JuMP.@constraint(pm.model,
    _PMs.var(pm, n, :q)[f_idx]*T2.get([n-1,n-1])
                    ==
                    sum(T3.get([n1-1,n2-1,n-1])*
                    g*(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2]) - (b)*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2])
                        for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                            )
end





function constraint_gp_power_branch_to_simplified(pm::AbstractACRModel, n::Int,f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
#simplified without shunt and taps
    vr_fr = Dict(nw => _PMs.var(pm, nw, :vr, f_bus) for nw in _PMs.nw_ids(pm))
    vr_to = Dict(nw => _PMs.var(pm, nw, :vr, t_bus) for nw in _PMs.nw_ids(pm))
    vi_fr = Dict(nw => _PMs.var(pm, nw, :vi, f_bus) for nw in _PMs.nw_ids(pm))
    vi_to = Dict(nw => _PMs.var(pm, nw, :vi, t_bus) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,
        _PMs.var(pm, n, :p)[t_idx]* T2.get([n-1,n-1])
        ==
        sum(T3.get([n1-1,n2-1,n-1]) *
        g*(vr_to[n1]*vr_fr[n2] + vi_to[n1]*vi_fr[n2]) + b*(vi_to[n1]*vr_fr[n2] - vr_to[n1]*vi_fr[n2])
        #g*(vr_to[n1]*vr_fr[n2] + vi_to[n1]*vi_fr[n2]) + b*(vi_to[n1]*vr_fr[n2] - vr_to[n1]*vi_fr[n2])
            for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
    )

    JuMP.@constraint(pm.model,
    _PMs.var(pm, n, :q)[t_idx]*T2.get([n-1,n-1])
                    ==
                    sum(T3.get([n1-1,n2-1,n-1])*
                    g*(vi_to[n1]*vr_fr[n2] - vr_to[n1]*vi_fr[n2]) - (b)*(vr_to[n1]*vr_fr[n2] + vi_to[n1]*vi_fr[n2])
                    #g*(vi_to[n1]*vr_fr[n2] - vr_to[n1]*vi_fr[n2]) - (b)*(vr_to[n1]*vr_to[n2] + vi_to[n1]*vi_fr[n2])
                        for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                            )
end




function constraint_bus_voltage_squared_cc_limit(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vs  = [_PMs.var(pm, n, :vs, i) for n in sorted_nw_ids(pm)]

    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vs, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vs, mop) <= vmax^2) #this copnstraint of constraint on p screws the problem

    #JuMP.@constraint(pm.model, vmin^2 <= sum([vs[n]*T2.get([n-1,n-1]) for n in _PMs.nw_ids(pm)]))
    #display(JuMP.@constraint(pm.model, sum([vs[n]*T2.get([n-1,n-1]) for n in _PMs.nw_ids(pm)])<=vmax^2))
    #JuMP.@constraint(pm.model, _PCE.var(vs, T2)>=0) #Tillmans paper has this
    #JuMP.@constraint(pm.model, _PCE.mean(vs, mop) >=0) #Tillmans code has this bound

    
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


function constraint_branch_series_current_squared_cc_limit(pm::AbstractACRModel, b, imax, λmax, T2, mop)
    css = [_PMs.var(pm, nw, :css, b) for nw in sorted_nw_ids(pm)]

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(css, mop) <= imax^2)
    #JuMP.@constraint(pm.model,  _PCE.var(css, T2) >= 0) #constraint in Tillemans
     # chance constraint bound
    JuMP.@constraint(pm.model,  _PCE.var(css,T2)
                                <=
                                ((imax^2 - _PCE.mean(css,mop)) / λmax)^2
    )
    JuMP.@constraint(pm.model,  _PCE.var(css,T2)
                                <=
                                (( _PCE.mean(css,mop)) / λmax)^2
                    )
            
end


function constraint_gen_power_real_cc_limit(pm::AbstractACRModel, g, pmin, pmax, λmin, λmax, T2, mop)
    pg  = [_PMs.var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)]
     # bounds on the expectation 
     #JuMP.@constraint(pm.model,  pmin <= _PCE.mean(pg, mop))
     #display(JuMP.@constraint(pm.model,  _PCE.mean(pg, mop) <= pmax))
     #display(JuMP.@constraint(pm.model,   sum([pg[n]*T2.get([n-1,n-1]) for n in _PMs.nw_ids(pm)]) >= pmin)) #given in Tillemans code
     JuMP.@constraint(pm.model,  _PCE.mean(pg, mop) <= pmax)
     # chance constraint bounds
    # JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
    #                             <=
    #                            ((_PCE.mean(pg, mop) - pmin) / λmin)^2
    #               )
     #JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
     #                            <=
     #                            ((pmax - _PCE.mean(pg, mop)) / λmax)^2
     #               )
 end



function constraint_gen_power_imaginary_cc_limit(pm::AbstractACRModel, g, qmin, qmax, λmin, λmax, T2, mop)
    qg  = [_PMs.var(pm, nw, :qg, g) for nw in sorted_nw_ids(pm)]
    # bounds on the expectation 
    JuMP.@constraint(pm.model,  qmin <= _PCE.mean(qg, mop))
    JuMP.@constraint(pm.model,  _PCE.mean(qg, mop) <= qmax)

    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((_PCE.mean(qg,mop) - qmin) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                               <=
                                ((qmax - _PCE.mean(qg,mop)) / λmax)^2
                    )
end

function constraint_power_balance(pm::AbstractACRModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)
    p    = _PMs.get(_PMs.var(pm, n),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = _PMs.get(_PMs.var(pm, n),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = _PMs.get(_PMs.var(pm, n),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = _PMs.get(_PMs.var(pm, n),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = _PMs.get(_PMs.var(pm, n),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = _PMs.get(_PMs.var(pm, n),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = _PMs.get(_PMs.var(pm, n),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = _PMs.get(_PMs.var(pm, n),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = _PMs.get(_PMs.var(pm, n), :p_dc, Dict()); _PMs._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = _PMs.get(_PMs.var(pm, n), :q_dc, Dict()); _PMs._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")


    JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        #+ sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        #+ sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg[g] for g in bus_gens)
        #- sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        #- sum(gs for gs in values(bus_gs))*(vr^2 + vi^2)
    )
    JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        #+ sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        #+ sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        #- sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        #+ sum(bs for bs in values(bus_bs))*(vr^2 + vi^2)
    )

end
