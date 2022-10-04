################################################################################
#  Copyright 2021, Arpan Koirala, Tom Van Acker                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# variables
## branch
""
function variable_branch_power(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_branch_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
function variable_branch_current(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_series_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"" 
function variable_branch_voltage_drop(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_voltage_drop_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_voltage_drop_img(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"variable: `vbdr[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_voltage_drop_real(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vbdr = _PM.var(pm, nw)[:vbdr] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_vbdr",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "vbdr_start", 0.01)
    )
    
    report && _PM.sol_component_value(pm, nw, :branch, :vbdr, _PM.ids(pm, nw, :branch), vbdr)
end
"variable: `vbdi[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_voltage_drop_img(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

    vbdi = _PM.var(pm, nw)[:vbdi] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_vbdi",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "vbdi_start", 0.01)
    )

    report && _PM.sol_component_value(pm, nw, :branch, :vbdi, _PM.ids(pm, nw, :branch), vbdi)
end

# general constraints
## bus
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

## branch
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
## branch
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

## branch
""
function constraint_gp_branch_series_current_magnitude_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    cmss = _PM.var(pm, n, :cmss, i)
    
    branch = _PM.ref(pm, n, :branch, i)
    g, b  = _PM.calc_branch_y(branch)
    
    vbdr = Dict(nw => _PM.var(pm, nw, :vbdr, i) for nw in _PM.nw_ids(pm))
    vbdi = Dict(nw => _PM.var(pm, nw, :vbdi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * cmss
                                ==
                                (g^2 + b^2) * sum(  T3.get([n1-1,n2-1,n-1]) * 
                                                    (vbdr[n1] * vbdr[n2] + vbdi[n1] * vbdi[n2])
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

# solution
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