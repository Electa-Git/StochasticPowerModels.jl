################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# general constraints
## reference
""
function constraint_bus_voltage_ref(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    constraint_bus_voltage_ref(pm, nw, i)
end

## bus
""
function constraint_current_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(_PM.con(pm, nw), :kcl_cr)
        _PM.con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PM.con(pm, nw), :kcl_ci)
        _PM.con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
end
""
function constraint_power_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus_arcs   = _PM.ref(pm, nw, :bus_arcs, i)
    bus_gens   = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads  = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance(pm, nw, i, bus_arcs, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end

# galerkin projection constraints
## bus
""
function constraint_gp_bus_voltage_magnitude_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_bus_voltage_magnitude_squared(pm, nw, i, T2, T3)
end

## branch
""
function constraint_gp_branch_series_current_magnitude_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_branch_series_current_magnitude_squared(pm, nw, b, T2, T3)
end
""
function constraint_gp_power_branch_to(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
end
""
function constraint_gp_power_branch_from(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
end

## generator
""
function constraint_gp_gen_power(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :gen, g, "gen_bus")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_gen_power_real(pm, nw, i, g, T2, T3)
    constraint_gp_gen_power_imaginary(pm, nw, i, g, T2, T3)
end

## load
""
function constraint_gp_load_power(pm::AbstractPowerModel, l::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :load, l, "load_bus") 

    pd  = _PM.ref(pm, nw, :load, l, "pd")
    qd  = _PM.ref(pm, nw, :load, l, "qd")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_load_power_real(pm, nw, i, l, pd, T2, T3)
    constraint_gp_load_power_imaginary(pm, nw, i, l, qd, T2, T3)
end

# chance constraint limit
## bus
""
function constraint_cc_bus_voltage_magnitude_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :bus, i, "vmin")
    vmax = _PM.ref(pm, nw, :bus, i, "vmax")
    
    λmin = _PM.ref(pm, nw, :bus, i, "λvmin")
    λmax = _PM.ref(pm, nw, :bus, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_bus_voltage_magnitude_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

## branch
""
function constraint_cc_branch_series_current_magnitude_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    cmax = _PM.ref(pm, nw, :branch, b, "cmax")
    λmax = _PM.ref(pm, nw, :branch, b, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_branch_series_current_magnitude_squared(pm, b, cmax, λmax, T2, mop)
end

## generator
""
function constraint_cc_gen_power(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    pmin = _PM.ref(pm, nw, :gen, g, "pmin")
    pmax = _PM.ref(pm, nw, :gen, g, "pmax")
    qmin = _PM.ref(pm, nw, :gen, g, "qmin")
    qmax = _PM.ref(pm, nw, :gen, g, "qmax")

    λpmin = _PM.ref(pm, nw, :gen, g, "λpmin")
    λpmax = _PM.ref(pm, nw, :gen, g, "λpmax")
    λqmin = _PM.ref(pm, nw, :gen, g, "λqmin")
    λqmax = _PM.ref(pm, nw, :gen, g, "λqmax")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_cc_gen_power_real(pm, g, pmin, pmax, λpmin, λpmax, T2, mop)
    constraint_cc_gen_power_imaginary(pm, g, qmin, qmax, λqmin, λqmax, T2, mop)
end