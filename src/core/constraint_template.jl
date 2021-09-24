################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# reference voltage 
""
function constraint_bus_voltage_ref(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    constraint_bus_voltage_ref(pm, nw, i)
end

function constraint_theta_ref(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = _PMs.ref(pm, nw, :bus, i)

    constraint_theta_ref(pm, nw, i, bus["vm"])
end

""
function constraint_current_limit(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "c_rating_a")
        constraint_current_limit(pm, nw, f_idx, branch["c_rating_a"])
    end
end

""
function constraint_voltage_magnitude_bounds(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = _PMs.ref(pm, nw, :bus, i)
    constraint_voltage_magnitude_bounds(pm, nw, i, bus["vmin"], bus["vmax"])
end

""
function constraint_voltage_setpoint(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = _PMs.ref(pm, nw, :bus, i)

    constraint_voltage_setpoint(pm, i, nw, bus["vm"],bus["bus_type"])
end

# current balance
""
function constraint_current_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(_PMs.con(pm, nw), :kcl_cr)
        _PMs.con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PMs.con(pm, nw), :kcl_ci)
        _PMs.con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PMs.ref(pm, nw, :bus, i)
    bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
end

# galerkin projection
""
function constraint_gp_bus_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_bus_voltage_squared(pm, nw, i, T2, T3)
end

""
function constraint_gp_branch_series_current_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_branch_series_current_squared(pm, nw, b, T2, T3)
end

""
function constraint_gp_branch_series_current_squared_new(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T4  = pm.data["T4"]

    constraint_gp_branch_series_current_squared_new(pm, nw, b, T2, T4)
end

""
function constraint_gp_current_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_current_squared(pm, nw, b, T2, T3)
end

""
function constraint_gp_current(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_branch_series_current_squared(pm, nw, b, T2, T3)
end

""
function constraint_gp_power_branch_to(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    tr, ti = _PMs.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
end

""
function constraint_gp_power_branch_from(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    tr, ti = _PMs.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
end

""
function constraint_gp_power_branch_to_simplified(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    #tr, ti = _PMs.calc_branch_t(branch)
    #g_to = branch["g_to"]
    #b_to = branch["b_to"]
    #tm = branch["tap"]

    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_to_simplified(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, T2, T3)
end

""
function constraint_gp_power_branch_from_simplified(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)   
    branch = _PMs.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PMs.calc_branch_y(branch)
    #tr, ti = _PMs.calc_branch_t(branch)
    #g_fr = branch["g_fr"]
    #b_fr = branch["b_fr"]
    #tm = branch["tap"]
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_from_simplified(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, T2, T3)
end
""
function constraint_gp_gen_power(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    i   = _PMs.ref(pm, nw, :gen, g, "gen_bus")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_gen_power_real(pm, nw, i, g, T2, T3)
    constraint_gp_gen_power_imaginary(pm, nw, i, g, T2, T3)
end

""
function constraint_gp_load_power(pm::AbstractPowerModel, l::Int; nw::Int=nw_id_default)
    i   = _PMs.ref(pm, nw, :load, l, "load_bus") 

    pd  = _PMs.ref(pm, nw, :load, l, "pd")
    qd  = _PMs.ref(pm, nw, :load, l, "qd")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_load_power_real(pm, nw, i, l, pd, T2, T3)
    constraint_gp_load_power_imaginary(pm, nw, i, l, qd, T2, T3)
end

# chance constraint limit
""
function constraint_bus_voltage_cc_limit(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PMs.ref(pm, nw, :bus, i, "vmin")
    vmax = _PMs.ref(pm, nw, :bus, i, "vmax")
    
    λmin = _PMs.ref(pm, nw, :bus, i, "λvmin")
    λmax = _PMs.ref(pm, nw, :bus, i, "λvmax")
    
    T2 = pm.data["T2"]
    T4 = pm.data["T4"]

    constraint_bus_voltage_cc_limit(pm, i, vmin, vmax, λmin, λmax, T2, T4)
end

""
function constraint_bus_voltage_squared_cc_limit(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PMs.ref(pm, nw, :bus, i, "vmin")
    vmax = _PMs.ref(pm, nw, :bus, i, "vmax")
    
    λmin = _PMs.ref(pm, nw, :bus, i, "λvmin")
    λmax = _PMs.ref(pm, nw, :bus, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_bus_voltage_squared_cc_limit(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

""
function constraint_branch_series_current_cc_limit(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    cmax = _PMs.ref(pm, nw, :branch, b, "cmax")
    λmax = _PMs.ref(pm, nw, :branch, b, "λcmax")

    branch = _PMs.ref(pm, nw, :branch, b)

    gs, bs = _PMs.calc_branch_y(branch)
    
    T2 = pm.data["T2"]
    T4 = pm.data["T4"]

    constraint_branch_series_current_cc_limit(pm, b, cmax, λmax, T2, T4,gs,bs)
end

""
function constraint_branch_series_current_squared_cc_limit(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    cmax = _PMs.ref(pm, nw, :branch, b, "cmax")
    λmax = _PMs.ref(pm, nw, :branch, b, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_branch_series_current_squared_cc_limit(pm, b, cmax, λmax, T2, mop)
end


""
function constraint_gen_power_cc_limit(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    pmin = _PMs.ref(pm, nw, :gen, g, "pmin")
    pmax = _PMs.ref(pm, nw, :gen, g, "pmax")
    qmin = _PMs.ref(pm, nw, :gen, g, "qmin")
    qmax = _PMs.ref(pm, nw, :gen, g, "qmax")

    λpmin = _PMs.ref(pm, nw, :gen, g, "λpmin")
    λpmax = _PMs.ref(pm, nw, :gen, g, "λpmax")
    λqmin = _PMs.ref(pm, nw, :gen, g, "λqmin")
    λqmax = _PMs.ref(pm, nw, :gen, g, "λqmax")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_gen_power_real_cc_limit(pm, g, pmin, pmax, λpmin, λpmax, T2, mop)
    constraint_gen_power_imaginary_cc_limit(pm, g, qmin, qmax, λqmin, λqmax, T2, mop)
end

""
function constraint_gen_power_cc_limit_without_aux(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    pmin = _PMs.ref(pm, nw, :gen, g, "pmin")
    pmax = _PMs.ref(pm, nw, :gen, g, "pmax")
    qmin = _PMs.ref(pm, nw, :gen, g, "qmin")
    qmax = _PMs.ref(pm, nw, :gen, g, "qmax")

    λpmin = _PMs.ref(pm, nw, :gen, g, "λpmin")
    λpmax = _PMs.ref(pm, nw, :gen, g, "λpmax")
    λqmin = _PMs.ref(pm, nw, :gen, g, "λqmin")
    λqmax = _PMs.ref(pm, nw, :gen, g, "λqmax")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_gen_power_real_cc_limit_wo_aux(pm, g, pmin, pmax, λpmin, λpmax, T2, mop)
    constraint_gen_power_imaginary_cc_limit_wo_aux(pm, g, qmin, qmax, λqmin, λqmax, T2, mop)
end


""
function constraint_power_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = _PMs.ref(pm, nw, :bus, i)
    bus_arcs = _PMs.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PMs.ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = _PMs.ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = _PMs.ref(pm, nw, :bus_gens, i)
    bus_loads = _PMs.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PMs.ref(pm, nw, :bus_shunts, i)
    bus_storage = _PMs.ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => _PMs.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PMs.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PMs.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end