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

function constraint_current_balance_with_RES(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(_PM.con(pm, nw), :kcl_cr)
         _PM.con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
     end
     if !haskey(_PM.con(pm, nw), :kcl_ci)
         _PM.con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
     end
 
     bus = _PM.ref(pm, nw, :bus, i)
     bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
     bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
     bus_gens = _PM.ref(pm, nw, :bus_gens, i)
     bus_convs_ac = _PM.ref(pm, nw, :bus_convs_ac, i)
     bus_loads = _PM.ref(pm, nw, :bus_loads, i)
     bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)
 
     bus_RES = _PM.ref(pm, nw, :bus_RES, i) 
 
 
     bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
     bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)
 
     constraint_current_balance_with_RES(pm, nw, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs, bus_RES)
 
 end

 function constraint_current_balance_dc(pm::_PM.AbstractIVRModel, i::Int; nw::Int=_PM.nw_id_default)
    bus_arcs_dcgrid = _PM.ref(pm, nw, :bus_arcs_dcgrid, i)
    bus_convs_dc = _PM.ref(pm, nw, :bus_convs_dc, i)
    pd = _PM.ref(pm, nw, :busdc, i)["Pdc"]
    constraint_current_balance_dc(pm, nw, bus_arcs_dcgrid, bus_convs_dc, pd)
end

function constraint_gp_ohms_dc_branch(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    branch = _PM.ref(pm, nw, :branchdc, b)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (b, f_bus, t_bus)
    t_idx = (b, t_bus, f_bus)

    p = _PM.ref(pm, nw, :dcpol)

    constraint_gp_ohms_dc_branch(pm, nw, b, T2, T3, f_bus, t_bus, f_idx, t_idx, branch["r"], p)
end

function constraint_ohms_dc_branch(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)

    branch = _PM.ref(pm, nw, :branchdc, b)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (b, f_bus, t_bus)
    t_idx = (b, t_bus, f_bus)

    p = _PM.ref(pm, nw, :dcpol)

    constraint_ohms_dc_branch(pm, nw, b, f_bus, t_bus, f_idx, t_idx, branch["r"], p)
end

function constraint_gp_filter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_filter_voltage_squared(pm, nw, i, T2, T3)
end

function constraint_gp_converter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_voltage_squared(pm, nw, i, T2, T3)
end

function constraint_gp_transformer_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_transformer_current_from_squared(pm, nw, i, T2, T3)
end

function constraint_gp_transformer_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_transformer_current_to_squared(pm, nw, i, T2, T3)
end

function constraint_gp_reactor_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_reactor_current_from_squared(pm, nw, i, T2, T3)
end

function constraint_gp_reactor_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_reactor_current_to_squared(pm, nw, i, T2, T3)
end

function constraint_gp_converter_current_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_current_squared(pm, nw, i, T2, T3)
end

function constraint_gp_iconv_lin_squared_1(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_iconv_lin_squared_1(pm, nw, i, T2, T3)
end

function constraint_gp_iconv_lin_squared_2(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_iconv_lin_squared_2(pm, nw, i, T2, T3)
end

function constraint_gp_converter_dc_power(pm::_PM.AbstractIVRModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    b_idx =   conv["busdc_i"]

    constraint_gp_converter_dc_power(pm, nw, i, T2, T3, b_idx)
end

function constraint_gp_converter_losses(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    
    conv = _PM.ref(pm, nw, :convdc, i)
    a = conv["LossA"]
    b = conv["LossB"]
    c = conv["LossCinv"]
    plmax = conv["LossA"] + conv["LossB"] * conv["Pacrated"] + conv["LossCinv"] * (conv["Pacrated"])^2
    constraint_gp_converter_losses(pm, nw, i, T2, T3, a, b, c, plmax)
end

function constraint_gp_converter_ac_power(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_converter_ac_power(pm, nw, i, T2, T3)
end

function constraint_conv_reactor(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    constraint_conv_reactor(pm, nw, i, conv["rc"], conv["xc"], Bool(conv["reactor"]))
end

#
function constraint_conv_filter(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    constraint_conv_filter(pm, nw, i, conv["bf"], Bool(conv["filter"]) )
end

#
function constraint_conv_transformer(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    constraint_conv_transformer(pm, nw, i, conv["rtf"], conv["xtf"], conv["busac_i"], conv["tm"], Bool(conv["transformer"]))
end

function constraint_gp_RES_power(pm::AbstractPowerModel, p::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :RES, p, "RES_bus") 

    pd  = _PM.ref(pm, nw, :RES, p, "pd")
    qd  = _PM.ref(pm, nw, :RES, p, "qd")

    # p_size= _PM.var(pm, 1, :p_size, p)
    # q_size= _PM.var(pm, 1, :q_size, p)

    p_size = _PM.ref(pm, nw, :RES, p, "p_size")
    q_size = _PM.ref(pm, nw, :RES, p, "q_size")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    # c = pm.data["curt"]

    constraint_gp_RES_power_real(pm, nw, i, p, pd, T2, T3, p_size)
    constraint_gp_RES_power_imaginary(pm, nw, i, p, qd, T2, T3, q_size)
end

function constraint_cc_filter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :convdc, i, "Vmmin")
    vmax = _PM.ref(pm, nw, :convdc, i, "Vmmax")
    
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_filter_voltage_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_converter_voltage_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :convdc, i, "Vmmin")
    vmax = _PM.ref(pm, nw, :convdc, i, "Vmmax")
    
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_converter_voltage_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_transformer_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_transformer_current_from_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_transformer_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_transformer_current_to_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_reactor_current_from_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_reactor_current_from_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_reactor_current_to_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_reactor_current_to_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_converter_current_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_converter_current_squared(pm, i, Imax, λmax, T2, mop)
end


function constraint_cc_dc_branch_current(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    vpu = 1;
    branch = _PM.ref(pm, nw, :branchdc, i)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    
    Imax = branch["rateA"]/vpu
    Imin = - branch["rateA"]/vpu
   
       
    λmax = _PM.ref(pm, nw, :branchdc, i, "λcmax")
    λmin = _PM.ref(pm, nw, :branchdc, i, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_dc_branch_current(pm, i, Imax, Imin, λmax, λmin, f_idx, t_idx, T2, mop)
end

function constraint_cc_iconv_lin_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Imax"]
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_iconv_lin_squared(pm, i, Imax, λmax, T2, mop)
end

function constraint_cc_iconv_lin(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Imax"]
    Imin = 0
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmax") # All λ values are equal.
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_iconv_lin(pm, i, Imax, Imin, λmax, λmin, T2, mop)
end

function constraint_cc_conv_ac_power(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Pacmin = _PM.ref(pm, nw, :convdc, i, "Pacmin")
    Pacmax = _PM.ref(pm, nw, :convdc, i, "Pacmax")
    Qacmin = _PM.ref(pm, nw, :convdc, i, "Qacmin")
    Qacmax = _PM.ref(pm, nw, :convdc, i, "Qacmax")

    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")

    #λqmin = _PM.ref(pm, nw, :convdc, g, "λqmin")
    #λqmax = _PM.ref(pm, nw, :convdc, g, "λqmax")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_cc_conv_ac_power_real(pm, i, Pacmin, Pacmax, λmin, λmax, T2, mop)
    constraint_cc_conv_ac_power_imaginary(pm, i, Qacmin, Qacmax, λmin, λmax, T2, mop)
end

function constraint_cc_conv_dc_power(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    Pdcmin = - _PM.ref(pm, nw, :convdc, i, "Pacrated")
    Pdcmax = _PM.ref(pm, nw, :convdc, i, "Pacrated")
    
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")


    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_cc_conv_dc_power(pm, i, Pdcmin, Pdcmax, λmin, λmax, T2, mop)
end


function constraint_cc_converter_dc_current(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vpu = 1;
    bigM = 1.2;
    conv = _PM.ref(pm, nw, :convdc, i)
    
    Imax = conv["Pacrated"]/vpu * bigM
    Imin = - conv["Pacrated"]/vpu * bigM
    
    λmax = _PM.ref(pm, nw, :convdc, i, "λvmax")
    λmin = _PM.ref(pm, nw, :convdc, i, "λvmin")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_converter_dc_current(pm, i, Imax, Imin, λmax, λmin, T2, mop)
end


function constraint_cc_conv_voltage_magnitude(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
 
    vmin = _PM.ref(pm, nw, :busdc, i, "Vdcmin")
    vmax = _PM.ref(pm, nw, :busdc, i, "Vdcmax")
    
    λmin = _PM.ref(pm, nw, :busdc, i, "λvmin")
    λmax = _PM.ref(pm, nw, :busdc, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_conv_voltage_magnitude(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

function constraint_cc_dc_branch_current_on_off(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    vpu = 1;
    branch = _PM.ref(pm, nw, :branchdc, i)
    f_bus = branch["fbusdc"]
    t_bus = branch["tbusdc"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    
    Imax = branch["rateA"]/vpu
    Imin = - branch["rateA"]/vpu
   
       
    λmax = _PM.ref(pm, nw, :branchdc, i, "λcmax")
    λmin = _PM.ref(pm, nw, :branchdc, i, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_dc_branch_current_on_off(pm, i, Imax, Imin, λmax, λmin, f_idx, t_idx, T2, mop)
end