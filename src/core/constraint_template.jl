################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# current balance
""
function constraint_current_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(con(pm, nw), :kcl_cr)
        con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(con(pm, nw), :kcl_ci)
        con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i, bus_arcs, bus_arcs_dc, bus_gens, bus_loads, bus_gs, bus_bs)
end

# galerkin projection
""
function constraint_gp_squared_voltage(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = ref(pm, :T2)
    T3  = ref(pm, :T3)

    constraint_gp_squared_voltage(pm, nw, i, T2, T3)
end

""
function constraint_gp_branch_series_current_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = ref(pm, :T2)
    T3  = ref(pm, :T3)

    constraint_gp_branch_series_current_squared(pm, nw, b, T2, T3)
end

""
function constraint_gp_gen_power(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    i   = ref(pm, nw, :gen, g, "gen_bus")

    T2  = ref(pm, :T2)
    T3  = ref(pm, :T3)

    constraint_gp_gen_power_real(pm, nw, i, g, T2, T3)
    constraint_gp_gen_power_imaginary(pm, nw, i, g, T2, T3)
end

""
function constraint_gp_load_power(pm::AbstractPowerModel, l::Int; nw::Int=nw_id_default)
    i   = ref(pm, nw, :load, l, "load_bus") 

    pd  = ref(pm, nw, :load, l, "pd")
    qd  = ref(pm, nw, :load, l, "qd")

    T2  = ref(pm, :T2)
    T3  = ref(pm, :T3)

    constraint_gp_load_power_real(pm, nw, i, l, pd, T2, T3)
    constraint_gp_load_power_imaginary(pm, nw, i, l, qd, T2, T3)
end

# chance constraint limit
""
function constraint_bus_voltage_squared_cc_limit(pm::AbstractPowerModel, i::Int)
    vmin = ref(pm, nw, :bus, i, "vmin")
    vmax = ref(pm, nw, :bus, i, "vmax")
    
    λ    = ref(pm, nw, :bus, i, "λ")
    
    T2  = ref(pm, :T2)
    mop = ref(pm, mop)

    constraint_bus_voltage_squared_cc_limit(pm, i, vmin, vmax, λ, T2, mop)
end

""
function constraint_branch_series_current_squared_cc_limit(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    imax = ref(pm, nw, :branch, b, "imax")
    
    λ   = ref(pm, nw, :branch, b, "λ")
    
    T2  = ref(pm, :T2)
    mop = ref(pm, mop)

    constraint_branch_series_current_squared_cc_limit(pm, b, imax, λ, T2, mop)
end

""
function constraint_gen_power_cc_limit(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    pgmin = ref(pm, nw, :gen, g, "pgmin")
    pgmax = ref(pm, nw, :gen, g, "pgmax")
    qgmin = ref(pm, nw, :gen, g, "qgmin")
    qgmax = ref(pm, nw, :gen, g, "qgmax")

    λ   = ref(pm, nw, :gen, g, "λ")

    T2  = ref(pm, :T2)
    mop = ref(pm, mop)

    constraint_gen_power_real_cc_limit(pm, g, pgmin, pgmax, λ, T2, mop)
    constraint_gen_power_imaginary_cc_limit(pm, g, qgmin, qgmax, λ, T2, mop)
end