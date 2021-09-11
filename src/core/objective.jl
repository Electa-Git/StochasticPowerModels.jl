################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

"expected cost of active power generation"
function objective_min_expected_generation_cost(pm::AbstractPowerModel; kwargs...)
    gen_cost = Dict()

    T2 = pm.data["T2"]

    for (g, gen) in _PMs.ref(pm, :gen, nw=1)
        pg = Dict(nw => _PMs.var(pm, nw, :pg, g) for nw in _PMs.nw_ids(pm))

        if length(gen["cost"]) == 1
            gen_cost[g] = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            gen_cost[g] = gen["cost"][1]*pg[1] + 
                          gen["cost"][2]
        elseif length(gen["cost"]) == 3
            gen_cost[g] = gen["cost"][1]*sum(T2.get([n-1,n-1]) * pg[n]^2 for n in _PMs.nw_ids(pm)) + 
                          gen["cost"][2]*pg[1] + 
                          gen["cost"][3]
        else
            gen_cost[g] = 0.0
        end
    end

    return JuMP.@objective(pm.model, Min,
            sum(gen_cost[g] for g in _PMs.ids(pm, :gen, nw=1))
    )
end

""
function objective_min_expected_fuel_cost(pm::AbstractPowerModel; kwargs...)
    model = _PMs.check_gen_cost_models(pm)

    if model == 1 
        #return objective_min_expected_fuel_cost_pwl(pm; kwargs...)
        Memento.error(_LOGGER, "pwl cost model not supported atm.")
    elseif model == 2
        return objective_min_expected_fuel_cost_polynomial(pm; kwargs...)
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end

""
function objective_min_fuel_cost_poly(pm::AbstractPowerModel; kwargs...)
    model = _PMs.check_gen_cost_models(pm)

    if model == 1 
        #return objective_min_expected_fuel_cost_pwl(pm; kwargs...)
        Memento.error(_LOGGER, "pwl cost model not supported atm.")
    elseif model == 2
        order = _PMs.calc_max_cost_index(pm.data)-1
        if order <=2
            return _objective_min_fuel_cost_poly_based(pm; kwargs...)
        end
    else
        Memento.error(_LOGGER, "Only cost models of types 1 and 2 are supported at this time, given cost model type of $(model)")
    end

end


""
function objective_min_expected_fuel_cost_polynomial(pm::AbstractPowerModel; kwargs...)
    mop = pm.data["mop"]
    order = _PMs.calc_max_cost_index(pm.data)-1

    if order <= 2
        return _objective_min_expected_fuel_cost_polynomial_linquad(pm, mop; kwargs...)
    else
        return _objective_min_expected_fuel_cost_polynomial_nl(pm, mop; kwargs...)
    end
end

""
function _objective_min_expected_fuel_cost_polynomial_linquad(pm::AbstractPowerModel, mop; report::Bool=true)
    gen_cost = Dict()
    
    for (g, gen) in _PMs.ref(pm, :gen, nw=1)
        exp_pg = _PCE.mean([_PMs.var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)], mop)        
        if length(gen["cost"]) == 1
            gen_cost[g] = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            gen_cost[g] = gen["cost"][1]*exp_pg + gen["cost"][2]
        elseif length(gen["cost"]) == 3
            gen_cost[g] = gen["cost"][1]*exp_pg^2 + gen["cost"][2]*exp_pg + gen["cost"][3]
        else
            gen_cost[g] = 0.0
        end
    end

    return JuMP.@objective(pm.model, Min,
            sum(gen_cost[g] for g in _PMs.ids(pm, :gen, nw=1))
    )
end

""
function _objective_min_fuel_cost_polynomial_nl(pm::AbstractPowerModel, mop; report::Bool=true)
    gen_cost = Dict()
    
    for (g, gen) in _PMs.ref(pm, :gen,nw=1)
        exp_pg = _PCE.mean([var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)], mop)    

        cost_rev = reverse(gen["cost"])
        if length(cost_rev) == 1
            gen_cost[g] = JuMP.@NLexpression(pm.model, cost_rev[1])
        elseif length(cost_rev) == 2
            gen_cost[g] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*exp_pg)
        elseif length(cost_rev) == 3
            gen_cost[g] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*exp_pg + cost_rev[3]*exp_pg^2)
        elseif length(cost_rev) >= 4
            cost_rev_nl = cost_rev[4:end]
            gen_cost[g] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*exp_pg + cost_rev[3]*exp_pg^2 + sum( v*exp_pg^(d+3) for (d,v) in enumerate(cost_rev_nl)) )
        else
            gen_cost[g] = JuMP.@NLexpression(pm.model, 0.0)
        end
    end

    return JuMP.@NLobjective(pm.model, Min,
        sum(gen_cost[g] for g in _PMs.ids(pm, :gen))
    )
end 

"only for 2nd degree/ square termed to be checked"
function _objective_min_fuel_cost_poly_based(pm::AbstractPowerModel; report::Bool=true)
    gen_cost = Dict()
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    
    for (g, gen) in _PMs.ref(pm, :gen,nw=1)
        #exp_pg = _PCE.mean([var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)], mop)    
        pg_sum= sum([_PMs.var(pm, n, :pg, g) * T2.get([n-1,n-1]) for n in sorted_nw_ids(pm)])
        pg_sum_square= sum([_PMs.var(pm, n, :pg, g)^2 * T2.get([n-1,n-1]) for n in sorted_nw_ids(pm)])
        cost_rev = reverse(gen["cost"])
        if length(cost_rev) == 1
            gen_cost[g] = JuMP.@expression(pm.model, cost_rev[1])
        elseif length(cost_rev) == 2
            gen_cost[g] = JuMP.@expression(pm.model, cost_rev[1] + cost_rev[2]*pg_sum)
        elseif length(cost_rev) == 3
            gen_cost[g] = JuMP.@expression(pm.model, cost_rev[1] + cost_rev[2]*pg_sum + cost_rev[3]*pg_sum_square)
        #elseif length(cost_rev) >= 4
        #    cost_rev_nl = cost_rev[4:end]
        #    gen_cost[g] = JuMP.@NLexpression(pm.model, cost_rev[1] + cost_rev[2]*pg_sum + cost_rev[3]*pg^2 + sum( v*exp_pg^(d+3) for (d,v) in enumerate(cost_rev_nl)) )
        else
            gen_cost[g] = JuMP.@expression(pm.model, 0.0)
        end
    end

    return JuMP.@objective(pm.model, Min,
        sum(gen_cost[g] for g in _PMs.ids(pm, 1, :gen))
    )
end