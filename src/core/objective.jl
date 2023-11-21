
"expected cost of active power generation"
function objective_min_expected_generation_cost(pm::AbstractPowerModel; kwargs...)
    gen_cost = Dict()

    T2 = pm.data["T2"]

    for (g, gen) in _PM.ref(pm, :gen, nw=1)
        pg = Dict(nw => _PM.var(pm, nw, :pg, g) for nw in _PM.nw_ids(pm))

        if length(gen["cost"]) == 1
            gen_cost[g] = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            gen_cost[g] = gen["cost"][1]*pg[1] + 
                          gen["cost"][2]
        elseif length(gen["cost"]) == 3
            gen_cost[g] = gen["cost"][1]*sum(T2.get([n-1,n-1]) * pg[n]^2 for n in _PM.nw_ids(pm)) + 
                          gen["cost"][2]*pg[1] + 
                          gen["cost"][3]
        else
            gen_cost[g] = 0.0
        end
    end

    return JuMP.@objective(pm.model, Min,
            sum(gen_cost[g] for g in _PM.ids(pm, :gen, nw=1))
    )
end

# function objective_min_expected_generation_cost_dim(pm::AbstractPowerModel; kwargs...)
#     gen_cost = Dict()

#     T2 = _FP.dim_meta(pm, :PCE_coeff, "T2")

#     for (n, network) in _PM.nws(pm) 

#         if _FP.is_first_id(pm,n,:PCE_coeff)

#             for (g, gen) in _PM.ref(pm, :gen, nw=n)            
#                 pg = Dict(nw => _PM.var(pm, nw, :pg, g) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

#                 if length(gen["cost"]) == 1
#                     gen_cost[g,n] = gen["cost"][1]
#                 elseif length(gen["cost"]) == 2
#                     gen_cost[g,n] = gen["cost"][1]*pg[n] + 
#                                   gen["cost"][2]
#                 elseif length(gen["cost"]) == 3
#                     gen_cost[g,n] = gen["cost"][1]*sum(T2.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n1, :PCE_coeff)-1]) * pg[n1]^2 for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))) + 
#                                   gen["cost"][2]*pg[n] + 
#                                   gen["cost"][3]
#                 else
#                     gen_cost[g,n] = 0.0
#                 end

#             end

#         end

        
#     end

#     return JuMP.@objective(pm.model, Min,
#                             sum(gen_cost[g,n] for (n, network) in _PM.nws(pm) for g in _PM.ids(pm, :gen, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
#                             )

# end

function objective_min_expected_generation_cost_dim(pm::AbstractPowerModel; kwargs...)
    gen_cost = Dict()
    pd_curt_cost = Dict()
    p_RES_curt_cost = Dict()

    T2 = _FP.dim_meta(pm, :PCE_coeff, "T2")
    curt_status = pm.data["curtailment"]

    VOLL = 1e8;

    RES_curt_penalty = 1e8;



    for (n, network) in _PM.nws(pm) 

        if _FP.is_first_id(pm,n,:PCE_coeff)

            for (g, gen) in _PM.ref(pm, :gen, nw=n)            
                pg = Dict(nw => _PM.var(pm, nw, :pg, g) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))

                if length(gen["cost"]) == 1
                    gen_cost[g,n] = gen["cost"][1]
                elseif length(gen["cost"]) == 2
                    gen_cost[g,n] = gen["cost"][1]*pg[n] + 
                                  gen["cost"][2]
                elseif length(gen["cost"]) == 3
                    gen_cost[g,n] = gen["cost"][1]*sum(T2.get([_FP.coord(pm, n1, :PCE_coeff)-1, _FP.coord(pm, n1, :PCE_coeff)-1]) * pg[n1]^2 for n1 in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff))) + 
                                  gen["cost"][2]*pg[n] + 
                                  gen["cost"][3]
                else
                    gen_cost[g,n] = 0.0
                end

            end
            
            if curt_status["Load Curtailment"] == true
                for (l, load) in _PM.ref(pm, :load, nw=n)

                    pd_curt = Dict(nw => _PM.var(pm, nw, :pd_curt, l) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    
                    pd_curt_cost[l,n] = pd_curt[n] * VOLL

                end
            end

            if curt_status["RES Curtailment"] == true
                for (p, RES) in _PM.ref(pm, :RES, nw=n)

                    p_RES_curt = Dict(nw => _PM.var(pm, nw, :p_RES_curt, p) for nw in _FP.similar_ids(pm, n; PCE_coeff=1:_FP.dim_length(pm, :PCE_coeff)))
                    
                    p_RES_curt_cost[p,n] = p_RES_curt[n] * RES_curt_penalty

                end
            end

        end

        
    end



    if (curt_status["Load Curtailment"] == false) && (curt_status["RES Curtailment"] == false) 
        return JuMP.@objective(pm.model, Min,
                            sum(gen_cost[g,n] for (n, network) in _PM.nws(pm) for g in _PM.ids(pm, :gen, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
                            )
    elseif (curt_status["Load Curtailment"] == true) && (curt_status["RES Curtailment"] == false) 
        return JuMP.@objective(pm.model, Min,
                            sum(gen_cost[g,n] for (n, network) in _PM.nws(pm) for g in _PM.ids(pm, :gen, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
                            + sum(pd_curt_cost[l,n] for (n, network) in _PM.nws(pm) for l in _PM.ids(pm, :load, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
                            )
    elseif (curt_status["Load Curtailment"] == false) && (curt_status["RES Curtailment"] == true) 
        return JuMP.@objective(pm.model, Min,
                            sum(gen_cost[g,n] for (n, network) in _PM.nws(pm) for g in _PM.ids(pm, :gen, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
                            + sum(p_RES_curt_cost[p,n] for (n, network) in _PM.nws(pm) for p in _PM.ids(pm, :RES, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
                            )
    elseif (curt_status["Load Curtailment"] == true) && (curt_status["RES Curtailment"] == true) 
        return JuMP.@objective(pm.model, Min,
                            sum(gen_cost[g,n] for (n, network) in _PM.nws(pm) for g in _PM.ids(pm, :gen, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
                            + sum(pd_curt_cost[l,n] for (n, network) in _PM.nws(pm) for l in _PM.ids(pm, :load, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
                            + sum(p_RES_curt_cost[p,n] for (n, network) in _PM.nws(pm) for p in _PM.ids(pm, :RES, nw=n) if _FP.is_first_id(pm,n,:PCE_coeff))
                            )
    end

end
