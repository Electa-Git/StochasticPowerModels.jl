"""
    StochasticPowerModels.extend_matlab_file(path::String)
"""
function extend_matlab_file(path::String)
    # data
    data = _PM.parse_file(path)

    # general data
    baseMVA = data["baseMVA"]

    # bus data 
    Nb   = length(data["bus"])
    μ, σ = zeros(Nb), zeros(Nb)
    for (l,load) in data["load"]
        bus = load["load_bus"]
        μ[bus] = load["pd"] * baseMVA
        σ[bus] = load["pd"] * baseMVA * 0.10
    end
    for (b,bus) in data["bus"]
        bus["dst_id"]   = 0
        bus["μ"]        = μ[parse(Int,b)]
        bus["σ"]        = σ[parse(Int,b)]
        bus["λvmin"]    = 1.03643
        bus["λvmax"]    = 1.03643
    end

    # generator data
    for (g,gen) in data["gen"]
        gen["λpmin"] = 1.03643
        gen["λpmax"] = 1.03643
        gen["λqmin"] = 1.03643
        gen["λqmax"] = 1.03643
    end

    # branch data
    for (l,branch) in data["branch"]
        branch["λcmax"] = 1.03643
    end

    # export file
    _PM.export_file(path[1:end-2] * "_spm.m", data)
end

function extend_matlab_file_AC(path::String)
    # data
    data = _PM.parse_file(path)

    α = 0.95; 

    λ_val = Distributions.quantile.(Distributions.Normal(), [1-α, α]);

    λ_val = round(λ_val[2], digits = 3);

    # λ_val = 1.65; #0.95 = 1.65, 0.90 = 1.285, 0.85 = 1.03643, 0.80 = 0.83

    # general data
    baseMVA = data["baseMVA"]

    # bus data 
    Nb   = length(data["bus"])
    μ, σ = zeros(Nb), zeros(Nb)
    for (l,load) in data["load"]

        bus = load["load_bus"]
        μ[bus] = load["pd"] * baseMVA
        # μ[bus] = 0.5 * load["pd"] * baseMVA
        σ[bus] = abs(load["pd"] * baseMVA * 0.10)
        #σ[bus] = 0
    end

    # for (l,load) in data["load"]

    #     bus = load["load_bus"]
    #     μ[bus] = load["pd"] * baseMVA

    #     if bus == 3
    #         σ[bus] = abs(load["pd"] * baseMVA * 0.10)

    #     elseif bus == 4
    #         σ[bus] = abs(load["pd"] * baseMVA * 0.15)
        
    #     else
    #         σ[bus] = 0

    #     end
    # end
    
    # for (b,bus) in data["bus"]

    #     if parse(Int,b) == 3
            
    #         bus["dst_id"]   = 1
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
        
    #     elseif parse(Int,b) == 4
            
    #         bus["dst_id"]   = 2
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
        
    #     else
    #         bus["dst_id"]   = 0
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
    #     end

    # end

    # for (b,bus) in data["bus"]

    #     if parse(Int,b) == 6
            
    #         bus["dst_id"]   = 1
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
        
    #     elseif parse(Int,b) == 20
            
    #         bus["dst_id"]   = 2
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
        
    #     else
    #         bus["dst_id"]   = 0
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
    #     end

    # end

    for (b,bus) in data["bus"]

        if parse(Int,b) in 1:67
            
            bus["dst_id"]   = 1
            bus["μ"]        = μ[parse(Int,b)]
            bus["σ"]        = σ[parse(Int,b)]
            bus["λvmin"]    = λ_val
            bus["λvmax"]    = λ_val
        
        else
            bus["dst_id"]   = 0
            bus["μ"]        = μ[parse(Int,b)]
            bus["σ"]        = σ[parse(Int,b)]
            bus["λvmin"]    = λ_val
            bus["λvmax"]    = λ_val
        end

    end
    
    
    # for (b,bus) in data["bus"]

    #     bus["dst_id"]   = 1
    #     bus["μ"]        = μ[parse(Int,b)]
    #     bus["σ"]        = σ[parse(Int,b)]
    #     bus["λvmin"]    = λ_val
    #     bus["λvmax"]    = λ_val


    # end
    
    # generator data
    for (g,gen) in data["gen"]
        gen["λpmin"] = λ_val
        gen["λpmax"] = λ_val
        gen["λqmin"] = λ_val
        gen["λqmax"] = λ_val
    end

    # branch data
    for (l,branch) in data["branch"]
        branch["λcmax"] = λ_val
    end

    # export file
    _PM.export_file(path[1:end-2] * "_SPM_$(round(Int,(100*α)))" * "cc_alpha.m", data)
end

function extend_matlab_file_ACDC(path::String)
    # data
    data = _PM.parse_file(path)

    α = 0.95; 

    λ_val = Distributions.quantile.(Distributions.Normal(), [1-α, α])

    λ_val = round(λ_val[2], digits = 3);

    # λ_val = 1.65; #0.95 = 1.65, 0.90 = 1.285, 0.85 = 1.03643, 0.80 = 0.83

    # general data
    baseMVA = data["baseMVA"]

    # bus data 
    Nb   = length(data["bus"])
    μ, σ = zeros(Nb), zeros(Nb)
    for (l,load) in data["load"]

        bus = load["load_bus"]
        μ[bus] = load["pd"] * baseMVA
        # μ[bus] = 0.5 * load["pd"] * baseMVA
        σ[bus] = abs(load["pd"] * baseMVA * 0.10)
        #σ[bus] = 0
    end

    # for (l,load) in data["load"]

    #     bus = load["load_bus"]
    #     μ[bus] = load["pd"] * baseMVA

    #     if bus == 3
    #         σ[bus] = abs(load["pd"] * baseMVA * 0.10)

    #     elseif bus == 4
    #         σ[bus] = abs(load["pd"] * baseMVA * 0.15)
        
    #     else
    #         σ[bus] = 0

    #     end
    # end
    
    for (b,bus) in data["bus"]

        if parse(Int,b) == 3
            
            bus["dst_id"]   = 1
            bus["μ"]        = μ[parse(Int,b)]
            bus["σ"]        = σ[parse(Int,b)]
            bus["λvmin"]    = λ_val
            bus["λvmax"]    = λ_val
        
        elseif parse(Int,b) == 4
            
            bus["dst_id"]   = 2
            bus["μ"]        = μ[parse(Int,b)]
            bus["σ"]        = σ[parse(Int,b)]
            bus["λvmin"]    = λ_val
            bus["λvmax"]    = λ_val
        
        else
            bus["dst_id"]   = 0
            bus["μ"]        = μ[parse(Int,b)]
            bus["σ"]        = σ[parse(Int,b)]
            bus["λvmin"]    = λ_val
            bus["λvmax"]    = λ_val
        end

    end

    # for (b,bus) in data["bus"]

    #     if parse(Int,b) == 6
            
    #         bus["dst_id"]   = 1
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
        
    #     elseif parse(Int,b) == 20
            
    #         bus["dst_id"]   = 2
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
        
    #     else
    #         bus["dst_id"]   = 0
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
    #     end

    # end

    # for (b,bus) in data["bus"]

    #     if parse(Int,b) in 1:67
            
    #         bus["dst_id"]   = 1
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
        
    #     else
    #         bus["dst_id"]   = 0
    #         bus["μ"]        = μ[parse(Int,b)]
    #         bus["σ"]        = σ[parse(Int,b)]
    #         bus["λvmin"]    = λ_val
    #         bus["λvmax"]    = λ_val
    #     end

    # end
        
    # for (b,bus) in data["bus"]

    #     bus["dst_id"]   = 1
    #     bus["μ"]        = μ[parse(Int,b)]
    #     bus["σ"]        = σ[parse(Int,b)]
    #     bus["λvmin"]    = λ_val
    #     bus["λvmax"]    = λ_val


    # end

    # dcbus data 
    for (b,busdc) in data["busdc"]
        busdc["λvmin"]    = λ_val
        busdc["λvmax"]    = λ_val
    end

    # convdc data 
    for (b,convdc) in data["convdc"]
        convdc["λvmin"]    = λ_val
        convdc["λvmax"]    = λ_val
    end


    # generator data
    for (g,gen) in data["gen"]
        gen["λpmin"] = λ_val
        gen["λpmax"] = λ_val
        gen["λqmin"] = λ_val
        gen["λqmax"] = λ_val
    end

    # branch data
    for (l,branch) in data["branch"]
        branch["λcmax"] = λ_val
    end

    # dc branch data
    for (l,branchdc) in data["branchdc"]
        branchdc["λcmax"] = λ_val
    end

    # export file
    _PM.export_file(path[1:end-2] * "_SPM_$(round(Int,(100*α)))" * "cc_alpha.m", data)
end
