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

