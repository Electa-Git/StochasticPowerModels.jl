################################################################################
#  Copyright 2021, Hakan Ergun                                              #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
function add_stochastic_data!(data)
    for (b, bus) in data["bus"]
        demand = 0
        for (l, load) = data["load"]
            if load["load_bus"] == parse(Int, b)
                demand =  load["pd"]
            end
        end
        bus["μ"] = demand
        bus["σ"] = bus["μ"] / 10
        bus["λvmax"] = 1.6
        bus["λvmin"] = 1.6
        bus["dst_id"] = mod(parse(Int, b), 3)
    end

    for (g, gen) in data["gen"]
        gen["λpmin"] = 1.6
        gen["λpmax"] = 1.6
        gen["λqmin"] = 1.6
        gen["λqmax"] = 1.6
    end

    for (br, branch) in data["branch"]
        branch["λcmax"] = 1.6
    end

    data["sdata"] = Dict{String, Any}()
    data["sdata"]["1"] = Dict{String, Any}()
    data["sdata"]["1"]["pa"] = 0.0
    data["sdata"]["1"]["pb"] = 0.0
    data["sdata"]["1"]["index"] = 1
    data["sdata"]["1"]["dst"] = "Uniform"
    data["sdata"]["1"]["source_id"] = ["sdata", 1]

    data["sdata"]["2"] = Dict{String, Any}()
    data["sdata"]["2"]["pa"] = 2.0
    data["sdata"]["2"]["pb"] = 2.0
    data["sdata"]["2"]["index"] = 2
    data["sdata"]["2"]["dst"] = "Beta"
    data["sdata"]["2"]["source_id"] = ["sdata", 2]

    data["sdata"]["3"] = Dict{String, Any}()
    data["sdata"]["3"]["pa"] = 0.0
    data["sdata"]["3"]["pb"] = 0.0
    data["sdata"]["3"]["index"] = 3
    data["sdata"]["3"]["dst"] = "Normal"
    data["sdata"]["3"]["source_id"] = ["sdata", 3]
    return data
end