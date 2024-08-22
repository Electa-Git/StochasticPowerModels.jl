###################################################################################
#  Copyright 2024, Kaan Yurtseven                                                 #
###################################################################################
# StochasticPowerModels.jl                                                        #
# An extention package of PowerModels.jl for Stochastic Power System Optimization #
# See http://github.com/Electa-Git/StochasticPowerModels.jl                       #
###################################################################################

function add_ref_RES!(ref::Dict{Symbol,<:Any}, data::Dict{String,<:Any})
    for (n, nw_ref) in ref[:it][:pm][:nw]
        if haskey(nw_ref, :RES)
            bus_RES = Dict([(i, []) for (i,bus) in nw_ref[:bus]])
                for (i,RES) in nw_ref[:RES]
                    push!(bus_RES[RES["RES_bus"]], i)
                end
            nw_ref[:bus_RES] = bus_RES

        else
            nw_ref[:bus_RES] = Dict{String, Any}()
        end
    end
end

