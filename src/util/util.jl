################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

""
var_min(var, λ) = var[1] - λ * sqrt(sum(var[2:end].^2))
var_max(var, λ) = var[1] + λ * sqrt(sum(var[2:end].^2))

""
is_constrained(pm, cmp, idx) = _PMs.ref(pm, 1, cmp, idx, "cstr")

""
function parse_dst(dst, pa, pb, deg)
    dst == "Beta"    && return _PCE.Beta01OrthoPoly(deg, pa, pb; Nrec=5*deg)
    dst == "Normal"  && return _PCE.GaussOrthoPoly(deg; Nrec=5*deg)
    dst == "Uniform" && return _PCE.Uniform01OrthoPoly(deg; Nrec=5*deg)
end

"""
    build_stochastic_data!(data::Dict{String,Any}, deg::Int)
"""
function build_stochastic_data(data, deg)
    # add maximum current
    for (nb, branch) in data["branch"]
        f_bus = branch["f_bus"]
        branch["cmax"] = branch["rate_a"] / data["bus"]["$f_bus"]["vmin"]
    end

    # build mop
    opq = [parse_dst(ns[2]["dst"], ns[2]["pa"], ns[2]["pb"], deg) for ns in data["sdata"]]
    mop = _PCE.MultiOrthoPoly(opq, deg)

    # build load matrix
    Nd, Npce = length(data["load"]), mop.dim
    pd, qd = zeros(Nd, Npce), zeros(Nd, Npce)
    for nd in 1:Nd 
        # reactive power
        qd[nd,1] = data["load"]["$nd"]["qd"]
        # active power
        nb = data["load"]["$nd"]["load_bus"]
        ni = data["bus"]["$nb"]["dst_id"]
        if ni == 0
            pd[nd,1] = data["load"]["$nd"]["pd"]
        else
            base = data["baseMVA"]
            μ, σ = data["bus"]["$nb"]["μ"] / base, data["bus"]["$nb"]["σ"] / base
            if mop.uni[ni] isa _PCE.GaussOrthoPoly
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            else
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni], kind="μσ")
            end
        end
    end

    # replicate the data
    data = _PMs.replicate(data, Npce)

    # add the stochastic data 
    data["T2"] = _PCE.Tensor(2,mop)
    data["T3"] = _PCE.Tensor(3,mop)
    data["mop"] = mop
    for nw in 1:Npce, nd in 1:Nd
        data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
        data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
    end

    return data
end

# ""
# function create_mop(data_dist)
#     m=Vector{Symbol}()
#     for id in sort!(collect(keys(data_dist)))
#         push!(m, Symbol(data_dist[id]["distribution"]))
#     end

#     m2=[]
#     for i=1:length(m)
#         dist_sym=m[i]
#         d=data_dist["$i"]
#         degree = d["deg"]
#         No_rec = (d["deg"] *d["Nrec"])
#         alpha= d["alpha"]
#         beta=d["beta"]
#         if dist_sym == :Beta01OrthoPoly
#             a = Meta.parse(string("$dist_sym($degree, $alpha, $beta; Nrec=$No_rec)"))
#         else
#             a = Meta.parse(string("$dist_sym($degree; Nrec=$No_rec)"))
#         end
#         push!(m2, a)
#     end
#     return m2
# end