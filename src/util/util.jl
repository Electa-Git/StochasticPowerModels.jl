################################################################################
#  Copyright 2021, Tom Van Acker, Arpan Koirala                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# input data
""
function parse_dst(dst, pa, pb, deg)
    dst == "Beta"    && return _PCE.Beta01OrthoPoly(deg, pa, pb; Nrec=5*deg)
    dst == "Normal"  && return _PCE.GaussOrthoPoly(deg; Nrec=5*deg)
    dst == "Uniform" && return _PCE.Uniform01OrthoPoly(deg; Nrec=5*deg)
end

"""
    StochasticPowerModels.build_stochastic_data(data::Dict{String,Any}, deg::Int)

Function to build the multi-network data representative of the polynomial chaos
expansion of a single-network data dictionary.
"""
function build_stochastic_data(data::Dict{String,Any}, deg::Int)
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
    data = _PM.replicate(data, Npce)

    # add the stochastic data 
    data["T2"] = _PCE.Tensor(2,mop)
    data["T3"] = _PCE.Tensor(3,mop)
    data["T4"] = _PCE.Tensor(4,mop)
    data["mop"] = mop
    for nw in 1:Npce, nd in 1:Nd
        data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
        data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
    end

    return data
end


function build_stochastic_data_ACDC_RES(data::Dict{String,Any}, deg::Int, p_size)

    # add maximum current
    for (nb, branch) in data["branch"]
        f_bus = branch["f_bus"]
        branch["cmax"] = branch["rate_a"] #/ data["bus"]["$f_bus"]["vmin"]
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


    # build RES matrix
    Nd_g = length(data["RES"])
    pd_g, qd_g = zeros(Nd_g, Npce), zeros(Nd_g, Npce)

    # Take the last entry of sdata as the irradiance.
    #data["sdata"]= s_dict 

    # data["RES"]=deepcopy(data["load"]); #this is temporary. Fix it to add PV in matpower format
    # [data["RES"][d]["μ"] = 0.8336 * 27.9230 for d in   keys(data["RES"])]
    # # [data["RES"][d]["σ"] = 0.09 for d in   keys(data["RES"])]
    # [data["RES"][d]["σ"] = 0.8336 * 27.9230 * 0.1079 for d in   keys(data["RES"])]
    
    [data["RES"][d]["μ"]=0 for d in   keys(data["RES"])]
    [data["RES"][d]["σ"]=1 for d in   keys(data["RES"])] # μ, σ parameters are support of Beta (a, b) not (α, β)
    [data["RES"][d]["p_size"] = p_size for d in keys(data["RES"])]
    [data["RES"][d]["q_size"] = p_size for d in keys(data["RES"])]   
    [data["RES"][d]["qd"] = 0 for d in keys(data["RES"])]
    


    for nd_g in 1:Nd_g 
        # reactive power
        qd_g[nd_g,1] = data["RES"]["$nd_g"]["qd"]
        # active power

        # np = length(opq) #Last distribution is for irradiance
        np_g = data["RES"]["$nd_g"]["dst_id"]

        base = data["baseMVA"]
        # μ, σ = data["RES"]["1"]["μ"] / base, data["RES"]["1"]["σ"] / base
        μ, σ = data["RES"]["1"]["μ"], data["RES"]["1"]["σ"]
        
            if mop.uni[np_g] isa _PCE.GaussOrthoPoly
                pd_g[nd_g,[1,np_g+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np_g])
                qd_g[nd_g,[1,np_g+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np_g])
            else
                pd_g[nd_g,[1,np_g+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np_g]) #Beta parameters without 'kind="μσ"'
                qd_g[nd_g,[1,np_g+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np_g]) # μ, σ parameters are support of Beta (a, b) not (α, β)
            end
       
    end

        # replicate the data
        data = _PM.replicate(data, Npce)

        # add the stochastic data 
        data["T2"] = _PCE.Tensor(2,mop)
        data["T3"] = _PCE.Tensor(3,mop)
        data["T4"] = _PCE.Tensor(4,mop)
        data["mop"] = mop
        for nw in 1:Npce, nd in 1:Nd
            data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
            data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
        end
    
        for nw in 1:Npce, nd_g in 1:Nd_g
    
            data["nw"]["$nw"]["RES"]["$nd_g"]["pd"] = pd_g[nd_g,nw]
            # data["nw"]["$nw"]["RES"]["$nd_g"]["qd"] = qd_g[nd_g,nw]
            data["nw"]["$nw"]["RES"]["$nd_g"]["qd"] = 0
        end
    
        return data


end




# output data
"""
    StochasticPowerModels.pce_coeff(result, element::String, id::Int, var::String)

Returns all polynomial chaos coefficients associated with the variable `var` of 
the `id`th element `element`.
"""
pce_coeff(result, element::String, id::Int, var::String) =
    [nw[2][element]["$id"][var] for nw in sort(collect(result["solution"]["nw"]), by=x->parse(Int,x[1]))]

"""
    StochasticPowerModels.sample(sdata, result, element::String, id::Int, var::String; sample_size::Int=1000)

Return an `sample_size` sample of the variable `var` of the `id`th element 
`element`.
"""
sample(result, element::String, id::Int, var::String; sample_size::Int=1000) =
    _PCE.samplePCE(sample_size, pce_coeff(result, element, id, var), result["mop"])

"""
    StochasticPowerModels.density(sdata, result, element::String, id::Int, var::String; sample_size::Int=1000)

Return an kernel density estimate of the variable `var` of the `id`th element 
`element`.
"""
density(result, element::String, id::Int, var::String; sample_size::Int=1000) =
    _KDE.kde(sample(result, element, id, var; sample_size=sample_size))

function print_summary(obj::Dict{String,<:Any}; kwargs...)
    if _IM.ismultinetwork(obj)
        for (n,nw) in obj["nw"]
            println("----------------")
            println("PCE index $n")
            _PM.summary(stdout, nw; kwargs...)
        end
    end
end




function dimension_management_PCE(sdata::Dict{String,Any})

    PCE_data = Dict(
        "T2" => sdata["T2"],
        "T3" => sdata["T3"],
        "T4" => sdata["T4"],
        "mop" => sdata["mop"],
    )
    num_of_coeff = PCE_data["mop"].dim

    delete!(sdata, "T2")
    delete!(sdata, "T3")
    delete!(sdata, "T4")
    delete!(sdata, "mop")

    _FP.add_dimension!(sdata, :PCE_coeff, num_of_coeff; metadata = PCE_data)

    return sdata
end



