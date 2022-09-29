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


""
function parse_dst_beta(dst, pa, pb, deg)
    dst == "Beta"    && return _PCE.Beta01OrthoPoly(deg, pa, pb; Nrec=5*deg)
    #dst == "Normal"  && return _PCE.GaussOrthoPoly(deg; Nrec=5*deg)
    #dst == "Uniform" && return _PCE.Uniform01OrthoPoly(deg; Nrec=5*deg)
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
        display(pd[nd,nw])
        data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
        data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
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



"Converts JSON file of three phase DN to single phase equivalent"
function build_mathematical_model_single_phase(dir, config_file_name;t_s=52, pd = 0.0, qd = 0.0, scale_factor = 1.0, curt=0.0, cross_area_fact=1.0)
#configuration = "star"

"""
Specify voltage and power base, voltage base should be the phase-to-ground voltage
of the feeder studied (kV), the power base can be arbitrairly chosen (MW)
"""
voltage_base = 0.230  # (kV)
power_base = 0.5  # (MW)
Z_base = voltage_base^2/power_base # (Ohm)
current_base = power_base/(voltage_base*1e-3) # (A)

mwpu = 1/power_base
kwpu = (1e-3)/power_base

network_model = Dict{String,Any}()
configuration_json_dict = Dict{Any,Any}()
device_df=CSV.read(dir*config_file_name[1:length(config_file_name)-19]*".csv", DataFrame)

dist_lv=CSV.read(dir*"beta_lm_2016_8_6"*".csv", DataFrame)
# dist_pv=CSV.read(dir*"beta_pm_2016_8_6"*".csv", DataFrame)
dist_pv=CSV.read(dir*"beta_pm_2022_181"*".csv", DataFrame)
dist_pv_ts= dist_pv[in([t_s]).(dist_pv.timeslot),:]
dist_lv_ts=dist_lv[in([t_s]).(dist_lv.timeslot),:]

dist_lv_ts_feeder = dist_lv_ts[in(unique(device_df.category)).(dist_lv_ts.cluster),:]
s_dict=Dict()
i=1
for dist in eachrow(dist_lv_ts_feeder)
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = dist["cluster"]
    s["pa"]= dist["alpha"]
    s["pb"]= dist["beta"]
    s["pc"]= dist["lower"]
    s["pd"]= dist["lower"]+dist["upper"]
    s_dict[string(i)] = s
    i=i+1
end


##add Irradiance if day time or there is some Irradiance
if dist_pv_ts.upper[1]>0
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = 55
    s["pa"]= dist_pv_ts[!,"alpha"][1]
    s["pb"]= dist_pv_ts[!,"beta"][1]
    s["pc"]= dist_pv_ts[!,"lower"][1]
    s["pd"]= dist_pv_ts[!,"lower"][1]+dist_pv_ts[!,"upper"][1]
    s_dict[string(i)] = s
end


#network_model["is_kron_reduced"] = true
network_model["p_factor"] = 0.95
network_model["q_factor"] = sqrt(1-0.95^2)
network_model["dcline"] = Dict{String,Any}()
network_model["switch"] = Dict{String,Any}()
#network_model["is_projected"] = true
network_model["per_unit"] = true
#network_model["data_model"] = MATHEMATICAL
network_model["shunt"] = Dict{String,Any}()
network_model["transformer"] = Dict{String,Any}()
network_model["bus"] = Dict{String,Any}()
network_model["map"] = Dict{String,Any}()
#network_model["conductors"] = 1
network_model["baseMVA"] =  power_base
network_model["basekv"] =  voltage_base
network_model["bus_lookup"] = Dict{Any,Int64}()
network_model["run_type"] = 1
network_model["load"] = Dict{String,Any}()
network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
"pg"            => 0.2,
"model"         => 2,
#"connections"   => [1, 2, 3],
"shutdown"      => 0.0,
"startup"       => 0.0,
#"configuration" => WYE,
"name"          => "virtual_generator",
"qg"            => 0.0,
"gen_bus"       => 1,
"vbase"         =>  voltage_base,
"source_id"     => Any["gen",1],
"index"         => 1,
"cost"          => [20000.0, 1400.0, 0.0],
"gen_status"    => 1,
"qmax"          => 1.275,
"qmin"          => -1.275,
"pmax"          => 1.5,
"pmin"          => -1.5,
"ncost"         => 3,
"λpmin"         => 1.65, #1.03643 ,
"λpmax"         => 1.65, #1.03643 ,
"λqmin"         => 1.65, #1.03643 ,
"λqmax"         => 1.65 #1.03643
))
network_model["settings"] = Dict{String,Any}(
"sbase_default"        => power_base,
"vbases_default"       => Dict{String,Any}(), #No default is specified for now, since default is never used
"voltage_scale_factor" => 1E3, #Voltages are thus expressed in kV
"sbase"                => power_base,
"power_scale_factor"   => 1E6, #Power is expressed in MW
"base_frequency"       => 50.0 #Hertz
)
network_model["branch"] = Dict{String,Any}()
network_model["storage"] = Dict{String,Any}()
open(dir * config_file_name,"r") do io
configuration_json_dict = JSON.parse(io)
end;
#voltage_base = configuration_json_dict["gridConfig"]["basekV"]
#power_base = configuration_json_dict["gridConfig"]["baseMVA"]
sub_dir=splitpath(config_file_name)[1]
configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
branches_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["branches_file"])[2]
buses_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["buses_file"])[2]
devices_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["devices_file"])[2]


open(dir * buses_file_name,"r") do io
buses_json_dict = JSON.parse(io)
    for bus in buses_json_dict
        id = bus["busId"] + 1 #Indexing starts at one in Julia
        id_s = string(id)
        network_model["bus_lookup"][id_s] = id
        network_model["settings"]["vbases_default"][id_s] =  voltage_base

        if id == 1 #Settings for slack bus
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => "slack",
                "bus_type"  => 3,
                ##"grounded"  => Bool[0, 0, 0],
                #"terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65, #1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      => 1.0,
                "vmax"      => 1,
                "va"        => 0.0,
                "vm"        => 1, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"run_type"  => 1
                )
        else
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => id_s,
                "bus_type"  => 1,
                #"grounded"  => Bool[0, 0, 0],
                #"terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65,#1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      =>  0.95,
                "vmax"      => 1.05, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "run_type"  => 1
                )
        end;
    end;
end;

#print(device_df)
open(dir * devices_file_name,"r") do io
devices_json_dict = JSON.parse(io)
  for device in devices_json_dict["LVcustomers"]
    id = device["deviceId"] + 1 #Indexing starts at one in Julia
    d=device_df[in(id-1).(device_df.dev_id),:]
    id_s = string(id)
    μ = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"lower"][1]
    σ  = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"upper"][1] 
    cons = convert(Float64,device["yearlyNetConsumption"])
    network_model["load"][id_s] = Dict{String,Any}(
        #"connections"   => vec(Int.(device["phases"])),
        "name"          => id_s*"-"*device["coded_ean"],
        "status"        => 1,
        "vbase"         =>  voltage_base,
        "vnom_kv"       => 1.0,
        "source_id"     => device["coded_ean"],
        "load_bus"      => device["busId"] + 1,
        "dispatchable"  => 0,
        "index"         => id,
        "yearlyNetConsumption" => cons,
        #"phases"        => device["phases"],
        "pd"            => max(0.1, μ)/1e3/power_base/ 3,
        "qd"            =>  max(0.01,μ)/1e3/ power_base/ 3/10,
        "p_inj"         => 0.0,
        "q_inj"         => 0.0,
        "conn_cap_kW"   => device["connectionCapacity"],
        "dst_id" => d[!,"category"][1],
        "cluster_id"  => findall(x->x==1,[s_dict["$i"]["dst_id"]==d[!,"category"][1] for i=1:length(s_dict)])[1],
        "μ"  => μ,
        "σ"  => σ 

    )
    end;
end;

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
impedance_dict = Dict{String,Any}(
"BT - Desconocido BT" => [0.21, 0.075],
"BT - MANGUERA" => [0.3586, 0.089], #[1.23, 0.08],
"BT - RV 0,6/1 KV 2*16 KAL" => [2.14, 0.09], #16 = 20 A
"BT - RV 0,6/1 KV 2*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => [0.2309, 0.085],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 4*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 4*50 KAL" => [0.71849, 0.093],
"BT - RV 0,6/1 KV 4*95 KAL" => [0.3586, 0.089],
"BT - RX 0,6/1 KV 2*16 Cu" => [1.23, 0.08],
"BT - RX 0,6/1 KV 2*2 Cu" => [9.9, 0.075],
"BT - RX 0,6/1 KV 2*4 Cu" => [4.95, 0.075],
"BT - RX 0,6/1 KV 2*6 Cu" => [3.3, 0.075],
"BT - RZ 0,6/1 KV 2*16 AL" => [2.14, 0.09],
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => [1.34, 0.097],
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => [0.9073, 0.095],
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => [0.718497, 0.093],
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => [0.4539, 0.091],
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => [0.3586, 0.089],
"BT - RZ 0,6/1 KV 4*16 AL" => [2.14, 0.09],
"aansluitkabel" => [1.15, 0.150]
)

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
currentmax_dict = Dict{String,Any}(
"BT - Desconocido BT" => 200,
"BT - MANGUERA" => 150, #200 certain  40.18#150
"BT - RV 0,6/1 KV 2*16 KAL" => 75,
"BT - RV 0,6/1 KV 2*25 KAL" => 100,
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 305,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
"BT - RV 0,6/1 KV 4*25 KAL" => 100,
"BT - RV 0,6/1 KV 4*50 KAL" => 150,
"BT - RV 0,6/1 KV 4*95 KAL" => 230,
"BT - RX 0,6/1 KV 2*16 Cu" => 75,#95,
"BT - RX 0,6/1 KV 2*2 Cu" => 40, #30,
"BT - RX 0,6/1 KV 2*4 Cu" => 60,#40,
"BT - RX 0,6/1 KV 2*6 Cu" => 80, #50,
"BT - RZ 0,6/1 KV 2*16 AL" => 20, #75,
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 305, #264,
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 305,#264,
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 100, #78.98,
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 150, #118.47,
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 230, # 182.21,
"BT - RZ 0,6/1 KV 4*16 AL" => 75,
"aansluitkabel" => 120 #200
)


    for branch in branches_json_dict
        id = branch["branchId"] +1
        id_s = string(id)
        network_model["branch"][id_s] = Dict{String,Any}(
            "shift"         => 0.0,
            #"f_connections" => [1, 2, 3],
            "name"          => id_s,
            "switch"        => false,
            "g_to"          => 0.0,
            "c_rating_a"    => currentmax_dict[branch["cableType"]],
            "vbase"         =>  voltage_base,
            "g_fr"          => 0.0,
            #"t_connections" => [1, 2, 3],
            "f_bus"         => branch["upBusId"]+1,
            "b_fr"          => 0.0,
            "c_rating_b"    => currentmax_dict[branch["cableType"]],
            "br_status"     => 1,
            "t_bus"         => branch["downBusId"]+1,
            "b_to"          => 0.0,
            "index"         => id,
            "angmin"        => -1.0472,
            "angmax"        => 1.0472,
            "transformer"   => false,
            "tap"           => 1.0,
            "c_rating_c"    => currentmax_dict[branch["cableType"]], 
            "λcmax"         => 1.65 #1.65 #2.5 #1.03643    
            )

        if haskey(impedance_dict,branch["cableType"])
            # network_model["branch"][id_s]["br_x"] = (cross_area_fact^(1/4)).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            # network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base
            network_model["branch"][id_s]["br_x"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base  
        end;
        

        if haskey(currentmax_dict,branch["cableType"])
            network_model["branch"][id_s]["rate_a"] = cross_area_fact.*((currentmax_dict[branch["cableType"]]*voltage_base)/1e3)/power_base
            
            #network_model["branch"][id_s]["I_rating"] = currentmax_dict[branch["cableType"]]/current_base

        end;
    end;
end;
end;


network_model["sdata"]= s_dict 
network_model["curt"]= curt

network_model["PV"]=deepcopy(network_model["load"]);
[network_model["PV"][d]["μ"]=s_dict[string(length(s_dict))]["pc"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["σ"]=s_dict[string(length(s_dict))]["pd"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["pd"]=s_dict[string(length(s_dict))]["pd"]/1e6/ power_base / 3 for d in   keys(network_model["PV"])]
return network_model
end;


"""
    StochasticPowerModels.build_stochastic_data_hc(data::Dict{String,Any}, deg::Int)

Function to build the multi-network data representative of the polynomial chaos
expansion of a single-network data dictionary for HC problem in DN.
"""
function build_stochastic_data_hc(data::Dict{String,Any}, deg::Int, t_s=50)
    # add maximum current
    curt=data["curt"]
    for (nb, branch) in data["branch"]
        f_bus = branch["f_bus"]
        branch["cmax"] = branch["rate_a"] / data["bus"]["$f_bus"]["vmin"]
    end


    # build mop
    opq = [parse_dst_beta(ns[2]["dst"], ns[2]["pa"], ns[2]["pb"], deg) for ns in sort(data["sdata"])]
    mop = _PCE.MultiOrthoPoly(opq, deg)

    # build load matrix
    Nd, Npce = length(data["load"]), mop.dim
    pd, qd = zeros(Nd, Npce), zeros(Nd, Npce)
    pd_g, qd_g = zeros(Nd, Npce), zeros(Nd, Npce)
    for nd in 1:Nd 
        # reactive power
        qd[nd,1] = data["load"]["$nd"]["qd"]
        # active power
        nb = data["load"]["$nd"]["load_bus"]
        ni = data["load"]["$nd"]["cluster_id"]
        if ni == 55
            pd[nd,1] = data["load"]["$nd"]["pd"]
        else
            base = data["baseMVA"]
            μ, σ = data["load"]["$nd"]["μ"] /1e3/ base/ 3, data["load"]["$nd"]["σ"] /1e3/ base/3
            if mop.uni[ni] isa _PCE.GaussOrthoPoly
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            else
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            end
        end
        np = length(opq)
        base = data["baseMVA"]
        μ, σ = data["PV"]["1"]["μ"]/1e6 / base / 3, data["PV"]["1"]["σ"] /1e6/ base / 3
        
            if mop.uni[np] isa _PCE.GaussOrthoPoly
                pd_g[nd,[1,np+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np])
            else
                pd_g[nd,[1,np+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np])
            end
        
    end

    # replicate the data
    data = _PM.replicate(data, Npce)

    # add the stochastic data 
    data["T2"] = _PCE.Tensor(2,mop)
    data["T3"] = _PCE.Tensor(3,mop)
    data["T4"] = _PCE.Tensor(4,mop)
    data["mop"] = mop
    data["curt"]= curt
    for nw in 1:Npce, nd in 1:Nd
       
        data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
        data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
    end
    
    for nw in 1:Npce, nd in 1:Nd

        data["nw"]["$nw"]["PV"]["$nd"]["pd"] = pd_g[nd,nw]
        data["nw"]["$nw"]["PV"]["$nd"]["qd"] = qd_g[nd,nw]
    end

    return data
end
