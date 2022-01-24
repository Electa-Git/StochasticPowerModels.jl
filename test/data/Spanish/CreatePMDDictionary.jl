# using JSON
# using DataFrames
# using CSV
# using PowerModelsDistribution
# using PowerModels
# using Ipopt
# using Statistics
# using Plots
"""
Parser for distribution network model in .JSON format to
PowerModelsDistribution MATHEMATICAL formultation (more info on this format:
https://lanl-ansi.github.io/PowerModelsDistribution.jl/latest/math-model/).
Time series analysis can be performed (which is currently not natively supported
in PowerDistributionModels, as of sept 2020)
.JSON files should be in a folder, the directory of the folder should be specified below
the configuration file of the feeder to be analyzed as well
Date: 16 sept 2020
Author: Alexander Hoogsteyn
"""
dir = "C:/Users/karpan/git/summer-job-Alexander/"
config_file_name = "POLA/65019_74469_configuration.json"


const _PMD = PowerModelsDistribution
time_steps=24
solver =Ipopt.Optimizer

network = build_mathematical_model(dir,feeder)
result=_PMD.solve_mc_pf(network,ACPUPowerModel, solver)
"""
Specify voltage and power base, voltage base should be the phase-to-ground voltage
of the feeder studied (kV), the power base can be arbitrairly chosen (MW)
"""
const voltage_base = 0.230  # (kV)
const power_base = 0.5  # (MW)
const Z_base = voltage_base^2/power_base # (Ohm)
const current_base = power_base/(voltage_base*1e-3) # (A)

const mwpu = 1/power_base
const kwpu = (1e-3)/power_base


"""
Function that builds a network model in the PowerModelsDistribution format
(mathematical) from a JSON file name & location. Three phase active & reactive power
for all devices in the network can be set using pd & qd respectively (in pu).
"""
function build_mathematical_model(dir, config_file_name; pd = [0.0, 0.0, 0.0], qd = [0.0, 0.0, 0.0], scale_factor = 1.0)
    configuration = "star"
    network_model = Dict{String,Any}()
    configuration_json_dict = Dict{Any,Any}()

    network_model["is_kron_reduced"] = true
    network_model["p_factor"] = 0.95
    network_model["q_factor"] = sqrt(1-0.95^2)
    network_model["dcline"] = Dict{String,Any}()
    network_model["switch"] = Dict{String,Any}()
    network_model["is_projected"] = true
    network_model["per_unit"] = true
    network_model["data_model"] = MATHEMATICAL
    network_model["shunt"] = Dict{String,Any}()
    network_model["transformer"] = Dict{String,Any}()
    network_model["bus"] = Dict{String,Any}()
    network_model["map"] = Dict{String,Any}()
    network_model["conductors"] = 3
    network_model["baseMVA"] =  power_base
    network_model["basekv"] =  voltage_base
    network_model["bus_lookup"] = Dict{Any,Int64}()
    network_model["run_type"] = 1
    network_model["load"] = Dict{String,Any}()
    network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
    "pg"            => [1.0, 1.0, 1.0],
    "model"         => 2,
    "connections"   => [1, 2, 3],
    "shutdown"      => 0.0,
    "startup"       => 0.0,
    "configuration" => WYE,
    "name"          => "virtual_generator",
    "qg"            => [1.0, 1.0, 1.0],
    "gen_bus"       => 1,
    "vbase"         =>  voltage_base,
    "source_id"     => "virtual_generator",
    "index"         => 1,
    "cost"          => [0.0, 0.0],
    "gen_status"    => 1,
    "qmax"          => [1.0, 1.0, 1.0],
    "qmin"          => [-1.0, -1.0, -1.0],
    "pmax"          => [1.0, 1.0, 1.0],
    "pmin"          => [-1.0, -1.0, -1.0],
    "ncost"         => 2
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
    configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
    branches_file_name = configuration_json_dict["gridConfig"]["branches_file"]
    buses_file_name = configuration_json_dict["gridConfig"]["buses_file"]
    devices_file_name = configuration_json_dict["gridConfig"]["devices_file"]


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
                    "grounded"  => Bool[0, 0, 0],
                    "terminals" => [1, 2, 3],
                    "vbase"     =>  voltage_base,
                    "index"     => id,
                    "bus_i"     => id,
                    "vmin"      => [0.0, 0.0, 0.0],
                    "vmax"      => [1.5, 1.5, 1.5],
                    "va"        => [0.0, 2.0944, -2.0944],
                    "vm"        => [1.0, 1.0, 1.0], 
                    "LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "run_type"  => 1
                    )
            else
                network_model["bus"][id_s] = Dict{String,Any}(
                    "name"      => id_s,
                    "bus_type"  => 1,
                    "grounded"  => Bool[0, 0, 0],
                    "terminals" => [1, 2, 3],
                    "vbase"     =>  voltage_base,
                    "index"     => id,
                    "bus_i"     => id,
                    "vmin"      => [0.0, 0.0, 0.0],
                    "vmax"      => [1.5, 1.5, 1.5], 
                    "LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                    "run_type"  => 1
                    )
            end;
        end;
    end;
    load_st_file = feeder[1:length(feeder)-19]*".csv"
    open(dir * devices_file_name,"r") do io
    devices_json_dict = JSON.parse(io)
      for device in devices_json_dict["LVcustomers"]
        id = device["deviceId"] + 1 #Indexing starts at one in Julia
        id_s = string(id)
        cons = convert(Float64,device["yearlyNetConsumption"])
        network_model["load"][id_s] = Dict{String,Any}(
            "model"         => POWER,
            "connections"   => vec(Int.(device["phases"])),
            "configuration" => configuration=="star" ? WYE : DELTA,
            "name"          => id_s*"-"*device["coded_ean"],
            "status"        => 1,
            "vbase"         =>  voltage_base,
            "vnom_kv"       => 1.0,
            "source_id"     => device["coded_ean"],
            "load_bus"      => device["busId"] + 1,
            "dispatchable"  => 0,
            "index"         => id,
            "yearlyNetConsumption" => cons,
            "phases"        => device["phases"],
            "pd"            => pd,
            "qd"            => qd,
            "p_inj"         => fill(0.0, 3),
            "q_inj"         => fill(0.0, 3),
            "conn_cap_kW"   => device["connectionCapacity"]
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
    "BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 264,
    "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
    "BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
    "BT - RV 0,6/1 KV 4*25 KAL" => 100,
    "BT - RV 0,6/1 KV 4*50 KAL" => 150,
    "BT - RV 0,6/1 KV 4*95 KAL" => 230,
    "BT - RX 0,6/1 KV 2*16 Cu" => 95,
    "BT - RX 0,6/1 KV 2*2 Cu" => 30,
    "BT - RX 0,6/1 KV 2*4 Cu" => 40,
    "BT - RX 0,6/1 KV 2*6 Cu" => 50,
    "BT - RZ 0,6/1 KV 2*16 AL" => 75,
    "BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 264,
    "BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 264,
    "BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 78.98,
    "BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
    "BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 118.47,
    "BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
    "BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 182.21,
    "BT - RZ 0,6/1 KV 4*16 AL" => 75,
    "aansluitkabel" => 200
    )


        for branch in branches_json_dict
            id = branch["branchId"] +1
            id_s = string(id)
            network_model["branch"][id_s] = Dict{String,Any}(
                "shift"         => [0.0, 0.0, 0.0],
                "f_connections" => [1, 2, 3],
                "name"          => id_s,
                "switch"        => false,
                "g_to"          => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
                "c_rating_a"    => [0.8, 0.8, 0.8],
                "vbase"         =>  voltage_base,
                "g_fr"          => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
                "t_connections" => [1, 2, 3],
                "f_bus"         => branch["upBusId"]+1,
                "b_fr"          => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
                "c_rating_b"    => [0.8, 0.8, 0.8],
                "br_status"     => 1,
                "t_bus"         => branch["downBusId"]+1,
                "b_to"          => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
                "index"         => id,
                "angmin"        => [-1.0472, -1.0472, -1.0472],
                "angmax"        => [1.0472, 1.0472, 1.0472],
                "transformer"   => false,
                "tap"           => [1.0, 1.0, 1.0],
                "c_rating_c"    => [0.8, 0.8, 0.8],              
                )

            if haskey(impedance_dict,branch["cableType"])
                network_model["branch"][id_s]["br_x"] = impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base .* [2.0 1.0 1.0; 1.0 2.0 1.0; 1.0 1.0 2.0] 
                network_model["branch"][id_s]["br_r"] = impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base .* [2.0 1.0 1.0; 1.0 2.0 1.0; 1.0 1.0 2.0]
            end;
            

            if haskey(currentmax_dict,branch["cableType"])
                network_model["branch"][id_s]["rate_a"] = fill(((currentmax_dict[branch["cableType"]]*voltage_base)/1e3)/power_base, 3)
                
                network_model["branch"][id_s]["I_rating"] = fill(currentmax_dict[branch["cableType"]]/current_base, 3)

            end;
        end;
    end;
end;
    return network_model
end;




##not required for Rickard for now
""
function build_mathematical_model_single_phase(dir, config_file_name; pd = 0.0, qd = 0.0, scale_factor = 1.0)
#configuration = "star"
network_model = Dict{String,Any}()
configuration_json_dict = Dict{Any,Any}()

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
network_model["weight"] = Dict{String,Any}("1" => Dict{String,Any}(
    "w1"    => 0,
    "w2"    => 0,
    "w3"    => 0,
    "w4"    => 0,
    "w5"    => 0,
    "w6"    => 1,
    "w7"    => 1,
    "flex_const" => 0
)
)
network_model["map"] = Dict{String,Any}()
#network_model["conductors"] = 1
network_model["baseMVA"] =  power_base
network_model["basekv"] =  voltage_base
network_model["bus_lookup"] = Dict{Any,Int64}()
network_model["run_type"] = 1
network_model["load"] = Dict{String,Any}()
network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
"pg"            => 1.0,
"model"         => 2,
#"connections"   => [1, 2, 3],
"shutdown"      => 0.0,
"startup"       => 0.0,
"configuration" => WYE,
"name"          => "virtual_generator",
"qg"            => 1.0,
"gen_bus"       => 1,
"vbase"         =>  voltage_base,
"source_id"     => "virtual_generator",
"index"         => 1,
"cost"          => [0.0, 0.0],
"gen_status"    => 1,
"qmax"          => 1.0,
"qmin"          => -1.0,
"pmax"          => 1.0,
"pmin"          => -1.0,
"ncost"         => 2
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
configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
branches_file_name = configuration_json_dict["gridConfig"]["branches_file"]
buses_file_name = configuration_json_dict["gridConfig"]["buses_file"]
devices_file_name = configuration_json_dict["gridConfig"]["devices_file"]


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
                "vmin"      => 0.0,
                "vmax"      => 1.5,
                "va"        => 0.0,
                "vm"        => 1.01, 
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
                "vmin"      =>  0.0,
                "vmax"      => 1.5, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "run_type"  => 1
                )
        end;
    end;
end;

open(dir * devices_file_name,"r") do io
devices_json_dict = JSON.parse(io)
  for device in devices_json_dict["LVcustomers"]
    id = device["deviceId"] + 1 #Indexing starts at one in Julia
    id_s = string(id)
    cons = convert(Float64,device["yearlyNetConsumption"])
    network_model["load"][id_s] = Dict{String,Any}(
        "model"         => POWER,
        #"connections"   => vec(Int.(device["phases"])),
        "configuration" => configuration=="star" ? WYE : DELTA,
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
        "pd"            => pd,
        "qd"            => qd,
        "p_inj"         => 0.0,
        "q_inj"         => 0.0,
        "conn_cap_kW"   => device["connectionCapacity"]
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
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 264,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
"BT - RV 0,6/1 KV 4*25 KAL" => 100,
"BT - RV 0,6/1 KV 4*50 KAL" => 150,
"BT - RV 0,6/1 KV 4*95 KAL" => 230,
"BT - RX 0,6/1 KV 2*16 Cu" => 95,
"BT - RX 0,6/1 KV 2*2 Cu" => 30,
"BT - RX 0,6/1 KV 2*4 Cu" => 40,
"BT - RX 0,6/1 KV 2*6 Cu" => 50,
"BT - RZ 0,6/1 KV 2*16 AL" => 75,
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 264,
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 264,
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 78.98,
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 118.47,
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 182.21,
"BT - RZ 0,6/1 KV 4*16 AL" => 75,
"aansluitkabel" => 200
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
            "c_rating_a"    => 0.8,
            "vbase"         =>  voltage_base,
            "g_fr"          => 0.0,
            #"t_connections" => [1, 2, 3],
            "f_bus"         => branch["upBusId"]+1,
            "b_fr"          => 0.0,
            "c_rating_b"    => 0.8,
            "br_status"     => 1,
            "t_bus"         => branch["downBusId"]+1,
            "b_to"          => 0.0,
            "index"         => id,
            "angmin"        => -1.0472,
            "angmax"        => 1.0472,
            "transformer"   => false,
            "tap"           => 1.0,
            "c_rating_c"    => 0.8,              
            )

        if haskey(impedance_dict,branch["cableType"])
            network_model["branch"][id_s]["br_x"] = impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            network_model["branch"][id_s]["br_r"] = impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
        end;
        

        if haskey(currentmax_dict,branch["cableType"])
            network_model["branch"][id_s]["rate_a"] = ((currentmax_dict[branch["cableType"]]*voltage_base)/1e3)/power_base
            
            network_model["branch"][id_s]["I_rating"] = currentmax_dict[branch["cableType"]]/current_base

        end;
    end;
end;
end;
return network_model
end;

"""
Returns an array containing an estimated load_profile based on the households
mean power consumption. The array whithin reference_profiles with the nearest
mean to target_mean is chosen. Units of target_mean and samples within
reference_profiles should be the same (ussualy kWh).
"""
function pick_load_profile(target_mean,reference_profiles)
    load_profile = zeros(Float64,length(reference_profiles[1]))
    smallest_error = Inf
    for i in reference_profiles
        error = abs(target_mean-mean(i))
        if error < smallest_error
            smallest_error = error
            load_profile  = i
        end;
    end;
    return load_profile
end;

"""
Reads out load profiles from CSV file
"""
function read_reference_profile_from_csv(csv_name,mean_power)
    csv_name_1 = dir*"Load profiles POLA/SM/profile_"*csv_name*"_2020.csv"
    if "profile_"*csv_name*"_2020.csv" in readdir(dir*"Load profiles POLA/SM/")
        df = DataFrame(CSV.File(csv_name_1,delim=",",header=0))
    else
        df=DataFrame(x=ones(481)*mean_power)
    end
        reference_profiles = []
    for i in names(df)
        push!(reference_profiles,df[!,i])
    end;
    return reference_profiles
end;


"""
Reads out load profiles from CSV file selecting representative days only.
"""
function read_reference_profile_from_csv_representative(csv_name, representative_days, day_time_steps)
    
    csv_path = dir*"Load profiles POLA/2020/profile_"*csv_name*"_2020.csv"
    yearlong_reference_profile = readdlm(csv_path)

    # if length(yearlong_reference_profile) == 8760
    #     println("Profile: one year, hourly time steps.")
    # elseif length(yearlong_reference_profile) == 8784
    #      println("Profile: leap year, hourly time steps.")
    # else
    #     println("Profile: undefined.")
    # end

    reference_profile = []

    for d in representative_days
        push!(reference_profile, yearlong_reference_profile[(day_time_steps*(d-1))+1:(day_time_steps*d)])
    end


    for xx in 1:length(reference_profile)
        for yy in 1:length(reference_profile[xx])
            if isnan(reference_profile[xx][yy]) == true
                println("Load profile $csv_name has Nan at time_step $xx -> replaced with previous time_step's load")
                reference_profile[xx][yy] = reference_profile[xx][yy-1]
            end
        end
    end

    return reference_profile
end;


findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))


"""
Reads out load profiles from CSV file selecting representative days only.
"""
function read_reference_profile_from_csv_representative_StROBe(device_name::String, conncap::Float64, profiles::Array, representative_days::Array, day_time_steps::Int64)
    
    tot_profiles = length(yearlong_profiles[1,:])
    max_yearlong_profiles = zeros(tot_profiles) #kW

    for i in 1:tot_profiles
        max_yearlong_profiles[i] = maximum(yearlong_profiles[:,i])
    end

    profile_num = findnearest(max_yearlong_profiles, conncap)[2]
  
    reference_profile = []

    for d in representative_days
        push!(reference_profile, yearlong_profiles[(day_time_steps*(d-1))+1:(day_time_steps*d), profile_num])
    end


    for xx in 1:length(reference_profile)
        for yy in 1:length(reference_profile[xx])
            if isnan(reference_profile[xx][yy]) == true
                println("Load profile $csv_name has Nan at time_step $xx -> replaced with previous time_step's load")
                reference_profile[xx][yy] = reference_profile[xx][yy-1]
            end
        end
    end  

    return reference_profile, profile_num
end;


"""
Reads out CSV files from dublin data set and returns it as reference profiles
that can be used to form a multinetwork model
"""
function read_dublin_data()
    reference_profiles = []
    for file in 2:2
        csv_name = dir*"DUBLIN/38_CER Electricity_Gas/File"*string(file)*".txt"
        df = DataFrame(CSV.File(csv_name,delim=" ",header=["id","time","power"]))
        data = unstack(df,:time,:id,:power)
        for i in propertynames(data)
            if findfirst(ismissing,data[1:17520,i]) == nothing
                push!(reference_profiles, data[1:17520,i])
            end;
        end;
    end;
    if length(reference_profiles) == 0
        print("NO REFERENCE PROFILES FOUND")
    end;
    return reference_profiles
end;

"""
Builds a multinetwork model in a format similar to the one used in PowerModels
but extended to be able to account for phase inbalance such as in PowerModelsDistribution
Additionaly, it needs reference load profiles to attach time series data to the customers
reference_profiles should by a Array that contains kWh use for each
half hour for an entire year. Optionally, a scaling factor for the profile can be given.
Optionally, the reference_profiles can be a mulidimentional array which contains a
series of load profile Arrays, then a profile that is closest to the customer
is selected (based on the total consumption)
"""
function build_mn_mathematical_model(dir,feeder,reference_profiles,time_steps,scale_factor=1.0,time_unit=0.5)
    if eltype(reference_profiles) == Float64    #If no multidimensional array is given
        reference_profiles = [reference_profiles]   #make it into one
    end;

    network_model = build_mathematical_model(dir,feeder)
    mn_model = replicate(network_model,time_steps,global_keys=Set{String}())
    mn_model["data_model"] = MATHEMATICAL
    for (id_s,device) in network_model["load"]
        mean_power = device["yearlyNetConsumption"]*time_unit/365/24 #Assuming yearlyNetConsumption contains total consumption of 365 days in kWh
        load_profile = pick_load_profile(mean_power,reference_profiles)   #Pick the best fit load profile
        load_profile = scale_factor*load_profile/1000/time_unit #scale from kWh to MW
        load_profile = load_profile/power_base #convert to per-uits
        print(load_profile)
        for step in 1:time_steps
            pd = load_profile[step]
            qd = pd/20
            if length(device["phases"]) == 3   #Three phase connection
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd .* [0.33, 0.33, 0.33]
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd .* [0.33, 0.33, 0.33]
            elseif device["phases"][1] == 1   #Connected to phase 1
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd .* [0.0, 1.0, 0.0]
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd .* [0.0, 1.0, 0.0]
            elseif device["phases"][1] == 2   #Connected to phase 2
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd .* [1.0, 0.0, 0.0]
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd .* [1.0, 0.0, 0.0]
            elseif device["phases"][1] == 3   #Connected to phase 3
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd .* [0.0, 0.0, 1.0]
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd .* [0.0, 0.0, 1.0]
            else
                print(device["phases"] * " is an unknown phase connection")
            end;
        end;
    end;
    return mn_model
end;

"""
This function to be used when loadprofile is in csv for each device
"""
function build_mn_mathematical_model_from_csv(dir,feeder,time_steps; scale_factor=1.0,time_unit=0.5)
    network_model = build_mathematical_model(dir,feeder)
    mn_model = PowerModels.replicate(network_model,time_steps,global_keys=Set{String}())
    mn_model["data_model"] = MATHEMATICAL
    for (id_s,device) in network_model["load"]
        mean_power = device["yearlyNetConsumption"]*time_unit/365/24 #Assuming yearlyNetConsumption contains total consumption of 365 days in kWh
        load_profile = read_reference_profile_from_csv(device["source_id"],mean_power)   #Pick the profile from csv
        load_profile = scale_factor*load_profile/1000/time_unit #scale from kWh to MW
        load_profile = load_profile/power_base #convert to per-uits

        for step in 1:time_steps
            pd = load_profile[1][step]
            qd = pd/20
            if length(device["phases"]) == 3   #Three phase connection
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd .* [0.33, 0.33, 0.33]
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd .* [0.33, 0.33, 0.33]
            elseif device["phases"][1] == 1   #Connected to phase 1
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd .* [0.0, 1.0, 0.0]
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd .* [0.0, 1.0, 0.0]
            elseif device["phases"][1] == 2   #Connected to phase 2
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd .* [1.0, 0.0, 0.0]
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd .* [1.0, 0.0, 0.0]
            elseif device["phases"][1] == 3   #Connected to phase 3
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd .* [0.0, 0.0, 1.0]
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd .* [0.0, 0.0, 1.0]
            else
                print(device["phases"] * " is an unknown phase connection")
            end;
        end;
    end;
    return mn_model
end;


"single phase eq"
function build_mn_single_phase_mathematical_model_from_csv(dir,feeder,time_steps; scale_factor=1.0,time_unit=0.5)
    network_model = build_mathematical_model_single_phase(dir,feeder)
    mn_model = PowerModels.replicate(network_model,time_steps,global_keys=Set{String}())
    #mn_model["data_model"] = MATHEMATICAL
    for (id_s,device) in network_model["load"]
        mean_power = device["yearlyNetConsumption"]*time_unit/365/24 #Assuming yearlyNetConsumption contains total consumption of 365 days in kWh
        load_profile = read_reference_profile_from_csv(device["source_id"],mean_power)   #Pick the profile from csv
        load_profile = scale_factor*load_profile/1000/time_unit #scale from kWh to MW
        load_profile = load_profile/power_base #convert to per-uits

        for step in 1:time_steps
            pd = load_profile[1][step]
            qd = pd/20
            #if length(device["phases"]) == 3   #Three phase connection
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd 
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd 
            
        end;
    end;
    return mn_model
end;


"single phase eq"
function build_mn_single_phase_mathematical_model_scenario(dir,feeder,time_steps;scenario=1, scale_factor=1.0,time_unit=1.0)
    network_model = build_mathematical_model_single_phase(dir,feeder)
    mn_model = PowerModels.replicate(network_model,time_steps,global_keys=Set{String}())
    #mn_model["data_model"] = MATHEMATICAL
    for (id_s,device) in network_model["load"]
        #mean_power = device["yearlyNetConsumption"]*time_unit/365/24 #Assuming yearlyNetConsumption contains total consumption of 365 days in kWh
        load_profile = read_reference_profile_from_scenarios(string(device["index"]-1),scenario=scenario)   #Pick the profile from csv
        load_profile = scale_factor*load_profile/1000/time_unit #scale from kWh to MW
        load_profile = load_profile/power_base #convert to per-uits
        pv_prof = read_reference_PV_profile_from_scenarios("SolarPV_normalized_profiles"; scenario=scenario)
        pv_size=[5,6,7,8]/1000/power_base
        st_time=[6,7,8,9]

        for step in 1:time_steps
            pd = load_profile[1][step]
            qd = pd/20
            #if length(device["phases"]) == 3   #Three phase connection
                mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd 
                mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd 
            
        end;


        
        ind=parse(Int,id_s)%4+1
        #[mn_model["nw"]["$j"]["load"]["$i"]["pd"]=mn_model["nw"]["$j"]["load"]["$i"]["pd"]/12  for j=st_time[ind]:time_steps-8]
        [mn_model["nw"]["$j"]["load"]["$id_s"]["pd"]=mn_model["nw"]["$j"]["load"]["$id_s"]["pd"]-pv_prof[1][j]*pv_size[ind]  for j=1:time_steps]

    end;
    return mn_model
end;

"""
Reads out load profiles from CSV file
"""
function read_reference_profile_from_scenarios(csv_name; scenario=1)
    csv_name_1 = "C:/Users/karpan/Documents/PhD/flexplan_PSCC/Scenarios_RobustSolution/"*"House"*csv_name*".csv"
   
    df = DataFrame(CSV.File(csv_name_1,delim=",",header=0))
    

    reference_profiles = []
    i = "Column"*string(scenario)
    push!(reference_profiles,df[!,i])
    
    return reference_profiles
end;


"""
Reads out PV profiles from CSV file for each scenario
"""
function read_reference_PV_profile_from_scenarios(csv_name; scenario=1)
    csv_name_1 = "C:/Users/karpan/Documents/PhD/flexplan_PSCC/Scenarios_RobustSolution/"*csv_name*".csv"
   
    df = DataFrame(CSV.File(csv_name_1,delim=",",header=0))
    

    reference_profiles = []
    i = "Column"*string(scenario)
    push!(reference_profiles,df[!,i])
    
    return reference_profiles
end;

"""
This function to be used for single phase equivalent when loadprofile is in csv for each device using also representative days.
"""
function build_mn_mathmodel_csvrepresent_withmarket(dir::String, feed_id::String, path_weather_data::String, path_market_data::String; day_time_steps::Int64=24, scale_factor::Float64=1.0, time_unit::Float64=1.0)

    feeder =  string(feed_id, "_configuration.json");

    network_model = build_mathematical_model(dir,feeder)
    
    days_csv = path_weather_data*"decision_variables_short.csv"
    days_data = readdlm(days_csv, ',')
    representative_days = days_data[2:end, 1]

    n_times = day_time_steps*length(representative_days)
    mn_model = PowerModels.replicate(network_model,n_times,global_keys=Set{String}())

    csv_path = "D:/1_03_20/_Genie Timeline/0/C/Users/u0126211/Documents/Power_Models/Finalizing/Clymans Wim/Load_Profiles_300_1hRes_StROBe.csv"
    yearlong_profiles = readdlm(csv_path, ';')[2:end, 2:end]./1e3 #kWh

    for (id_s,device) in network_model["load"]

        # load_profile = read_reference_profile_from_csv_representative(device["source_id"], representative_days, day_time_steps)   #Pick the profile from csv

        load_profile, profile_num = read_reference_profile_from_csv_representative_StROBe(device["source_id"], device["connectionCapacity"], yearlong_profiles, representative_days, day_time_steps)
        yearlong_profiles = yearlong_profiles[:, 1:end .!= profile_num]

        load_profile = load_profile.*(scale_factor/1000/time_unit) #scale from kWh to MW
        load_profile = load_profile./network_model["baseMVA"] #convert to per-units

        for day in 1:length(load_profile)
            for step in 1:length(load_profile[day])
                pd = load_profile[day][step]*0.99
                qd = load_profile[day][step]*sqrt(1-0.99^2)
                if length(device["phases"]) == 3   #Three phase connection
                    mn_model["nw"][string(24*(day-1)+step)]["load"][id_s]["pd"] = pd .* [0.33, 0.33, 0.33]
                    mn_model["nw"][string(24*(day-1)+step)]["load"][id_s]["qd"] = qd .* [0.33, 0.33, 0.33]
                elseif device["phases"][1] == 1   #Connected to phase 1
                    mn_model["nw"][string(24*(day-1)+step)]["load"][id_s]["pd"] = pd .* [1.0, 0.0, 0.0]
                    mn_model["nw"][string(24*(day-1)+step)]["load"][id_s]["qd"] = qd .* [1.0, 0.0, 0.0]
                elseif device["phases"][1] == 2   #Connected to phase 2
                    mn_model["nw"][string(24*(day-1)+step)]["load"][id_s]["pd"] = pd .* [0.0, 1.0, 0.0]
                    mn_model["nw"][string(24*(day-1)+step)]["load"][id_s]["qd"] = qd .* [0.0, 1.0, 0.0]
                elseif device["phases"][1] == 3   #Connected to phase 3
                    mn_model["nw"][string(24*(day-1)+step)]["load"][id_s]["pd"] = pd .* [0.0, 0.0, 1.0]
                    mn_model["nw"][string(24*(day-1)+step)]["load"][id_s]["qd"] = qd .* [0.0, 0.0, 1.0]
                else
                    print(device["phases"] * " is an unknown phase connection")
                end;
            end;
        end;


    end;

    mn_model = incorporate_market_data(path_weather_data, path_market_data, mn_model)

    return mn_model
end;




"""
Alias for when a single array is passed instead of a multidimentional array
"""
function build_mn_mathematical_model(dir,feeder,reference_profiles::Array{Float64,1},kwargs...)
    return build_mn_mathematical_model(dir,feeder,reference_profiles=[reference_profiles],kwargs...)
end;

"""
Variant on build_mc_pf() from PowerModelsDistribution that allows you to build
multinetworks, similar to how it is implemented in build_mn_pf() from PowerModels
"""
function build_mn_mc_pf(pm::AbstractPowerModel)
    for (n, network) in nws(pm)
        variable_mc_bus_voltage(pm, nw=n, bounded=false)
        variable_mc_branch_power(pm,nw=n,bounded=false)
        variable_mc_transformer_power(pm, nw=n, bounded=false)
        variable_mc_gen_power_setpoint(pm, nw=n, bounded=false)
        variable_mc_load_setpoint(pm, nw=n, bounded=false)
        variable_mc_storage_power(pm, nw=n, bounded=false)
        constraint_mc_model_voltage(pm,nw=n)
        for (i,bus) in ref(pm, :ref_buses, nw=n)
            @assert bus["bus_type"] == 3
            constraint_mc_theta_ref(pm, i, nw=n)
            constraint_mc_voltage_magnitude_only(pm, i, nw=n)
        end;
        # gens should be constrained before KCL, or Pd/Qd undefined
        for id in ids(pm, :gen, nw=n)
        constraint_mc_gen_setpoint(pm, id, nw=n)
        end;
        # loads should be constrained before KCL, or Pd/Qd undefined
        for id in ids(pm, :load, nw=n)
            constraint_mc_load_setpoint(pm, id, nw=n)
        end;
        for (i,bus) in ref(pm, :bus, nw=n)
            constraint_mc_load_power_balance(pm, i, nw=n)
            # PV Bus Constraints
            if length(ref(pm, :bus_gens, i, nw=n)) > 0 && !(i in ids(pm,:ref_buses, nw=n))
                # this assumes inactive generators are filtered out of bus_gens
                @assert bus["bus_type"] == 2

                constraint_mc_voltage_magnitude_only(pm, i, nw=n)
                for j in ref(pm, :bus_gens, i, nw=n)
                    constraint_mc_gen_power_setpoint_real(pm, j, nw=n)
                end;
            end;
        end;
        for i in ids(pm, :storage, nw=n)
            constraint_storage_state(pm, i, nw=n)
            constraint_storage_complementarity_nl(pm, i, nw=n)
            constraint_mc_storage_losses(pm, i, nw=n)
            constraint_mc_storage_thermal_limit(pm, i, nw=n)
        end;
        for i in ids(pm, :branch, nw=n)
            constraint_mc_ohms_yt_from(pm, i, nw=n)
            constraint_mc_ohms_yt_to(pm, i, nw=n)
        end;
        for i in ids(pm, :transformer, nw=n)
            constraint_mc_transformer_power(pm, i, nw=n)
        end;
    end;
end;
"Alias of build_mn_mc_pf() for when IVR formultation of pf is used"
function build_mn_mc_pf(pm::AbstractIVRModel)
    for (n, network) in nws(pm)
        # Variables
        variable_mc_bus_voltage(pm, nw=n, bounded = false)
        variable_mc_branch_current(pm, nw=n, bounded = false)
        variable_mc_transformer_current(pm, nw=n, bounded = false)
        variable_mc_gen_power_setpoint(pm, nw=n, bounded = false)
        variable_mc_load_setpoint(pm, nw=n, bounded = false)
        # Constraints
        for (i,bus) in ref(pm, :ref_buses, nw=n)
            @assert bus["bus_type"] == 3
            constraint_mc_theta_ref(pm, i, nw=n)
            constraint_mc_voltage_magnitude_only(pm, i, nw=n)
        end;
        # gens should be constrained before KCL, or Pd/Qd undefined
        for id in ids(pm, :gen, nw=n)
            constraint_mc_gen_setpoint(pm, id, nw=n)
        end;
        # loads should be constrained before KCL, or Pd/Qd undefined
        for id in ids(pm, :load, nw=n)
            constraint_mc_load_setpoint(pm, id, nw=n)
        end;
        for (i,bus) in ref(pm, :bus, nw=n)
            constraint_mc_load_current_balance(pm, i, nw=n)
            # PV Bus Constraints
            if length(ref(pm, :bus_gens, i, nw=n)) > 0 && !(i in ids(pm,:ref_buses, nw=n))
                # this assumes inactive generators are filtered out of bus_gens
                @assert bus["bus_type"] == 2
                constraint_mc_voltage_magnitude_only(pm, i, nw=n)
                for j in ref(pm, :bus_gens, i, nw=n)
                    constraint_mc_gen_power_setpoint_real(pm, j, nw=n)
                end;
            end;
        end;
        for i in ids(pm, :branch, nw=n)
            constraint_mc_current_from(pm, i, nw=n)
            constraint_mc_current_to(pm, i, nw=n)

            constraint_mc_bus_voltage_drop(pm, i, nw=n)
        end;

        for i in ids(pm, :transformer)
            constraint_mc_transformer_power(pm, i, nw=n)
        end;
    end;
end;

"Alias to run multinetwork powerflow using run_mc_model() from PowerModelsDistribution"
function run_mn_mc_pf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mn_mc_pf, multiconductor=true, multinetwork=true, kwargs...)
end;

"Run time series analysis of the specified feeder in POLA using Dublin data load profiles as reference"
function time_series_analysis()
    reference_profiles = read_dublin_data()
    multinetwork_model = build_mn_mathematical_model(dir,feeder,reference_profiles,48)
    solver = with_optimizer(Ipopt.Optimizer, print_level = 3, tol=1e-9)
    return run_mn_mc_pf(multinetwork_model,ACPPowerModel,solver)
end

"Run Time series analysis by simply rerunning PowerModelsDistribution for each time step "
function naive_time_series_analysis()
    #reference_profiles = read_dublin_data() #To pick one of the profile from DUBLIN LOADprofile data
    reference_profile=
    multinetwork_model = build_mn_mathematical_model(dir,feeder,reference_profiles,48)
    solver = with_optimizer(Ipopt.Optimizer, print_level = 1, tol=1e-9)
    result = Dict{String,Any}()
    for (n,network) in multinetwork_model["nw"]
        network["per_unit"] = true
        result[n] = run_mc_pf(network,ACPPowerModel, solver)
    end;
    return result
end;

"TOD: This function is yet not working."
function naive_time_series_analysis_from_csv()
    #reference_profiles = read_dublin_data() #To pick one of the profile from DUBLIN LOADprofile data
    multinetwork_model = build_mn_mathematical_model_from_csv(dir,feeder,48)
    solver = with_optimizer(Ipopt.Optimizer, print_level = 1, tol=1e-9)
    result = Dict{String,Any}()
    for (n,network) in multinetwork_model["nw"]
        network["per_unit"] = true
        result[n] = run_mc_pf(network,ACPPowerModel, solver)
    end;
    return result
end;



"""
This function to be used for flexplan_energy for each device
"""
function build_mn_mathematical_model_from_scenarios_with_pp_plus(dir,feeder; time_steps=24, scale_factor=1.0, time_unit=1.0, scenario=1)
    mn_model = build_mn_single_phase_mathematical_model_scenario(dir, feeder,time_steps; scenario=scenario, scale_factor=1.0,time_unit=1.0)
    #mn_model = replicate(network_model,time_steps)#,global_keys=Set{String}())
    
    
    solver = Ipopt.Optimizer
    len_nodes=length(mn_model["nw"]["1"]["bus"])
    result = Dict{String,Any}()
    result_pf = Dict{String,Any}()
    voltage_5= zeros(time_steps,len_nodes)
    power_1=zeros(time_steps,len_nodes-1)
    loss_1=zeros(time_steps,len_nodes-1)
    for (id_s,device) in mn_model["nw"]["1"]["load"]
        for step in 1:time_steps
            mn_model["nw"]["$(step)"]["load"][id_s]["pd_s"] = mn_model["nw"]["$(step)"]["load"][id_s]["pd"]
            mn_model["nw"]["$(step)"]["load"][id_s]["qd_s"] = mn_model["nw"]["$(step)"]["load"][id_s]["qd"]
        end;
    end
   
        for (id_s,bus) in mn_model["nw"]["1"]["bus"]
            for step in 1:time_steps 
                mn_model["nw"]["$(step)"]["bus"][id_s]["pp_plus"] = 1
                mn_model["nw"]["$(step)"]["bus"][id_s]["pp_minus"] = -1
                mn_model["nw"]["$(step)"]["bus"][id_s]["pq_plus"] = 1
                mn_model["nw"]["$(step)"]["bus"][id_s]["pq_minus"] = -1
                mn_model["nw"]["$(step)"]["bus"][id_s]["pc_gen"] = 0
                mn_model["nw"]["$(step)"]["bus"][id_s]["pc_load"] = 0
            end;
        end;
    return mn_model
end;
