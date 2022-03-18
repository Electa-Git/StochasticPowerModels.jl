using Pkg
Pkg.activate(".")
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels
using PowerModelsDistribution
using JSON
using DataFrames
using CSV
using Statistics
using Plots
# constants 
const PM = PowerModels
const SPM = StochasticPowerModels

# solvers
ipopt_solver = Ipopt.Optimizer

# input
deg  = 2
aux  = true
red  = false

feeder = "POLA/1076069_1274129_mod_configuration.json" 


#feeder = "All_feeder/65025_80123_configuration.json"#1076069_1274125_configuration.json"

# data
file  = joinpath(BASE_DIR, "test/data/Spanish/")

data  = SPM.build_mathematical_model_single_phase(file, feeder, t_s= 59)
[data["PV"]["$i"]["p_size"]=1 for  i=1:length(data["load"])] 
result_pf= SPM.run_pf_deterministic(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red, stochastic=false)

"""
s2 = Dict("output" => Dict("duals" => true))
result_hc_2= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
e=1;
m=1
for i=1:length(data["bus"])
    if -result_hc_2["solution"]["nw"]["1"]["bus"]["$i"]["dual_voltage_max"]>500
        l_old=data["bus"]["$i"]["λvmax"]
        m=sample(result_hc_2, "bus", i, "vm"; sample_size=100000)
        if quantile(m,[0.95])[1]>1.001*data["bus"]["$i"]["vmax"]
            data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]-mean(m))/std(m)+0.3
        elseif quantile(m,[0.95])[1]>1.0001*data["bus"]["$i"]["vmax"]
                data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]-mean(m))/std(m)+0.1
        else
            data["bus"]["$i"]["λvmax"]= 2*l_old-(data["bus"]["$i"]["vmax"]-mean(m))/std(m)
        end
    end
end

[data["branch"]["$i"]["λcmax"]= 3 for i=1:18]
data["branch"]["13"]["λcmax"]= 2
result_hc_1= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)



while e>0.005
    mean_voltage=[result_hc["solution"]["nw"]["1"]["bus"]["$j"]["vm"] for j=1:length(data["bus"])]
    i=argmax(mean_voltage)
    samp = sample(result_hc, "bus", i, "vm"; sample_size=10000)
    if sum([a>data["bus"]["$i"]["vmax"] for a in samp])/10000 > 0.055
        if sum([a>data["bus"]["$i"]["vmax"] for a in samp])/10000 > 0.08
            [data["bus"]["$i"]["λvmax"]= data["bus"]["$i"]["λvmax"]+0.1 for i=45:52]
        else 
            [data["bus"]["$i"]["λvmax"]= data["bus"]["$i"]["λvmax"]+0.05 for i=1:length(data["bus"])]
        end
    elseif sum([a>data["bus"]["$i"]["vmax"] for a in samp])/10000 < 0.045
        if sum([a>data["bus"]["$i"]["vmax"] for a in samp])/10000 < 0.02
            [data["bus"]["$i"]["λvmax"]= data["bus"]["$i"]["λvmax"]-0.1 for i=1:length(data["bus"])]
        else 
            [data["bus"]["$i"]["λvmax"]= data["bus"]["$i"]["λvmax"]-0.05 for i=1:length(data["bus"])]
        end
    end

        result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
        mean_voltage=[result_hc["solution"]["nw"]["1"]["bus"]["$j"]["vm"] for j=1:length(data["bus"])]
        i=argmax(mean_voltage)
        e =  abs(0.05-sum([a>1.05 for a in sample(result_hc, "bus", i, "vm"; sample_size=10000)])/10000)
    print("iteration:$m")
    print(e)
    m=m+1
end

[data["bus"]["$i"]["λvmax"]=2.4 for i=1:63]
[data["branch"]["$i"]["λcmax"]=1.5 for i=1:63]

result_hc2_4= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)

a=[result_hc2_4["solution"]["nw"]["1"]["PV"]["$i"]["p_size"] for i=1:52]
hc2= Dict()
for i=30:80
 data  = SPM.build_mathematical_model_single_phase(file, feeder, t_s= i)
 [data["bus"]["$i"]["λvmax"]=2.4 for i=1:63]
 [data["branch"]["$i"]["λcmax"]=1.5 for i=1:63]
 result_hc2 = run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
 PV_HC = result_hc2["objective"]
 hc2[i]=PV_HC
end
#sdata = build_stochastic_data_hc(data, deg)
#remove the existing generators and keep onl;y in slack bus
#data["gen"] = Dict(k => v for (k, v) in data["gen"] if data["bus"][string(data["gen"][k]["gen_bus"])]["bus_type"] == 3)
#[data["bus"][k]["bus_type"]=1 for (k,v) in data["bus"] if data["bus"][k]["bus_type"]==2]

#-----------------------------------
# run the convenience functions for stochastic OPF for IVR and ACR
#result_ivr = run_sopf_iv(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
result_hc2 = run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)

@assert result_hc2["termination_status"] == PM.LOCALLY_SOLVED
#@assert result_hc["termination_status"] == PM.LOCALLY_SOLVED
#the optimal objective values (expectation) are 
obj_ivr = result_hc2["objective"] 
#bj_acr = result_hc["objective"] 

# print variables for all polynomial indices k
SPM.print_summary(result_hc2["solution"])

# print variables for a specific index k
k=1
SPM.print_summary(result_hc2["solution"]["nw"]["$k"]["PV"])

# get polynomial chaos coefficients for specific component
crd_pv = pce_coeff(result_hc2, "PV", 1, "crd_pv") 

# obtain 10 samples of the generator active power output variable
crd_sample = sample(result_hc2, "PV", 1, "crd_pv"; sample_size=10) 

# obtain an kernel density estimate of the PV active power output variable
crd_density = density(result_ivr, "PV", 1, "crd_pv"; sample_size=10) 

#-----------------------------------
# alternatively, you can first read in PowerModels dict, 
# from a file with stochastic data extensions:
"""