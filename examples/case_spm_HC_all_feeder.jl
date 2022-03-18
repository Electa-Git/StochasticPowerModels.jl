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
dir="C://Users//karpan//.julia//dev//StochasticPowerModels//test//data//Spanish//All_feeder"
all_feeder=CSV.read(dir*"/ts_all_feeder.csv",DataFrame,header=["conf", "ts","HC1","HC2"])
all_feeder[!,"HC1"]=zeros(160)
all_feeder[!,"HC2"]=zeros(160)
hc1=[]
hc2=[]
for b in eachrow(all_feeder)
    feeder="All_feeder/"*b.conf
    file  = joinpath(BASE_DIR, "test/data/Spanish/")
    data  = SPM.build_mathematical_model_single_phase(file, feeder, t_s= b.ts)

    s2 = Dict("output" => Dict("duals" => true))
    result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
    e=1;
    m=1

    for i=1:length(data["bus"])
        if -result_hc["solution"]["nw"]["1"]["bus"]["$i"]["dual_voltage_max"]>500
            l_old=data["bus"]["$i"]["λvmax"]
            m=sample(result_hc, "bus", i, "vm"; sample_size=100000)
            if quantile(m,[0.95])[1]>1.001*data["bus"]["$i"]["vmax"]
                data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]-mean(m))/std(m)+0.3
            elseif quantile(m,[0.95])[1]>1.0001*data["bus"]["$i"]["vmax"]
                    data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]-mean(m))/std(m)+0.1
            else
                data["bus"]["$i"]["λvmax"]= 2*l_old-(data["bus"]["$i"]["vmax"]-mean(m))/std(m)
            end
        end
    end


    result_hc_1= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)

    if result_hc_1["termination_status"]== PM.LOCALLY_SOLVED
        push!(hc1,result_hc["objective"])
        push!(hc2,result_hc_1["objective"])
    else
        push!(hc1,-1)
        push!(hc2,-1)
    end
end
all_feeder[!,"HC1"]=hc1
all_feeder[!,"HC2"]=hc2


