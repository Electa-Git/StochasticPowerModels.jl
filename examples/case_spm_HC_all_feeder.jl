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
r  = false
dir="C://Users//karpan//.julia//dev//StochasticPowerModels//test//data//Spanish//All_feeder"
all_feeder=CSV.read(dir*"/ts_all_feeder.csv",DataFrame,header=["conf", "ts","HC1","HC2"])
all_feeder[!,"HC1"]=zeros(nrow(all_feeder))
all_feeder[!,"HC2"]=zeros(nrow(all_feeder))
all_feeder[!,"HC3"]=zeros(nrow(all_feeder))
hc1=[]
hc2=[]
hc3=[]
t_cc=[]
t_opf=[]
global i=0
for b in eachrow(all_feeder)
    feeder="All_feeder/"*b.conf
    i=i+1
    print("Feeder no: $i \n")
    #feeder="All_feeder/"*all_feeder[1,"conf"]
    file  = joinpath(BASE_DIR, "test/data/Spanish/")
    data  = SPM.build_mathematical_model_single_phase(file, feeder, t_s= 59)
    [data["bus"]["$i"]["vmin"]=0.9 for i=1:length(data["bus"])]
    s2 = Dict("output" => Dict("duals" => true))
    result_hc_2= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2)
    e=1;
    m=1

    if result_hc_2["termination_status"]== PM.LOCALLY_SOLVED
    for i=1:length(data["bus"])
        if -result_hc_2["solution"]["nw"]["1"]["bus"]["$i"]["dual_voltage_max"]>500
            l_old=data["bus"]["$i"]["λvmax"]
            m=sample(result_hc_2, "bus", i, "vs"; sample_size=100000)
            if quantile(m,[0.95])[1]>1.001*data["bus"]["$i"]["vmax"]^2
                data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]^2-mean(m))/std(m)+0.3
            elseif quantile(m,[0.95])[1]>1.0001*data["bus"]["$i"]["vmax"]^2
                    data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]^2-mean(m))/std(m)+0.1
            elseif quantile(m,[0.95])[1]< 0.99*data["bus"]["$i"]["vmax"]^2
                data["bus"]["$i"]["λvmax"]= l_old-0.8
            elseif quantile(m,[0.95])[1]< 1*data["bus"]["$i"]["vmax"]^2
                data["bus"]["$i"]["λvmax"]= l_old-0.2
            end
        end
    end

    result_hc_1= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2)

    
        push!(hc1,result_hc_2["objective"])
        push!(hc2,result_hc_1["objective"])
        push!(t_cc, result_hc_1["solve_time"])
    else
        push!(hc1,-1)
        push!(hc2,-1)
        push!(t_cc, -1)
    end
    result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2, stochastic=false)
    
    if result_hc["termination_status"]== PM.LOCALLY_SOLVED
        push!(hc3,result_hc["objective"])
        push!(t_opf, result_hc["solve_time"])
    else
        push!(hc3,-1)
        push!(t_opf, result_hc["solve_time"])
    end
end
all_feeder[!,"HC1"]=hc1
all_feeder[!,"HC2"]=hc2
all_feeder[!,"HC3"]=hc3
all_feeder[!,"t_opf"]=t_opf
all_feeder[!,"t_cc"]=t_cc
CSV.write("PV_HC_new_lower_vmin.csv",all_feeder)
"""
#deterministic

for b in eachrow(all_feeder)
    feeder="All_feeder/"*b.conf
    file  = joinpath(BASE_DIR, "test/data/Spanish/")
    data  = SPM.build_mathematical_model_single_phase(file, feeder, t_s= b.ts)

    s2 = Dict("output" => Dict("duals" => true))
    result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2, stochastic=false)
    e=1;
    if result_hc["termination_status"]== PM.LOCALLY_SOLVED
        push!(hc3,result_hc["objective"])
        #push!(hc2,result_hc_1["objective"])
    else
        push!(hc3,-1)
        #push!(hc2,-1)
    end
end

all_feeder[!,"HC3"]=hc3

"""