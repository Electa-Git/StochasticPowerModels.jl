using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels


# constants 
const PM = PowerModels
const SPM = StochasticPowerModels

# solvers
ipopt_solver = Ipopt.Optimizer

# input
deg  = 2
aux  = true
red  = false
case = "case5_spm.m"

# data
file  = joinpath(BASE_DIR, "test/data/matpower", case)
data  = PM.parse_file(file)
#-----------------------------------
# run the convenience functions for stochastic OPF for IVR and ACR
result_ivr = run_sopf_iv(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
result_hc = run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)

@assert result_ivr["termination_status"] == PM.LOCALLY_SOLVED
@assert result_hc["termination_status"] == PM.LOCALLY_SOLVED
#the optimal objective values (expectation) are 
obj_ivr = result_ivr["objective"] 
obj_acr = result_hc["objective"] 

# print variables for all polynomial indices k
SPM.print_summary(result_ivr["solution"])

# print variables for a specific index k
k=1
SPM.print_summary(result_ivr["solution"]["nw"]["$k"])

# get polynomial chaos coefficients for specific component
pg_coeff = pce_coeff(result_ivr, "gen", 1, "pg") 

# obtain 10 samples of the generator active power output variable
pg_sample = sample(result_ivr, "gen", 1, "pg"; sample_size=10) 

# obtain an kernel density estimate of the generator active power output variable
pg_density = density(result_ivr, "gen", 1, "pg"; sample_size=10) 

#-----------------------------------
# alternatively, you can first read in PowerModels dict, 
# from a file with stochastic data extensions:
