using Pkg
Pkg.activate(".")
# Pkg.instantiate()

using JuMP
using Ipopt
using PowerModels
using Memento
using InfrastructureModels
using StochasticPowerModels
using PowerModelsACDC
using FlexPlan

# using HSL_jll


#Constants
const _PM = PowerModels
const _SPM = StochasticPowerModels
const _PMACDC = PowerModelsACDC
const _FP = FlexPlan

Memento.setlevel!(Memento.getlogger(StochasticPowerModels), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsACDC), "error")

#Solver inputs
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>1000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-4,
                                                               # "hsllib" => HSL_jll.libhsl_path,
                                                               # "linear_solver" => "ma27",
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

#gPC degree input
deg  = 2

#RES input
pen_level = 0.4

#Case file and data reading
case = "case5_ACDC_SPM_95cc.m"
file  = joinpath(BASE_DIR, "test/data/matpower", case)
data = _PM.parse_file(file)
_PMACDC.process_additional_data!(data)

s = Dict("RES Curtailment" => true,
         "Load Curtailment" => false,
         )


total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus

sdata = _SPM.build_stochastic_data_ACDC_RES(data, deg, p_size)
_SPM.dimension_management_PCE(sdata)
sdata["curtailment"] = s


println("\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_spm_VaR = _SPM.solve_sopf_iv_acdc_VaR(sdata, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm_VaR["primal_status"]), ")")


# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm_VaR["primal_status"])
print("Objective: ")
print(result_spm_VaR["objective"])
print("\nSolve Time: ")
print(result_spm_VaR["solve_time"])