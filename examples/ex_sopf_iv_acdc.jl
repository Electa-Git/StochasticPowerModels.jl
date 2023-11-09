using Pkg
Pkg.activate(".")

using JuMP
using Ipopt
using PowerModels
using Memento
using InfrastructureModels
using StochasticPowerModels
using PowerModelsACDC

#Constants
const _PM = PowerModels
const _SPM = StochasticPowerModels
const _PMACDC = PowerModelsACDC

Memento.setlevel!(Memento.getlogger(StochasticPowerModels), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsACDC), "error")

#Solver inputs
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>3000, 
                                                               "sb"=>"yes", 
                                                            #    "fixed_variable_treatment" => "relax_bounds",
                                                )

#gPC degree input
deg  = 2

#RES input
pen_level = 0.4

#Case file and data reading
case = "case5_ACDC_SPM_95cc_RES.m"
file  = joinpath(BASE_DIR, "test/data/matpower", case)
data = _PM.parse_file(file)


println("\nPenetration Level = $pen_level")

total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus

println("   Solution progress: Solving...")
result_spm = _SPM.solve_sopf_iv_acdc(file, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm["primal_status"]), ")")

# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm["primal_status"])
print("Objective: ")
print(result_spm["objective"])
print("\nSolve Time: ")
print(result_spm["solve_time"])       