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
using Plots

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
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, 
                                                               "max_iter"=>2000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-4,
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

#gPC degree input
deg  = 2

#RES input
pen_level = 0.3

#Case file and data reading
case = "case5_ACDC_mod_SPMACDC_95cc.m"
file  = joinpath(BASE_DIR, "test/data/matpower", case)
data = _PM.parse_file(file)
_PMACDC.process_additional_data!(data)

s = Dict("RES Curtailment" => false,
         "Load Curtailment" => false,
         )


total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus

sdata = _SPM.build_stochastic_data_ACDC_RES(data, deg, p_size)

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

sdata["curtailment"] = s


println("\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_spm = _SPM.solve_sopf_iv_acdc(sdata, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm["primal_status"]), ")")




data_ots = _PM.deepcopy(data)



for (b, branch) in data_ots["branch"]

   if parse(Int,b) == 2

      data_ots["branch"][b]["br_status_initial"] = 1

   else

      data_ots["branch"][b]["br_status_initial"] = 1

   end

end

for (b, branch) in data_ots["branchdc"]

   data_ots["branchdc"][b]["br_status_initial"] = 1

end



sdata = _SPM.build_stochastic_data_ACDC_RES(data_ots, deg, p_size)

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

sdata["curtailment"] = s


# println("\nPenetration Level = $pen_level")
# println("   Solution progress: Solving...")
result_spm_ots = _SPM.solve_sots_iv_acdc(sdata, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm_ots["primal_status"]), ")")




# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm["primal_status"])
print("Objective: ")
print(result_spm["objective"])
print("\nSolve Time: ")
print(result_spm["solve_time"])


# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm_ots["primal_status"])
print("Objective: ")
print(result_spm_ots["objective"])
print("\nSolve Time: ")
print(result_spm_ots["solve_time"])













# bus = 5

# sample1 = _SPM.sample(result_spm, "gen", bus, "pg"; sample_size=100000); 
# sample2 = _SPM.sample(result_spm, "RES", bus, "p_RES_curt"; sample_size=100000); 


# Plots.histogram(sample1)
# # Plots.histogram!(sample2)

# Plots.histogram(sample2)
# Plots.histogram!(sample1)



# # for x in sample1
# #    if x>0 && x<1e-2
# #       print("\n$x")
# #    end
# # end

# minimum(sample1)
# minimum(sample2)