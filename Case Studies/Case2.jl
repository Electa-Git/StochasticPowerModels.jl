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
using XLSX

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
                                                               "max_iter"=>10000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-4,
                                                               #"fixed_variable_treatment" => "relax_bounds",
                                                )

#gPC degree input
deg  = 2

#RES input
pen_level = 0.75

#Case file and data reading
case = "case67acdc_SPM_85cc.m"
file  = joinpath(BASE_DIR, "test/data/matpower", case)
data = _PM.parse_file(file)
_PMACDC.process_additional_data!(data)

s = Dict("RES Curtailment" => true,
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



# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm["primal_status"])
print("Objective: ")
print(result_spm["objective"])
print("\nSolve Time: ")
print(result_spm["solve_time"])





for i = 1:num_of_coeff
   result_spm["solution"]["nw"]["$i"]["RES"]["4"] = Dict()
   result_spm["solution"]["nw"]["$i"]["RES"]["4"]["p_RES_curt"] = sum([result_spm["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:3])
   result_spm["solution"]["nw"]["$i"]["RES"]["4"]["p_RES"] = sum([result_spm["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES"] for j=1:3])

end

# @suppress begin
   global sample1 = Dict()
   global sample2 = Dict()
   for idx = 1:4
      sample1[idx] = _SPM.sample(result_spm, "RES", idx, "p_RES_curt"; sample_size=10000);
      sample2[idx] = _SPM.sample(result_spm, "RES", idx, "p_RES"; sample_size=10000);
   end
# end

file_name1 = "Results\\Results Curtailment $pen_level - Beta 0.85.xlsx"
fid    = XLSX.openxlsx(file_name1, mode="w")
header = ["RES 1 Curt";"RES 2 Curt";"RES 3 Curt";"RES Tot Curt";
         "RES 1";"RES 2";"RES 3";"RES Tot";]
data   = [[sample1[1]];[sample1[2]];[sample1[3]];[sample1[4]];
         [sample2[1]];[sample2[2]];[sample2[3]];[sample2[4]];]
XLSX.writetable(file_name1, data,header)




Plots.histogram(sample1[3])



# using Statistics
# mean(sample2[1])









# curt_tot = sum([result_spm["solution"]["nw"]["1"]["RES"]["$i"]["p_RES_curt"] for i=1:3])
# curt_cost = 5e6 * curt_tot
# operational_cost = result_spm["objective"] - curt_cost

# sum([result_spm["solution"]["nw"]["1"]["gen"]["$i"]["pg"] for i=1:17])
# bus = 3

# sample1 = _SPM.sample(result_spm, "RES", bus, "p_RES_curt"; sample_size=100000); 

# # sample2 = _SPM.sample(result_spm, "RES", bus, "p_RES_curt"; sample_size=100000); 


# Plots.histogram(sample1)
# # Plots.histogram!(sample2)

# Plots.histogram(sample2)
# Plots.histogram!(sample1)



# for x in sample1
#    if x>0 && x<1e-2
#       print("\n$x")
#    end
# end

# minimum(sample1)
# minimum(sample2)