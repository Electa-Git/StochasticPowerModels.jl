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
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>2000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-6,
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

#gPC degree input
deg  = 2

#RES input
pen_level = 0.4

#Case file and data reading
case = "case5_ACDC_SPM_95cc_cost.m"
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


# println("\nPenetration Level = $pen_level")
# println("   Solution progress: Solving...")
# result_spm_variance = _SPM.solve_sopf_iv_acdc_variance(sdata, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
# println("   Solution progress: Solved! (", string(result_spm_variance["primal_status"]), ")")

println("\nPenetration Level = $pen_level")
println("   Solution progress: Solving...")
result_spm_VaR = _SPM.solve_sopf_iv_acdc_VaR(sdata, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
println("   Solution progress: Solved! (", string(result_spm_VaR["primal_status"]), ")")








for i = 1:num_of_coeff
   result_spm["solution"]["nw"]["$i"]["gen"]["6"] = Dict()
   result_spm["solution"]["nw"]["$i"]["gen"]["6"]["pg"] = sum([result_spm["solution"]["nw"]["$i"]["gen"]["$j"]["pg"] for j=1:5])

   result_spm_variance["solution"]["nw"]["$i"]["gen"]["6"] = Dict()
   result_spm_variance["solution"]["nw"]["$i"]["gen"]["6"]["pg"] = sum([result_spm_variance["solution"]["nw"]["$i"]["gen"]["$j"]["pg"] for j=1:5])

   result_spm_VaR["solution"]["nw"]["$i"]["gen"]["6"] = Dict()
   result_spm_VaR["solution"]["nw"]["$i"]["gen"]["6"]["pg"] = sum([result_spm_VaR["solution"]["nw"]["$i"]["gen"]["$j"]["pg"] for j=1:5])


   # result_spm["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
   # result_spm["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_spm["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])

   # result_spm_variance["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
   # result_spm_variance["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_spm_variance["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])

   # result_spm_VaR["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
   # result_spm_VaR["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_spm_VaR["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])
end






# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm["primal_status"])
print("Objective: ")
print(result_spm["objective"])
print("\nSolve Time: ")
print(result_spm["solve_time"])

# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm_variance["primal_status"])
print("Objective: ")
print(result_spm_variance["objective"])
print("\nSolve Time: ")
print(result_spm_variance["solve_time"])

# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_spm_VaR["primal_status"])
print("Objective: ")
print(result_spm_VaR["objective"])
print("\nSolve Time: ")
print(result_spm_VaR["solve_time"])

bus = 6

sample1 = _SPM.sample(result_spm, "gen", 6, "pg"; sample_size=10000);    
sample2 = _SPM.sample(result_spm_variance, "gen", 6, "pg"; sample_size=10000); 
sample3 = _SPM.sample(result_spm_VaR, "gen", 6, "pg"; sample_size=10000); 

using Statistics

Statistics.mean(sample1)
Statistics.var(sample1)

Statistics.mean(sample2)
Statistics.var(sample2)

Statistics.mean(sample3)
Statistics.var(sample3)

Plots.histogram(sample1, bin = 100)
# Plots.histogram!(sample2, bin = 100)
Plots.histogram!(sample3, bin = 100)



# sample1 = _SPM.sample(result_spm, "gen", 6, "pg"; sample_size=10000);    
# sample2 = _SPM.sample(result_spm_variance, "gen", 6, "pg"; sample_size=10000); 
# sample3 = _SPM.sample(result_spm_VaR, "gen", 6, "pg"; sample_size=10000); 

# sample4 = _SPM.sample(result_spm, "RES", 6, "p_RES_curt"; sample_size=10000); 
# sample5 = _SPM.sample(result_spm_variance, "RES", 6, "p_RES_curt"; sample_size=10000); 
# sample6 = _SPM.sample(result_spm_VaR, "RES", 6, "p_RES_curt"; sample_size=10000); 


# using XLSX

# file_name1 = "Results\\Risk\\Results - Risk SOPF - beta = 0.25 - alpha = 1.285 - $case ($pen_level).xlsx"
# fid    = XLSX.openxlsx(file_name1, mode="w")
# header = ["Pg"; "Pg_variance"; "Pg_VaR"; "RES_curt"; "RES_curt_variance"; "RES_curt_VaR"]
# data   = [[sample1]; [sample2]; [sample3]; [sample4]; [sample5]; [sample6]]
# XLSX.writetable(file_name1, data,header)






# sample1 = _SPM.sample(result_SOPF[pen_level_start], "branch", 2, "cr_fr"; sample_size=10000);    
# sample2 = _SPM.sample(result_SOPF[pen_level_start], "branch", 3, "cr_fr"; sample_size=10000); 
# sample3 = _SPM.sample(result_SOPF[pen_level_start], "branch", 14, "cr_fr"; sample_size=10000); 
# sample4 = _SPM.sample(result_SOPF[pen_level_start], "branch", 15, "cr_fr"; sample_size=10000); 
# sample5 = _SPM.sample(result_SOPF[pen_level_start], "branch", 16, "cr_fr"; sample_size=10000); 
# sample6 = _SPM.sample(result_SOPF[pen_level_start], "branch", 17, "cr_fr"; sample_size=10000); 

# sample7 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 1, "pf"; sample_size=10000); 
# sample8 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 2, "pf"; sample_size=10000); 
# sample9 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 3, "pf"; sample_size=10000);    
# sample10 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 4, "pf"; sample_size=10000); 
# sample11 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 6, "pf"; sample_size=10000); 
# sample12 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 7, "pf"; sample_size=10000); 
# sample13 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 8, "pf"; sample_size=10000); 
# sample14 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 9, "pf"; sample_size=10000); 
# sample15 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 11, "pf"; sample_size=10000); 
# sample16 = _SPM.sample(result_SOPF[pen_level_start], "branchdc", 12, "pf"; sample_size=10000); 


# file_name1 = "Results\\Results - SOPF (2e-6) - (atol = 1e-6) - $case ($pen_level_start).xlsx"
# fid    = XLSX.openxlsx(file_name1, mode="w")
# header = ["1"; "2"; "3"; "4"; "5"; "6"; "7"; "8"; "9"; "10";
#             "11"; "12"; "13"; "14"; "15";
#             "16"]
# data   = [[sample1]; [sample2]; [sample3]; [sample4]; [sample5];
#             [sample6]; [sample7]; [sample8]; [sample9]; [sample10];
#             [sample11]; [sample12]; [sample13]; [sample14]; [sample15];
#             [sample16]]
# XLSX.writetable(file_name1, data,header)














# Plots.histogram(sample2)
# Plots.histogram!(sample1)



# # for x in sample1
# #    if x>0 && x<1e-2
# #       print("\n$x")
# #    end
# # end

# minimum(sample1)
# minimum(sample2)