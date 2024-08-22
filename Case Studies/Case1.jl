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
using Distributions
using PolyChaos

#Constants
const _PM = PowerModels
const _SPM = StochasticPowerModels
const _PMACDC = PowerModelsACDC
const _FP = FlexPlan

using Suppressor

Memento.setlevel!(Memento.getlogger(StochasticPowerModels), "error")
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModels), "error")
Memento.setlevel!(Memento.getlogger(PowerModelsACDC), "error")
Memento.setlevel!(Memento.getlogger(Distributions), "error")
Memento.setlevel!(Memento.getlogger(PolyChaos), "error")

#Solver inputs
ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, 
                                                               "max_iter"=>2000, 
                                                               "sb"=>"yes", 
                                                               "tol" => 1e-4,
                                                               #"dependency_detector" => "mumps"
                                                               # "fixed_variable_treatment" => "relax_bounds",
                                                )

#gPC degree input
deg  = 2

#RES inputs
pen_level_start = 1.30
pen_level_step = 0.05
pen_level_end = 1.30

#Report decision
Report_obj = false
Report_curt = false

#Necessary initializations
result_case1 = Dict()
result_case2 = Dict()
obj_case1 = Dict()
obj_case2 = Dict()
stat_case1 = Dict()
stat_case2 = Dict()
time_case1 = Dict()
time_case2 = Dict()
sample1 = Dict()
sample2 = Dict()
p_size_dict = Dict()

global feas_ctr1 = 0
global feas_ctr2 = 0

global solve_case1 = false
global solve_case2 = true

#Case file and data reading
case1 = "case5_AC_SPM_95cc.m"
file1  = joinpath(BASE_DIR, "test/data/matpower/", case1)
data1 = _PM.parse_file(file1)
_PMACDC.process_additional_data!(data1)

case2 = "case5_ACDC_SPM_95cc.m"
file2  = joinpath(BASE_DIR, "test/data/matpower/", case2)
data2 = _PM.parse_file(file2)
_PMACDC.process_additional_data!(data2)

s = Dict("RES Curtailment" => true,
         "Load Curtailment" => false,
         )


#Penetration level iteration
for pen_level = pen_level_start:pen_level_step:pen_level_end

   if solve_case1 || solve_case2
      println("Penetration Level = $pen_level")

      global total_load = sum([data2["load"]["$i"]["pd"] for i=1:length(data2["load"])])
      global p_size = pen_level * total_load / length(data2["RES"]) #Calculate RES size for each RES bus
      p_size_dict[pen_level] = p_size

      if solve_case1

         sdata1 = _SPM.build_stochastic_data_ACDC_RES(data1, deg, p_size)

         PCE_data1 = Dict(
               "T2" => sdata1["T2"],
               "T3" => sdata1["T3"],
               "T4" => sdata1["T4"],
               "mop" => sdata1["mop"],
            )
         local num_of_coeff = PCE_data1["mop"].dim

         delete!(sdata1, "T2")
         delete!(sdata1, "T3")
         delete!(sdata1, "T4")
         delete!(sdata1, "mop")

         _FP.add_dimension!(sdata1, :PCE_coeff, num_of_coeff; metadata = PCE_data1)

         sdata1["curtailment"] = s

         println("   Case 1 - Solution progress: Solving...")
         global result_case1[pen_level] = _SPM.solve_sopf_iv_acdc(sdata1, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);

         #Store necessary values for reporting
         obj_case1[pen_level] = result_case1[pen_level]["objective"]
         stat_case1[pen_level] = string(result_case1[pen_level]["primal_status"])
         time_case1[pen_level] = result_case1[pen_level]["solve_time"]
         
         if sdata1["curtailment"]["RES Curtailment"] == true
            for i = 1:num_of_coeff
               result_case1[pen_level]["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
               result_case1[pen_level]["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_case1[pen_level]["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])
            end

            @suppress begin
               global sample1[pen_level] = Dict()
               for idx = 1:6
                  sample1[pen_level][idx] = _SPM.sample(result_case1[pen_level], "RES", idx, "p_RES_curt"; sample_size=10000);
               end
            end
         end
         if Report_curt
            file_name1 = "Results\\Results Curtailment $pen_level - $case1.xlsx"
            fid    = XLSX.openxlsx(file_name1, mode="w")
            header = ["RES 1";"RES 2";"RES 3";"RES 4";"RES 5";"RES Tot";]
            data   = [[sample1[pen_level][1]];[sample1[pen_level][2]];[sample1[pen_level][3]];[sample1[pen_level][4]];[sample1[pen_level][5]];[sample1[pen_level][6]];]
            XLSX.writetable(file_name1, data,header)
         end

         if string(result_case1[pen_level]["primal_status"]) != "FEASIBLE_POINT" && string(result_case1[pen_level]["primal_status"]) !="NEARLY_FEASIBLE_POINT"
            global feas_ctr1 += 1
         else
            global feas_ctr1 = 0
         end

         println("   Case 1 - Solution progress: Solved! (", string(result_case1[pen_level]["primal_status"]), ")")

         if feas_ctr1 >= 2
            global solve_case1 = false
            println("   Case 1 reached infeasibility on penetration level of $pen_level.")
         end

      end

      if solve_case2

         sdata2 = _SPM.build_stochastic_data_ACDC_RES(data2, deg, p_size)

         PCE_data2 = Dict(
               "T2" => sdata2["T2"],
               "T3" => sdata2["T3"],
               "T4" => sdata2["T4"],
               "mop" => sdata2["mop"],
            )
         local num_of_coeff = PCE_data2["mop"].dim

         delete!(sdata2, "T2")
         delete!(sdata2, "T3")
         delete!(sdata2, "T4")
         delete!(sdata2, "mop")

         _FP.add_dimension!(sdata2, :PCE_coeff, num_of_coeff; metadata = PCE_data2)

         sdata2["curtailment"] = s

         println("   Case 2 - Solution progress: Solving...")
         global result_case2[pen_level] = _SPM.solve_sopf_iv_acdc(sdata2, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);

         #Store necessary values for reporting
         obj_case2[pen_level] = result_case2[pen_level]["objective"]
         stat_case2[pen_level] = string(result_case2[pen_level]["primal_status"])
         time_case2[pen_level] = result_case2[pen_level]["solve_time"]
         
         if sdata2["curtailment"]["RES Curtailment"] == true
            for i = 1:num_of_coeff
               result_case2[pen_level]["solution"]["nw"]["$i"]["RES"]["6"] = Dict()
               result_case2[pen_level]["solution"]["nw"]["$i"]["RES"]["6"]["p_RES_curt"] = sum([result_case2[pen_level]["solution"]["nw"]["$i"]["RES"]["$j"]["p_RES_curt"] for j=1:5])
            end
            @suppress begin
               global sample2[pen_level] = Dict()
               for idx = 1:6
                  sample2[pen_level][idx] = _SPM.sample(result_case2[pen_level], "RES", idx, "p_RES_curt"; sample_size=10000);
               end
            end
         end
         if Report_curt
            file_name2 = "Results\\Results Curtailment $pen_level - $case2.xlsx"
            fid    = XLSX.openxlsx(file_name2, mode="w")
            header = ["RES 1";"RES 2";"RES 3";"RES 4";"RES 5";"RES Tot";]
            data   = [[sample2[pen_level][1]];[sample2[pen_level][2]];[sample2[pen_level][3]];[sample2[pen_level][4]];[sample2[pen_level][5]];[sample2[pen_level][6]];]
            XLSX.writetable(file_name2, data,header)
         end

         if string(result_case2[pen_level]["primal_status"]) != "FEASIBLE_POINT" && string(result_case2[pen_level]["primal_status"]) !="NEARLY_FEASIBLE_POINT"
            global feas_ctr2 += 1
         else
            global feas_ctr2 = 0
         end

         println("   Case 2 - Solution progress: Solved! (", string(result_case2[pen_level]["primal_status"]), ")")

         if feas_ctr2 >= 2
            global solve_case2 = false
            println("   Case 2 reached infeasibility on penetration level of $pen_level.")
         end

      end

   end

end


# Reporting
if Report_obj
   if solve_case1
      file_name1 = "Results\\Results - $case1 ($pen_level_start - $pen_level_end).xlsx"
      fid    = XLSX.openxlsx(file_name1, mode="w")
      header = ["Penetration Level";"Objective Value";"Status";"Time"]
      data   = [[collect(keys(obj_case1))];[collect(values(obj_case1))];[collect(values(stat_case1))];[collect(values(time_case1))]]
      XLSX.writetable(file_name1, data,header)
   end
   if solve_case2
      file_name2 = "Results\\Results - $case2 ($pen_level_start - $pen_level_end).xlsx"
      fid    = XLSX.openxlsx(file_name2, mode="w")
      header = ["Penetration Level";"Objective Value";"Status";"Time"]
      data   = [[collect(keys(obj_case2))];[collect(values(obj_case2))];[collect(values(stat_case2))];[collect(values(time_case2))]]
      XLSX.writetable(file_name2, data,header)
   end
end




# #Case file and data reading
# case = "case5_AC_SPM_95cc.m"
# file  = joinpath(BASE_DIR, "test/data/matpower", case)
# data = _PM.parse_file(file)
# _PMACDC.process_additional_data!(data)

# s = Dict("RES Curtailment" => true,
#          "Load Curtailment" => false,
#          )


# total_load = sum([data["load"]["$i"]["pd"] for i=1:length(data["load"])])
# p_size = pen_level * total_load / length(data["RES"]) #Calculate RES size for each RES bus

# sdata = _SPM.build_stochastic_data_ACDC_RES(data, deg, p_size)

# PCE_data = Dict(
#       "T2" => sdata["T2"],
#       "T3" => sdata["T3"],
#       "T4" => sdata["T4"],
#       "mop" => sdata["mop"],
#    )
# num_of_coeff = PCE_data["mop"].dim

# delete!(sdata, "T2")
# delete!(sdata, "T3")
# delete!(sdata, "T4")
# delete!(sdata, "mop")

# _FP.add_dimension!(sdata, :PCE_coeff, num_of_coeff; metadata = PCE_data)

# sdata["curtailment"] = s


# println("\nPenetration Level = $pen_level")
# println("   Solution progress: Solving...")
# result_spm = _SPM.solve_sopf_iv_acdc(sdata, _PM.IVRPowerModel, ipopt_solver, deg=deg, p_size=p_size);
# println("   Solution progress: Solved! (", string(result_spm["primal_status"]), ")")



# Show results on the terminal
println("\n\n>>> SPMACDC Results >>>")
println(result_case2[pen_level_start]["primal_status"])
print("Objective: ")
print(result_case2[pen_level_start]["objective"])
print("\nSolve Time: ")
print(result_case2[pen_level_start]["solve_time"])

bus = 5

sample1 = _SPM.sample(result_case2[pen_level_start], "RES", bus, "p_RES_curt"; sample_size=100000); 

# sample2 = _SPM.sample(result_spm, "RES", bus, "p_RES_curt"; sample_size=100000); 


Plots.histogram(sample1)
# using Statistics
# mean(sample1)
# # # Plots.histogram!(sample2)

# # Plots.histogram(sample2)
# # Plots.histogram!(sample1)



# # for x in sample1
# #    if x>0 && x<1e-2
# #       print("\n$x")
# #    end
# # end

# # minimum(sample1)
# # minimum(sample2)