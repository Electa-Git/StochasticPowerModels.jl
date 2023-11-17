################################################################################
#  Copyright 2023, Kaan Yurtseven                                              #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/Electa-Git/StochasticPowerModels.jl                    #
################################################################################

""
# function solve_sopf_iv_acdc_dim(file::String, model_constructor, optimizer; deg::Int=1, p_size=0, solution_processors=[sol_data_model!], kwargs...)
#     data = _PM.parse_file(file)
#     _PMACDC.process_additional_data!(data)
    
#     return solve_sopf_iv_acdc_dim(data, model_constructor, optimizer; deg=deg, p_size=p_size, ref_extensions = [_PMACDC.add_ref_dcgrid!, _SPM.add_ref_RES!], solution_processors=solution_processors, kwargs...)
# end

""
function solve_sopf_iv_ac(data::Dict, model_constructor, optimizer; deg::Int=1, p_size=0, solution_processors=[sol_data_model!], kwargs...)
    # @assert _IM.ismultinetwork(data) == false "The data supplied is multinetwork, it should be single-network"
    @assert model_constructor <: _PM.AbstractIVRModel "This problem type only supports the IVRModel"
    
    # sdata = build_stochastic_data_ACDC_RES(data, deg, p_size)
    result = _PM.solve_model(data, model_constructor, optimizer, build_sopf_iv_ac; multinetwork=true, ref_extensions = [_PMACDC.add_ref_dcgrid!, _SPM.add_ref_RES!], solution_processors=solution_processors, kwargs...)
    result["mop"] = _FP.dim_meta(data, :PCE_coeff, "mop")
    
    return result
end


""
function build_sopf_iv_ac(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        
        if _FP.is_first_id(pm,n,:PCE_coeff)
            bounded = true
        else 
            bounded = false
        end
        
        # if n == 1
        #     bounded = true
        # else 
        #     bounded = false
        # end

        _PM.variable_dcline_current(pm, nw=n)

        variable_bus_voltage(pm, nw=n, bounded=false)

        variable_branch_current(pm, nw=n, bounded=false)

        variable_gen_power(pm, nw=n, bounded=bounded)
        variable_gen_current(pm, nw=n, bounded=bounded)
        variable_load_current(pm, nw=n, bounded=bounded)

        variable_RES_current(pm, nw=n, bounded=bounded)

    end
    
    # objective_min_expected_generation_cost(pm)
    objective_min_expected_generation_cost_dim(pm)

    for (n, network) in _PM.nws(pm)

    # for n in sorted_nw_ids(pm)        
        
        first_stage_model_ac!(pm, n)

    end

end

function first_stage_model_ac!(pm, n)
    for i in _PM.ids(pm, :ref_buses, nw=n)
        constraint_bus_voltage_ref(pm, i, nw=n) #dim
    end

    for i in _PM.ids(pm, :bus, nw=n)
        constraint_current_balance_with_RES_ac(pm, i, nw=n) #no need dim

        constraint_gp_bus_voltage_magnitude_squared(pm, i, nw=n) #dim
    end

    for b in _PM.ids(pm, :branch, nw=n)
        _PM.constraint_voltage_drop(pm, b, nw=n)
        constraint_gp_branch_series_current_magnitude_squared(pm, b, nw=n) #dim
    end

    for g in _PM.ids(pm, :gen, nw=n)
        constraint_gp_gen_power(pm, g, nw=n) #dim
    end

    for l in _PM.ids(pm, :load, nw=n)
        constraint_gp_load_power(pm, l, nw=n) #dim
    end



    for p in _PM.ids(pm, :RES, nw=n)
        constraint_gp_RES_power(pm, p, nw=n) #dim

    end

    

    if _FP.is_first_id(pm,n,:PCE_coeff)


        for i in _PM.ids(pm, :bus, nw=n)
            constraint_cc_bus_voltage_magnitude_squared(pm, i, nw=n) #dim
        end
    
        for b in _PM.ids(pm, :branch, nw=n)
            constraint_cc_branch_series_current_magnitude_squared(pm, b, nw=n) #dim
        end
    
        for g in _PM.ids(pm, :gen, nw=n)
            constraint_cc_gen_power(pm, g, nw=n) #dim
        end
    

    end


    

end

