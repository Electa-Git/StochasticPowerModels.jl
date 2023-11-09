################################################################################
#  Copyright 2023, Kaan Yurtseven                                              #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/Electa-Git/StochasticPowerModels.jl                    #
################################################################################

""
function solve_sopf_iv_acdc(file::String, model_constructor, optimizer; deg::Int=1, p_size=0, solution_processors=[sol_data_model!], kwargs...)
    data = _PM.parse_file(file)
    _PMACDC.process_additional_data!(data)
    
    return solve_sopf_iv_acdc(data, model_constructor, optimizer; deg=deg, p_size=p_size, ref_extensions = [_PMACDC.add_ref_dcgrid!, _SPM.add_ref_RES!], solution_processors=solution_processors, kwargs...)
end

""
function solve_sopf_iv_acdc(data::Dict, model_constructor, optimizer; deg::Int=1, p_size=0, solution_processors=[sol_data_model!], kwargs...)
    @assert _IM.ismultinetwork(data) == false "The data supplied is multinetwork, it should be single-network"
    @assert model_constructor <: _PM.AbstractIVRModel "This problem type only supports the IVRModel"
    
    sdata = build_stochastic_data_ACDC_RES(data, deg, p_size)
    result = _PM.solve_model(sdata, model_constructor, optimizer, build_sopf_iv_acdc; multinetwork=true, ref_extensions = [_PMACDC.add_ref_dcgrid!, _SPM.add_ref_RES!], solution_processors=solution_processors, kwargs...)
    result["mop"] = sdata["mop"]
    
    return result
end


""
function build_sopf_iv_acdc(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        if n == 1
            bounded = true
        else 
            bounded = false
        end

        _PM.variable_dcline_current(pm, nw=n)

        variable_bus_voltage(pm, nw=n, bounded=false)

        variable_branch_current(pm, nw=n, bounded=false)

        variable_gen_power(pm, nw=n, bounded=bounded)
        variable_gen_current(pm, nw=n, bounded=bounded)
        variable_load_current(pm, nw=n, bounded=bounded)

        #DC grid variables
        variable_active_dcbranch_flow(pm, nw=n, bounded=bounded)
        variable_dcbranch_current(pm, nw=n, bounded=bounded)
        variable_dcgrid_voltage_magnitude(pm, nw=n, bounded=bounded)
        #DC converter variables
        variable_dc_converter(pm, nw=n, bounded=bounded)
        variable_dc_converter_squared(pm, nw=n, bounded=bounded)

        variable_RES_current(pm, nw=n, bounded=bounded)

    end
    
    objective_min_expected_generation_cost(pm)


    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n) 
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance_with_RES(pm, i, nw=n)

            constraint_gp_bus_voltage_magnitude_squared(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_voltage_drop(pm, b, nw=n)
            constraint_gp_branch_series_current_magnitude_squared(pm, b, nw=n)
        end

        for g in _PM.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end


        #DC Constraints        
        for i in _PM.ids(pm, :busdc, nw=n)
            constraint_current_balance_dc(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :branchdc, nw=n)
            constraint_gp_ohms_dc_branch(pm, i, nw=n) 
            constraint_ohms_dc_branch(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :convdc, nw=n)
            # constraint_gp_filter_voltage_squared(pm, i, nw=n)
            # constraint_gp_converter_voltage_squared(pm, i, nw=n)
            # constraint_gp_transformer_current_from_squared(pm, i, nw=n)
            # constraint_gp_transformer_current_to_squared(pm, i, nw=n)


            # constraint_gp_reactor_current_from_squared(pm, i, nw=n)
            # constraint_gp_reactor_current_to_squared(pm, i, nw=n)

            constraint_gp_iconv_lin_squared_1(pm, i, nw=n)
            constraint_gp_iconv_lin_squared_2(pm, i, nw=n)

            constraint_gp_converter_dc_power(pm, i, nw=n)
            constraint_gp_converter_losses(pm, i, nw=n)
            constraint_gp_converter_ac_power(pm, i, nw=n)

            constraint_conv_transformer(pm, i, nw=n)
            constraint_conv_reactor(pm, i, nw=n)
            constraint_conv_filter(pm, i, nw=n)

        end

        for p in _PM.ids(pm, :RES, nw=n)
            constraint_gp_RES_power(pm, p, nw=n)

        end

    end

    for i in _PM.ids(pm, :bus, nw=1)
        constraint_cc_bus_voltage_magnitude_squared(pm, i, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_cc_branch_series_current_magnitude_squared(pm, b, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_cc_gen_power(pm, g, nw=1)
    end

    for i in _PM.ids(pm, :busdc, nw=1)
        constraint_cc_conv_voltage_magnitude(pm, i, nw=1)
    end

    for i in _PM.ids(pm, :branchdc, nw=1)
        constraint_cc_dc_branch_current(pm, i, nw=1)
    end

    for i in _PM.ids(pm, :convdc, nw=1)

        # constraint_cc_iconv_lin_squared(pm, i, nw=1)
        constraint_cc_iconv_lin(pm, i, nw=1)

        constraint_cc_conv_ac_power(pm, i, nw=1)
        constraint_cc_conv_dc_power(pm, i, nw=1)
        constraint_cc_converter_dc_current(pm, i, nw=1)
        
        
        # constraint_cc_filter_voltage_squared(pm, i, nw=1)
        # constraint_cc_converter_voltage_squared(pm, i, nw=1)
        # constraint_cc_transformer_current_from_squared(pm, i, nw=1)
        # constraint_cc_transformer_current_to_squared(pm, i, nw=1)


        # constraint_cc_reactor_current_from_squared(pm, i, nw=1)
        # constraint_cc_reactor_current_to_squared(pm, i, nw=1)

    end

end