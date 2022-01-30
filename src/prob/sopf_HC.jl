################################################################################
#  Copyright 2021, Tom Van Acker, Frederik Geth                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# NOTE: dc lines are omitted from the current formulation                      #
################################################################################

""
function run_sopf_hc(data::Dict, model_constructor, optimizer; aux::Bool=true, deg::Int=1, red::Bool=false, solution_processors=[sol_data_model!], kwargs...)
    @assert _IM.ismultinetwork(data) == false "The data supplied is multinetwork, it should be single-network"
    @assert model_constructor <: _PM.AbstractIVRModel "This problem type only supports the IVRModel"
    sdata = build_stochastic_data_hc(data, deg)
    if aux && !red #only with aux and without reduced variable
        result = _PM.run_model(sdata, model_constructor, optimizer, build_sopf_hc; multinetwork=true, solution_processors=solution_processors, kwargs...)
    end
    result["mop"] = sdata["mop"]
    return result
end

""
function run_sopf_hc(file::String, model_constructor, optimizer; aux::Bool=true, deg::Int=1, red::Bool=false, solution_processors=[sol_data_model!], kwargs...)
    data = _PM.parse_file(file)
    return run_sopf_hc(data, model_constructor, optimizer; aux=aux, deg=deg, red=red, solution_processors=solution_processors, kwargs...)
end

""
function build_sopf_hc(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        variable_bus_voltage(pm, nw=n, aux=true)

        variable_branch_current(pm, nw=n, aux=true)

        variable_gen_power(pm, nw=n, bounded=false)                             # enforcing bounds alters the objective 
        variable_gen_current(pm, nw=n, bounded=false)                           # enforcing bounds makes problem infeasible
        variable_load_current(pm, nw=n)
        variable_PV_current(pm, nw=n)
        if n==1
        variable_PV_size(pm, nw=n, bounded=true)
        end
    end

    for i in _PM.ids(pm, :bus, nw=1)
        constraint_bus_voltage_squared_cc_limit(pm, i, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_gen_power_cc_limit(pm, g, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_branch_series_current_squared_cc_limit(pm, b, nw=1)
    end

    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_current_balance_with_PV(pm, i, nw=n)
            constraint_gp_bus_voltage_squared(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_current_from(pm, b, nw=n)
            _PM.constraint_current_to(pm, b, nw=n)

            _PM.constraint_voltage_drop(pm, b, nw=n)

            constraint_gp_branch_series_current_squared(pm, b, nw=n)
        end

        for g in _PM.ids(pm, :gen, nw=n)
            constraint_gp_gen_power(pm, g, nw=n)
        end

        for l in _PM.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end

        for p in _PM.ids(pm, :PV, nw=n)
            constraint_gp_pv_power(pm, p, nw=n)
        end
    end

    objective_max_PV(pm)
end