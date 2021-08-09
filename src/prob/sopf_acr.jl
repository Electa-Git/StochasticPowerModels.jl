################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# NOTE: dc lines are omitted from the current formulation                      #
################################################################################

""
function run_sopf_acr(file, model_constructor::Type, optimizer; kwargs...)
    return _PMs.run_model(file, model_constructor, optimizer, build_sopf_acr; multinetwork=true, kwargs...)
end

""
function build_sopf_acr(pm::AbstractPowerModel)
    for (n, network) in _PMs.nws(pm) 
        variable_bus_voltage(pm, nw=n, bounded=false)
        variable_gen_power(pm, nw=n, bounded=false)
        variable_branch_power(pm, nw=n, bounded=false)
        variable_branch_current(pm, nw=n, bounded=false)   
        _PMs.variable_dcline_power(pm, nw=n, bounded=false)   
    end

    for i in _PMs.ids(pm, :bus,nw=1)
        constraint_bus_voltage_squared_cc_limit(pm, i,nw=1)
    end

    for g in _PMs.ids(pm, :gen,nw=1)
        #if g==1 || g==2 
        constraint_gen_power_cc_limit(pm, g,nw=1)
        #end
    end

    for b in _PMs.ids(pm, :branch, nw=1)
       constraint_branch_series_current_squared_cc_limit(pm, b,nw=1)
    end

    for (n, network) in _PMs.nws(pm)

        for i in _PMs.ids(pm, :ref_buses, nw=n)
           constraint_theta_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_power_balance(pm, i, nw=n)
            constraint_gp_bus_voltage_squared(pm, i, nw=n)
            #constraint_voltage_magnitude_bounds(pm, i, nw=n)
            #constraint_voltage_setpoint(pm, i, nw=n)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            constraint_gp_power_branch_to(pm, b, nw=n)
            constraint_gp_power_branch_from(pm, b, nw=n)
            
            constraint_branch_voltage(pm, b, nw=n)
            constraint_gp_current_squared(pm, b, nw=n) 
           
            
                
            #following are simplified only g and based in Tillemans paper; but apparantly not working
                #constraint_gp_power_branch_to_simplified(pm, b, nw=n)
                #constraint_gp_power_branch_from_simplified(pm, b, nw=n)
            
        end



        #for d in _PMs.ids(pm, :dcline, nw=n)
        #   _PMs.constraint_dcline_power_losses(pm, d, nw=n)
        #end

    end
    
   
    # for d in _PMs.ids(pm, :dcline)                                                 # needs to be implemented, similar to constraint_branch_series_current_squared_cc_limit
    #     constraint_dcline_current_squared_cc_limit(pm, d)
    # end

    objective_min_expected_fuel_cost(pm) 
    #objective_min_fuel_cost_poly(pm)                                      # needs to be implemented, based on final polynomial.
end