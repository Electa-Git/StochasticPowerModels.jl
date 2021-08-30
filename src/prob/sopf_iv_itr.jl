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
function run_sopf_iv_itr(data, model_constructor, optimizer; deg::Int=1, max_iter::Int=10, time_limit::Float64=3600.0, kwargs...)
    Memento.info(_LOGGER, "maximum iterations set to value of $max_iter")

    # build the stochastic data
    sdata = build_stochastic_data(data, deg)

    # initialize
    iter = 1
    start_time = time()
    global violated = false
    cnstr_gen, cnstr_bus, cnstr_branch = String[], String[], String[]
    
    # instantiate and solve the deterministic opf
    dtr_pm = _PMs.instantiate_model(data, model_constructor, _PMs.build_opf_iv)
    dtr_result = _PMs.optimize_model!(dtr_pm, optimizer=optimizer)

    # instantiate the unconstrained stochastic opf
    stc_pm = _PMs.instantiate_model(sdata, model_constructor, build_sopf_iv_unc)

    # add the initial bounds based on the deterministic solution
    for (ng, gen) in data["gen"] if !in(ng, cnstr_gen)
        gen_id = parse(Int, ng)
        gen_res = dtr_result["solution"]["gen"][ng]

        pg, qg = gen_res["pg"], gen_res["qg"]
        pmin, pmax = gen["pmin"], gen["pmax"]
        qmin, qmax = gen["qmin"], gen["qmax"]

        if isapprox(pmin, pg, rtol=1e-6) || isapprox(pg, pmax, rtol=1e-6) ||
           isapprox(qmin, qg, rtol=1e-6) || isapprox(qg, qmax, rtol=1e-6)
            
            violated = true
            push!(cnstr_gen, ng)

            constraint_gen_power_cc_limit(stc_pm, gen_id, nw=1)
        end
    end end
    for (nb, bus) in data["bus"] if !in(nb, cnstr_bus)
        bus_id = parse(Int, nb)
        bus_res = dtr_result["solution"]["bus"][nb]

        vr, vi = bus_res["vr"], bus_res["vi"]
        vmin, vmax = bus["vmin"], bus["vmax"]

        vm = sqrt(vr^2 + vi^2)

        vs  = [_PMs.var(stc_pm, nw, :vs, bus_id) for nw in _PMs.nw_ids(stc_pm)]

        if isapprox(vmin, vm, rtol=1e-6) || isapprox(vm, vmax, rtol=1e-6) 
            violated = true
            push!(cnstr_bus, nb)

            JuMP.unfix.(vs)
            for (n, network) in _PMs.nws(stc_pm)
                constraint_gp_bus_voltage_squared(stc_pm, bus_id, nw=n)
            end
            constraint_bus_voltage_squared_cc_limit(stc_pm, bus_id, nw=1)
        end
    end end
    for (nb, branch) in data["branch"] if !in(nb, cnstr_branch)
        branch_id = parse(Int, nb)
        branch_res = dtr_result["solution"]["branch"][nb]

        csr, csi = branch_res["csr_fr"], branch_res["csi_fr"]
        cmax = branch["cmax"]

        csm = sqrt(csr^2 + csi^2)

        css  = Dict(nw => _PMs.var(stc_pm, nw, :css, branch_id) for nw in _PMs.nw_ids(stc_pm))

        if isapprox(csm, cmax, rtol=1e-6) 
            violated = true
            push!(cnstr_branch, nb)

            JuMP.unfix.(css)
            for (n, network) in _PMs.nws(stc_pm)
                constraint_gp_branch_series_current_squared(stc_pm, branch_id, nw=n)
            end
            constraint_branch_series_current_squared_cc_limit(stc_pm, branch_id, nw=1)
        end
    end end

    while violated && iter <= max_iter && (time() - start_time) < time_limit
        violated = false

        # solve the stochastic opf
        @time global stc_result = _PMs.optimize_model!(stc_pm, optimizer=optimizer)

        # add the necessary bounds based on the last stochastic solution
        for (ng, gen) in data["gen"] if !in(ng, cnstr_gen)
            gen_id = parse(Int, ng)

            pg = [stc_result["solution"]["nw"]["$nw"]["gen"][ng]["pg"] for nw in sorted_nw_ids(stc_pm)]
            qg = [stc_result["solution"]["nw"]["$nw"]["gen"][ng]["qg"] for nw in sorted_nw_ids(stc_pm)]
            
            pmin, pmax = gen["pmin"], gen["pmax"]
            qmin, qmax = gen["qmin"], gen["qmax"]
    
            λpmin, λpmax = gen["λpmin"], gen["λpmax"]
            λqmin, λqmax = gen["λqmin"], gen["λqmax"]
    
            if pmin > var_min(pg, λpmin) || var_max(pg, λpmax) > pmax || 
               qmin > var_min(qg, λqmin) || var_max(qg, λqmax) > qmax
                
                violated = true
                push!(cnstr_gen, ng)
    
                constraint_gen_power_cc_limit(stc_pm, gen_id, nw=1)
            end
        end end
        for (nb, bus) in data["bus"] if !in(nb, cnstr_bus)
            bus_id = parse(Int, nb)

            vr = [stc_result["solution"]["nw"]["$nw"]["bus"][nb]["vr"] for nw in sorted_nw_ids(stc_pm)]
            vi = [stc_result["solution"]["nw"]["$nw"]["bus"][nb]["vi"] for nw in sorted_nw_ids(stc_pm)]
    
            vmin, vmax = bus["vmin"], bus["vmax"]
    
            λvmin, λvmax = bus["λvmin"], bus["λvmax"]
    
            vm = sqrt.(vr.^2 + vi.^2)
    
            vs  = [_PMs.var(stc_pm, nw, :vs, bus_id) for nw in _PMs.nw_ids(stc_pm)]
    
            if vmin > var_min(vm, λvmin) || var_max(vm, λvmax) > vmax
                violated = true
                push!(cnstr_bus, nb)
    
                JuMP.unfix.(vs)
                for (n, network) in _PMs.nws(stc_pm)
                    constraint_gp_bus_voltage_squared(stc_pm, bus_id, nw=n)
                end
                constraint_bus_voltage_squared_cc_limit(stc_pm, bus_id, nw=1)
            end
        end end
        for (nb, branch) in data["branch"] if !in(nb, cnstr_branch)
            branch_id = parse(Int, nb)

            csr = [stc_result["solution"]["nw"]["$nw"]["branch"][nb]["csr_fr"] for nw in sorted_nw_ids(stc_pm)]
            csi = [stc_result["solution"]["nw"]["$nw"]["branch"][nb]["csi_fr"] for nw in sorted_nw_ids(stc_pm)]
            
            cmax = branch["cmax"]
    
            λcmax = branch["λcmax"]
    
            csm = sqrt.(csr.^2 + csi.^2)

            css = [_PMs.var(stc_pm, nw, :css, branch_id) for nw in _PMs.nw_ids(stc_pm)]
    
            if var_max(csm, λcmax) > cmax 
                violated = true
                push!(cnstr_branch, nb)
    
                JuMP.unfix.(css)
                for (n, network) in _PMs.nws(stc_pm)
                    constraint_gp_branch_series_current_squared(stc_pm, branch_id, nw=n)
                end
                constraint_branch_series_current_squared_cc_limit(stc_pm, branch_id, nw=1)
            end
        end end

        iter += 1
    end

    stc_result["cnstr"] = Dict("gen" => cnstr_gen, 
                               "bus" => cnstr_bus, 
                               "branch" => cnstr_branch)
    stc_result["solve_time"] = time() - start_time
    stc_result["iterations"] = iter

    return dtr_result, stc_result
end

""
function build_sopf_iv_unc(pm::AbstractPowerModel)
    for (n, network) in _PMs.nws(pm) 
        variable_bus_voltage(pm, nw=n, bounded=false, aux_fix=true)
        variable_branch_current(pm, nw=n, aux_fix=true)

        variable_gen_power(pm, nw=n, bounded=false)
        variable_gen_current(pm, nw=n, bounded=false)
        variable_load_current(pm, nw=n, bounded=false)
    end

    for (n, network) in _PMs.nws(pm)
        for i in _PMs.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PMs.ids(pm, :bus, nw=n)
            constraint_current_balance(pm, i, nw=n)
        end

        for b in _PMs.ids(pm, :branch, nw=n)
            _PMs.constraint_current_from(pm, b, nw=n)
            _PMs.constraint_current_to(pm, b, nw=n)

            _PMs.constraint_voltage_drop(pm, b, nw=n)
        end

        for g in _PMs.ids(pm, :gen, nw=n) 
            constraint_gp_gen_power(pm, g, nw=n) 
        end

        for l in _PMs.ids(pm, :load, nw=n)
            constraint_gp_load_power(pm, l, nw=n)
        end
    end

    objective_min_expected_fuel_cost(pm)
end

