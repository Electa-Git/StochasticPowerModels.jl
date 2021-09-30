function run_case_118(path_to_result::String)
    # case
    case = "case 118"

    # dataframe
    df = _DFs.DataFrame(case=String[],
                        form=String[], 
                        prob=String[],
                        deg=Int64[], 
                        aux=String[],
                        limit=Float64[],
                        std=Float64[], 
                        solve_time=Float64[],
                        termination_status=String[], 
                        objective=Float64[])

    # data
    path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case118_spm.m")
    data = _PMs.parse_file(path)

    # solver
    solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, 
                                            "print_level" => 0, 
                                            "max_cpu_time" => 3600.0,
                                            "linear_solver" => "ma27")

    # test loop
    for deg in [1,2], aux in [true,false], ε in [0.05, 0.10, 0.15], σ in [0.10, 0.15] if !(deg == 2 && !aux)
        # aux string
        aux_bool = ifelse(aux, "true", "false")

        # adjust the λ
        λ = _DST.quantile(_DST.Normal(0.0,1.0), 1.0 - ε)
        for (b, bus) in data["bus"]
            bus["λvmin"] = λ
            bus["λvmax"] = λ
        end 
        for (g, gen) in data["gen"]
            gen["λpmin"] = λ
            gen["λpmax"] = λ
            gen["λqmin"] = λ
            gen["λqmax"] = λ
        end
        for (br, branch) in data["branch"]
            branch["λcmax"] = λ
        end

        # add load mean and variance
        for (l, load) in data["load"]
            bus_id = load["load_bus"]
            data["bus"]["$(bus_id)"]["μ"] = load["pd"]
            data["bus"]["$(bus_id)"]["σ"] = load["pd"] * σ
        end

        # solve problem
        println("acr full")
        result_acr_stc = _SPM.run_sopf_acr(data, _PMs.ACRPowerModel, solver, aux=aux, deg=deg)
        push!(df, [case, "acr", "full", deg, aux_bool, ε, σ, result_acr_stc["solve_time"], string(result_acr_stc["termination_status"]), result_acr_stc["objective"]])

        # solve ivr problem
        println("ivr full")
        result_ivr_stc = run_sopf_iv(data, _PMs.IVRPowerModel, solver, aux=aux, deg=deg)
        push!(df, [case, "ivr", "full", deg, aux_bool, ε, σ, result_ivr_stc["solve_time"], string(result_ivr_stc["termination_status"]), result_ivr_stc["objective"]])

        # solve reduced ivr problem
        println("ivr reduced")
        result_ivr_red = run_sopf_iv_reduced(data, _PMs.IVRPowerModel, solver, aux=aux, deg=deg)
        push!(df, [case, "ivr", "reduced", deg, aux_bool, ε, σ, result_ivr_red["solve_time"], string(result_ivr_red["termination_status"]), result_ivr_red["objective"]])

    end end

    # write away results
    CSV.write(path_to_result, df)

end