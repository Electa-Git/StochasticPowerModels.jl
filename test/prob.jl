################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

@testset "Problem Formulations" begin

    @testset "Convenience functions for IVR stochastic OPF" begin
        deg  = 1
        aux  = true
        red  = true
        case = "case5_spm.m"

        file  = joinpath(_SPM.BASE_DIR, "test/data/matpower", case)
        result_ivr = _SPM.run_sopf_iv(file, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
        @test result_ivr["termination_status"] == LOCALLY_SOLVED
        obj1 = result_ivr["objective"]

        data  = _PM.parse_file(file)
        result_ivr2 = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver; aux=aux, deg=deg, red=red)
        @test result_ivr2["termination_status"] == LOCALLY_SOLVED
        obj2 = result_ivr2["objective"]

        sdata = _SPM.build_stochastic_data(data, deg)
        result_ivr3 = _PM.run_model(sdata, _PM.IVRPowerModel, ipopt_solver, _SPM.build_sopf_iv_reduced_with_aux; multinetwork=true, solution_processors=[_PM.sol_data_model!])
        @test result_ivr3["termination_status"] == LOCALLY_SOLVED
        obj3 = result_ivr3["objective"]

        @test isapprox(obj1, obj2)
        @test isapprox(obj1, obj3)
    end

    @testset "Convenience functions for ACR stochastic OPF" begin
        deg  = 1
        aux  = true
        case = "case5_spm.m"

        file  = joinpath("../test/data/matpower", case)
        result_acr = _SPM.run_sopf_acr(file, _PM.ACRPowerModel, ipopt_solver, aux=aux, deg=deg)
        @test result_acr["termination_status"] == LOCALLY_SOLVED
        obj1 = result_acr["objective"]

        data  = _PM.parse_file(file)
        result_acr2 = _SPM.run_sopf_acr(data, _PM.ACRPowerModel, ipopt_solver; aux=aux, deg=deg)
        @test result_acr2["termination_status"] == LOCALLY_SOLVED
        obj2 = result_acr2["objective"]

        sdata = _SPM.build_stochastic_data(data, deg)
        result_acr3 = _PM.run_model(sdata, _PM.ACRPowerModel, ipopt_solver, _SPM.build_sopf_acr_with_aux; multinetwork=true, solution_processors=[_PM.sol_data_model!])
        @test result_acr3["termination_status"] == LOCALLY_SOLVED
        obj3 = result_acr3["objective"]

        @test isapprox(obj1, obj2)
        @test isapprox(obj1, obj3)
    end


    @testset "IVR vs ACR - deg = 1, aux = true, case = 5-bus" begin
        # input
        deg  = 1
        aux  = true
        red  = false
        case = "case5_spm.m"

        # data
        path  = joinpath("../test/data/matpower", case)
        data  = _PM.parse_file(path)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_ivr = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
        result_acr = _SPM.run_sopf_acr(data, _PM.ACRPowerModel, ipopt_solver, aux=aux, deg=deg)

        @test result_ivr["termination_status"] == LOCALLY_SOLVED
        @test result_acr["termination_status"] == LOCALLY_SOLVED

        # test for objective
        @test isapprox(result_ivr["objective"], result_acr["objective"], rtol=1e-8)

        # test for real bus voltage
        bus_vr_ivr = [_SPM.pce_coeff(result_ivr, "bus", nb, "vr") for nb in 1:Nb]
        bus_vr_acr = [_SPM.pce_coeff(result_acr, "bus", nb, "vr") for nb in 1:Nb]
        for nb in 1:Nb 
            @test all(isapprox.(bus_vr_ivr[nb], bus_vr_acr[nb], atol=1e-6)) 
        end

        # test for imaginary bus voltage
        bus_vi_ivr = [_SPM.pce_coeff(result_ivr, "bus", nb, "vi") for nb in 1:Nb]
        bus_vi_acr = [_SPM.pce_coeff(result_acr, "bus", nb, "vi") for nb in 1:Nb]
        for nb in 1:Nb
            @test all(isapprox.(bus_vi_ivr[nb], bus_vi_acr[nb], atol=1e-6)) 
        end

        # test for active gen power
        bus_pg_ivr = [_SPM.pce_coeff(result_ivr, "gen", ng, "pg") for ng in 1:Ng]
        bus_pg_acr = [_SPM.pce_coeff(result_acr, "gen", ng, "pg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_pg_ivr[ng], bus_pg_acr[ng], atol=1e-6)) 
        end

        # test for reactive gen power
        bus_qg_ivr = [_SPM.pce_coeff(result_ivr, "gen", ng, "qg") for ng in 1:Ng]
        bus_qg_acr = [_SPM.pce_coeff(result_acr, "gen", ng, "qg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_qg_ivr[ng], bus_qg_acr[ng], atol=1e-6)) 
        end
    end

    @testset "IVR vs ACR - deg = 2, aux = true, case = 5-bus" begin
        # input
        deg  = 2
        aux  = true
        red  = false
        case = "case5_spm.m"

        # data
        path  = joinpath("../test/data/matpower", case)
        data  = _PM.parse_file(path)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_ivr = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
        result_acr = _SPM.run_sopf_acr(data, _PM.ACRPowerModel, ipopt_solver, aux=aux, deg=deg)

        @test result_ivr["termination_status"] == LOCALLY_SOLVED
        @test result_acr["termination_status"] == LOCALLY_SOLVED
        
        # test for objective
        @test isapprox(result_ivr["objective"], result_acr["objective"], rtol=1e-8)

        # test for real bus voltage
        bus_vr_ivr = [_SPM.pce_coeff(result_ivr, "bus", nb, "vr") for nb in 1:Nb]
        bus_vr_acr = [_SPM.pce_coeff(result_acr, "bus", nb, "vr") for nb in 1:Nb]
        for nb in 1:Nb 
            @test all(isapprox.(bus_vr_ivr[nb], bus_vr_acr[nb], atol=1e-6)) 
        end

        # test for imaginary bus voltage
        bus_vi_ivr = [_SPM.pce_coeff(result_ivr, "bus", nb, "vi") for nb in 1:Nb]
        bus_vi_acr = [_SPM.pce_coeff(result_acr, "bus", nb, "vi") for nb in 1:Nb]
        for nb in 1:Nb
            @test all(isapprox.(bus_vi_ivr[nb], bus_vi_acr[nb], atol=1e-6)) 
        end

        # test for active gen power
        bus_pg_ivr = [_SPM.pce_coeff(result_ivr, "gen", ng, "pg") for ng in 1:Ng]
        bus_pg_acr = [_SPM.pce_coeff(result_acr, "gen", ng, "pg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_pg_ivr[ng], bus_pg_acr[ng], atol=1e-6)) 
        end

        # test for reactive gen power
        bus_qg_ivr = [_SPM.pce_coeff(result_ivr, "gen", ng, "qg") for ng in 1:Ng]
        bus_qg_acr = [_SPM.pce_coeff(result_acr, "gen", ng, "qg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_qg_ivr[ng], bus_qg_acr[ng], atol=1e-6)) 
        end
    end

    @testset "IVR vs ACR - deg = 1, aux = false, case = 5-bus" begin
        # input
        deg  = 1
        aux  = false
        red  = false
        case = "case5_spm.m"

        # data
        path  = joinpath("../test/data/matpower", case)
        data  = _PM.parse_file(path)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_ivr = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
        result_acr = _SPM.run_sopf_acr(data, _PM.ACRPowerModel, ipopt_solver, aux=aux, deg=deg)

        @test result_ivr["termination_status"] == LOCALLY_SOLVED
        @test result_acr["termination_status"] == LOCALLY_SOLVED
        
        # test for objective
        @test isapprox(result_ivr["objective"], result_acr["objective"], rtol=1e-8)

        # test for real bus voltage
        bus_vr_ivr = [_SPM.pce_coeff(result_ivr, "bus", nb, "vr") for nb in 1:Nb]
        bus_vr_acr = [_SPM.pce_coeff(result_acr, "bus", nb, "vr") for nb in 1:Nb]
        for nb in 1:Nb 
            @test all(isapprox.(bus_vr_ivr[nb], bus_vr_acr[nb], atol=1e-6)) 
        end

        # test for imaginary bus voltage
        bus_vi_ivr = [_SPM.pce_coeff(result_ivr, "bus", nb, "vi") for nb in 1:Nb]
        bus_vi_acr = [_SPM.pce_coeff(result_acr, "bus", nb, "vi") for nb in 1:Nb]
        for nb in 1:Nb
            @test all(isapprox.(bus_vi_ivr[nb], bus_vi_acr[nb], atol=1e-6)) 
        end

        # test for active gen power
        bus_pg_ivr = [_SPM.pce_coeff(result_ivr, "gen", ng, "pg") for ng in 1:Ng]
        bus_pg_acr = [_SPM.pce_coeff(result_acr, "gen", ng, "pg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_pg_ivr[ng], bus_pg_acr[ng], atol=1e-6)) 
        end

        # test for reactive gen power
        bus_qg_ivr = [_SPM.pce_coeff(result_ivr, "gen", ng, "qg") for ng in 1:Ng]
        bus_qg_acr = [_SPM.pce_coeff(result_acr, "gen", ng, "qg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_qg_ivr[ng], bus_qg_acr[ng], atol=1e-6)) 
        end
    end

    @testset "IVR: w vs w/o aux - deg = 1, case = 5-bus" begin
        # input
        deg  = 1
        red  = false
        case = "case5_spm.m"

        # data
        path  = joinpath("../test/data/matpower", case)
        data  = _PM.parse_file(path)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_w_aux  = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=true, deg=deg, red=red)
        result_wo_aux = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=false, deg=deg, red=red)

        @test result_w_aux["termination_status"] == LOCALLY_SOLVED
        @test result_wo_aux["termination_status"] == LOCALLY_SOLVED

        
        # test for objective
        @test isapprox(result_w_aux["objective"], result_wo_aux["objective"], rtol=1e-8)

        # test for real bus voltage
        bus_vr_w_aux  = [_SPM.pce_coeff(result_w_aux, "bus", nb, "vr") for nb in 1:Nb]
        bus_vr_wo_aux = [_SPM.pce_coeff(result_wo_aux, "bus", nb, "vr") for nb in 1:Nb]
        for nb in 1:Nb 
            @test all(isapprox.(bus_vr_w_aux[nb], bus_vr_wo_aux[nb], atol=1e-6)) 
        end

        # test for imaginary bus voltage
        bus_vi_w_aux  = [_SPM.pce_coeff(result_w_aux, "bus", nb, "vi") for nb in 1:Nb]
        bus_vi_wo_aux = [_SPM.pce_coeff(result_wo_aux, "bus", nb, "vi") for nb in 1:Nb]
        for nb in 1:Nb
            @test all(isapprox.(bus_vi_w_aux[nb], bus_vi_wo_aux[nb], atol=1e-6)) 
        end

        # test for active gen power
        bus_pg_w_aux  = [_SPM.pce_coeff(result_w_aux, "gen", ng, "pg") for ng in 1:Ng]
        bus_pg_wo_aux = [_SPM.pce_coeff(result_wo_aux, "gen", ng, "pg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_pg_w_aux[ng], bus_pg_wo_aux[ng], atol=1e-6)) 
        end

        # test for reactive gen power
        bus_qg_w_aux  = [_SPM.pce_coeff(result_w_aux, "gen", ng, "qg") for ng in 1:Ng]
        bus_qg_wo_aux = [_SPM.pce_coeff(result_wo_aux, "gen", ng, "qg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_qg_w_aux[ng], bus_qg_wo_aux[ng], atol=1e-6)) 
        end
    end

    @testset "IVR-RED: true vs false - deg = 1, aux = true, case = 5-bus" begin
        # input
        deg  = 1
        aux  = true
        case = "case5_spm.m"

        # data
        path  = joinpath("../test/data/matpower", case)
        data  = _PM.parse_file(path)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_ivr = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=false)
        result_red = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=true)

        @test result_ivr["termination_status"] == LOCALLY_SOLVED
        @test result_red["termination_status"] == LOCALLY_SOLVED
        
        # test for objective
        @test isapprox(result_ivr["objective"], result_red["objective"], rtol=1e-8)

        # test for real bus voltage
        bus_vr_ivr = [_SPM.pce_coeff(result_ivr, "bus", nb, "vr") for nb in 1:Nb]
        bus_vr_red = [_SPM.pce_coeff(result_red, "bus", nb, "vr") for nb in 1:Nb]
        for nb in 1:Nb 
            @test all(isapprox.(bus_vr_ivr[nb], bus_vr_red[nb], atol=1e-6)) 
        end

        # test for imaginary bus voltage
        bus_vi_ivr = [_SPM.pce_coeff(result_ivr, "bus", nb, "vi") for nb in 1:Nb]
        bus_vi_red = [_SPM.pce_coeff(result_red, "bus", nb, "vi") for nb in 1:Nb]
        for nb in 1:Nb
            @test all(isapprox.(bus_vi_ivr[nb], bus_vi_red[nb], atol=1e-6)) 
        end

        # test for active gen power
        bus_pg_ivr = [_SPM.pce_coeff(result_ivr, "gen", ng, "pg") for ng in 1:Ng]
        bus_pg_red = [_SPM.pce_coeff(result_red, "gen", ng, "pg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_pg_ivr[ng], bus_pg_red[ng], atol=1e-6)) 
        end

        # test for reactive gen power
        bus_qg_ivr = [_SPM.pce_coeff(result_ivr, "gen", ng, "qg") for ng in 1:Ng]
        bus_qg_red = [_SPM.pce_coeff(result_red, "gen", ng, "qg") for ng in 1:Ng]
        for ng in 1:Ng
            @test all(isapprox.(bus_qg_ivr[ng], bus_qg_red[ng], atol=1e-6)) 
        end
    end

    @testset "IVR: case30_spm_muhlpfordt.m, σ = 0.15 * μ, ε = 0.15" begin
        # input
        aux  = true
        red  = true
        case = "case30_spm_muhlpfordt.m"

        # data
        path  = joinpath("../test/data/matpower", case)
        data  = _PM.parse_file(path)
    
        # solve problem
        result_one = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=1, red=red)
        result_two = _SPM.run_sopf_iv(data, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=2, red=red)

        @test result_one["termination_status"] == LOCALLY_SOLVED
        @test result_two["termination_status"] == LOCALLY_SOLVED
        
        # test for objective
        @test isapprox(result_one["objective"], 599.35, atol=1e-2)
        @test isapprox(result_two["objective"], 599.35, atol=1e-2)

        # sample against chance constraint
        sample_size = 10000
        for ns in [("gen", 3, "pg", "pmax"), ("gen", 4, "pg", "pmax")]
            elm, id, var, lim = ns 
            limit = data[elm]["$id"][lim]
            ε_one = count(x -> x <= limit, sample(result_one, elm, id, var, sample_size=sample_size)) / sample_size
            ε_two = count(x -> x <= limit, sample(result_two, elm, id, var, sample_size=sample_size)) / sample_size
            @test isapprox(ε_one, 1.0 - 0.15, atol = 1.25e-2)
            @test isapprox(ε_two, 1.0 - 0.15, atol = 1.25e-2)
        end
    end
end