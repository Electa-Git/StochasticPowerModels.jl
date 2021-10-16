################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

@testset "Problem Formulations" begin

    @testset "IVR vs ACR - deg = 1, aux = true, case = 5-bus" begin
        # input
        deg  = 1
        aux  = true
        red  = false
        case = "matpower/case5_spm.m"

        # data
        path  = joinpath(_SPM.BASE_DIR,"test/data/$case")
        data  = _PMs.parse_file(path)
        sdata = _SPM.build_stochastic_data(data, deg)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_ivr = _SPM.run_sopf_iv(sdata, _PMs.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
        result_acr = _SPM.run_sopf_acr(sdata, _PMs.ACRPowerModel, ipopt_solver, aux=aux, deg=deg)
        
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
        case = "matpower/case5_spm.m"

        # data
        path  = joinpath(_SPM.BASE_DIR,"test/data/$case")
        data  = _PMs.parse_file(path)
        sdata = _SPM.build_stochastic_data(data, deg)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_ivr = _SPM.run_sopf_iv(sdata, _PMs.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
        result_acr = _SPM.run_sopf_acr(sdata, _PMs.ACRPowerModel, ipopt_solver, aux=aux, deg=deg)
        
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
        case = "matpower/case5_spm.m"

        # data
        path  = joinpath(_SPM.BASE_DIR,"test/data/$case")
        data  = _PMs.parse_file(path)
        sdata = _SPM.build_stochastic_data(data, deg)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_ivr = _SPM.run_sopf_iv(sdata, _PMs.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)
        result_acr = _SPM.run_sopf_acr(sdata, _PMs.ACRPowerModel, ipopt_solver, aux=aux, deg=deg)
        
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

    @testset "IVR-AUX: w vs w/o - deg = 1, case = 5-bus" begin
        # input
        deg  = 1
        red  = false
        case = "matpower/case5_spm.m"

        # data
        path  = joinpath(_SPM.BASE_DIR,"test/data/$case")
        data  = _PMs.parse_file(path)
        sdata = _SPM.build_stochastic_data(data, deg)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_w_aux  = _SPM.run_sopf_iv(sdata, _PMs.IVRPowerModel, ipopt_solver, aux=true, deg=deg, red=red)
        result_wo_aux = _SPM.run_sopf_iv(sdata, _PMs.IVRPowerModel, ipopt_solver, aux=false, deg=deg, red=red)
        
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
        case = "matpower/case5_spm.m"

        # data
        path  = joinpath(_SPM.BASE_DIR,"test/data/$case")
        data  = _PMs.parse_file(path)
        sdata = _SPM.build_stochastic_data(data, deg)

        # element cardinality
        Nb = length(data["bus"])
        Ng = length(data["gen"])
    
        # solve problem
        result_ivr = _SPM.run_sopf_iv(sdata, _PMs.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=false)
        result_red = _SPM.run_sopf_iv(sdata, _PMs.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=true)
        
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
end