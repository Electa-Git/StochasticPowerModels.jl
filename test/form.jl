################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

@testset "Formulations" begin 

    @testset "ACR vs IVR (deg=1)" begin
        # Example from: Chance-constrained ac optimal power flow: A polynomial 
        # chaos approach (p. 8)

        # data
        path = joinpath(_SPM.BASE_DIR, "test/data/matpower/case30_spm_muhlpfordt.m")
        data = _PMs.parse_file(path)

        # solve problem
        result_ivr = run_sopf_iv(data, _PMs.IVRPowerModel, ipopt_solver, deg=1)
        result_acr = run_sopf_acr(data, _PMs.ACRPowerModel, ipopt_solver, deg=1)

        # test
        @test isapprox(result_ivr["objective"], 599.35, rtol=1e-4)
        @test isapprox(result_acr["objective"], 599.35, rtol=1e-4)
    end

    @testset "ACR vs IVR (deg=2)" begin
        # Example from: Chance-constrained ac optimal power flow: A polynomial 
        # chaos approach (p. 8)

        # data
        path = joinpath(_SPM.BASE_DIR, "test/data/matpower/case30_spm_muhlpfordt.m")
        data = _PMs.parse_file(path)

        # solve problem
        result_ivr = run_sopf_iv(data, _PMs.IVRPowerModel, ipopt_solver, deg=2)
        result_acr = run_sopf_acr(data, _PMs.ACRPowerModel, ipopt_solver, deg=2)

        # test
        @test isapprox(result_ivr["objective"], 599.35, rtol=1e-4)
        @test isapprox(result_acr["objective"], 599.35, rtol=1e-4)
    end

    @testset "IVR vs reduced IVR" begin
        # data
        path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case3_spm.m")
        data = _PMs.parse_file(path)
    
        # solve problem
        result_ivr = run_sopf_iv(data, _PMs.IVRPowerModel, ipopt_solver, deg = 1)
        result_red = run_sopf_iv_reduced(data, _PMs.IVRPowerModel, ipopt_solver, deg = 1)

        sol_ivr = result_ivr["solution"]["nw"]
        sol_red = result_red["solution"]["nw"]
        bus_vs_ivr = [[sol_ivr["$nw"]["bus"]["$nb"]["vs"] for nw in 1:4] for nb in 1:3]
        bus_vs_red = [[sol_red["$nw"]["bus"]["$nb"]["vs"] for nw in 1:4] for nb in 1:3]

        @test all(isapprox.(bus_vs_ivr[1], bus_vs_red[1], atol=1e-6))
        @test all(isapprox.(bus_vs_ivr[2], bus_vs_red[2], atol=1e-6))
        @test all(isapprox.(bus_vs_ivr[3], bus_vs_red[3], atol=1e-6))

        @test isapprox(result_ivr["objective"], result_red["objective"], rtol=1e-6)

    end
    
end