################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

@testset "Formulations" begin 

    @testset "ACR vs IVR" begin
        # data
        path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case3_spm.m")
        data = _PMs.parse_file(path)
    
        # solve problem
        result_acr = run_sopf_acr(data, _PMs.ACRPowerModel, ipopt_solver, deg = 1)
        result_ivr = run_sopf_iv(data, _PMs.IVRPowerModel, ipopt_solver, deg = 1)

        sol_acr = result_acr["solution"]["nw"]
        sol_ivr = result_ivr["solution"]["nw"]
        bus_vs_acr = [[sol_acr["$nw"]["bus"]["$nb"]["vs"] for nw in 1:4] for nb in 1:3]
        bus_vs_ivr = [[sol_ivr["$nw"]["bus"]["$nb"]["vs"] for nw in 1:4] for nb in 1:3]

        @test all(isapprox.(bus_vs_acr[1], bus_vs_ivr[1], atol=1e-6))
        @test all(isapprox.(bus_vs_acr[2], bus_vs_ivr[2], atol=1e-6))
        @test all(isapprox.(bus_vs_acr[3], bus_vs_ivr[3], atol=1e-6))

        @test isapprox(result_acr["objective"], result_ivr["objective"], rtol=1e-6)

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