################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

@testset "Problem Formulations" begin

    @testset "IVR vs iterative IVR"
        # data
        path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case14_spm.m")
        data = _PMs.parse_file(path)
    
        # solve problem
        result_stc = run_sopf_iv(data, _PMs.IVRPowerModel, solver, deg = 1)
    
        # solve problem iteratively
        result_dtr, result_itr = run_sopf_iv_itr(data, _PMs.IVRPowerModel, solver, deg = 1);
    
        # assert
        @test isapprox(result_stc["objective"], result_itr["objective"], rtol=1e-6)
    end

end