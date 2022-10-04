
@testset "Utilities" begin
    deg  = 1
    case = "case5_spm.m"

    file  = joinpath("../test/data/matpower", case)

    # solve problem
    result = _SPM.solve_sopf_iv(file, _PM.IVRPowerModel, ipopt_solver, deg=deg)
    @test result["termination_status"] == LOCALLY_SOLVED

    pg_coeff = pce_coeff(result, "gen", 1, "pg") 
    @test isapprox(pg_coeff, [1.5000000149998045;    0.000459749295155349;   0.00010187050082378998])

    vr_coeff = pce_coeff(result, "bus", 2, "vr") 
    @test isapprox(vr_coeff, [0.9902405059792688;  -0.0015349621281835933;  -0.0002292009456770499  ])


    pg_sample = sample(result, "gen", 1, "pg"; sample_size=10) 
    @test length(pg_sample) == 10

    pg_density = density(result, "gen", 1, "pg"; sample_size=10) 
    # @test TODO: what to test here? do we re-export KernelDensity.UnivariateKDE?

end

