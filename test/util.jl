
@testset "Utilities" begin
    deg  = 1
    aux  = true
    red  = false
    case = "case5_spm.m"

    # data
    path  = joinpath("../test/data/matpower",case)
    data  = _PM.parse_file(path)
    sdata = _SPM.build_stochastic_data(data, deg)

    # element cardinality
    Nb = length(data["bus"])
    Ng = length(data["gen"])

    # solve problem
    result = _SPM.run_sopf_iv(sdata, _PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red)

    @test result["termination_status"] == LOCALLY_SOLVED


    pg_coeff = pce_coeff(result, "gen", 1, "pg") 
    @test isapprox(pg_coeff, [1.5000000149998045;    0.000459749295155349;   0.00010187050082378998])

    vr_coeff = pce_coeff(result, "bus", 2, "vr") 
    @test isapprox(vr_coeff, [  0.9902405059792688;  -0.0015349621281835933;  -0.0002292009456770499  ])


    pg_sample = sample(sdata, result, "gen", 1, "pg"; sample_size=10) 
    @test length(pg_sample) == 10

    pg_density = density(sdata, result, "gen", 1, "pg"; sample_size=10) 
    # @test TODO: what to test here? do we re-export KernelDensity.UnivariateKDE?

    _SPM.print_summary(result["solution"])
end

