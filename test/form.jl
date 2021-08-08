################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

@testset "Formulations" begin 

    @testset "ACR vs IVR" begin
        path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case3.m")
        data = _PMs.parse_file(path)
        
        # set-up stochastic input
        deg  = 1
        opq  = [Uniform01OrthoPoly(deg; Nrec=5*deg), 
                Uniform01OrthoPoly(deg; Nrec=5*deg),
                Uniform01OrthoPoly(deg; Nrec=5*deg)]
        mop  = MultiOrthoPoly(opq, deg)
        Npce = mop.dim
        Nd = length(data["load"])
        pd, qd = zeros(Nd,Npce), zeros(Nd,Npce)
        for nd in 1:Nd
            pd[nd,[1,nd+1]] = convert2affinePCE(data["load"]["$nd"]["pd"], 0.10, mop.uni[1], kind="μσ")
            qd[nd,[1,nd+1]] = convert2affinePCE(data["load"]["$nd"]["qd"], 0.05, mop.uni[2], kind="μσ")
        end
        for bus in data["bus"]
            bus[2]["λvmin"], bus[2]["λvmax"] = 1.6, 1.6
        end
        for gen in data["gen"]
            gen[2]["pmin"] = 0.0
            gen[2]["λpmin"], gen[2]["λpmax"] = 1.6, 1.6
            gen[2]["λqmin"], gen[2]["λqmax"] = 1.6, 2.5
        end
        for branch in data["branch"]
            branch[2]["imax"] = branch[2]["rate_a"]/0.9
            branch[2]["λimax"] = 1.6
        end

        # replicate data
        data = _PMs.replicate(data, Npce)

        # add stochastic input to data
        for nw in 1:Npce, nd in 1:Nd
            data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
            data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
        end
        data["T2"] = Tensor(2,mop)
        data["T3"] = Tensor(3,mop)
        data["mop"] = mop
        
        res_acr = run_sopf_acr(data, _PMs.ACRPowerModel, ipopt_solver)
        res_ivr = run_sopf_iv(data, _PMs.IVRPowerModel, ipopt_solver)

        @test isapprox(res_acr["objective"], res_ivr["objective"], rtol=1e-6)
    end

end