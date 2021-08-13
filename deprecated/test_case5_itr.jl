# using pkgs
using Ipopt
using PolyChaos
using PowerModels
using StochasticPowerModels

# constants 
const _PMs = PowerModels
const _SPM = StochasticPowerModels

# solver
solver = Ipopt.Optimizer

# deterministic data
path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case5.m")
data = _PMs.parse_file(path)

res_dtr_ivr = _PMs.run_opf_iv(data, _PMs.IVRPowerModel, solver)

check_deterministic_solution!(data, res_dtr_ivr["solution"])

# stochastic data
deg  = 1
opq  = [Beta01OrthoPoly(deg, 1.2, 1.2; Nrec=5*deg), 
        Beta01OrthoPoly(deg, 1.2, 1.2; Nrec=5*deg), 
        GaussOrthoPoly(deg; Nrec=5*deg),
        GaussOrthoPoly(deg; Nrec=5*deg)]
mop  = MultiOrthoPoly(opq, deg)
Npce = mop.dim
Nopq = length(opq)

# build load matrix
Nd = length(data["load"])
lib = Dict(1 => 1, 2 => 2, 3 => 3)
pd, qd = zeros(Nd,Npce), zeros(Nd,Npce)

for nd in 1:Nd
    Pd, Qd = data["load"]["$nd"]["pd"], data["load"]["$nd"]["qd"]
    if nd in keys(lib)
        ns = lib[nd]
        if mop.uni[ns] isa GaussOrthoPoly
            pd[nd,[1,ns+1]] = convert2affinePCE(Pd, 0.15*Pd, mop.uni[ns])
        else
            pd[nd,[1,ns+1]] = convert2affinePCE(Pd, 0.15*Pd, mop.uni[ns], kind="μσ")
        end
    else
        pd[nd,1] = Pd
    end
    qd[nd,1] = Qd
end

# add the λs
for bus in data["bus"]
    bus[2]["λvmin"], bus[2]["λvmax"] = 1.0364, 1.0364 
    # λ = quantile(Normal(0.0,1.0),1.0-0.15)
end
for gen in data["gen"]
    gen[2]["λpmin"], gen[2]["λpmax"] = 1.0364, 1.0364
    gen[2]["λqmin"], gen[2]["λqmax"] = 1.0364, 1.0364
end
for branch in data["branch"]
    ## why no sqrt(3) in powermodel, but in the paper of line it seems to be there 29 / sqrt(3) ≈ 16 as in the paper
    branch[2]["λimax"] = 1.0364
end

# replicated data
sdata = _PMs.replicate(data, Npce)

# adding the stochastic data
for nw in 1:Npce, nd in 1:Nd
    sdata["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
    sdata["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
end

# add stochastic base data
sdata["T2"] = Tensor(2,mop)
sdata["T3"] = Tensor(3,mop)
sdata["mop"] = mop

# solve 
res_stc_ivr, solve_time = run_sopf_iv_itr(sdata, _PMs.IVRPowerModel, solver, res_dtr_ivr["solve_time"]);
