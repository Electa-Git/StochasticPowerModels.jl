## TODO
## 1) enable load matrix builder for degrees higher than 1.
## 2) construct script to run deterministic opf, only activate violated constraints.

# using pkgs
using Ipopt
using PolyChaos
using PowerModels
using StochasticPowerModels

# constants 
const _PMs = PowerModels
const _SPM = StochasticPowerModels

# data path
path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case30.m")

# build uncertainty data
deg  = 1
opq  = [Beta01OrthoPoly(deg, 2.0, 2.0; Nrec=5*deg), 
        Beta01OrthoPoly(deg, 2.0, 5.0; Nrec=5*deg), 
        GaussOrthoPoly(deg; Nrec=5*deg),
        GaussOrthoPoly(deg; Nrec=5*deg)]
mop  = MultiOrthoPoly(opq, deg)
Npce = mop.dim
Nopq = length(opq)

# parse matpower data 
data = _PMs.parse_file(path)

# reset the line rates
data["branch"]["30"]["rate_a"] *= 11/16
data["branch"]["35"]["rate_a"] *= 12/16

# build load matrix
Nd = length(data["load"])
lib = Dict(2 => 1, 3 => 1, 4 => 2, 10 => 4, 21 => 4, 24 => 3)
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
    gen[2]["pmin"] = 0.0
    gen[2]["λpmin"], gen[2]["λpmax"] = 1.0364, 1.0364
    gen[2]["λqmin"], gen[2]["λqmax"] = 1.0364, 1.0364
end
for branch in data["branch"]
    f_bus = branch[2]["f_bus"]
    branch[2]["imax"] = branch[2]["rate_a"]/data["bus"]["$f_bus"]["vmin"] 
    ## why no sqrt(3) in powermodel, but in the paper of line it seems to be there 29 / sqrt(3) ≈ 16 as in the paper
    branch[2]["λimax"] = 1.0364
end

# replicated data
data = _PMs.replicate(data, Npce)

# adding the stochastic data
for nw in 1:Npce, nd in 1:Nd
    data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
    data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
end

# add stochastic base data
data["T2"] = Tensor(2,mop)
data["T3"] = Tensor(3,mop)
data["mop"] = mop

# solve 
solver = Ipopt.Optimizer
result = _SPM.run_sopf_iv(data, _PMs.IVRPowerModel, solver)


_SPM.plotHist_volt(result, "vs", mop, 1000, busno=[6,7,8,9] ,pdf=true) #cdf=true to plot CDF #PDF=

_SPM.plotHist_gen(result, "qg", mop, 1000, pdf=true)
_SPM.plotHist_gen(result, "pg", mop, 1000, pdf=true)
_SPM.plotHist_branch(result, "css", mop, 1000, pdf=true)
