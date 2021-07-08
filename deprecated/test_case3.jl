# using pkgs
using Ipopt
using PolyChaos
using PowerModels
using StochasticPowerModels

# constants 
const _PMs = PowerModels
const _SPM = StochasticPowerModels

# data path
path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case3.m")

# build uncertainty data
deg  = 2
opq  = [Uniform01OrthoPoly(deg; Nrec=5*deg), 
        Uniform01OrthoPoly(deg; Nrec=5*deg),
        Uniform01OrthoPoly(deg; Nrec=5*deg)]
mop  = MultiOrthoPoly(opq, deg)
Npce = mop.dim

# parse matpower data 
data = _PMs.parse_file(path)

# build load matrix
Nd = length(data["load"])
pd, qd = zeros(Nd,Npce), zeros(Nd,Npce)
for nd in 1:Nd
    pd[nd,[1,nd+1]] = convert2affinePCE(data["load"]["$nd"]["pd"], 0.10, mop.uni[1], kind="μσ")
    qd[nd,[1,nd+1]] = convert2affinePCE(data["load"]["$nd"]["qd"], 0.05, mop.uni[2], kind="μσ")
end

# add the λs
for bus in data["bus"]
    bus[2]["λvmin"], bus[2]["λvmax"] = 1.6, 1.6
end
for gen in data["gen"]
    gen[2]["pmin"] = 0.0
    gen[2]["λpmin"], gen[2]["λpmax"] = 1.6, 1.6
    gen[2]["λqmin"], gen[2]["λqmax"] = 1.6, 1.6
end
for branch in data["branch"]
    branch[2]["imax"] = branch[2]["rate_a"]/0.9
    branch[2]["λimax"] = 1.6
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
#data["T4"] = Tensor(4,mop)
data["mop"] = mop

# solve 
solver = Ipopt.Optimizer
result = _SPM.run_sopf_acr(data, _PMs.ACRPowerModel, solver)

#result2 = run_sopf_iv(data, _PMs.IVRPowerModel, solver)

#a1=[([(result["solution"]["nw"]["$i"]["bus"]["$j"]["vs"]) for i in 1:4]) for j in 1:3]
""

#_SPM.plotHist_volt(result, "vs", mop, 1000) #cdf=true to plot CDF #PDF=

#_SPM.plotHist_gen(result, "pg", mop, 1000)



_SPM.plotHist_volt(result, "vs", mop, 1000, pdf=true) #cdf=true to plot CDF #PDF=

#_SPM.plotHist_gen(result, "qg", mop, 1000, pdf=true)
_SPM.plotHist_gen(result, "pg", mop, 1000, pdf=true)
_SPM.plotHist_branch(result, "css", mop, 1000, pdf=true)

if haskey(result["solution"]["nw"]["1"], "load")
    _SPM.plotHist_load(result, "cid", mop, 1000, pdf=true)
    _SPM.plotHist_load(result, "cid", mop, 1000)
end