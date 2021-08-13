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

# data path
path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case3.m")
data = _PMs.parse_file(path)

res_dtr_ivr = _PMs.run_opf_iv(data, _PMs.IVRPowerModel, solver)

check_deterministic_solution!(data, res_dtr_ivr["solution"])

# stochastic data
deg  = 1
opq  = [Uniform01OrthoPoly(deg; Nrec=5*deg), 
        Uniform01OrthoPoly(deg; Nrec=5*deg),
        Uniform01OrthoPoly(deg; Nrec=5*deg)]
mop  = MultiOrthoPoly(opq, deg)
Npce = mop.dim
Nopq = length(opq)

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
    gen[2]["λqmin"], gen[2]["λqmax"] = 1.6, 2.5
end
for branch in data["branch"]
    branch[2]["λimax"] = 1.6
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
res_stc_ivr, solve_time = run_sopf_iv_itr(sdata, _PMs.IVRPowerModel, solver, res_dtr_ivr["solve_time"])