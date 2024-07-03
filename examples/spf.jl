################################################################################
#  Copyright 2024, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# using pkgs 
using JuMP, Ipopt, PolyChaos
using PowerModels
using StochasticPowerModels

# constants 
const PM  = PowerModels
const PCE = PolyChaos 
const SPM = StochasticPowerModels 

# solvers
ipopt_solver = Ipopt.Optimizer

# input 
deg     = 2
case    = "case_spf.m"

# data
data    = PM.parse_file(joinpath(SPM.BASE_DIR, "test/data/matpower", case))

# build stochastic data
mom     = [ ( 0.800,1.100) # |U|₁ 
            (-10.00,10.00) # ∠U₁
            ( 0.800,1.100) # |U|₂
            (-10.00,10.00) # ∠U₂
            ( 0.000,100.0) # P₁
            (-0.484,0.484) # tanϕ₁
            ( 0.000,100.0) # P₂
            (-0.484,0.484) # tanϕ₂
            ( 0.000,100.0) # P₃
            (-0.484,0.484) # tanϕ₃
            ( 0.000,100.0) # P₄
            (-0.484,0.484) # tanϕ₄
        ]

opq     = [ PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # |U|₁
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # ∠U₁
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # |U|₂
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # ∠U₂
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # P₁
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # tanϕ₁
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # P₂
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # tanϕ₂
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # P₃
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # tanϕ₃
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # P₄
            PCE.Uniform01OrthoPoly(deg; Nrec=5*deg) # tanϕ₄
        ]

mop     = PCE.MultiOrthoPoly(opq, deg)

sdata   = PM.replicate(data, mop.dim)

sdata["T2"] = PCE.Tensor(2,mop)
sdata["T3"] = PCE.Tensor(3,mop)

sdata["mop"] = mop

# set the polynomial chaos expansion coefficients
input = zeros(length(opq), mop.dim)
for nd in 1:length(opq)
    input[nd,[1,nd+1]] = PCE.convert2affinePCE(mom[nd][1], mom[nd][2], mop.uni[nd])
end

for (nw, ntw) in sdata["nw"]
    nw = parse(Int, nw)
    # ref_bus 1
    ntw["bus"]["1"]["vr_ref"] = input[1,nw] * cosd(input[2,nw])
    ntw["bus"]["1"]["vi_ref"] = input[1,nw] * sind(input[2,nw])
    # ref_bus 2
    ntw["bus"]["2"]["vr_ref"] = input[3,nw] * cosd(input[4,nw])
    ntw["bus"]["2"]["vi_ref"] = input[3,nw] * sind(input[4,nw])
    # load 1
    ntw["load"]["1"]["pd"] = input[5,nw] / 100.0
    ntw["load"]["1"]["qd"] = input[5,nw] * input[6,nw] / 100.0
    # load 2
    ntw["load"]["2"]["pd"] = input[7,nw] / 100.0
    ntw["load"]["2"]["qd"] = input[7,nw] * input[8,nw] / 100.0
    # load 3
    ntw["load"]["3"]["pd"] = input[9,nw] / 100.0
    ntw["load"]["3"]["qd"] = input[9,nw] * input[10,nw] / 100.0
    # load 4
    ntw["load"]["4"]["pd"] = input[11,nw] / 100.0
    ntw["load"]["4"]["qd"] = input[11,nw] * input[12,nw] / 100.0
end

result = solve_spf_iv(sdata, PM.IVRPowerModel, ipopt_solver)