using PowerModels



# Load the necessary Pkgs
using JuMP
using LinearAlgebra
using MosekTools
using PolyChaos

# Parameters
Nn  = 2
Nb  = 1
Ng  = 1
Nd  = 1

A   = [-1   1]
Bb  = diagm(0 => [-0.5])
Ψ   = [ zeros(Nb) -Bb*A[:,2:end]*inv(A[:,2:end]'*Bb*A[:,2:end]) ]

Cg  = [ 1
        0]
Cd  = [ 0
        1]
cg  = 4.0
λg, λb = 1.6*ones(Ng), 1.6*ones(Nb)                                             # Chance constraint lambda'sets
pgmin, pgmax = zeros(Ng), 10*ones(Ng)                                           # Generator limits
pbmin, pbmax = -10*ones(Nb), 10*ones(Nb) 

# PolyChaos definition of the stochastic variables
deg = 1
opq = [Uniform01OrthoPoly(deg, Nrec=5*deg)]
mop = MultiOrthoPoly(opq, deg)
K   = mop.dim

pd   = zeros(Nd,K)
pd[1,[1,2]] = convert2affinePCE(1.0, 0.1, mop.uni[1], kind="μσ")

# Second-order cone definition
function buildSOC(x::Vector, mop::MultiOrthoPoly)
    t = [sqrt(Tensor(2,mop).get([i,i])) for i in 0:mop.dim-1 ]
    (t.*x)[2:end]
end

# JuMP model
"""
pdᵢ ~ U(μᵢ,σᵢ)                                      ∀ i ∈ Nd
pb = Ψ(Cg*pg + Cd*pd)                   

min.    ∑ᵢ cgᵢ ⋅ 𝔼(pgᵢ)
s.t.    ∑ᵢ pgᵢ - ∑ᵢ pdᵢ = 0
        p̲gᵢ ≤ 𝔼(pgᵢ) ± λgᵢ ⋅ √(𝕍(pgᵢ)) ≤ p̄gᵢ,       ∀ i ∈ Ng
        p̲bᵢ ≤ 𝔼(pbᵢ) ± λbᵢ ⋅ √(𝕍(pbᵢ)) ≤ p̄bᵢ,       ∀ i ∈ Nb
"""
m = Model(with_optimizer(Mosek.Optimizer))

@variable(m, pg[i in 1:Ng,j in 1:K], base_name = "pg")                          # Generator power variable ∀ generators, PCE nodes
                                                                                # Power Balance
@constraint(m, [j in 1:K], sum(pg[i,j] for i in 1:Ng) - sum(pd[i,j] for i in 1:Nd) == 0.0)
                                                                                # Generator power limits
@constraint(m, [i in 1:Ng], [1/λg[i]*(pgmax[i] - mean(pg[i,:],mop)); buildSOC(pg[i,:],mop)] in SecondOrderCone())
@constraint(m, [i in 1:Ng], [1/λg[i]*(mean(pg[i,:],mop) - pgmin[i]); buildSOC(pg[i,:],mop)] in SecondOrderCone())
                                                                                # Branch power limits
pb = Ψ*(Cg*pg + Cd*pd)
@constraint(m, [i in 1:Nb], [1/λb[i]*(pbmax[i] - mean(pb[1,:],mop)); buildSOC(pb[i,:],mop)] in SecondOrderCone())
@constraint(m, [i in 1:Nb], [1/λb[i]*(mean(pb[1,:],mop) - pbmin[i]); buildSOC(pb[i,:],mop)] in SecondOrderCone())
                                                                                # Minimize generator costs
@objective(m, Min, sum(mean(pg[i,:],mop)*cg[i] for i in 1:Ng))

# Solve model
@time optimize!(m)

Pg1 = value.(pg)
Pb1 = Ψ*(Cg*Pg1 + Cd*pd)

Ki = [[1],[]]
Ko = [[],[1]]
Kd = [[],[1]]
Kg = [[1],[]]

B = [(1,1,2)]

# JuMP model
"""
pdᵢ ~ U(μᵢ,σᵢ)                                      ∀ i ∈ Nd         

min.    ∑ᵢ cgᵢ ⋅ 𝔼(pgᵢ)
s.t.    ∑ᵢ pgᵢ - ∑ᵢ pdᵢ + ∑ᵢ pb = 0                      ∀ i ∈ N
        p̲gᵢ ≤ 𝔼(pgᵢ) ± λgᵢ ⋅ √(𝕍(pgᵢ)) ≤ p̄gᵢ,       ∀ i ∈ Ng
        p̲bᵢ ≤ 𝔼(pbᵢ) ± λbᵢ ⋅ √(𝕍(pbᵢ)) ≤ p̄bᵢ,       ∀ i ∈ Nb
        
"""
n = Model(with_optimizer(Mosek.Optimizer))

@variable(n, an[i in 1:Nn, j in 1:K], base_name = "an")                         # Nodal voltage angle ∀ nodes, PCE nodes
@variable(n, pg[i in 1:Ng, j in 1:K], base_name = "pg")                         # Generator power variable ∀ generators, PCE nodes
@variable(n, pb[i in 1:Nb, j in 1:K], base_name = "pb")                         # Branch power variables ∀ branches, PCE nodes    
                                                                                # Power Balance
@constraint(n, [i in 1:Nn, j in 1:K], sum(pg[k,j] for k in Kg[i]) + sum(pb[k,j] for k in Ki[i]) ==
                                      sum(pd[k,j] for k in Kd[i]) + sum(pb[k,j] for k in Ko[i]))
                                                                                # Branch power flow
@constraint(n, [x in B, j in 1:K], pb[x[1],j] == Bb[x[1],x[1]]*(an[x[3],j]-an[x[2],j]))
                                                                                # Generator power limits
@constraint(n, [i in 1:Ng, j in 1:K], [1/λg[i]*(pgmax[i] - mean(pg[i,k],mop)); buildSOC(pg[i,k],mop)] in SecondOrderCone())
@constraint(n, [i in 1:Ng], [1/λg[i]*(mean(pg[i,:],mop) - pgmin[i]); buildSOC(pg[i,:],mop)] in SecondOrderCone())
                                                                                # Branch power limits
@constraint(n, [i in 1:Nb], [1/λb[i]*(pbmax[i] - mean(pb[1,:],mop)); buildSOC(pb[i,:],mop)] in SecondOrderCone())
@constraint(n, [i in 1:Nb], [1/λb[i]*(mean(pb[1,:],mop) - pbmin[i]); buildSOC(pb[i,:],mop)] in SecondOrderCone())
                                                                                # Minimize generator costs
@objective(n, Min, sum(mean(pg[i,:],mop)*cg[i] for i in 1:Ng))

# Solve model
@time optimize!(n)

Pg2 = value.(pg)
Pb2 = value.(pb)

@assert isapprox(Pg1,Pg2, atol = 1e-6)
@assert isapprox(Pb1,Pb2, atol = 1e-6)