# An uncertainty file needs to be created, similar to what has been done for
# PMS. However, in contrast to that file 

unc_id,dst,param
1,"Normal", [5.0,0.5]
2,"Uniform",[5.0,7.5]

# This file need to be parsed into an uncertainty array
unc = [ Normal(5.0,0.5),
        Uniform(5.0,7.5)]

# The maximal degree
deg = 2

# Create the polynomials
gpc = [ GaussOrthoPoly(deg) for n in 1:60]


# Create the affine coefficients
pce = [ convert2affinePCE(unc[1].μ, unc[1].σ, gpc[1]),
        convert2affinePCE(unc[2].a, unc[2].b, gpc[2])]


mop = MultiOrthoPoly(gpc, deg)

# Feed the mop through the power flow algorithm
# xxx

# Only remember the necessary data?
# ?

# Resample
@time ξ = sampleMeasure(10000, mop);
@time samples = evaluatePCE(rand(1891),ξ,mop);


deg = 4;
ns  = 10000;
op = GaussOrthoPoly(deg);
pc = convert2affinePCE(2.0, 0.2, op);

@time ξ = sampleMeasure(ns,op);
@time samples = evaluatePCE(pc,ξ,op);