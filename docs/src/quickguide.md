# Quick Start Guide

Once StochasticPowerModels is installed, Ipopt is installed, and a network data file (e.g. `"case5_spm.m"` ) has been acquired, an stochatic AC Optimal Power Flow can be executed with,

```julia
using StochasticPowerModels
using Ipopt

result = run_ac_sopf("matpower/case5_spm.m", Ipopt.Optimizer)
```

