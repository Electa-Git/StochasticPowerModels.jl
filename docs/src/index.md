# HarmonicPowerModels.jl Documentation

```@meta
CurrentModule = HarmonicPowerModels
```

## Overview

HarmonicPowerModels.jl is a research-grade Julia/JuMP package for experimentation with Steady-State Power Network Optimization with Harmonics, extending PowerModels.jl.


## Installation

For now, HarmonicPowerModels is unregistered. Nevertheless, you can install it through

```
] add https://github.com/timmyfaraday/HarmonicPowerModels.jl.git
```

At least one solver is required for running HarmonicPowerModels.  The open-source solver Ipopt is recommended, as it is fast, scaleable and can be used to solve a wide variety of the problems and network formulations provided in PowerModels.  The Ipopt solver can be installed via the package manager with

```julia
] add Ipopt
```

Test that the package works by running

```julia
] test HarmonicPowerModels
```
