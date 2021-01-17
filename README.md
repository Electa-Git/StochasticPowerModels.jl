# StochasticPowerModels

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://timmyfaraday.github.io/StochasticPowerModels.jl/dev)
[![Build Status](https://travis-ci.com/timmyfaraday/StochasticPowerModels.jl.svg?branch=master)](https://travis-ci.com/timmyfaraday/StochasticPowerModels.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/timmyfaraday/StochasticPowerModels.jl?svg=true)](https://ci.appveyor.com/project/timmyfaraday/StochasticPowerModels-jl)
[![Codecov](https://codecov.io/gh/timmyfaraday/StochasticPowerModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/timmyfaraday/StochasticPowerModels.jl)

StochasticPowerModels.jl is an extension package of PowerModels(Distribution).jl
for Stochastic (Optimal) Power Flow. It is designed to enable inclusion of 
uncertainty in Steady-State Power Network Optimization. 

## Core Problem Specification

- Stochastic Power Flow (sPF)
- Stochastic Optimal Power Flow (sOPF)

## Core Network Formulation

- PowerModels.jl Formulation
    - IVR
- PowerModelsDistribution.jl Formulation
    - IVR

## Core Stochastic Specification

- Polynomial Chaos Expansion

## Network Data Formats

- Matpower ".m" files
- PTI ".raw" files (PSS(R)E v33 specfication)
- OpenDSS ".dss" files in the PowerModelsDistribution format
- CSV ".csv" file with uncertainty data

## Installation

The latest stable release of StochasticPowerModels can be installed using the 
Julia package manager:

```
] add https://github.com/timmyfaraday/StochasticPowerModels.jl.git
```

In order to test whether the package works, run:

```
] test StochasticPowerModels
```

## Acknowledgements

This code has been developed at KU Leuven (University of Leuven). The primary
developers are Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday))
and Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala))

## License

This code is provided under a BSD license.