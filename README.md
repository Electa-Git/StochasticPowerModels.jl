# StochasticPowerModels

<a href="https://github.com/timmyfaraday/StochasticPowerModels.jl/actions?query=workflow%3ACI"><img src="https://github.com/timmyfaraday/StochasticPowerModels.jl/workflows/CI/badge.svg"></img></a>
<a href="https://codecov.io/gh/timmyfaraday/StochasticPowerModels.jl"><img src="https://img.shields.io/codecov/c/github/timmyfaraday/StochasticPowerModels.jl?logo=Codecov"></img></a>

StochasticPowerModels.jl is an extension package of PowerModels.jl for 
Stochastic (Optimal) Power Flow. It is designed to enable inclusion of 
uncertainty in Steady-State Power Network Optimization. 

## Core Problem Specification

- Stochastic Power Flow (sPF)
- Stochastic Optimal Power Flow (sOPF)

## Core Network Formulation

- PowerModels.jl Formulation
    - IVR
    - ACR

## Core Stochastic Specification

- Polynomial Chaos Expansion

## Network Data Formats

- Matpower ".m" files
- PTI ".raw" files (PSS(R)E v33 specfication)
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
developer is Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)), 
with support from the following contributors:
- Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven, ACR formulation

## License

This code is provided under a BSD license.
