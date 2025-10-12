# StochasticPowerModels

<a href="https://github.com/timmyfaraday/StochasticPowerModels.jl/actions?query=workflow%3ACI"><img src="https://github.com/timmyfaraday/StochasticPowerModels.jl/workflows/CI/badge.svg"></img></a>
<a href="https://codecov.io/gh/timmyfaraday/StochasticPowerModels.jl"><img src="https://img.shields.io/codecov/c/github/timmyfaraday/StochasticPowerModels.jl?logo=Codecov"></img></a>


StochasticPowerModels.jl is an extension package of PowerModels.jl for 
Stochastic Power System Optimization.

Note that development is ongoing, and changes can be breaking without notice. We plan to register the package once we feel comfortable with the state of the implementation.

## Core Problem Specification

* Stochastic Optimal Power Flow (SOPF) for AC and AC/DC grids with RES
* Risk-based SOPF for AC and AC/DC grids with RES
* Stochastic Optimal Transmission Switching (SOTS) for AC and AC/DC grids with RES

## Core Network Formulation

- Exact
    - ACR
    - IVR 

## Core Stochastic Specification
For now, we only support Polynomial Chaos Expansion. We may add alternative stochastic optimization methods at a later stage.

## Network Data with Stochastic Data Extension

- Matpower ".m" files, extended to include:
    - stochastic germ: `mpc.sdata`,
    - stochastic bus data: `mpc.bus_sdata`, including: `dst_id`, `μ`, `σ`, `λvmin` and `λvmax`,
    - stochastic gen data: `mpc.gen_sdata`, including: `λpmin`, `λpmax`, `λqmin` and `λqmax`, and
    - stochastic branch data: `mpc.branch_sdata`, including: `λcmax`.

For an example, the user is referred to `/test/data/matpower/case5_spm.m`

## Installation

The latest stable release of StochasticPowerModels can be installed using the 
Julia package manager:

```
] add https://github.com/Electa-Git/StochasticPowerModels.jl.git
```

To test whether the package works, run:

```
] test StochasticPowerModels
```

## Acknowledgements

- Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)) _(Main Developer)_ , KU Leuven
    - SOPF for AC grids, IVR formulation   
- Kaan Yurtseven ([@kaanyurtseven](https://github.com/kaanyurtseven)), KU Leuven  
  - SOPF for AC/DC grids with RES, IVR formulation  
  - Risk-based SOPF for AC/DC grids with RES, IVR formulation  
  - SOTS for AC/DC grids with RES, IVR formulation
- Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven  
  - SOPF for AC Grids, ACR formulation
- Frederik Geth ([@frederikgeth](https://github.com/frederikgeth)), CSIRO  
  - SOPF for AC grids, reduced IVR formulation

## Citing StochasticPowerModels

If you find StochasticPowerModels useful in your work, we kindly request that you cite the following publications:

#### SOPF for hybrid AC/DC Grids with RES _(IVR Formulation)_:
- [Stochastic optimal power flow for hybrid AC/DC grids considering continuous non-Gaussian uncertainty](https://doi.org/10.1016/j.ijepes.2025.110828) 

#### SOPF for AC Grids _(IVR Formulation)_:
- [General polynomial chaos in the current–voltage formulation of the optimal power flow problem](https://doi.org/10.1016/j.epsr.2022.108472)

#### Risk-based SOPF for hybrid AC/DC Grids with RES _(IVR Formulation)_:
- [Risk-based Stochastic Optimal Power Flow for AC/DC Grids Using Polynomial Chaos Expansion](https://doi.org/10.1109/ISGTEUROPE62998.2024.10863605)

#### Security-constrained SOPF for hybrid AC/DC Grids with RES _(DC Formulation)_: 
- [Security-Constrained Stochastic Optimal Power Flow For Hybrid AC/DC Grids Considering Non-Gaussian Uncertainty](https://doi.org/10.1109/PowerTech59965.2025.11180443)

## License

This code is provided under a BSD license.
