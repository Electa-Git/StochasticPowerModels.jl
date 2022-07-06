# LVDS data in StochasticPowerModels format

<a href="https://github.com/timmyfaraday/StochasticPowerModels.jl/actions?query=workflow%3ACI"><img src="https://github.com/timmyfaraday/StochasticPowerModels.jl/workflows/CI/badge.svg"></img></a>
<a href="https://codecov.io/gh/timmyfaraday/StochasticPowerModels.jl"><img src="https://img.shields.io/codecov/c/github/timmyfaraday/StochasticPowerModels.jl?logo=Codecov"></img></a>
<a href="https://timmyfaraday.github.io/StochasticPowerModels.jl/"><img src="https://github.com/timmyfaraday/StochasticPowerModels.jl/workflows/Documentation/badge.svg"></img></a>


StochasticPowerModels.jl is an extension package of PowerModels.jl for 
Stochastic (Optimal) Power Flow. It is designed to enable inclusion of 
uncertainty in Steady-State Power Network Optimization. 

Note that development is ongoing, and changes can be breaking without notice. We plan to register the package once we feel comfortable with the state of the implementation.

## Formulation description

- Stochastic Optimal Power Flow (sOPF)

## HC Formulations
Using IVR power flow
- Deterministic
    OPF-HC: Assumes high irradiance and 0.1 kW load for each consumer
-Stochastic
	gPC-CC-OPF: Assumes uncertainty in Irradiance and load and chance constraints in nodal voltage and current

	

## Core Stochastic Specification
For now, we only support Polynomial Chaos Expansion. We may add alternative stochastic optimization methods at a later stage.

- Polynomial Chaos Expansion
    - with auxiliary variables/constraints

## Network Data with Stochastic Data Extension
The Folder Spanish has Spanish network on folder All_feeder in JSON format. 

The file `CreatePMDDictionary` is the parser file to convert the JSON file into `PowerModels.jl` format.

The original dataset consists of a full Low voltage network of sub-urban region with 30 transformers, 160 feeders, 10290 nodes and 8087 consumers, with load profiles of 20 days from actual smart-meter.
The paper is available in https://www.sciencedirect.com/science/article/pii/S0142061519318836

and the link for full data-set: https://data.mendeley.com/datasets/685vgp64sm/1

For the original networks, the line impedance is specified 4x4 matrice without mutual impedance and the load from smart meter data for 20 days.

However, in this part the load and irradiance are defined as Beta distribution in folders `beta_lm_2016_8_6.csv` and `beta_pm_2016_8_6.csv` respectively for a high irradiance day in spring. 

For each feeder there are 4 JSON file describing the feeder topology:
	- <feeder_name>_configuration.json, 
	- <feeder_name>_branches.json, 
	- <feeder_name>_buses.json, and 
	- <feeder_name>_devices.json 
and 1 csv file:
  	_<feeder_name>.csv- which is the linking file between devices and the load uncertainty. 

- Matpower ".m" files, extended to include:
    - stochastic germ: `mpc.sdata`,
    - stochastic bus data: `mpc.bus_sdata`, including: `dst_id`, `μ`, `σ`, `λvmin` and `λvmax`,
    - stochastic gen data: `mpc.gen_sdata`, including: `λpmin`, `λpmax`, `λqmin` and `λqmax`, and
    - stochastic branch data: `mpc.branch_sdata`, including: `λcmax`.


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

The primary developer is Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)), 
with support from the following contributors:
- Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven, ACR formulation
- Frederik Geth ([@frederikgeth](https://github.com/frederikgeth)), CSIRO, reduced IVR formulation

## License

This code is provided under a BSD license.
