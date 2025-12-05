# Network Formulations

## Type Hierarchy
We begin with the top of the hierarchy, where we can distinguish between current-voltage and power-voltage formulations
```julia
AbstractACPModel <: AbstractPowerModel
AbstractIVRModel <: AbstractPowerModel
```

## Power Models
Each of these forms can be used as the model parameter for a PowerModel:
```julia
ACPPowerModel <: AbstractACPForm
IVRPowerModel <: AbstractIVRModel
```