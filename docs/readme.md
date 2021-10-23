# Building the Documentation for StochasticPowerModels.jl

## Installation
We rely on [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). To install it, run the following command in a julia session:

```julia
Pkg.add("Documenter")
```

## Building the Docs
To preview the html output of the documents, run the following command:

```julia
julia --color=yes make.jl
```

You can then view the documents in `build/index.html`.

To generate the pdf, run 

```julia
julia --pdf make.jl
```

You will need to have an active docker installation on your local machine. [Install guide](https://docs.docker.com/desktop/windows/install/).

**Warning**: Do not `git commit` the contents of build (or any other content generated by Documenter) to your repository's master branch. This helps to avoid including unnessesary changes for anyone reviewing commits that happen to include documentation changes.

For further details, please read the [documentation for Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/).