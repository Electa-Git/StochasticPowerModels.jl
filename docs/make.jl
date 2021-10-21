push!(LOAD_PATH,"../src/")
using Documenter

makedocs(sitename="StochasticPowerModels Documentation")

# A flag to check if we are running in a GitHub action.
const _IS_GITHUB_ACTIONS = get(ENV, "GITHUB_ACTIONS", "false") == "true"

# Pass --pdf to build the PDF. On GitHub actions, we always build the PDF.
const _PDF = findfirst(isequal("--pdf"), ARGS) !== nothing || _IS_GITHUB_ACTIONS
# const _PDF = true

if _PDF
    latex_platform = _IS_GITHUB_ACTIONS ? "docker" : "native"
    @time Documenter.makedocs(
        sitename = "StochasticPowerModels",
        authors = "The StochasticPowerModels core developers and contributors",
        format = Documenter.LaTeX(platform = latex_platform),
        build = "latex_build",
    )
    # Hack for deploying: copy the pdf (and only the PDF) into the HTML build
    # directory! We don't want to copy everything in `latex_build` because it
    # includes lots of extraneous LaTeX files.
    cp(
        joinpath(@__DIR__, "latex_build", "StochasticPowerModels.pdf"),
        joinpath(@__DIR__, "build", "StochasticPowerModels.pdf"),
    )
end

Documenter.deploydocs(
    repo = "github.com/timmyfaraday/StochasticPowerModels.jl.git",
    push_preview = true,
)