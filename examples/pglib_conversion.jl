using PowerModels
const PM = PowerModels

pglibpath = "C:\\Users\\Frederik Geth\\Downloads\\pglibsmall"

function add_stochastic_data!(data)
    for (i,gen) in data["gen"]
        gen["λpmin"] = 1.03643
        gen["λpmax"] = 1.03643
        gen["λqmin"] = 1.03643
        gen["λqmax"] = 1.03643
    end

    for (l,branch) in data["branch"]
        branch["λcmax"] = 1.03643
    end

    data["sdata"] = Dict()

    data["sdata"]["1"] = Dict()
    data["sdata"]["1"]["dst"] = "Beta" 
    data["sdata"]["1"]["pa"] = 2.0
    data["sdata"]["1"]["pb"] = 2.0 

    data["sdata"]["2"] = Dict()
    data["sdata"]["2"]["dst"] = "Normal" 
    data["sdata"]["2"]["pa"] = 0.0
    data["sdata"]["2"]["pb"] = 0.0 
    
end

for file in readdir(pglibpath)
    data = PM.parse_file(joinpath(pglibpath, file))

    add_stochastic_data!(data)

    newfilename = file[1:end-2]*"_spm.m"
    PM.export_file(joinpath(pglibpath, newfilename), data)

    #TODO add writer for sdata
end

    