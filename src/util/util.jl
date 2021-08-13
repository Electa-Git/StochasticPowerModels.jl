################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

sorted_nw_ids(pm) = sort(collect(_PMs.nw_ids(pm)))

is_constrained(pm, cmp, idx) = _PMs.ref(pm, 1, cmp, idx, "cstr")

function check_deterministic_solution!(file, sol)
    # bus
    for nb in keys(file["bus"])
        file["bus"][nb]["cstr"] = 
            check_deterministic_voltage_magnitude_bounds(file["bus"][nb], sol["bus"][nb])
    end
    # gen
    for ng in keys(file["gen"])
        file["gen"][ng]["cstr"] =
            check_deterministic_generator_power_bounds(file["gen"][ng], sol["gen"][ng])
    end
    # branch 
    for nb in keys(file["branch"])
        fbus = file["branch"][nb]["f_bus"]
        file["branch"][nb]["imax"] = file["branch"][nb]["rate_a"] / file["bus"]["$fbus"]["vmin"]
        file["branch"][nb]["cstr"] =
            check_deterministic_branch_current_bounds(file["branch"][nb], sol["branch"][nb])
    end
end

function check_stochastic_solution!(file, sol)
    cstr = false

    # bus
    for nb in keys(file["nw"]["1"]["bus"])
        nb_cstr = check_stochastic_voltage_magnitude_bounds(file["nw"]["1"]["bus"][nb], sol, nb)
        if nb_cstr != file["nw"]["1"]["bus"][nb]["cstr"]  
            cstr = true
        end
        file["nw"]["1"]["bus"][nb]["cstr"] = nb_cstr            
    end
    # gen
    for ng in keys(file["nw"]["1"]["gen"])
        ng_cstr = check_stochastic_generator_power_bounds(file["nw"]["1"]["gen"][ng], sol, ng)
        if ng_cstr != file["nw"]["1"]["gen"][ng]["cstr"] 
            cstr = true
        end
        file["nw"]["1"]["gen"][ng]["cstr"] = ng_cstr
    end
    # branch 
    for nb in keys(file["nw"]["1"]["branch"])
        nb_cstr = check_stochastic_branch_current_bounds(file["nw"]["1"]["branch"][nb], sol, nb)
        if nb_cstr != file["nw"]["1"]["branch"][nb]["cstr"] 
            cstr = true
        end
        file["nw"]["1"]["branch"][nb]["cstr"] = nb_cstr
    end

    return cstr
end

function check_deterministic_voltage_magnitude_bounds(data, sol)
    vr, vi = sol["vr"], sol["vi"]
    vmin, vmax = data["vmin"], data["vmax"]

    vm = sqrt(vr^2 + vi^2)

    return isapprox(vm, vmin, rtol=1e-6) || isapprox(vm, vmax, rtol=1e-6)
end
function check_deterministic_generator_power_bounds(data, sol)
    pg, qg = sol["pg"], sol["qg"]
    pmin, pmax = data["pmin"], data["pmax"]
    qmin, qmax = data["qmin"], data["qmax"]

    return isapprox(pg, pmin, rtol=1e-6) || isapprox(pg, pmax, rtol=1e-6) ||
           isapprox(qg, qmin, rtol=1e-6) || isapprox(qg, qmax, rtol=1e-6)
end
function check_deterministic_branch_current_bounds(data, sol)
    csr, csi = sol["csr_fr"], sol["csi_fr"]
    cmax = data["imax"]

    cm = sqrt(csr^2 + csi^2)
    
    return isapprox(cm, cmax, rtol=1e-6)
end

function check_stochastic_voltage_magnitude_bounds(data, sol, nb)
    vs = [sol["nw"][nw]["bus"][nb]["vs"] for nw in sort(collect(keys(sol["nw"])))]
    
    vmin, vmax = data["vmin"], data["vmax"]
    λvmin, λvmax = data["λvmin"], data["λvmax"]

    vsmin = vs[1] - λvmin * sqrt(sum(vs[2:end].^2))
    vsmax = vs[1] + λvmax * sqrt(sum(vs[2:end].^2))

    return vsmin <= vmin^2 || vmax^2 <= vsmax
end
function check_stochastic_generator_power_bounds(data, sol, ng)
    pg = [sol["nw"][nw]["gen"][ng]["pg"] for nw in sort(collect(keys(sol["nw"])))]
    qg = [sol["nw"][nw]["gen"][ng]["qg"] for nw in sort(collect(keys(sol["nw"])))]

    pmin, pmax = data["pmin"], data["pmax"]
    qmin, qmax = data["qmin"], data["qmax"]
    λpmin, λpmax = data["λpmin"], data["λpmax"]
    λqmin, λqmax = data["λqmin"], data["λqmax"]

    pgmin = pg[1] - λpmin * sqrt(sum(pg[2:end].^2))
    pgmax = pg[1] + λpmax * sqrt(sum(pg[2:end].^2))
    qgmin = qg[1] - λqmin * sqrt(sum(qg[2:end].^2))
    qgmax = qg[1] + λqmax * sqrt(sum(qg[2:end].^2))

    return pgmin <= pmin || pmax <= pgmax || qgmin <= qmin || qmax <= qgmax
end
function check_stochastic_branch_current_bounds(data, sol, nb)
    css = [sol["nw"][nw]["branch"][nb]["css"] for nw in sort(collect(keys(sol["nw"])))]

    cmax = data["imax"]
    λimax = data["λimax"]

    csss = css[1] - λimax * sqrt(sum(css[2:end].^2))

    return cmax^2 <= csss
end