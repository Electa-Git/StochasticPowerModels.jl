################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
""
function variable_unit_current(pm::AbstractIVRModel; nw::Int=pm.cnw, 
                               bounded::Bool=true, report::Bool=true, kwargs...)
    variable_unit_current_real(pm, nw=nw, bounded=bounded, report=report, kwargs...)
    variable_unit_current_imaginary(pm, nw=nw, bounded=bounded, report=report, kwargs...)

    # store active and reactive power expressions for use in objective + post processing
    ntw = nws(pm)
    pu, qu = Dict(), Dict()
    for (i,unit) in ref(pm, nw, :unit)
        id = gen["unit_bus"]
        vr = [var(pm, n, :vr, id) for n in nws(pm)]
        vi = [var(pm, n, :vi, id) for n in nws(pm)]
        cr = [var(pm, n, :cru, i) for n in nws(pm)]
        ci = [var(pm, n, :ciu, i) for n in nws(pm)]
        pu[i] = _JMP.@NLexpression(pm.model, 
                    sum(vr[n1]*cr[n2] + vi[n1]*ci[n2] for n1 in ntw, n2 in ntw))
        qu[i] = _JMP.@NLexpression(pm.model, 
                    sum(vi[n1]*cr[n2] - vr[n1]*ci[n2] for n1 in ntw, n2 in ntw))
    end
    var(pm, nw)[:pu] = pu
    var(pm, nw)[:qu] = qu
    report && _IM.sol_component_value(pm, nw, :unit, :pu, ids(pm, nw, :unit), pu)
    report && _IM.sol_component_value(pm, nw, :unit, :qu, ids(pm, nw, :unit), qu)

end
