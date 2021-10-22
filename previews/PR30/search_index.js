var documenterSearchIndex = {"docs":
[{"location":"formulations/#Network-Formulations","page":"Network Formulations","title":"Network Formulations","text":"","category":"section"},{"location":"formulations/#Type-Hierarchy","page":"Network Formulations","title":"Type Hierarchy","text":"","category":"section"},{"location":"formulations/","page":"Network Formulations","title":"Network Formulations","text":"We begin with the top of the hierarchy, where we can distinguish between current-voltage and power-voltage formulations","category":"page"},{"location":"formulations/","page":"Network Formulations","title":"Network Formulations","text":"AbstractACPModel <: AbstractPowerModel\nAbstractIVRModel <: AbstractPowerModel","category":"page"},{"location":"formulations/#Power-Models","page":"Network Formulations","title":"Power Models","text":"","category":"section"},{"location":"formulations/","page":"Network Formulations","title":"Network Formulations","text":"Each of these forms can be used as the model parameter for a PowerModel:","category":"page"},{"location":"formulations/","page":"Network Formulations","title":"Network Formulations","text":"ACPPowerModel <: AbstractACPForm\nIVRPowerModel <: AbstractIVRModel","category":"page"},{"location":"quickguide/#Quick-Start-Guide","page":"Getting Started","title":"Quick Start Guide","text":"","category":"section"},{"location":"quickguide/","page":"Getting Started","title":"Getting Started","text":"Once StochasticPowerModels is installed, Ipopt is installed, and a network data file (e.g. \"case5_spm.m\" ) has been acquired, an stochatic AC Optimal Power Flow can be executed with,","category":"page"},{"location":"quickguide/","page":"Getting Started","title":"Getting Started","text":"using StochasticPowerModels\nusing Ipopt\n\nresult = run_ac_sopf(\"matpower/case5_spm.m\", Ipopt.Optimizer)","category":"page"},{"location":"math-model/#The-StochasticPowerModels-Mathematical-Model","page":"Mathematical Model","title":"The StochasticPowerModels Mathematical Model","text":"","category":"section"},{"location":"math-model/#Sets-and-Parameters","page":"Mathematical Model","title":"Sets and Parameters","text":"","category":"section"},{"location":"math-model/","page":"Mathematical Model","title":"Mathematical Model","text":"StochasticPowerModels implements a generalized polynomial chaos expansion version of the AC Optimal Power Flow problem from Matpower.   The core generalizations of the deterministic OPF problem are,","category":"page"},{"location":"math-model/","page":"Mathematical Model","title":"Mathematical Model","text":"Support for multiple load (S^d_k) and shunt (Y^s_k) components on each bus i\nLine charging that supports a conductance and asymmetrical values (Y^c_ij Y^c_ji)","category":"page"},{"location":"math-model/","page":"Mathematical Model","title":"Mathematical Model","text":"beginalign\n\nmboxsets  nonumber \n N mbox - busesnonumber \n R mbox - reference busesnonumber \n E E^R mbox - branches forward and reverse orientation nonumber \n G G_i mbox - generators and generators at bus i nonumber \n L L_i mbox - loads and loads at bus i nonumber \n S S_i mbox - shunts and shunts at bus i nonumber \n\nmboxdata  nonumber \n S^gl_k S^gu_k  forall k in G nonumber mbox - generator complex power bounds\n c_2k c_1k c_0k  forall k in G nonumber  mbox - generator cost components\n v^l_i v^u_i  forall i in N nonumber mbox - voltage bounds\n S^d_k  forall k in L nonumber mbox - load complex power consumption\n Y^s_k  forall k in S nonumber mbox - bus shunt admittance\n Y_ij Y^c_ij Y^c_ji  forall (ij) in E nonumber mbox - branch pi-section parameters\n T_ij  forall (ij) in E nonumber mbox - branch complex transformation ratio\n s^u_ij   forall (ij) in E nonumber mbox - branch apparent power limit\n i^u_ij   forall (ij) in E nonumber mbox - branch current limit\n theta^Delta l_ij theta^Delta u_ij  forall (ij) in E nonumber mbox - branch voltage angle difference bounds\n\nendalign","category":"page"},{"location":"math-model/#Stochastic-Optimal-Power-Flow-in-Current-Voltage-Variables","page":"Mathematical Model","title":"Stochastic Optimal Power Flow in Current-Voltage Variables","text":"","category":"section"},{"location":"math-model/","page":"Mathematical Model","title":"Mathematical Model","text":"A variable I^s_ij, representing the current in the direction i to j, through the series part of the pi-section, is used. The mathematical structure for a current-voltage formulation is conceived as:","category":"page"},{"location":"math-model/","page":"Mathematical Model","title":"Mathematical Model","text":"beginalign\n\nmboxvariables   nonumber \n I^g_k  forall kin G nonumber \n V_i  forall iin N nonumber \n I^s_ij  forall (ij) in E cup E^R  mbox - branch complex (series) current\n I_ij  forall (ij) in E cup E^R  mbox - branch complex (total) current labelvar_total_current\n\nmboxminimize   sum_k in G c_2k (Re(S^g_k))^2 + c_1kRe(S^g_k) + c_0k nonumber\n\nmboxsubject to   nonumber \n angle V_r = 0   forall r in R nonumber \n S^gl_k leq Re(V_i (I^g_k)^*) + j Im(V_i (I^g_k)^*) leq S^gu_k  forall k in G   labeleq_complex_power_definition_gen\n v^l_i leq V_i leq v^u_i  forall i in N nonumber\n sum_substackk in G_i I^g_k - sum_substackk in L_i (S^d_kV_i)^* - sum_substackk in S_i Y^s_k V_i = sum_substack(ij)in E_i cup E_i^R I_ij  forall iin N  labeleq_kcl_current \n I_ij =  fracI^s_ijT_ij^* + Y^c_ij fracV_iT_ij^2   forall (ij)in E labeleq_current_from \n I_ji = -I^s_ij + Y^c_ji V_j   forall (ij)in E labeleq_current_to \n fracV_iT_ij = V_j + z_ij I^s_ij   forall (ij) in E labeleq_ohms_iv \n S_ij = V_i I_ij leq s^u_ij  forall (ij) in E cup E^R nonumber\n I_ij leq i^u_ij  forall (ij) in E cup E^R nonumber\n theta^Delta l_ij leq angle (V_i V^*_j) leq theta^Delta u_ij  forall (ij) in E nonumber\n\nendalign","category":"page"},{"location":"math-model/#Stochastic-Optimal-Power-Flow-in-Power-Voltage-Variables","page":"Mathematical Model","title":"Stochastic Optimal Power Flow in Power-Voltage Variables","text":"","category":"section"},{"location":"math-model/","page":"Mathematical Model","title":"Mathematical Model","text":"A complete mathematical model is as follows,","category":"page"},{"location":"math-model/","page":"Mathematical Model","title":"Mathematical Model","text":"beginalign\n\nmboxvariables   nonumber \n S^g_k  forall kin G mbox - generator complex power dispatch labelvar_generation\n V_i  forall iin N labelvar_voltage mbox - bus complex voltage\n S_ij  forall (ij) in E cup E^R  labelvar_complex_power mbox - branch complex power flow\n\nmboxminimize   sum_k in G c_2k (Re(S^g_k))^2 + c_1kRe(S^g_k) + c_0k labeleq_objective\n\nmboxsubject to   nonumber \n angle V_r = 0   forall r in R labeleq_ref_bus\n S^gl_k leq S^g_k leq S^gu_k  forall k in G  labeleq_gen_bounds\n v^l_i leq V_i leq v^u_i  forall i in N labeleq_voltage_bounds\n sum_substackk in G_i S^g_k - sum_substackk in L_i S^d_k - sum_substackk in S_i (Y^s_k)^* V_i^2 = sum_substack(ij)in E_i cup E_i^R S_ij  forall iin N labeleq_kcl_shunt \n S_ij = left( Y_ij + Y^c_ijright)^* fracV_i^2T_ij^2 - Y^*_ij fracV_i V^*_jT_ij  forall (ij)in E labeleq_power_from\n S_ji = left( Y_ij + Y^c_ji right)^* V_j^2 - Y^*_ij fracV^*_i V_jT^*_ij  forall (ij)in E labeleq_power_to\n S_ij leq s^u_ij  forall (ij) in E cup E^R labeleq_thermal_limit\n I_ij leq i^u_ij  forall (ij) in E cup E^R labeleq_current_limit\n theta^Delta l_ij leq angle (V_i V^*_j) leq theta^Delta u_ij  forall (ij) in E labeleq_angle_difference\n\nendalign","category":"page"},{"location":"#StochasticPowerModels.jl-Documentation","page":"Home","title":"StochasticPowerModels.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = StochasticPowerModels","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"StochasticPowerModels.jl is a research-grade Julia/JuMP package for experimentation with Steady-State Power Network Optimization under uncertainty, extending PowerModels.jl.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For now, StochasticPowerModels is unregistered. Nevertheless, you can install it through","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add https://github.com/timmyfaraday/StochasticPowerModels.jl.git","category":"page"},{"location":"","page":"Home","title":"Home","text":"At least one solver is required for running StochasticPowerModels.  The open-source solver Ipopt is recommended, as it is fast, scaleable and can be used to solve a wide variety of the problems and network formulations provided in PowerModels.  The Ipopt solver can be installed via the package manager with","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add Ipopt","category":"page"},{"location":"","page":"Home","title":"Home","text":"Test that the package works by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"] test StochasticPowerModels","category":"page"}]
}
