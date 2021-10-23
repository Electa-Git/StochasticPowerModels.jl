# The StochasticPowerModels Mathematical Model


## Sets and Parameters

StochasticPowerModels implements a generalized polynomial chaos expansion version of the AC Optimal Power Flow problem from [Matpower](http://www.pserc.cornell.edu/matpower/).  
The core generalizations of the deterministic OPF problem are,

- Support for multiple load ($S^d_k$) and shunt ($Y^s_{k}$) components on each bus $i$
- Line charging that supports a conductance and asymmetrical values ($Y^c_{ij}, Y^c_{ji}$)


```math
\begin{align}
%
\mbox{sets:} & \nonumber \\
& I \mbox{ - buses}\nonumber \\
& R \mbox{ - reference buses}\nonumber \\
& E^F, E^R \mbox{ - branches, forward and reverse orientation} \nonumber \\
& G, G_i \mbox{ - generators and generators at bus $i$} \nonumber \\
& L, L_i \mbox{ - loads and loads at bus $i$} \nonumber \\
& S, S_i \mbox{ - shunts and shunts at bus $i$} \nonumber \\
& K \mbox{ - polynomial chaos basis} \\
%
\mbox{data:} & \nonumber \\
& S^{gl}_k, S^{gu}_k \;\; \forall k \in G \nonumber \mbox{ - generator complex power bounds}\\
& c_{2k}, c_{1k}, c_{0k} \;\; \forall k \in G \nonumber  \mbox{ - generator cost components}\\
& v^l_i, v^u_i \;\; \forall i \in N \nonumber \mbox{ - voltage bounds}\\
& S^d_k \;\; \forall k \in L \nonumber \mbox{ - load complex power consumption}\\
& Y^s_{k} \;\; \forall k \in S \nonumber \mbox{ - bus shunt admittance}\\
& Y_{ij}, Y^c_{ij}, Y^c_{ji} \;\; \forall (i,j) \in E \nonumber \mbox{ - branch pi-section parameters}\\
& {T}_{ij} \;\; \forall (i,j) \in E \nonumber \mbox{ - branch complex transformation ratio}\\
& s^u_{ij}  \;\; \forall (i,j) \in E \nonumber \mbox{ - branch apparent power limit}\\
& i^u_{ij}  \;\; \forall (i,j) \in E \nonumber \mbox{ - branch current limit}\\
& \theta^{\Delta l}_{ij}, \theta^{\Delta u}_{ij} \;\; \forall (i,j) \in E \nonumber \mbox{ - branch voltage angle difference bounds}\\
%
\end{align}
```

## Stochastic Optimal Power Flow in Current-Voltage Variables
A variable $I^{s}_{ij}$, representing the current in the direction $i$ to $j$, through the series part of the pi-section, is used.
The mathematical structure for a current-voltage formulation is conceived as:

```math
\begin{align}
%
\mbox{variables: } & \nonumber \\
& I_{g,k} \;\; \forall g \in G, k \in K 
\mbox{ - generator complex current} 
\label{var_gen_current_ivr} \nonumber \\
& V_{i,k} \;\; \forall i \in I, k \in K
\mbox{ - bus complex voltage} 
\label{var_bus_voltage_ivr} \nonumber \\
& I^{s}_{ij,k} \;\; \forall (i,j) \in E^F \cup E^R, k \in K 
\mbox{ - branch complex (series) current} 
\label{var_series_branch_current_ivr} \nonumber \\
& I_{ij,k} \;\; \forall (i,j) \in E^F \cup E^R, k \in K 
\mbox{ - branch complex (total) current} 
\label{var_total_branch_current_ivr} \nonumbed \\
%
\mbox{minimize: } & 
\sum_{g \in G} \mathbb{E} \left[ 
    c^2_{g} (\Re(\bar{S_g}))^2 + c^1_{g} \Re(\bar{S_g}) + c^0_{g} 
\right]    
\nonumber\\
%
\mbox{subject to: } & \nonumber \\
& \mathbb{E} \[\bar{\angle V_{r}} \] = 0 \;\; 
\forall r \in R 
\label{eq_expected_reference_bus_voltage_ivr} \\
& |V_{r,k}| = 0 \;\;
\forall r \in R, k \in K_{0} 
\label{eq_reference_bus_voltage_ivr} \\
& S_{u,k} = \sum_{ \in K} M_{,k} \Re(V_{i,} (I_{u,})^*) + j \Im(V_{i,} (I_{u,})^*) \;\; 
\forall u \in G \cup L, k \in K: u \to i in I   
\label{eq_complex_gen_power_definition_ivr} \\
& \mathbb{P} \[ S^{min}_{g} \leq S_g \] \geq 1 - \varepsilon \;\;
\forall g \in G
\label{eq_complex_gen_power_lb_ivr} \\
& \mathbb{P} \[ S_g \leq S^{max}_{g} \] \geq 1 - \varepsilon \;\;
\forall g \in G
\label{eq_complex_gen_power_ub_ivr} \\
& \mathbb{P} \[ v^{min}_i \leq |V_i| \] \geq 1 - \varepsilon \;\; 
\forall i \in I 
\label{eq_voltage_bus_lb_ivr} \\
& \mathbb{P} \[ |V_i| \leq v^{max}_i \] \geq 1 - \varepsilon \;\; 
\forall i \in I 
\label{eq_bus_voltage_ub_ivr} \\
&   \sum_{\substack{g \in G_i}} I_{g,k} - 
    \sum_{\substack{l \in L_i}} I_{l,k} - 
    \sum_{\substack{s \in S_i}} Y_s V_{i,k} 
    = 
    \sum_{\substack{(i,j) \in E_i^F \cup E_i^R}} I_{ij,k} \;\; 
\forall i \in I, k \in K
\label{eq_kcl_current_ivr} \\
& I_{ij,k} = \frac{I^{s}_{ij,k}}{T_{ij}^*} + Y^c_{ij} \frac{V_{i,k}}{|T_{ij}|^2} \;\; \forall (i,j) \in E^F, k \in K 
\label{eq_current_from_ivr} \\
& I_{ji,k} = -I^{s}_{ij,k} + Y^c_{ji} V_{j,k}  \;\; 
\forall (j,i) \in E^R, k \in K 
\label{eq_current_to_ivr} \\
& \frac{V_{i,k}}{{T}_{ij}} = V_{j,k} + z_{ij} I^{s}_{ij,k}  \;\; 
\forall (i,j) \in E^F, k \in K 
\label{eq_ohms_ivr} \\
& \mathbb{P} \[ |I_{ij}| \leq i^u_{ij} \] \geq 1 - \varepsilon \;\; 
\forall (i,j) \in E^F \cup E^R 
\label{eq_branch_current_ub_ivr} \\
%& \theta^{\Delta l}_{ij} \leq \angle (V_i V^*_j) \leq \theta^{\Delta u}_{ij} \;\; \forall (i,j) \in E \nonumber --> is this necessary?
%& |S_{ij,k}| = |V_{i}| |I_{ij}| \leq s^u_{ij} \;\; \forall (i,j) \in E \cup E^R \nonumber\\ --> is this necessary?
%
\end{align}
```

##  Stochastic Optimal Power Flow in Power-Voltage Variables

A complete mathematical model is as follows,

```math
\begin{align}
%
\mbox{variables: } & \nonumber \\
& S^g_k \;\; \forall k\in G \mbox{ - generator complex power dispatch} \label{var_generation}\\
& V_i \;\; \forall i\in N \label{var_voltage} \mbox{ - bus complex voltage}\\
& S_{ij} \;\; \forall (i,j) \in E \cup E^R  \label{var_complex_power} \mbox{ - branch complex power flow}\\
%
\mbox{minimize: } & \sum_{k \in G} c_{2k} (\Re(S^g_k))^2 + c_{1k}\Re(S^g_k) + c_{0k} \label{eq_objective}\\
%
\mbox{subject to: } & \nonumber \\
& \angle V_{r} = 0  \;\; \forall r \in R \label{eq_ref_bus}\\
& S^{gl}_k \leq S^g_k \leq S^{gu}_k \;\; \forall k \in G  \label{eq_gen_bounds}\\
& v^l_i \leq |V_i| \leq v^u_i \;\; \forall i \in N \label{eq_voltage_bounds}\\
& \sum_{\substack{k \in G_i}} S^g_k - \sum_{\substack{k \in L_i}} S^d_k - \sum_{\substack{k \in S_i}} (Y^s_k)^* |V_i|^2 = \sum_{\substack{(i,j)\in E_i \cup E_i^R}} S_{ij} \;\; \forall i\in N \label{eq_kcl_shunt} \\
& S_{ij} = \left( Y_{ij} + Y^c_{ij}\right)^* \frac{|V_i|^2}{|{T}_{ij}|^2} - Y^*_{ij} \frac{V_i V^*_j}{{T}_{ij}} \;\; \forall (i,j)\in E \label{eq_power_from}\\
& S_{ji} = \left( Y_{ij} + Y^c_{ji} \right)^* |V_j|^2 - Y^*_{ij} \frac{V^*_i V_j}{{T}^*_{ij}} \;\; \forall (i,j)\in E \label{eq_power_to}\\
& |S_{ij}| \leq s^u_{ij} \;\; \forall (i,j) \in E \cup E^R \label{eq_thermal_limit}\\
& |I_{ij}| \leq i^u_{ij} \;\; \forall (i,j) \in E \cup E^R \label{eq_current_limit}\\
& \theta^{\Delta l}_{ij} \leq \angle (V_i V^*_j) \leq \theta^{\Delta u}_{ij} \;\; \forall (i,j) \in E \label{eq_angle_difference}
%
\end{align}
```
