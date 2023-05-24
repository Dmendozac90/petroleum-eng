"""
Pressure Traverse Computation

@author: Daniel Mendoza
@email: dmendozac90@gmail.com
"""

from . import multiphase_flow_params
from . import beggs_brill_params
from . import pressure_traverse_utilities

def pressure_traverse(
        P_inlet, P_outlet, T_inlet, dT_dz, q_o_sc, q_w_sc, q_g_sc, R_p, API, 
        gamma_w, gamma_g, P_b, R_sb, R_s, beta_o, beta_w, z, beta_g, mu_o, 
        mu_w, mu_g, sigma_o, sigma_w, c_o, TDS, d, theta, delta, EOT, nodes, 
        outlet_depth=0, dP_dz=0.002, method="beggs-brill", P_sep=None, 
        T_sep=None, tol=1e-2, depth=None, pressure=None
    ):
    """
    Recursive implementation of pressure loss calculation across a segment of 
    a wellbore. 

    Parameters
    ----------
    P_inlet : inlet pressure (psi)
    P_outlet : outlet pressure (psi)
    T_inlet : inlet temperature (°F)
    dT_dz : temperature gradient (°F/ft)
    q_o_sc : surface oil rate (STB/D)
    q_w_sc : surface water rate (STB/D)
    q_g_sc : surface gas rate (Mscf/D)
    R_p : producing gas-oil ratio (scf/STB)
    API : oil gravity (°API)
    gamma_w : water gravity (unitless)
    gamma_g : gas gravity (unitless)
    P_b : bubble-point pressure (psia)
    R_sb : solution gas-oil ratio at bubble-point pressure (scf/STB)
    R_s : solution gas-oil ratio (scf/STB)
    beta_o : oil formation volume factor (bbl/STB)
    beta_w : water formation volume factor (bbl/STB)
    z : compressibility factor (unitless)
    beta_g : gas formation volume factor (ft**3/STB)
    mu_o : oil viscosity (cp)
    mu_w : water viscosity (cp)
    mu_g : gas viscosity (cp)
    sigma_o : oil surface tension (dyne/cm)
    sigma_w : water surface tension (dyne/cm)
    c_o : oil compressibility (psi**-1)
    TDS : total dissolved solids (ppm)
    d : pipe inner diameter (in)
    theta : inclination (°)
    delta : pipe roughness (ft)
    EOT : end-of-tubing depth (ft)
    nodes : nodes
    outlet_depth : outlet depth (ft)
    dP_dz : pressure gradient for initial outlet depth (psi/ft)
    method : method for pressure drop calculation 
    P_sep : seperator pressure (only needed for vasquez-beggs correlations)
    T_sep : seperator temperature (only needed for vasquez-beggs correlations)
    tol : error tolerance
    depth : pressure traverse depth
    pressure : pressure traverse
    
    Returns
    -------
    depth : node depth (ft)
    pressure : node pressure (psi)
    """
    
    if depth is None:
        depth = []
    if pressure is None:
        pressure = []

    dz = EOT / nodes
    outlet_depth += dz
    T_outlet =  T_inlet + dT_dz * dz
    #terminal condition
    if outlet_depth > EOT:
        return None, None
    #implement recursion
    else:
        P_avg = (P_inlet + P_outlet) / 2
        T_avg = (T_inlet + T_outlet) / 2
        #check if reservoir fluid is saturated
        saturated = P_avg <= P_b
        #saturated reservoir fluid properties
        if saturated:
            R_s_avg, beta_o_avg, beta_w_avg, z_avg, beta_g_avg, mu_o_avg, \
                mu_w_avg, mu_g_avg = \
                    pressure_traverse_utilities.saturated_reservoir_params(
                        P_avg, T_avg, API, gamma_g, TDS, P_sep, T_sep, R_s, \
                        beta_o, beta_w, z, beta_g, mu_o, mu_w, mu_g
                    )
        #undersaturated reservoir fluid properties
        else:
            R_s_avg, beta_o_avg, beta_w_avg, beta_g_avg, mu_o_avg, mu_w_avg, \
                mu_g_avg = \
                pressure_traverse_utilities.undersaturated_reservoir_params(
                P_avg, T_avg, API, gamma_g, TDS, P_sep, T_sep, R_sb, P_b, \
                beta_o, beta_w, mu_o, mu_w, c_o
                )

        q_o, q_w, q_l, q_g = pressure_traverse_utilities.volumetric_rates(
            q_o_sc, q_w_sc, q_g_sc, beta_o_avg, beta_w_avg, beta_g_avg, 
            R_s_avg
        )

        rho_o, rho_w, rho_g = pressure_traverse_utilities.insitu_densities(
            q_o, q_w, q_g, q_o_sc, q_w_sc, API, gamma_w, gamma_g, R_s_avg, R_p
        )

        lambda_l, f_o, f_w = pressure_traverse_utilities.fluid_fractions(
            q_l, q_o, q_w, q_g
        )
        
        v_sl, v_sg, v_m = pressure_traverse_utilities.fluid_velocities(
            q_l, q_g, d
        )

        rho_l, mu_l, sigma_l = pressure_traverse_utilities.liquid_params(
            rho_o, rho_w, mu_o_avg, mu_w_avg, sigma_o, sigma_w, f_o, f_w
        )

        N_Lv, N_gv, N_d, N_L = pressure_traverse_utilities.duns_ross_params(
            rho_l, mu_l, sigma_l, v_sl, v_sg, d
        )        
        
        rho_n, mu_n = pressure_traverse_utilities.no_slip_params(
            rho_l, rho_g, mu_l, mu_g_avg, lambda_l
        )
        
        if method == 'beggs-brill':
            H_l_theta = pressure_traverse_utilities.beggs_brill_liquid_holdup(
                N_Lv, v_m, lambda_l, d, theta
            )

            f_two_phase = pressure_traverse_utilities.beggs_brill_friction(
                rho_n, mu_n, lambda_l, v_m, H_l_theta, d, delta
            )

            rho_s = multiphase_flow_params.insitu_slip_mixture_density(
                rho_l, rho_g, H_l_theta
            )

        dp_dz = multiphase_flow_params.pressure_loss_gradient(
            rho_n, rho_s, v_m, f_two_phase, d, theta
        )

        P_calc = P_inlet + dp_dz * dz
        error = abs(P_calc - P_outlet) / P_calc  
        #evaluate if converged
        while error > tol:
            #recurse if not converged
            outlet_depth -= dz
            P_outlet = P_calc

            return pressure_traverse(
                P_inlet, P_outlet, T_inlet, dT_dz, q_o_sc, q_w_sc, q_g_sc, 
                R_p, API, gamma_w, gamma_g, P_b, R_sb, R_s, beta_o, beta_w, 
                z, beta_g, mu_o, mu_w, mu_g, sigma_o, sigma_w, c_o, TDS, d, 
                theta, delta, EOT, nodes, outlet_depth, depth=depth, 
                pressure=pressure
            ) 
        #update when converged
        P_inlet = P_calc
        P_outlet = P_inlet + dP_dz * dz 
        T_inlet = T_outlet
        depth.append(outlet_depth)
        pressure.append(P_inlet)

        if outlet_depth == EOT:
            return depth, pressure
    
    return pressure_traverse(
        P_inlet, P_outlet, T_inlet, dT_dz, q_o_sc, q_w_sc, q_g_sc, R_p, API, 
        gamma_w, gamma_g, P_b, R_sb, R_s, beta_o, beta_w, z, beta_g, mu_o, 
        mu_w, mu_g, sigma_o, sigma_w, c_o, TDS, d, theta, delta, EOT, nodes, 
        outlet_depth, depth=depth, pressure=pressure
    ) 

