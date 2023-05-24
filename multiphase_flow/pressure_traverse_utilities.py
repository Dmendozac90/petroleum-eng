"""
Pressure Traverse Helper Functions

@author: Daniel Mendoza
@email: dmendozac90@gmail.com
"""

from . import reservoir_correlations
from . import multiphase_flow_params
from . import beggs_brill_params


def saturated_reservoir_params(
        P, T, API, gamma_g, TDS, P_sep, T_sep, R_s, beta_o, beta_w, z, 
        beta_g, mu_o, mu_w, mu_g
    ):
    """
    Helper function that returns saturated fluid properties for pressure 
    traverse computation. 
    """

    if R_s == 'standing':
        R_s = reservoir_correlations.standing_r_s(
            P, T, API, gamma_g
        )
    elif R_s == 'vasquez-beggs':
        R_s= reservoir_correlations.vasquez_beggs_r_s(
            P, T, P_sep, T_sep, API, gamma_g
        )
    elif R_s == 'almarhoun':
        R_s = reservoir_correlations.al_marhoun_r_s(
            P, T, API, gamma_g
        )
    elif R_s == 'glaso':
        R_s = reservoir_correlations.glaso_r_s(
            P, T, API, gamma_g
        )
    elif R_s == 'petrosky-farshad':
        R_s = reservoir_correlations.petrosky_farshad_r_s(
            P, T, API, gamma_g
        )

    if beta_o == 'standing':
        beta_o = reservoir_correlations.standing_beta_o(
            T, API, gamma_g, R_s
        )
    elif beta_o == 'vasquez-beggs':
        beta_o = reservoir_correlations.vasquez_beggs_beta_o(
            T, P_sep, T_sep, API, gamma_g, R_s
        )
    elif beta_o == 'almarhoun':
        beta_o =reservoir_correlations.al_marhoun_beta_o(
            T, API, gamma_g, R_s
        )
    elif beta_o == 'glaso':
        beta_o = reservoir_correlations.glaso_beta_o(
            T, API, gamma_g, R_s
        )  
    elif beta_o == 'petrosky-farshad':
        beta_o = reservoir_correlations.petrosky_farshad_beta_o(
            T, API, gamma_g, R_s
        )

    if beta_w == 'analytic':
        beta_w = reservoir_correlations.beta_w(
            P, T
        )
    
    if z == 'analytic':
        z = reservoir_correlations.hall_yarb_compressibility(
            P, T, gamma_g
        )

    if beta_g == 'analytic':
        beta_g = reservoir_correlations.beta_g(
            P, T, z
        )

    if mu_o == 'standing':
        mu_o = reservoir_correlations.standing_mu_ob(
            T, API, R_s
        )
    elif mu_o == 'beggs-robinson':
        mu_o = reservoir_correlations.beggs_robinson_mu_ob(
            T, API, R_s
        )

    if mu_w == 'matthews-russel':
        mu_w = reservoir_correlations.matthews_russel_mu_w(
            P, T, TDS
        ) 

    if mu_g == 'carr':
        mu_g = reservoir_correlations.carr_gas_viscosity(
            P, T, gamma_g
        )

    return R_s, beta_o, beta_w, z, beta_g, mu_o, mu_w, mu_g

def undersaturated_reservoir_params(
        P, T, API, gamma_g, TDS, P_sep, T_sep, R_s, P_b, beta_o, beta_w, mu_o, mu_w, c_o
    ):
    """
    Helper function that returns undersaturated fluid properties for pressure 
    traverse computation. 
    """

    R_s = R_s

    if c_o == 'spivey':
        c_o = reservoir_correlations.spivey_oil_compressibility(
            P, P_b, T, API, gamma_g, R_s
        )
    elif c_o == 'vasquez-beggs':
        c_o = \
            reservoir_correlations.vasquez_beggs_oil_compressibility(
            P, T, P_sep, T_sep, API, gamma_g, R_s
            )
    elif c_o == 'petrosky-farshad':
        c_o = reservoir_correlations.petrosky_farshad_oil_compressibility(
            P, T, API, gamma_g, R_s
        )
    
    if beta_o == 'standing':
        beta_ob = reservoir_correlations.standing_beta_o(
            T, API, gamma_g, R_s
        )
    elif beta_o == 'vasquez-beggs':
        beta_ob = reservoir_correlations.vasquez_beggs_beta_o(
            T, P_sep, T_sep, API, gamma_g, R_s
        )
    elif beta_o == 'almarhoun':
        beta_ob =reservoir_correlations.al_marhoun_beta_o(
            T, API, gamma_g, R_s
        )
    elif beta_o == 'glaso':
        beta_ob = reservoir_correlations.glaso_beta_o(
            T, API, gamma_g, R_s
        )  
    elif beta_o == 'petrosky-farshad':
        beta_ob = reservoir_correlations.petrosky_farshad_beta_o(
            T, API, gamma_g, R_s
        )
    beta_o = reservoir_correlations.undersaturated_beta_o(
                P, P_b, beta_ob, c_o
            )
    
    if beta_w == 'analytic':
        beta_w = reservoir_correlations.beta_w(
            P, T
        )

    beta_g = 0.

    if type(mu_o) == str:
        mu_o = reservoir_correlations.bergman_sutton_mu_o(
            P, P_b, T, API, R_s
        )

    if mu_w == 'matthews-russel':
        mu_w = reservoir_correlations.matthews_russel_mu_w(
            P, T, TDS
        ) 

    mu_g = 0.

    return R_s, beta_o, beta_w, beta_g, mu_o, mu_w, mu_g

def volumetric_rates(q_o_sc, q_w_sc, q_g_sc, beta_o, beta_w, beta_g, R_s):
    """
    
    """
    q_l_sc = q_o_sc + q_w_sc
    q_o = multiphase_flow_params.insitu_q_o(q_o_sc, beta_o)
    q_w = multiphase_flow_params.insitu_q_w(q_w_sc, beta_w)
    q_l = multiphase_flow_params.insitu_q_l(q_o, q_w)
    q_g = multiphase_flow_params.insitu_q_g(q_g_sc, q_l_sc, R_s, beta_g)

    return q_o, q_w, q_l, q_g

def mass_flux_rates(q_o_sc, q_w_sc, API, gamma_w, gamma_g, R_s, R_p):
    """
    
    """

    w_o = multiphase_flow_params.insitu_oil_mass_flux(
        q_o_sc, API, gamma_g, R_s
    )

    w_w = multiphase_flow_params.insitu_water_mass_flux(
        q_w_sc, gamma_w
    )

    w_g = multiphase_flow_params.insitu_gas_mass_flux(
        R_p, R_s, gamma_g, q_o_sc
    )

    return w_o, w_w, w_g

def insitu_densities(
        q_o, q_w, q_g, q_o_sc, q_w_sc, API, gamma_w, gamma_g, R_s, R_p
):
    """
    
    """
    w_o, w_w, w_g = mass_flux_rates(
        q_o_sc, q_w_sc, API, gamma_w, gamma_g, R_s, R_p
    )
    rho_o = multiphase_flow_params.insitu_oil_density(
        w_o, q_o
    )
    rho_w = multiphase_flow_params.insitu_water_density(
        w_w, q_w
    )
    rho_g = multiphase_flow_params.insitu_gas_density(
        w_g, q_g
    )

    return rho_o, rho_w, rho_g

def fluid_fractions(q_l, q_o, q_w, q_g):
    """
    
    """
    
    lambda_l = multiphase_flow_params.no_slip_liquid_volume_fraction(
        q_l, q_g
    )
    f_o = multiphase_flow_params.oil_cut(
        q_o, q_w
    )
    f_w = multiphase_flow_params.water_cut(
        q_o, q_w
    )

    return lambda_l, f_o, f_w

def fluid_velocities(q_l, q_g, d):
    """
    
    """

    v_sl = multiphase_flow_params.superficial_liquid_velocity(
        q_l, d
    )
    v_sg = multiphase_flow_params.superficial_gas_velocity(
        q_g, d
    )
    v_m = multiphase_flow_params.mixture_velocity(
        v_sl, v_sg
    )

    return v_sl, v_sg, v_m

def liquid_params(rho_o, rho_w, mu_o, mu_w, sigma_o, sigma_w, f_o, f_w):
    """
    
    """

    rho_l = multiphase_flow_params.insitu_liquid_density(
        rho_o, rho_w, f_o, f_w
    )
    mu_l = multiphase_flow_params.insitu_liquid_viscosity(
        mu_o, mu_w, f_o, f_w
    )
    sigma_l = multiphase_flow_params.liquid_surface_tension(
        sigma_o, sigma_w, f_o, f_w
    )

    return rho_l, mu_l, sigma_l

def duns_ross_params(rho_l, mu_l, sigma_l, v_sl, v_sg, d):
    """
    
    """

    N_Lv = multiphase_flow_params.liquid_velocity_number(
        v_sl, rho_l, sigma_l
    )
    N_gv = multiphase_flow_params.gas_velocity_number(
        v_sg, rho_l, sigma_l
    )
    N_d = multiphase_flow_params.pipe_diameter_number(
        rho_l, sigma_l, d
    )
    N_L = multiphase_flow_params.liquid_viscosity_number(
        mu_l, rho_l, sigma_l
    )

    return N_Lv, N_gv, N_d, N_L

def no_slip_params(rho_l, rho_g, mu_l, mu_g, lambda_l):
    """
    
    """

    rho_n = multiphase_flow_params.insitu_no_slip_mixture_density(
        rho_l, rho_g, lambda_l
    )
    mu_n = multiphase_flow_params.insitu_no_slip_mixture_viscosity(
        mu_l, mu_g, lambda_l
    )

    return rho_n, mu_n

def beggs_brill_liquid_holdup(N_Lv, v_m, lambda_l, d, theta):
    """
    
    """

    N_fr = beggs_brill_params.froude_number(
        v_m, d
    )
    if theta > 0:
        H_l_theta = beggs_brill_params.corrected_liquid_holdup_uphill(
            N_fr, N_Lv, lambda_l, theta
        )
    elif theta < 0:
        H_l_theta = beggs_brill_params.corrected_liquid_holdup_downhill(
            N_fr, N_Lv, lambda_l, theta
        )
    else:
        H_l_theta = beggs_brill_params.horizontal_liquid_holdup(
            N_fr, lambda_l
        )

    return H_l_theta

def beggs_brill_friction(rho_n, mu_n, lambda_l, v_m, H_l_theta, d, delta):
    """
    
    """

    N_re = multiphase_flow_params.reynolds_number(
        rho_n, v_m, d, mu_n
    )
    y = multiphase_flow_params.parameter_y(
        lambda_l, H_l_theta
    )
    s = multiphase_flow_params.parameter_s(
        y
    )
    f_fn = multiphase_flow_params.two_phase_friction_ratio(
        s
    )
    f = multiphase_flow_params.friction_factor(
        0.01, N_re, delta, d
    )
    f_two_phase = multiphase_flow_params.two_phase_friction_factor(
        f, f_fn
    )

    return f_two_phase