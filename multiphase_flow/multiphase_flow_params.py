"""
Multiphase Flow Parameters

@author: Daniel Mendoza
@email: dmendozac90@gmail.com
"""

from . import utilities
from numpy import log, log10, exp, sin

G = 32.17 #ft/sec**2
G_c = 32.17 #lbm ft/lbf sec**2
RHO_W = 62.37 #lbm/ft**3

def insitu_q_o(q_o_sc, beta_o):
    """
    Insitu oil rate computation. Volume is adjusted for pressure and 
    temperature effects via the oil formation volume factor. The conversion 
    from days to seconds and barrels to cubic feet can often be encountred 
    as a single constant equivalent to 6.498356e-05.

    Parameters
    ----------
    q_os_c : surface oil rate (STB/D)
    beta_o : oil formation volume factor (ft**3/STB)

    Returns
    -------
    q_o : insitu oil rate (ft**3/sec)
    """

    #convert barrels to cubic feet
    q_o_sc = utilities.volume_bbl_to_cubic_ft(q_o_sc)
    #convert daily rate to seconds
    q_o_sc = utilities.time_days_to_sec(q_o_sc)

    q_o = q_o_sc * beta_o

    return q_o

def insitu_q_w(q_w_sc, beta_w):
    """
    Insitu water rate computation. Volume is adjusted for pressure and 
    temperature effects via the water formation volume factor. The conversion 
    from days to seconds and barrels to cubic feet can often be encountred 
    as a single constant equivalent to 6.498356e-05.

    Parameters
    ----------
    q_w_sc : surface water rate (STB/D)
    beta_w : water formation volume factor (ft**3/STB)

    Returns
    -------
    q_w : insitu water rate (ft**3/sec)
    """

    #convert barrels to cubic feet
    q_w_sc = utilities.volume_bbl_to_cubic_ft(q_w_sc)
    #convert daily rate to seconds
    q_w_sc = utilities.time_days_to_sec(q_w_sc)

    q_w = q_w_sc * beta_w

    return q_w

def insitu_q_g(q_g_sc, q_l_sc, R_s, beta_g): 
    """
    Insitu gas rate computation. This quantity is the total surface gas less 
    the quantity evolved into solution from associated pressure temperature 
    effects. The conversion from days to seconds coupled with the leading 
    coefficient of the gas formation volume factor can be encountered as a 
    single constant equivalent to 3.263889e-07.

    Parameters
    ----------
    q_g_sc : surface gas rate (Mscf/D)
    q_l_sc : surface liquid rate (STB/D)
    R_s : solution gas-oil ratio (scf/STB)
    beta_g : (ft**3/scf)

    Returns
    -------
    q_g : insitu gas rate (ft**3/sec)
    """

    #convert Mscf to scf
    q_g_sc = utilities.volume_mscf_to_scf(q_g_sc)
    #insitu gas rate is surface gas rate less gas evolved back into solution
    q_g = (q_g_sc - q_l_sc * R_s) * beta_g
    #convert to from days to seconds
    q_g = utilities.time_days_to_sec(q_g)

    return q_g

def insitu_q_l(q_o, q_w):
    """
    Insitu liquid rate computation.

    Parameters
    ----------
    q_o : insitu oil rate (ft**3/sec)
    q_w : insitu water rate (ft**3/sec)

    Returns
    -------
    q_l : insitu liquid rate (ft**3/sec)
    """
    
    q_l = q_o + q_w

    return q_l

def insitu_oil_mass_flux(q_o_sc, API, gamma_g, R_s):
    """
    Insitu oil mass flux computation. This parameter is required to determine 
    the insitu density.

    Parameters
    ----------
    q_o_sc : surface oil rate (STB/D)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns 
    -------
    w_o : oil mass flux rate (lbm/s)
    """
    
    #convert density from lbm/ft**3 to lbm/STB
    rho_w_sc = utilities.density_per_cubic_ft_to_per_bbl(RHO_W)
    gamma_o = utilities.density_api_to_specific_gravity(API)
    #convert surface oil rate from STB/D to STB/sec
    q_o_sc = utilities.time_days_to_sec(q_o_sc)

    w_o = rho_w_sc * gamma_o * q_o_sc + R_s * q_o_sc * \
        utilities.air_density() * gamma_g

    return w_o

def insitu_water_mass_flux(q_w_sc, gamma_w):
    """
    Insitu water mass flux computation. This parameter is required to 
    determine the insitu density.

    Parameters
    ----------
    q_w_sc : surface water rate (STB/D)
    gamma_w : water specific gravity (unitless)

    Returns
    --------
    w_w : water mass flux rate (lmb/s)
    """ 

    #convert density from lbm/ft**3 to lbm/STB
    rho_w_sc = utilities.density_per_cubic_ft_to_per_bbl(RHO_W)
    #convert surface oil rate from STB/D to STB/sec
    q_w_sc = utilities.time_days_to_sec(q_w_sc)

    w_w = rho_w_sc * gamma_w * q_w_sc

    return w_w

def insitu_gas_mass_flux(R_p, R_s, gamma_g, q_o_sc):
    """
    Insitu gas mass flux computation. This parameter is required to determine 
    the insitu density.

    Parameters
    ----------
    R_p : producing gas-oil ratio (scf/STB)
    R_s : solution gas-oil ratio (scf/STB)
    gamma_g : gas specific gravity (unitless)
    q_o_sc : surface oil rate (STB/D)

    Returns
    --------
    w_g : gas mass flux rate (lmb/s)
    """

    #convert density from lbm/ft**3 to lbm/STB
    q_o_sc = utilities.time_days_to_sec(q_o_sc)
    #The difference in producing and solution gas-oil ratio corresponds to 
    #the total surface gas less the quantity that evolves into solution as 
    #pressure and temperature change
    w_g = (R_p - R_s) * utilities.air_density() * gamma_g * q_o_sc
    
    if w_g < 0:
        w_g = 0.

    return w_g

def insitu_oil_density(w_o, q_o):
    """
    Insitu oil density computation.

    Parameters
    ----------
    w_o : oil mass flux rate (lbm/s)
    q_o : insitu oil rate (ft**3/s)

    Returns
    -------
    rho_o : insitu oil density (lbm/ft**3)
    """

    rho_o = w_o / q_o

    return rho_o

def insitu_water_density(w_w, q_w):
    """
    Insitu water density computation.

    Parameters
    ----------
    w_w : water mass flux rate (lbm/s)
    q_w : insitu water rate (ft**3/s)

    Returns
    -------
    rho_w : insitu water density (lbm/ft**3)
    """
    
    if q_w != 0:
        rho_w = w_w / q_w
    #case where there is no water production
    else:
        rho_w = 0.

    return rho_w

def insitu_gas_density(w_g, q_g):
    """
    Insitu gas density computation.

    Parameters
    ----------
    w_g : gas mass flux rate (lbm/s)
    q_g : insitu gas rate (ft**3/s)

    Returns
    -------
    rho_g : insitu gas density (lbm/ft**3)
    """

    if q_g == 0.:
        rho_g = 0.
    else:
        rho_g = w_g / q_g

    return rho_g

def insitu_liquid_density(rho_o, rho_w, f_o, f_w):
    """
    Insitu liquid density computation.

    Parameters
    ----------
    rho_o : insitu oil density (lbm/ft**3)
    rho_w : insitu water density (lbm/ft**3)
    f_o : oil cut (unitless)
    f_w : water cut (unitless)

    Returns
    -------
    rho_l : insitu liquid density (lbm/ft**3)
    """

    rho_l = rho_o * f_o + rho_w * f_w

    return rho_l

def no_slip_liquid_volume_fraction(q_l, q_g):
    """
    Liquid volume fraction computation.

    Parameters
    ----------
    q_l : insitu liquid rate (ft**3/s)
    q_g : insitu gas rate (ft**3/s)

    Returns
    -------
    lambda_l : no-slip liquid volume fraction (unitless)
    """

    lambda_l = q_l / (q_l + q_g)

    return lambda_l

def oil_cut(q_o, q_w):
    """
    Oil fraction.

    Parameters
    ----------
    q_o : insitu oil rate (ft**3/s)
    q_w : insitu water rate (ft**3/s)

    Returns
    -------
    f_o : oil fraction (unitless)
    """

    f_o = q_o / (q_o + q_w)
    
    return f_o

def water_cut(q_o, q_w):
    """
    Water fraction.

    Parameters
    ----------
    q_o : insitu oil rate (ft**3/s)
    q_w : insitu water rate (ft**3/s)

    Returns
    -------
    f_w : water fraction (unitless)
    """

    f_w = q_w / (q_o + q_w)

    return f_w

def superficial_liquid_velocity(q_l, d):
    """
    Superficial liquid velocity computation. This parameters assumes that the 
    liquid phase occupies the entire pipe area.

    Parameters
    ----------
    q_l : insitu liquid rate (ft**3/s)
    d : pipe inner diameter (in)

    Returns
    -------
    v_sl : superficial liquid velocity (ft/sec)
    """

    v_sl = q_l / utilities.tubing_cross_sectional_area(d)

    return v_sl

def superficial_gas_velocity(q_g, d):
    """
    Superficial gas velocity computation. This parameters assumes that the 
    gas phase occupies the entire pipe area.

    Parameters
    ----------
    q_g : insitu gas rate (ft**3/s)
    d : pipe inner diameter (in)

    Returns
    -------
    v_sg : superficial gas velocity (ft/s)
    """

    v_sg = q_g / utilities.tubing_cross_sectional_area(d)

    return v_sg

def mixture_velocity(v_sl, v_sg):
    """
    Mixture velocity computation.

    Parameters
    ----------
    v_sl : superficial liquid velocity (ft/s)
    v_sg : superficial gas velocity (ft/s)

    Returns
    -------
    v_m : superficial mixture velocity (ft/s)
    """

    v_m = v_sl + v_sg

    return v_m

def insitu_liquid_viscosity(mu_o, mu_w, f_o, f_w):
    """
    Insitu liquid viscosity computation.

    Parameters
    ----------
    mu_o : insitu oil viscosity (cp)
    mu_w : insitu water viscosity (cp)
    f_o : oil cut (unitless)
    f_w : water cut (unitless)

    Returns
    -------
    mu_l : insitu liquid viscosity (cp)
    """

    mu_l = mu_o * f_o + mu_w * f_w

    return mu_l

def insitu_no_slip_mixture_viscosity(mu_l, mu_g, lambda_l):
    """
    Insitu no-slip mixture viscosity computation. 

    Parameters
    ----------
    mu_l : insitu liquid viscosity (cp)
    mu_g : insitu gas viscosity (cp)
    lambda_l : no-slip liquid volume fraction (unitless)

    Returns
    -------
    mu_n : insitu no-slip mixture viscosity (cp)
    """

    mu_n = mu_l * lambda_l + mu_g * (1 - lambda_l)

    return mu_n

def insitu_slip_mixture_viscosity(mu_l, mu_g, H_l):
    """
    Insitu slip mixture viscosity computation. 

    Parameters
    ----------
    mu_l : insitu liquid viscosity (cp)
    mu_g : insitu gas viscosity (cp)
    H_l : liquid holdup (unitless)

    Returns
    -------
    mu_s : insitu slip mixture viscosity (cp)
    """

    mu_s = mu_l * H_l + mu_g * (1 - H_l)

    return mu_s

def insitu_no_slip_mixture_density(rho_l, rho_g, lambda_l):
    """
    Insitu no-slip mixture density computation.

    Parameters
    ----------
    rho_l : insitu liquid density (lbm/ft**3)
    rho_g : insitu gas density (lbm/ft**3)
    lambda_l : no-slip liquid volume fraction (unitless)

    Returns
    -------
    rho_n : insitu no-slip mixture density (lbm/ft**3)
    """

    rho_n = rho_l * lambda_l + rho_g * (1 - lambda_l)

    return rho_n

def insitu_slip_mixture_density(rho_l, rho_g, H_l):
    """
    Insitu slip mixture density computation. 

    Parameters
    ----------
    rho_l : insitu liquid density (lbm/ft**3)
    rho_g : insitu gas density (lbm/ft**3)
    H_l : liquid holdup (unitless)

    Returns
    -------
    rho_s : insitu slip mixture density (lbm/ft**3)
    """

    rho_s = rho_l * H_l + rho_g * (1 - H_l)

    return rho_s

def liquid_surface_tension(sigma_o, sigma_w, f_o, f_w):
    """
    Insitu liquid surface tension computation.

    Parameters
    ----------
    sigma_o : oil surface tension (dyne/cm)
    sigma_w : water surface tension (dyne/cm)
    f_o : oil cut (unitless)
    f_w : water cut (unitless)

    Returns
    -------
    sigma_l : liquid surface tension(dyne)
    """

    sigma_l = sigma_o * f_o + sigma_w * f_w

    return sigma_l

def liquid_velocity_number(v_sl, rho_l, sigma_l):
    """
    Liquid velocity number computation.

    Parameters
    ----------
    v_sl : superficial liquid velocity (ft/s)
    rho_l : insitu liquid density (lbm/ft**3)
    sigma_l : insitu liquid surface tension (dyne/cm)

    Returns
    -------
    N_Lv : liquid velocity number (unitless)
    """

    N_Lv = 1.938 * v_sl * (rho_l / sigma_l) ** 0.25

    return N_Lv

def gas_velocity_number(v_sg, rho_l, sigma_l):
    """
    Gas velocity number computation.

    Parameters
    ----------
    v_sl : superficial liquid velocity (ft/s)
    rho_l : insitu liquid density (lbm/ft**3)
    sigma_l : insitu liquid surface tension (dyne/cm)

    Returns
    -------
    N_gv : gas velocity number (unitless)
    """

    N_Lv = 1.938 * v_sg * (rho_l / sigma_l) ** 0.25

    return N_Lv

def pipe_diameter_number(rho_l, sigma_l, d):
    """
    Pipe diameter number computation.

    Parameters
    ----------
    rho_l : insitu liquid density (lbm/ft**3)
    sigma_l : insitu liquid surface tension (dyne/cm)
    d : pipe inner diameter (in)

    Returns
    -------
    N_d : pipe diameter number (unitless)
    """

    d = utilities.length_inches_to_ft(d)

    N_d = 120.872 * d * (rho_l / sigma_l) ** 0.5

    return N_d

def liquid_viscosity_number(mu_l, rho_l, sigma_l):
    """
    Liquid viscosity number computation.

    Parameters
    ----------
    mu_l : liquid viscosity (cp)
    rho_l : liquid density (cp)
    sigma_l : liquid surface tension (dyne/cm)

    Returns 
    -------
    N_L : liquid viscosity number (unitless)
    """

    N_L = 0.15726 *  mu_l * (1 / rho_l / sigma_l ** 3) ** 0.25 

    return N_L

def reynolds_number(rho, v, d, mu):
    """
    Reynolds number computation. The Reynolds number is defined as follows:

    N_re = rho * v * d / mu

    When using english units one must ensure that the units cancel to yield a 
    dimensionless number. Under the specified parameter units, this equation  
    requires the unit conversion constant g_c (lbm*ft)/(lbf*s**2) thus:

    N_re = rho * v * d / (g_c * mu)

    Additionally, because the units of centipoise are defined in metric 
    units, one must use a conversion constant to yield a unitless parameter. 
    This can be derived as follows:

    10 p = 1 Pa * s 
    1 cp = 0.001 Pa * s
    1 cp = 0.001 (N/m**2) * s

    The following conversion constants can then be substituted

    1 lbf = 4.448 N
    1 ft = 0.3048 m

    to yield

    1 cp = 0.001 * (N/m**2) * (0.3048 m / 1 ft) ** 2 * (1 lbf/4.448 N) * s

    1 cp = 0.001 * (0.3048) ** 2 (lbf/ft**2) * (1 / 4.448) * s

    1 cp = 0.001 * 0.020886475 (lbf/ft**2) * s

    1 cp = 2.0886475e-05 (lbf/ft**2) * s
    
    Parameters
    ----------
    rho : density (lbm/ft**3)
    v : velocity (ft/sec)
    d : pipe inner diameter (in)
    mu : viscosiy (cp)

    Returns
    -------
    N_re : Reynolds number (unitless)
    """

    #convert cp to (lbf / ft**2) * s
    mu = utilities.viscosity_cp_to_field_units(mu)
    #convert inches to ft
    d = utilities.length_inches_to_ft(d)

    N_re = rho * v * d / mu / G_c

    return N_re

def friction_factor(f, N_re, delta, d, tol=1e-6):
    """
    This function computes the Moody friction factor using the modified 
    Colebrook-White equation. The non-modified Colebrook-White equation 
    yields the Fanning friction factor which is 4 times less than the  
    Moody friction factor.

    Parameters
    ----------
    f : initial guess of Moody friction factor (unitless)
    N_re : Reynolds number (unitless)
    delta : tubing roughness (in)
    d : pipe inner diameter (in)
    tol : error tolerance 

    Returns
    -------
    f : Moody friction factor (unitless)
    """

    if N_re < 2000:
        f = 64 / N_re
    
    else:
        epsilon = delta / utilities.length_inches_to_ft(d)
        lhs = f ** -0.5
        rhs = -2 * log10(epsilon / 3.7 + 2.51 / (N_re * f ** 0.5))
        error = abs(lhs - rhs)
        while error > tol:
            f = rhs ** -2
            return friction_factor(f, N_re, delta, d)
        
    return f

def parameter_y(lambda_l, H_l):
    """
    Parameter y computation. This parameter is required for the two-phase 
    friction factor.

    Parameters
    ----------
    lambda_l : no-slip liquid volume fraction (unitless)
    H_l : liquid holdup (unitless)

    Returns
    -------
    y : parameter y (unitless)
    """

    y = lambda_l / H_l ** 2

    return y

def parameter_s(y):
    """
    Parameter s computation. This parameter is required for the two-phase 
    friction factor.

    Parameters
    ----------
    y  : parameter y (unitless)

    Returns
    -------
    s : parameter s (unitless)
    """
    if y > 1 and y < 1.2:
        s = log(2.2 * y -1.2)
    
    else:
        s = log(y) / (-0.0532 + 3.182 * log(y) - 0.8725 * \
                         log(y) ** 2 + 0.01853 * log(y) ** 4)
    
    return s

def two_phase_friction_ratio(s):
    """
    Two-phase friction ratio calculation.

    Parameters
    ----------
    s : parameter s (unitless)

    Returns
    -------
    f_fn : two-phase friction factor ratio (unitless)
    """

    f_fn = exp(s)

    return f_fn

def two_phase_friction_factor(f, f_fn):
    """
    Two-phase friction calculation.

    Parameters
    ----------
    f : friction factor (unitless)
    f_fn : two-phase friction ratio (unitless)

    Returns
    -------
    f_two_phase : two-phase friction factor (unitless)
    """

    f_two_phase = f * f_fn

    return f_two_phase

def friction_pressure_loss_gradient(rho_n, v_m, f, d):
    """
    Friction pressure loss calculation.

    Parameters
    ----------
    rho_n : insitu no-slip mixture density (lbm/ft**3)
    v_m : superficial mixture velocity (ft/s)
    f : two-phase friction factor (unitless)
    d : pipe inner diameter (in)

    Returns
    ------- 
    dp_dz_f : friction pressure loss gradient (lbf/ft**3)
    """
    
    d = utilities.length_inches_to_ft(d)

    dp_dz_f = (f * rho_n * v_m ** 2) / (2 *  d * G_c)

    return dp_dz_f

def elevation_pressure_loss_gradient(rho_s, theta):
    """
    Elevation pressure loss calculation.

    Parameters
    ----------
    rho_s : slip mixture density (lbm/ft**3)

    Returns
    -------
    dp_dz_ele : elevation pressure loss gradient (lbf/ft**3)
    """

    theta = utilities.angle_degrees_to_radians(theta)

    #G/G_c is unity; however, the units do not cancel out leaving (lbf/lbm) 
    #and thus converting density from (lbm/ft**3) to (lbf/ft**3)
    dp_dz_ele = rho_s * G / G_c * sin(theta)

    return dp_dz_ele

def pressure_loss_gradient(rho_n, rho_s, v_m, f, d, theta):
    """
    Pressure loss gradient calculation.

    Parameters
    ----------
    rho_n : insitu no-slip mixture density (lbm/ft**3)
    rho_s : slip mixture density (lbm/ft**3)
    v_m : superficial mixture velocity (ft/s)
    f : two-phase friction factor (unitless)
    d : pipe inner diameter (in)

    Returns
    -------
    dp_dz : pressure loss gradient (psi/ft)
    """

    dp_dz_f = friction_pressure_loss_gradient(rho_n, v_m, f, d)
    dp_dz_ele = elevation_pressure_loss_gradient(rho_s, theta)
    dp_dz = dp_dz_f + dp_dz_ele
    dp_dz = utilities.density_to_psi_per_foot(dp_dz)

    return dp_dz