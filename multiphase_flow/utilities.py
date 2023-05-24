"""
Petroleum Engineering Utilites

@author: Daniel Mendoza
@email: dmendozac90@gmail.com
"""

from numpy import pi

G = 32.17 #ft/sec**2
G_c = 32.17 #lbm ft/lbf sec**2
RHO_W_SC = 62.37 #lbm/ft**3
M_AIR = 28.966 #lbm/lb mol
P_SC = 14.696 #psia
R = 10.732 #psia ft**3 / lb-mole °R
T_SC = 60 + 459.67 #°R
FT = 0.3048 #meters
LBF = 4.448 #Newtons  
CP = 0.001 #Pa * sec

def length_inches_to_ft(l):
    """
    Length conversion from inches to feet.
    """

    l /= 12

    return l

def angle_degrees_to_radians(theta):
    """
    Angle converson from degrees to radians.
    """

    theta = theta * pi / 180

    return theta

def area_ft_to_inches(a):
    """
    Area conversion from squared feet to squared inches.
    """

    a /= 144

    return a

def volume_mscf_to_scf(v):
    """
    Volume conversion from Mscf to scf.
    """

    v *= 1000 #scf/Mscf

    return v

def volume_bbl_to_cubic_ft(v):
    """
    Volume conversion from barrels to cubic feet.
    """

    v *= 5.61458

    return v

def density_per_cubic_ft_to_per_bbl(rho):
    """
    Density conversion from per cubic feet to per barrel.
    """

    rho *= 5.61458

    return rho

def density_api_to_specific_gravity(API):
    """
    Density conversion from °API to specific gravity.
    """

    gamma_o = 141.5 / (131.5 + API)

    return gamma_o

def density_to_psi_per_foot(rho):
    """
    Density conversion from lbm/ft**3 to psi/ft. This conversion requires 
    that lbm be converted into lbf using the gravitational acceleration and 
    conversion constants.
    """

    #G/G_c is unity; however, the units do not cancel out leaving (lbf/lbm) 
    #and thus converting density from (lbm/ft**3) to (lbf/ft**3)
    rho = rho * G / G_c / 144

    return rho

def time_days_to_sec(t):
    """
    Time conversion from days to seconds.
    """

    t /= 86400 #seconds/day

    return t

def pressure_to_psi_per_foot(p):
    """
    Pressure conversion from lbf/ft**3 to psi/ft.
    """

    p /= 144

    return p

def viscosity_cp_to_field_units(mu):
    """
    Viscosity conversion from centipoise to lbf / ft**2 * s. This can be 
    derived as follows:

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
    """

    mu = mu * CP * FT ** 2 / LBF 

    return mu

def tubing_cross_sectional_area(d):
    """
    Cross-sectional area computation for cylindrical conduits.

    Parameters
    ----------
    d : inner diameter (in)

    Returns
    -------
    a : cross-sectional area (ft**2)
    """

    d = length_inches_to_ft(d)

    a = pi / 4 * d ** 2

    return a

def air_density(P=P_SC, T=T_SC, M=M_AIR, R=R):
    """
    Air density at the following condtions:

    P = 14.696 psia
    T = 519.67 Rankine
    """

    rho_air = (P_SC * M_AIR) / (R * T_SC)

    return rho_air