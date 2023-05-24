"""
Beggs-Brill Parameters

@author: Daniel Mendoza
@email: dmendozac90@gmail.com
"""

from . import utilities
from numpy import log, sin

G = 32.17 #ft/sec**2
G_c = 32.17 #lbm ft/lbf sec**2

def froude_number(v_m, d):
    """
    Froude number computation.

    Parameters
    ----------
    v_m : mixture velocity (ft/s)
    d : pipe inner diameter (in)

    Returns
    -------
    N_fr : Froude number (unitless)
    """
    d = utilities.length_inches_to_ft(d)
    N_fr = v_m ** 2 / G / d

    return N_fr

def mod_flow_pattern_param_1(lambda_l):
    """
    L_1 computation.

    Parameters
    ----------
    lambda_l : no-slip liquid volume fraction (unitless)

    Returns 
    -------
    L_1 : modified flow parameter (unitless)
    """

    L_1 = 316 * lambda_l ** 0.302

    return L_1

def mod_flow_pattern_param_2(lambda_l):
    """
    L_2 computation.

    Parameters
    ----------
    lambda_l : no-slip liquid volume fraction (unitless)

    Returns 
    -------
    L_2 : modified flow parameter (unitless)
    """

    L_2 = 9.25e-4 * lambda_l ** -2.468

    return L_2

def mod_flow_pattern_param_3(lambda_l):
    """
    L_3 computation.

    Parameters
    ----------
    lambda_l : no-slip liquid volume fraction (unitless)

    Returns 
    -------
    L_3 : modified flow parameter (unitless)
    """

    L_3 = 0.1 * lambda_l ** -1.452

    return L_3

def mod_flow_pattern_param_4(lambda_l):
    """
    L_4 computation.

    Parameters
    ----------
    lambda_l : no-slip liquid volume fraction (unitless)

    Returns 
    -------
    L_4 : modified flow parameter (unitless)
    """

    L_4 = 0.5 * lambda_l ** -6.738

    return L_4

def flow_pattern_prediction(N_fr, lambda_l):
    """
    Flow pattern prediction. This function returns a tuple that indicates the 
    preicted flow pattern as follows:

    (segregated, transition, intermittent, distributed)

    Parameters
    ----------
    N_fr : Fround number (unitless)
    lambda_l : no-slip liquid volume fraction (unitless)

    Returns 
    -------
    tuple : modified flow parameter prediction 
    """

    L_1 = mod_flow_pattern_param_1(lambda_l)
    L_2 = mod_flow_pattern_param_2(lambda_l)
    L_3 = mod_flow_pattern_param_3(lambda_l)
    L_4 = mod_flow_pattern_param_4(lambda_l)

    #segregated
    if (lambda_l < 0.01 and N_fr < L_1) or (lambda_l >= 0.01 and N_fr < L_2):
        return (True, False, False, False)
    #transition
    elif (lambda_l >= 0.01 and L_2 <= N_fr <= L_3):
        return (False, True, False, False)
    #intermittent
    elif (0.01 <= lambda_l and L_3 < N_fr <= L_1) or (lambda_l >= 0.4 and L_3 < N_fr <= L_4):
        return (False, False, True, False)
    #distributed
    elif (lambda_l < 0.4 and N_fr >= L_1) or (lambda_l >= 0.4 and N_fr > L_4):
        return (False, False, False, True)

def horizontal_liquid_holdup(N_fr, lambda_l):
    """
    Horizontal liquid holdup computation. This function returns a tuple to 
    handle the transition flow pattern in which the horizontal liquid holdup 
    is a function of the segregated and intermittent horizontal liquid holdup 
    values. 

    Parameters
    ----------
    N_fr : Fround number (unitless)
    lambda_l : no-slip liquid volume fraction (unitless)

    Returns
    -------
    H_l_0 : horizontal liquid holdup (unitless)
    """

    segregated, transition, intermittent, distributed = \
        flow_pattern_prediction(N_fr, lambda_l)
    
    if segregated:
        #horizontal liquid holdup coefficients
        a, b, c = 0.98, 0.4846, 0.0868
        
        H_l_0 = a * lambda_l ** b / (N_fr ** c)

        return H_l_0, 0

    elif intermittent:
        #horizontal liquid holdup coefficients
        a, b, c = 0.845, 0.5351, 0.0173

        H_l_0 = a * lambda_l ** b / (N_fr ** c)

        return H_l_0, 0

    elif distributed:
        #horizontal liquid holdup coefficients
        a, b, c = 1.065, 0.5824, 0.0609

        H_l_0 = a * lambda_l ** b / (N_fr ** c)

        return H_l_0, 0
    
    elif transition:
        #horizontal liquid holdup coefficients
        a_seg, b_seg, c_seg = 0.98, 0.4846, 0.0868
        a_int, b_int, c_int = 0.845, 0.5351, 0.0173

        H_l_0_seg = a_seg * lambda_l ** b_seg / (N_fr ** c_seg)
        H_l_0_int = a_int * lambda_l ** b_int / (N_fr ** c_int)

        return H_l_0_seg, H_l_0_int
    
def corrected_liquid_holdup_uphill(N_fr, N_Lv, lambda_l, theta):
    """
    Corrected horizontal liquid holdup computation for uphill conditions. 

    Parameters
    ----------
    N_fr : Fround number (unitless)
    lambda_l : no-slip liquid volume fraction (unitless)
    
    Returns
    -------
    H_l_theta : corrected liquid holdup (unitless)
    """

    theta = utilities.angle_degrees_to_radians(theta)

    segregated, transition, intermittent, distributed = \
        flow_pattern_prediction(N_fr, lambda_l)
    
    A = 1

    if segregated:
        #C-coefficient coefficients
        e, f, g, h = 0.011, -3.768, 3.539, -1.614
        C = (1 - lambda_l) * log(e * lambda_l ** f * N_Lv ** g * N_fr ** h)
        if C < 0:
            C = 0
        psi = 1 + C * (sin(1.8 * theta) - 0.333 * sin(1.8 * theta) ** 3)

    elif intermittent:
        e, f, g, h = 2.96, 0.305, -0.4473, 0.0978
        C = (1 - lambda_l) * log(e * lambda_l ** f * N_Lv ** g * N_fr ** h)
        if C < 0:
            C = 0
        psi = 1 + C * (sin(1.8 * theta) - 0.333 * sin(1.8 * theta) ** 3)
    
    elif distributed:
        C, psi = 0, 1

    elif transition:
        L_2 = mod_flow_pattern_param_2(lambda_l)
        L_3 = mod_flow_pattern_param_3(lambda_l)
        A = (L_3 - N_fr) / (L_3 - L_2)

    H_l_0_1, H_l_0_2 = horizontal_liquid_holdup(N_fr, lambda_l)
    H_l_theta = (A * H_l_0_1 + (1 - A) * H_l_0_2) * psi
    #Payne modification to adjust for underprection of friction factors
    H_l_theta *= 0.924 
    #undo modification if the corrected liquid holdup is less than the 
    #no-slip liquid volume fraction
    if H_l_theta < lambda_l:
        H_l_theta /= 0.924

    return H_l_theta

def corrected_liquid_holdup_downhill(N_fr, N_Lv, lambda_l, theta):
    """
    Corrected horizontal liquid holdup computation for downhill conditions. 

    Parameters
    ----------
    N_fr : Fround number (unitless)
    lambda_l : no-slip liquid volume fraction (unitless)
    
    Returns
    -------
    H_l_theta : corrected liquid holdup (unitless)
    """

    _, transition, _, _ = \
        flow_pattern_prediction(N_fr, lambda_l)
    
    A = 1
    
    e, f, g, h = 4.7, -0.3692, 0.1244, -0.5056
    C = (1 - lambda_l) * log(e * lambda_l ** f * N_Lv ** g * N_fr ** h)
    psi = 1 + C * (sin(1.8 * theta) - 0.333 * sin(1.8 * theta) ** 3)
    
    if transition:
        L_2 = mod_flow_pattern_param_2(lambda_l)
        L_3 = mod_flow_pattern_param_3(lambda_l)
        A = (L_3 - N_fr) / (L_3 - L_2)

    H_l_0_1, H_l_0_2 = horizontal_liquid_holdup(N_fr, lambda_l)
    H_l_theta = (A * H_l_0_1 + (1 - A) * H_l_0_2) * psi
    #Payne modification to adjust for underprection of friction factors
    H_l_theta *= 0.685

    return H_l_theta
    