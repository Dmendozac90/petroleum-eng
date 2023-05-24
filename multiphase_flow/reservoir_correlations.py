"""
Petroleum Fluids Empirical Correlations

@author: Daniel Mendoza
@email: dmendozac90@gmail.com
"""

from numpy import log10, log, exp, sqrt

def sutton_pseudo_critical_temperature(gamma_g, y_CO2=None, y_H2S=None):
    """
    Sutton pseudo-critical temperature correlation. This equation is valid 
    over the following specific gas gravities: 0.57 < gamma_g < 1.68. The
    correlation uses the correction factors defined by Wichert and Aziz to 
    improve the accuracy in the Standing and Katz compressibility factor when
    significant quantities of acidic gases are present in the gas sample. It 
    is applicable to concentrations of CO2 < 54.4 mol% and H2S < 73.8 mol%. 
    Wichert and Azis reported an average absolute error of 0.97% over the 
    following ranges of data: 154 psia < p < 7026 psia and 40°F < T < 300°F.
    
    Parameters
    ----------
    gamma_g : gas specific gravity (unitless)
    y_CO2 : carbon dioxide mole fraction (unitless)
    y_H2S : hydrogen sulfide mole fraction (unitless)

    Returns
    -------
    T_pc : pseudo-critical temperature (°R)

    References
    ----------
    [1]: Towler, Brian F., *Gas Properties.* General Engineering vol 1., 
    edited by John R. Fanchi, Society of Petroleum Engineers, 2006, 
    pp. 226-229.
    """

    #calculate pseudo-critical temperature
    T_pc = 169.2 + 349.5 * gamma_g - 74 * gamma_g ** 2
    #correct for acidic gases 
    if y_CO2 is not None or y_H2S is not None:
        if y_CO2 is None:
            y_CO2 = 0
        if y_H2S is None:
            y_H2S = 0
        A = y_CO2 + y_H2S
        B = y_H2S
        epsilon = 120 * (A ** 0.9 - A ** 1.6) + 15 * (B ** 0.5 - B ** 4)
        T_pc -= epsilon
    
    return T_pc

def sutton_pseudo_critical_pressure(gamma_g, y_CO2=None, y_H2S=None):
    """
    Sutton pseudo-critical pressure correlation. This equation is valid over 
    the following specific gas gravities: 0.57 < gamma_g < 1.68. The 
    correlation uses the correction factors defined by Wichert and Aziz to 
    improve the accuracy in the Standing and Katz compressibility factor when
    significant quantities of acidic gases are present in the gas sample. It 
    is applicable to concentrations of CO2 < 54.4 mol% and H2S < 73.8 mol%. 
    Wichert and Azis reported an average absolute error of 0.97% over the 
    following ranges of data: 154 psia < p < 7026 psia and 40°F < T < 300°F.
    
    Parameters
    ----------  
    gamma_g : gas specific gravity (unitless)
    y_CO2 : carbon dioxide mole fraction (unitless)
    y_H2S : hydrogen sulfide mole fraction (unitless)

    Returns
    -------
    P_pc : pseudo-critical pressure (psia)

    References
    ----------
    [1]: Towler, Brian F. *Gas Properties.* General Engineering vol 1., 
    edited by John R. Fanchi, Society of Petroleum Engineers, 2006, 
    pp. 226-229.
    """

    #calculate pseudo-critical pressure
    P_pc = 756.8 - 131.07 * gamma_g - 3.6 * gamma_g ** 2
    #correct for acidic gases 
    if y_CO2 is not None or y_H2S is not None:
        if y_CO2 is None:
            y_CO2 = 0
        if y_H2S is None:
            y_H2S = 0
        A = y_CO2 + y_H2S
        B = y_H2S
        epsilon = 120 * (A ** 0.9 - A ** 1.6) + 15 * (B ** 0.5 - B ** 4)
        T_pc = sutton_pseudo_critical_temperature(gamma_g, 0, 0)
        T_pc_adj = sutton_pseudo_critical_temperature(gamma_g, y_CO2, y_H2S)
        P_pc = P_pc * T_pc_adj / (T_pc - y_H2S * (1 - y_H2S) * epsilon)

    return P_pc    

def carr_gas_viscosity(P, T, gamma_g, y_N2=None, y_CO2=None, y_H2S=None):
    """
    Car gas viscosity correlation. 
    
    Parameters
    ----------
    P : pressure (psia)
    T : temperature (°F)
    gamma_g : gas specific gravity (unitless)
    y_N2 : nitrogen mole fraction (unitless)
    y_CO2 : carbon dioxide mole fraction (unitless)
    y_H2S : hydrogen sulfide mole fraction (unitless)

    Returns
    -------
    mu_g : gas viscosity (cp)

    References
    ----------
    [2]: Danesh, Ali, *PVT Behaviour of Petroleum Reservoir Fluids*, 
    Elsevier Science & Technology Books, 1998, pp. 83-85
    """
    
    T += 459.67
    P_pr = P / sutton_pseudo_critical_pressure(gamma_g, y_CO2, y_H2S)
    T_pr = T / sutton_pseudo_critical_temperature(gamma_g, y_CO2, y_H2S)
    
    a_0 = -2.46211820
    a_1 = 2.97054714
    a_2 = -0.286264054
    a_3 = 0.00805420522
    a_4 = 2.80860949
    a_5 = -3.49803305
    a_6 = 0.36037302
    a_7 = -0.0104432413
    a_8 = -0.793385684
    a_9 = 1.39643306
    a_10 = -0.149144925
    a_11 = 0.00441015512
    a_12 = 0.0839387178
    a_13 = -0.186408848
    a_14 = 0.0203367881
    a_15 = -0.000609579263
    
    #compute gas viscosity for impurities at 1 atm. 
    if y_N2 is None:
        mu_1N2 = 0
    else:
        mu_1N2 = (9.59e-3 + 8.48e-3 * log10(gamma_g)) * y_N2
    if y_CO2 is None:
        mu_1CO2 = 0
    else:
        mu_1CO2 = (6.24e-3 + 9.08e-3 * log10(gamma_g)) * y_CO2
    if y_H2S is None:
        mu_1H2S = 0
    else:
        mu_1H2S = (3.73e-3 + 8.49e-3 * log10(gamma_g)) * y_H2S
    
    #compute gas viscosity for hydrocarbon mixture at 1 atm.
    mu_1HC = 8.188e-3 - 6.15e-3 * log10(gamma_g) + \
             (1.709e-5 - 2.062e-6 * gamma_g) * (T - 459.67)
    #sum hydrocarbon and non-hydrocarbon gas viscosities
    mu_1 = mu_1HC + mu_1N2 + mu_1CO2 + mu_1H2S 
    #compute pressure-corrected natural log viscosity ratio 
    mu_r = a_0 + a_1 * P_pr + a_2 * P_pr ** 2 + a_3 * P_pr ** 3 + \
           T_pr * (a_4 + a_5 * P_pr + a_6 * P_pr ** 2 + a_7 * P_pr ** 3) + \
           T_pr ** 2 * (a_8 + a_9 * P_pr + a_10 * P_pr ** 2 + a_11 * \
                        P_pr ** 3) + \
           T_pr ** 3 * (a_12 + a_13 * P_pr + a_14 * P_pr ** 2 + a_15 * \
                        P_pr ** 3)
    #compute temperature-corrected viscosity
    mu_g = (mu_1 / T_pr) * exp(mu_r)

    return mu_g

def hall_yarb_compressibility(
        P, T, gamma_g, y_CO2=None, y_H2S=None, tol=1e-6):
    """
    Hall-Yarborough compressibility factor correlation. The average difference
    between the Standing-Katz correlation chart and this method is -0.158% and
    the average absolute difference is 0.518%. 
    
    Parameters
    ----------
    P : pressure (psia)
    T : temperature (°F)
    gamma_g : gas specific gravity (unitless)
    y_CO2 : carbon dioxide mole fraction (unitless)
    y_H2S : hydrogen sulfide mole fraction (unitless)
    tol : error tolerence (unitless)

    Returns
    -------
    z : gas compressibility factor (unitless)

    References
    ----------
    [3]: Dake, L.P. *Fundamentals of Reservoir Engineering.* 17th ed.,
    Elsevier Science B.V., 1978, pp. 19-20
    """
    
    T = T + 459.67
    P_pr = P / sutton_pseudo_critical_pressure(gamma_g, y_CO2, y_H2S)
    T_pr = T / sutton_pseudo_critical_temperature(gamma_g, y_CO2, y_H2S)
    t_r = T_pr ** -1
    D = 2.18 + 2.82 * t_r
    C = t_r * (90.7 -242.2 * t_r + 42.4 * t_r ** 2)
    B = t_r * (14.76 - 9.76 * t_r + 4.58 * t_r ** 2)
    A = 0.06125 * t_r * exp(-1.2 * (1 - t_r) ** 2)
    #define parameters to initialize loop
    f_Y = 1
    Y = 0.001
    #Newton-Raphson method
    while abs(f_Y) > tol: 
        f_Y = ((Y + Y ** 2 + Y ** 3 - Y ** 4) / (1 - Y) ** 3) - A * P_pr - \
              B * Y ** 2 + C * Y ** D
        df_dY = ((1 + 4 * Y + 4 * Y ** 2 - 4 * Y ** 3 + Y ** 4) / \
                 (1 - Y) ** 4) - 2 * B * Y + C * D * Y ** (D - 1)
        Y = Y - f_Y / df_dY
    z = A * P_pr / Y

    return z

def api(gamma_o):
    """
    Conversion from oil gravity to API gravity.

    Parameters
    ----------
    gamma_o : oil gravity (unitless)

    Returns 
    -------
    api : api gravity (°API)
    """

    return 141.5 / gamma_o - 131.5

def oil_gravity(API):
    """
    Conversion from API gravity to to oil gravity.

    Parameters
    ----------
    API : api gravity (°API)

    Returns
    -------
    gamma_o : oil gravity (unitless)
    """

    return 141.5 / (131.5 + API)

def normalized_gas_gravity(P_sep, T_sep, API, gamma_g):
    """
    Normalized gas gravity at the reference separator pressure and 
    temperature. Vazques-Beggs normalize the gas gravity to a separator 
    pressure of 100 psig to account for the dependence of gas gravity on
    seperator conditions.

    Parameters
    ----------
    P_sep : separator pressure (psia)
    T_sep : separator temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)

    Returns
    -------
    gamma_gn : normalized gas gravity (unitless)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    gamma_gn = gamma_g * (1 + 5.912e-5 * API * T_sep * \
                          log10(P_sep / 114.7))
    
    return gamma_gn

def spivey_oil_compressibility(P, P_b, T, API, gamma_gsp, R_sb):
    """
    Spivey isothermal oil compressibility correlation for pressures at or 
    above bubblepoint pressure. The average relative error and average 
    absolute relative error are 0.1% and 7.1%, respectively. The correlation 
    is applicable within the following limits:

    11.6 <= API <= 57.7°API
    0.561 <= gamma_gsp <= 1.798
    120.7 <= Pb <= 6658.7 psia
    414.7 <= P <= 8014.7 psia
    12 <= R_sb <= 1808 SCF/STB
    70.7 <= T <= 320°F
    
    Parameters
    ----------
    P : pressure (psia)
    P_b : bubble-point pressure (psia)
    T : temperature (°F)
    API : oil gravity (°API)
    gamma_gsp : seperator gas specific gravity (unitless)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)

    Returns
    -------
    c_o : oil compressibility (psi**-1)

    References
    ----------
    [5]: Ezekew, Nnaemeka. *Petroleum Reservoir Engineering Practice.* 
    Pearson Education, 2010, pp. 96-97
    """
    
    Z_1 = 3.011 - 2.6254 * log(API) + 0.497 * log(API) ** 2
    Z_2 = -0.0835 - 0.259 * log(gamma_gsp) + 0.382 * log(gamma_gsp) ** 2
    Z_3 = 3.51 - 0.0289 * log(P_b) - 0.0584 * log(P_b) ** 2
    Z_4 = 0.327 - 0.608 * log(P / P_b) + 0.0911 * log(P / P_b) ** 2
    Z_5 = -1.918 - 0.642 * log(R_sb) + 0.154 * log(R_sb) ** 2
    Z_6 = 2.52 - 2.73 * log(T) + 0.429 * log(T) ** 2
    
    Z = Z_1 + Z_2 + Z_3 + Z_4 + Z_5 + Z_6
    
    ln_c_ob = 2.434 + 0.475 * Z + 0.048 * Z ** 2
    #equation returns in units of microsips (1e-6 psi**-1) thus needs to be 
    #multiplied by 1e-6 to obtain in units of psi**-1
    return exp(ln_c_ob) * 1e-6

def mccain_oil_compressibility(P, P_b, T, API, R_sb):
    """
    McCain oil compressibility correlation at pressures below the bubblepoint 
    pressure. This correlation is valid up to 5300 psia and 330°F.

    Parameters
    ----------
    P : pressure (psia)
    P_b : bubble-point pressure (psia)
    T : temperature (°F)
    API : oil gravity (°API)
    R_sb: bubble-point solution gas-oil ratio (SCF/STB)

    Returns
    -------
    c_o : oil compressibility (psi**-1)

    References
    ----------
    [6]: McCain, William Jr. *The Properties of Petroleum Fluids.* 2nd ed., 
    PennWell Books, 1989, pp. 523
    """

    T += 459.67
    ln_c_o = -7.573 - 1.45 * log(P) - 0.383 * log(P_b) + 1.402 * \
        log(T) + 0.256 * log(API) + 0.449 * log(R_sb)

    return exp(ln_c_o)

def vasquez_beggs_oil_compressibility(P, T, P_sep, T_sep, API, gamma_g, R_sb):
    """
    Vasquez-Beggs oil compressibility correlation for pressures at or above 
    bubblepoint pressure. The correlation is applicable within the following 
    limits:

    15 <= Pb <= 6055 psia
    75 <= T <= 294°F
    0 <= Rsb <= 2199 SCF/STB
    15.3 <= API <= 59.3 °API
    0.51 <= gamma_g <= 1.35

    Parameters
    ----------
    P : pressure (psia)
    T : temperature (°F)
    P_sep : separator pressure (psia)
    T_sep : separator temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)
 
    Returns
    -------
    c_o : oil compressibility (psi**-1)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    gamma_gn = normalized_gas_gravity(P_sep, T_sep, API, gamma_g)
    c_o = (-1433 + 5 * R_sb + 17.2 * T - 1180 * gamma_gn + 12.61 * API) / (1e5 * P)

    return c_o

def petrosky_farshad_oil_compressibility(P, T, API, gamma_g, R_sb):
    """
    Petrosky-Farshad oil compressibility correlation for pressures at or above
    bubblepoint pressure. The correlation is applicable within the following
    limits:

    1574 <= P_b <= 4640 psia
    114 <= T <= 288°F
    217 <= R_sb <= 1406 SCF/STB
    16.3 <= API <= 45 °API
    0.58 <= gamma_g <= 0.86

    Parameters
    ----------
    P : pressure (psia)
    T : temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)

    Returns
    -------
    c_o : oil compressibility (psi**-1)

    References
    ----------
    [1]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    c_o = 1.705e-7 * R_sb ** 0.69357 * gamma_g ** 0.1885 * API ** 0.3272 * \
          T ** 0.6729 * P ** (-0.5906)
    
    return c_o

def standing_p_b(T_res, API, gamma_g, R_sb):
    """
    Standing bubblepoiont pressure correlation. The correlation is applicable
    within the following limits:

    130 <= P_b <= 7000 psia
    100 <= T <= 258°F
    20 <= R_sb <= 1425 SCF/STB
    16.5 <= API <= 63.8 °API
    0.59 <= gamma_g <= 0.95

    Parameters
    ----------
    T_res : reservoir temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas gravity (unitless)
    R_sb : solution gas-oil ratio at bubble-point pressure (SCF/STB)

    Returns
    -------
    P_b : bubble-point pressure (psia)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """
    
    C_pb = (R_sb / gamma_g) ** 0.83 * 10 ** (0.00091 * T_res - 0.0125 * API)
    P_b = 18.2 * (C_pb - 1.4)

    return P_b

def valko_mccain_p_b(T_res, API, gamma_g, R_sb):
    """
    Valko and McCain bubblepoint pressure correlation. The average relative
    error and the absolute relative error are 0.0% and 10.9%, respectively. 
    The correlation is applicable within the following limits:

    31.7 <= P_b <= 7127 psia
    74 <= T <= 341.6°F
    6.0 <= R_sb <= 3298.6 SCF/STB
    6 <= API <= 63.7 °API
    0.51 <= gamma_g <= 3.44    
    
    Parameters
    ----------
    T_res : reservoir temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas gravity (unitless)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)

    Returns
    -------
    P_b : bubble-point pressure (psia)

    References
    ----------
    [5]: Ezekew, Nnaemeka. *Petroleum Reservoir Engineering Practice.* 
    Pearson Education, 2010, pp. 93-95
    """
    
    Z_1 = -5.48 - 0.0378 * log(R_sb) + 0.281 * log(R_sb) ** 2 -\
          0.0206 * log(R_sb) ** 3
    Z_2 = 1.27 - 0.0449 * API + 4.36e-4 * API ** 2 - 4.76e-6 * API ** 3
    Z_3 = 4.51 - 10.84 * gamma_g + 8.39 * gamma_g ** 2 - 2.34 * gamma_g ** 3
    Z_4 = -0.7835 + 6.23e-3 * T_res -1.22e-5 * T_res ** 2 + 1.03e-8 * \
        T_res ** 3

    Z = Z_1 + Z_2 + Z_3 + Z_4

    ln_P_b = 7.475 + 0.7132 * Z + 0.0075 * Z ** 2

    return exp(ln_P_b)

def vasquez_beggs_p_b(T_res, P_sep, T_sep, API, gamma_g, R_sb):
    """
    Vasquez Beggs bubblepoint pressure correlation. The correlation is 
    applicable within the following limits:

    15 <= P_b <= 6055 psia
    75 <= T <= 294°F
    0 <= R_sb <= 2199 SCF/STB
    15.3 <= API <=59.3 °API
    0.51 <= gamma_g <= 1.53   

    Parameters
    ----------
    T_res : reservoir temperature (°F)
    P_sep : seperator pressure (psia)
    T_sep : seprator temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas gravity (unitless)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)

    Returns
    -------
    P_b : bubble-point pressure (psia)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    if API > 30:
        A, B, C = 56.06, -10.393, 0.8425
    else:
        A, B, C  = 27.64, -11.172, 0.9143

    gamma_gn = normalized_gas_gravity(P_sep, T_sep, API, gamma_g)
    a = B * API / (T_res + 459.67)
    P_b = (A * (R_sb / gamma_gn) * 10 ** a) ** C

    return P_b

def al_marhoun_p_b(T_res, API, gamma_g, R_sb):
    """
    Al-Marhoun bubblepoint pressure correlation. The correlation is applicable
    within the following limits:

    20 <= P_b <= 3573 psia
    75 <= T <= 240°F
    24 <= R_sb <= 1901 SCF/STB
    14.3 <= API <= 44.6 °API
    0.752 <= gamma_g <= 1.367   

    Parameters
    ----------
    T_res : reservoir temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas gravity (unitless)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)

    Returns
    -------
    P_b : bubble-point pressure (psia)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    T_res += 459.67
    a = 0.722569
    b = -1.879109
    c = 3.04659
    d = 1.302347

    gamma_o = oil_gravity(API)
    X = R_sb ** a * gamma_g ** b * gamma_o ** c * T_res ** d
    P_b = -64.13891 + 7.02362e-3 * X - 2.278475e-9 * X ** 2

    return P_b

def glaso_p_b(T_res, API, gamma_g, R_sb):
    """
    Glaso bubblepoint pressure correlation. The correlation is applicable
    within the following limits:

    165 <= P_b <= 7142 psia
    80 <= T <= 280°F
    90 <= R_sb <= 2637 SCF/STB
    22.3 <= API <= 48.1 °API
    0.65 <= gamma_g <= 1.28   

    Parameters
    ----------
    T_res : reservoir temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas gravity (unitless)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)

    Returns
    -------
    P_b : bubble-point pressure (psia)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    X = (R_sb / gamma_g) ** 0.816 * (T_res ** 0.172 / API ** 0.989)
    log_P_b = 1.7669 + 1.7447 * log10(X) - 0.30218 * log10(X) ** 2

    return 10 ** log_P_b

def petrosky_farshad_p_b(T_res, API, gamma_g, R_sb):
    """
    Glasso bubblepoint pressure correlation. The correlation is applicable
    within the following limits:

    1574 <= P_b <= 4640 psia
    114 <= T <= 288°F
    217 <= R_sb <= 1406 SCF/STB
    16.3 <= API <= 45 °API
    0.58 <= gamma_g <= 0.86

    Parameters
    ----------
    T_res : reservoir temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas gravity (unitless)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns
    -------
    P_b : bubble-point pressure (psia)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    a = 7.916e-4 * API ** 1.541 - 4.561e-5 * T_res ** 1.3911
    P_b = 112.727 * (R_sb ** 0.577421 / (gamma_g ** 0.8439 * 10 ** a) - 12.34)

    return P_b

def standing_r_s(P, T, API, gamma_g):
    """
    Standing solution gas-oil ratio correlation. The correlation is applicable
    within the following limits:

    130 <= P_b <= 7000 psia
    100 <= T <= 258°F
    20 <= R_sb <= 1425 SCF/STB
    16.5 <= API <= 63.8 °API
    0.59 <= gamma_g <= 0.95
    
    Parameters
    ----------
    P : pressure (psia)
    T : temeperature(°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)

    Returns
    -------
    R_s : solution gas-oil ratio (SCF/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """
    
    R_s = gamma_g * (((P / 18.2) + 1.4) * (10 ** (0.0125 * API)) / \
                     (10 ** (0.00091 * T))) ** 1.2048
    
    return R_s

def vasquez_beggs_r_s(P, T, P_sep, T_sep, API, gamma_g):
    """
    Vasquez-Beggs solution gas-oil ratio correlation. The correlation is 
    applicable within the following limits:

    15 <= P_b <= 6055 psia
    75 <= T <= 294°F
    0 <= R_sb <= 2199 SCF/STB
    15.3 <= API <=59.3 °API
    0.51 <= gamma_g <= 1.53 
    
    Parameters
    ----------
    P : pressure (psia)
    T : temeperature(°F)
    P_sep : seperator pressure (psia)
    T_sep : seprator temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas gravity (unitless)

    Returns
    -------
    R_s : solution gas-oil ratio (SCF/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    if API > 30:
        A, B, C = 0.0178, 1.1870, 23.931
    else:
        A, B, C = 0.0362, 1.0937, 25.7240
    
    T += 459.67
    gamma_gn = normalized_gas_gravity(P_sep, T_sep, API, gamma_g)
    R_s = A * gamma_gn * P ** B * exp(C * (API / T))

    return R_s

def al_marhoun_r_s(P, T, API, gamma_g):
    """
    Al-Marhoun solution gas-oil ratio correlation. The correlation is 
    applicable within the following limits:

    20 <= P_b <= 3573 psia
    75 <= T <= 240°F
    24 <= R_sb <= 1901 SCF/STB
    14.3 <= API <= 44.6 °API
    0.752 <= gamma_g <= 1.367   

    Parameters
    ----------
    P : pressure (psia)
    T : temeperature(°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)

    Returns
    -------
    R_s : solution gas-oil ratio (SCF/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    T += 459.67
    a = -2.278475e-9
    b = 7.02362e-3
    c = -64.13891 - P
    d = -1.879109
    e = 3.04659
    f = 1.302347
    g = 0.722569

    gamma_o = oil_gravity(API)
    x = (-b + sqrt(b ** 2 - 4 * a * c)) / (2 * a)

    R_s = (x / (gamma_g ** d * gamma_o ** e * T ** f)) ** (1 / g)

    return R_s

def glaso_r_s(P, T, API, gamma_g):
    """
    Glaso solution gas-oil ratio correlation. The correlation is applicable
    within the following limits:

    165 <= P_b <= 7142 psia
    80 <= T <= 280°F
    90 <= R_sb <= 2637 SCF/STB
    22.3 <= API <= 48.1 °API
    0.65 <= gamma_g <= 1.28   
    
    Parameters
    ----------
    P : pressure (psia)
    T : temeperature(°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)

    Returns
    -------
    R_s : solution gas-oil ratio (SCF/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    a = -0.30218
    b = 1.7447
    c = 1.7669 - log10(P)

    X = (-b + sqrt(b ** 2 - 4 * a * c)) / (2 * a)

    R_s = gamma_g * (10 ** X * API ** 0.989 / T ** 0.172) ** (1 / 0.816)

    return R_s

def petrosky_farshad_r_s(P, T, API, gamma_g):
    """
    Petrosky-Farshad solution gas-oil ratio correlation. The correlation is 
    applicable within the following limits:

    1574 <= P_b <= 4640 psia
    114 <= T <= 288°F
    217 <= R_sb <= 1406 SCF/STB
    16.3 <= API <= 45 °API
    0.58 <= gamma_g <= 0.86
    
    Parameters
    ----------
    P : pressure (psia)
    T : temeperature(°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)

    Returns
    -------
    R_s : solution gas-oil ratio (SCF/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    a = 7.916e-4 * API ** 1.541 - 4.561e-5 * T ** 1.3911
    R_s = ((P / 112.727 + 12.34) * gamma_g ** 0.8439 * 10 ** a) ** 1.73184

    return R_s

def standing_beta_o(T, API, gamma_g, R_s):
    """
    Standing oil formation volume factor correlation. The correlation is 
    applicable within the following limits:

    130 <= P_b <= 7000 psia
    100 <= T <= 258°F
    20 <= R_sb <= 1425 SCF/STB
    16.5 <= API <= 63.8 °API
    0.59 <= gamma_g <= 0.95
    
    Parameters
    ----------
    T : temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns
    -------
    beta_o : oil formation volume factor (RB/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """
    
    gamma_o = oil_gravity(API)
    beta_o = 0.972 + 1.47e-4 * (R_s * (gamma_g / gamma_o) ** 0.5 + 1.25 * T) \
        ** 1.175
    
    return beta_o

def vasquez_beggs_beta_o(T, P_sep, T_sep, API, gamma_g, R_s):
    """
    Vasquez-Beggs oil formation volume factor correlation. The correlation is
    applicable within the following limits:

    15 <= P_b <= 6055 psia
    75 <= T <= 294°F
    0 <= R_sb <= 2199 SCF/STB
    15.3 <= API <=59.3 °API
    0.51 <= gamma_g <= 1.53

    Parameters
    ----------
    T : temperature (°F)
    P_sep : seperator pressure (psia)
    T_sep : seprator temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns
    -------
    beta_o : oil formation volume factor (RB/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    if API > 30:
        A, B, C = 4.67e-4, 1.1e-5, 1.337e-9
    else:
        A, B, C = 4.677e-4, 1.751e-5, -1.8106e-8

    gamma_gn = gamma_g * (1 + 5.912e-5 * API * T_sep * \
                          log10(P_sep / 114.7))
    beta_o = 1 + A * R_s + (T - 60) * (API / gamma_gn) * (B + C * R_s)

    return beta_o

def al_marhoun_beta_o(T, API, gamma_g, R_s):
    """
    Al-Marhoun oil formation volume factor correlation. The correlation is 
    applicable within the following limits:

    20 <= P_b <= 3573 psia
    75 <= T <= 240°F
    24 <= R_sb <= 1901 SCF/STB
    14.3 <= API <= 44.6 °API
    0.752 <= gamma_g <= 1.367   

    Parameters
    ----------
    T : temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns
    -------
    beta_o : oil formation volume factor (RB/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    T += 459.67
    gamma_o = oil_gravity(API)
    X = R_s ** 0.74239 * gamma_g ** 0.323294 * gamma_o ** -1.20204

    beta_o = 0.497069 + 8.62963e-4 * T + 1.82594e-3 * X + 3.18099e-6 * X ** 2

    return beta_o

def glaso_beta_o(T, API, gamma_g, R_s):
    """
    Glaso oil formation volume factor correlation. The correlation is 
    applicable within the following limits:

    165 <= P_b <= 7142 psia
    80 <= T <= 280°F
    90 <= R_sb <= 2637 SCF/STB
    22.3 <= API <= 48.1 °API
    0.65 <= gamma_g <= 1.28   

    Parameters
    ----------
    T : temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns
    -------
    beta_o : oil formation volume factor (RB/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    gamma_o = oil_gravity(API)
    X = R_s * (gamma_g / gamma_o) ** 0.526 +0.968 * T

    beta_o = 1 + 10 ** (-6.58511 + 2.91329 * log10(X) - 0.27683 * \
                        log10(X) ** 2)

    return beta_o

def petrosky_farshad_beta_o(T, API, gamma_g, R_s):
    """
    Petrosky-Farshad oil formation volume factor correlation. The correlation 
    is applicable within the following limits:

    1574 <= P_b <= 4640 psia
    114 <= T <= 288°F
    217 <= R_sb <= 1406 SCF/STB
    16.3 <= API <= 45 °API
    0.58 <= gamma_g <= 0.86
    
    
    Parameters
    ---------- 
    T : temperature (°F)
    API : oil gravity (°API)
    gamma_g : gas specific gravity (unitless)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns
    -------
    beta_o : oil formation volume factor (RB/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    gamma_o = oil_gravity(API)

    beta_o = 1.0113 + 7.2046e-5 * (R_s ** 0.3738 * (gamma_g ** 0.2914 / \
                gamma_o ** 0.6265) + 0.24626 * T ** 0.5371) ** 3.0936

    return beta_o

def undersaturated_beta_o(P, P_b, beta_ob, c_o):
    """
    Undersaturated calculation of oil formation volume factor.

    Parameters
    ----------
    P : pressure (psia)
    P_b : bubblepoint pressure (psia)
    beta_ob : formation oil volume factor at bubblepoint pressure (RB/STB)
    c_o : undersaturated oil compressibility (psi**-1)

    Returns
    -------
    beta_o : oil formation volume factor (RB/STB)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 276.
    """

    beta_o = beta_ob * exp(-c_o * (P - P_b))

    return beta_o

def standing_mu_od(T, API):
    """
    Standing dead oil viscosity correlation. The correlation limits are not
    listed.
    
    Parameters
    ----------
    T : temperature (°F)
    API : oil gravity (°API)

    Returns
    -------
    mu_od : dead oil viscosity (cp)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """
    
    A = 10 ** (0.42 + (8.33 / API))
    mu_od = (0.32 + (1.8e7 / API ** 4.53)) * (360 / (T + 200)) ** A
    
    return mu_od

def standing_mu_ob(T, API, R_s):
    """
    Standing saturated oil viscosity correlation. The correlation limits are not
    given.
    
    Parameters
    ----------
    T : temperature (°F)
    API : oil gravity (°API)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns
    -------
    mu_ob : saturated oil viscosity (cp)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """
    
    e = 3.74e-3 * R_s
    d = 1.1e-3 * R_s
    c = 8.68e-5 * R_s
    b = (0.68 / 10 ** c) + (0.25 / 10 ** d) + (0.062 / 10 ** e)
    a = R_s * (2.2e-7 * R_s - 7.4e-4)
    mu_od = standing_mu_od(API, T)
    mu_ob = 10 ** a * mu_od ** b

    return mu_ob

def ng_egbogah_mu_od(T, API):
    """
    Ng and Egbogah dead oil viscosity correlation. The correlation 
    is applicable within the following limits:

    59 <= T <= 176°F
    5 <= API <= 58 °API
    
    Parameters
    ----------
    
    T : temperature (°F)
    API : oil gravity (°API)

    Returns
    -------
    mu_od : dead oil viscosity (cp)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """
    
    log_10_mu_od = 1.8653 - 2.5086e-2 * API - 0.5644 * log10(T)
    mu_od = 10 ** 10 ** log_10_mu_od - 1
    
    return mu_od

def beggs_robinson_mu_od(T, API):
    """
    Ng and Egbogah dead oil viscosity correlation. The correlation 
    is applicable within the following limits:

    70 <= T <= 295°F
    16 <= API <= 58 °API
    
    Parameters
    ----------
    T : temperature (°F)
    API : oil gravity (°API)

    Returns
    -------
    mu_od : dead oil viscosity (cp)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    X = 10 ** (3.0324 - 0.02023 * API) * T ** -1.163
    mu_od = 10 ** X - 1

    return mu_od

def beggs_robinson_mu_ob(T, API, R_s):
    """
    Beggs and Robinson saturated oil viscosity correlation. The correlation 
    is applicable within the following limits:

    70 <= T <= 295°F
    16 <= API <= 58 °API
    20 <= R_s <= 2070 SCF/STB 
    
    Parameters
    ----------
    T : temperature (°F)
    API : oil gravity (°API)
    R_s : solution gas-oil ratio (SCF/STB)

    Returns
    -------
    mu_ob: saturated oil viscosity (cp)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """
    
    a = 10.715 * (R_s + 100) ** -0.515
    b = 5.44 * (R_s + 150) ** -0.338
    mu_od = beggs_robinson_mu_od(T, API)
    mu_ob = a * mu_od ** b

    return mu_ob

def bergman_sutton_mu_o(P, P_b, T, API, R_sb):
    """
    Bergman and Sutton undersaturated oil viscosity correlation. The 
    correlation is applicable within the following limits:

    0 <= P <= 5265 psia
    70 <= T <= 295°F
    16 <= API <= 58 °API
    20 <= R_sb <= SCF/STB
    
    Parameters
    ----------
    P : pressure (psia)
    P_b : bubblepoint pressure (psia)
    T : temperature (°F)
    API : oil gravity (°API)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)  

    Returns
    -------
    mu_o : undersaturated oil viscosity (cp)

    References
    ----------
    [5]: Ezekew, Nnaemeka. *Petroleum Reservoir Engineering Practice.* 
    Pearson Education, 2010, pp. 99
    """
    
    mu_ob = beggs_robinson_mu_ob(T, API, R_sb)
    alpha = 6.5698e-7 * log(mu_ob) ** 2 - 1.48211e-5 * log(mu_ob) + \
            2.27877e-4
    beta = 2.24623e-2 * log(mu_ob) + 0.873204
    mu_o = mu_ob * exp(alpha * (P - P_b) ** beta)

    return mu_o

def vasquez_beggs_mu_o(P, P_b, T, API, R_sb):
    """
    Vasquez-Beggs undersaturated-oil viscosity correlation. The correlation
    limits are not given.

    Parameters
    ----------
    P : pressure (psia)
    P_b : bubblepoint pressure (psia)
    T : temperature (°F)
    API : oil gravity (°API)
    R_sb : solution gas-oil ratio at bubblepoint pressure (SCF/STB)  

    Returns
    -------
    mu_o : undersaturated oil viscosity (cp)

    References
    ----------
    [6]: McCain, William Jr. *The Properties of Petroleum Fluids.* 2nd ed., 
    PennWell Books, 1989, pp. 524
    """
    
    B = 2.6 * P ** 1.187 * exp(-11.513 - 8.98e-5 * P)
    mu_ob = beggs_robinson_mu_ob(T, API, R_sb)
    mu_o = mu_ob * (P / P_b) ** B

    return mu_o

def glaso_mu_od(T, API):
    """
    Glaso dead-oil viscosity correlation. The correlation is applicable within
    the following limits:

    50 <= T <= 300°F
    20.1 <= API <= 48.1 °API 

    Parameters
    ----------
    T : temperature (°F)
    API : oil gravity (°API)

    Returns
    -------
    mu_od: dead oil viscosity (cp)

    References
    ----------
    [4]: Sutton, Robert P. *Oil System Correlations.* 
    General Engineering vol 1., edited by John R. Fanchi, 
    Society of Petroleum Engineers, 2006, pp. 310-331.
    """

    A = 10.313 * log10(T) - 36.447
    mu_od = 3.141e10 * T ** -3.444 * log10(API) ** A

    return mu_od

def WOR_calculation(WC, q_l):
    """
    Compute the water-oil ratio from the water-cut and liquid production
    rate.
    
    Recall that the water-cut is defined as:
    
                                WC = q_w / q_l
                                
    Thus q_w is
    
                                q_w = WC * q_l
    
    WOR can then be caluclated
    
                                WOR = WC * q_l / (q_l - WC * q_l)

    Parameters
    ----------
    WC : water-cut (unitless)
    q_l : liquid flow rate (STB/D)

    Returns
    -------
    WOR : Water-oil ratio
    """
    
    WOR = WC * q_l / (q_l - (WC * q_l))
    
    return WOR

def dens_g_per_cc(rho_w):
    """
    Conversion from lbm/ft**3 to g/cc
    
    Parameters
    ----------
    rho_w : density (lbm/ft**3)

    Returns
    -------
    rho_w : water density (g/cc)
    """
    
    rho_w = rho_w * 453.592 / (2.54 * 12) ** 3
    
    return rho_w

def ppm_to_weight_percent_solids(ppm):
    """
    Conversion from TDS to weight percent solids. 
    
    Parameters
    ----------
    ppm : total dissolved solids (ppm)

    Returns
    -------
    C_w : weight percent solids (%)
    """
    
    C_w = ppm * 1e-4
    
    return C_w

def ppm_to_mg_per_l(ppm, rho_w):
    """
    Conversion form ppm to mg/l.
    
    Parameters
    ----------
    ppm : salinity (ppm)
    rho_w : brine density (lbf/ft**3)

    Returns
    -------
    mgl : total dissolved solids (mg/l)
    """
    
    rho_w = dens_g_per_cc(rho_w)
    mgl = ppm * rho_w
    
    return mgl

def brine_density(S):
    """
    Calculate brine density at standard condtions from salinity.
    
    Parameters
    ----------
    S : salinity in weight percent solids (unitless)

    Returns 
    -------
    rho_w : brine density (lbm/ft**3)
    """
    
    rho_w = 62.368 + 0.438603 * S +1.60074 * 1e-3 * S ** 2 
    
    return rho_w

def mccain_mu_w(T, TDS):
    """
    McCain atmospheric water viscosity correlation. The user must be aware
    that the TDS parameter is in (ppm). The solids present in oilfield 
    waters are reported in several different ways (ppm, mg/l, wt% solids)
    and therefore must be converted into ppm for this correlation to give
    correct atmospheric water viscosity estimates. Reported error is 
    at 5%.
    
    Parameters
    ----------
    T : temeperature (°F)
    TDS : total dissolved solids (ppm)

    Returns 
    -------
    mu_w_atm : water viscosity at atmospheric pressure (cp)

    References
    ----------
    [6]: McCain, William Jr. *The Properties of Petroleum Fluids.* 2nd ed., 
    PennWell Books, 1989, pp. 527
    """
    
    S = TDS / 10e6 * 10e2
    A = 109.574 - 8.40564 * S + 0.313314 * S ** 2 + 8.72213e-3 * S ** 3
    B = 1.12166 - 2.63951e-2 * S + 6.79461e-4 * S ** 2 + \
        5.47119e-5 * S ** 3 - 1.55586e-6 * S ** 4
    mu_w_atm = A * T ** -B

    return mu_w_atm

def matthews_russel_mu_w(P, T, TDS):
    """
    Matthews-Russel pressure-adjusted water viscosity correlation.
    The user must be aware that the TDS parameter is in (ppm). The 
    solids present in oilfield waters are reported in several different 
    ways (ppm, mg/l, wt% solids) and therefore must be converted into 
    ppm for this correlation to give correct pressure-adjusted water 
    viscosity estimates.
    
    Parameters
    ----------
    P : pressure (psia)
    T : temperature (°F)
    TDS : total dissolved solids (ppm)

    Returns 
    -------
    mu_w : water viscosity reservoir pressure (cp)

    References
    ----------
    [6]: McCain, William Jr. *The Properties of Petroleum Fluids.* 2nd ed., 
    PennWell Books, 1989, pp. 527
    """
    
    mu_w_atm = mccain_mu_w(T, TDS)
    mu_w = mu_w_atm * (0.9994 + 4.0295e-5 * P + 3.1062e-9 * P ** 2)

    return mu_w

def beta_w(P, T):
    """
    Water formation volume factor correlation. 

    Parameters
    ----------
    P : pressure (psia)
    T : temperature (°F)

    Returns
    -------
    beta_w : water formation volume factor (RB/STB)

    References
    ----------
    [6]: McCain, William Jr. *The Properties of Petroleum Fluids.* 2nd ed., 
    PennWell Books, 1989, pp. 525
    """

    delta_wt = -1.0001 * 1e-2 + 1.33391 * 1e-4 * T + 5.50654 * 1e-7 * T ** 2
    delta_wp = -1.95301 * 1e-9 * P * T - 1.72834 * 1e-13 * P ** 2 * T - \
        3.58922 * 1e-7 * P - 2.25341 * 1e-10 * P ** 2
    
    beta_w = (1 + delta_wp) * (1 + delta_wt)

    return beta_w

def beta_g(P, T, z):
    """
    
    """

    P_sc = 14.696
    T_sc = 80 + 459.67
    z_sc = 1

    beta_g = (z * T * P_sc) / (z_sc * T_sc * P)

    return beta_g