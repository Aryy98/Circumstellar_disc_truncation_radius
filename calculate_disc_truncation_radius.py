# ========= Script to compute the outer disc truncation radius for S-type binaries =================
# Main developer: A.Nigioni. Co-developer: J.Venturini (September 2025).
# It uses the fit of Manara et al. 2019, where the truncation radius is calculated based in the binary parameters and Reynolds number.
# Instead of Reynolds number, we set the dependence on the disc alhpa viscosity parameter and use planet formation simulations to convert to Reynolds number.
# (see details in Venturini et al. 2025).


import numpy as np
import pandas as pd


# DEFINE BINARY PARAMETERS
# Gamma-Cephei as an example:
M1 = 1.27 #primary star mass in Msun
M2 = 0.328 #secondary star mass in Msun
abin = 19.56 #binary separation in au
ebin = 0.41 #binary eccentricity
Stype_truncation_secondary = False #Should be set to 'True' if one wants the truncation radius of the circumsecondary disc.


# DISC'S VISCOSITY
alpha = 0.001


# --------------FUNCTION DEFINITION-------------------
# Definition of the function to calculate the truncation radius as in Manara+2019 (their Appendix C)
#
def calculate_truncation(Mstar,Mstar2,e_star2,alpha_disk, a_star2,Stype_truncation_secondary):
    H_over_r_from_populations = 0.0983 
    log10Re = [np.log10(1e4), np.log10(1e5), np.log10(1e6)] 
    mu = [0.1, 0.2, 0.3, 0.4, 0.5]   
    b_Re_1e4_circumprimary = [-0.66,-0.72,-0.76,-0.77,-0.78]  
    b_Re_1e5_circumprimary = [-0.75,-0.78,-0.80,-0.81,-0.81] 
    b_Re_1e6_circumprimary = [-0.78,-0.80,-0.81,-0.82,-0.82]
    c_Re_1e4_circumprimary = [0.84,0.88,0.92,0.95,0.94]
    c_Re_1e5_circumprimary = [0.68,0.72,0.75,0.78,0.78]
    c_Re_1e6_circumprimary = [0.56,0.60,0.63,0.66,0.66]
    b_Re_1e4_circumsecondary = [-0.81,-0.81,-0.79,-0.80,-0.79]
    b_Re_1e5_circumsecondary = [-0.81,-0.82,-0.82,-0.82,-0.81]
    b_Re_1e6_circumsecondary = [-0.83,-0.83,-0.83,-0.83,-0.82]
    c_Re_1e4_circumsecondary = [0.98,0.99,0.97,0.98,0.95]
    c_Re_1e5_circumsecondary = [0.80,0.82,0.81,0.80,0.78]
    c_Re_1e6_circumsecondary = [0.69,0.70,0.69,0.68,0.66]

    b_ = []
    c_ = []

    mu_Manara = Mstar2 / (Mstar + Mstar2)
    
    # Eggleton term
    if (Stype_truncation_secondary == False):
        q13 = (Mstar / Mstar2) ** (1.0 / 3.0)
        q23 = q13 ** 2
        f_Eggleton = (0.49 * q23) / (0.6 * q23 + np.log(1.0 + q13)) 
        
    else:
        q13 = (Mstar2 / Mstar) ** (1.0 / 3.0)
        q23 = q13 ** 2
        f_Eggleton = (0.49 * q23) / (0.6 * q23 + np.log(1.0 + q13)) 
        
    # calculate Reynolds number
    if (e_star2 > 0.0):
        Reynolds = (alpha_disk**-1)*(H_over_r_from_populations**-2)
        if ((Reynolds < 1e4)):
            Reynolds = 1e4
        if ((Reynolds > 1e6)):
            Reynolds = 1e6  
        logReynolds = np.log10(Reynolds)

    # calculate truncation radius
    if (Stype_truncation_secondary == False):
        if (e_star2>0.0):
            for i in range(len(mu)):
                DEGREE = 2
                X = log10Re
                Y = [b_Re_1e4_circumprimary[i], b_Re_1e5_circumprimary[i], b_Re_1e6_circumprimary[i]]
                coeff_b_Re = np.polyfit(X, Y, deg = DEGREE)
                Y = [c_Re_1e4_circumprimary[i], c_Re_1e5_circumprimary[i], c_Re_1e6_circumprimary[i]]
                coeff_c_Re = np.polyfit(X, Y, deg = DEGREE)

                b_.append(np.polyval(coeff_b_Re, logReynolds))
                c_.append(np.polyval(coeff_c_Re, logReynolds))

            DEGREE = 4
            coeff_b_mu = np.polyfit(mu, b_, deg = DEGREE)
            coeff_c_mu = np.polyfit(mu, c_, deg = DEGREE)

            b = np.polyval(coeff_b_mu, mu_Manara)
            c = np.polyval(coeff_c_mu, mu_Manara)  

            f_Manara = f_Eggleton*(b*(e_star2**c) + 0.88*((mu_Manara)**0.01))
        else:
            f_Manara = f_Eggleton*0.88*((mu_Manara)**0.01)
    else:
        if (e_star2>0.0):
            for i in range(len(mu)):
                DEGREE = 2
                X = log10Re
                Y = [b_Re_1e4_circumsecondary[i], b_Re_1e5_circumsecondary[i], b_Re_1e6_circumsecondary[i]]
                coeff_b_Re = np.polyfit(X, Y, deg = DEGREE)
                Y = [c_Re_1e4_circumsecondary[i], c_Re_1e5_circumsecondary[i], c_Re_1e6_circumsecondary[i]]
                coeff_c_Re = np.polyfit(X, Y, deg = DEGREE)

                b_.append(np.polyval(coeff_b_Re, logReynolds))
                c_.append(np.polyval(coeff_c_Re, logReynolds))

            DEGREE = 4
            coeff_b_mu = np.polyfit(mu, b_, deg = DEGREE)
            coeff_c_mu = np.polyfit(mu, c_, deg = DEGREE)

            b = np.polyval(coeff_b_mu, mu_Manara)
            c = np.polyval(coeff_c_mu, mu_Manara)  

            f_Manara = f_Eggleton*(b*(e_star2**c) + 0.88*((mu_Manara)**0.01))
        else:
            f_Manara = f_Eggleton*0.88*((mu_Manara)**0.01)       

    return a_star2 * f_Manara
    
# --------------END FUNCTION SECTION-------------------


# Print the truncation radius for the defined binary parameters in Rt:
Rt = calculate_truncation(M1,M2,ebin,alpha,abin,Stype_truncation_secondary)

if Stype_truncation_secondary == False:
    print('For the defined binary parameters and disc viscosity, the outer trancation radius around the primary star [in au] is:')
    print(Rt)
else:
    print('For the defined binary parameters and disc viscosity, the outer trancation radius around the secondary star [in au] is:')
    print(Rt)
