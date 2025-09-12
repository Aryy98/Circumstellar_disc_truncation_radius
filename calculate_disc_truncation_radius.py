#%%
import numpy as np
import pandas as pd
from pathlib import Path as path
import matplotlib.pyplot as plt

# -------------- READ PARAMETERS FROM FILE-------------------

data = pd.read_csv('input.csv', sep=',')

alpha = data['alpha'] #viscosity
M1 = data['Mstar1'] #primary star mass in Msun
M2 = data['Mstar2'] #secondary star mass in Msun  
abin =  data['a_bin'] #binary separation in au   
ebin = data['e_bin'] #binary eccentricity

# --------------END PARAMETERS SECTION-------------------

# --------------PHTSICAL CONSTANTS-------------------
# Costanti fisiche in CGS
Lsun = 3.828e33        # erg/s
Msun = 1.989e33        # g
au = 1.496e13          # cm
sigma = 5.670374419e-5 # erg/cm^2/s/K^4
pi = np.pi
# --------------END PHTSICAL CONSTANTS SECTION-------------------

# --------------FUNCTION DEFINITION-------------------
# Define function to calculate the truncation radius as in Manara+2019 (Appendix C)

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

# -------------- CALL FUNCTION TO CALCULATE THE TRUNCATION RADIUS -------------------

Rt1 = [] #au
Rt2 = [] #au

for i in range(len(data['alpha'])):
    # calculate disc truncation around primary (last argument of function needs to be False)
    Rt1.append(calculate_truncation(M1[i],M2[i],ebin[i],alpha[i],abin[i],False))

    # calculate disc truncation around secondary (last argument of function needs to be True)
    Rt2.append(calculate_truncation(M1[i],M2[i],ebin[i],alpha[i],abin[i],True))

# --------------CALCULATE IRRADIATION FROM SECONDARY-------------------

Tirr2_at_Rt1 = []  
Tirr1_at_Rt2 = []

for i in range(len(data)):

    Mstar = M2[i] * Msun          
    Lstar = Lsun * (Mstar / Msun)**3.5  

    a_bin = abin[i] * au    
    e_bin = ebin[i]      

    Rtab_au = np.linspace(1, 1000, 10000) 
    Rtab_cm = Rtab_au * au

    Nbin = 1000
    M_array = 2 * pi * np.arange(Nbin) / (Nbin - 1)
    E_array = M_array.copy()
    for _ in range(10):
        E_array = M_array + e_bin * np.sin(E_array)

    rbin_array = a_bin * (1 - e_bin * np.cos(E_array))  # cm

    Tirr2 = np.zeros_like(Rtab_cm)

    for j, R in enumerate(Rtab_cm):
        d_array = np.abs(R - rbin_array)
        d_array[d_array < 1e10] = 1e10    
        inv_d2_mean = np.mean(d_array**2)
        Tirr2[j] = (Lstar / (16 * pi * sigma * inv_d2_mean))**0.25

    # Calcolo Tirr2 al luogo di Rt1
    Rt1_cm = Rt1[i] * au
    Tirr_interp = np.interp(Rt1_cm, Rtab_cm, Tirr2)
    Tirr2_at_Rt1.append(Tirr_interp)
    
    Mstar = M1[i] * Msun          
    Lstar = Lsun * (Mstar / Msun)**3.5  

    Nbin = 1000
    M_array = 2 * pi * np.arange(Nbin) / (Nbin - 1)
    E_array = M_array.copy()
    for _ in range(10):
        E_array = M_array + e_bin * np.sin(E_array)

    rbin_array = a_bin * (1 - e_bin * np.cos(E_array))  # cm

    Tirr1 = np.zeros_like(Rtab_cm)

    for j, R in enumerate(Rtab_cm):
        d_array = np.abs(R - rbin_array)
        d_array[d_array < 1e10] = 1e10    
        inv_d2_mean = np.mean(d_array**2)
        Tirr1[j] = (Lstar / (16 * pi * sigma * inv_d2_mean))**0.25

    # Calcolo Tirr2 al luogo di Rt1
    Rt2_cm = Rt2[i] * au
    Tirr_interp = np.interp(Rt2_cm, Rtab_cm, Tirr1)
    Tirr1_at_Rt2.append(Tirr_interp)

    # Plot per ogni sistema (opzionale)
    plt.figure()
    plt.plot(Rtab_au, Tirr2, label="Tirr from secondary on primary")
    plt.vlines(a_bin/au, 0, max(Tirr2), 'red', 'dashed', label="a_bin")
    plt.vlines(Rt1[i], 0, max(Tirr2), 'blue', 'dashed', label="Rt1")
    plt.plot(Rtab_au, Tirr1, label="Tirr from primary on secondary")
    plt.vlines(Rt2[i], 0, max(Tirr1), 'cyan', 'dashed', label="Rt2")
    plt.xlabel("Distance [au]")
    plt.ylabel("Irradiation temperature [K]")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(False)
    plt.legend()
    plt.show()
    
    
# -------------------------

for i in range(len(data)):

    Mstar = M2[i] * Msun          
    Lstar = Lsun * (Mstar / Msun)**3.5  

    a_bin = abin[i] * au    
    e_bin = ebin[i]      

    Rtab_au = np.linspace(0.1, 1000, 10000) 
    Rtab_cm = Rtab_au * au

    Nbin = 1000
    M_array = 2 * pi * np.arange(Nbin) / (Nbin - 1)
    E_array = M_array.copy()
    for _ in range(10):
        E_array = M_array + e_bin * np.sin(E_array)

    rbin_array = a_bin * (1 - e_bin * np.cos(E_array))  # cm

    Tirr2 = np.zeros_like(Rtab_cm)

    for j, R in enumerate(Rtab_cm):
        d_array = np.abs(R - rbin_array)
        d_array[d_array < 1e10] = 1e10    
        inv_d2_mean = np.mean(d_array**2)
        Tirr2[j] = (Lstar / (16 * pi * sigma * inv_d2_mean))**0.25
    
    Mstar = M1[i] * Msun          
    Lstar = Lsun * (Mstar / Msun)**3.5  
    Tirr1 = (Lstar / (16 * pi * sigma * (Rtab_cm**2)))**0.25


    # Plot per ogni sistema (opzionale)
    plt.figure()
    plt.plot(Rtab_au, Tirr1, label="Tirr1")
    plt.gca().invert_xaxis()
    plt.gca().invert_xaxis() 
    plt.vlines(a_bin/au, 0, max(Tirr2), 'red', 'dashed', label="a_bin")
    plt.vlines(Rt1[i], 0, max(Tirr2), 'blue', 'dashed', label="Rt1")
    plt.plot(Rtab_au, Tirr2, label="Tirr2")
    plt.xlabel("Distance [au]")
    plt.ylabel("Irradiation temperature [K]")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(False)
    plt.legend()
    plt.show()    

# --------------END CALCULATE IRRADIATION FROM SECONDARY SECTION-------------------

# -------------- PRINT ON OUTPUT FILE AND SAVE -------------------

file_name = 'output.csv'
data_out = {'alpha': alpha, 
                'Mstar1': M1,
                'Mstar2': M2,
                'a_bin': abin,
                'e_bin': ebin,
                'Rt1': Rt1,
                'Rt2': Rt2,
                'Tirr2_at_Rt1': Tirr2_at_Rt1,
                'Tirr1_at_Rt2': Tirr1_at_Rt2
        }

df_out = pd.DataFrame(data_out)
print(df_out)
df_out.to_csv(file_name, header=True, index=None, sep=',')

# %%

import numpy as np
import pandas as pd
from pathlib import Path as path
import matplotlib.pyplot as plt
data = pd.read_csv('Teff_MG.csv', sep=',', header=None)

coeff = np.polyfit(data[0],data[1],5)
T = np.linspace(min(data[0]),max(data[0]),10000)
M = np.polyval(coeff,T)
#plt.plot(T,M)
#plt.xlim(10000,3000)
#plt.ylim(12,0)

x = np.polyval(coeff,3710)
x_up = np.polyval(coeff,3710+180)
x_low = np.polyval(coeff,3710-180)
Gmag2 = x+5*np.log10(114.046/10)
Gmag2_low = Gmag2 - (x_up+5*np.log10(114.046/10) )
Gmag2_up = - Gmag2 + (x_low+5*np.log10(114.046/10))
print(x,Gmag2,Gmag2_up,Gmag2_low)

Gmag1 = 10.9976	
Gmag_err = 0.0003


f1f2 = 10**((-Gmag1+Gmag2)/2.5)
f2f1 = 10**((Gmag1-Gmag2)/2.5)

print(f1f2,f2f1)
# %%
