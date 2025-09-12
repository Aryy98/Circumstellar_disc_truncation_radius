# Circumstellar Disc Truncation Radius Calculator

This Python script computes the truncation radius of a circumstellar protoplanetary disc in an binary system using the recipe from *Manara et al. (2019)* as described in *Venturini et al. (submitted)*

---

## Authors

- **Main developer:** Arianna Nigioni  
- **Co-developer:** Julia Venturini  
- **Date:** September 2025  

---

## Overview

The script calculates the truncation radius of a protoplanetary disc around a star in a binary system using:

- Binary parameters: primary ($M_1$) and secondary ($M_2$) star masses (Note: $M_1$ > $M_2$), binary separation, and binary eccentricity.  
- Disc viscosity: parameterized by the viscosity $\alpha$
- Option to compute the truncation radius for either the **circumprimary** or **circumsecondary** disc.

The truncation radius is computed following the method outlined in *Manara et al. (2019, Appendix C)*, with adjustments described in *Venturini et al. (submitted)*.

---

## Requirements
 
- Python packages:  
  - `numpy`  
  - `pandas`
 
## Usage

- Define the binary system parameters (Example provided: Gamma-Cephei system):
  - `M1 = 1.27      # Mass of the primary star [Msun]`
  - M2 = 0.328     # Mass of the secondary star [Msun]
  - abin = 19.56   # Binary separation [au]
  - ebin = 0.41    # Binary eccentricity
  - Stype_truncation_secondary = False  # True for circumsecondary disc
  - alpha = 0.001  # Disc alpha viscosity


- Call the function:
  - Rt = calculate_truncation(M1, M2, ebin, alpha, abin, Stype_truncation_secondary)
  The function returns the outer truncation radius in astronomical units (au).

- Resulting truncation radius value is printed on screen

## Function description 

The function `calculate_truncation` computes the truncation radius of a circumstellar protoplanetary disc in an binary system using the recipe from *Manara et al. (2019)* as described in *Venturini et al. (submitted)*. 
The truncation radius is expressed as:

$R_{\rm trunc}(M_1, M_2, e_{\rm bin}, a_{\rm bin}) = R_{\rm Egg} \times (b \, e_{\rm bin}^c + h \, \mu^k)$

where \(b, c, h, k\) are fitting parameters, and  

$\frac{R_{\rm Egg}}{a_{\rm bin}} = \frac{0.49 \, q^{-2/3}}{0.6 \, q^{-2/3} + \ln(1+q^{-1/3})}$

with \(q = M_2/M_1\) (\(M_2 < M_1\)) and \(\mu = M_2/(M_1 + M_2)\) *(Eggleton 1983)*. The parameters \(h = 0.88\) and \(k = 0.01\) are found by fitting the data from *Papaloizou & Pringle (1977)*.  

The coefficients \(b\) and \(c\) depend on the mass parameter μ and the Reynolds number \(\mathcal{R}\), with separate tables for circumprimary and circumsecondary discs (see Table C.1 in *Manara et al. 2019*, which fits the data from *Artymowicz & Lubow (1994)*). 
To compute \(b\) and \(c\) for arbitrary \(\mathcal{R}\) and μ, the function:

1. Computes the Reynolds number from the disc viscosity α and a fixed disc aspect ratio \(H/r = 0.0983\), derived from synthetic populations of single-star discs simulated with the Bern Model *(Emsenhuber et al. 2021)*. The aspect ratio is measured beyond 10 au and averaged across 1000 discs for α = 10⁻³ and 10⁻⁴.  
2. Performs a **two-step interpolation**:
   - Fits second-degree polynomials in \(\log_{10}(\mathcal{R})\) to the tabulated \(b\) and \(c\) values for each μ.  
   - Evaluates these polynomials at the computed Reynolds number and then fits fourth-degree polynomials in μ to obtain the final coefficients at the desired mass parameter.  

This procedure allows the model to calculate the truncation radius for either the **circumprimary** or **circumsecondary** disc, depending on the choice of the user.

- Input parameters:
  - Mstar: Mass of the primary star [Msun]
  - Mstar2: Mass of the secondary star [Msun]
  - e_star2: Binary eccentricity
  - alpha_disk: Disc alpha viscosity
  - a_star2: Binary separation [au]
  - Stype_truncation_secondary: Boolean, True for circumsecondary disc, False for circumprimary disc
- Output: Truncation radius in au.

## References
- Artymowicz, P. & Lubow, S. H. 1994, ApJ, 421, 651
- Eggleton, P. P. 1983, ApJ, 268, 368
- Emsenhuber, A., Mordasini, C., Burn, R., et al. 2021, A&A, 656, A69
- Manara, C. F., Tazzari, M., Long, F., et al. 2019, A&A, 628, A95
- Papaloizou, J. & Pringle, J. E. 1977, MNRAS, 181, 441
- Venturini, J., Nigioni, A., t al. submitted.
