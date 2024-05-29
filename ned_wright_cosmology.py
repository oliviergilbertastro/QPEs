#!/usr/bin/env python
  
import sys
from math import *

def calculate_cosmo(z, H0=None, Omega_m=None, Omega_v=None, verbose=False):
    """
    Based on Ned Wright's online cosmology calculator: https://astro.ucla.edu/~wright/CosmoCalc.html

    Calculates a bunch of parameters from a redshift and an input cosmology

    z: Redshift [0,+inf]
    H0: Hubble constant (km/s)
    Omega_m: Fraction of matter [0,1]
    Omega_v: Fraction of vacuum [0,1], Omega_m+Omega_v <= 1
    verbose: Boolean to decide if parameters are printed out or not

    Returns: dictionnary of parameters, to check what each parameter is, use verbose=True
    """
    # if no values, assume standard cosmology Model, input is z
    if H0==None and Omega_m==None and Omega_v==None:
        z = z                         # redshift
        H0 = 67                         # Hubble constant
        WM = 0.31                       # Omega(matter)
        WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda

    # if one value, assume Benchmark Model with given Ho
    elif Omega_m==None and Omega_v==None:
        z = z                         # redshift
        H0 = H0                         # Hubble constant
        WM = 0.31                       # Omega(matter)
        WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda

    # if Univ is Open, use Ho, Wm and set Wv to 0.
    elif Omega_v==None:
        z = z                         # redshift
        H0 = H0                         # Hubble constant
        WM = Omega_m                    # Omega(matter)
        WV = 0.0                        # Omega(vacuum) or lambda

    # if Univ is General, use Ho, Wm and given Wv
    else:
        z = z                         # redshift
        H0 = H0                         # Hubble constant
        WM = Omega_m                    # Omega(matter)
        WV = Omega_v                    # Omega(vacuum) or lambda


    # initialize constants

    WR = 0.        # Omega(radiation)
    WK = 0.        # Omega curvaturve = 1-Omega(total)
    c = 299792.458 # velocity of light in km/sec
    Tyr = 977.8    # coefficent for converting 1/H into Gyr
    DTT = 0.5      # time from z to now in units of 1/H0
    DTT_Gyr = 0.0  # value of DTT in Gyr
    age = 0.5      # age of Universe in units of 1/H0
    age_Gyr = 0.0  # value of age in Gyr
    zage = 0.1     # age of Universe at redshift z in units of 1/H0
    zage_Gyr = 0.0 # value of zage in Gyr
    DCMR = 0.0     # comoving radial distance in units of c/H0
    DCMR_Mpc = 0.0 
    DCMR_Gyr = 0.0
    DA = 0.0       # angular size distance
    DA_Mpc = 0.0
    DA_Gyr = 0.0
    kpc_DA = 0.0
    DL = 0.0       # luminosity distance
    DL_Mpc = 0.0
    DL_Gyr = 0.0   # DL in units of billions of light years
    V_Gpc = 0.0
    a = 1.0        # 1/(1+z), the scale factor of the Universe
    az = 0.5       # 1/(1+z(object))

    h = H0/100.
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0.
    n=1000         # number of points in integrals
    for i in range(n):
        a = az*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot

    zage = az*age/n
    zage_Gyr = (Tyr/H0)*zage
    DTT = 0.0
    DCMR = 0.0

    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)

    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage
    age_Gyr = age*(Tyr/H0)
    DTT_Gyr = (Tyr/H0)*DTT
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR

    # tangential comoving distance

    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio =  0.5*(exp(x)-exp(-x))/x 
        else:
            ratio = sin(x)/x
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/6. + y*y/120.
    DCMT = ratio*DCMR
    DA = az*DCMT
    DA_Mpc = (c/H0)*DA
    kpc_DA = DA_Mpc/206.264806
    DA_Gyr = (Tyr/H0)*DA
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    DL_Gyr = (Tyr/H0)*DL

    # comoving volume computation

    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
        else:
            ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/5. + (2./105.)*y*y
    VCM = ratio*DCMR*DCMR*DCMR/3.
    V_Gpc = 4.*pi*((0.001*c/H0)**3)*VCM

    if verbose:
        print('For H_o = ' + '%1.1f' % H0 + ', Omega_M = ' + '%1.2f' % WM + ', Omega_vac = ',)
        print('%1.2f' % WV + ', z = ' + '%1.3f' % z)
        print('It is now ' + '%1.3f' % age_Gyr + ' Gyr since the Big Bang.')
        print('The age at redshift z was ' + '%1.3f' % zage_Gyr + ' Gyr.')
        print('The light travel time was ' + '%1.3f' % DTT_Gyr + ' Gyr.')
        print('The comoving radial distance, which goes into Hubbles law, is',)
        print('%1.3f' % DCMR_Mpc + ' Mpc or ' + '%1.3f' % DCMR_Gyr + ' Gly.')
        print('The comoving volume within redshift z is ' + '%1.3f' % V_Gpc + ' Gpc^3.')
        print('The angular size distance D_A is ' + '%1.3f' % DA_Mpc + ' Mpc or',)
        print('%1.3f' % DA_Gyr + ' Gly.')
        print('This gives a scale of ' + '%.3f' % kpc_DA + ' kpc/".')
        print('The luminosity distance D_L is ' + '%1.3f' % DL_Mpc + ' Mpc or ' + '%1.3f' % DL_Gyr + ' Gly.')
        print('The distance modulus, m-M, is '+'%1.3f' % (5*log10(DL_Mpc*1e6)-5))
    return {"age_Gyr": age_Gyr,
            "zage_Gyr": zage_Gyr,
            "DTT_Gyr": DTT_Gyr,
            "DCMR_Mpc": DCMR_Mpc,
            "DCMR_Gyr": DCMR_Gyr,
            "V_Gpc": V_Gpc,
            "DA_Mpc": DA_Mpc,
            "DA_Gyr": DA_Gyr,
            "kpc_DA": kpc_DA,
            "DL_Mpc": DL_Mpc,
            "DL_Gyr": DL_Gyr,
            "D_mod": (5*log10(DL_Mpc*1e6)-5)
            }


#age_Gyr, zage_Gyr, DTT_Gyr, DCMR_Mpc, DCMR_Gyr, V_Gpc, DA_Mpc, DA_Gyr, kpc_DA, DL_Mpc, DL_Gyr, (5*log10(DL_Mpc*1e6)-5)

if __name__ == "__main__":
  res = calculate_cosmo(z=0.246, H0=67, Omega_m=0.31, Omega_v=0.69, verbose=True)