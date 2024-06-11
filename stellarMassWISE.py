"""
Functions to estimate stellar mass from W1 and W2 magnitudes (WISE survey)

These functions come from Jarrett, Cluver et al.: https://arxiv.org/pdf/2301.05952
"""
import numpy as np
import matplotlib.pyplot as plt

def getMW1(L_W1):
    """
    Get W1 filter absolute magnitude (in vega mags) from the L_W1 luminosity (in solar luminosities)
    """
    M_SUN = 3.24
    M_W1 = M_SUN-2.5*np.log10(L_W1)
    return M_W1

def getLW1(M_W1):
    """
    Get L_W1 luminosity (in solar luminosities) from the W1 filter absolute magnitude (in vega mags)
    """
    M_SUN = 3.24
    L_W1 = 10**(-0.4*(M_W1-M_SUN))
    return L_W1

def getStarMass_W1(M_W1, returnLog=False):
    """
    Get Stellar Mass (in solar Masses) from the W1 filter absolute magnitude (in vega mags)
    """
    A = (-12.62185, 5.00155, -0.43857, 0.01593)
    L_W1 = getLW1(M_W1)
    logStarMass = A[0]+A[1]*np.log10(L_W1)+A[2]*np.log10(L_W1)**2+A[3]*np.log10(L_W1)**3
    if returnLog:
        return logStarMass
    return 10**(logStarMass)

def getStarMass_W1_W2(M_W1, M_W2, returnLog=False):
    """
    Get Stellar Mass (in solar Masses) from the W1 and W2 filters absolute magnitudes (in vega mags)
    """
    L_W1 = getLW1(M_W1)
    C12 = M_W1-M_W2
    A = (-0.376, -1.053)
    logMassRatio = A[0]+A[1]*C12
    logStarMass = logMassRatio + np.log10(L_W1)
    if returnLog:
        return logStarMass
    return 10**(logStarMass)



if __name__ == "__main__":
    lw1s = 10**(np.linspace(7,12,100))
    mw1s = getMW1(lw1s)
    plt.plot(np.log10(lw1s), getStarMass_W1(mw1s, returnLog=True), "--")
    plt.xlabel(r"$\log(L_\mathrm{W1})$ [$L_\odot$]", fontsize=16)
    plt.ylabel(r"$\log(M_\star)$ [$M_\odot$]", fontsize=16)
    plt.show()