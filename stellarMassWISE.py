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

def getLW1(M_W1, sigma_M_W1=None):
    """
    Get L_W1 luminosity (in solar luminosities) from the W1 filter absolute magnitude (in vega mags)

    M_W1: W1 filter VEGA magnitude
    sigma_M_W1: (optional) uncertainty on VEGA magnitude
    """
    M_SUN = 3.24
    L_W1 = 10**(-0.4*(M_W1-M_SUN))
    if sigma_M_W1 != None:
        sigma_L_W1 = L_W1*np.abs(0.4*np.log(10)*sigma_M_W1)
        return L_W1, sigma_M_W1
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

def getStarMass_W1_W2(M_W1, M_W2, sigma_M_W1=None, sigma_M_W2=None, returnLog=False):
    """
    Get Stellar Mass (in solar Masses) from the W1 and W2 filters absolute magnitudes (in vega mags)
    """
    L_W1 = getLW1(M_W1, sigma_M_W1=sigma_M_W1)
    if sigma_M_W1 != None:
        #Unpack the uncertainty
        L_W1, sigma_L_W1 = L_W1
    C12 = M_W1-M_W2
    A = (-0.376, -1.053)
    logMassRatio = A[0]+A[1]*C12
    logStarMass = logMassRatio + np.log10(L_W1)
    if sigma_M_W1 != None and sigma_M_W2 != None:
        #Calculate uncertainties
        sigma_C12 = np.sqrt(sigma_M_W1**2+sigma_M_W2**2)
        sigma_logMassRatio = np.abs(A[1])*sigma_C12
        sigma_log_L_W1 = np.abs(sigma_L_W1/(L_W1*np.log(10)))
        sigma_logStarMass = np.sqrt(sigma_logMassRatio**2+sigma_log_L_W1**2)
        if returnLog:
            return logStarMass, sigma_logStarMass
        return 10**(logStarMass), 10**(logStarMass)*np.log(10)*sigma_logStarMass

    if returnLog:
        return logStarMass
    return 10**(logStarMass)





#WISE DATA FOR EACH OBJECT:
from download_data import objects_names
from utils import print_table

#Dictionnary of magnitudes (and their uncertainties)
wise_mags = {
    # QPE hosts
    "GSN 069":          {"W1":12.838,"W2":12.828},
    "RX J1301.9+2747":  {"W1":12.580,"W2":12.558},
    "eRO-QPE1":         {"W1":15.327,"W2":15.266},
    "eRO-QPE2":         {"W1":13.674,"W2":13.518},
    "AT 2019vcb":       {"W1":15.778,"W2":15.525},
    "2MASX J0249":      {"W1":13.753,"W2":13.700},
    "eRO-QPE3":         {"W1":15.075,"W2":15.393},
    "eRO-QPE4":         {"W1":13.989,"W2":13.922},
    "AT 2019qiz":       {"W1":12.481,"W2":12.521},

    # TDE hosts
    "ASASSN-14ae":      {"W1":14.296,"W2":14.150},
    "ASASSN-14li":      {"W1":13.061,"W2":13.049},
    "PTF-09ge":         {"W1":14.423,"W2":14.295},
    "RBS 1032":         {"W1":14.471,"W2":14.382},
    "SDSS J1323":       {"W1":14.803,"W2":14.729},
    "SDSS J0748":       {"W1":13.505,"W2":12.728},
    "SDSS J1342":       {"W1":13.428,"W2":12.354},
    "SDSS J1350":       {"W1":13.925,"W2":13.062},
    "SDSS J0952":       {"W1":13.666,"W2":12.638},
    "SDSS J1201":       {"W1":14.895,"W2":14.677},
    }

sigma_wise_mags = {
    # QPE hosts
    "GSN 069":          {"W1":0.023,"W2":0.026},
    "RX J1301.9+2747":  {"W1":0.023,"W2":0.024},
    "eRO-QPE1":         {"W1":0.036,"W2":0.083},
    "eRO-QPE2":         {"W1":0.025,"W2":0.028},
    "AT 2019vcb":       {"W1":0.046,"W2":0.108},
    "2MASX J0249":      {"W1":0.026,"W2":0.033},
    "eRO-QPE3":         {"W1":0.036,"W2":0.096},
    "eRO-QPE4":         {"W1":0.026,"W2":0.036},
    "AT 2019qiz":       {"W1":0.024,"W2":0.024},

    # TDE hosts
    "ASASSN-14ae":      {"W1":0.027,"W2":0.046},
    "ASASSN-14li":      {"W1":0.024,"W2":0.027},
    "PTF-09ge":         {"W1":0.026,"W2":0.036},
    "RBS 1032":         {"W1":0.028,"W2":0.046},
    "SDSS J1323":       {"W1":0.031,"W2":0.053},
    "SDSS J0748":       {"W1":0.025,"W2":0.026},
    "SDSS J1342":       {"W1":0.025,"W2":0.024},
    "SDSS J1350":       {"W1":0.025,"W2":0.024},
    "SDSS J0952":       {"W1":0.026,"W2":0.026},
    "SDSS J1201":       {"W1":0.033,"W2":0.056},
    }

    

#Get apparent magnitudes:
w1_mags = []
w2_mags = []
w1_mags_unc = []
w2_mags_unc = []
for i in range(len(objects_names)):
    w1_mags.append(wise_mags[objects_names[i]]["W1"])
    w2_mags.append(wise_mags[objects_names[i]]["W2"])
    w1_mags_unc.append(sigma_wise_mags[objects_names[i]]["W1"])
    w2_mags_unc.append(sigma_wise_mags[objects_names[i]]["W2"])
#w1_mags, w2_mags, w1_mags_unc, w2_mags_unc = np.array(w1_mags), np.array(w2_mags), np.array(w1_mags_unc), np.array(w2_mags_unc)

#Get absolute magnitudes:
from ned_wright_cosmology import calculate_cosmo
from paper_data import QPE_redshifts
distances = []
w1_abs_mags = []
w2_abs_mags = []

for i in range(len(objects_names)):
    distances.append(calculate_cosmo(QPE_redshifts[i])["DL_Mpc"])
    w1_abs_mags.append(w1_mags[i]-2.5*np.log10((distances[i]*1E6/10)**2))
    w2_abs_mags.append(w2_mags[i]-2.5*np.log10((distances[i]*1E6/10)**2))


#Calculate stellar masses:
QPE_stellarMasses_WISE = []
for i in range(len(objects_names)):
    sm = getStarMass_W1_W2(w1_abs_mags[i], w2_abs_mags[i], w1_mags_unc[i], w2_mags_unc[i])
    QPE_stellarMasses_WISE.append((sm[0], sm[1], sm[1])) #Setting 0 uncertainty for the moment
QPE_stellarMasses_WISE = np.array(QPE_stellarMasses_WISE)


#Calculate stellar masses:
TDE_stellarMasses_WISE = []
for i in range(len(objects_names)):
    sm = getStarMass_W1_W2(w1_abs_mags[i], w2_abs_mags[i], w1_mags_unc[i], w2_mags_unc[i])
    TDE_stellarMasses_WISE.append((sm[0], sm[1], sm[1])) #Setting 0 uncertainty for the moment
TDE_stellarMasses_WISE = np.array(TDE_stellarMasses_WISE)

if __name__ == "__main__":
    print_table(np.array([objects_names, np.around(distances, 2), np.around(w1_abs_mags, 2), w1_mags_unc, np.around(w2_abs_mags, 2), w2_mags_unc, np.around(np.log10(QPE_stellarMasses_WISE[:,0]), 3)]).T, ["Name", "Distance (Mpc)", "W1", "+/-", "W2", "+/-", "log Stellar Mass (M_sun)"], title="QPE hosts WISE properties", borders=2)
    if False:
        lw1s = 10**(np.linspace(7,12,100))
        mw1s = getMW1(lw1s)
        print(mw1s)
        plt.plot(np.log10(lw1s), getStarMass_W1(mw1s, returnLog=True), "--")
        plt.xlabel(r"$\log(L_\mathrm{W1})$ [$L_\odot$]", fontsize=16)
        plt.ylabel(r"$\log(M_\star)$ [$M_\odot$]", fontsize=16)
        plt.show()