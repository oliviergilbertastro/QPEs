"""
Program to compare QPE hosts and TDE hosts simply, but this compares the QPE hosts properties measuredfrom uncalibrated files,
a new version of this program will be made once the griz bands have all been fitted with galight
"""
import pickle
import numpy as np
from galight_modif.tools.plot_tools import total_compare
from download_data import objects
from ned_wright_cosmology import calculate_cosmo
from utils import print_table
import matplotlib.pyplot as plt
import corner
from scipy.stats import norm
from paper_data import *
from download_data import *#objects, comparisons, objects_names, comparisons_names, objects_types, QPE_bands_list
from stellarMassblackHoleMass import stellarMass_mBH, log10_stellarMass_mBH
from stellarMassWISE import TDE_stellarMasses_WISE, QPE_stellarMasses_WISE
from compareModels import plot_sersicIndex_mBH, plot_surfaceStellarMassDensity_mBH, plot_sersicIndex_surfaceStellarMassDensity

def get_QPE_n_and_r50(ra_dec, model="None", band="i"):
    picklename = f'ra{str(ra_dec[0])}_dec{str(ra_dec[1])}_{model}_{band}.pkl'
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        raise NameError("There is no picklerun with this name.")

    #Calculate the Sersic index + uncertainties
    chain = fitting_run_result.samples_mcmc
    lo, mid, hi = np.percentile(chain[:, 1],16), np.percentile(chain[:, 1],50), np.percentile(chain[:, 1],84)
    plus, minus = (hi-mid), (mid-lo)
    sersic_index_data = [mid, minus, plus]

    #Calculate the Sersic half-light radius + uncertainties:
    lo, mid, hi = np.percentile(chain[:, 0],16), np.percentile(chain[:, 0],50), np.percentile(chain[:, 0],84)
    plus, minus = (hi-mid), (mid-lo)
    r50_data = [mid, minus, plus]



    return sersic_index_data, r50_data


def chooseStellarMassMethod():
    sMass_option = input("Which stellar mass method do you want to use?\n    1. Litterature\n    2. Litterature + SDSS prospector\n    3. BH mass relation\n    4. WISE mass relation\n")
    try:
        sMass_option = int(sMass_option)
    except:
        sMass_option = 1
    
    if sMass_option == 1:
        QPE_stellar_masses = QPE_stellar_masses_litterature
        TDE_stellar_masses = TDE_stellar_masses_litterature
    elif sMass_option == 2:
        QPE_stellar_masses = QPE_stellar_masses_litterature_sdssProspector
        TDE_stellar_masses = TDE_stellar_masses_litterature
    elif sMass_option == 3:
        QPE_stellar_masses = stellarMass_mBH(QPE_mBH)
        QPE_stellar_masses[:,1:] = np.zeros_like(QPE_stellar_masses[:,1:]) #Make all uncertainties zero as they are currently not calculated properly
        QPE_stellar_masses = list(QPE_stellar_masses)
        for i in range(len(QPE_stellar_masses)):
            QPE_stellar_masses[i] = tuple(QPE_stellar_masses[i])
        TDE_stellar_masses = stellarMass_mBH(TDE_mBH)
        TDE_stellar_masses = add_0_uncertainties(TDE_stellar_masses)
        TDE_stellar_masses = list(TDE_stellar_masses)
        for i in range(len(TDE_stellar_masses)):
            TDE_stellar_masses[i] = tuple(TDE_stellar_masses[i])
    elif sMass_option == 4:
        QPE_stellar_masses = QPE_stellarMasses_WISE
        QPE_stellar_masses = list(QPE_stellar_masses)
        for i in range(len(QPE_stellar_masses)):
            QPE_stellar_masses[i] = tuple(QPE_stellar_masses[i])
        TDE_stellar_masses = TDE_stellarMasses_WISE
        TDE_stellar_masses = list(TDE_stellar_masses)
        for i in range(len(TDE_stellar_masses)):
            TDE_stellar_masses[i] = tuple(TDE_stellar_masses[i])
    else:
        raise ValueError("This is not a valide option.")
    return QPE_stellar_masses, TDE_stellar_masses

def add_0_uncertainties(a):
    a = np.array(a)
    placeholder = np.zeros((a.shape[0],3))
    placeholder[:,0] = a
    a = placeholder
    return a

def stellarMassDensity(M_star, r50):
    '''
    Calculate the stellar surface mass density \Sigma_{M_\ast} from the total stellar mass and the half-light radius

    M_star: stellar mass in solar masses -> tuple/list/array such as (mass, errlo, errhi)
    r50: half-light radius in kpc -> tuple/list/array such as (radius, errlo, errhi)

    returns [value, errlo, errhi]
    '''
    try:
        res = M_star[0]/r50[0]**2
        errlo = res*np.sqrt((M_star[1]/M_star[0])**2+(2*r50[2]/r50[0])**2)
        errhi = res*np.sqrt((M_star[2]/M_star[0])**2+(2*r50[1]/r50[0])**2)
        return [res, errlo, errhi]
    except:
        #In the exception where the user did not input uncertainties, the value will still be calculated.
        return M_star/r50[0]**2#/(2*np.pi) #is there a constant 2pi we need to divide by???



#filter bands used to fit each object to simplify the user input, some other bands were also fitted and/or could be fitted, but these ones work well



if __name__ == "__main__":

    #First, load the TDE sersic indices and half-light radii into arrays or list, idc:
    QPE_sersicIndices = []
    QPE_r50s = []
    for i in range(len(objects)):
        n, r50 = get_QPE_n_and_r50(objects[i], objects_types[i], QPE_bands_list[i])
        QPE_sersicIndices.append(n)
        QPE_r50s.append(r50)
    QPE_stellar_masses, TDE_stellar_masses = chooseStellarMassMethod()
    #Transform lists into arrays
    QPE_sersicIndices, QPE_r50s, QPE_stellar_masses, QPE_mBH, TDE_stellar_masses, TDE_r50s, TDE_mBH = np.array(QPE_sersicIndices), np.array(QPE_r50s), np.array(QPE_stellar_masses), np.array(QPE_mBH), np.array(TDE_stellar_masses), np.array(TDE_r50s), np.array(TDE_mBH)
    #Transfer the QPE r50s from arcseconds to kpc:
    for i in range(len(objects)):
        cosmology_params = calculate_cosmo(QPE_redshifts[i], H0=67, Omega_m=0.31, Omega_v=0.69)
        kpc_per_arcsec = cosmology_params["kpc_DA"]
        QPE_r50s[i] = QPE_r50s[i]*kpc_per_arcsec
    #Add +/-0 uncertainties to TDE r50s
    TDE_r50s = add_0_uncertainties(TDE_r50s)
    #Create upper-limit list for TDE stellar surface densities:
    TDE_uplims = [False, False, False, False, False, False, False, False, True, False]
    #Calculate stellar mass surface densities
    QPE_stellarDensities = []
    for i in range(len(objects)):
        QPE_stellarDensities.append(stellarMassDensity(QPE_stellar_masses[i], QPE_r50s[i]))
        QPE_stellarDensities[i] = (np.log10(QPE_stellarDensities[i][0]), QPE_stellarDensities[i][1]/(QPE_stellarDensities[i][0]*np.log(10)), QPE_stellarDensities[i][2]/(QPE_stellarDensities[i][0]*np.log(10)))
    QPE_stellarDensities = np.array(QPE_stellarDensities)

    TDE_stellarDensities = []
    for i in range(len(TDE_names)):
        TDE_stellarDensities.append(stellarMassDensity(TDE_stellar_masses[i], TDE_r50s[i]))
        TDE_stellarDensities[i] = (np.log10(TDE_stellarDensities[i][0]), TDE_stellarDensities[i][1]/(TDE_stellarDensities[i][0]*np.log(10)), TDE_stellarDensities[i][2]/(TDE_stellarDensities[i][0]*np.log(10)))
    TDE_stellarDensities = np.array(TDE_stellarDensities)

    

    # QPE hosts properties
    print_table(np.array([objects_names, np.around(QPE_sersicIndices[:,0], 4), np.around(QPE_r50s[:,0], 4), np.around(np.log10(QPE_stellar_masses[:,0]), 4), np.around(np.log10(QPE_mBH[:,0]), 4), np.around(QPE_stellarDensities[:,0], 4)]).T,
                    header=["Name", "Sérsic index", "r50 (kpc)", "log Stellar mass", "log M_BH", "log Surface density"],
                    title="QPE hosts properties",
                    space_between_columns=4,
                    space_between_rows=0,
                    borders=2)
    
    # TDE hosts properties
    print_table(np.array([TDE_names, np.around(TDE_sersicIndices, 4), np.around(TDE_r50s[:,0], 4), np.around(np.log10(TDE_stellar_masses[:,0]), 4), np.around(np.log10(TDE_mBH), 4), np.around(TDE_stellarDensities[:,0], 4)]).T,
                    header=["Name", "Sérsic index", "r50 (kpc)", "log Stellar mass", "log M_BH", "log Surface density"],
                    title="TDE hosts properties",
                    space_between_columns=4,
                    space_between_rows=0,
                    borders=2)
    

    #Plots:
    plot_sersicIndex_mBH(QPE_mBH, QPE_sersicIndices, TDE_mBH, TDE_sersicIndices)
    plot_surfaceStellarMassDensity_mBH(QPE_mBH, QPE_stellarDensities, TDE_mBH, TDE_stellarDensities, uplims=TDE_uplims)
    plot_surfaceStellarMassDensity_mBH(QPE_stellar_masses, QPE_stellarDensities, TDE_stellar_masses, TDE_stellarDensities, uplims=TDE_uplims)
    plot_sersicIndex_surfaceStellarMassDensity(QPE_sersicIndices, QPE_stellarDensities, TDE_sersicIndices, TDE_stellarDensities)

    