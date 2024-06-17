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
from compareModels import plot_sersicIndex_mBH, plot_surfaceStellarMassDensity_mBH, plot_sersicIndex_surfaceStellarMassDensity, plot_surfaceStellarMassDensity_stellarMass
import sys
from prospector.prospector_myFits import QPE_stellar_masses_desiProspector

def get_QPE_n_and_r50(ra_dec, model="None", band="i", survey="DESI", redshift=0):
    objID = objects.index(ra_dec)
    picklename = f"{objects_names[objID]}_{band}-band_{model}_{survey}.pkl"
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        raise NameError("There is no picklerun with this name.")

    #Calculate the Sersic index + uncertainties
    chain = fitting_run_result.samples_mcmc
    params = fitting_run_result.param_mcmc
    if "R_sersic" in params and "n_sersic" in params:
        lo, mid, hi = np.percentile(chain[:, 1],16), np.percentile(chain[:, 1],50), np.percentile(chain[:, 1],84)
        plus, minus = (hi-mid), (mid-lo)
        sersic_index_data = [mid, minus, plus]

        #Calculate the Sersic half-light radius + uncertainties:
        lo, mid, hi = np.percentile(chain[:, 0],16), np.percentile(chain[:, 0],50), np.percentile(chain[:, 0],84)
        plus, minus = (hi-mid), (mid-lo)
        r50_data = [mid, minus, plus]
    else:
        sersic_index_data = [fitting_run_result.final_result_galaxy[0]["n_sersic"], 0, 0]
        r50_data = [fitting_run_result.final_result_galaxy[0]["R_sersic"], 0, 0]
    #Convert r50 from arcsec to kpc
    cosmology_params = calculate_cosmo(redshift, H0=67, Omega_m=0.31, Omega_v=0.69)
    kpc_per_arcsec = cosmology_params["kpc_DA"]
    r50_data = np.array(r50_data)*kpc_per_arcsec

    magnitude = fitting_run_result.final_result_galaxy[0]['magnitude']

    return sersic_index_data, r50_data, magnitude


def chooseStellarMassMethod():
    sMass_option = input("Which stellar mass method do you want to use?\n    1. Litterature\n    2. Litterature + SDSS prospector\n    3. BH mass relation\n    4. WISE mass relation\n    5. DESI prospector\n")
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
    elif sMass_option == 5:
        QPE_stellar_masses = QPE_stellar_masses_desiProspector
        QPE_stellar_masses = list(QPE_stellar_masses)
        for i in range(len(QPE_stellar_masses)):
            QPE_stellar_masses[i] = tuple(QPE_stellar_masses[i])
        TDE_stellar_masses = TDE_stellar_masses_litterature
    else:
        raise ValueError("This is not a valide option.")
    return QPE_stellar_masses, TDE_stellar_masses

def add_0_uncertainties(a):
    a = np.array(a)
    placeholder = np.zeros((a.shape[0],3))
    placeholder[:,0] = a
    a = placeholder
    return a

def stellarMassDensity(M_star, r50, returnLog=False):
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
        if returnLog:
            return [np.log10(res), errlo/(res*np.log(10)), errhi/(res*np.log(10))]
        return [res, errlo, errhi]
    except:
        #In the exception where the user did not input uncertainties, the value will still be calculated.
        if returnLog:
            return np.log10(M_star/r50[0]**2)
        return M_star/r50[0]**2#/(2*np.pi) #is there a constant 2pi we need to divide by???


def checkWhichFiltersWork(list_of_dicts):
    working = {"g":[], "r":[], "i":[], "z":[]}
    for i in range(len(list_of_dicts)):
        for band in "griz":
            try:
                idc = list_of_dicts[i][band]
                working[band].append("\x1b[32mY\x1b[0m")
            except:
                working[band].append("\x1b[31mN\x1b[0m")
    print_table(np.array([objects_names, working["g"], working["r"], working["i"], working["z"]]).T,
                header=["Name", "g", "r", "i", "z"],
                title="Working filters",
                borders=2,
                override_length=[15, 1, 1, 1, 1],
                )
    return

def printPropertyAcrossFilters(list_of_dicts, name_of_property="Name of property", round_to_n_decimals=2):
    properties = {"g":[], "r":[], "i":[], "z":[]}
    for i in range(len(list_of_dicts)):
        for band in "griz":
            try:
                properties[band].append(f"{list_of_dicts[i][band][0]:0.2f}")
            except:
                try:
                    properties[band].append(f"{list_of_dicts[i][band]:0.2f}")
                except:
                    properties[band].append("-")
    print_table(np.array([objects_names, properties["g"], properties["r"], properties["i"], properties["z"]]).T,
                header=["Name", "g", "r", "i", "z"],
                title=name_of_property,
                borders=2,
                )
    return


def makeFinalCornerPlot(QPE_mBH, QPE_sersicIndices, QPE_stellarSurfaceDensities, TDE_mBH, TDE_sersicIndices, TDE_stellarSurfaceDensities):
    import seaborn as sns
    QPE_mBH = np.log10(QPE_mBH)
    TDE_mBH = np.log10(TDE_mBH)

    QPE_data = np.array([QPE_mBH, QPE_sersicIndices, QPE_stellarSurfaceDensities])
    TDE_data = np.array([TDE_mBH, TDE_sersicIndices, TDE_stellarSurfaceDensities])

    if False:
        #Make gaussian distributions of all parameters:
        mus, stds = [], []
        QPE_xaxis, TDE_xaxis = [], []
        QPE_pdfs, TDE_pdfs = [], []
        for i in range(len(QPE_data)):

            mu1, std1 = norm.fit(QPEsersicIndices[:,0])
            mu2, std2 = norm.fit(TDEsersicIndices)
            alldata = np.concatenate((QPEsersicIndices[:,0], TDEsersicIndices))
            mini, maxi = np.min(alldata), np.max(alldata)
            # Plot the histogram.
            #plt.hist(QPEsersicIndices[:,0], bins=25, density=False, alpha=0.5, color='b', label='QPE')
            #plt.hist(TDEsersicIndices, bins=25, density=False, alpha=0.5, color='r', label='TDE')
            plt.hist(QPEsersicIndices[:,0], range=(mini, maxi), bins=50, density=False, alpha=0.4, color='b')
            plt.hist(TDEsersicIndices, range=(mini, maxi), bins=50, density=False, alpha=0.4, color='r')
            # Plot the PDF.
            xmin, xmax = plt.xlim()
            x = np.linspace(xmin, xmax, 100)
            p1 = norm.pdf(x, mu1, std1)

    QPE_data, TDE_data = QPE_data.T, TDE_data.T
    fig = corner.corner(QPE_data, color="blue")
    corner.corner(TDE_data, fig=fig, color="red")
    plt.show()
    
    #sns.pairplot(pd.DataFrame([QPE_data, TDE_data]), hue="Hosts")
    return


















survey = "DESI"
if __name__ == "__main__":
    #survey = input("Which survey? [default=DESI]\n")
    survey = "DESI" if survey == "" else survey
#Do this part so other programs can load it (especially the magnitudes from prospector)
#First, load the TDE sersic indices and half-light radii into arrays or list, idc:
QPE_sersicIndices = []
QPE_r50s = []
QPE_magnitudes = []
QPE_unreddenedMagnitudes = []
for i in range(len(objects)):
    QPE_sersicIndices.append({})
    QPE_r50s.append({})
    QPE_magnitudes.append({})
    QPE_unreddenedMagnitudes.append({})
    for band in "griz":
        try:
            n, r50, mag = get_QPE_n_and_r50(objects[i], objects_types[i], redshift=QPE_redshifts[i], band=band, survey=survey)
            QPE_sersicIndices[-1][band] = n
            QPE_r50s[-1][band] = r50
            QPE_magnitudes[-1][band] = mag
            QPE_unreddenedMagnitudes[-1][band] = mag - QPE_extinction[objects_names[i]][band]
        except:
            pass

if __name__ == "__main__":
    #print(objects[i], objects_types[i], band)
    checkWhichFiltersWork(QPE_sersicIndices)
    printPropertyAcrossFilters(QPE_sersicIndices, "Sérsic Index")
    printPropertyAcrossFilters(QPE_r50s, "Sérsic half-light radius (kpc)")
    printPropertyAcrossFilters(QPE_magnitudes, "Magnitude")
    printPropertyAcrossFilters(QPE_unreddenedMagnitudes, "Dereddened magnitude")

    #From now on, keep only the r-band properties:
    QPE_placeholder_properties = {"n_sersic":[], "r_50":[]}
    for i in range(len(objects)):
        QPE_placeholder_properties["n_sersic"].append(QPE_sersicIndices[i]["r"])
        QPE_placeholder_properties["r_50"].append(QPE_r50s[i]["r"])
    QPE_sersicIndices = QPE_placeholder_properties["n_sersic"]
    QPE_r50s = QPE_placeholder_properties["r_50"]
    
    QPE_stellar_masses, TDE_stellar_masses = chooseStellarMassMethod()
    
    #Calculate stellar surface densities:
    QPE_stellar_surface_densities = []
    for i in range(len(objects)):
        QPE_stellar_surface_densities.append(stellarMassDensity(QPE_stellar_masses[i], QPE_r50s[i], returnLog=True))
    
    TDE_stellar_surface_densities = []
    TDE_r50s = add_0_uncertainties(TDE_r50s)
    for i in range(len(TDE_stellar_masses)):
        TDE_stellar_surface_densities.append(stellarMassDensity(TDE_stellar_masses[i], TDE_r50s[i], returnLog=True))
    
    #Transform lists into arrays
    QPE_mBH = np.array(QPE_mBH)
    QPE_sersicIndices = np.array(QPE_sersicIndices)
    QPE_stellar_surface_densities = np.array(QPE_stellar_surface_densities)
    TDE_mBH = np.array(TDE_mBH)
    TDE_sersicIndices = np.array(TDE_sersicIndices)
    TDE_stellar_surface_densities = np.array(TDE_stellar_surface_densities)

    for i in range(len(QPE_stellar_masses)):
        QPE_stellar_masses[i] = tuple(QPE_stellar_masses[i])
    for i in range(len(TDE_stellar_masses)):
        TDE_stellar_masses[i] = tuple(TDE_stellar_masses[i])
    QPE_stellar_masses = np.array(QPE_stellar_masses)
    TDE_stellar_masses = np.array(TDE_stellar_masses)


    plot_surfaceStellarMassDensity_stellarMass(QPE_stellar_masses, QPE_stellar_surface_densities, TDE_stellar_masses, TDE_stellar_surface_densities)
    makeFinalCornerPlot(QPE_mBH[:,0], QPE_sersicIndices[:,0], QPE_stellar_surface_densities[:,0], TDE_mBH, TDE_sersicIndices, TDE_stellar_surface_densities[:,0])
    sys.exit()

    # QPE hosts properties
    print_table(np.array([objects_names, np.around(QPE_sersicIndices[:,0], 4), np.around(QPE_r50s[:,0], 4), np.around(np.log10(QPE_stellar_masses[:,0]), 4), np.around(np.log10(QPE_mBH[:,0]), 4)]).T,
                    header=["Name", "Sérsic index", "r50 (kpc)", "log Stellar mass", "log M_BH"],
                    title="QPE hosts properties",
                    space_between_columns=4,
                    space_between_rows=0,
                    borders=2)
    
    # TDE hosts properties
    print_table(np.array([TDE_names, np.around(TDE_sersicIndices, 4), np.around(TDE_r50s[:,0], 4), np.around(np.log10(TDE_stellar_masses[:,0]), 4), np.around(np.log10(TDE_mBH), 4)]).T,
                    header=["Name", "Sérsic index", "r50 (kpc)", "log Stellar mass", "log M_BH"],
                    title="TDE hosts properties",
                    space_between_columns=4,
                    space_between_rows=0,
                    borders=2)

    