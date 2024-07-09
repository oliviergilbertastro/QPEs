"""
Save magnitudes to numpy array text file so they are faster to load and don't need to have original files after running this once.
"""
import pickle
import numpy as np
from ned_wright_cosmology import calculate_cosmo
from utils import print_table, myCornerPlot, toLog, myFinalPlot
import matplotlib.pyplot as plt
from paper_data import *
from download_data import *
import sys
from prospector.prospector_myFits import QPE_stellar_masses_desiProspector, TDE_stellar_masses_desiProspector
#from bulgeRatio import QPE_bulgeRatios, TDE_bulgeRatios

QPE_bulgeRatios, TDE_bulgeRatios = np.loadtxt("QPE_bulgeRatios.txt"), np.loadtxt("TDE_bulgeRatios.txt")

QPE_and_TDEs = [4,5,8]

def get_n_and_r50(objID, model="None", band="i", survey="DESI", redshift=0, qpe_oder_tde="QPE"):
    if qpe_oder_tde == "QPE":
        picklename = f"{objects_names[objID]}_{band}-band_{model}_{survey}.pkl"
    else:
        picklename = f"{TDE_names[objID]}_{band}-band_{model}_{survey}.pkl"
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/big_fits/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))
    #print(picklename)
    #Calculate the Sersic index + uncertainties
    chain = fitting_run_result.samples_mcmc
    params = fitting_run_result.param_mcmc
    if "R_sersic_lens_light1" in params and "n_sersic_lens_light1" in params:
        lo, mid, hi = np.percentile(chain[:, 7],16), np.percentile(chain[:, 7],50), np.percentile(chain[:, 7],84)
        plus, minus = (hi-mid), (mid-lo)
        sersic_index_data = [mid, minus, plus]

        #Calculate the Sersic half-light radius + uncertainties:
        lo, mid, hi = np.percentile(chain[:, 6],16), np.percentile(chain[:, 6],50), np.percentile(chain[:, 6],84)
        plus, minus = (hi-mid), (mid-lo)
        r50_data = [mid, minus, plus]
    elif "R_sersic_lens_light0" in params and "n_sersic_lens_light0" in params:
        lo, mid, hi = np.percentile(chain[:, 1],16), np.percentile(chain[:, 1],50), np.percentile(chain[:, 1],84)
        plus, minus = (hi-mid), (mid-lo)
        sersic_index_data = [mid, minus, plus]

        #Calculate the Sersic half-light radius + uncertainties:
        lo, mid, hi = np.percentile(chain[:, 0],16), np.percentile(chain[:, 0],50), np.percentile(chain[:, 0],84)
        plus, minus = (hi-mid), (mid-lo)
        r50_data = [mid, minus, plus]
    else:
        #print(qpe_oder_tde)
        sersic_index_data = [fitting_run_result.final_result_galaxy[0]["n_sersic"], 0, 0]
        r50_data = [fitting_run_result.final_result_galaxy[0]["R_sersic"], 0, 0]
    #Convert r50 from arcsec to kpc
    cosmology_params = calculate_cosmo(redshift, H0=67, Omega_m=0.31, Omega_v=0.69)
    kpc_per_arcsec = cosmology_params["kpc_DA"]
    r50_data = np.array(r50_data)*kpc_per_arcsec

    magnitude = fitting_run_result.final_result_galaxy[0]['magnitude']
    # Add the magnitude of the disk if there is a bulge+disk decomposition
    try:
        magnitude = -2.5*np.log10(10**(-0.4*magnitude)+10**(-0.4*fitting_run_result.final_result_galaxy[1]['magnitude']))
    except:
        pass

    return sersic_index_data, r50_data, magnitude



def printPropertyAcrossFilters(list_of_dicts, name_of_property="Name of property", round_to_n_decimals=2, qpe_oder_tde="QPE"):
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
    names = objects_names if qpe_oder_tde=="QPE" else TDE_names
    print_table(np.array([names, properties["g"], properties["r"], properties["i"], properties["z"]]).T,
                header=["Name", "g", "r", "i", "z"],
                title=name_of_property,
                borders=2,
                )
    return










# We don't use the magnitudes calculated from the single Sérsic profile fit, so we recalculate them using the bulge models:
QPE_magnitudes = []
QPE_unreddenedMagnitudes = []
for i in range(len(objects)):
    QPE_magnitudes.append({})
    QPE_unreddenedMagnitudes.append({})
    for band in "griz":
        try:
            _, _, mag = get_n_and_r50(i, "Bulge", redshift=QPE_redshifts[i], band=band, survey="DESI_PSF_FINAL2")
            QPE_magnitudes[-1][band] = mag
            QPE_unreddenedMagnitudes[-1][band] = mag - QPE_extinction[objects_names[i]][band]
        except:
            pass




# We don't use the magnitudes calculated from the single Sérsic profile fit, so we recalculate them using the bulge models:
TDE_magnitudes = []
TDE_unreddenedMagnitudes = []
for i in range(len(TDE_coords)):
    TDE_magnitudes.append({})
    TDE_unreddenedMagnitudes.append({})
    for band in "griz":
        try:
            _, _, mag = get_n_and_r50(i, "Bulge", redshift=TDE_redshifts[i], band=band, survey="DESI_PSF_FINAL2", qpe_oder_tde="TDE")
            TDE_magnitudes[-1][band] = mag
            TDE_unreddenedMagnitudes[-1][band] = mag - TDE_extinction[TDE_names[i]][band]
        except:
            pass


QPE_unreddenedMagnitudes = np.array(QPE_unreddenedMagnitudes)
TDE_unreddenedMagnitudes = np.array(TDE_unreddenedMagnitudes)


if __name__ == "__main__":

    printPropertyAcrossFilters(QPE_magnitudes, "Magnitude")
    printPropertyAcrossFilters(QPE_unreddenedMagnitudes, "Dereddened magnitude")

    printPropertyAcrossFilters(TDE_magnitudes, "Magnitude", qpe_oder_tde="TDE")
    printPropertyAcrossFilters(TDE_unreddenedMagnitudes, "Dereddened magnitude", qpe_oder_tde="TDE")

    if input("Overwrite the pickle files used in prospector_myFits? [y/n]") == "y":
        with open('QPE_magnitudes.pkl', 'wb') as f:
            pickle.dump(QPE_unreddenedMagnitudes, f)
        with open('TDE_magnitudes.pkl', 'wb') as f:
            pickle.dump(TDE_unreddenedMagnitudes, f)