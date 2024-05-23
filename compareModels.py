#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
import numpy as np
from galight_modif.tools.plot_tools import total_compare
from download_data import objects

def compareModels(ra_dec, models=["None", "AGN", "Bulge", "Bulge+AGN"], band="i", stellar_mass=5E9, verbose=True):
        
    if band in ["g", "G"]:
        band = "g"
    elif band in ["r", "R"]:
        band = "r"
    elif band in ["i", "I"]:
        band = "i"
    elif band in ["z", "Z"]:
        band = "z"
    else:
        raise ValueError(f"band {band} is not a supported filter band")
    
    fitting_run_results = []
    types = []
    for type in models:
        picklename = f'ra{str(ra_dec[0])}_dec{str(ra_dec[1])}_{type}_{band}.pkl'
        try:
            fitting_run_results.append(pickle.load(open("galight_fitruns/"+picklename,'rb')))  #fitting_run_result is actually the fit_run in galightFitting.py.
            types.append(type)
        except:
            if verbose:
                print(f"Model {type} has not been computed yet. Skipping it.")

    bics = []
    if verbose:
        print("BIC:")
        print("-------------------------------------")
    for i in range(len(fitting_run_results)):
        bics.append(fitting_run_results[i].fitting_seq.bic)
        if verbose:
            print(types[i], ":", fitting_run_results[i].fitting_seq.bic)
    best_model_index = bics.index(np.min(bics))
    if verbose:
        print("-------------------------------------")
        print("Best model:", types[best_model_index])
    #Calculate the Sersic index + uncertainties
    chain = fitting_run_results[best_model_index].samples_mcmc
    lo, mid, hi = np.percentile(chain[:, 1],16), np.percentile(chain[:, 1],50), np.percentile(chain[:, 1],84)
    plus, minus = (hi-mid), (mid-lo)
    sersic_index_str = (f"{mid:.2f}_"+r"{-"+f"{minus:.2f}"+r"}^{+"+f"{plus:.2f}"+r"}")

    #Calculate the Sersic half-light radius + uncertainties:
    lo, mid, hi = np.percentile(chain[:, 0],16), np.percentile(chain[:, 0],50), np.percentile(chain[:, 0],84)
    plus, minus = (hi-mid), (mid-lo)

    stellar_density_str = "-"
    if stellar_mass != None:
        smd = np.array(stellarMassDensity(stellar_mass, [mid, minus, plus]))/1E9
        stellar_density_str = (f"{smd[0]:.2f}_"+r"{-"+f"{smd[1]:.2f}"+r"}^{+"+f"{smd[2]:.2f}"+r"}")

    return types[bics.index(np.min(bics))] + " " + "?"*(len(models)-len(bics)) + f' n = {sersic_index_str}      \Sigma_star = {stellar_density_str}'

def stellarMassDensity(M_star, r50):
    '''
    Calculate the stellar surface mass density \Sigma_{M_\ast} from the total stellar mass and the half-light radius

    M_star: stellar mass tuple/list/array such as (mass, errlo, errhi)
    r50: half-light radius tuple/list/array such as (radius, errlo, errhi)

    returns [value, errlo, errhi]
    '''
    try:
        res = M_star[0]/r50[0]**2
        errlo = res*np.sqrt((M_star[1]/M_star[0])**2+(2*r50[2]/r50[0])**2)
        errhi = res*np.sqrt((M_star[2]/M_star[0])**2+(2*r50[1]/r50[0])**2)
        return [res, errlo, errhi]
    except:
        #In the exception where the user did not input uncertainties, the value will still be calculated.
        return M_star/r50**2
    

from download_data import objects
if input("Compare specific object? [y/n]") == "y":
    objID = int(input("Enter the object ID you want to load [0-7]:\n"))
    band = input("Enter the filter band you want to load [g,r,i,z]:\n")
    compareModels(objects[objID], band=band, verbose=True)

if input("See best model for all objects? [y/n]") == "y":
    bands_list = ["i", "i", "i", "i", "z", "i", "i", "r"]
    stellar_masses = [
                (None),           # GSN 069
                (None),           # RX J1301.9+2747
                (3.8E9, 1.9E9, 0.4E9),           # eRO-QPE1
                (1.01E9, 0.5E9, 0.01E9),           # eRO-QPE2
                (None),           # AT 2019vcb
                (None),           # 2MASX J0249
                (2.56E9, 1.40E9, 0.24E9),           # eRO-QPE3
                (1.6E10, 0.6E10, 0.7E10),           # eRO-QPE4
                ]
    print("Best models:")
    print("-------------------------------------------------")
    for i in range(8):
        print(f"Object {i}: {compareModels(objects[i], band=bands_list[i], stellar_mass=stellar_masses[i], models=['None', 'AGN'], verbose=False)}")
    print("-------------------------------------------------")

