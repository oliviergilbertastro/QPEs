#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
import numpy as np
from galight_modif.tools.plot_tools import total_compare
from download_data import objects

def compareModels(ra_dec, models=["None", "AGN", "Bulge", "Bulge+AGN"], band="i", verbose=True):
        
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
    print(f"{mid:.2f}_"+r"{-"+f"{minus:.2f}"+r"}^{+"+f"{plus:.2f}"+r"}")
    return types[bics.index(np.min(bics))] + " " + "?"*(len(models)-len(bics)) + f" n = {mid} +{plus}/{minus}"


from download_data import objects
if input("Compare specific object? [y/n]") == "y":
    objID = int(input("Enter the object ID you want to load [0-7]:\n"))
    band = input("Enter the filter band you want to load [g,r,i,z]:\n")
    compareModels(objects[objID], band=band, verbose=True)

if input("See best model for all objects? [y/n]") == "y":
    bands_list = ["i", "i", "i", "i", "z", "i", "i", "r"]
    print("Best models:")
    print("-------------------------------------------------")
    for i in range(8):
        print(f"Object {i}: {compareModels(objects[i], band=bands_list[i], models=['None', 'AGN'], verbose=False)}")
    print("-------------------------------------------------")

