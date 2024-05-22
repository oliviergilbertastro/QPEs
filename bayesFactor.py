#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
import numpy as np
from galight_modif.tools.plot_tools import total_compare
from download_data import objects
import corner
import matplotlib.pyplot as plt

def getPosteriors(ra_dec, models=["None", "AGN", "Bulge", "Bulge+AGN"], band="i", verbose=True):
        
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
    for i in range(len(fitting_run_results)):
        if fitting_run_results[i].fitting_kwargs_list[-1][0] == 'MCMC':
            samples_mcmc = fitting_run_results[i].samples_mcmc
            n, num_param = np.shape(samples_mcmc)
            if verbose:
                print(n, num_param)
                print(fitting_run_results[i].param_mcmc)
            plot = corner.corner(samples_mcmc, labels=fitting_run_results[i].param_mcmc, show_titles=True)
            plt.show()
        

    return samples_mcmc


from download_data import objects
if input("Compare specific object? [y/n]") == "y":
    objID = int(input("Enter the object ID you want to load [0-7]:\n"))
    band = input("Enter the filter band you want to load [g,r,i,z]:\n")
    getPosteriors(objects[objID], band=band, models=['None', 'AGN'], verbose=True)
