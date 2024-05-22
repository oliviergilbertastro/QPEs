#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
import numpy as np
from galight_modif.tools.plot_tools import total_compare
from download_data import objects
import corner
import matplotlib.pyplot as plt
from copy import copy

def getPosteriors(ra_dec, models=["None", "AGN", "Bulge", "Bulge+AGN"], band="i", verbose=True, if_plot=False):
        
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
    
    #Fetch the MCMC posteriors
    posteriors = []
    for i in range(len(fitting_run_results)):
        if fitting_run_results[i].fitting_kwargs_list[-1][0] == 'MCMC':
            samples_mcmc = fitting_run_results[i].samples_mcmc
            n, num_param = np.shape(samples_mcmc)
            amp_params = int(fitting_run_results[i].fitting_seq.param_class.num_param_linear())
            total_params = num_param + amp_params
            if verbose:
                print(n, num_param)
                print("Total amplitude params:", fitting_run_results[i].fitting_seq.param_class.num_param_linear())
                print("LensLight amp params:", fitting_run_results[i].fitting_seq.param_class.lensLightParams.num_param_linear())
                print("PointSource amp params:", fitting_run_results[i].fitting_seq.param_class.pointSourceParams.num_param_linear())
                print("PointSource model amp kwargs list:", fitting_run_results[i].fitting_seq.param_class.pointSourceParams.kwargs_fixed)
                print("LensLight model amp kwargs list:", fitting_run_results[i].fitting_seq.param_class.lensLightParams.kwargs_fixed)
                print("LensLight model amp number of params:", fitting_run_results[i].fitting_seq.param_class.lensLightParams._lightModel.num_param_linear(kwargs_list=fitting_run_results[i].fitting_seq.param_class.lensLightParams.kwargs_fixed, list_return=True))
                kwargs_result = fitting_run_results[i].fitting_seq.best_fit()
                ps_result = kwargs_result['kwargs_ps']
                source_result = kwargs_result['kwargs_lens_light']
                print(ps_result)
                print(source_result)
            #plot = corner.corner(samples_mcmc, labels=fitting_run_results[i].param_mcmc, show_titles=True)
            #plt.show()
            #plt.plot(samples_mcmc[:,1])
            #plt.ylabel("Sérsic index")
            #plt.show()
            histograms_data = np.empty((n, total_params))
            histograms_data[:,:num_param] = np.copy(samples_mcmc)
            amplitudes = np.array(fitting_run_results[i].amplitudes_list)[:,0]
            labels = copy(fitting_run_results[i].param_mcmc)
            for k in range(amp_params):
                histograms_data[:,num_param+k] = amplitudes[:,k]
                if k < fitting_run_results[i].fitting_seq.param_class.lensLightParams.num_param_linear():
                    labels.append("host amplitude")
                else:
                    labels.append("agn amplitude")
            if False:
                fig, axes = plt.subplots(nrows=1, ncols=total_params)
                for k in range(total_params):
                    axes[k].hist(histograms_data[:,k], bins=20)
                    axes[k].set_xlabel(labels[k], fontsize=15)
                plt.show()
            if if_plot:
                plot = corner.corner(histograms_data, labels=labels, show_titles=True)
                plt.show()
        posteriors.append(histograms_data)
    return posteriors


from download_data import objects
if input("Compare specific object? [y/n]") == "y":
    objID = int(input("Enter the object ID you want to load [0-7]:\n"))
    band = input("Enter the filter band you want to load [g,r,i,z]:\n")
    getPosteriors(objects[objID], band=band, models=['None', 'AGN'], verbose=False)
