import pickle
import numpy as np
from utils import print_table, myCornerPlot
import matplotlib.pyplot as plt
from download_data import *
import sys



def get_agn_flux(picklename):
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        fitting_run_result = pickle.load(open("galight_fitruns/big_fits/"+picklename,'rb'))
    # Get the BIC
    flux = fitting_run_result.mcmc_flux_list
    plt.plot(flux)
    plt.show()
    return flux



if __name__ == "__main__":

    survey = "DESI"
    band = "r"

    for model in ["AGN"]:
        bics.append([])
        aics.append([])
        for name in objects_names:
            picklename = f"{name}_{band}-band_{model}_{survey}.pkl"
            try:
                bic, aic = get_bic(picklename)
            except:
                bic, aic = "-", "-"
            bics[-1].append(bic)
            aics[-1].append(aic)
    data = [objects_names]
    data.extend(bics)
    print_table(np.array(data).T, header=["Names", "None ","AGN","Bulge"], title="QPE BICs", borders=2)
    data = [objects_names]
    data.extend(aics)
    print_table(np.array(data).T, header=["Names", "None ","AGN","Bulge"], title="QPE AICs", borders=2)
    #print("BIC:",make_list_of_best_models(bics))
    #print("AIC:",make_list_of_best_models(aics))
