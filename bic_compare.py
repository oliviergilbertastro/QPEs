import pickle
import numpy as np
from utils import print_table, myCornerPlot
import matplotlib.pyplot as plt
from download_data import *
import sys



def get_bic(picklename):
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        fitting_run_result = pickle.load(open("galight_fitruns/big_fits/"+picklename,'rb'))
    # Get the BIC
    bic = fitting_run_result.fitting_seq.bic
    return bic



if __name__ == "__main__":

    survey = "DESI"
    band = "r"
    bics = []
    for model in ["None","AGN","Bulge"]:
        bics.append([])
        for name in objects_names:
            picklename = f"{name}_{band}-band_{model}_{survey}.pkl"
            try:
                bic = get_bic(picklename)
            except:
                bic = "-"
            bics[-1].append(bic)
    data = [objects_names]
    data.extend(bics)
    print_table(np.array(data).T, header=["Names", "None ","AGN","Bulge"], title="QPE BICs", borders=2)

    survey = "DESI_PSF"
    band = "r"
    bics = []
    for model in ["None","AGN","Bulge"]:
        bics.append([])
        for name in TDE_names:
            picklename = f"{name}_{band}-band_{model}_{survey}.pkl"
            try:
                bic = get_bic(picklename)
            except:
                bic = "-"
            bics[-1].append(bic)
    data = [TDE_names]
    data.extend(bics)
    print_table(np.array(data).T, header=["Names", "None ","AGN","Bulge"], title="TDE BICs", borders=2)