import pickle
import numpy as np
from utils import print_table, myCornerPlot
import matplotlib.pyplot as plt
from download_data import *
import sys



def get_n_and_r50(objID, model="None", band="i", survey="DESI", redshift=0, qpe_oder_tde="QPE"):
    if qpe_oder_tde == "QPE":
        picklename = f"{objects_names[objID]}_{band}-band_{model}_{survey}.pkl"
    else:
        picklename = f"{TDE_names[objID]}_{band}-band_{model}_{survey}_PSF.pkl"
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        fitting_run_result = pickle.load(open("galight_fitruns/big_fits/"+picklename,'rb'))

    # Get the BIC
    bic = fitting_run_result.fitting_seq.bic

    return bic