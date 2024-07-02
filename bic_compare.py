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
    aic = fitting_run_result.fitting_seq.aic
    return bic, aic

def get_n(picklename):
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        fitting_run_result = pickle.load(open("galight_fitruns/big_fits/"+picklename,'rb'))
    # Get the sersic index
    try:
        n = fitting_run_result.final_result_galaxy[0]["n_sersic"], fitting_run_result.final_result_galaxy[1]["n_sersic"]
    except:
        n = fitting_run_result.final_result_galaxy[0]["n_sersic"]
    return n

def make_list_of_best_models(ic):
    """ic is an array of information criterions"""
    ic = np.array(ic).T
    models = []
    choices = ["None","AGN","Bulge","Bulge_fixed"]
    for row in ic:
        models.append(choices[list(row).index(np.min(row))])
    return models


def printSersicAcrossModels(properties, name_of_property="Name of property", round_to_n_decimals=2, qpe_oder_tde="QPE"):
    names = objects_names if qpe_oder_tde=="QPE" else TDE_names
    properties["None"] = np.around(properties["None"], round_to_n_decimals)
    properties["AGN"] = np.around(properties["AGN"], round_to_n_decimals)
    properties["Bulge (bulge)"] = np.around(np.array(properties["Bulge"])[:,0], round_to_n_decimals)
    properties["Bulge (disk)"] = np.around(np.array(properties["Bulge"])[:,1], round_to_n_decimals)

    print_table(np.array([names, properties["None"], properties["AGN"], properties["Bulge (bulge)"], properties["Bulge (disk)"]]).T,
                header=["Name","None ","AGN","Bulge (bulge)","Bulge (disk)"],
                title=name_of_property,
                borders=2,
                )
    return


if __name__ == "__main__":

    survey = "DESI_PSF"
    band = "r"
    bics = []
    aics = []
    n_sersics = {}
    for model in ["None","AGN","Bulge"]:
        bics.append([])
        aics.append([])
        n_sersics[model] = []
        for name in objects_names:
            picklename = f"{name}_{band}-band_{model}_{survey}.pkl"
            try:
                bic, aic = get_bic(picklename)
                n = get_n(picklename)
            except:
                bic, aic = "-", "-"
                n = "-"
            bics[-1].append(bic)
            aics[-1].append(aic)
            n_sersics[model].append(n)
    data = [objects_names]
    data.extend(bics)
    print_table(np.array(data).T, header=["Names", "None ","AGN","Bulge"], title="QPE BICs", borders=2)
    data = [objects_names]
    data.extend(aics)
    print_table(np.array(data).T, header=["Names", "None ","AGN","Bulge"], title="QPE AICs", borders=2)
    
    printSersicAcrossModels(n_sersics, "QPE sersic indices")

    survey = "DESI_PSF"
    band = "r"
    bics = []
    aics = []
    n_sersics = {}
    for model in ["None","AGN","Bulge","Bulge_fixed"]:
        bics.append([])
        aics.append([])
        n_sersics[model] = []
        for name in TDE_names:
            picklename = f"{name}_{band}-band_{model}_{survey}.pkl"
            try:
                bic, aic = get_bic(picklename)
                n = get_n(picklename)
            except:
                bic, aic = "-", "-"
                n = "-"
            bics[-1].append(bic)
            aics[-1].append(aic)
            n_sersics[model].append(n)
    data = [TDE_names]
    data.extend(bics)
    print_table(np.array(data).T, header=["Names", "None ", "AGN", "Bulge", "Bulge_fixed"], title="TDE BICs", borders=2)
    data = [TDE_names]
    data.extend(aics)
    print_table(np.array(data).T, header=["Names", "None ", "AGN", "Bulge", "Bulge_fixed"], title="TDE AICs", borders=2)
    print("BIC:",make_list_of_best_models(bics))
    print("AIC:",make_list_of_best_models(aics))
    printSersicAcrossModels(n_sersics, "TDE sersic indices", qpe_oder_tde="TDE")
    plt.plot(n_sersics["None"], label="Sérsic")
    plt.plot(n_sersics["AGN"], label="Sérsic+AGN")
    plt.ylabel("Sérsic index", fontsize=16)
    plt.legend()
    print(np.mean(n_sersics["AGN"])-np.mean(n_sersics["None"]))
    plt.show()