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

def make_list_of_best_models(ic, choices=["None","AGN","Bulge","Bulge_fixed","Bulge+AGN"]):
    """ic is an array of information criterions"""
    ic = np.array(ic).T
    models = []
    for row in ic:
        models.append(choices[list(row).index(np.min(row))])
    return models


def printSersicAcrossModels(properties, name_of_property="Name of property", round_to_n_decimals=2, qpe_oder_tde="QPE"):
    names = objects_names if qpe_oder_tde=="QPE" else TDE_names
    properties["None"] = np.around(properties["None"], round_to_n_decimals)
    properties["AGN"] = np.around(properties["AGN"], round_to_n_decimals)
    properties["Bulge (bulge)"] = np.around(np.array(properties["Bulge"])[:,0], round_to_n_decimals)
    properties["Bulge (disk)"] = np.around(np.array(properties["Bulge"])[:,1], round_to_n_decimals)
    properties["Bulge+AGN (bulge)"] = np.around(np.array(properties["Bulge+AGN"])[:,0], round_to_n_decimals)
    properties["Bulge+AGN (disk)"] = np.around(np.array(properties["Bulge+AGN"])[:,1], round_to_n_decimals)

    print_table(np.array([names, properties["None"], properties["AGN"], properties["Bulge (bulge)"], properties["Bulge (disk)"], properties["Bulge+AGN (bulge)"], properties["Bulge+AGN (disk)"]]).T,
                header=["Name","None ","AGN","Bulge (bulge)","Bulge (disk)","Bulge+AGN (bulge)","Bulge+AGN (disk)"],
                title=name_of_property,
                borders=2,
                )
    return




survey = "DESI_PSF"
band = "r"
bics = []
aics = []
n_sersics = {}
for model in ["None","AGN","Bulge","Bulge+AGN"]:
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
if __name__ == "__main__":
    print_table(np.array(data).T, header=["Names", "None ","AGN","Bulge","Bulge+AGN"], title="QPE BICs", borders=2)
data = [objects_names]
data.extend(aics)

QPE_best_models = make_list_of_best_models(bics, choices=["None","AGN","Bulge","Bulge+AGN"])
if __name__ == "__main__":
    print_table(np.array(data).T, header=["Names", "None ","AGN","Bulge","Bulge+AGN"], title="QPE AICs", borders=2)
    print("BIC:",make_list_of_best_models(bics, choices=["None","AGN","Bulge","Bulge+AGN"]))
    print("AIC:",make_list_of_best_models(aics, choices=["None","AGN","Bulge","Bulge+AGN"]))
    printSersicAcrossModels(n_sersics, "QPE sersic indices")

survey = "DESI_PSF"
band = "r"
bics = []
aics = []
n_sersics = {}
for model in ["None","AGN","Bulge","Bulge_fixed","Bulge+AGN"]:
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
if __name__ == "__main__":
    print_table(np.array(data).T, header=["Names", "None ", "AGN", "Bulge", "Bulge_fixed","Bulge+AGN"], title="TDE BICs", borders=2)
data = [TDE_names]
data.extend(aics)

TDE_best_models = make_list_of_best_models(bics)
if __name__ == "__main__":
    print_table(np.array(data).T, header=["Names", "None ", "AGN", "Bulge", "Bulge_fixed","Bulge+AGN"], title="TDE AICs", borders=2)
    print("BIC:",make_list_of_best_models(bics))
    print("AIC:",make_list_of_best_models(aics))
    printSersicAcrossModels(n_sersics, "TDE sersic indices", qpe_oder_tde="TDE")
    plt.plot(n_sersics["None"], 'o', label="Sérsic")
    plt.plot(n_sersics["AGN"], 'o', label="Sérsic+AGN")
    plt.ylabel("Sérsic index", fontsize=16)
    plt.legend()
    print(np.mean(n_sersics["AGN"])-np.mean(n_sersics["None"]))
    plt.show()