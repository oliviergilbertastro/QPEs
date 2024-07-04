import pickle
import numpy as np
from utils import print_table, myCornerPlot
import matplotlib.pyplot as plt
from download_data import *
#from bic_compare import QPE_best_models, TDE_best_models
import sys


def get_fluxes(picklename):
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/big_fits/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
        #print("bigfit")
    except:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))
    # Get the sersic index
    try:
        fluxes = fitting_run_result.mcmc_flux_list
        labels = fitting_run_result.labels_flux
    except:
        pass
    return fluxes, labels




survey = "DESI_PSF"
band = "g"
QPE_bulgeRatios = []
for i in range(len(objects_names)):
    fluxes, labels = get_fluxes(f"{objects_names[i]}_{band}-band_{'Bulge'}_{survey}.pkl")
    bulge_flux = fluxes[:,labels.index("Galaxy_0 flux")] # Bulge flux
    bulge_flux = (np.median(bulge_flux), np.median(bulge_flux)-np.quantile(bulge_flux, 0.16), np.quantile(bulge_flux, 0.84)-np.median(bulge_flux))
    disk_flux = fluxes[:,labels.index("Galaxy_1 flux")] # Disk flux
    disk_flux = (np.median(disk_flux), np.median(disk_flux)-np.quantile(disk_flux, 0.16), np.quantile(disk_flux, 0.84)-np.median(disk_flux))
    total_flux = (bulge_flux[0]+disk_flux[0], np.sqrt(bulge_flux[1]**2+disk_flux[1]**2), np.sqrt(bulge_flux[2]**2+disk_flux[2]**2))
    bulge_ratio = bulge_flux[0]/total_flux[0]
    bulge_ratio = (bulge_ratio, bulge_ratio*np.sqrt((bulge_flux[1]/bulge_flux[0])**2+(total_flux[1]/total_flux[0])**2), bulge_ratio*np.sqrt((bulge_flux[2]/bulge_flux[0])**2+(total_flux[2]/total_flux[0])**2))
    if bulge_ratio[0] != bulge_ratio[0]:
        bulge_ratio = (0,0,0)
    QPE_bulgeRatios.append(bulge_ratio)
QPE_bulgeRatios = np.array(QPE_bulgeRatios)
TDE_bulgeRatios = []
for i in range(len(TDE_names)):
    fluxes, labels = get_fluxes(f"{TDE_names[i]}_{band}-band_{'Bulge'}_{survey}.pkl")
    bulge_flux = fluxes[:,labels.index("Galaxy_0 flux")] # Bulge flux
    bulge_flux = (np.median(bulge_flux), np.median(bulge_flux)-np.quantile(bulge_flux, 0.16), np.quantile(bulge_flux, 0.84)-np.median(bulge_flux))
    disk_flux = fluxes[:,labels.index("Galaxy_1 flux")] # Disk flux
    disk_flux = (np.median(disk_flux), np.median(disk_flux)-np.quantile(disk_flux, 0.16), np.quantile(disk_flux, 0.84)-np.median(disk_flux))
    total_flux = (bulge_flux[0]+disk_flux[0], np.sqrt(bulge_flux[1]**2+disk_flux[1]**2), np.sqrt(bulge_flux[2]**2+disk_flux[2]**2))
    bulge_ratio = bulge_flux[0]/total_flux[0]
    bulge_ratio = (bulge_ratio, bulge_ratio*np.sqrt((bulge_flux[1]/bulge_flux[0])**2+(total_flux[1]/total_flux[0])**2), bulge_ratio*np.sqrt((bulge_flux[2]/bulge_flux[0])**2+(total_flux[2]/total_flux[0])**2))
    TDE_bulgeRatios.append(bulge_ratio)
TDE_bulgeRatios = np.array(TDE_bulgeRatios)

if __name__ == "__main__":
    plt.errorbar(range(9), QPE_bulgeRatios[:,0], [QPE_bulgeRatios[:,1], QPE_bulgeRatios[:,2]], fmt="o", label="QPEs")
    plt.errorbar(range(10), TDE_bulgeRatios[:,0], [TDE_bulgeRatios[:,1], TDE_bulgeRatios[:,2]], fmt="o", label="TDEs")
    plt.legend()
    plt.ylabel(r"$(B/T)_g$", fontsize=16)
    plt.xlabel("Index", fontsize=16)
    plt.show()