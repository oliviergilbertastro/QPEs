"""
Program to compare QPE hosts and TDE hosts "apples to apples" with the LEGACY DESI DECam survey
"""
import pickle
import numpy as np
from ned_wright_cosmology import calculate_cosmo
from utils import print_table, myCornerPlot, toLog, myFinalPlot, myCombinedFinalPlot, recombine_arrays, add_0_uncertainties, makeLatexTable
import matplotlib.pyplot as plt
from paper_data import *
from download_data import *
import sys
import pandas as pd


QPE_and_TDEs = [4,5,8] # indices of QPE+TDE hosts which are currently only in the QPE arrays


if __name__ == "__main__":

    # load QPE and TDE data
    qpe_load_0 = np.loadtxt("QPE_allRelevantData_0.txt")
    qpe_load_1 = np.loadtxt("QPE_allRelevantData_1.txt")
    qpe_load_2 = np.loadtxt("QPE_allRelevantData_2.txt")
    QPE_fullData = recombine_arrays(qpe_load_0,qpe_load_1,qpe_load_2)

    tde_load_0 = np.loadtxt("TDE_allRelevantData_0.txt")
    tde_load_1 = np.loadtxt("TDE_allRelevantData_1.txt")
    tde_load_2 = np.loadtxt("TDE_allRelevantData_2.txt")
    TDE_fullData = recombine_arrays(tde_load_0,tde_load_1,tde_load_2)
    TDE_fullData = TDE_fullData[:10,:,:]

    # load reference catalog
    refCat = np.loadtxt("referenceCatalog_modif2.txt")
    fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
    refCat = pd.read_csv("referenceCatalog_modif2.txt", delimiter=" ", header=None, names=fieldnames)

    QPE_r50s, TDE_r50s = QPE_fullData[:,5,:], TDE_fullData[:,5,:]
    QPE_sersicIndices, TDE_sersicIndices = QPE_fullData[:,1,:], TDE_fullData[:,1,:]
    QPE_bulgeRatios, TDE_bulgeRatios = QPE_fullData[:,0,:], TDE_fullData[:,0,:]
    QPE_SMSDs, TDE_SMSDs = QPE_fullData[:,2,:], TDE_fullData[:,2,:]
    QPE_stellar_masses, TDE_stellar_masses = QPE_fullData[:,3,:], TDE_fullData[:,3,:]
    QPE_mBH, TDE_mBH = QPE_fullData[:,4,:], TDE_fullData[:,4,:]

    makeLatexTable(np.concatenate((objects_names, TDE_names)),
                   np.concatenate((QPE_redshifts,TDE_redshifts)),
                   np.concatenate((QPE_r50s,TDE_r50s)),
                   np.concatenate((QPE_sersicIndices,TDE_sersicIndices)),
                   np.concatenate((QPE_bulgeRatios,TDE_bulgeRatios)),
                   np.concatenate((QPE_SMSDs,TDE_SMSDs)),
                   np.concatenate((QPE_stellar_masses,TDE_stellar_masses)),
                   references="abccefddghijklmnnno"
                   )

    QPE_data = np.array([QPE_sersicIndices, QPE_bulgeRatios, QPE_SMSDs, QPE_stellar_masses, QPE_mBH])
    TDE_data = np.array([np.concatenate((TDE_sersicIndices, QPE_sersicIndices[QPE_and_TDEs])), np.concatenate((TDE_bulgeRatios, QPE_bulgeRatios[QPE_and_TDEs])), np.concatenate((TDE_SMSDs, QPE_SMSDs[QPE_and_TDEs])), np.concatenate((TDE_stellar_masses, QPE_stellar_masses[QPE_and_TDEs])), np.concatenate((TDE_mBH, QPE_mBH[QPE_and_TDEs]))])
    myCombinedFinalPlot([QPE_data, TDE_data], referenceCatalogData=refCat, columns_compare=((60,12,68),63,67), save_plot="combined_final", fontsize=16, markersize=9, levels=[0.3,0.5,0.89,0.9,1])

    # Make big plot
    QPE_data  = np.array([QPE_mBH[:,0], QPE_stellar_masses[:,0], QPE_bulgeRatios[:,0], QPE_r50s[:,0], QPE_sersicIndices[:,0], QPE_SMSDs[:,0]])
    TDE_data = np.array([TDE_mBH[:,0], TDE_stellar_masses[:,0], TDE_bulgeRatios[:,0], TDE_r50s[:,0], TDE_sersicIndices[:,0], TDE_SMSDs[:,0]])
    double_hosts_data = QPE_data[:,QPE_and_TDEs]
    TDE_data = np.vstack((TDE_data.T, double_hosts_data.T)).T
    myCornerPlot(
        [QPE_data,TDE_data,double_hosts_data],
        labels=["$\log(M_\mathrm{BH})$", "$\log(M_\star)$", "$(B/T)_g$", "$r_{50}$", "$n_\mathrm{Sérsic}$", "$\log(\Sigma_{M_\star})$"],
        units=["$M_\odot$", "$M_\odot$", " ", "$\mathrm{kpc}$", " ", "$M_\odot/\mathrm{kpc}^2$"],
        smoothness=6,
        refCat=refCat,
        columns_compare=[67,63,12,59,60,68],
        save_plot="corner_plot"
        )
