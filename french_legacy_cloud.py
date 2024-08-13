"""
Program to compare QPE hosts and TDE hosts "apples to apples" with the LEGACY DESI DECam survey
"""
import pickle
import numpy as np
from ned_wright_cosmology import calculate_cosmo
from utils import print_table, myCornerPlot, toLog, myFinalPlot, myCombinedFinalPlot, recombine_arrays, add_0_uncertainties, makeLatexTable, redshiftMass, cut_from_catalog
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

    french_tde_load_0 = np.loadtxt("french_TDE_allRelevantData_0.txt")
    french_tde_load_1 = np.loadtxt("french_TDE_allRelevantData_1.txt")
    french_tde_load_2 = np.loadtxt("french_TDE_allRelevantData_2.txt")
    french_TDE_fullData = recombine_arrays(french_tde_load_0,french_tde_load_1,french_tde_load_2)
    french_TDE_fullData = french_TDE_fullData[:,:,:]

    offset = input("Assuming ____ was wrong. Choices: ['Me','Mendel','no one (Enter)']")
    m_star_offsets = [0,0]
    referenceCata_extension = ""
    if offset.lower() == "me":
        QPE_fullData[:,3,0] = QPE_fullData[:,3,0]+0.2370666
        TDE_fullData[:,3,0] = TDE_fullData[:,3,0]+0.2370666
        french_TDE_fullData[:,3,0] = french_TDE_fullData[:,3,0]+0.2370666
        referenceCata_extension = "_prospector"
        m_star_offsets = [0.2370666, 0.2370666]
    elif offset.lower() == "mendel":
        referenceCata_extension = "_mendel"

    # load reference catalog
    refCat = np.loadtxt(f"referenceCatalog_final{referenceCata_extension}.txt")
    fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
    refCat = pd.read_csv(f"referenceCatalog_final{referenceCata_extension}.txt", delimiter=" ", header=None, names=fieldnames)

    QPE_r50s, TDE_r50s, french_TDE_r50s = QPE_fullData[:,5,:], TDE_fullData[:,5,:], french_TDE_fullData[:,5,:]
    QPE_sersicIndices, TDE_sersicIndices, french_TDE_sersicIndices = QPE_fullData[:,1,:], TDE_fullData[:,1,:], french_TDE_fullData[:,1,:]
    QPE_bulgeRatios, TDE_bulgeRatios, french_TDE_bulgeRatios = QPE_fullData[:,0,:], TDE_fullData[:,0,:], french_TDE_fullData[:,0,:]
    QPE_SMSDs, TDE_SMSDs, french_TDE_SMSDs = QPE_fullData[:,2,:], TDE_fullData[:,2,:], french_TDE_fullData[:,2,:]
    QPE_stellar_masses, TDE_stellar_masses, french_TDE_stellar_masses = QPE_fullData[:,3,:], TDE_fullData[:,3,:], french_TDE_fullData[:,3,:]
    QPE_mBH, TDE_mBH, french_TDE_mBH = QPE_fullData[:,4,:], TDE_fullData[:,4,:], french_TDE_fullData[:,4,:]

    # Concatenate data from Law-Smith and French
    TDE_r50s = np.concatenate((TDE_r50s,french_TDE_r50s))
    TDE_sersicIndices = np.concatenate((TDE_sersicIndices, french_TDE_sersicIndices))
    TDE_bulgeRatios = np.concatenate((TDE_bulgeRatios, french_TDE_bulgeRatios))
    TDE_SMSDs = np.concatenate((TDE_SMSDs, french_TDE_SMSDs))
    TDE_stellar_masses = np.concatenate((TDE_stellar_masses, french_TDE_stellar_masses))
    TDE_mBH = np.concatenate((TDE_mBH, french_TDE_mBH))
    TDE_names = np.concatenate((TDE_names, french_TDE_names))
    TDE_redshifts = np.concatenate((TDE_redshifts, french_TDE_redshifts))

    if input("Show only French TDEs? [y/n]") == "y":
        TDE_r50s = french_TDE_r50s
        TDE_sersicIndices = french_TDE_sersicIndices
        TDE_bulgeRatios = french_TDE_bulgeRatios
        TDE_SMSDs = french_TDE_SMSDs
        TDE_stellar_masses = french_TDE_stellar_masses
        TDE_mBH = french_TDE_mBH
        TDE_names = french_TDE_names
        TDE_redshifts = french_TDE_redshifts
        QPE_and_TDEs = []

    try:
        makeLatexTable(np.concatenate((objects_names, TDE_names)),
                    np.concatenate((QPE_redshifts,TDE_redshifts)),
                    np.concatenate((QPE_r50s,TDE_r50s)),
                    np.concatenate((QPE_sersicIndices,TDE_sersicIndices)),
                    np.concatenate((QPE_bulgeRatios,TDE_bulgeRatios)),
                    np.concatenate((QPE_SMSDs,TDE_SMSDs)),
                    np.concatenate((QPE_stellar_masses,TDE_stellar_masses)),
                    references="",#"abccefddghijklmnnno",
                    filename="french_latexTable.txt"
                    )
    except:
        print("\x1b[31mCouldn't make a latex table :(\x1b[0m")

    QPE_redshifts = add_0_uncertainties(QPE_redshifts)
    TDE_redshifts = add_0_uncertainties(TDE_redshifts)

    QPE_data = np.array([QPE_redshifts, QPE_stellar_masses, QPE_mBH])
    TDE_data = np.array([np.concatenate((TDE_redshifts, QPE_redshifts[QPE_and_TDEs])), np.concatenate((TDE_stellar_masses, QPE_stellar_masses[QPE_and_TDEs])), np.concatenate((TDE_mBH, QPE_mBH[QPE_and_TDEs]))])

    redshiftMass([QPE_data, TDE_data], referenceCatalogData=refCat, columns_compare=(1,63,67), save_plot="french_redshift_distribution", fontsize=16, markersize=10,
                        levels=[0.5,0.7,0.9,1],
                        smoothness=7,
                        referenceSmoothness=7,
                        bins=10,
                        kernelDensitiesReference=False,
                        extremums={"param": (0.01,0.1),
                                   "m_star": (9.2+m_star_offsets[0],10.4+m_star_offsets[0]),
                                   }
                )

    QPE_data = np.array([QPE_sersicIndices, QPE_bulgeRatios, QPE_SMSDs, QPE_stellar_masses, QPE_mBH])
    TDE_data = np.array([np.concatenate((TDE_sersicIndices, QPE_sersicIndices[QPE_and_TDEs])), np.concatenate((TDE_bulgeRatios, QPE_bulgeRatios[QPE_and_TDEs])), np.concatenate((TDE_SMSDs, QPE_SMSDs[QPE_and_TDEs])), np.concatenate((TDE_stellar_masses, QPE_stellar_masses[QPE_and_TDEs])), np.concatenate((TDE_mBH, QPE_mBH[QPE_and_TDEs]))])
    myCombinedFinalPlot([QPE_data, TDE_data], referenceCatalogData=refCat, columns_compare=((60,12,68),63,67), save_plot="french_combined_final", fontsize=16, markersize=10,
                        levels=[0.5,0.7,0.9,1],
                        smoothness=12,
                        referenceSmoothness=20,
                        bins=10,
                        kernelDensitiesReference=False,
                        extremums={"n_sersic": (0,5.5),
                                   "bt_ratio": (-0.15,1.05),
                                   "ssmd": (8.2,10.6),
                                   "m_star": (9+m_star_offsets[0],11.25+m_star_offsets[0]),
                                   "m_bh": (4.5,9),
                                   }
                        )
    # Make big plot
    QPE_data  = np.array([QPE_mBH[:,0], QPE_stellar_masses[:,0], QPE_redshifts[:,0], QPE_r50s[:,0], QPE_bulgeRatios[:,0], QPE_sersicIndices[:,0], QPE_SMSDs[:,0]])
    TDE_data = np.array([TDE_mBH[:,0], TDE_stellar_masses[:,0], TDE_redshifts[:,0], TDE_r50s[:,0], TDE_bulgeRatios[:,0], TDE_sersicIndices[:,0], TDE_SMSDs[:,0]])
    double_hosts_data = QPE_data[:,QPE_and_TDEs]
    TDE_data = np.vstack((TDE_data.T, double_hosts_data.T)).T
    myCornerPlot(
        [QPE_data,TDE_data,double_hosts_data],
        labels=["$\log(M_\mathrm{BH})$", "$\log(M_\star)$", "$z$", "$r_{50}$", "$(B/T)_g$", "$n_\mathrm{Sérsic}$", "$\log(\Sigma_{M_\star})$"],
        units=["$[M_\odot]$", "$[M_\odot]$", " ", "$[\mathrm{kpc}]$", " ", " ", "$[M_\odot/\mathrm{kpc}^2]$"],
        smoothness=6,
        markersize=10,
        levels=[0.5,0.7,0.9,1],
        refCat=refCat,
        columns_compare=[67,63,1,59,12,60,68],
        save_plot="french_corner_plot",
        extremums={"$\log(M_\mathrm{BH})$": (4.5,9),
                   "$\log(M_\star)$": (9,11.25),
                   "$(B/T)_g$": (-0.15,1.05),
                   "$r_{50}$": (0,10.5),
                   "$n_\mathrm{Sérsic}$": (0,5.5),
                   "$\log(\Sigma_{M_\star})$": (8.2,10.6) 
                   }
        )
    # Make big plot
    QPE_data  = np.array([QPE_mBH[:,0], QPE_stellar_masses[:,0], QPE_bulgeRatios[:,0], QPE_sersicIndices[:,0], QPE_SMSDs[:,0]])
    TDE_data = np.array([TDE_mBH[:,0], TDE_stellar_masses[:,0], TDE_bulgeRatios[:,0], TDE_sersicIndices[:,0], TDE_SMSDs[:,0]])
    double_hosts_data = QPE_data[:,QPE_and_TDEs]
    TDE_data = np.vstack((TDE_data.T, double_hosts_data.T)).T
    myCornerPlot(
        [QPE_data,TDE_data,double_hosts_data],
        labels=["$\log(M_\mathrm{BH})$", "$\log(M_\star)$", "$(B/T)_g$", "$n_\mathrm{Sérsic}$", "$\log(\Sigma_{M_\star})$"],
        units=["$[M_\odot]$", "$[M_\odot]$", " ", " ", "$[M_\odot/\mathrm{kpc}^2]$"],
        smoothness=6,
        markersize=10,
        levels=[0.5,0.7,0.9,1],
        refCat=refCat,
        columns_compare=[67,63,12,60,68],
        save_plot="french_corner_plot",
        extremums={"$\log(M_\mathrm{BH})$": (4.5,9),
                   "$\log(M_\star)$": (9,11.25),
                   "$(B/T)_g$": (-0.15,1.05),
                   "$r_{50}$": (0,10.5),
                   "$n_\mathrm{Sérsic}$": (0,5.5),
                   "$\log(\Sigma_{M_\star})$": (8.2,10.6) 
                   }
        )
