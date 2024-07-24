"""
Program to compare QPE hosts and TDE hosts "apples to apples" with the LEGACY DESI DECam survey
"""
import pickle
import numpy as np
from ned_wright_cosmology import calculate_cosmo
from utils import print_table, myCornerPlot, toLog, myFinalPlot, myCombinedFinalPlot, recombine_arrays, add_0_uncertainties
import matplotlib.pyplot as plt
from paper_data import *
from download_data import *
import sys
import pandas as pd


QPE_and_TDEs = [4,5,8]


def add_0_uncertainties(a):
    a = np.array(a)
    placeholder = np.zeros((a.shape[0],3))
    placeholder[:,0] = a
    a = placeholder
    return a



def makeLatexTable(names, redshifts, r50s, n_sersics, bt_ratios, ssmds, references=None, filename="latexTable.txt", verbose=False):
    """
    columns are: NAME, REDSHIFT, R_50, N_SERSIC, BT_RATIO, SSMD

    makes a .txt file
    """
    total_string = ""
    each_lines = []
    if references is None:
        references = ["" for i in range(len(names))]

    assert len(names) == len(redshifts) == len(r50s[:,0]) == len(n_sersics[:,0]) == len(bt_ratios[:,0]) == len(ssmds[:,0]) == len(references)
    length = len(names)

    def LatexUncertainty(value):
        return f"${round(value[0],2)}_{r'{-'+str(round(value[1],2))+r'}'}^{r'{+'+str(round(value[2],2))+r'}'}$"

    for i in range(length):
        each_lines.append(names[i]+r"$^{\rm "+references[i]+"}$" + " & " + "$" + str(round(redshifts[i],3)) + "$" + " & " + LatexUncertainty(r50s[i]) + " & " + LatexUncertainty(n_sersics[i]) + " & " + LatexUncertainty(bt_ratios[i]) + " & " + LatexUncertainty(ssmds[i]) + r" \\" + "\n")        

    # swap lines around here easily
    onlyQPEs = np.array(each_lines)[[0,1,2,3,6,7]]
    QPE_TDEs = np.array(each_lines)[[4,5,8]]
    onlyTDEs = np.array(each_lines)[[9,10,11,12,13,14,15,16,17,18]]
    
    total_string += r"\multicolumn{6}{c}{\emph{QPE host galaxies}}\\"+"\n\hline\n\hline\n"
    for line in onlyQPEs:
        total_string += line
        if line != onlyQPEs[-1]:
            total_string += r"\vspace{2pt}"+"\n"
    total_string += "\hline\n"+r"\multicolumn{6}{c}{\emph{TDE+QPE host galaxies}}\\"+"\n\hline\n\hline\n"
    for line in QPE_TDEs:
        total_string += line
        if line != QPE_TDEs[-1]:
            total_string += r"\vspace{2pt}"+"\n"
    total_string += "\hline\n"+r"\multicolumn{6}{c}{\emph{TDE host galaxies}}\\"+"\n\hline\n\hline\n"
    for line in onlyTDEs:
        total_string += line
        if line != onlyTDEs[-1]:
            total_string += r"\vspace{2pt}"+"\n"

    with open(filename, "w") as text_file:
        text_file.write(total_string)
    if verbose:
        print(total_string)





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

    # load reference catalog
    refCat = np.loadtxt("referenceCatalog_modif2.txt")
    fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
    refCat = pd.read_csv("referenceCatalog_modif2.txt", delimiter=" ", header=None, names=fieldnames)

    QPE_r50s, TDE_r50s = QPE_fullData[:,:,5], TDE_fullData[:,:,5]
    QPE_sersicIndices, TDE_sersicIndices = QPE_fullData[:,:,1], TDE_fullData[:,:,1]
    QPE_bulgeRatios, TDE_bulgeRatios = QPE_fullData[:,:,0], TDE_fullData[:,:,0]
    QPE_SMSDs, TDE_SMSDs = QPE_fullData[:,:,2], TDE_fullData[:,:,2]
    QPE_stellar_masses, TDE_stellar_masses = QPE_fullData[:,:,3], TDE_fullData[:,:,3]
    QPE_mBH, TDE_mBH = QPE_fullData[:,:,4], TDE_fullData[:,:,4]

    makeLatexTable(np.concatenate((objects_names, TDE_names)),
                   np.concatenate((QPE_redshifts,TDE_redshifts)),
                   np.concatenate((QPE_r50s,TDE_r50s)),
                   np.concatenate((QPE_sersicIndices,TDE_sersicIndices)),
                   np.concatenate((QPE_bulgeRatios,TDE_bulgeRatios)),
                   np.concatenate((QPE_SMSDs,TDE_SMSDs)),
                   references="abccefddghijklmnnno"
                   )

    QPE_data = np.array([QPE_sersicIndices, QPE_bulgeRatios, QPE_SMSDs, QPE_stellar_masses, QPE_mBH])
    TDE_data = np.array([np.concatenate((TDE_sersicIndices, QPE_sersicIndices[QPE_and_TDEs])), np.concatenate((TDE_bulgeRatios, QPE_bulgeRatios[QPE_and_TDEs])), np.concatenate((TDE_SMSDs, QPE_SMSDs[QPE_and_TDEs])), np.concatenate((TDE_stellar_masses, QPE_stellar_masses[QPE_and_TDEs])), np.concatenate((add_0_uncertainties(TDE_mBH), QPE_mBH[QPE_and_TDEs]))])
    myCombinedFinalPlot([QPE_data, TDE_data], referenceCatalogData=refCat, columns_compare=((60,12,68),63,67), save_plot="combined_final", fontsize=16, markersize=9, levels=[0.3,0.5,0.89,0.9,1])

    # Make big plot
    QPE_data  = np.array([QPE_mBH[:,0], QPE_stellar_masses[:,0], QPE_bulgeRatios[:,0], QPE_r50s[:,0], QPE_sersicIndices[:,0], QPE_SMSDs[:,0]])
    TDE_data = np.array([TDE_mBH, TDE_stellar_masses[:,0], TDE_bulgeRatios[:,0], TDE_r50s[:,0], TDE_sersicIndices[:,0], TDE_SMSDs[:,0]])
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






    sys.exit()

    #From now on, keep only the r-band properties:
    QPE_placeholder_properties = {"n_sersic":[], "r_50":[]}
    for i in range(len(objects)):
        QPE_placeholder_properties["n_sersic"].append(QPE_sersicIndices[i]["r"])
        QPE_placeholder_properties["r_50"].append(QPE_r50s[i]["r"])
    QPE_sersicIndices = QPE_placeholder_properties["n_sersic"]
    QPE_r50s = QPE_placeholder_properties["r_50"]
    QPE_stellar_masses = QPE_stellar_masses_desiProspector
    QPE_stellar_masses = list(QPE_stellar_masses)
    for i in range(len(QPE_stellar_masses)):
        QPE_stellar_masses[i] = tuple(QPE_stellar_masses[i])
    TDE_stellar_masses = TDE_stellar_masses_desiProspector
    
    #Calculate stellar surface densities:
    QPE_stellar_surface_densities = []
    for i in range(len(objects)):
        QPE_stellar_surface_densities.append(stellarMassDensity(QPE_stellar_masses[i], QPE_r50s[i], returnLog=True))
    
    TDE_stellar_surface_densities = []
    TDE_r50s = add_0_uncertainties(TDE_r50s)
    for i in range(len(TDE_stellar_masses)):
        TDE_stellar_surface_densities.append(stellarMassDensity(TDE_stellar_masses[i], TDE_r50s[i], returnLog=True))
    
    #Transform lists into arrays
    QPE_mBH = np.array(QPE_mBH)
    QPE_sersicIndices = np.array(QPE_sersicIndices)
    QPE_stellar_surface_densities = np.array(QPE_stellar_surface_densities)
    QPE_r50s = np.array(QPE_r50s)
    TDE_mBH = np.array(TDE_mBH)
    TDE_sersicIndices = np.array(TDE_sersicIndices)
    TDE_stellar_surface_densities = np.array(TDE_stellar_surface_densities)
    TDE_r50s = np.array(TDE_r50s)

    for i in range(len(QPE_stellar_masses)):
        QPE_stellar_masses[i] = tuple(QPE_stellar_masses[i])
    for i in range(len(TDE_stellar_masses)):
        TDE_stellar_masses[i] = tuple(TDE_stellar_masses[i])
    QPE_stellar_masses = np.array(QPE_stellar_masses)
    TDE_stellar_masses = np.array(TDE_stellar_masses)
    #plot_surfaceStellarMassDensity_stellarMass(QPE_stellar_masses, QPE_stellar_surface_densities, TDE_stellar_masses, TDE_stellar_surface_densities)
    #Make a kind of corner plot, but with fitted gaussians to better illustrate the distributions:
    QPE_data  = np.array([np.log10(QPE_mBH[:,0]), np.log10(QPE_stellar_masses[:,0]), QPE_r50s[:,0], QPE_sersicIndices[:,0], QPE_stellar_surface_densities[:,0]])
    TDE_data = np.array([np.log10(TDE_mBH), np.log10(TDE_stellar_masses[:,0]), TDE_r50s[:,0], TDE_sersicIndices, TDE_stellar_surface_densities[:,0]])
    double_hosts_data = QPE_data[:,QPE_and_TDEs]
