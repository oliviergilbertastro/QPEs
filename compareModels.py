#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
import numpy as np
from galight_modif.tools.plot_tools import total_compare
from download_data import objects
from ned_wright_cosmology import calculate_cosmo

def compareModels(ra_dec, models=["None", "AGN", "Bulge", "Bulge+AGN"], band="i", stellar_mass=5E9, z=0, verbose=True, returnData=False):
        
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

    bics = []
    if verbose:
        print("BIC:")
        print("-------------------------------------")
    for i in range(len(fitting_run_results)):
        bics.append(fitting_run_results[i].fitting_seq.bic)
        if verbose:
            print(types[i], ":", fitting_run_results[i].fitting_seq.bic)
    best_model_index = bics.index(np.min(bics))
    if verbose:
        print("-------------------------------------")
        print("Best model:", types[best_model_index])
    #Calculate the Sersic index + uncertainties
    chain = fitting_run_results[best_model_index].samples_mcmc
    lo, mid, hi = np.percentile(chain[:, 1],16), np.percentile(chain[:, 1],50), np.percentile(chain[:, 1],84)
    plus, minus = (hi-mid), (mid-lo)
    sersic_index_data = [mid, minus, plus]
    sersic_index_str = (f"{mid:.2f}_"+r"{-"+f"{minus:.2f}"+r"}^{+"+f"{plus:.2f}"+r"}")

    #Calculate the Sersic half-light radius + uncertainties:
    lo, mid, hi = np.percentile(chain[:, 0],16), np.percentile(chain[:, 0],50), np.percentile(chain[:, 0],84)
    plus, minus = (hi-mid), (mid-lo)

    theta_rad = np.array([mid, minus, plus])*np.pi/(180*3600)
    cosmology_params = calculate_cosmo(z, H0=67, Omega_m=0.31, Omega_v=0.69)
    kpc_per_arcsec = cosmology_params["kpc_DA"]
    half_light_radius = np.array([mid, minus, plus])*kpc_per_arcsec
    try:
        index = objects.index(ra_dec)
    except:
        index = comparisons.index(ra_dec)
    #print(index)

    stellar_density_str = "-"
    smd = 0
    smd_data = 0
    if stellar_mass != None:
        #print("r50", [mid, minus, plus])
        smd = np.array(stellarMassDensity(stellar_mass, half_light_radius))
        smd_data = smd.copy()
        smd = np.log10([smd[0], smd[0]-smd[1], smd[0]+smd[2]])
        smd = np.array([smd[0], smd[0]-smd[1], smd[2]-smd[0]])
        stellar_density_str = (f"{smd[0]:.2f}_"+r"{-"+f"{smd[1]:.2f}"+r"}^{+"+f"{smd[2]:.2f}"+r"}")
    
    if returnData:
        return sersic_index_data, smd

    return types[bics.index(np.min(bics))] + " " + "?"*(len(models)-len(bics)) + f' n = {sersic_index_str}      log(\Sigma_star) = {stellar_density_str}'


def stellarMassDensity(M_star, r50):
    '''
    Calculate the stellar surface mass density \Sigma_{M_\ast} from the total stellar mass and the half-light radius

    M_star: stellar mass in solar masses -> tuple/list/array such as (mass, errlo, errhi)
    r50: half-light radius in kpc -> tuple/list/array such as (radius, errlo, errhi)

    returns [value, errlo, errhi]
    '''
    try:
        res = M_star[0]/r50[0]**2
        errlo = res*np.sqrt((M_star[1]/M_star[0])**2+(2*r50[2]/r50[0])**2)
        errhi = res*np.sqrt((M_star[2]/M_star[0])**2+(2*r50[1]/r50[0])**2)
        return [res, errlo, errhi]
    except:
        #In the exception where the user did not input uncertainties, the value will still be calculated.
        return M_star/r50[0]**2


import matplotlib.pyplot as plt
import corner
from scipy.stats import norm

def plot_sersicIndex_mBH(QPEmBH, QPEsersicIndices, TDEmBH, TDEsersicIndices, verbose=False):
    '''
    Plots the sersic indices as a function of the black hole masses
    '''
    ax1 = plt.subplot(111)
    ax1.errorbar(QPEmBH[:,0], QPEsersicIndices[:,0], yerr=[QPEsersicIndices[:,1],QPEsersicIndices[:,2]], xerr=[QPEmBH[:,1],QPEmBH[:,2]], fmt='D', color='green', label='QPE hosts')
    ax1.scatter(TDEmBH, TDEsersicIndices, marker='o', color='orange', label='TDE hosts')
    ax1.set_xscale("log")
    ax1.yaxis.set_tick_params(labelsize=15)
    ax1.xaxis.set_tick_params(labelsize=15)
    plt.xlabel(r'$M_\mathrm{BH}$ [$M_\odot$]', size=17)
    plt.ylabel(r'Galaxy Sérsic Index', size=17)
    plt.legend(fontsize=15)
    plt.show()


    #Plot Sérsic difference:
    # Fit a normal distribution to the data:
    mu1, std1 = norm.fit(QPEsersicIndices[:,0])
    mu2, std2 = norm.fit(TDEsersicIndices)
    alldata = np.concatenate((QPEsersicIndices[:,0], TDEsersicIndices))
    mini, maxi = np.min(alldata), np.max(alldata)
    mini, maxi = 0.5, 8.5
    # Plot the histogram.
    #plt.hist(QPEsersicIndices[:,0], bins=25, density=False, alpha=0.5, color='b', label='QPE')
    #plt.hist(TDEsersicIndices, bins=25, density=False, alpha=0.5, color='r', label='TDE')
    plt.hist(QPEsersicIndices[:,0], range=(mini, maxi), bins=50, density=False, alpha=0.4, color='b')
    plt.hist(TDEsersicIndices, range=(mini, maxi), bins=50, density=False, alpha=0.4, color='r')
    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p1 = norm.pdf(x, mu1, std1)
    plt.plot(x, p1, color='b', linewidth=2, linestyle="dashed", label='QPE hosts')
    p2 = norm.pdf(x, mu2, std2)
    plt.plot(x, p2, color='r', linewidth=2, linestyle="dashed", label='TDE hosts')
    if True: #ADD SDSS CONTROL
        p3 = norm.pdf(x, 1.21, 1.29)
        plt.plot(x, p3, color='black', linewidth=2, linestyle="dashed", label="SDSS controls")
    plt.xlabel("Sérsic index", fontsize=17)
    plt.ylabel("Number of hosts", fontsize=17)
    plt.legend(fontsize=14)
    if verbose:
        print("mu1:", mu1)
        print("mu2:", mu2)
        print("std1:", std1)
        print("std2:", std2)
    plt.show()

    #fig1 = plt.figure(1)
    #plot1 = corner.corner(np.array([np.log10(QPEmBH[:,0]), QPEsersicIndices[:,0]]).T, labels=["Black hole mass", "Sérsic index"], show_titles=True, fig=fig1)
    #plt.suptitle("QPE hosts")
    #fig2 = plt.figure(2)
    #plot2 = corner.corner(np.array([np.log10(TDEmBH), TDEsersicIndices]).T, labels=["Black hole mass", "Sérsic index"], show_titles=True, fig=fig2)
    #plt.suptitle("TDE hosts")
    #plt.show()
    return

def plot_surfaceStellarMassDensity_mBH(QPEmBH, QPEstellarDensities, TDEmBH, TDEstellarDensities, verbose=False):
    uplims = TDEstellarDensities[:,1] == "upper limit"
    #TDE_stellarDensities[np.array(TDEstellarDensities)[:,1] == "upper limit"][:,1:] = (0, 0)
    for i in range(len(uplims)):
        if uplims[i]:
            TDEstellarDensities[i,1:] = (0.1,0)
    TDEstellarDensities = np.array(TDEstellarDensities, dtype=float)
    #plt.errorbar(Edd_ratio, a_ox, [a_ox_err_lo, a_ox_err_hi], [Edd_ratio_err_lo, Edd_ratio_err_hi], fmt='o', markeredgecolor='black', markeredgewidth=2, uplims=uplims, label=r'$\alpha_\mathrm{OX}$')
    ax1 = plt.subplot(111)
    ax1.errorbar(QPEmBH[:,0], QPEstellarDensities[:,0], yerr=[QPEstellarDensities[:,1],QPEstellarDensities[:,2]], xerr=[QPEmBH[:,1],QPEmBH[:,2]], fmt='D', color='green', label='QPE hosts')
    ax1.errorbar(TDEmBH, TDEstellarDensities[:,0], yerr=[TDEstellarDensities[:,1],TDEstellarDensities[:,2]], uplims=uplims, markeredgewidth=2, fmt='o', color='orange', label='TDE hosts')
    ax1.set_xscale("log")
    ax1.yaxis.set_tick_params(labelsize=15)
    ax1.xaxis.set_tick_params(labelsize=15)
    plt.xlabel(r'$M_\mathrm{BH}$ [$M_\odot$]', size=17)
    plt.ylabel(r'$\log(\Sigma_{M_\star})$', size=17)
    plt.legend(fontsize=15)
    plt.show()

    
    if input("Include the upper limits as data points in histogram? [y/n]") != "y":
        TDEstellarDensities_withoutUplims = []
        for i in range(len(uplims)):
            if not uplims[i]:
                TDEstellarDensities_withoutUplims.append(TDEstellarDensities[i])
        TDEstellarDensities = np.array(TDEstellarDensities_withoutUplims)

    #Plot Sérsic difference:
    # Fit a normal distribution to the data:
    mu1, std1 = norm.fit(QPEstellarDensities[:,0])
    mu2, std2 = norm.fit(TDEstellarDensities[:,0])
    alldata = np.concatenate((QPEstellarDensities[:,0], TDEstellarDensities[:,0]))
    mini, maxi = np.min(alldata), np.max(alldata)
    # Plot the histogram.
    #plt.hist(QPEsersicIndices[:,0], bins=25, density=False, alpha=0.5, color='b', label='QPE')
    #plt.hist(TDEsersicIndices, bins=25, density=False, alpha=0.5, color='r', label='TDE')
    plt.hist(QPEstellarDensities[:,0], range=(mini, maxi), bins=50, density=False, alpha=0.4, color='b')
    plt.hist(TDEstellarDensities[:,0], range=(mini, maxi), bins=50, density=False, alpha=0.4, color='r')
    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p1 = norm.pdf(x, mu1, std1)
    plt.plot(x, p1, color='b', linewidth=2, linestyle="dashed", label='QPE hosts')
    p2 = norm.pdf(x, mu2, std2)
    plt.plot(x, p2, color='r', linewidth=2, linestyle="dashed", label='TDE hosts')
    plt.xlabel(r'$\log(\Sigma_{M_\star})$', fontsize=17)
    plt.ylabel("Number of hosts", fontsize=17)
    plt.legend(fontsize=14)
    if verbose:
        print("mu1:", mu1)
        print("mu2:", mu2)
        print("std1:", std1)
        print("std2:", std2)
    plt.show()
    return

def plot_sersicIndex_surfaceStellarMassDensity(QPEsersicIndices, QPEstellarDensities, TDEsersicIndices, TDEstellarDensities):
    uplims = TDEstellarDensities[:,1] == "upper limit"
    #TDE_stellarDensities[np.array(TDEstellarDensities)[:,1] == "upper limit"][:,1:] = (0, 0)
    for i in range(len(uplims)):
        if uplims[i]:
            TDEstellarDensities[i,1:] = (0.1,0)
    TDEstellarDensities = np.array(TDEstellarDensities, dtype=float)
    #plt.errorbar(Edd_ratio, a_ox, [a_ox_err_lo, a_ox_err_hi], [Edd_ratio_err_lo, Edd_ratio_err_hi], fmt='o', markeredgecolor='black', markeredgewidth=2, uplims=uplims, label=r'$\alpha_\mathrm{OX}$')
    ax1 = plt.subplot(111)
    ax1.errorbar(QPEsersicIndices[:,0], QPEstellarDensities[:,0], yerr=[QPEstellarDensities[:,1],QPEstellarDensities[:,2]], xerr=[QPEsersicIndices[:,1],QPEsersicIndices[:,2]], fmt='D', color='green', label='QPE hosts')
    ax1.errorbar(TDEsersicIndices, TDEstellarDensities[:,0], yerr=[TDEstellarDensities[:,1],TDEstellarDensities[:,2]], uplims=uplims, markeredgewidth=2, fmt='o', color='orange', label='TDE hosts')
    #ax1.set_xscale("log")
    ax1.yaxis.set_tick_params(labelsize=15)
    ax1.xaxis.set_tick_params(labelsize=15)
    plt.xlabel(r'Sérsic index', size=17)
    plt.ylabel(r'$\log(\Sigma_{M_\star})$', size=17)
    plt.legend(fontsize=15)
    plt.show()



#Import all the data from another file for convenience
from paper_data import TDE_mBH, TDE_sersicIndices, TDE_stellarDensities, QPE_redshifts, TDE_redshifts, QPE_mBH, QPE_stellar_masses, TDE_stellar_masses


#filter bands used to fit each object to simplify the user input, some other bands were also fitted and/or could be fitted, but these ones work well
QPE_bands_list = ["i", "i", "i", "i", "z", "i", "i", "r", "r"]
TDE_bands_list = ["r", "r", "g",]

from download_data import objects, comparisons, objects_names, comparisons_names

if __name__ == "__main__":
    if input("See options to compare QPE and TDE fits? [y/n]") == "y":
        if input("Compare fits for a specific QPE host galaxy? [y/n]") == "y":
            objID = int(input(f"Enter the object ID you want to load [0-{len(objects)-1}]:\n"))
            band = input("Enter the filter band you want to load [g,r,i,z]:\n")
            compareModels(objects[objID], band=band, stellar_mass=QPE_stellar_masses[objID], z=QPE_redshifts[objID], verbose=True)

        if input("See best model for all QPE host galaxies? [y/n]") == "y":
            print("Best models:")
            print("-------------------------------------------------")
            for i in range(len(objects)):
                print(f"{objects_names[i]}: {compareModels(objects[i], band=QPE_bands_list[i], stellar_mass=QPE_stellar_masses[i], z=QPE_redshifts[i], models=['None', 'AGN'], verbose=False)}")
            print("-------------------------------------------------")
            #The Sérsic index of AT 2019qiz is fitted to be 3. something, but according to https://arxiv.org/pdf/2006.0245, it is 5.2 when fitted with a point source AGN

        if input("Compare fits for a specific TDE host galaxy? [y/n]") == "y":
            objID = int(input(f"Enter the object ID you want to load [0-{len(comparisons)-1}]:\n"))
            band = input("Enter the filter band you want to load [g,r,i,z]:\n")
            compareModels(comparisons[objID], band=band, verbose=True, stellar_mass=TDE_stellar_masses[objID], z=TDE_redshifts[objID])

        if input("See best model for all TDE host galaxies? [y/n]") == "y":
            print("Best models:")
            print("-------------------------------------------------")
            for i in range(len(comparisons)):
                print(f"{comparisons_names[i]}: {compareModels(comparisons[i], band=TDE_bands_list[i], stellar_mass=TDE_stellar_masses[i], z=TDE_redshifts[i], models=['None'], verbose=False)}")
            print("-------------------------------------------------")

    if input("Plot log(mBH)-Sérsic index for host galaxies comparison? [y/n]") == "y":
        QPE_mBH = np.array(QPE_mBH)
        QPE_sersicIndices = []
        for i in range(len(objects)):
            sersic, stellarDensity = compareModels(objects[i], band=QPE_bands_list[i], stellar_mass=QPE_stellar_masses[i], z=QPE_redshifts[i], models=['None', 'AGN'], verbose=False, returnData=True)
            QPE_sersicIndices.append(sersic)
        QPE_sersicIndices = np.array(QPE_sersicIndices)
        plot_sersicIndex_mBH(QPE_mBH, QPE_sersicIndices, TDE_mBH, TDE_sersicIndices)

    if input("Plot log(mBH)-Stellar mass surface density for host galaxies comparison? [y/n]") == "y":
        QPE_mBH = np.array(QPE_mBH)
        used_QPE_mBH = [] #only temporary until all QPE_stellarMasses are found, then this will be the same as QPE_mBH
        QPE_stellarDensities = []
        for i in range(len(objects)):
            sersic, stellarDensity = compareModels(objects[i], band=QPE_bands_list[i], stellar_mass=QPE_stellar_masses[i], z=QPE_redshifts[i], models=['None', 'AGN'], verbose=False, returnData=True)
            if QPE_stellar_masses[i] != None:
                QPE_stellarDensities.append(stellarDensity)
                used_QPE_mBH.append(QPE_mBH[i])
        used_QPE_mBH = np.array(used_QPE_mBH)
        QPE_stellarDensities = np.array(QPE_stellarDensities)
        TDE_mBH = np.array(TDE_mBH)
        TDE_stellarDensities = np.array(TDE_stellarDensities)
        plot_surfaceStellarMassDensity_mBH(used_QPE_mBH, QPE_stellarDensities, TDE_mBH, TDE_stellarDensities)

    if input("Plot Sérsic index - Stellar mass surface density for host galaxies comparison? [y/n]") == "y":
        QPE_mBH = np.array(QPE_mBH)
        used_QPE_mBH = [] #only temporary until all QPE_stellarMasses are found, then this will be the same as QPE_mBH
        QPE_sersicIndices = []
        QPE_stellarDensities = []
        for i in range(len(objects)):
            sersic, stellarDensity = compareModels(objects[i], band=QPE_bands_list[i], stellar_mass=QPE_stellar_masses[i], z=QPE_redshifts[i], models=['None', 'AGN'], verbose=False, returnData=True)
            if QPE_stellar_masses[i] != None:
                QPE_stellarDensities.append(stellarDensity)
                QPE_sersicIndices.append(sersic)
                used_QPE_mBH.append(QPE_mBH[i])
        used_QPE_mBH = np.array(used_QPE_mBH)
        QPE_sersicIndices = np.array(QPE_sersicIndices)
        QPE_stellarDensities = np.array(QPE_stellarDensities)
        TDE_mBH = np.array(TDE_mBH)
        TDE_stellarDensities = np.array(TDE_stellarDensities)
        plot_sersicIndex_surfaceStellarMassDensity(QPE_sersicIndices, QPE_stellarDensities, TDE_sersicIndices, TDE_stellarDensities)