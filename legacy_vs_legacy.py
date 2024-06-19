"""
Program to compare QPE hosts and TDE hosts "apples to apples" with the LEGACY DESI DECam survey
"""
import pickle
import numpy as np
from ned_wright_cosmology import calculate_cosmo
from utils import print_table
import matplotlib.pyplot as plt
from paper_data import *
from download_data import *
import sys
from prospector.prospector_myFits import QPE_stellar_masses_desiProspector

def get_n_and_r50(objID, model="None", band="i", survey="DESI", redshift=0, qpe_oder_tde="QPE"):
    if qpe_oder_tde == "QPE":
        picklename = f"{objects_names[objID]}_{band}-band_{model}_{survey}.pkl"
    else:
        picklename = f"{TDE_names[objID]}_{band}-band_{model}_{survey}.pkl"
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        raise NameError("There is no picklerun with this name.")

    #Calculate the Sersic index + uncertainties
    chain = fitting_run_result.samples_mcmc
    params = fitting_run_result.param_mcmc
    if "R_sersic" in params and "n_sersic" in params:
        lo, mid, hi = np.percentile(chain[:, 1],16), np.percentile(chain[:, 1],50), np.percentile(chain[:, 1],84)
        plus, minus = (hi-mid), (mid-lo)
        sersic_index_data = [mid, minus, plus]

        #Calculate the Sersic half-light radius + uncertainties:
        lo, mid, hi = np.percentile(chain[:, 0],16), np.percentile(chain[:, 0],50), np.percentile(chain[:, 0],84)
        plus, minus = (hi-mid), (mid-lo)
        r50_data = [mid, minus, plus]
    else:
        sersic_index_data = [fitting_run_result.final_result_galaxy[0]["n_sersic"], 0, 0]
        r50_data = [fitting_run_result.final_result_galaxy[0]["R_sersic"], 0, 0]
    #Convert r50 from arcsec to kpc
    cosmology_params = calculate_cosmo(redshift, H0=67, Omega_m=0.31, Omega_v=0.69)
    kpc_per_arcsec = cosmology_params["kpc_DA"]
    r50_data = np.array(r50_data)*kpc_per_arcsec

    magnitude = fitting_run_result.final_result_galaxy[0]['magnitude']

    return sersic_index_data, r50_data, magnitude


def add_0_uncertainties(a):
    a = np.array(a)
    placeholder = np.zeros((a.shape[0],3))
    placeholder[:,0] = a
    a = placeholder
    return a

def stellarMassDensity(M_star, r50, returnLog=False):
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
        if returnLog:
            return [np.log10(res), errlo/(res*np.log(10)), errhi/(res*np.log(10))]
        return [res, errlo, errhi]
    except:
        #In the exception where the user did not input uncertainties, the value will still be calculated.
        if returnLog:
            return np.log10(M_star/r50[0]**2)
        return M_star/r50[0]**2#/(2*np.pi) #is there a constant 2pi we need to divide by???

def checkWhichFiltersWork(list_of_dicts, qpe_oder_tde="QPE"):
    working = {"g":[], "r":[], "i":[], "z":[]}
    for i in range(len(list_of_dicts)):
        for band in "griz":
            try:
                idc = list_of_dicts[i][band]
                working[band].append("\x1b[32mY\x1b[0m")
            except:
                working[band].append("\x1b[31mN\x1b[0m")
    names = objects_names if qpe_oder_tde=="QPE" else TDE_names
    print_table(np.array([names, working["g"], working["r"], working["i"], working["z"]]).T,
                header=["Name", "g", "r", "i", "z"],
                title="Working filters",
                borders=2,
                override_length=[15, 1, 1, 1, 1],
                )
    return

def printPropertyAcrossFilters(list_of_dicts, name_of_property="Name of property", round_to_n_decimals=2, qpe_oder_tde="QPE"):
    properties = {"g":[], "r":[], "i":[], "z":[]}
    for i in range(len(list_of_dicts)):
        for band in "griz":
            try:
                properties[band].append(f"{list_of_dicts[i][band][0]:0.2f}")
            except:
                try:
                    properties[band].append(f"{list_of_dicts[i][band]:0.2f}")
                except:
                    properties[band].append("-")
    names = objects_names if qpe_oder_tde=="QPE" else TDE_names
    print_table(np.array([names, properties["g"], properties["r"], properties["i"], properties["z"]]).T,
                header=["Name", "g", "r", "i", "z"],
                title=name_of_property,
                borders=2,
                )
    return




def myCornerPlot(data, labels=None, fontsize=15):
    """
    data should be [data_set1, data_set2, ...] each containing multiple parameters
    """
    for i in range(len(data)-1):
        assert len(data[i]) == len(data[i+1])
    # Create the plot axes:
    fig = plt.figure(figsize=(10,8))
    plot_size = len(data[0])
    hist_axes = []
    corner_axes = []
    for i in range(plot_size):
        hist_axes.append(plt.subplot(plot_size,plot_size,i*plot_size+i+1))
        corner_axes.append([])
        for k in range(plot_size-(i+1)):
            if i == 0:
                corner_axes[i].append(plt.subplot(plot_size,plot_size,(i+k+1)*plot_size+(i+1),sharex=hist_axes[i]))
            else:
                corner_axes[i].append(plt.subplot(plot_size,plot_size,(i+k+1)*plot_size+(i+1),sharex=hist_axes[i],sharey=corner_axes[i-1][k+1]))
            if k != plot_size-(i+1)-1:
                corner_axes[i][k].get_xaxis().set_visible(False)
            if i != 0:
                corner_axes[i][k].get_yaxis().set_visible(False)
            corner_axes[i][k].xaxis.set_tick_params(labelsize=fontsize-2)
            corner_axes[i][k].yaxis.set_tick_params(labelsize=fontsize-2)
        if i == plot_size-1:
            hist_axes[i].get_yaxis().set_visible(False)
            hist_axes[i].xaxis.set_tick_params(labelsize=fontsize-2)
        else:
            hist_axes[i].get_xaxis().set_visible(False)
            hist_axes[i].get_yaxis().set_visible(False)

    # Show data in each plot:

    from sklearn.neighbors import KernelDensity
    #Plot kernel histograms:
    for i in range(plot_size):
        if labels is not None:
            hist_axes[i].set_title(labels[i], fontsize=fontsize)
        x_min, x_max = np.min(data[0][i]), np.max(data[0][i])
        for j in range(len(data)):
            x_min = np.min(data[j][i]) if x_min > np.min(data[j][i]) else x_min
            x_max = np.max(data[j][i]) if x_max < np.max(data[j][i]) else x_max
        for j in range(len(data)-1):
            X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
            kde = KernelDensity(kernel="gaussian", bandwidth=0.75).fit(data[j][i][:,np.newaxis])
            log_dens = kde.score_samples(X_plot)
            hist_axes[i].fill_between(X_plot[:, 0], np.exp(log_dens), fc=["blue","red","orange"][j%3], alpha=[0.4,0.4][j])

 
    for i in range(plot_size):
        for k in range(len(corner_axes[i])):
            for j in range(len(data)):
                corner_axes[i][k].plot(data[j][i], data[j][i+k+1], ["o","d","*"][j%3], color=["blue","red","orange"][j%3])
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.97, top=0.94, wspace=0, hspace=0)
    plt.show()
    return
















survey = "DESI"
if __name__ == "__main__":
    #survey = input("Which survey? [default=DESI]\n")
    survey = "DESI" if survey == "" else survey
#Do this part so other programs can load it (especially the magnitudes from prospector)
#First, load the TDE sersic indices and half-light radii into arrays or list, idc:
QPE_sersicIndices = []
QPE_r50s = []
QPE_magnitudes = []
QPE_unreddenedMagnitudes = []
for i in range(len(objects)):
    QPE_sersicIndices.append({})
    QPE_r50s.append({})
    QPE_magnitudes.append({})
    QPE_unreddenedMagnitudes.append({})
    for band in "griz":
        try:
            n, r50, mag = get_n_and_r50(i, objects_types[i], redshift=QPE_redshifts[i], band=band, survey=survey)
            QPE_sersicIndices[-1][band] = n
            QPE_r50s[-1][band] = r50
            QPE_magnitudes[-1][band] = mag
            QPE_unreddenedMagnitudes[-1][band] = mag - QPE_extinction[objects_names[i]][band]
        except:
            pass

TDE_sersicIndices = []
TDE_r50s = []
TDE_magnitudes = []
TDE_unreddenedMagnitudes = []
for i in range(len(TDE_coords)):
    TDE_sersicIndices.append({})
    TDE_r50s.append({})
    TDE_magnitudes.append({})
    TDE_unreddenedMagnitudes.append({})
    for band in "griz":
        try:
            n, r50, mag = get_n_and_r50(i, "None", redshift=TDE_redshifts[i], band=band, survey=survey, qpe_oder_tde="TDE")
            TDE_sersicIndices[-1][band] = n
            TDE_r50s[-1][band] = r50
            TDE_magnitudes[-1][band] = mag
            TDE_unreddenedMagnitudes[-1][band] = mag - QPE_extinction[objects_names[i]][band] #Change this
        except:
            pass

if __name__ == "__main__":
    #print(objects[i], objects_types[i], band)
    checkWhichFiltersWork(QPE_sersicIndices)
    printPropertyAcrossFilters(QPE_sersicIndices, "Sérsic Index")
    printPropertyAcrossFilters(QPE_r50s, "Sérsic half-light radius (kpc)")
    printPropertyAcrossFilters(QPE_magnitudes, "Magnitude")
    printPropertyAcrossFilters(QPE_unreddenedMagnitudes, "Dereddened magnitude")

    checkWhichFiltersWork(TDE_sersicIndices, qpe_oder_tde="TDE")
    printPropertyAcrossFilters(TDE_sersicIndices, "Sérsic Index", qpe_oder_tde="TDE")

    #From now on, keep only the r-band properties:
    QPE_placeholder_properties = {"n_sersic":[], "r_50":[]}
    TDE_placeholder_properties = {"n_sersic":[], "r_50":[]}
    for i in range(len(objects)):
        QPE_placeholder_properties["n_sersic"].append(QPE_sersicIndices[i]["r"])
        QPE_placeholder_properties["r_50"].append(QPE_r50s[i]["r"])
        TDE_placeholder_properties["n_sersic"].append(TDE_sersicIndices[i]["r"])
        TDE_placeholder_properties["r_50"].append(TDE_r50s[i]["r"])
    QPE_sersicIndices = QPE_placeholder_properties["n_sersic"]
    QPE_r50s = QPE_placeholder_properties["r_50"]
    TDE_sersicIndices = TDE_placeholder_properties["n_sersic"]
    TDE_r50s = TDE_placeholder_properties["r_50"]

    #Transform lists into arrays
    QPE_sersicIndices = np.array(QPE_sersicIndices)
    QPE_r50s = np.array(QPE_r50s)
    TDE_sersicIndices = np.array(TDE_sersicIndices)
    TDE_r50s = np.array(TDE_r50s)

    QPE_data  = np.array([QPE_r50s[:,0], QPE_sersicIndices[:,0]])
    TDE_data = np.array([TDE_r50s[:,0], TDE_sersicIndices[:,0]])
    double_hosts_data = QPE_data[:,[4,8]]
    TDE_data = np.vstack((TDE_data.T, double_hosts_data.T)).T
    myCornerPlot([
        QPE_data,TDE_data,double_hosts_data
        ],
                 labels=["$r_{50}$","$n_\mathrm{Sérsic}$"])

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
    double_hosts_data = QPE_data[:,[4,8]]
