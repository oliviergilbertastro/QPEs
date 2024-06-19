'''
Get SDSS photometry data of objects, fit galaxy models to it using nested sampling to obtain the total stellar mass in the galaxy
'''


import os
import sys
import copy
import sys
import parse
import astropy.units
#Let's just go in the parent directory since all of our imports are from there (much less complicated)
parent_dir = parse.parse("{}/prospector", os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(parent_dir)
import download_data, paper_data
os.environ["SPS_HOME"] = "/Users/oliviergilbert/Desktop/QPEs/fsps-master" #Need to download fsps from https://github.com/cconroy20/fsps and put it locally, change the path in .bash_profile to point to it
import astropy.table
import fsps
import dynesty
import sedpy
import h5py, astropy
import numpy as np
import astroquery
from sedpy.observate import load_filters
from prospect.utils.obsutils import fix_obs
from download_data import objects_names, objects
if __name__ != "prospector.prospector_myFits":
    from cleanEverything_v2 import QPE_magnitudes, QPE_unreddenedMagnitudes
from paper_data import QPE_redshifts
from prospect.sources import CSPSpecBasis
from prospect.io import read_results as reader
import matplotlib.pyplot as plt
from prospect.plotting import corner
from prospect.plotting.utils import best_sample

def fit_SED(objID, bands="griz", redshift=0, magnitudes_dict=None):
    """
    Fits SED using prospector
    """
    
    filters = load_filters([f"sdss_{b}0" for b in bands])
    working_bands = ""
    for b in bands:
        try:
            idc = magnitudes_dict[0][f"cModelMag_{b}"]
            working_bands += b
        except:
            pass
    filters = load_filters([f"sdss_{b}0" for b in working_bands])
    print(working_bands)

    cat = copy.copy(magnitudes_dict)
    for b in working_bands:
        print(b, cat[0][f"cModelMag_{b}"])
    maggies = np.array([10**(-0.4 * cat[0][f"cModelMag_{b}"]) for b in working_bands])
    magerr = np.array([cat[0][f"cModelMagErr_{b}"] for b in working_bands])
    magerr = np.hypot(magerr, 0.05)

    obs = dict(wavelength=None, spectrum=None, unc=None, redshift=redshift,
            maggies=maggies, maggies_unc=magerr * maggies / 1.086, filters=filters)
    obs = fix_obs(obs)
    #print(obs)


    from prospect.models.templates import TemplateLibrary
    from prospect.models import SpecModel
    model_params = TemplateLibrary["parametric_sfh"]
    #model_params.update(TemplateLibrary["nebular"]) #This is to add nebular emission lines (I don't care about that)
    #model_params.update(TemplateLibrary["continuity_psb_sfh"]) #This is to add additional parameters
    model_params["zred"]["init"] = obs["redshift"]
    #Let's fix some of the parameters:
    #model_params["dust2"]["isfree"] = False
    #model_params["tage"]["isfree"] = False
    #model_params["tau"]["isfree"] = False

    #Fit directly in SURVIVING stellar masses
    model_params["mass_units"]=dict(init="mstar", isfree=False, N=1)


    print("-------------------------------------")
    for k in model_params.keys():
        string = "\x1b[32mfree\x1b[0m" if model_params[k]["isfree"] else "\x1b[31mfixed\x1b[0m"
        print(f"{k}: {string}")
    print("-------------------------------------")
    


    model = SpecModel(model_params)
    print("Number of free parameters:", len(model.free_params))
    #assert len(model.free_params) == 5
    #print(model)

    noise_model = (None, None)

    
    sps = CSPSpecBasis(zcontinuous=1)
    #print(sps.ssp.libraries)

    current_parameters = ",".join([f"{p}={v}" for p, v in zip(model.free_params, model.theta)])
    print(current_parameters)
    spec, phot, mfrac = model.predict(model.theta, obs=obs, sps=sps)
    print("mfrac:",mfrac)
    surviving_mass = np.sum(model.params["mass"]) * mfrac
    print("surviving_mass:", surviving_mass)
    #print(phot / obs["maggies"])
    #with open(f"prospector/data/{objects_names[objID]}_LEGACY_mfrac.txt", "w") as text_file:
    #    text_file.write(str(mfrac))

    from prospect.fitting import lnprobfn, fit_model
    fitting_kwargs = dict(nlive_init=400, nested_method="rwalk", nested_target_n_effective=1000, nested_dlogz_init=0.05)
    output = fit_model(obs, model, sps, optimize=False, dynesty=True, lnprobfn=lnprobfn, noise=noise_model, **fitting_kwargs)
    result, duration = output["sampling"]
    from prospect.io import write_results as writer
    hfile = f"prospector/fits_data/{objects_names[objID]}_survMass_prospector_SED.h5"
    writer.write_hdf5(hfile, {}, model, obs,
                    output["sampling"][0], None,
                    sps=sps,
                    tsample=output["sampling"][1],
                    toptimize=0.0)
    #uncomment the rest to run the python file only once to make all files
    #objID = objID + 1
    #fit_SED(objID, bands="griz", redshift=QPE_redshifts[objID], magnitudes_dict=magnitudes_dicts[objID])

def read_SED(objID):
    try:
        # Directly fitting the surviving mass
        hfile = f"prospector/fits_data/{objects_names[objID]}_survMass_prospector_SED.h5"
        out, out_obs, out_model = reader.results_from(hfile)
    except:
        # Taking the fitted total formed mass and then multiplying it by the surviving fraction approximation
        hfile = f"prospector/fits_data/{objects_names[objID]}_prospector_SED.h5"
        out, out_obs, out_model = reader.results_from(hfile)
    for k in out.keys():
        print(k, ":", out[k])
    #Plot the posterior's corner plot
    nsamples, ndim = out["chain"].shape
    cfig, axes = plt.subplots(ndim, ndim, figsize=(10,9))
    axes = corner.allcorner(out["chain"].T, out["theta_labels"], axes, weights=out["weights"], color="royalblue", show_titles=True)
    pbest = best_sample(out)
    corner.scatter(pbest[:, None], axes, color="firebrick", marker="o")
    plt.show()

    #Making my own version with log(Mass):
    data = copy.copy(out)
    mass = data["chain"][:,0]
    logMass = np.log10(mass)
    data["chain"][:,0] = logMass
    data["theta_labels"][0] = r"$\log(M_\star)$"
    data["theta_labels"][1] = r"$\log(Z_\star/Z_\odot)$"
    cfig, axes = plt.subplots(ndim, ndim, figsize=(10,9))
    axes = corner.allcorner(data["chain"].T, out["theta_labels"], axes, weights=out["weights"], color="royalblue", show_titles=True)
    pbest = best_sample(data)
    print(pbest)
    corner.scatter(pbest[:, None], axes, color="firebrick", marker="o")
    plt.show()

    #Plot the observed SED of the spectrum and SED of the highest probability posterior sample
    sfig, saxes = plt.subplots(2, 1, gridspec_kw=dict(height_ratios=[1, 4]), sharex=True)
    ax = saxes[1]
    pwave = np.array([f.wave_effective for f in out_obs["filters"]])

    # plot the data
    ax.plot(pwave, out_obs["maggies"], linestyle="", marker="o", color="k")
    ax.errorbar(pwave,  out_obs["maggies"], out_obs["maggies_unc"], linestyle="", color="k", zorder=10)
    ax.set_ylabel(r"$f_\nu$ (maggies)")
    ax.set_xlabel(r"$\lambda$ (AA)")
    ax.set_xlim(3e3, 1e4)
    ax.set_ylim(out_obs["maggies"].min() * 0.1, out_obs["maggies"].max() * 5)
    ax.set_yscale("log")

    # get the best-fit SED
    bsed = out["bestfit"]
    ax.plot(bsed["restframe_wavelengths"] * (1+out_obs["redshift"]), bsed["spectrum"], color="firebrick", label="MAP sample")
    ax.plot(pwave, bsed["photometry"], linestyle="", marker="s", markersize=10, mec="orange", mew=3, mfc="none")
    ax = saxes[0]
    chi = (out_obs["maggies"] - bsed["photometry"]) / out_obs["maggies_unc"]
    ax.plot(pwave, chi, linestyle="", marker="o", color="k")
    ax.axhline(0, color="k", linestyle=":")
    ax.set_ylim(-2, 2)
    ax.set_ylabel(r"$\chi_{\rm best}$")
    plt.show()


def getStellarMass(objID):
    if False:
        # Directly fitting the surviving mass
        hfile = f"prospector/fits_data/{objects_names[objID]}_survMass_prospector_SED.h5"
        out, out_obs, out_model = reader.results_from(hfile)
        mass = out["chain"][:,0]
    else:
        # Taking the fitted total formed mass and then multiplying it by the surviving fraction approximation
        hfile = f"prospector/fits_data/{objects_names[objID]}_prospector_SED.h5"
        out, out_obs, out_model = reader.results_from(hfile)
        mass = out["chain"][:,0]
        mfrac = np.float64(open(f"prospector/data/{objects_names[objID]}_LEGACY_mfrac.txt","r").read())
        mass = mass*mfrac


    return (np.quantile(mass, 0.5), (np.quantile(mass, 0.5)-np.quantile(mass, 0.16)), (np.quantile(mass, 0.84)-np.quantile(mass, 0.5)))


def makeAstropyTableFromDictionnary(dict):
    if dict == None:
        return None
    names = np.array([key for key in dict.keys()])
    data = np.array([val for val in dict.values()])
    dtypes = [type(val) for val in dict.values()]
    return astropy.table.table.Table(data=data, names=names, dtype=dtypes)




if __name__ == "__main__":

    magnitudes_dicts = []
    for i in range(len(objects_names)):
        magnitudes_dicts.append({"name":objects_names[i]})
        for band in "griz":
            try:
                magnitudes_dicts[i][f"cModelMag_{band}"] = QPE_unreddenedMagnitudes[i][band]
                magnitudes_dicts[i][f"cModelMagErr_{band}"] = 0
            except:
                pass

    for i in range(len(magnitudes_dicts)):
        magnitudes_dicts[i] = makeAstropyTableFromDictionnary(magnitudes_dicts[i])

    if input("Fit objects? [y/n]") == "y":
        objID = int(input(f"Input object ID you want to fit [0-{len(objects_names)-1}]:\n"))
        fit_SED(objID, bands="griz", redshift=QPE_redshifts[objID], magnitudes_dict=magnitudes_dicts[objID])

    if input("Read object? [y/n]") == "y":
        objID = int(input(f"Input object ID you want to read [0-{len(objects)-1}]:\n"))
        read_SED(objID)

else:
    QPE_stellar_masses_desiProspector = []
    for i in range(len(objects_names)):
        QPE_stellar_masses_desiProspector.append(getStellarMass(i))