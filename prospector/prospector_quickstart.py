'''
Get SDSS photometry data of objects, fit galaxy models to it using nested sampling to obtain the total stellar mass in the galaxy
'''


import os
import sys
import copy
import importlib.util
import sys
spec = importlib.util.spec_from_file_location("download_data", "/Users/oliviergilbert/Desktop/QPEs/QPEs/download_data.py")
download_data = importlib.util.module_from_spec(spec)
sys.modules["download_data"] = download_data
spec.loader.exec_module(download_data)
spec = importlib.util.spec_from_file_location("paper_data", "/Users/oliviergilbert/Desktop/QPEs/QPEs/paper_data.py")
paper_data = importlib.util.module_from_spec(spec)
sys.modules["paper_data"] = paper_data
spec.loader.exec_module(paper_data)


#Little code to ensure we are in the user directory running this code so the env variable works correctly  
home = os.getcwd()
if home != "/Users/oliviergilbert":
    os.system("cd ../../../ \nsource .bash_profile \npython3 /Users/oliviergilbert/Desktop/QPEs/QPEs/prospector/prospector_quickstart.py")
    sys.exit()

#Need to download fsps from https://github.com/cconroy20/fsps and put it locally, change the path in .bash_profile to point to it

import astropy.table
import fsps
import dynesty
import sedpy
import h5py, astropy
import numpy as np
import astroquery

def fit_SED(pos, bands="ugriz", redshift=0, magnitudes_dict=None):
    """
    pos = (ra,dec)
    """
    ra, dec = pos
    print(ra, dec)
    from astroquery.sdss import SDSS
    from astropy.coordinates import SkyCoord
    mcol = [f"cModelMag_{b}" for b in bands]
    ecol = [f"cModelMagErr_{b}" for b in bands]
    data_release = int(input("Which data release do you want to use? [1-18]"))
    cat = SDSS.query_crossid(SkyCoord(ra=ra, dec=dec, unit="deg"),
                            data_release=data_release,
                            photoobj_fields=mcol + ecol + ["specObjID"])
    #This code is only if you don't have a redshift for the source and find the redshift in the spectra
    #shdus = SDSS.get_spectra(plate=2101, mjd=53858, fiberID=220)[0]
    #assert int(shdus[2].data["SpecObjID"][0]) == cat[0]["specObjID"]
    #redshift = shdus[2].data[0]["z"]

    from sedpy.observate import load_filters
    from prospect.utils.obsutils import fix_obs

    filters = load_filters([f"sdss_{b}0" for b in bands])
    working_filters = ""
    for b in bands:
        try:
            cat[0][f"cModelMag_{b}"]
            working_filters += b
        except:
            pass
    print("working_filters:", working_filters)
    print("--------------------------------------------------------")
    print(type(cat))
    print(cat)
    print("--------------------------------------------------------")
    if magnitudes_dict != None:
        cat = copy.copy(magnitudes_dict)

    maggies = np.array([10**(-0.4 * cat[0][f"cModelMag_{b}"]) for b in bands])
    magerr = np.array([cat[0][f"cModelMagErr_{b}"] for b in bands])
    magerr = np.hypot(magerr, 0.05)

    obs = dict(wavelength=None, spectrum=None, unc=None, redshift=redshift,
            maggies=maggies, maggies_unc=magerr * maggies / 1.086, filters=filters)
    obs = fix_obs(obs)
    #print(obs)


    from prospect.models.templates import TemplateLibrary
    from prospect.models import SpecModel
    model_params = TemplateLibrary["parametric_sfh"]
    model_params.update(TemplateLibrary["nebular"])
    model_params["zred"]["init"] = obs["redshift"]

    model = SpecModel(model_params)
    assert len(model.free_params) == 5
    #print(model)

    noise_model = (None, None)

    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=1)
    #print(sps.ssp.libraries)

    current_parameters = ",".join([f"{p}={v}" for p, v in zip(model.free_params, model.theta)])
    print(current_parameters)
    spec, phot, mfrac = model.predict(model.theta, obs=obs, sps=sps)
    print(phot / obs["maggies"])

    from prospect.fitting import lnprobfn, fit_model
    fitting_kwargs = dict(nlive_init=400, nested_method="rwalk", nested_target_n_effective=1000, nested_dlogz_init=0.05)
    output = fit_model(obs, model, sps, optimize=False, dynesty=True, lnprobfn=lnprobfn, noise=noise_model, **fitting_kwargs)
    result, duration = output["sampling"]
    from prospect.io import write_results as writer
    hfile = f"/Users/oliviergilbert/Desktop/QPEs/QPEs/prospector/fits_data/quickstart_dynesty_mcmc_{ra}_{dec}_{data_release}.h5"
    writer.write_hdf5(hfile, {}, model, obs,
                    output["sampling"][0], None,
                    sps=sps,
                    tsample=output["sampling"][1],
                    toptimize=0.0)

def read_SED(pos, data_release=16):
    ra, dec = pos
    hfile = f"/Users/oliviergilbert/Desktop/QPEs/QPEs/prospector/fits_data/quickstart_dynesty_mcmc_{ra}_{dec}_{data_release}.h5"
    from prospect.io import read_results as reader
    out, out_obs, out_model = reader.results_from(hfile)
    for k in out.keys():
        print(k, ":", out[k])

    import matplotlib.pyplot as plt
    #Plot the posterior's corner plot
    from prospect.plotting import corner
    nsamples, ndim = out["chain"].shape
    cfig, axes = plt.subplots(ndim, ndim, figsize=(10,9))
    axes = corner.allcorner(out["chain"].T, out["theta_labels"], axes, weights=out["weights"], color="royalblue", show_titles=True)

    from prospect.plotting.utils import best_sample
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



if __name__ == "__main__":
    objects = download_data.objects
    QPE_redshifts = paper_data.QPE_redshifts
    bands_for_each_obj = ["ugriz", 
                          "ugriz", 
                          "ugriz", 
                          "ugriz", 
                          "ugriz", 
                          "ugriz", 
                          "ugriz", 
                          "ugriz", 
                          "ugriz",
                          ]
    magnitudes_dicts = [
        astropy.table.table.Table(np.array(["obj_0",1237667444048855291,16.92195,15.41549,14.9153,14.63878,14.35094,0.009043842,0.002720573,0.002580761,0.00274345,0.004231831,2523306138166913024,"GALAXY"]), names=np.array(["name","objID","cModelMag_u","cModelMag_g","cModelMag_r","cModelMag_i","cModelMag_z","cModelMagErr_u","cModelMagErr_g","cModelMagErr_r","cModelMagErr_i","cModelMagErr_z","specObjID","type"]))
    ]
    #print(magnitudes_dicts[0]["name"])
    if input("Fit objects? [y/n]") == "y":
        #fit_thingamabob((204.46376, 35.79883), redshift=0.07260209)
        objID = int(input(f"Input object ID you want to fit [0-{len(objects)-1}]:\n"))
        fit_SED(objects[objID], bands=bands_for_each_obj[objID], redshift=QPE_redshifts[objID])

    if input("Read object? [y/n]") == "y":
        #read_thingamabob((204.46376, 35.79883))
        objID = int(input(f"Input object ID you want to read [0-{len(objects)-1}]:\n"))
        data_release = input("Which data release?")
        data_release = 16 if data_release == "" else int(data_release)
        read_SED(objects[objID], data_release=data_release)