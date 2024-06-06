import sys
import parse
import os
import copy
import numpy as np
import matplotlib.pyplot as plt

parent_dir = parse.parse("{}/prospector", os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(parent_dir)
import download_data, paper_data

def read_SED(pos, data_release=16):
    ra, dec = pos
    hfile = f"prospector/fits_data/quickstart_dynesty_mcmc_{ra}_{dec}_{data_release}.h5"
    from prospect.io import read_results as reader
    out, out_obs, out_model = reader.results_from(hfile)
    for k in out.keys():
        print(k, ":", out[k])

    
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

    #This is mainly for GSN069 and its data comes from https://simbad.u-strasbg.fr/simbad/sim-id?protocol=html&Ident=GSN_069&bibdisplay=none
    magnitudes_dicts = [
        {"name":"GSN 069",
         "objID": "idk",
         "cModelMag_J": 13.579,
         "cModelMag_H": 12.888,
         "cModelMag_Ks": 12.748,
         "cModelMagErr_J": 0.056,
         "cModelMagErr_H": 0.069,
         "cModelMagErr_Ks": 0.123,
         "specObjID": "idk",
         "type": "GALAXY"},
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ]

    if input("Read object? [y/n]") == "y":
        #read_thingamabob((204.46376, 35.79883))
        objID = int(input(f"Input object ID you want to read [0-{len(download_data.objects)-1}]:\n"))
        if magnitudes_dicts[objID] != None:
            data_release = "custom" #Not really a data release, just to read file name
        else:
            data_release = input("Which data release? [default=18]")
            data_release = 18 if data_release == "" else int(data_release)
        read_SED(download_data.objects[objID], data_release=data_release)