import os
import sys

#Little code to ensure we are in the user directory running this code so the env variable works correctly  
home = os.getcwd()
if home != "/Users/oliviergilbert":
    os.system("cd ../../../ \nsource .bash_profile \npython3 /Users/oliviergilbert/Desktop/QPEs/QPEs/prospector/prospector_quickstart.py")
    sys.exit()

#Need to download fsps from https://github.com/cconroy20/fsps and put it locally, change the path in .bash_profile to point to it

import fsps
import dynesty
import sedpy
import h5py, astropy
import numpy as np
import astroquery

def fit_thingamabob(pos):
    """
    pos = (ra,dec)
    """
    ra, dec = pos
    from astroquery.sdss import SDSS
    from astropy.coordinates import SkyCoord
    bands = "ugriz"
    mcol = [f"cModelMag_{b}" for b in bands]
    ecol = [f"cModelMagErr_{b}" for b in bands]
    cat = SDSS.query_crossid(SkyCoord(ra=ra, dec=dec, unit="deg"),
                            data_release=16,
                            photoobj_fields=mcol + ecol + ["specObjID"])
    shdus = SDSS.get_spectra(plate=2101, mjd=53858, fiberID=220)[0]
    assert int(shdus[2].data["SpecObjID"][0]) == cat[0]["specObjID"]

    from sedpy.observate import load_filters
    from prospect.utils.obsutils import fix_obs

    filters = load_filters([f"sdss_{b}0" for b in bands])
    maggies = np.array([10**(-0.4 * cat[0][f"cModelMag_{b}"]) for b in bands])
    magerr = np.array([cat[0][f"cModelMagErr_{b}"] for b in bands])
    magerr = np.hypot(magerr, 0.05)

    obs = dict(wavelength=None, spectrum=None, unc=None, redshift=shdus[2].data[0]["z"],
            maggies=maggies, maggies_unc=magerr * maggies / 1.086, filters=filters)
    obs = fix_obs(obs)
    print(obs)


    from prospect.models.templates import TemplateLibrary
    from prospect.models import SpecModel
    model_params = TemplateLibrary["parametric_sfh"]
    model_params.update(TemplateLibrary["nebular"])
    model_params["zred"]["init"] = obs["redshift"]

    model = SpecModel(model_params)
    assert len(model.free_params) == 5
    print(model)

    noise_model = (None, None)

    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=1)
    print(sps.ssp.libraries)

    current_parameters = ",".join([f"{p}={v}" for p, v in zip(model.free_params, model.theta)])
    print(current_parameters)
    spec, phot, mfrac = model.predict(model.theta, obs=obs, sps=sps)
    print(phot / obs["maggies"])

    from prospect.fitting import lnprobfn, fit_model
    fitting_kwargs = dict(nlive_init=400, nested_method="rwalk", nested_target_n_effective=1000, nested_dlogz_init=0.05)
    output = fit_model(obs, model, sps, optimize=False, dynesty=True, lnprobfn=lnprobfn, noise=noise_model, **fitting_kwargs)
    result, duration = output["sampling"]
    print("This is still running... line74")
    from prospect.io import write_results as writer
    hfile = f"/Users/oliviergilbert/Desktop/QPEs/QPEs/prospector/fits_data/quickstart_dynesty_mcmc_{ra}_{dec}.h5"
    writer.write_hdf5(hfile, {}, model, obs,
                    output["sampling"][0], None,
                    sps=sps,
                    tsample=output["sampling"][1],
                    toptimize=0.0)

def read_thingamabob(pos):
    ra, dec = pos
    hfile = f"/Users/oliviergilbert/Desktop/QPEs/QPEs/prospector/fits_data/quickstart_dynesty_mcmc_{ra}_{dec}.h5"
    from prospect.io import read_results as reader
    out, out_obs, out_model = reader.results_from(hfile)

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
    if input("Fit objects? [y/n]") == "y":
        fit_thingamabob((204.46376, 35.79883))

    if input("Read object? [y/n]") == "y":
        read_thingamabob((204.46376, 35.79883))