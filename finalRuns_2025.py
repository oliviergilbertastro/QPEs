from multiprocessing import Process

from download_data import TDE_coords, TDE_names, objects, objects_names, TDE_types


#ignore warnings:
import shutup
shutup.please()



import numpy as np
import astropy.io.fits as pyfits
import galight.tools.astro_tools as astro_tools
from galight.data_process import DataProcess
from galight.fitting_specify import FittingSpecify
from galight.fitting_process import FittingProcess
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mpl_interactions.ipyplot as iplt
from matplotlib.widgets import Slider
from download_data import objects
import copy
import lenstronomy.Util.param_util as param_util
from matplotlib.colors import LogNorm
import pickle

SUBTRACT_NOISE = False

def galight_fit_short(ra_dec, img_path, oow_path=None, exp_path=None, psf_path=None, type="AGN", pixel_scale=0.262, PSF_pos_list=None, band="i", nsigma=15, radius=60, exp_sz_multiplier=1, npixels=5, survey="DESI", savename=None, threshold=5, fitting_level="deep", fixed_n_list=None):
    
    
    if type in ["AGN", "agn", "Agn"]:
        type = "AGN"
        number_of_ps = 1
        bulge = False
    elif type in ["Bulge", "BULGE", "bulge"]:
        type = "Bulge"
        number_of_ps = 0
        bulge = True
    elif type in ["Bulge+AGN", "Bulge+Agn", "BULGE+AGN", "bulge+agn", "AGN+Bulge", "Agn+Bulge", "AGN+BULGE", "agn+bulge"]:
        type = "Bulge+AGN"
        number_of_ps = 1
        bulge = True
    elif type in ["None", "Sersic", "none", "sersic", "single", ""]:
        type = "None"
        number_of_ps = 0
        bulge = False
    elif type in ["Bulge_fixed", "BULGE_FIXED", "bulge_fixed"]:
        type = "Bulge_fixed"
        number_of_ps = 0
        bulge = True
    else:
        raise ValueError(f"type {type} is not a supported fitting type")
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
    if PSF_pos_list == None:
        PSF_pos_list = None
    elif len(PSF_pos_list) == 0:
        PSF_pos_list = None

    img = pyfits.open(img_path)
    if oow_path != None:
        wht_img = pyfits.open(oow_path)
    fov_noise_map = None
    zp = 22.5
    if survey == "PANSTARRS":
        fov_image = img[0].data
        header = img[0].header
        exp =  astro_tools.read_fits_exp(img[0].header)  #Read the exposure time 
        if exp_path != None:
            exp_img = pyfits.open(exp_path)
            exp_map = exp_img[0].data
        band_index = ["g","r","i","z"].index(band)
        zp = [24.41,24.68,24.56,24.22][band_index]

    elif survey == "COADDED_DESI":
        band_index = ["g","r","i","z"].index(band)
        fov_image = (img[0].data)#[:,:,band_index]
        header = img[0].header
        exp =  1  #Read the exposure time
        exp_map = exp
        wht = wht_img[1].data
        mean_wht = exp * (pixel_scale)**2  #The drizzle information is used to derive the mean WHT value.
        exp_map = exp * wht/mean_wht  #Derive the exposure time map for each pixel
        fov_noise_map = 1/np.sqrt(wht)


    elif survey == "PS":
        header = img[1].header
        #header.set('RADESYS','FK5') #Correction so target actually shows up in the image
        fov_image = (img[1].data)
        exp =  1  #Read the exposure time
        exp_map = pyfits.open(exp_path)[1].data
        band_index = ["g","r","i","z"].index(band)
        zp = [24.41,24.68,24.56,24.22][band_index]

    plt.imshow(np.log10(fov_image))
    plt.show()
    #data_process = DataProcess(fov_image = fov_image, target_pos = [1432., 966.], pos_type = 'pixel', header = header,
    #                        rm_bkglight = False, exptime = exp_map, if_plot=True, zp = 22.5)  #zp use 27.0 for convinence.
    data_process = DataProcess(fov_image = fov_image, target_pos = [ra_dec[0], ra_dec[1]], pos_type = 'wcs', header = header,
                            rm_bkglight = True, exptime = exp_map, if_plot=True, zp = zp, fov_noise_map=fov_noise_map)  #zp use 27.0 for convinence.

    data_process.generate_target_materials(radius=radius, create_mask = True, nsigma=nsigma,
                                        exp_sz= exp_sz_multiplier, npixels = npixels, if_plot=True, show_materials=False)

    #Create a dictionnary of all informations that could be useful for us
    coolinfos = {}
    coolinfos["cutout_radius"] = radius
    coolinfos["survey"] = survey
    coolinfos["pix_scale"] = pixel_scale
    print('---------------DATA PROCESS PARAMETERS-------------')
    print('target_pos:', data_process.target_pos)
    print('zero point:', data_process.zp) #zp is in the AB system and should be 22.5: https://www.legacysurvey.org/svtips/
    print('---------------------------------------------------')

    if psf_path == None:
        if PSF_pos_list == None and survey == "COADDED_DESI":
            data_process.find_PSF(radius = 30, user_option = True, threshold=threshold)  #Try this line out!
        else:
            data_process.find_PSF(radius = 30, PSF_pos_list = PSF_pos_list, pos_type="wcs", user_option=True, threshold=threshold)
        #Select which PSF id you want to use to make the fitting.
        data_process.psf_id_for_fitting = int(input('Use which PSF? Input a number.\n'))
    else:
        psf_img = pyfits.open(psf_path)[0].data
        fwhm = data_process.use_custom_psf(psf_img, if_plot=False)
        print(f'FWHM = {np.around(fwhm*pixel_scale, decimals=3)}"')

    #Check if all the materials is given, if so to pass to the next step.
    data_process.checkout()

    if bulge:
        apertures = copy.deepcopy(data_process.apertures)
        comp_id = 0
        add_aperture0 = copy.deepcopy(apertures[comp_id])
        #This setting assigns comp0 as 'bulge' and comp1 as 'disk'
        add_aperture0.a, add_aperture0.b = add_aperture0.a/2, add_aperture0.b/2
        apertures = apertures[:comp_id] + [add_aperture0] + apertures[comp_id:]
        data_process.apertures = apertures #Pass apertures to the data_process

        #Adding a prior so that 1)the size of the bulge is within a range to the disk size. 2) disk have more ellipticity.
        def condition_bulgedisk(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction):
            logL = 0
            #note that the Comp[0] is the bulge and the Comp[1] is the disk.
            phi0, q0 = param_util.ellipticity2phi_q(kwargs_lens_light[0]['e1'], kwargs_lens_light[0]['e2'])
            phi1, q1 = param_util.ellipticity2phi_q(kwargs_lens_light[1]['e1'], kwargs_lens_light[1]['e2'])
            cond_0 = (kwargs_lens_light[0]['R_sersic'] > kwargs_lens_light[1]['R_sersic'] * 0.9)
            cond_1 = (kwargs_lens_light[0]['R_sersic'] < kwargs_lens_light[1]['R_sersic']*0.15)
            cond_2 = (q0 < q1)
            if cond_0 or cond_1 or cond_2:
                logL -= 10**15
            return logL
    else:
        condition_bulgedisk = None

    if False:
        def condition_bulgedisk(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction, kwargs_tracer_source):
                logL = 0
                cond_0 = (kwargs_lens_light[0]['n_sersic'] < 2)
                if cond_0:
                    logL -= 10**15
                return logL


    #PREPARE THE FITTING
    fit_sepc = FittingSpecify(data_process)
    

    #Prepare the fitting sequence, keywords see notes above.
    fit_sepc.prepare_fitting_seq(point_source_num = number_of_ps,
                                fix_Re_list=None,
                                fix_n_list=fixed_n_list, #To fix the Sérsic index at 2.09: [[0,2.09]]
                                condition=condition_bulgedisk)


    #Build up and to pass to the next step.
    fit_sepc.build_fitting_seq()


    #Setting the fitting method and run.

    #Pass fit_sepc to FittingProcess,
    if savename == None:
        savename = f'ra{str(ra_dec[0])}_dec{str(ra_dec[1])}_{type}_{band}'
    fit_run = FittingProcess(fit_sepc, savename = savename, fitting_level=fitting_level) 

    #Setting the fitting approach and Run:
    fit_run.run(algorithm_list = ['PSO', 'MCMC'], setting_list = None)
    #fit_run.run(algorithm_list = ['PSO', 'MCMC'], setting_list = [None, {'n_burn': 200, 'n_run': 1000, 'walkerRatio': 10, 'sigma_scale': .1}])
    fit_run.mcmc_result_range()
    if "I dont want to plot these as Im leaving the lab" == "y":
        # Plot all the fitting results:
        if type == "None" or type == "Bulge":
            fit_run.plot_final_galaxy_fit(target_ID=f'{str(ra_dec[0])+str(ra_dec[1])}-{band}')
        else:
            fit_run.plot_final_qso_fit(target_ID=f'{str(ra_dec[0])+str(ra_dec[1])}-{band}')
    fit_run.coolinfos = coolinfos
    #Save the fitting class as pickle format:
    fit_run.dump_result()



def fit_single_object(qpe_oder_tde="QPE", objID=0, bands="r", types=["None"], psf_band="r", fixed_n_list=None):
    if qpe_oder_tde == "TDE":
        coords = TDE_coords
        path_section = "tde"
        names = TDE_names
    else:
        coords = objects
        path_section = "qpe"
        names = objects_names

    for current_type in types:
        for band in bands:
            picklename = f"{names[objID]}_r-band_None_DESI_PSF.pkl"
            print(f"\x1b[34m{names[objID]}\x1b[0m")

            fixed_n_list = None

            img_path = f"data/images/{path_section}{objID}_{band}.fits"
            oow_path = f"data/images/{path_section}{objID}_{band}.fits"
            #Use the co-add PSF model from the survey
            psf_path = f"data/images/{path_section}{objID}_{psf_band}_PSF.fits"
            galight_fit_short(
                coords[objID],
                img_path,
                oow_path,
                None,
                psf_path,
                current_type,
                0.262,
                None,
                band,
                15,
                60,
                1,
                5,
                "COADDED_DESI",
                f"{names[objID]}_{band}-band_{current_type}_DESI_PSF_FINAL_2025",
                5,
                "mega_deep",
                fixed_n_list,
                )
            print("\x1b[32mThis one did work\x1b[0m")
            
            try:
                pass
            except Exception as e:
                print(e)
                print("\x1b[31mThis one didn't work\x1b[0m")
                #psf_bands = "griz"
                #fit_single_object(qpe_oder_tde=qpe_oder_tde, objID=objID, bands=bands, types=types, psf_band=psf_bands[(psf_bands.index(psf_band)+1)%len(psf_bands)])


def galight_fit_2_comp(id:int, save:bool=True, if_plot:bool=False, verbose:bool=False, radius:int=50, custom_mask:list=None, exp_sz:float= 1.5) -> None:
    """
    Fit an image with Point-Source + Sersic

    Parameter
        --------
        id : int, cweb-id for BLAGN
        save : bool, flag to choose whether to save plots and pickled results
        if_plot : bool, flag to choose whether to show the plots while running
        verbose : bool, flag to choose whether to print all information to terminal while running
        radius : int, radius of cutout image
        custom_mask : list, optional list that can be provided to keep more than 1 region in the fitting (only useful if bright contaminant nearby)
        exp_sz : float, multiplier of the source regions, which often creates a better mask
    """
    for filter in ["r"]:
        try:
            
            fov_image, exp, wht, sig_map, header, full_tile = get_data_for_fitting(id, filter, cutout_size=(10,10))
            pix_per_arcsec = 1.0 / WCS(header).proj_plane_pixel_scales()[0].to('deg').value / 3600.0  # pixels per arcsec
            pixel_scale = 1/pix_per_arcsec # arcsec/pixel
            mean_wht = exp * (pixel_scale)**2  #The drizzle information is used to derive the mean WHT value.
            exp_map = exp * wht/mean_wht  #Derive the exposure time map for each pixel
            fov_noise_map = sig_map#1/np.sqrt(wht)
            fov_image_data = np.array(fov_image, dtype=np.float32) # Needs to be converted to array
            fov_image_data = np.nan_to_num(fov_image_data, nan=0) # Remove NaNs from fov_image
            fov_noise_map = np.nan_to_num(fov_noise_map, nan=np.nanmedian(fov_noise_map)) # Remove NaNs from fov_noise_map
            data_process = DataProcess(fov_image = fov_image_data, target_pos = coords, pos_type = 'wcs', header = header,
                                    rm_bkglight = True, exptime = exp_map, if_plot=if_plot, zp = zp, fov_noise_map=fov_noise_map, verbose=verbose)
            data_process.generate_target_materials(radius=radius, psf_only=False, create_mask = True, custom_mask=custom_mask, nsigma=1, exp_sz= exp_sz, npixels = 10, if_plot=if_plot, show_materials=if_plot)
            if verbose:
                print('---------------DATA PROCESS PARAMETERS-------------')
                print('target_pos:', data_process.target_pos)
                print('zero point:', data_process.zp) #zp is in the AB system
                print('---------------------------------------------------')
            psf_img = load_Jackie_psf(filter)
            fwhm = data_process.use_custom_psf(psf_img, if_plot=if_plot) # Custom function to use our own PSF that we constructed instead of building it with Galight
            if verbose: print(f'FWHM = {np.around(fwhm*pixel_scale, decimals=3)}"')
            #Check if all the materials is given, if so to pass to the next step.
            data_process.checkout()

            # Add condition that flux of PS must be greater than flux of host:
            def condition_ps_brighter_than_sersic(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps,kwargs_special, kwargs_extinction):
                """
                Forces the point-source's flux to be greater than the Sérsic's flux by applying a likelihood penalty if violated.
                """
                logL = 0
                # --- Extract fluxes ---
                ps_flux = kwargs_ps[0].get('flux', kwargs_ps[0].get('point_amp'))
                # Sérsic flux (sum over all lens-light components)
                sersic_flux = 0
                for comp in kwargs_lens_light:
                    sersic_flux += comp['flux']
                if not (ps_flux > sersic_flux):
                    logL -= 1e15
                return logL

            #PREPARE THE FITTING
            fit_sepc = FittingSpecify(data_process)
            #Prepare the fitting sequence, keywords see notes above.
            fit_sepc.prepare_fitting_seq(point_source_num = 1, # Only one Point-source
                                        fix_Re_list=None,
                                        fix_n_list=None,
                                        extend_source_model=None, # We don't want to fit additional Sérsic component(s) (for the moment, at least...)
                                        condition=condition_ps_brighter_than_sersic,
                                        ps_pix_center_list=[[0, 0],],
                                        )
            #Build up and to pass to the next step.
            fit_sepc.build_fitting_seq()
            if if_plot:
                fit_sepc.plot_fitting_sets()
            #Pass fit_sepc to FittingProcess,
            fit_run = FittingProcess(fit_sepc, savename = f"/quasar/ogilbert/dataproducts/galight_fits/{id}_{filter}_2comp", fitting_level="norm") 
            #Setting the fitting approach and Run:
            fit_run.run(algorithm_list = ['PSO','MCMC'], setting_list = None)
            # Plot all the fitting results:
            if if_plot:
                fit_run.plot_all(target_ID=f'{id}_{filter}', show_plot=if_plot, save_plot=False)
            fit_run.plot_final_qso_fit(target_ID=f'{id}_{filter}', save_plot=save, savename=f"/quasar/ogilbert/figures/galight_fits_2comp/{id}_{filter}")
            #Save the fitting class as pickle format:
            if save:
                fit_run.dump_result()
            plt.close('all')
        except:
            traceback.print_exc()
            print(f"Filter {filter} failed for ID {id}.")
    pass

