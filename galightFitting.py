import numpy as np
import astropy.io.fits as pyfits
import galight.tools.astro_tools as astro_tools
from galight_modif.data_process import DataProcess
from galight_modif.fitting_specify import FittingSpecify
from galight_modif.fitting_process import FittingProcess
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

SUBTRACT_NOISE = True

def galight_fit(ra_dec, img_path, oow_path, type="AGN", median_noise=0, PSF_pos_list=None, band="i", nsigma=15, radius=60):
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

    if len(PSF_pos_list) == 0:
        PSF_pos_list = None

    img = pyfits.open(img_path)
    wht_img = pyfits.open(oow_path)

    #Showing the image in the selected band
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, sharex=ax1, sharey=ax1)
    ax1.imshow(copy.copy(wht_img[1].data).T, norm=LogNorm())
    ax2.imshow(copy.copy(img[1].data).T, norm=LogNorm())
    ax1.invert_xaxis()
    plt.show()

    fov_image = img[1].data
    header = img[1].header
    exp =  astro_tools.read_fits_exp(img[0].header)  #Read the exposure time 
    wht = wht_img[1].data
    mean_wht = exp * (0.262)**2  #The drizzle information is used to derive the mean WHT value.
    exp_map = exp * wht/mean_wht  #Derive the exposure time map for each pixel

    #data_process = DataProcess(fov_image = fov_image, target_pos = [1432., 966.], pos_type = 'pixel', header = header,
    #                        rm_bkglight = False, exptime = exp_map, if_plot=True, zp = 22.5)  #zp use 27.0 for convinence.
    data_process = DataProcess(fov_image = fov_image, target_pos = [ra_dec[0], ra_dec[1]], pos_type = 'wcs', header = header,
                            rm_bkglight = False, exptime = exp_map, if_plot=False, zp = 22.5)  #zp use 27.0 for convinence.

    data_process.generate_target_materials(radius=60, create_mask = False, nsigma=15,
                                        exp_sz= 1, npixels = 5, if_plot=False)
    
    if SUBTRACT_NOISE:
        #To get a good fit, we remove the median noise from the initial image and restart the procedure
        fov_image = fov_image-np.ones(fov_image.shape)*data_process.bkg_mid
        data_process = DataProcess(fov_image = fov_image, target_pos = [ra_dec[0], ra_dec[1]], pos_type = 'wcs', header = header,
                            rm_bkglight = False, exptime = exp_map, if_plot=True, zp = 22.5)
        data_process.generate_target_materials(radius=radius, create_mask = True, nsigma=nsigma,
                                        exp_sz= 1, npixels = 5, if_plot=True)

    print('---------------DATA PROCESS PARAMETERS-------------')
    print('target_pos:', data_process.target_pos)
    print('zero point:', data_process.zp) #zp is in the AB system and should be 22.5: https://www.legacysurvey.org/svtips/
    print('kwargs: ', data_process.arguments)
    print('---------------------------------------------------')


    #data_process.find_PSF(radius = 30, user_option = True)  #Try this line out! 
    data_process.find_PSF(radius = 30, PSF_pos_list = PSF_pos_list, pos_type="wcs", user_option=True)

    #Plot the FOV image and label the position of the target and the PSF
    data_process.plot_overview(label = 'Example', target_label = None)

    # Compare the 1D profile of all the components.
    data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False,
                    prf_name_list = (['target'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_list))]) )

    #Select which PSF id you want to use to make the fitting.
    data_process.psf_id_for_fitting = int(input('Use which PSF? Input a number.\n'))

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
        def condition_bulgedisk(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction, kwargs_tracer_source):
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




    #PREPARE THE FITTING
    fit_sepc = FittingSpecify(data_process)

    #Prepare the fitting sequence, keywords see notes above.
    fit_sepc.prepare_fitting_seq(point_source_num = number_of_ps, fix_n_list= None, condition=condition_bulgedisk)

    #Plot the initial settings for fittings.
    fit_sepc.plot_fitting_sets()

    #Build up and to pass to the next step.
    fit_sepc.build_fitting_seq()


    #Setting the fitting method and run.

    #Pass fit_sepc to FittingProcess,
    fit_run = FittingProcess(fit_sepc, savename = f'ra{str(ra_dec[0])}_dec{str(ra_dec[1])}_{type}_{band}', fitting_level='deep') 

    #Setting the fitting approach and Run:
    fit_run.run(algorithm_list = ['PSO', 'MCMC'], setting_list = None)
    #fit_run.run(algorithm_list = ['PSO', 'MCMC'], setting_list = [None, {'n_burn': 200, 'n_run': 1000, 'walkerRatio': 10, 'sigma_scale': .1}])

    # Plot all the fitting results:
    fit_run.plot_all(target_ID=f'{str(ra_dec[0])+str(ra_dec[1])}-{band}')

    #Save the fitting class as pickle format:
    fit_run.dump_result()