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

#THIS IS FOR OBJECT0 IN THE i-BAND ONLY
img_path = "/Users/oliviergilbert/Downloads/image-decam-498255-N24-i.fits.gz"
oow_path = "/Users/oliviergilbert/Downloads/iv-decam-498255-N24-i.fits.gz"

img = pyfits.open(img_path)
wht_img = pyfits.open(oow_path)

if True:
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
    ax1.imshow(wht_img[1].data)
    ax2.imshow(img[1].data)
    plt.show()

fov_image = img[1].data-np.ones(img[1].data.shape)*863.3368 #trying to remove the noise
header = img[0].header
exp =  astro_tools.read_fits_exp(header)  #Read the exposure time 
wht = wht_img[1].data
mean_wht = exp * (0.262)**2  #The drizzle information is used to derive the mean WHT value.
exp_map = exp * wht/mean_wht  #Derive the exposure time map for each pixel


data_process = DataProcess(fov_image = fov_image, target_pos = [1432., 966.], pos_type = 'pixel', header = header,
                          rm_bkglight = False, exptime = exp_map, if_plot=True, zp = 22.5)  #zp use 27.0 for convinence.


data_process.generate_target_materials(radius=60, create_mask = True, nsigma=15,
                                      exp_sz= 1, npixels = 5, if_plot=True)

print('---------------DATA PROCESS PARAMETERS-------------')
print('target_pos:', data_process.target_pos)
print('zero point:', data_process.zp) #zp is in the AB system and should be 22.5: https://www.legacysurvey.org/svtips/
print('kwargs: ', data_process.arguments)
print('---------------------------------------------------')


#data_process.find_PSF(radius = 30, user_option = True)  #Try this line out! 
data_process.find_PSF(radius = 30, PSF_pos_list = [[ 1617., 1027.], [1502., 667.], [1150., 1437.]], user_option=True)

#Plot the FOV image and label the position of the target and the PSF
data_process.plot_overview(label = 'Example', target_label = None)

# Compare the 1D profile of all the components.
data_process.profiles_compare(norm_pix = 5, if_annuli=False, y_log = False,
                  prf_name_list = (['target'] + ['PSF{0}'.format(i) for i in range(len(data_process.PSF_list))]) )


#Select which PSF id you want to use to make the fitting.
data_process.psf_id_for_fitting = int(input('Use which PSF? Input a number.\n'))

#Check if all the materials is given, if so to pass to the next step.
data_process.checkout()



#Add bulge

#Modify the fitting of the component (i.e., galaxy) id = 0 into to components (i.e., bulge + disk)
import copy
apertures = copy.deepcopy(data_process.apertures)
comp_id = 0
add_aperture0 = copy.deepcopy(apertures[comp_id])
#This setting assigns comp0 as 'bulge' and comp1 as 'disk'
add_aperture0.a, add_aperture0.b = add_aperture0.a/2, add_aperture0.b/2
apertures = apertures[:comp_id] + [add_aperture0] + apertures[comp_id:]
data_process.apertures = apertures #Pass apertures to the data_process

#Adding a prior so that 1)the size of the bulge is within a range to the disk size. 2) disk have more ellipticity.
import lenstronomy.Util.param_util as param_util
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


#Start to produce the class and params for fitting.
#For more details, see notebook galight_HSC_QSO.ipynb
fit_sepc = FittingSpecify(data_process)
#The 'fix_n_list' will fix Sersic_n as 4 for the comp0, and as 1 for the comp1.
fit_sepc.prepare_fitting_seq(point_source_num = 0, fix_n_list= [[0,1], [1,1.14]],condition=condition_bulgedisk)

#The settings of the parameters is the dict defined by fit_sepc.kwargs_params. One can modify the values herein to change the default setting manually. 

#Plot the initial settings for fittings. 
fit_sepc.build_fitting_seq()
fit_sepc.plot_fitting_sets()




#Setting the fitting method and run.


#Pass fit_sepc to FittingProcess,
# savename: The name of the saved files.    
fit_run = FittingProcess(fit_sepc, savename = 'ra19.7861250_dec-34.1916944_bulge', fitting_level='deep') 
#For the fitting_level, you can also put ['norm', 'deep'] for the later ['PSO', 'MCMC'] corresplingly.

#Setting the fitting approach and Run: 
#     algorithm_list: The fitting approaches that would be used: e.g. ['PSO', 'PSO', 'MCMC']
#     setting_list: The detailed setting for the fitting:
#     -for PSO:
#         input template: {'sigma_scale': 0.8, 'n_particles': 50, 'n_iterations': 50}
#     -for MCMC:
#         input template: {'n_burn': 50, 'n_run': 100, 'walkerRatio': 10, 'sigma_scale': .1}
#     if setting_list = [None, None, None], default values would be given 
fit_run.run(algorithm_list = ['PSO', 'MCMC'], setting_list = None)
#fit_run.run(algorithm_list = ['PSO', 'MCMC'], setting_list = [None, {'n_burn': 200, 'n_run': 1000, 'walkerRatio': 10, 'sigma_scale': .1}])

# Plot all the fitting results, including:
#         run_diag() : The convergence of the chains.
#         model_plot(): The model plot (by lenstronomy)
#         plot_params_corner(): The mcmc corner for all the chains (MCMC should be peformed) 
#         plot_flux_corner(): The flux corner for all the component (MCMC should be peformed)
#         plot_final_qso_fit() or plot_final_galaxy_fit(): Plot the overall plot (data, model, data-ps, resudal, 1D profile)
fit_run.plot_all(target_ID='GSN 069')

#Save the fitting class as pickle format:
#     Note, if you use python3 (or 2), load with python3 (or 2)
fit_run.dump_result()
fitting_run_result = fit_run
from galight_modif.tools.plot_tools import total_compare
data = fitting_run_result.fitting_specify_class.kwargs_data['image_data']
noise = fitting_run_result.fitting_specify_class.kwargs_data['noise_map']
galaxy_list = fitting_run_result.image_host_list
galaxy_total_image = np.zeros_like(galaxy_list[0])
for i in range(len(galaxy_list)):
    galaxy_total_image = galaxy_total_image+galaxy_list[i]
model = galaxy_total_image
norm_residual = (data - model)/noise
flux_list_2d = [data, model, norm_residual]
label_list_2d = ['data', 'model', 'normalized residual']
flux_list_1d = [data, model, galaxy_list[0], galaxy_list[1], -model]
label_list_1d = ['data', 'model', 'bulge', 'disk']
total_compare(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d, deltaPix = fitting_run_result.fitting_specify_class.deltaPix,
                      zp=fitting_run_result.zp, if_annuli=False, arrows= False, show_plot = True,
                      target_ID = 'target_ID', sum_rest = True)