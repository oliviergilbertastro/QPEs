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


#PREPARE THE FITTING


fit_sepc = FittingSpecify(data_process)


#For the extended source model, we can use any of the following:
'''
    "GAUSSIAN",
    "GAUSSIAN_ELLIPSE",
    "ELLIPSOID",
    "MULTI_GAUSSIAN",
    "MULTI_GAUSSIAN_ELLIPSE",
    "SERSIC",
    "SERSIC_ELLIPSE",
    "SERSIC_ELLIPSE_Q_PHI",
    "CORE_SERSIC",
    "SHAPELETS",
    "SHAPELETS_POLAR",
    "SHAPELETS_POLAR_EXP",
    "SHAPELETS_ELLIPSE",
    "HERNQUIST",
    "HERNQUIST_ELLIPSE",
    "PJAFFE",
    "PJAFFE_ELLIPSE",
    "UNIFORM",
    "POWER_LAW",
    "NIE",
    "CHAMELEON",
    "DOUBLE_CHAMELEON",
    "TRIPLE_CHAMELEON",
    "INTERPOL",
    "SLIT_STARLETS",
    "SLIT_STARLETS_GEN2",
    "LINEAR",
    "LINEAR_ELLIPSE",
'''

source_params = [
    [{'R_sersic': 2, 'n_sersic': 1, 'e1': 0.020985783580179947, 'e2': 0.020877686531223398, 'center_x': -0.03150875524892178, 'center_y': 0.006710267576362281}],   #init
    [{'n_sersic': 2, 'R_sersic': 0.01469085455157253, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.07924000173807065, 'center_y': 0.07924000173807065}],                     #sigma
    [{'amp': 1, 'e1': 0.4999904134855958, 'e2': 0.058201039792618724, 'center_x': -0.01088466290826861, 'center_y': 0.009910087886402434}],                                                                                                                                           #fixed
    [{'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.001981000043451766, 'n_sersic': 0.2, 'center_x': -0.42770876393927504, 'center_y': -0.389489741113991}],                #lower bound
    [{'e1': 0.5, 'e2': 0.5, 'R_sersic': 305.2036281827358795, 'n_sersic': 9.0, 'center_x': 0.3646912534414315, 'center_y': 0.40291027626671555}]                      #upper bound
          ]

#Prepare the fitting sequence, keywords see notes above.
fit_sepc.prepare_fitting_seq(point_source_num = 1, fix_n_list= None, fix_center_list = None, 
                            extend_source_model=["SERSIC_ELLIPSE"], source_params = None, ps_params = None)

#Plot the initial settings for fittings. 
print("plot_fitting_sets")
fit_sepc.plot_fitting_sets()

#Build up and to pass to the next step.
print("build_fitting_seq")
fit_sepc.build_fitting_seq()




#Setting the fitting method and run.


#Pass fit_sepc to FittingProcess,
# savename: The name of the saved files.    
fit_run = FittingProcess(fit_sepc, savename = 'ra19.7861250_dec-34.1916944', fitting_level='deep') 
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