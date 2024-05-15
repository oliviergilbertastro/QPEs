#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
picklename = 'ra19.7861250_dec-34.1916944_AGN.pkl'
fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galaxyFitting.py.
#fitting_run_result.plot_final_qso_fit()
print(fitting_run_result.final_result_galaxy)
import numpy as np

from galight_modif.tools.plot_tools import total_compare
data = fitting_run_result.fitting_specify_class.kwargs_data['image_data']
noise = fitting_run_result.fitting_specify_class.kwargs_data['noise_map']
galaxy_list = fitting_run_result.image_host_list
ps_list = fitting_run_result.image_ps_list
galaxy_total_image = np.zeros_like(galaxy_list[0])
for i in range(len(galaxy_list)):
    galaxy_total_image = galaxy_total_image+galaxy_list[i]
for i in range(len(ps_list)):
    galaxy_total_image = galaxy_total_image+ps_list[i]
model = galaxy_total_image
norm_residual = (data - model)/noise
flux_list_2d = [data, model, norm_residual]
label_list_2d = ['data', 'model', 'normalized residual']
flux_list_1d = [data, model, ps_list[0], galaxy_list[0], -model]
label_list_1d = ['data', 'model', 'AGN', 'galaxy']
total_compare(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d, deltaPix = fitting_run_result.fitting_specify_class.deltaPix,
                      zp=fitting_run_result.zp, if_annuli=False, arrows= False, show_plot = True,
                      target_ID = 'target_ID', sum_rest = True)

fitting_run_result.fitting_specify_class.plot_fitting_sets()
obj_id = int(input('Which obj to measure using statmorph?\n'))
morph = fitting_run_result.cal_statmorph(obj_id=obj_id, segm=fitting_run_result.fitting_specify_class.segm_deblend , if_plot = True)


from statmorph.utils.image_diagnostics import make_figure
import matplotlib.pyplot as plt
fig = make_figure(morph)
plt.show()
print('xc_asymmetry =', morph.xc_asymmetry)
print('yc_asymmetry =', morph.yc_asymmetry)
print('ellipticity_asymmetry =', morph.ellipticity_asymmetry)
print('elongation_asymmetry =', morph.elongation_asymmetry)
print('orientation_asymmetry =', morph.orientation_asymmetry)
print('C =', morph.concentration)
print('A =', morph.asymmetry)
print('S =', morph.smoothness)