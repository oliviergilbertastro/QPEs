#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
picklename = 'ra19.7861250_dec-34.1916944.pkl'
fitting_run_result = pickle.load(open(picklename,'rb'))  #fitting_run_result is actually the fit_run in previous box.
#fitting_run_result.plot_final_qso_fit()
print(fitting_run_result.final_result_galaxy)


fitting_run_result.fitting_specify_class.plot_fitting_sets()

morph = fitting_run_result.cal_statmorph(obj_id=0, segm=fitting_run_result.fitting_specify_class.segm_deblend , if_plot = True)


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