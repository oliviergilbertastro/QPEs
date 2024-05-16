#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
import numpy as np
from galight_modif.tools.plot_tools import total_compare

def loadRun(ra_dec, type="AGN"):
    if type in ["AGN", "agn", "Agn"]:
        type = "AGN"
    elif type in ["Bulge", "BULGE", "bulge"]:
        type = "Bulge"
    elif type in ["Bulge+AGN", "Bulge+Agn", "BULGE+AGN", "bulge+agn", "AGN+Bulge", "Agn+Bulge", "AGN+BULGE", "agn+bulge"]:
        type = "Bulge+AGN"
    else:
        raise ValueError(f"type {type} is not a supported fitting type")

    picklename = f'ra{str(ra_dec[0])}_dec{str(ra_dec[1])}_{type}.pkl'
    fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    fitting_run_result.run_diag()
    fitting_run_result.model_plot()
    if fitting_run_result.fitting_kwargs_list[-1][0] == 'MCMC':
        fitting_run_result.plot_params_corner()
        fitting_run_result.plot_flux_corner()  

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
    if type == "AGN":
        flux_list_1d = [data, model, ps_list[0], galaxy_list[0], -model]
        label_list_1d = ['data', 'model', 'AGN', 'galaxy']
    elif type == "Bulge":
        flux_list_1d = [data, model, galaxy_list[0], galaxy_list[1], -model]
        label_list_1d = ['data', 'model', 'bulge', 'disk']
    else:
        flux_list_1d = [data, model, ps_list[0], galaxy_list[0], galaxy_list[1], -model]
        label_list_1d = ['data', 'model', 'AGN', 'bulge', 'disk']
    total_compare(flux_list_2d, label_list_2d, flux_list_1d, label_list_1d, deltaPix = fitting_run_result.fitting_specify_class.deltaPix,
                        zp=fitting_run_result.zp, if_annuli=False, arrows= False, show_plot = True, mask_image = fitting_run_result.fitting_specify_class.kwargs_likelihood['image_likelihood_mask_list'][0],
                        target_ID = f'{str(ra_dec[0])+str(ra_dec[1])}', sum_rest = True)


    fitting_run_result.fitting_specify_class.plot_fitting_sets()
    obj_id = int(input('Which component to measure using statmorph?\n'))
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


from download_data import objects
objID = int(input("Enter the object ID you want to load [0-7]:\n"))
type = input("Enter the type of fitting you want to load:\n")
loadRun(objects[objID], type=type)