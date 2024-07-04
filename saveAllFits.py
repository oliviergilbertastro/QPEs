#Load the saved fitting class, the fitting_run_result would be the loaded as fit_run() in previous fittings.
import pickle
import numpy as np
from galight_modif.tools.plot_tools import total_compare
from download_data import objects
import matplotlib.pyplot as plt
import copy
from matplotlib.colors import LogNorm
import matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('black')
import parse
import os
from datetime import datetime

save_folder = f"/Users/oliviergilbert/Desktop/QPEs/fits/auto_fits/"

def saveFit(picklename=None, savename=None):
    try:
        fitting_run_result = pickle.load(open("galight_fitruns/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.
    except:
        fitting_run_result = pickle.load(open("galight_fitruns/big_fits/"+picklename,'rb'))  #fitting_run_result is actually the fit_run in galightFitting.py.



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
    elif type == "Bulge+AGN":
        flux_list_1d = [data, model, ps_list[0], galaxy_list[0], galaxy_list[1], -model]
        label_list_1d = ['data', 'model', 'AGN', 'bulge', 'disk']
    elif type == "None":
        flux_list_1d = [data, model, galaxy_list[0], -model]
        label_list_1d = ['data', 'model', 'Sérsic']

    from galight_modif.fitting_process import ModelPlot
    data = fitting_run_result.fitting_specify_class.kwargs_data['image_data']
    if 'psf_error_map' in fitting_run_result.fitting_specify_class.kwargs_psf.keys():
        _modelPlot = ModelPlot(fitting_run_result.fitting_specify_class.kwargs_data_joint['multi_band_list'],
                                fitting_run_result.fitting_specify_class.kwargs_model, fitting_run_result.kwargs_result,
                                arrow_size=0.02, cmap_string="gist_heat", 
                                image_likelihood_mask_list=fitting_run_result.fitting_specify_class.kwargs_likelihood['image_likelihood_mask_list'] )    
        _, psf_error_map, _, _ = _modelPlot._imageModel.image_linear_solve(inv_bool=True, **fitting_run_result.kwargs_result)
        noise = np.sqrt(fitting_run_result.fitting_specify_class.kwargs_data['noise_map']**2+np.abs(psf_error_map[0]))
    else:
        noise = fitting_run_result.fitting_specify_class.kwargs_data['noise_map']
    
    ps_list = fitting_run_result.image_ps_list

    ps_image = np.zeros_like(data)

    #target_ID = f'{str(ra_dec[0])+symbol+str(ra_dec[1])}-{band}'
    target_ID = parse.parse("{}_{}.pkl", picklename)[0]
    print(target_ID)

    for i in range(len(ps_list)):
        ps_image = ps_image+ps_list[i]
    galaxy_list = fitting_run_result.image_host_list
    galaxy_image = np.zeros_like(data)
    for i in range(len(galaxy_list)):
        galaxy_image = galaxy_image+galaxy_list[i]
    model = ps_image + galaxy_image
    data_removePSF = data - ps_image
    norm_residual = (data - model)/noise
    if len(ps_list) != 0: #If there is an AGN
        flux_dict_2d = {'data':data, 'model':model, 'normalized residual':norm_residual}
        flux_dict_1d = {'data':data, 'model':model, 'AGN':ps_image, 'Sérsic':galaxy_image}
        fig = total_compare(list(flux_dict_2d.values()), list(flux_dict_2d.keys()), list(flux_dict_1d.values()), list(flux_dict_1d.keys()), deltaPix = fitting_run_result.fitting_specify_class.deltaPix,
                        zp=fitting_run_result.zp, if_annuli=False,
                        mask_image = fitting_run_result.fitting_specify_class.kwargs_likelihood['image_likelihood_mask_list'][0],
                        target_ID = target_ID, cmap=my_cmap, center_pos= [-fitting_run_result.final_result_ps[0]['ra_image'][0]/fitting_run_result.fitting_specify_class.deltaPix, 
                                                                        fitting_run_result.final_result_ps[0]['dec_image'][0]/fitting_run_result.fitting_specify_class.deltaPix], figsize=(13,4.5), show_plot=False )
    else: #If there is no AGN
        flux_dict_2d = {'data':data, 'model':model, 'normalized residual':norm_residual}
        flux_dict_1d = {'data':data, 'model':model, 'Sérsic':galaxy_image}
        if len(galaxy_list) == 2:
            flux_dict_1d = {'data':data, 'model':model, 'Bulge':galaxy_list[0], 'Disk':galaxy_list[1]}
        fig = total_compare(list(flux_dict_2d.values()), list(flux_dict_2d.keys()), list(flux_dict_1d.values()), list(flux_dict_1d.keys()), deltaPix = fitting_run_result.fitting_specify_class.deltaPix,
                        zp=fitting_run_result.zp, if_annuli=False,
                        mask_image = fitting_run_result.fitting_specify_class.kwargs_likelihood['image_likelihood_mask_list'][0],
                        target_ID = target_ID, cmap=my_cmap, figsize=(13,4.5), show_plot=False)
    #flux_dict_2d['data-point source'] = flux_dict_2d.pop('data$-$point source')
    fitting_run_result.flux_2d_out = flux_dict_2d
    fitting_run_result.flux_1d_out = flux_dict_1d
    if savename is None:
        savename = f"{picklename}_fit"
    
    
    fig.savefig(f"{save_folder}{savename}.pdf")
    return

from download_data import objects, comparisons, objects_names, objects_types, TDE_names, TDE_coords, TDE_types



if input("Save CO-ADDED SURVEY_PSF QPE hosts? [y/n]") == "y":
    time_dir = str(datetime.now()).replace(' ', '_').replace(':', '_')
    os.mkdir(f"{save_folder}/{time_dir}")
    for band in "r":
        for i in range(len(objects_names)):
            picklename=f"{objects_names[i]}_{band}-band_{'None'}_DESI_PSF.pkl"
            saveFit(picklename, savename=f"{time_dir}/sersicfit_{objects_names[i]}")

elif input("Save CO-ADDED SURVEY_PSF TDE hosts? [y/n]") == "y":
    time_dir = str(datetime.now()).replace(' ', '_').replace(':', '_')
    os.mkdir(f"{save_folder}/{time_dir}")
    for band in "r":
        for i in range(len(TDE_names)):
            picklename=f"{TDE_names[i]}_{band}-band_{'None'}_DESI_PSF.pkl"
            saveFit(picklename, savename=f"{time_dir}/sersicfit_{TDE_names[i]}")

elif input("Save CO-ADDED QPE hosts? [y/n]") == "y":
    time_dir = str(datetime.now()).replace(' ', '_').replace(':', '_')
    os.mkdir(f"{save_folder}/{time_dir}")
    for band in "r":
        for i in range(len(objects_names)):
            picklename=f"{objects_names[i]}_{band}-band_{objects_types[i]}_DESI.pkl"
            saveFit(picklename, savename=f"{time_dir}/sersicfit_{objects_names[i]}")