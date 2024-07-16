"""
Needed to make a lot of tables in the terminal and there is no support for it, so I made this program.
"""

import numpy as np
from math import floor, ceil
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.neighbors import KernelDensity
import seaborn as sns

def print_table(a, header=None, title=None, space_between_columns=2, space_between_rows=0, borders=1, header_color="yellow", border_color="grey", override_length=None):
    """
    Nicely print out a table

    a: array to print
    header: either an array of column names or a boolean. True -> first row will be header
    title: string, title that will be centered above the table
    space_between_columns: int, self-explanatory
    space_between_rows: int, self-explanatory
    """
    a = np.array(a, dtype=str)
    if type(header) == str or (type(header) == bool and header):
        pass
    elif (type(header) == bool and header == False) or header == None:
        header = None
    else:
        a = np.vstack((header, a))

    #Initialize the ascii characters to create the table depending on the borders parameter:
    if borders == None or borders == 0 or borders == False or borders == "none":
        characters = [" "," "," "," "," "," "]
    elif borders == "bold" or borders == 2:
        characters = ["═","║","╔","╗","╚","╝"]
    elif borders == 1 or borders == True or borders == "normal":
        characters = ["─","│","┌","┐","└","┘"]
    else:
        if type(borders) == str and len(borders) == 1:
            characters = [*(borders*6)]
        else:
            raise ValueError(f"Border style '{borders}' does not exist, use the keyword 'none', 'normal' or 'bold'.")
    
    possible_colors = ["black","red","green","yellow","blue","magenta","cyan","white"]
    #Initialize the colors:
    #Header color
    if header_color == None or header_color == "grey":
        header_color = "0"
    elif type(header_color) == str:
        header_color = header_color.lower()
        if header_color in possible_colors:
            header_color = str(possible_colors.index(header_color)+30)
        else:
            print(f"Color '{header_color}' not implemented, defaulting to grey.\nPossible colors are: {['grey']+possible_colors}")
            header_color = "0"
    else:
        raise ValueError(f"Parameter 'header_color' needs to be a string.")
    #Borders color
    if border_color == None or border_color == "grey":
        border_color = "0"
    elif type(border_color) == str:
        border_color = border_color.lower()
        if border_color in possible_colors:
            border_color = str(possible_colors.index(border_color)+30)
        else:
            print(f"Color '{border_color}' not implemented, defaulting to grey.\nPossible colors are: {['grey']+possible_colors}")
            border_color = "0"
    else:
        raise ValueError("Parameter 'border_color' needs to be a string.")

    for i in range(len(characters)):
        characters[i] = f"\x1b[{border_color}m{characters[i]}\x1b[0m"

    #Replace (None) elements with "-":
    a[a == "None"] = "-"

    #Get longest string in each column:
    column_maxes = []
    vfunc = np.vectorize(lambda x: len(x))
    a_lens = vfunc(a)
    for i in range(a.shape[1]):
        column_maxes.append(np.max(a_lens[:,i]))
    if override_length != None:
        column_maxes = override_length
    total_length = np.sum(column_maxes)+(len(column_maxes)-1)*space_between_columns #To include spaces between each column

    #Actually start printing table:
    top_and_bottom_bounds = (characters[2]+characters[0]*(total_length+2)+characters[3], characters[4]+characters[0]*(total_length+2)+characters[5])
    print()
    usable_length = total_length+4
    if title != None:
        title = floor((usable_length-len(title))/2)*" "+title+ceil((usable_length-len(title))/2)*" "
        print(f"\x1b[{header_color}m{title}\x1b[0m")
    print(top_and_bottom_bounds[0])
    #Print each row:
    for row in range(a.shape[0]):
        row_string = ""
        for column in range(a.shape[1]):
            row_string += a[row, column] + " "*(column_maxes[column]-a_lens[row,column])
            if column < a.shape[1]-1:
                row_string += " "*space_between_columns
        if row == 0 and header != None:
            row_string = f"\x1b[{header_color}m{row_string}\x1b[0m"
        row_string = f"{characters[1]} {row_string} {characters[1]}"
        if row != (a.shape[0]-1):
            row_string += f"\n{characters[1]} {' '*(total_length)} {characters[1]}"*space_between_rows
            if row == 0 and header != None:
                row_string += f"\n{characters[1]} {' '*(total_length)} {characters[1]}"
        print(row_string)
    print(top_and_bottom_bounds[1])
    print()








def myCornerPlot(data, labels=None, fontsize=15, smoothness=6):
    """
    data should be [data_set1, data_set2, ...] each containing multiple parameters
    """
    for i in range(len(data)-1):
        assert len(data[i]) == len(data[i+1])
    # Create the plot axes:
    fig = plt.figure(figsize=(10,8))
    plot_size = len(data[0])
    hist_axes = []
    corner_axes = []
    for i in range(plot_size):
        hist_axes.append(plt.subplot(plot_size,plot_size,i*plot_size+i+1))
        corner_axes.append([])
        for k in range(plot_size-(i+1)):
            if i == 0:
                corner_axes[i].append(plt.subplot(plot_size,plot_size,(i+k+1)*plot_size+(i+1),sharex=hist_axes[i]))
            else:
                corner_axes[i].append(plt.subplot(plot_size,plot_size,(i+k+1)*plot_size+(i+1),sharex=hist_axes[i],sharey=corner_axes[i-1][k+1]))
            if k != plot_size-(i+1)-1:
                corner_axes[i][k].get_xaxis().set_visible(False)
            if i != 0:
                corner_axes[i][k].get_yaxis().set_visible(False)
            corner_axes[i][k].xaxis.set_tick_params(labelsize=fontsize-2)
            corner_axes[i][k].yaxis.set_tick_params(labelsize=fontsize-2)
        if i == plot_size-1:
            hist_axes[i].get_yaxis().set_visible(False)
            hist_axes[i].xaxis.set_tick_params(labelsize=fontsize-2)
        else:
            hist_axes[i].get_xaxis().set_visible(False)
            hist_axes[i].get_yaxis().set_visible(False)

    # Show data in each plot:

    
    #Plot kernel histograms:
    for i in range(plot_size):
        if labels is not None:
            hist_axes[i].set_title(labels[i], fontsize=fontsize)
        x_min, x_max = np.min(data[0][i]), np.max(data[0][i])
        for j in range(len(data)):
            x_min = np.min(data[j][i]) if x_min > np.min(data[j][i]) else x_min
            x_max = np.max(data[j][i]) if x_max < np.max(data[j][i]) else x_max
        for j in range(len(data)-1):
            X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
            bandwidth = np.abs(x_max-x_min)/smoothness
            kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(data[j][i][:,np.newaxis])
            log_dens = kde.score_samples(X_plot)
            hist_axes[i].fill_between(X_plot[:, 0], np.exp(log_dens), fc=["blue","red","orange"][j%3], alpha=[0.4,0.4][j])

 
    for i in range(plot_size):
        for k in range(len(corner_axes[i])):
            for j in range(len(data)):
                corner_axes[i][k].plot(data[j][i], data[j][i+k+1], ["o","*","*"][j%3], color=["blue","red","red"][j%3], markersize=[8,7,7][j%3])
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.97, top=0.94, wspace=0, hspace=0)
    plt.show()
    return



if __name__ == "__main__":
    data = np.array([["potato", 5178, 13095, 3151],
            ["123", None, 1023, 51515],
            ["potato", 5178, 13012495, 51515],
            ["123", 123, 1024443, 51515],
            ["potaawddwto", 5178, 13095, 51515],
            ["123", "something", 1023, 51515],
            ["potato", 5178, 13095, 51515],
            ["123", None, 1023, 51515]])
    print_table(data,
                space_between_columns=4,
                space_between_rows=0,
                header=[f"Column{i}" for i in range(data.shape[1])],
                title="My Table",
                borders="bold",
                header_color="yellow",
                border_color="blue",
                )
    

def myFinalPlot(data, main_property=r"Sérsic index", referenceCatalogData=None, columns_compare=None, fontsize=15, smoothness=6, save_plot=True):
    """
    Originally made for the Sérsic index, but tweaked so it can accomodate the Bulge/Total light ratio
    """
    for i in range(len(data)-1):
        assert len(data[i]) == len(data[i+1])
    # data should be in the shape of [QPE data, TDE data]
    QPE_data, TDE_data = data
    QPE_data, TDE_data = np.array(QPE_data), np.array(TDE_data)

    # QPE and TDE data should be in the shape of np.array([n_sersic, m_star, m_BH]), with each one containing their uncertainties

    # Create the plot axes:
    fig = plt.figure(figsize=(6,6))
    gs = mpl.gridspec.GridSpec(10, 10, wspace=0.0, hspace=0.0)    
    n_hist_ax = fig.add_subplot(gs[0:2, 0:8]) # Sérsic index histogram
    mS_hist_ax =fig.add_subplot(gs[2:6, 8:10]) # Stellar mass histogram
    mBH_hist_ax = fig.add_subplot(gs[6:10, 8:10]) # Black hole mass histogram
    mS_ax = fig.add_subplot(gs[2:6, 0:8], sharex=n_hist_ax, sharey=mS_hist_ax) #  n vs log(M_star)
    mBH_ax = fig.add_subplot(gs[6:10, 0:8], sharex=n_hist_ax, sharey=mBH_hist_ax) # n vs log(M_BH)
    # Hide some of the axes ticks
    plt.setp(n_hist_ax.get_xticklabels(), visible=False)
    plt.setp(mS_hist_ax.get_xticklabels(), visible=False)
    plt.setp(mBH_hist_ax.get_xticklabels(), visible=False)
    plt.setp(mS_ax.get_xticklabels(), visible=False)
    plt.setp(n_hist_ax.get_yticklabels(), visible=False)
    plt.setp(mS_hist_ax.get_yticklabels(), visible=False)
    plt.setp(mBH_hist_ax.get_yticklabels(), visible=False)

    # fancy colors
    cmap = mpl.colormaps.get_cmap("viridis")
    naxes = len(fig.axes)

    #----------Make the histograms---------

    # n_sersic
    x_min = np.min(QPE_data[0,:,0]) if np.min(TDE_data[0,:,0]) > np.min(QPE_data[0,:,0]) else np.min(TDE_data[0,:,0])
    x_max = np.max(QPE_data[0,:,0]) if np.max(TDE_data[0,:,0]) < np.max(QPE_data[0,:,0]) else np.max(TDE_data[0,:,0])
    if referenceCatalogData is not None:
        x_min = x_min if np.min(referenceCatalogData[f"col_{columns_compare[0]}"]) > x_min else np.min(referenceCatalogData[f"col_{columns_compare[0]}"])
        x_max = x_max if np.max(referenceCatalogData[f"col_{columns_compare[0]}"]) < x_max else np.max(referenceCatalogData[f"col_{columns_compare[0]}"])
    X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(x_max-x_min)/smoothness
    if referenceCatalogData is not None:
        kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(np.array(referenceCatalogData[f"col_{columns_compare[0]}"])[:,np.newaxis])
        log_dens = kde.score_samples(X_plot)
        n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[0,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[0,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    # m_star
    y_min = np.min(QPE_data[1,:,0]) if np.min(TDE_data[1,:,0]) > np.min(QPE_data[1,:,0]) else np.min(TDE_data[1,:,0])
    y_max = np.max(QPE_data[1,:,0]) if np.max(TDE_data[1,:,0]) < np.max(QPE_data[1,:,0]) else np.max(TDE_data[1,:,0])
    if referenceCatalogData is not None:
        y_min = y_min if np.min(referenceCatalogData[f"col_{columns_compare[1]}"]) > y_min else np.min(referenceCatalogData[f"col_{columns_compare[1]}"])
        y_max = y_max if np.max(referenceCatalogData[f"col_{columns_compare[1]}"]) < y_max else np.max(referenceCatalogData[f"col_{columns_compare[1]}"])
    Y_plot = np.linspace(y_min, y_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(y_max-y_min)/smoothness
    if referenceCatalogData is not None:
        kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(np.array(referenceCatalogData[f"col_{columns_compare[1]}"])[:,np.newaxis])
        log_dens = kde.score_samples(Y_plot)
        mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[1,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[1,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    # m_bh
    y_min = np.min(QPE_data[2,:,0]) if np.min(TDE_data[2,:,0]) > np.min(QPE_data[2,:,0]) else np.min(TDE_data[2,:,0])
    y_max = np.max(QPE_data[2,:,0]) if np.max(TDE_data[2,:,0]) < np.max(QPE_data[2,:,0]) else np.max(TDE_data[2,:,0])
    if referenceCatalogData is not None:
        y_min = y_min if np.min(referenceCatalogData[f"col_{columns_compare[2]}"]) > y_min else np.min(referenceCatalogData[f"col_{columns_compare[2]}"])
        y_max = y_max if np.max(referenceCatalogData[f"col_{columns_compare[2]}"]) < y_max else np.max(referenceCatalogData[f"col_{columns_compare[2]}"])
    Y_plot = np.linspace(y_min, y_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(y_max-y_min)/smoothness
    if referenceCatalogData is not None:
        kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(np.array(referenceCatalogData[f"col_{columns_compare[2]}"])[:,np.newaxis])
        log_dens = kde.score_samples(Y_plot)
        mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[2,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[2,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    #----------Make the plots---------

    # n_sersic vs m_star
    if referenceCatalogData is not None:
        sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[1]}", fill=True, color="black", ax=mS_ax, label="Reference")
    mS_ax.errorbar(QPE_data[0,:,0], QPE_data[1,:,0], yerr=[QPE_data[1,:,1],QPE_data[1,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=8, label="QPE")
    mS_ax.errorbar(TDE_data[0,:,0], TDE_data[1,:,0], yerr=[TDE_data[1,:,1],TDE_data[1,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=7, mec="white", mew=0.5, label="TDE")
    mS_ax.legend(loc="upper left", fontsize=fontsize-3)
    # n_sersic vs m_bh
    if referenceCatalogData is not None:
        sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[2]}", fill=True, color="black", ax=mBH_ax)
    mBH_ax.errorbar(QPE_data[0,:,0], QPE_data[2,:,0], yerr=[QPE_data[2,:,1],QPE_data[2,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=8)
    mBH_ax.errorbar(TDE_data[0,:,0], TDE_data[2,:,0], yerr=[TDE_data[2,:,1],TDE_data[2,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", mec="white", mew=0.5, markersize=7)

    #----------Make the labels---------
    mBH_ax.set_xlabel(main_property, fontsize=fontsize)
    mS_ax.set_ylabel(r"$\log(M_\star)$", fontsize=fontsize)
    mBH_ax.set_ylabel(r"$\log(M_\mathrm{BH})$", fontsize=fontsize)
    mBH_ax.xaxis.set_tick_params(labelsize=fontsize-2)
    mS_ax.yaxis.set_tick_params(labelsize=fontsize-2)
    mBH_ax.yaxis.set_tick_params(labelsize=fontsize-2)
    #plt.subplots_adjust(left=0.06, bottom=0.06, right=0.97, top=0.94, wspace=0, hspace=0)
    if save_plot:
        plt.savefig(f"finalPlot_{columns_compare[0]}.pdf")
    plt.show()




def toLog(a):
    """Convert array of data and uncertainties in log base"""
    a = np.array(a)
    data, lo, hi = a[:,0], a[:,1], a[:,2]
    lo = np.abs(lo/(data*np.log(10)))
    hi = np.abs(hi/(data*np.log(10)))
    data = np.log10(data)
    return np.array([data, hi, lo]).T

def SDSS_objid_to_values(objid):

    # Determined from http://skyserver.sdss.org/dr7/en/help/docs/algorithm.asp?key=objID

    bin_objid = bin(objid)
    bin_objid = bin_objid[2:len(bin_objid)]
    bin_objid = bin_objid.zfill(64)

    empty = int( '0b' + bin_objid[0], base=0)
    skyVersion = int( '0b' + bin_objid[1:4+1], base=0)
    rerun = int( '0b' + bin_objid[5:15+1], base=0)
    run = int( '0b' + bin_objid[16:31+1], base=0)
    camcol = int( '0b' + bin_objid[32:34+1], base=0)
    firstField = int( '0b' + bin_objid[35+1], base=0)
    field = int( '0b' + bin_objid[36:47+1], base=0)
    object_num = int( '0b' + bin_objid[48:63+1], base=0)

    return skyVersion, rerun, run, camcol, field, object_num

from astropy.coordinates import SkyCoord
import astropy.units as u
def get_smallest_sep(pos, ras, decs):
    c1 = SkyCoord(pos[0]*u.deg, pos[1]*u.deg)
    c2 = SkyCoord(ras*u.deg, decs*u.deg)
    sep = (c1.separation(c2)).arcsec
    smallest_sep = min(sep)
    index = list(sep).index(smallest_sep)
    return index, smallest_sep


def get_smallest_sep_v2(pos, ras, decs):
    c1 = SkyCoord(pos[0]*u.deg, pos[1]*u.deg)
    c2 = SkyCoord(ras*u.deg, decs*u.deg)
    sep = (c1.separation(c2)).arcsec
    idx = (sep).argmin()
    return idx, sep[idx]



def print_color(message, color="yellow"):
    possible_colors = ["black","red","green","yellow","blue","magenta","cyan","white"]
    #Initialize the colors:
    #Header color
    if color == None or color == "grey":
        color = "0"
    elif type(color) == str:
        color = color.lower()
        if color in possible_colors:
            color = str(possible_colors.index(color)+30)
        else:
            print(f"Color '{color}' not implemented, defaulting to grey.\nPossible colors are: {['grey']+possible_colors}")
            color = "0"
    else:
        raise ValueError(f"Parameter 'header_color' needs to be a string.")
    print(f"\x1b[{color}m{message}\x1b[0m")
















import time
if __name__ == "__main__":

    start_time = time.time()
    print(get_smallest_sep([2.36286871e+02, -5.18003056e-01], [2.36247096e+02,2.36286871e+02,2.36336501e+02,2.15640297e+02], [-4.75263889e-01,-5.18003056e-01,-4.89098889e-01,1.06889111e+00]))
    print("\x1b[33m: --- %s seconds ---\x1b[0m" % (time.time() - start_time))
    start_time = time.time()
    print(get_smallest_sep_v2([2.36286871e+02, -5.18003056e-01], [2.36247096e+02,2.36286871e+02,2.36336501e+02,2.15640297e+02], [-4.75263889e-01,-5.18003056e-01,-4.89098889e-01,1.06889111e+00]))
    print("\x1b[33m: --- %s seconds ---\x1b[0m" % (time.time() - start_time))


    from paper_data import TDE_sersicIndices, TDE_stellar_masses_litterature, TDE_mBH
    from legacy_vs_legacy import add_0_uncertainties
    TDE_data = np.array([add_0_uncertainties(TDE_sersicIndices), toLog(TDE_stellar_masses_litterature), toLog(add_0_uncertainties(TDE_mBH))])
    data = np.array([TDE_data/2, TDE_data])
    #myFinalPlot(data, save_plot=False)