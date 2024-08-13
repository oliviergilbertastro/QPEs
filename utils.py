"""
Needed to make a lot of tables in the terminal and there is no support for it, so I made this program.
There are now various useful general functions in here, and some functions made specifically for the QPE vs TDE paper
"""

import numpy as np
from math import floor, ceil
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.neighbors import KernelDensity
import seaborn as sns
from tqdm import tqdm

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








def myCornerPlot(data, labels=None, units=None, fontsize=15, smoothness=6, linewidth=3, extremums=None, background_colors=["#f0f0f0","#969696","#252525"], levels=4, markersize=8, refCat=None, columns_compare=None, save_plot=None):
    """
    data should be [data_set1, data_set2, ...] each containing multiple parameters
    """
    for i in range(len(data)-1):
        assert len(data[i]) == len(data[i+1])

    # labels are now required (just put empty ones if you don't want them to appear)
    if extremums is None:
        extremums = {}
    for label in labels:
        if label not in extremums:
            extremums[label] = None
    
    
    # Create the plot axes:
    fig = plt.figure(figsize=(10,8))
    plot_size = len(data[0])
    if refCat is not None:
        assert plot_size == len(columns_compare)
    if labels is not None:
        assert plot_size == len(labels)
    if units is not None:
        assert plot_size == len(units)
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
        if refCat is not None:
            x_min = x_min if x_min < np.min(refCat[f"col_{columns_compare[i]}"]) else np.min(refCat[f"col_{columns_compare[i]}"])
            x_max = x_max if x_max > np.max(refCat[f"col_{columns_compare[i]}"]) else np.max(refCat[f"col_{columns_compare[i]}"])
        if extremums[labels[i]] is not None:
            x_min, x_max = extremums[labels[i]]
        for j in range(len(data)-1):
            X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
            bandwidth = np.abs(x_max-x_min)/smoothness
            if (refCat is not None) and j == 0:
                kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(np.array(refCat[f"col_{columns_compare[i]}"])[:,np.newaxis])
                log_dens = kde.score_samples(X_plot)
                #hist_axes[i].fill_between(X_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
                hist_axes[i].plot(X_plot[:, 0], np.exp(log_dens), color="black", linewidth=linewidth)
            kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(data[j][i][:,np.newaxis])
            log_dens = kde.score_samples(X_plot)
            hist_axes[i].fill_between(X_plot[:, 0], np.exp(log_dens), fc=["blue","red","orange"][j%3], alpha=[0.4,0.4][j])


 
    for i in range(plot_size):
        for k in range(len(corner_axes[i])):
            for j in range(len(data)):
                corner_axes[i][k].plot(data[j][i], data[j][i+k+1], ["o","*","*"][j%3], color=["blue","red","red"][j%3], markersize=[markersize,markersize-1,markersize-1][j%3])
            if refCat is not None:
                sns.kdeplot(refCat, x=f"col_{columns_compare[i]}", y=f"col_{columns_compare[i+k+1]}", fill=True, levels=levels, colors=background_colors, ax=corner_axes[i][k])
                #sns.kdeplot(refCat, x=f"col_{columns_compare[i]}", y=f"col_{columns_compare[i+k+1]}", fill=False, levels=levels, linewidths=0.5, color="black", ax=corner_axes[i][k])
        print_color(f"{labels[i]} :")
        print(f"QPE: {(np.min(data[0][i]), np.median(data[0][i]), np.max(data[0][i]))}")
        print(f"TDE: {(np.min(data[1][i]), np.median(data[1][i]), np.max(data[1][i]))}")
        print(f"ref: {(np.min(refCat[f'col_{columns_compare[i]}']), np.median(refCat[f'col_{columns_compare[i]}']), np.max(refCat[f'col_{columns_compare[i]}']))}")
    
    # Make units labels and set axis limits:
    if units is None:
        units = [" " for i in range(plot_size)]
    for i in range(plot_size):
        if i < plot_size-1:
            if i > 0:
                corner_axes[0][i-1].set_ylabel(labels[i]+"\n"+units[i], fontsize=fontsize-1)
                if extremums[labels[i]] is not None:
                    corner_axes[0][i-1].set_ylim(*extremums[labels[i]])
            corner_axes[i][-1].set_xlabel(labels[i]+"\n"+units[i], fontsize=fontsize-1)
            if extremums[labels[i]] is not None:
                corner_axes[i][-1].set_xlim(*extremums[labels[i]])
        else:
            corner_axes[0][i-1].set_ylabel(labels[i]+"\n"+units[i], fontsize=fontsize-1)
            if extremums[labels[i]] is not None:
                corner_axes[0][i-1].set_ylim(*extremums[labels[i]])
            hist_axes[i].set_xlabel(labels[i]+"\n"+units[i], fontsize=fontsize-1)
            if extremums[labels[i]] is not None:
                hist_axes[i].set_xlim(*extremums[labels[i]])
    
    plt.subplots_adjust(left=0.095, bottom=0.1, right=0.99, top=0.955, wspace=0, hspace=0)
    if save_plot is not None:
        plt.savefig(f"{save_plot}.pdf")
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
    

def myFinalPlot(data, main_property=r"Sérsic index", referenceCatalogData=None, columns_compare=None, fontsize=15, smoothness=6, save_plot=None, markersize=9, extremums=None, background_colors = ["#f0f0f0","#969696","#252525"], linewidth=2, levels = 4):
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
    fig = plt.figure(figsize=(4,4))
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

    # param
    x_min = np.min(QPE_data[0,:,0]) if np.min(TDE_data[0,:,0]) > np.min(QPE_data[0,:,0]) else np.min(TDE_data[0,:,0])
    x_max = np.max(QPE_data[0,:,0]) if np.max(TDE_data[0,:,0]) < np.max(QPE_data[0,:,0]) else np.max(TDE_data[0,:,0])
    if referenceCatalogData is not None:
        x_min = x_min if np.min(referenceCatalogData[f"col_{columns_compare[0]}"]) > x_min else np.min(referenceCatalogData[f"col_{columns_compare[0]}"])
        x_max = x_max if np.max(referenceCatalogData[f"col_{columns_compare[0]}"]) < x_max else np.max(referenceCatalogData[f"col_{columns_compare[0]}"])
    if extremums["param"] is not None:
        x_min, x_max = extremums["param"]
    X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(x_max-x_min)/smoothness
    if referenceCatalogData is not None:
        kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(np.array(referenceCatalogData[f"col_{columns_compare[0]}"])[:,np.newaxis])
        log_dens = kde.score_samples(X_plot)
        #n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
        n_hist_ax.plot(X_plot[:, 0], np.exp(log_dens), color="black", linewidth=linewidth)
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
    if extremums["m_star"] is not None:
        y_min, y_max = extremums["m_star"]
    Y_plot = np.linspace(y_min, y_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(y_max-y_min)/smoothness
    if referenceCatalogData is not None:
        kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(np.array(referenceCatalogData[f"col_{columns_compare[1]}"])[:,np.newaxis])
        log_dens = kde.score_samples(Y_plot)
        #mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
        mS_hist_ax.plot(np.exp(log_dens), Y_plot[:, 0], color="black", linewidth=linewidth)
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
    if extremums["m_bh"] is not None:
        y_min, y_max = extremums["m_bh"]
    Y_plot = np.linspace(y_min, y_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(y_max-y_min)/smoothness
    if referenceCatalogData is not None:
        kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(np.array(referenceCatalogData[f"col_{columns_compare[2]}"])[:,np.newaxis])
        log_dens = kde.score_samples(Y_plot)
        #mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
        mBH_hist_ax.plot(np.exp(log_dens), Y_plot[:, 0], color="black", linewidth=linewidth)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[2,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[2,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    #----------Make the plots---------

    # param vs m_star
    if referenceCatalogData is not None:
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[1]}", fill=True, levels=levels, color="black", ax=mS_ax, label="Reference")
        sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[1]}", fill=True, colors=background_colors, levels=levels, ax=mS_ax, label="Reference")
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[1]}", fill=False, levels=7, linewidths=0.5, color="black", ax=mS_ax, label="Reference")
    mS_ax.errorbar(QPE_data[0,:,0], QPE_data[1,:,0], yerr=[QPE_data[1,:,1],QPE_data[1,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize, label="QPE")
    mS_ax.errorbar(TDE_data[0,:,0], TDE_data[1,:,0], yerr=[TDE_data[1,:,1],TDE_data[1,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1, mec="white", mew=0.5, label="TDE")
    mS_ax.legend(loc="lower right", fontsize=fontsize-3)
    # param vs m_bh
    if referenceCatalogData is not None:
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[2]}", fill=True, levels=7, color="black", ax=mBH_ax)
        sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[2]}", fill=True, colors=background_colors, levels=levels, ax=mBH_ax, label="Reference")
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[2]}", fill=False, levels=7, linewidths=0.5, color="black", ax=mBH_ax)
    mBH_ax.errorbar(QPE_data[0,:,0], QPE_data[2,:,0], yerr=[QPE_data[2,:,1],QPE_data[2,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize)
    mBH_ax.errorbar(TDE_data[0,:,0], TDE_data[2,:,0], yerr=[TDE_data[2,:,1],TDE_data[2,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", mec="white", mew=0.5, markersize=markersize-1)

    #----------Make the labels---------
    mBH_ax.set_xlabel(main_property, fontsize=fontsize)
    mS_ax.set_ylabel(r"$\log(M_\star)$", fontsize=fontsize)
    mBH_ax.set_ylabel(r"$\log(M_\mathrm{BH})$", fontsize=fontsize)
    mBH_ax.xaxis.set_tick_params(labelsize=fontsize-2)
    mS_ax.yaxis.set_tick_params(labelsize=fontsize-2)
    mBH_ax.yaxis.set_tick_params(labelsize=fontsize-2)

    if extremums["param"] is not None:
        n_hist_ax.set_xlim(*extremums["param"])
    if extremums["m_star"] is not None:
        mS_hist_ax.set_ylim(*extremums["m_star"])
    if extremums["m_bh"] is not None:
        mBH_hist_ax.set_ylim(*extremums["m_bh"])


    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.975, top=0.975, wspace=0, hspace=0)
    if save_plot is not None:
        plt.savefig(f"{save_plot}.pdf")
    plt.show()


def redshiftMass(data, referenceCatalogData=None, columns_compare=None, kernelDensitiesReference=False, fontsize=15, smoothness=6, referenceSmoothness=12, save_plot=None, markersize=9, extremums=None, bins=50, background_colors = ["#f0f0f0","#969696","#252525"], linewidth=2, levels = 4):
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
    mS_hist_ax =fig.add_subplot(gs[2:10, 8:10]) # Stellar mass histogram
    mS_ax = fig.add_subplot(gs[2:10, 0:8], sharex=n_hist_ax, sharey=mS_hist_ax) #  n vs log(M_star)
    # Hide some of the axes ticks
    plt.setp(n_hist_ax.get_xticklabels(), visible=False)
    plt.setp(mS_hist_ax.get_xticklabels(), visible=False)
    plt.setp(n_hist_ax.get_yticklabels(), visible=False)
    plt.setp(mS_hist_ax.get_yticklabels(), visible=False)



    #----------Make the histograms---------

    # param
    x_min = np.min(QPE_data[0,:,0]) if np.min(TDE_data[0,:,0]) > np.min(QPE_data[0,:,0]) else np.min(TDE_data[0,:,0])
    x_max = np.max(QPE_data[0,:,0]) if np.max(TDE_data[0,:,0]) < np.max(QPE_data[0,:,0]) else np.max(TDE_data[0,:,0])
    if referenceCatalogData is not None:
        x_min = x_min if np.min(referenceCatalogData[f"col_{columns_compare[0]}"]) > x_min else np.min(referenceCatalogData[f"col_{columns_compare[0]}"])
        x_max = x_max if np.max(referenceCatalogData[f"col_{columns_compare[0]}"]) < x_max else np.max(referenceCatalogData[f"col_{columns_compare[0]}"])
    if extremums["param"] is not None:
        x_min, x_max = extremums["param"]
    X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(x_max-x_min)/smoothness
    if referenceCatalogData is not None:
        if kernelDensitiesReference:
            kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(x_max-x_min)/referenceSmoothness).fit(np.array(referenceCatalogData[f"col_{columns_compare[0]}"])[:,np.newaxis])
            log_dens = kde.score_samples(X_plot)
            #n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
            n_hist_ax.plot(X_plot[:, 0], np.exp(log_dens), color="black", linewidth=linewidth)
            print(np.sum(np.exp(log_dens)))
        else:
            heights, bin = np.histogram(referenceCatalogData[f"col_{columns_compare[0]}"], bins=bins, range=(x_min,x_max), density=True)
            n_hist_ax.plot(np.linspace(x_min,x_max,bins), heights, color="black", linewidth=2)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[0,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    print(np.sum(np.exp(log_dens)))
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[0,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)
    print(np.sum(np.exp(log_dens)))

    # m_star
    y_min = np.min(QPE_data[1,:,0]) if np.min(TDE_data[1,:,0]) > np.min(QPE_data[1,:,0]) else np.min(TDE_data[1,:,0])
    y_max = np.max(QPE_data[1,:,0]) if np.max(TDE_data[1,:,0]) < np.max(QPE_data[1,:,0]) else np.max(TDE_data[1,:,0])
    if referenceCatalogData is not None:
        y_min = y_min if np.min(referenceCatalogData[f"col_{columns_compare[1]}"]) > y_min else np.min(referenceCatalogData[f"col_{columns_compare[1]}"])
        y_max = y_max if np.max(referenceCatalogData[f"col_{columns_compare[1]}"]) < y_max else np.max(referenceCatalogData[f"col_{columns_compare[1]}"])
    if extremums["m_star"] is not None:
        y_min, y_max = extremums["m_star"]
    Y_plot = np.linspace(y_min, y_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(y_max-y_min)/smoothness
    if referenceCatalogData is not None:
        if kernelDensitiesReference:
            kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(y_max-y_min)/referenceSmoothness).fit(np.array(referenceCatalogData[f"col_{columns_compare[1]}"])[:,np.newaxis])
            log_dens = kde.score_samples(Y_plot)
            #mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
            mS_hist_ax.plot(np.exp(log_dens), Y_plot[:, 0], color="black", linewidth=linewidth)
            # We don't want to smooth twice:
        else:
            heights, bin = np.histogram(referenceCatalogData[f"col_{columns_compare[1]}"], bins=bins, range=(y_min,y_max), density=True)
            mS_hist_ax.plot(heights, np.linspace(y_min,y_max,bins), color="black", linewidth=2)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[1,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[1,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)


    #----------Make the plots---------

    # param vs m_star
    if referenceCatalogData is not None:
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[1]}", fill=True, levels=levels, color="black", ax=mS_ax, label="Reference")
        sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[1]}", fill=True, colors=background_colors, levels=levels, ax=mS_ax, label="Reference")
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0]}", y=f"col_{columns_compare[1]}", fill=False, levels=7, linewidths=0.5, color="black", ax=mS_ax, label="Reference")
    mS_ax.errorbar(QPE_data[0,:,0], QPE_data[1,:,0], yerr=[QPE_data[1,:,1],QPE_data[1,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize, label="QPE")
    mS_ax.errorbar(TDE_data[0,:,0], TDE_data[1,:,0], yerr=[TDE_data[1,:,1],TDE_data[1,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1, mec="white", mew=0.5, label="TDE")
    mS_ax.legend(loc="lower right", fontsize=fontsize-3)

    #----------Make the labels---------
    mS_ax.set_xlabel(r"Redshift $z$", fontsize=fontsize)
    mS_ax.set_ylabel(r"$\log(M_\star)$", fontsize=fontsize)
    mS_ax.xaxis.set_tick_params(labelsize=fontsize-2)
    mS_ax.yaxis.set_tick_params(labelsize=fontsize-2)

    if extremums["param"] is not None:
        n_hist_ax.set_xlim(*extremums["param"])
    if extremums["m_star"] is not None:
        mS_hist_ax.set_ylim(*extremums["m_star"])


    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.975, top=0.975, wspace=0, hspace=0)
    if save_plot is not None:
        plt.savefig(f"{save_plot}.pdf")
    plt.show()


def myCombinedFinalPlot(data, referenceCatalogData=None, columns_compare=None, kernelDensitiesReference=False, extremums=None, fontsize=15, markersize=8, smoothness=6, referenceSmoothness=12, save_plot=None, bins=50, background_colors = ["#f0f0f0","#969696","#252525"], linewidth=2, levels = 4):
    """
    Originally made for the Sérsic index, but tweaked so it can accomodate the Bulge/Total light ratio
    """
    for i in range(len(data)-1):
        assert len(data[i]) == len(data[i+1])
    if extremums is None:
        extremums = {}
    for param in ["n_sersic", "bt_ratio", "ssmd", "m_star", "m_bh"]:
        if param not in extremums:
            extremums[param] = None
    # data should be in the shape of [QPE data, TDE data]
    QPE_data, TDE_data = data
    QPE_data, TDE_data = np.array(QPE_data), np.array(TDE_data)

    # QPE and TDE data should be in the shape of np.array(n_sersic, bt_ratio, ssmd, m_star, m_BH]), with each one containing their uncertainties
    # Columns compare should be [[n_sersic, bt_ratio, ssmd], m_star, m_BH]

    # Create the plot axes:
    fig = plt.figure(figsize=(10.5,7.5))
    gs = mpl.gridspec.GridSpec(10, 14, wspace=0.0, hspace=0.0)    
    n_hist_ax = fig.add_subplot(gs[0:2, 0:4]) # Sérsic index histogram
    bt_hist_ax = fig.add_subplot(gs[0:2, 4:8]) # B/T histogram
    ssmd_hist_ax = fig.add_subplot(gs[0:2, 8:12]) # SSMD histogram
    mS_hist_ax =fig.add_subplot(gs[2:6, 12:14]) # Stellar mass histogram
    mBH_hist_ax = fig.add_subplot(gs[6:10, 12:14]) # Black hole mass histogram
    n_mS_ax = fig.add_subplot(gs[2:6, 0:4], sharex=n_hist_ax, sharey=mS_hist_ax) #  n vs log(M_star)
    n_mBH_ax = fig.add_subplot(gs[6:10, 0:4], sharex=n_hist_ax, sharey=mBH_hist_ax) # n vs log(M_BH)
    bt_mS_ax = fig.add_subplot(gs[2:6, 4:8], sharex=bt_hist_ax, sharey=mS_hist_ax) #  bt vs log(M_star)
    bt_mBH_ax = fig.add_subplot(gs[6:10, 4:8], sharex=bt_hist_ax, sharey=mBH_hist_ax) # bt vs log(M_BH)
    ssmd_mS_ax = fig.add_subplot(gs[2:6, 8:12], sharex=ssmd_hist_ax, sharey=mS_hist_ax) #  ssmd vs log(M_star)
    ssmd_mBH_ax = fig.add_subplot(gs[6:10, 8:12], sharex=ssmd_hist_ax, sharey=mBH_hist_ax) # ssmd vs log(M_BH)
    # Hide some of the axes ticks
    plt.setp(n_hist_ax.get_xticklabels(), visible=False)
    plt.setp(bt_hist_ax.get_xticklabels(), visible=False)
    plt.setp(ssmd_hist_ax.get_xticklabels(), visible=False)
    plt.setp(mS_hist_ax.get_xticklabels(), visible=False)
    plt.setp(mBH_hist_ax.get_xticklabels(), visible=False)
    plt.setp(n_mS_ax.get_xticklabels(), visible=False)
    plt.setp(bt_mS_ax.get_xticklabels(), visible=False)
    plt.setp(ssmd_mS_ax.get_xticklabels(), visible=False)

    plt.setp(n_hist_ax.get_yticklabels(), visible=False)
    plt.setp(bt_hist_ax.get_yticklabels(), visible=False)
    plt.setp(ssmd_hist_ax.get_yticklabels(), visible=False)
    plt.setp(mS_hist_ax.get_yticklabels(), visible=False)
    plt.setp(mBH_hist_ax.get_yticklabels(), visible=False)
    plt.setp(bt_mS_ax.get_yticklabels(), visible=False)
    plt.setp(ssmd_mS_ax.get_yticklabels(), visible=False)
    plt.setp(bt_mBH_ax.get_yticklabels(), visible=False)
    plt.setp(ssmd_mBH_ax.get_yticklabels(), visible=False)


    #----------Make the histograms---------

    # n_sersic
    x_min = np.min(QPE_data[0,:,0]) if np.min(TDE_data[0,:,0]) > np.min(QPE_data[0,:,0]) else np.min(TDE_data[0,:,0])
    x_max = np.max(QPE_data[0,:,0]) if np.max(TDE_data[0,:,0]) < np.max(QPE_data[0,:,0]) else np.max(TDE_data[0,:,0])
    if referenceCatalogData is not None:
        x_min = x_min if np.min(referenceCatalogData[f"col_{columns_compare[0][0]}"]) > x_min else np.min(referenceCatalogData[f"col_{columns_compare[0][0]}"])
        x_max = x_max if np.max(referenceCatalogData[f"col_{columns_compare[0][0]}"]) < x_max else np.max(referenceCatalogData[f"col_{columns_compare[0][0]}"])
    if extremums["n_sersic"] is not None:
        x_min, x_max = extremums["n_sersic"]
    X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(x_max-x_min)/smoothness
    if referenceCatalogData is not None:
        if kernelDensitiesReference:
            kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(x_max-x_min)/referenceSmoothness).fit(np.array(referenceCatalogData[f"col_{columns_compare[0][0]}"])[:,np.newaxis])
            log_dens = kde.score_samples(X_plot)
            #n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
            n_hist_ax.plot(X_plot[:, 0], np.exp(log_dens), color="black", linewidth=linewidth)
        else:
            heights, bin = np.histogram(referenceCatalogData[f"col_{columns_compare[0][0]}"], bins=bins, range=(x_min, x_max), density=True)
            n_hist_ax.plot(np.linspace(x_min,x_max,bins), heights, color="black", linewidth=2)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[0,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[0,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    n_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    # B/T ratio
    x_min = np.min(QPE_data[1,:,0]) if np.min(TDE_data[1,:,0]) > np.min(QPE_data[1,:,0]) else np.min(TDE_data[1,:,0])
    x_max = np.max(QPE_data[1,:,0]) if np.max(TDE_data[1,:,0]) < np.max(QPE_data[1,:,0]) else np.max(TDE_data[1,:,0])
    if referenceCatalogData is not None:
        x_min = x_min if np.min(referenceCatalogData[f"col_{columns_compare[0][1]}"]) > x_min else np.min(referenceCatalogData[f"col_{columns_compare[0][1]}"])
        x_max = x_max if np.max(referenceCatalogData[f"col_{columns_compare[0][1]}"]) < x_max else np.max(referenceCatalogData[f"col_{columns_compare[0][1]}"])
    if extremums["bt_ratio"] is not None:
        x_min, x_max = extremums["bt_ratio"]
    X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(x_max-x_min)/smoothness
    if referenceCatalogData is not None:
        if kernelDensitiesReference:
            kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(x_max-x_min)/referenceSmoothness).fit(np.array(referenceCatalogData[f"col_{columns_compare[0][1]}"])[:,np.newaxis])
            log_dens = kde.score_samples(X_plot)
            #bt_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
            bt_hist_ax.plot(X_plot[:, 0], np.exp(log_dens), color="black", linewidth=linewidth)
        else:
            heights, bin = np.histogram(referenceCatalogData[f"col_{columns_compare[0][1]}"], bins=bins, range=(x_min, x_max), density=True)
            bt_hist_ax.plot(np.linspace(x_min,x_max,bins), heights, color="black", linewidth=2)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[1,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    bt_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[1,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    bt_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    # SSMD
    x_min = np.min(QPE_data[2,:,0]) if np.min(TDE_data[2,:,0]) > np.min(QPE_data[2,:,0]) else np.min(TDE_data[2,:,0])
    x_max = np.max(QPE_data[2,:,0]) if np.max(TDE_data[2,:,0]) < np.max(QPE_data[2,:,0]) else np.max(TDE_data[2,:,0])
    if referenceCatalogData is not None:
        x_min = x_min if np.min(referenceCatalogData[f"col_{columns_compare[0][2]}"]) > x_min else np.min(referenceCatalogData[f"col_{columns_compare[0][2]}"])
        x_max = x_max if np.max(referenceCatalogData[f"col_{columns_compare[0][2]}"]) < x_max else np.max(referenceCatalogData[f"col_{columns_compare[0][2]}"])
    if extremums["ssmd"] is not None:
        x_min, x_max = extremums["ssmd"]
    X_plot = np.linspace(x_min, x_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(x_max-x_min)/smoothness
    if referenceCatalogData is not None:
        if kernelDensitiesReference:
            kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(x_max-x_min)/referenceSmoothness).fit(np.array(referenceCatalogData[f"col_{columns_compare[0][2]}"])[:,np.newaxis])
            log_dens = kde.score_samples(X_plot)
            #ssmd_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
            ssmd_hist_ax.plot(X_plot[:, 0], np.exp(log_dens), color="black", linewidth=linewidth)
        else:
            heights, bin = np.histogram(referenceCatalogData[f"col_{columns_compare[0][2]}"], bins=bins, range=(x_min, x_max), density=True)
            ssmd_hist_ax.plot(np.linspace(x_min,x_max,bins), heights, color="black", linewidth=2)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[2,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    ssmd_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[2,:,0][:,np.newaxis])
    log_dens = kde.score_samples(X_plot)
    ssmd_hist_ax.fill_between(X_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    # m_star
    y_min = np.min(QPE_data[3,:,0]) if np.min(TDE_data[3,:,0]) > np.min(QPE_data[3,:,0]) else np.min(TDE_data[3,:,0])
    y_max = np.max(QPE_data[3,:,0]) if np.max(TDE_data[3,:,0]) < np.max(QPE_data[3,:,0]) else np.max(TDE_data[3,:,0])
    if referenceCatalogData is not None:
        y_min = y_min if np.min(referenceCatalogData[f"col_{columns_compare[1]}"]) > y_min else np.min(referenceCatalogData[f"col_{columns_compare[1]}"])
        y_max = y_max if np.max(referenceCatalogData[f"col_{columns_compare[1]}"]) < y_max else np.max(referenceCatalogData[f"col_{columns_compare[1]}"])
    if extremums["m_star"] is not None:
        y_min, y_max = extremums["m_star"]
    Y_plot = np.linspace(y_min, y_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(y_max-y_min)/smoothness
    if referenceCatalogData is not None:
        if kernelDensitiesReference:
            kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(y_max-y_min)/referenceSmoothness).fit(np.array(referenceCatalogData[f"col_{columns_compare[1]}"])[:,np.newaxis])
            log_dens = kde.score_samples(Y_plot)
            #mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
            mS_hist_ax.plot(np.exp(log_dens), Y_plot[:, 0], color="black", linewidth=linewidth)
        else:
            heights, bin = np.histogram(referenceCatalogData[f"col_{columns_compare[1]}"], bins=bins, range=(y_min, y_max), density=True)
            mS_hist_ax.plot(heights, np.linspace(y_min,y_max,bins), color="black", linewidth=2)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[3,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[3,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mS_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    # m_bh
    y_min = np.min(QPE_data[4,:,0]) if np.min(TDE_data[4,:,0]) > np.min(QPE_data[4,:,0]) else np.min(TDE_data[4,:,0])
    y_max = np.max(QPE_data[4,:,0]) if np.max(TDE_data[4,:,0]) < np.max(QPE_data[4,:,0]) else np.max(TDE_data[4,:,0])
    if referenceCatalogData is not None:
        y_min = y_min if np.min(referenceCatalogData[f"col_{columns_compare[2]}"]) > y_min else np.min(referenceCatalogData[f"col_{columns_compare[2]}"])
        y_max = y_max if np.max(referenceCatalogData[f"col_{columns_compare[2]}"]) < y_max else np.max(referenceCatalogData[f"col_{columns_compare[2]}"])
    if extremums["m_bh"] is not None:
        y_min, y_max = extremums["m_bh"]
    Y_plot = np.linspace(y_min, y_max, 1000)[:,np.newaxis]
    bandwidth = np.abs(y_max-y_min)/smoothness
    if referenceCatalogData is not None:
        if kernelDensitiesReference:
            kde = KernelDensity(kernel="gaussian", bandwidth=np.abs(y_max-y_min)/referenceSmoothness).fit(np.array(referenceCatalogData[f"col_{columns_compare[2]}"])[:,np.newaxis])
            log_dens = kde.score_samples(Y_plot)
            #mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="grey", alpha=0.4)
            mBH_hist_ax.plot(np.exp(log_dens), Y_plot[:, 0], color="black", linewidth=linewidth)
        else:
            heights, bin = np.histogram(referenceCatalogData[f"col_{columns_compare[2]}"], bins=bins, range=(y_min, y_max), density=True)
            mBH_hist_ax.plot(heights, np.linspace(y_min,y_max,bins), color="black", linewidth=2)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(QPE_data[4,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="blue", alpha=0.4)
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(TDE_data[4,:,0][:,np.newaxis])
    log_dens = kde.score_samples(Y_plot)
    mBH_hist_ax.fill_betweenx(Y_plot[:, 0], np.exp(log_dens), fc="red", alpha=0.4)

    #----------Make the plots---------


    
    # n_sersic vs m_star
    if referenceCatalogData is not None:
        kde00 = sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][0]}", y=f"col_{columns_compare[1]}", fill=True, colors=background_colors, levels=levels, ax=n_mS_ax, label="Reference")
        for contour, color in zip(kde00.collections, background_colors):
            contour.set_facecolor(color)
            #contour.set_color(color)
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][0]}", y=f"col_{columns_compare[1]}", fill=False, levels=levels, linewidths=0.5, color="black", ax=n_mS_ax, label="Reference")
    n_mS_ax.errorbar(QPE_data[0,:,0], QPE_data[3,:,0], yerr=[QPE_data[3,:,1],QPE_data[3,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize, label="QPE")
    n_mS_ax.errorbar(TDE_data[0,:,0], TDE_data[3,:,0], yerr=[TDE_data[3,:,1],TDE_data[3,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1, label="TDE")
    n_mS_ax.legend(loc="lower right", fontsize=fontsize-3)
    # n_sersic vs m_bh
    if referenceCatalogData is not None:
        kde01 = sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][0]}", y=f"col_{columns_compare[2]}", fill=True, colors=background_colors, levels=levels, ax=n_mBH_ax)
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][0]}", y=f"col_{columns_compare[2]}", fill=False, levels=levels, linewidths=0.5, color="black", ax=n_mBH_ax)
    n_mBH_ax.errorbar(QPE_data[0,:,0], QPE_data[4,:,0], yerr=[QPE_data[4,:,1],QPE_data[4,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize)
    n_mBH_ax.errorbar(TDE_data[0,:,0], TDE_data[4,:,0], yerr=[TDE_data[4,:,1],TDE_data[4,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1)

    # bt vs m_star
    if referenceCatalogData is not None:
        kde10 = sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][1]}", y=f"col_{columns_compare[1]}", fill=True, colors=background_colors, levels=levels, ax=bt_mS_ax, label="Reference")
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][1]}", y=f"col_{columns_compare[1]}", fill=False, levels=levels, linewidths=0.5, color="black", ax=bt_mS_ax, label="Reference")
    bt_mS_ax.errorbar(QPE_data[1,:,0], QPE_data[3,:,0], yerr=[QPE_data[3,:,1],QPE_data[3,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize, label="QPE")
    bt_mS_ax.errorbar(TDE_data[1,:,0], TDE_data[3,:,0], yerr=[TDE_data[3,:,1],TDE_data[3,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1, label="TDE")
    # n_sersic vs m_bh
    if referenceCatalogData is not None:
        kde11 = sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][1]}", y=f"col_{columns_compare[2]}", fill=True, colors=background_colors, levels=levels, ax=bt_mBH_ax)
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][1]}", y=f"col_{columns_compare[2]}", fill=False, levels=levels, linewidths=0.5, color="black", ax=bt_mBH_ax)
    bt_mBH_ax.errorbar(QPE_data[1,:,0], QPE_data[4,:,0], yerr=[QPE_data[4,:,1],QPE_data[4,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize)
    bt_mBH_ax.errorbar(TDE_data[1,:,0], TDE_data[4,:,0], yerr=[TDE_data[4,:,1],TDE_data[4,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1)


    # ssmd vs m_star
    if referenceCatalogData is not None:
        kde20 = sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][2]}", y=f"col_{columns_compare[1]}", fill=True, colors=background_colors, levels=levels, ax=ssmd_mS_ax, label="Reference")
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][2]}", y=f"col_{columns_compare[1]}", fill=False, levels=3, linewidths=0.5, color="black", ax=ssmd_mS_ax, label="Reference")
    ssmd_mS_ax.errorbar(QPE_data[2,:,0], QPE_data[3,:,0], yerr=[QPE_data[3,:,1],QPE_data[3,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize, label="QPE")
    ssmd_mS_ax.errorbar(TDE_data[2,:,0], TDE_data[3,:,0], yerr=[TDE_data[3,:,1],TDE_data[3,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1, label="TDE")
    # n_sersic vs m_bh
    if referenceCatalogData is not None:
        kde21 = sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][2]}", y=f"col_{columns_compare[2]}", fill=True, colors=background_colors, levels=levels, ax=ssmd_mBH_ax)
        #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][2]}", y=f"col_{columns_compare[2]}", fill=False, levels=3, linewidths=0.5, color="black", ax=ssmd_mBH_ax)
    ssmd_mBH_ax.errorbar(QPE_data[2,:,0], QPE_data[4,:,0], yerr=[QPE_data[4,:,1],QPE_data[4,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize)
    ssmd_mBH_ax.errorbar(TDE_data[2,:,0], TDE_data[4,:,0], yerr=[TDE_data[4,:,1],TDE_data[4,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1)



    #----------Make the labels---------
    n_mBH_ax.set_xlabel(r"$\text{Sérsic index } n$", fontsize=fontsize)
    bt_mBH_ax.set_xlabel(r"$(B/T)_g$", fontsize=fontsize)
    ssmd_mBH_ax.set_xlabel(r"$\log(\Sigma_{M_\star})$ $[M_\odot ~\mathrm{kpc}^{-2}]$", fontsize=fontsize)
    n_mS_ax.set_ylabel(r"$\log(M_\star)$ $[M_\odot]$", fontsize=fontsize)
    n_mBH_ax.set_ylabel(r"$\log(M_\mathrm{BH})$ $[M_\odot]$", fontsize=fontsize)
    n_mBH_ax.xaxis.set_tick_params(labelsize=fontsize-2)
    bt_mBH_ax.xaxis.set_tick_params(labelsize=fontsize-2)
    ssmd_mBH_ax.xaxis.set_tick_params(labelsize=fontsize-2)
    n_mS_ax.yaxis.set_tick_params(labelsize=fontsize-2)
    n_mBH_ax.yaxis.set_tick_params(labelsize=fontsize-2)

    if extremums["n_sersic"] is not None:
        n_hist_ax.set_xlim(*extremums["n_sersic"])
    if extremums["bt_ratio"] is not None:
        bt_hist_ax.set_xlim(*extremums["bt_ratio"])
    if extremums["ssmd"] is not None:
        ssmd_hist_ax.set_xlim(*extremums["ssmd"])
    if extremums["m_star"] is not None:
        mS_hist_ax.set_ylim(*extremums["m_star"])
    if extremums["m_bh"] is not None:
        mBH_hist_ax.set_ylim(*extremums["m_bh"])

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.99, top=0.975, wspace=0, hspace=0)
    if save_plot is not None:
        plt.savefig(f"{save_plot}.pdf")
    plt.show()


def makeLatexTable(names, redshifts, r50s, n_sersics, bt_ratios, ssmds, M_stars, references=None, filename="latexTable.txt", verbose=False):
    """
    columns are: NAME, REDSHIFT, R_50, N_SERSIC, BT_RATIO, SSMD, M_stars

    makes a .txt file ready to be copy&pasted into overleaf or your favorite LaTex editor near you.
    """
    total_string = ""
    each_lines = []
    if references is None:
        references = ["" for i in range(len(names))]
    if len(references) < len(names):
        references = references + "?"*(len(names)-len(references))
    print(len(names),len(redshifts), len(r50s[:,0]), len(n_sersics[:,0]), len(bt_ratios[:,0]), len(ssmds[:,0]), len(M_stars[:,0]), len(references))
    assert len(names) == len(redshifts) == len(r50s[:,0]) == len(n_sersics[:,0]) == len(bt_ratios[:,0]) == len(ssmds[:,0]) == len(M_stars[:,0]) == len(references)
    length = len(names)

    def LatexUncertainty(value):
        return f"${round(value[0],2)}_{r'{-'+str(round(value[1],2))+r'}'}^{r'{+'+str(round(value[2],2))+r'}'}$"

    for i in range(length):
        each_lines.append(names[i]+r"$^{\rm "+references[i]+"}$" + " & " + "$" + str(round(redshifts[i],3)) + "$" + " & " + LatexUncertainty(r50s[i]) + " & " + LatexUncertainty(n_sersics[i]) + " & " + LatexUncertainty(bt_ratios[i]) + " & " + LatexUncertainty(ssmds[i]) + " & " + LatexUncertainty(M_stars[i]) + r" \\" + "\n")        

    # swap lines around here easily
    onlyQPEs = np.array(each_lines)[[0,1,2,3,6,7]]
    QPE_TDEs = np.array(each_lines)[[4,5,8]]
    onlyTDEs = np.array(each_lines)[[9,10,11,12,13,14,15,16,17,18]]
    
    total_string += r"\multicolumn{7}{c}{\emph{QPE host galaxies}}\\"+"\n\hline\n\hline\n"
    for line in onlyQPEs:
        total_string += line
        if line != onlyQPEs[-1]:
            total_string += r"\vspace{2pt}"+"\n"
    total_string += "\hline\n"+r"\multicolumn{7}{c}{\emph{TDE+QPE host galaxies}}\\"+"\n\hline\n\hline\n"
    for line in QPE_TDEs:
        total_string += line
        if line != QPE_TDEs[-1]:
            total_string += r"\vspace{2pt}"+"\n"
    total_string += "\hline\n"+r"\multicolumn{7}{c}{\emph{TDE host galaxies}}\\"+"\n\hline\n\hline\n"
    for line in onlyTDEs:
        total_string += line
        if line != onlyTDEs[-1]:
            total_string += r"\vspace{2pt}"+"\n"

    with open(filename, "w") as text_file:
        text_file.write(total_string)
    if verbose:
        print(total_string)




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



def print_color(message, color="yellow", **kwargs):
    """print(), but with a color option"""
    possible_colors = ["black","red","green","yellow","blue","magenta","cyan","white"]
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
    print(f"\x1b[{color}m{message}\x1b[0m", **kwargs)



import copy

def recombine_arrays(data, lo, hi):
    """
    Recombine 2D arrays into a single 3D array
    data: array containing the data
    lo: array containing the lower uncertainties
    hi: array containing the high uncertainties
    """

    assert data.shape == lo.shape == hi.shape

    new_array = np.zeros((*(data.shape),3))
    new_array[:,:,0] = data
    new_array[:,:,1] = lo
    new_array[:,:,2] = hi

    return new_array

def cut_from_catalog(catalog, index, bounds, verbose=False):
    """
    catalog: numpy array
    index: int of index of parameter column which is under study here
    bounds: tuple of bound which we want to keep

    returns: numpy array with only remaining objects
    """
    catalog = np.array(catalog)
    if bounds[0] == None:
        good_indices_lo = catalog[:,index] == catalog[:,index]
    else:
        good_indices_lo = catalog[:,index] >= bounds[0]
    if bounds[1] == None:
        good_indices_hi = catalog[:,index] == catalog[:,index]
    else:
        good_indices_hi = catalog[:,index] <= bounds[1]
    good_indices = []
    for i in range(len(good_indices_lo)):
        good_indices.append(good_indices_lo[i] and good_indices_hi[i])
    cut_catalog = catalog[good_indices]
    if verbose:
        print(f"\x1b[31m{catalog.shape[0]-cut_catalog.shape[0]} objects cut\x1b[0m")
        print(f"\x1b[32m{cut_catalog.shape[0]} objects remaining\x1b[0m")
    return cut_catalog

def mergeCatalogs_withObjIDs(cat1,cat2,columnsToAdd=[0]):
    """
    Merge catalog1 with catalog2 assuming their first columns are the objIDs
    """
    good_indices = []
    properties_toAdd = []
    for i in range(len(columnsToAdd)):
        properties_toAdd.append([])
    for i in tqdm(range(len(cat1[:,0]))):
        try:
            index = list(cat2[:,0]).index(cat1[i,0])
            good_indices.append(i)
            for k in range(len(columnsToAdd)):
                properties_toAdd[k].append(cat2[index,columnsToAdd[k]])
            #print(f"{cat1[i,0]} vs {cat2[index,0]}")
        except:
            pass
    cat1 = cat1[[good_indices],:][0]
    for i in range(len(columnsToAdd)):
        cat1 = np.vstack((cat1.T, np.array(properties_toAdd[i]))).T
    return cat1

def cut_from_array(a, indices):
    indices = [(i if i >= 0 else len(a)+i) for i in indices]
    return a[[(i not in indices) for i in range(len(a))]]


if __name__ == "__main__":
    data = np.loadtxt("TDE_allRelevantData_0.txt")
    lo = np.loadtxt("TDE_allRelevantData_1.txt")
    hi = np.loadtxt("TDE_allRelevantData_2.txt")
    print(recombine_arrays(data, lo, hi)[:,0,:])

def add_0_uncertainties(a):
    a = np.array(a)
    placeholder = np.zeros((a.shape[0],3))
    placeholder[:,0] = a
    a = placeholder
    return a

import time
if __name__ == "__main__":

    start_time = time.time()
    print(get_smallest_sep([2.36286871e+02, -5.18003056e-01], [2.36247096e+02,2.36286871e+02,2.36336501e+02,2.15640297e+02], [-4.75263889e-01,-5.18003056e-01,-4.89098889e-01,1.06889111e+00]))
    print("\x1b[33m: --- %s seconds ---\x1b[0m" % (time.time() - start_time))
    start_time = time.time()
    print(get_smallest_sep_v2([2.36286871e+02, -5.18003056e-01], [2.36247096e+02,2.36286871e+02,2.36336501e+02,2.15640297e+02], [-4.75263889e-01,-5.18003056e-01,-4.89098889e-01,1.06889111e+00]))
    print("\x1b[33m: --- %s seconds ---\x1b[0m" % (time.time() - start_time))

    from paper_data import TDE_sersicIndices, TDE_stellar_masses_litterature, TDE_mBH
    TDE_data = np.array([add_0_uncertainties(TDE_sersicIndices), toLog(TDE_stellar_masses_litterature), toLog(add_0_uncertainties(TDE_mBH))])
    data = np.array([TDE_data/2, TDE_data])
    #myFinalPlot(data, save_plot=False)