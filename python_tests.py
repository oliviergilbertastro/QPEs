
if False:
    func = lambda x1, x2: 4*x1+x2

    print(func(2,3))
    print((lambda x1, x2: 4*x1+x2)(4,5))

    class EmptyClass():
        def __init__(self, var):
            self.var = var

    test = EmptyClass({"potato": 1})
    var = test.var
    res = var["potato"]
    var2 = var

    test.var["potato"] = 2
    print(var)
    print(var2)
    print(res)

    print("\x1b[33mYELLOW\x1b[0m")


    import matplotlib.pyplot as plt
    import astropy.io.fits as pyfits

    image = pyfits.open(f"data/science/tde0/rings.v3.skycell.2044.052.stk.r.unconv.fits")

    for i in image:
        print(i)
        print(i.header)
        print(i.data)
    plt.imshow(image[1].data)
    plt.show()

if False:
    import matplotlib.pyplot as plt
    import numpy as np
    # example data
    x = np.arange(0.1, 4, 0.5)
    y = np.exp(-x)

    # example error bar values that vary with x-position
    error = 0.1 + 0.2 * x

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
    ax0.errorbar(x, y, yerr=error, fmt='-o')
    ax0.set_title('variable, symmetric error')

    # error bar values w/ different -/+ errors that
    # also vary with the x-position
    lower_error = 0.4 * error
    upper_error = error
    asymmetric_error = [lower_error, upper_error]

    print(asymmetric_error)
    ax1.errorbar(x, y, xerr=asymmetric_error, fmt='o')
    ax1.set_title('variable, asymmetric error')
    ax1.set_yscale('log')
    plt.show()


if False:
    import matplotlib.pyplot as plt
    import numpy as np

    import matplotlib.tri as tri

    np.random.seed(19680801)
    npts = 200
    ngridx = 200
    ngridy = 200
    x = np.random.uniform(-2, 2, npts)
    y = np.random.uniform(-2, 2, npts)
    z = x * np.exp(-x**2 - y**2)

    fig, (ax1, ax2) = plt.subplots(nrows=2)

    # -----------------------
    # Interpolation on a grid
    # -----------------------
    # A contour plot of irregularly spaced data coordinates
    # via interpolation on a grid.

    # Create grid values first.
    xi = np.linspace(-2.1, 2.1, ngridx)
    yi = np.linspace(-2.1, 2.1, ngridy)

    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    # Note that scipy.interpolate provides means to interpolate data on a grid
    # as well. The following would be an alternative to the four lines above:
    # from scipy.interpolate import griddata
    # zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear')

    ax1.contour(xi, yi, zi, levels=14, linewidths=0.5, colors='k')
    cntr1 = ax1.contourf(xi, yi, zi, levels=14, cmap="RdBu_r")

    fig.colorbar(cntr1, ax=ax1)
    ax1.plot(x, y, 'ko', ms=3)
    ax1.set(xlim=(-2, 2), ylim=(-2, 2))
    ax1.set_title('grid and contour (%d points, %d grid points)' %
                (npts, ngridx * ngridy))

    # ----------
    # Tricontour
    # ----------
    # Directly supply the unordered, irregularly spaced coordinates
    # to tricontour.

    ax2.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
    cntr2 = ax2.tricontourf(x, y, z, levels=14, cmap="RdBu_r")

    fig.colorbar(cntr2, ax=ax2)
    ax2.plot(x, y, 'ko', ms=3)
    ax2.set(xlim=(-2, 2), ylim=(-2, 2))
    ax2.set_title('tricontour (%d points)' % npts)

    plt.subplots_adjust(hspace=0.5)
    plt.show()

if False:
    import matplotlib.pyplot as plt
    import plotly.express as px
    import seaborn as sns
    import pandas as pd
    import numpy as np
    df = px.data.tips()

    #fig = px.density_contour(df, x="total_bill", y="tip")
    #fig.update_traces(contours_coloring="fill", contours_showlabels = True)
    #fig.show()

    #Import reference catalog:
    refCat = np.loadtxt("referenceCatalog.txt")
    fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
    refCat = pd.read_csv("referenceCatalog.txt", delimiter=" ", header=None, names=fieldnames)
    sns.kdeplot(refCat, x="col_60", y="col_63", fill=True, color="black")
    plt.ylabel("$M_\star$")
    plt.xlabel("n_sersic")
    plt.show()

    sns.kdeplot(refCat, x="col_67", y="col_63", fill=True, color="black")
    plt.ylabel("$M_\star$")
    plt.xlabel("$M_{BH}$")
    plt.show()

if False:
    import numpy as np
    from download_data import TDE_names
    from paper_data import TDE_stellarDensities
    import matplotlib.pyplot as plt

    tde_stellarMassesDensities = np.loadtxt("TDE_distribution.txt")[:,2]
    for i, name in enumerate(TDE_names):
        print(f"{name} : {np.around((tde_stellarMassesDensities[i]), decimals=2)} vs {TDE_stellarDensities[i][0]} -> Delta = {np.around(np.abs(tde_stellarMassesDensities[i]-TDE_stellarDensities[i][0]), decimals=2)} dex")
    plt.plot(tde_stellarMassesDensities[:-3], "o", label="Mine")
    plt.plot([i[0] for i in TDE_stellarDensities], "o", label="Graur")
    plt.ylabel(r"$\Sigma_{M_\star}$", fontsize=17)
    plt.legend()
    plt.show()

if False:
    # kde plot try to change filled contours color
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from utils import recombine_arrays
    import pandas as pd

    # load QPE and TDE data
    qpe_load_0 = np.loadtxt("QPE_allRelevantData_0.txt")
    qpe_load_1 = np.loadtxt("QPE_allRelevantData_1.txt")
    qpe_load_2 = np.loadtxt("QPE_allRelevantData_2.txt")
    QPE_fullData = recombine_arrays(qpe_load_0,qpe_load_1,qpe_load_2)

    tde_load_0 = np.loadtxt("TDE_allRelevantData_0.txt")
    tde_load_1 = np.loadtxt("TDE_allRelevantData_1.txt")
    tde_load_2 = np.loadtxt("TDE_allRelevantData_2.txt")
    TDE_fullData = recombine_arrays(tde_load_0,tde_load_1,tde_load_2)
    TDE_fullData = TDE_fullData[:,:,:]
    QPE_r50s, TDE_r50s = QPE_fullData[:,5,:], TDE_fullData[:,5,:]
    QPE_sersicIndices, TDE_sersicIndices = QPE_fullData[:,1,:], TDE_fullData[:,1,:]
    QPE_bulgeRatios, TDE_bulgeRatios = QPE_fullData[:,0,:], TDE_fullData[:,0,:]
    QPE_SMSDs, TDE_SMSDs = QPE_fullData[:,2,:], TDE_fullData[:,2,:]
    QPE_stellar_masses, TDE_stellar_masses = QPE_fullData[:,3,:], TDE_fullData[:,3,:]
    QPE_mBH, TDE_mBH = QPE_fullData[:,4,:], TDE_fullData[:,4,:]

    QPE_and_TDEs = [4,5,8]
    QPE_data = np.array([QPE_sersicIndices, QPE_bulgeRatios, QPE_SMSDs, QPE_stellar_masses, QPE_mBH])
    TDE_data = np.array([np.concatenate((TDE_sersicIndices, QPE_sersicIndices[QPE_and_TDEs])), np.concatenate((TDE_bulgeRatios, QPE_bulgeRatios[QPE_and_TDEs])), np.concatenate((TDE_SMSDs, QPE_SMSDs[QPE_and_TDEs])), np.concatenate((TDE_stellar_masses, QPE_stellar_masses[QPE_and_TDEs])), np.concatenate((TDE_mBH, QPE_mBH[QPE_and_TDEs]))])

    # load reference catalog
    refCat = np.loadtxt("referenceCatalog_modif2.txt")[:100,:]
    fieldnames = [f"col_{i}" for i in range(refCat.shape[1])]
    refCat = pd.read_csv("referenceCatalog_modif2.txt", delimiter=" ", header=None, names=fieldnames)
    columns_compare = ((60,12,68),63,67)
    markersize=9
    levels=[0.2,0.3,0.5,0.9,1]
    cmap = "Greys"

    n_mBH_ax = plt.subplot(111)
    kde = sns.kdeplot(refCat, x=f"col_{columns_compare[0][0]}", y=f"col_{columns_compare[2]}", fill=True, levels=levels, ax=n_mBH_ax, colors=["black", "red", "blue", "yellow"])
    #sns.kdeplot(referenceCatalogData, x=f"col_{columns_compare[0][0]}", y=f"col_{columns_compare[2]}", fill=False, levels=levels, linewidths=0.5, color="black", ax=n_mBH_ax)
    n_mBH_ax.errorbar(QPE_data[0,:,0], QPE_data[4,:,0], yerr=[QPE_data[4,:,1],QPE_data[4,:,2]], xerr=[QPE_data[0,:,1],QPE_data[0,:,2]], fmt="o", color="blue", markersize=markersize)
    n_mBH_ax.errorbar(TDE_data[0,:,0], TDE_data[4,:,0], yerr=[TDE_data[4,:,1],TDE_data[4,:,2]], xerr=[TDE_data[0,:,1],TDE_data[0,:,2]], fmt="*", color="red", markersize=markersize-1)
    plt.xlabel(r"$\text{Sérsic index } n$", fontsize=17)
    plt.ylabel("$\log(M_\mathrm{BH})$", fontsize=17)
    plt.show()