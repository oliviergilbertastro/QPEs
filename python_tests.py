
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