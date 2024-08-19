import numpy as np
from matplotlib import pyplot as plt

from astroML.density_estimation import XDGMM
from astroML.plotting.tools import draw_ellipse



if "example" == True:
    # Sample the dataset. 
    # Here we use sample size = 400 in the example, 
    # which converges in shorter time, and gives reasonable result.
    N = 400
    np.random.seed(0)

    # generate the true data
    x_true = (1.4 + 2 * np.random.random(N)) ** 2
    y_true = 0.1 * x_true ** 2

    # add scatter to "true" distribution
    dx = 0.1 + 4. / x_true ** 2
    dy = 0.1 + 10. / x_true ** 2

    x_true += np.random.normal(0, dx, N)
    y_true += np.random.normal(0, dy, N)

    # define a function to plot all distributions in the same format
    def plot_distribution(text, sample_x, sample_y):
        plt.figure(figsize=(5, 3.75))
        plt.scatter(sample_x, sample_y, s=4,lw=0,c='k')
        plt.xlim(-1, 13)
        plt.ylim(-6, 16)
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.title(text,fontsize=10)

    # plot true distribution
    plot_distribution('True Distribution', x_true, y_true)

    plt.show()



    # add noise to get the "observed" distribution
    dx = 0.2 + 0.5 * np.random.random(N)
    dy = 0.2 + 0.5 * np.random.random(N)

    x = x_true + np.random.normal(0, dx)
    y = y_true + np.random.normal(0, dy)

    # plot noisy distribution
    plot_distribution('Noisy Distribution', x, y)

    plt.show()

    # stack the results for computation
    X = np.vstack([x, y]).T
    Xerr = np.zeros(X.shape + X.shape[-1:])
    diag = np.arange(X.shape[-1])
    Xerr[:, diag, diag] = np.vstack([dx ** 2, dy ** 2]).T

    clf = XDGMM(n_components=10, max_iter=200)

    clf.fit(X, Xerr)
    sample = clf.sample(N)

    # plot noisy distribution
    plot_distribution('Extreme Deconvolution Resampling', sample[:, 0], sample[:, 1])
    plt.show()

    # Plot the results
    fig = plt.figure(figsize=(5, 3.75))
    fig.subplots_adjust(left=0.1, right=0.95,
                        bottom=0.1, top=0.95,
                        wspace=0.02, hspace=0.02)

    ax1 = fig.add_subplot(221)
    ax1.scatter(x_true, y_true, s=4, lw=0, c='k')

    ax2 = fig.add_subplot(222)
    ax2.scatter(x, y, s=4, lw=0, c='k')

    ax3 = fig.add_subplot(223)
    ax3.scatter(sample[:, 0], sample[:, 1], s=4, lw=0, c='k')



    titles = ["True Distribution", "Noisy Distribution",
            "Extreme Deconvolution\n  resampling",
            "Extreme Deconvolution\n  cluster locations"]

    ax = [ax1, ax2, ax3]

    for i in range(len(ax)):
        ax[i].set_xlim(-1, 13)
        ax[i].set_ylim(-6, 16)

        ax[i].xaxis.set_major_locator(plt.MultipleLocator(4))
        ax[i].yaxis.set_major_locator(plt.MultipleLocator(5))

        ax[i].text(0.05, 0.95, titles[i],
                ha='left', va='top', transform=ax[i].transAxes)

        if i in (0, 1):
            ax[i].xaxis.set_major_formatter(plt.NullFormatter())
        else:
            ax[i].set_xlabel('$x$')

        if i in (1, 3):
            ax[i].yaxis.set_major_formatter(plt.NullFormatter())
        else:
            ax[i].set_ylabel('$y$')

    plt.show()