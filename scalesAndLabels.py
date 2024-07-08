"""
Opens jpeg images of QPE and TDE hosts, adds their names and a scale bar and saves them as pdfs.
"""

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from download_data import *
import matplotlib as matt
matt.rcParams['font.family'] = 'STIXGeneral'

def scale_bar(ax, d, dist=1/0.13, text='1"', color='black', flipped=False, fontsize=15):
    if flipped:
        p0 = d - d / 15. - dist
        p1 = d / 15.
        ax.plot([p0, p0 + dist], [p1, p1], linewidth=3, color=color)
        ax.text(p0 + dist / 2., p1 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')
    else:
        p0 = d / 15.
        ax.plot([p0, p0 + dist], [p0, p0], linewidth=3, color=color)
        ax.text(p0 + dist / 2., p0 + 0.02 * d, text, fontsize=fontsize, color=color, ha='center')



fontsize = 35

for i in range(len(objects_names)):
    img = np.asarray(Image.open(r'data/images/'+f'{objects_names[i]}'+r'.jpeg'))
    ax1 = plt.subplot(111)
    frame_size=256
    fig = plt.gcf()
    fig.set_size_inches(4,4)
    plt.imshow(img, origin="lower")
    scale_bar(ax1, frame_size, dist=1/0.262, text='1"', color = 'white', fontsize=fontsize+10)
    ax1.text(frame_size*0.05, frame_size*0.87, objects_names[i],fontsize=fontsize, weight='bold', color='white')
    ax1.set_xticks([])
    ax1.set_yticks([])
    plt.subplots_adjust(0.0,0.0,1.,1.,0.,0.)
    plt.show()
    fig.savefig(f"/Users/oliviergilbert/Desktop/QPEs/figures/cutouts/qpe/{objects_names[i]}.pdf")

for i in range(len(TDE_names)):
    img = np.asarray(Image.open(r'data/images/'+f'{TDE_names[i]}'+r'.jpeg'))
    ax1 = plt.subplot(111)
    frame_size=256
    fig = plt.gcf()
    fig.set_size_inches(4,4)
    plt.imshow(img, origin="lower")
    scale_bar(ax1, frame_size, dist=1/0.262, text='1"', color = 'white', fontsize=fontsize+10)
    ax1.text(frame_size*0.05, frame_size*0.87, TDE_names[i],fontsize=fontsize, weight='bold', color='white')
    ax1.set_xticks([])
    ax1.set_yticks([])
    plt.subplots_adjust(0.0,0.0,1.,1.,0.,0.)
    plt.show()
    fig.savefig(f"/Users/oliviergilbert/Desktop/QPEs/figures/cutouts/tde/{TDE_names[i]}.pdf")