import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_utr_heatmap(data, title=None, cmap=cm.OrRd, orientation='vertical'):
    utr = data.copy()
    max_val = np.max(utr)
    min_val = np.min(utr)
    midpoint = (max_val + min_val) / 2
    fig, ax = plt.subplots(figsize=(18,18))
    im = ax.imshow(utr,cmap=cmap)
    if title:
        ax.set_title(title, fontsize=20)
    else:
        ax.set_title("Read depth at 5'UTR", fontsize=20)
    ax.set_xlabel('Position from TSS', fontsize=20)
    ax.get_yaxis().set_ticks([])
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, ticks=[min_val,0,max_val], orientation=orientation)
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(labelsize=20)
    #cbar.ax.set_xticklabels(['0', '{0}'.format(int(midpoint)), '{0}'.format(int(max_val))], fontsize=20)
    plt.show()


