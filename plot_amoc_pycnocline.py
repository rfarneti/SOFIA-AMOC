
import numpy as np
import xarray as xr
import cf_xarray as cfxr
import momsofia as ms
#import pint_xarray
import cftime
import datetime
import nc_time_axis
#import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colorbar
import cmocean
from ocean_basins import basin_mask_ht_full
from pathlib import Path
import matplotlib.colors as mcolors
import scipy.stats as stats
import my_style_latex


AMOC_PATH  = Path('./DATA/MOC/')
DEPTH_PATH = Path('./DATA/')
DIR_OUT = Path('./Figures_transports/')
PRINTING = True


FORCINGS = [
         #('0.1Sv',     cmocean.cm.haline_r, 'o')
         ('1.0Sv',     cmocean.cm.haline_r, 'o')
         #('1.0Sv_60S', cmocean.cm.haline_r, 'o') 
]

LATITUDES = ['30','50']

SIGMA = 000.0
name_var = 'sigma0.0' if SIGMA == 0.0 else 'sigma2000.0'
label = r'$\Psi$ $\propto$ D^{-1/2}$'

figs = []
axes = []
for i, lat in enumerate(LATITUDES):
    fig, ax = plt.subplots(figsize=(6,5), num=i+1)
    figs.append(fig)
    axes.append(ax)

#cmap = plt.get_cmap('GnBu')
#cmap2 = cmocean.cm.haline
#cmap = [cmap1, cmap2]
for i, latitude in enumerate(LATITUDES):
    ax = axes[i]
    fig = figs[i]
    plt.figure(fig.number)

    for forcing, cmap, marker in FORCINGS:
        ds_depth = xr.open_dataset(DEPTH_PATH /f'pycnocline_depth_{name_var}_{forcing}.nc')
        depth = ds_depth['Depth']
        time = ds_depth['time']
        ds_amoc  = xr.open_dataset(AMOC_PATH /f'AMOC_index_{latitude}N_{forcing}.nc')
        amoc = ds_amoc['AMOC']

        #get a scaling 
        #C = depth[0] * amoc[0]**(1/2)
        #amoc_ext = np.arange(0,25,1)
        #scaling = C.values * (1/amoc_ext)**(1/2)

        sc = ax.scatter(depth, amoc, c=time, marker=marker, cmap=cmap, s=40, edgecolors='0.1')
        #ax.scatter(depth[-1], amoc[-1], color='k', marker=marker, s=100)
        #ax.plot(scaling, amoc_ext, color='black', linestyle='--', linewidth=1.0, label=label)

    bounds = np.linspace(0, 150.0, 16)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    cbar = plt.colorbar(sc, ticks=bounds)
    cbar.set_label('Time [Years]', fontsize=14)

    ax.set_xlim(750,900) if SIGMA == 2000.0 else ax.set_xlim(800,940)
    ax.set_yticks(np.arange(800,830,5)) if forcing == '0.1Sv' else ax.set_xticks(np.arange(750,925,25))
    ax.set_ylim(4,20)
    plt.ylabel(rf'$\Psi_z$ [Sv]', fontsize=18)
    plt.xlabel(rf'D [m]', fontsize=18)
    plt.title(rf'a)', loc='left', fontsize=16)
    plt.title(rf'AMOC vs Pycnocline Depth ({forcing})', loc='center', fontsize=16)
    #plt.title(rf'({forcing})',   loc='right', fontsize=14)
    plt.grid(alpha=0.3)
    ax.tick_params(labelsize=14)
 
    if PRINTING:
      plt.savefig(DIR_OUT / f'AMOC_vs_pycnocline_depth_{latitude}N_sigma{SIGMA}_{forcing}.png', bbox_inches='tight',dpi=300)

plt.show()


#Perform a linear regression on the log-transformed data
#amoc_arr = np.array(amoc)
#depth_arr = np.array(depth)

#log_psi = np.log(amoc_arr)
#log_d = np.log(depth_arr)
#slope, intercept, r_value, p_value, std_err = stats.linregress(log_psi, log_d)

#psi_range = np.linspace(min(amoc_arr), max(amoc_arr), 100)
#fit_line = np.exp(intercept) * (psi_range**slope)
#plt.plot(fit_line, psi_range, 'k--', label=f'Fit: $D \propto \Psi^{{{slope:.2f}}}$')

#print(f"Calculated scaling exponent: {slope:.2f}")
#print(f"Correlation (R^2): {r_value**2:.3f}")
