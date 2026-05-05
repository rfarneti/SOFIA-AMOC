""" Computes the maximum of the streamfunction between [lat_min ; lat_max] 
    and plots a hovmoller as a function of depth
"""
import numpy as np
import xarray as xr
import cf_xarray as cfxr
import cftime
import datetime
import nc_time_axis 
#import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colorbar
import cmocean
import cartopy.feature as cfeature
from pathlib import Path
import my_style_latex

PRINTING = True
DIR_OUT = Path('./Figures_transports/')
BASE_PATH = Path('/Volumes/Antwater2/DATA/MOC/')

EXPERIMENTS = [
   ("0.1Sv", "a)", "(0.1Sv)"),
   ("1.0Sv", "b)", "(1.0Sv)")
]
 
OFFSET_YEAR = 1600
DAYS_IN_YEAR = 365

labelsize = 12
titlesize = 12

LATITUDE = ['-32','30','50']
start_time,end_time,dt = 0, 150, 10

# --- Get the control AMOC first
#amoc_c = xr.open_dataset(BASE_PATH / f"MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv.nc", decode_times=False)['AMOC']
amoc_c = xr.open_dataset(BASE_PATH / f"MOC_FAFANTWATER_ideal_runoff_pv60_rest_ice_0.0Sv.nc", decode_times=False)['AMOC']


# --- Plotting
for i, lat in enumerate(LATITUDE):
    latt = int(lat)
    hem = 'N' if latt >= 0 else 'S'
    fig, axes = plt.subplots(2,1, figsize = (6,6), sharex=True, constrained_layout=True)
    vmin, vmax, ci, cmap = -20.0,20.0,2.0, cmocean.cm.delta
    levels = np.arange(vmin,vmax+ci,ci)
    
    for ii, (exp_id, label, name) in enumerate(EXPERIMENTS):
        ax = axes[ii]
        file_path = BASE_PATH / f"MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_{exp_id}.nc" 
        amoc = xr.open_dataset(file_path, decode_times=False)['AMOC']
        time_years = amoc.TIME / DAYS_IN_YEAR - OFFSET_YEAR
        amoc = amoc.assign_coords(TIME=time_years)
        cf = amoc.sel(YU_OCEAN=lat).plot.contourf(ax=ax, x='TIME',
                      yincrease=False, cmap=cmap, levels=levels, extend='both', add_colorbar=False)
        amoc.sel(YU_OCEAN=lat).plot.contour(ax=ax, x='TIME', 
                 levels=[0.0], colors='chocolate', linewidths=0.5)
        ax.set_title('')
        ax.set_title(f"{name}", loc='right', size=titlesize)
        ax.set_title(f"{label} $\Psi_z(y=${abs(latt)}$^\circ${hem}$,z,t)$", loc='left', size=titlesize)
        ax.set_ylim([5000,0])
        ax.set_xlim([start_time,end_time])
        ax.set_xticks(np.arange(start_time,end_time+dt,dt))
        ax.set_yticks(np.arange(0,5000,1000))
        ax.set_ylabel('Depth [m]', fontsize=labelsize)
        ax.tick_params(axis='both', which='major', labelsize=10)
        if ax==axes[-1]:
           ax.set_xlabel('Time [years]', fontsize=labelsize)
           ax.set_xticklabels(np.arange(start_time,end_time+dt,dt))
        else:
           ax.set_xlabel('')
           ax.set_xticklabels([])
    cb = fig.colorbar(cf, ax=axes[-1], orientation='horizontal',
                      fraction=0.06, pad=0.05, shrink=0.6, aspect=30)
    cb.set_label('[Sv]', fontsize=10)

    if PRINTING:
       fig.savefig(DIR_OUT / f'AMOC_hovmoller_{abs(latt)}{hem}.png', bbox_inches='tight',dpi=300)


plt.show()





