"""
	Plots the anomalous Global and Atlantic MOC, in depth space,
	for different anomalous forcing
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
import my_style


PRINTING = True
DIR_OUT = Path('./Figures_transports/')
BASE_PATH = Path('/Volumes/Antwater2/DATA/MOC/')

EXPERIMENTS = [
   ("0.1Sv",     "a)", "0.1Sv"),
   ("1.0Sv",     "b)", "1.0Sv"),
   ("1.0Sv_60S", "c)", "1.0Sv_60S"),
]
 
START_YEAR = '1741' ; END_YEAR = '1750' 

OFFSET_YEAR = 1600
DAYS_IN_YEAR = 365

t_start, t_end = f"{START_YEAR}-07-02", f"{END_YEAR}-07-02"

titlesize = 18
labelsize = 14

vmin, vmax, ci, cmap = -25,25,2.5, cmocean.cm.delta
levels = np.arange(vmin,vmax+ci,ci)
dvmin, dvmax, dci = -12,13,1
dlevels = np.arange(dvmin,dvmax,dci)


# --- Compute time-mean from the Control
ds_ctrl  = xr.open_dataset(BASE_PATH / f'MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv.nc')
MOCc  = ds_ctrl['MOC'].sel(TIME=slice(t_start,t_end)).cf.mean(dim='TIME') 
AMOCc = ds_ctrl['AMOC'].sel(TIME=slice(t_start,t_end)).cf.mean(dim='TIME')


# --- GMOC
fig, axes = plt.subplots(3, 1, figsize=(8,10), sharex=True, constrained_layout=True)
for ii, (exp_id, label, name) in enumerate(EXPERIMENTS):
    ax=axes[ii]

    ds   = xr.open_dataset(BASE_PATH / f"MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_{exp_id}.nc")
    MOCa  = ds['MOC'].sel(TIME=slice(t_start,t_end)).cf.mean(dim='TIME')
    diff  = (MOCa - MOCc).squeeze()

    cf = diff.plot.contourf(ax=ax, levels=dlevels,extend='both', cmap=cmap, add_colorbar=False)
    MOCc.plot.contour(ax=ax, levels=levels, colors='k', linewidths=0.5, alpha=0.5)
    
    ax.set_title(f'{label} $\Delta\Psi$', loc='left', size=titlesize);
    ax.set_title(f'({name})', loc='right', size=labelsize);
    ax.set_facecolor([0.7, 0.7,0.7])
    ax.set_yticks([1000, 2000, 3000, 4000, 5000])
    ax.set_ylabel('Depth [m]', fontsize=labelsize)
    ax.invert_yaxis()
    ax.set_xticks(np.arange(-90,91,30))
    ax.tick_params(axis='both', which='major', labelsize=labelsize)
    if ax==axes[-1]:
       ax.set_xlabel('Latitude [$^{\circ}$]', fontsize=labelsize)
       ax.set_xticklabels(('90S','60S','30S','0','30N','60N','90N'))
    else:
       ax.set_xlabel('')
       ax.set_xticklabels([])
cb = fig.colorbar(cf, ax=axes[-1], orientation='horizontal', 
                  fraction=0.06, pad=0.05, shrink=0.6, aspect=30)
cb.set_label('[Sv]', fontsize=labelsize)

if PRINTING:       
   fig.savefig(DIR_OUT / f'GMOC_anomalies.png',bbox_inches='tight',dpi=300)


# --- AMOC
fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True, constrained_layout=True)
for ii, (exp_id, label, name) in enumerate(EXPERIMENTS):
    ax=axes[ii]
    
    ds = xr.open_dataset(BASE_PATH / f"MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_{exp_id}.nc")
    AMOCa  = ds['AMOC'].sel(TIME=slice(t_start,t_end)).mean(dim='TIME')
    diff  = (AMOCa - AMOCc).squeeze()

    cf = diff.plot.contourf(ax=ax, levels=dlevels,extend='both', cmap=cmap, add_colorbar=False)
    AMOCc.plot.contour(ax=ax, levels=levels, colors='k', linewidths=0.5, alpha=0.5)
    ax.set_title(f'{label} $\Delta\Psi$', loc='left', size=titlesize);
    ax.set_title(f'({name})', loc='right', size=labelsize);
    ax.set_facecolor([0.7, 0.7,0.7])
    ax.set_yticks([1000, 2000, 3000, 4000, 5000])
    ax.set_ylabel('Depth [m]', fontsize=labelsize)
    ax.invert_yaxis()
    ax.set_xlim([-32,90]) 
    ax.set_xticks(np.arange(-30,91,30))
    ax.tick_params(axis='both', which='major', labelsize=labelsize)
    if ax==axes[-1]:
       ax.set_xlabel('Latitude [$^{\circ}$]', fontsize=labelsize)
       ax.set_xticklabels(('30S','0','30N','60N','90N'))
    else:
       ax.set_xlabel('')
       ax.set_xticklabels([])
cb = fig.colorbar(cf, ax=axes[-1], orientation='horizontal', 
                  fraction=0.06, pad=0.05, shrink=0.6, aspect=30)
cb.set_label('[Sv]', fontsize=labelsize)

if PRINTING:
   fig.savefig(DIR_OUT / f'AMOC_anomalies.png', bbox_inches='tight',dpi=300)


plt.show()





