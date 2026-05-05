'''
   Compute the buoyancy frequency 

   N^2 = g * (alpha * dT/dz - beta dS/dz) 

   - calculate alpha
   - calculate beta
   - calculate N2

'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
import momsofia as ms
#import cftime
#import datetime
#import nc_time_axis
#import netCDF4 as nc
import cmocean as cm
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colorbar
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from ocean_basins import basin_mask_ht_full
from pathlib import Path
import my_style

# --- Choose basin, strength of forcing and period
OCEAN = 'Atlantic'
EXPT = ['0.0Sv','0.1Sv','1.0Sv','1.0Sv_60S']
DECADE = '91'
PRINTING = True


BASE_PATH = Path('/Volumes/Antwater2/DATA/')
DIR_OUT = Path('./Figures_N2_PV/')


# Create dictionaries for all experiments
psi_args = {
    "0.0Sv": {
#        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv",
        "dirin": "FAFANTWATER_ideal_runoff_pv60_rest_ice_0.0Sv",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
    },
    "0.1Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.1Sv",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
    },
    "1.0Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
    },
    "1.0Sv_60S": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv_60S",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
    },
}


# --- Generate Basin Mask for the chosen sector
basin_mask = ms.ocean_basins.basin_mask_ht_full(OCEAN)


# Select the time period
#start_time = psi_args[expt[-1]]["start_time"]
#end_time   = psi_args[expt[-1]]["end_time"]
if DECADE == '91':
   decade = 'year_91'
   years = '(Years 91-100)'
   start_time = '1691-07-02'
   end_time   = '1700-07-02'
elif DECADE == '141':
   decade = 'year_141'
   years = '(Years 141-150)'
   decade = 'year_141'
   start_time = '1741-07-02'
   end_time   = '1750-07-02'
else:
   print('Wrong selection of decade')
   exit()



# --- The experiments
for i in range(len(EXPT)):
    thetao = xr.open_dataset(BASE_PATH / psi_args[EXPT[i]]["dirin"] / 'pp' / psi_args[EXPT[i]]["t"])['temp']
    so     = xr.open_dataset(BASE_PATH / psi_args[EXPT[i]]["dirin"] / 'pp' / psi_args[EXPT[i]]["s"])['salt']
    thetao = thetao.sel(time=slice(start_time, end_time))
    so     = so.sel(time=slice(start_time, end_time))
    print('--- Computing alpha, beta, and N^2 for experiment {EXPT}'.format(EXPT=EXPT[i]))
    n2, alpha, beta = ms.derived.calc_n2(thetao * basin_mask, so * basin_mask, eos="Wright", gravity=-9.8,
                      patm=101325.0, zcoord="st_ocean", interfaces=None,
                      adjust_negative=True)

    # Plot N^2 in the chosen sector
    labelsize = 14
    titlesize = 14
    cmap = cm.cm.deep
    fig = plt.figure(i+1, figsize=(8,5))
    ax1 = fig.add_subplot(1,1,1)
    
    maxn2, delta = 3.0, 0.2
    levels = np.arange(0.0, maxn2+delta, delta)
    n2 *= 1e4
    n2.sel(st_ocean=slice(0,2200)).mean(['xt_ocean','time']).plot.contourf(ax=ax1, 
           levels=levels, cmap=cmap, add_colorbar=False)
    n2.sel(st_ocean=slice(0,2200)).mean(['xt_ocean','time']).plot.contour(ax=ax1, 
           colors='k', levels=[0.0,0.2,0.4,0.6,0.8,1.0,2.0,3.0], linewidths=0.5)
    plt.gca().invert_yaxis()
    ax1.set_title(f'{OCEAN} N$^{2}$', size=titlesize, loc='left')
    ax1.set_title(f'{EXPT[i]}', size=titlesize, loc='center')
    ax1.set_title(f'{years}', size=titlesize, loc='right')
    ax1.set_ylabel(r'Depth [m]', fontsize=labelsize)
    ax1.set_xlabel(r'Latitude [$^\circ$]', fontsize=labelsize)
    ax1.set_facecolor([0.7, 0.7,0.7])
    norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
    ax_cb = plt.gcf().add_axes((.92,.2,.02,.6))
    cb = matplotlib.colorbar.ColorbarBase(ax=ax_cb, cmap=cmap, norm=norm, boundaries=levels,
         extend='max', orientation='vertical', label=r'10$^{-4}~[$s$^{-2}$]')
    cb.set_ticks(np.linspace(0.0,maxn2,11))
    
    if PRINTING:
       plt.savefig(DIR_OUT / f'N2_{OCEAN}_{EXPT[i]}_{decade}.png', bbox_inches='tight',dpi=300)

plt.show()



