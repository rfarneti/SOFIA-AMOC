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
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from ocean_basins import basin_mask_ht_full
from pathlib import Path
from matplotlib.patches import Polygon
import my_style

# --- Choose basin, strength of forcing and period
OCEAN = 'Atlantic'
EXPT = ['0.1Sv','1.0Sv','0.0Sv']
DECADE = '141' #41 91 141
PRINTING = True

D_LIMS = (0, 2500)     #vertical limits

BASE_PATH = Path('/Volumes/Antwater2/DATA/')
DIR_OUT = Path('./Figures_N2_PV/')
NAMEPLOT = ['N2_anom_atl', 'N2_anom_pac', 'N2_anom_ind']


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
if OCEAN == 'Atlantic':
   nameplot = NAMEPLOT[0]
elif OCEAN == 'Pacific':
   nameplot = NAMEPLOT[1]
elif OCEAN == 'Indian':
   nameplot = NAMEPLOT[2]
else:
   print('basin Sector not present!')
   exit()
basin_mask = ms.ocean_basins.basin_mask_ht_full(OCEAN)


# Select the time period
#start_time = psi_args[expt[-1]]["start_time"]
#end_time   = psi_args[expt[-1]]["end_time"]
if DECADE == '41':
   decade = 'year_41'
   years = '(Years 41-50)'
   start_time = '1631-07-02'
   end_time   = '1640-07-02'
elif DECADE == '91':
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


# --- The Control 0.0Sv
thetao = xr.open_dataset(BASE_PATH / psi_args[EXPT[-1]]["dirin"] / 'pp' / psi_args[EXPT[-1]]["t"])['temp']
so     = xr.open_dataset(BASE_PATH / psi_args[EXPT[-1]]["dirin"] / 'pp' / psi_args[EXPT[-1]]["s"])['salt']
thetao = thetao.sel(time=slice(start_time, end_time))
so     = so.sel(time=slice(start_time, end_time))

sigma2c = xr.open_dataset(BASE_PATH / psi_args[EXPT[-1]]["dirin"] / 'pp' / psi_args[EXPT[-1]]["potrho2"])['pot_rho_2']
sigma2c = sigma2c.sel(time=slice(start_time, end_time)) * basin_mask

print('--- Computing alpha, beta, and N^2 for experiment {EXPT}'.format(EXPT=EXPT[-1]))
n2c, alpha, beta = ms.derived.calc_n2(thetao * basin_mask, so * basin_mask, eos="Wright", gravity=-9.8,
                  patm=101325.0, zcoord="st_ocean", interfaces=None,
                  adjust_negative=True)
n2c *= 1e5

# --- The experiments
labelsize = 16
titlesize = 16
cmap = cm.cm.delta
sigma_levels = [36.4,36.6,36.8,36.9]
maxn2, deltan2 = 5.0, 0.5
levels = np.arange(-maxn2, maxn2+deltan2, deltan2)

for i in range(len(EXPT)-1):
    if DECADE == '41':
       letter = "a)" if EXPT[i] == '0.1Sv' else "d)"
    elif DECADE == '91':
       letter = "b)" if EXPT[i] == '0.1Sv' else "e)"
    elif DECADE == '141': 
       letter = "c)" if EXPT[i] == '0.1Sv' else "f)"
    else:
       print('Wrong selection of decade')
       exit()      

    thetao = xr.open_dataset(BASE_PATH / psi_args[EXPT[i]]["dirin"] / 'pp' / psi_args[EXPT[i]]["t"])['temp']
    so     = xr.open_dataset(BASE_PATH / psi_args[EXPT[i]]["dirin"] / 'pp' / psi_args[EXPT[i]]["s"])['salt']
    thetao = thetao.sel(time=slice(start_time, end_time))
    so     = so.sel(time=slice(start_time, end_time))
    sigma2p = xr.open_dataset(BASE_PATH / psi_args[EXPT[i]]["dirin"] / 'pp' / psi_args[EXPT[i]]["potrho2"])['pot_rho_2']
    sigma2p = sigma2p.sel(time=slice(start_time, end_time)) * basin_mask
    
    print('--- Computing alpha, beta, and N^2 for experiment {EXPT}'.format(EXPT=EXPT[i]))
    n2p, alphap, betap = ms.derived.calc_n2(thetao * basin_mask, so * basin_mask, eos="Wright", gravity=-9.8,
                      patm=101325.0, zcoord="st_ocean", interfaces=None,
                      adjust_negative=True)
    n2p *= 1e5

    # Plot N^2 anomalies in the chosen sector
    fig = plt.figure(i+1, figsize=(8,5))
    ax1 = fig.add_subplot(1,1,1)
    
    (n2p - n2c).sel(st_ocean=slice(*D_LIMS)).mean(['xt_ocean','time']).plot.contourf(ax=ax1, 
           levels=levels, cmap=cmap, add_colorbar=False)
    #(n2p - n2c).sel(st_ocean=slice(*D_LIMS)).mean(['xt_ocean','time']).plot.contour(ax=ax1, 
    #       colors='k', levels=[-2.0,0.0,0.5,1.0,2.0], linewidths=0.5)

    # Add the sigma_2 contours
    csc=(sigma2c-1000).sel(st_ocean=slice(*D_LIMS)).mean(['xt_ocean','time']).plot.contour(ax=ax1,
          colors='black', levels=sigma_levels,
          linewidths=0.5, linestyles='-')
    csp=(sigma2p-1000).sel(st_ocean=slice(*D_LIMS)).mean(['xt_ocean','time']).plot.contour(ax=ax1,
          colors='red', levels=sigma_levels,
          linewidths=1.0, linestyles='-')
    ax1.clabel(csp, inline=False, fmt='%1.1f', rightside_up=False, fontsize=12)
 
    plt.gca().invert_yaxis()
    ax1.set_title(f'{OCEAN} N$^{2}$', size=titlesize, loc='center')
    ax1.set_title(f'{letter} {EXPT[i]}', size=titlesize, loc='left')
    ax1.set_title(f'{years}', size=titlesize, loc='right')
    ax1.set_ylabel(r'Depth [m]', fontsize=labelsize)
    ax1.set_xlabel(r'Latitude [$^\circ$]', fontsize=labelsize)
    ax1.set_xlim([-80,80])
    ax1.set_xticks(np.arange(-80,81,20))
    ax1.set_xticklabels(('80S','60S','40S','20S','0','20N','40N','60N','80N'))
    ax1.tick_params(labelsize=12)
    ax1.set_facecolor([0.7, 0.7,0.7])
    norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
    ax_cb = plt.gcf().add_axes((.92,.15,.02,.7))
    cb = matplotlib.colorbar.ColorbarBase(ax=ax_cb, cmap=cmap, norm=norm, boundaries=levels,
         extend='both', orientation='vertical', label=r'10$^{-5} [$s$^{-2}$]')
    cb.set_ticks(np.linspace(-maxn2,maxn2,11))
   
    # Add grey area for ACC channel
    #vertsACC = [(-64, 3000), (-55, 3000), (-55, 0), (-65, 0)]
    #polyACC = Polygon(vertsACC, facecolor='0.8', edgecolor='1.0')
    #ax1.add_patch(polyACC)
 
    if PRINTING:
       plt.savefig(DIR_OUT / f'N2_anom_{OCEAN}_{EXPT[i]}_{decade}.png', bbox_inches='tight',dpi=300)

plt.show()



