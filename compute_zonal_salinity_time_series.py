'''

 Plot Atlantic zonal-mean pot_rho_0 for different forcings

 One can choose 
 the basin: Atlantic, Pacific, Indian
 the field: so, thetao, pot_rho_0
 the forcing: 0.1Sv, ..., 1.0Sv

'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
import momsofia as ms
#import cf_xarray.units
#import pint_xarray
import cftime
import datetime
import nc_time_axis 
#import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colorbar
import cmocean
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ocean_basins import basin_mask_ht_full

printing = True

## Choose basin, variable, strength of forcing
ocean    = 'Atlantic'
forcing  = ["0.0Sv", "0.1Sv", "1.0Sv", "1.0Sv_60S"]

## Set path, plot name, latitude limits 
path = '/Volumes/Antwater/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_'

dirout = './Figures_Boeira/'
nameplot = ['zonal_mean_Atlantic_rho2_time_series']

## limits for \sigma_n and \sigma_s
'''         30S-60N ; 40N-60N  Bonan et al. 2022
            34S-40N ; 40N-60N  vanWesten et al. 2025
        50S-30S     ; 50N-80N  Swingedouw et al. 2009
        56S-34S     ; 43N-65N  Vanderborght et al. 2025 
'''
sn1 =  50; sn2 = 80
ss1 = -30; ss2 = 50

axissize  = 14
#labelsize = 8
titlesize= 14
year = 365; offset = 1600

##=========================================================================
## Generate Basin Mask for the chosen sector
if ocean=='Atlantic':
   nameplot=nameplot[0]
elif ocean=='Pacific':
   nameplot=nameplot[1]
elif ocean=='Indian':
   nameplot=nameplot[2]
else:
   print('basin Sector not present!')
   exit() 
print('Ocean is =',    ocean)
basin_mask = ms.ocean_basins.basin_mask_ht_full(ocean, check=False)


## Load the area of the cell
area = xr.open_dataset(path+forcing[0]+'/pp/'+'areacello_t.nc')['area_t']


## Compute the volume-weighted zonal-mean
def zonal_mean_salt(ii, area, basin_mask):
    var = xr.open_dataset(path+forcing[ii]+'/pp/'+'so'+'.nc', decode_times = False)['salt']
    dzt = xr.open_dataset(path+forcing[ii]+'/pp/'+'dht.nc', decode_times = False)['dht']
    volume = ms.derived.calc_volume(area, dzt)
    volume = volume * basin_mask
    var  = var  * basin_mask
    var_z  = (var * volume).cf.mean(dim='longitude') / volume.cf.mean('longitude')
    return var_z


for ii in range(len(forcing)):
    fig = plt.figure(ii+1, figsize=(8,6))
    ax1 = fig.add_subplot(1,1,1)
    so_z = zonal_mean_salt(ii, area, basin_mask)
    so_n = so_z.sel(yt_ocean=(slice(sn1,sn2))).cf.mean(dim=['latitude','vertical'])
    so_s = so_z.sel(yt_ocean=(slice(ss1,ss2))).cf.mean(dim=['latitude','vertical'])
    so_n['time'] = so_n.time/year - offset
    so_s['time'] = so_s.time/year - offset

    ax1.plot(so_s.time, so_s - so_s[0], 
             color='k', linestyle='--', linewidth=2.0, label=r'$\sigma_b$')
    ax1.plot(so_n.time, so_n - so_n[0], 
             color='k', linestyle=':', linewidth=2.0, label=r'$\sigma_n$')
    ax1.plot(so_s.time, (so_n - so_n[0]) - (so_s - so_s[0]), 
             color='k', linestyle='-', linewidth=2.0, label=r'$\sigma_n$-$\sigma_b$')

    ax1.legend(loc='upper left', ncol=1, fontsize=titlesize, frameon=False)
    for ax in [ax1]:
        ax.set_title(forcing[ii], fontsize=titlesize, loc='right')
        ax.set_ylabel('[kg m$^{-3}]$', fontsize=titlesize)
        ax.set_xlabel("Time [Years]",fontsize=titlesize)
        ax.set_yticks(np.arange(-1.2,1.6,0.40), fontsize=axissize)
        ax.set_ylim([-1.2,1.2])
        ax.set_xlim([0,150])
        ax.set_xticks(np.arange(0,160,20), fontsize=axissize)
        ax.tick_params(labelsize=axissize)
        ax.axhline(y=0.0, linewidth=0.5, color='grey', linestyle="solid")
        #ax.grid()
    if ii==0:
      ax1.set_title(r'$\bf{a}$', loc='left', size=axissize)
    if ii==1:
      ax1.set_title(r'$\bf{b}$', loc='left', size=axissize)
    if ii==2:
      ax1.set_title(r'$\bf{c}$', loc='left', size=axissize)
    if ii==3:
      ax1.set_title(r'$\bf{d}$', loc='left', size=axissize)

    if printing:
       plt.savefig('{dirout}{nameplot}_{forcing}.png'.format(dirout=dirout, 
                   nameplot=nameplot, forcing=forcing[ii]),
                   format='png', transparent = True, bbox_inches='tight',dpi=300)
 

plt.show()


