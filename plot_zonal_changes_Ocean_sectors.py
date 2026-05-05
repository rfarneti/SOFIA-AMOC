'''

 Plot zonal-mean anomalies in T, S, rho_0 for different Oceans

 One can choose 
 the basin: Atlantic, Pacific, Indian
 the field: so, thetao, pot_rho_0
 the forcing: 0.1Sv, ..., 1.0Sv

'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
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
from pathlib import Path
import my_style


##====================================
## Choose basin, variable, strength of forcing
OCEAN    = 'Atlantic' # Atlantic, Pacific
VARIABLE = 'thetao' # thetao, so, pot_rho_0, pot_rho_2
FORCING  = '0.1Sv'  # 0.1Sv, 1.0Sv

PRINTING = True
##====================================

BASE_PATH = Path('/Volumes/Antwater2/DATA/')
PREFIX = 'FAFANTWATER_ideal_runoff_Boeira_rest_ice_'

DIR_OUT = Path('./Figures_tracers/')
NAMEPLOTS = ['zonal_anomalies_Atlantic',
             'zonal_anomalies_Pacific',
             'zonal_anomalies_Indian']

#VARFILE = ['thetao',     'so',         'pot_rho_0',           'pot_rho_2']
ANOM    = ['$\Delta$ T', '$\Delta$ S', '$\Delta$ $\sigma_0$', '$\Delta$ $\sigma_2$']

period = ['(Years 41-50)', '(Years 91-100)', '(Years 141-150)']
date41  = "1641-07-02"
date50  = '1650-07-02'
date91  = "1691-07-02"
date100 = '1700-07-02'
date141 = '1741-07-02'
date150 = '1750-07-02'


if OCEAN=='Atlantic':
   nameplot=NAMEPLOTS[0]
elif OCEAN=='Pacific':
   nameplot=NAMEPLOTS[1]
elif OCEAN=='Indian':
   nameplot=NAMEPLOTS[2]
else:
   print('basin Sector not present!')
   exit() 

if VARIABLE=='thetao':
   anomin    = ANOM[0]
elif VARIABLE=='so':
   anomin    = ANOM[1]
elif VARIABLE=='pot_rho_0':
   anomin    = ANOM[2]
elif VARIABLE=='pot_rho_2':
   anomin    = ANOM[3]
else:
   print('variable not present!')
   exit()

print(f'Variable is = {VARIABLE}')
print(f'Ocean is = {OCEAN}')
print(f'Forcing is = {FORCING}')



# Generate Basin Mask for the chosen sector
basin_mask = basin_mask_ht_full(OCEAN)


# --- Plotting
axissize=18
labelsize=14

fig = plt.figure(1, figsize=(8,10))
ax1 = fig.add_subplot(3,1,1)
ax2 = fig.add_subplot(3,1,2)
ax3 = fig.add_subplot(3,1,3)

if VARIABLE == 'so' or VARIABLE == 'pot_rho_0' or VARIABLE == 'pot_rho_2':
   vmin, vmax, ci, cmap = -2.0,2.0,0.1, cmocean.cm.curl
elif VARIABLE == 'thetao':
     #vmin, vmax, ci, cmap = -5.0,5.0,0.5, cmocean.cm.curl
     vmin, vmax, ci, cmap = -1.5,1.5,0.25, cmocean.cm.curl
levels = np.arange(vmin,vmax+ci,ci)

forcing0 = '0.0Sv'
if VARIABLE=='so':
   S  = xr.open_dataset(BASE_PATH / f"{PREFIX}{FORCING}" / 'pp' / f'{VARIABLE}.nc')['salt']
   Sc = xr.open_dataset(BASE_PATH / f"{PREFIX}{forcing0}" / 'pp' / f'{VARIABLE}.nc')['salt']
elif VARIABLE=='thetao':
   S  = xr.open_dataset(BASE_PATH / f"{PREFIX}{FORCING}" / 'pp' / f'{VARIABLE}.nc')['temp']
   Sc = xr.open_dataset(BASE_PATH / f"{PREFIX}{forcing0}" / 'pp' / f'{VARIABLE}.nc')['temp']
elif VARIABLE=='pot_rho_0':
   S  = xr.open_dataset(BASE_PATH / f"{PREFIX}{FORCING}" / 'pp' / f'{VARIABLE}.nc')['pot_rho_0']
   Sc = xr.open_dataset(BASE_PATH / f"{PREFIX}{forcing0}" / 'pp' / f'{VARIABLE}.nc')['pot_rho_0']
elif VARIABLE=='pot_rho_2':
   S  = xr.open_dataset(BASE_PATH / f"{PREFIX}{FORCING}" / 'pp' / f'{VARIABLE}.nc')['pot_rho_2']
   Sc = xr.open_dataset(BASE_PATH / f"{PREFIX}{forcing0}" / 'pp' / f'{VARIABLE}.nc')['pot_rho_2']
else:
   print("variable not available!")
   exit()

anom41 =  ((S.sel(time=slice(date41,date50))   * basin_mask).cf.mean(dim=['time','longitude']) - 
           (Sc.sel(time=slice(date41,date50))  * basin_mask).cf.mean(dim=['time','longitude'])).squeeze()
anom91 =  ((S.sel(time=slice(date91,date100))  * basin_mask).cf.mean(dim=['time','longitude']) - 
           (Sc.sel(time=slice(date91,date100)) * basin_mask).cf.mean(dim=['time','longitude'])).squeeze()
anom141 = ((S.sel(time=slice(date141,date150)) * basin_mask).cf.mean(dim=['time','longitude']) - 
           (Sc.sel(time=slice(date141,date150))* basin_mask).cf.mean(dim=['time','longitude'])).squeeze()

anom41.plot.contourf(ax=ax1, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
anom41.plot.contour(ax=ax1, colors='k', levels=levels,linewidths=0.5)
anom91.plot.contourf(ax=ax2, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
anom91.plot.contour(ax=ax2, colors='k', levels=levels,linewidths=0.5)
anom141.plot.contourf(ax=ax3, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
anom141.plot.contour(ax=ax3, colors='k', levels=levels,linewidths=0.5)

ax1.set_title(period[0], loc='right', size=labelsize)
ax2.set_title(period[1], loc='right', size=labelsize)
ax3.set_title(period[2], loc='right', size=labelsize)
if FORCING == '0.1Sv':
   ax1.set_title(f'a) {FORCING}', loc='center', size=axissize)
   ax2.set_title(f'b) {FORCING}', loc='center', size=axissize)
   ax3.set_title(f'c) {FORCING}', loc='center', size=axissize)
else:
   ax1.set_title(f'd) {FORCING}', loc='center', size=axissize)
   ax2.set_title(f'e) {FORCING}', loc='center', size=axissize)
   ax3.set_title(f'f) {FORCING}', loc='center', size=axissize)

for ax in [ax1, ax2, ax3]:
    ax.set_facecolor([0.5, 0.5,0.5])
    ax.set_ylim([5500,0])
    ax.set_xlim([-80,80])
    ax.set_xticks(np.arange(-80,81,20))
    ax.set_xticklabels(('80S','60S','40S','20S','0','20N','40N','60N','80N'))
    ##ax.set_title(FORCING, loc='center', size=axissize)
    if FORCING == '0.1Sv':
       ax.set_ylabel('Depth [m]', fontsize=axissize)
    else:
       ax.set_ylabel('Depth [m]', fontsize=axissize)
    if ax==ax3:
       ax3.set_xlabel('Latitude', fontsize=axissize)
    else:
       ax.set_xlabel('')
       ax.set_xticklabels([])

norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
ax = plt.gcf().add_axes((.3,.02,.4,.01))
cb1 = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels,
                                      orientation='horizontal', extend='both',)

for cb in [cb1]:
    cb.ax.set_title(anomin, size=14)
    cb.ax.tick_params(labelsize=8)

if PRINTING:
   plt.savefig(DIR_OUT / f'{nameplot}_{VARIABLE}_{FORCING}.png', bbox_inches='tight',dpi=300)
 

plt.show()


