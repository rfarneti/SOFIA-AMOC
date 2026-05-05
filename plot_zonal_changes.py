'''

 Plot GLOBAL zonal-mean anomalies in T, S and \sigma_0
 for forcings of 0.1Sv and 1.0Sv

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
from pathlib import Path
import my_style


PRINTING = True
DIR_OUT = Path('./Figures_tracers/')
NAMEPLOT = ['zonal_anomalies_Global']

FORCINGS  = ['0.1Sv',
            '0.1Sv_tas_corrected',
            '1.0Sv',
           ]

PATH = [
        #'/Volumes/Antwater2/DATA/FAFANTWATER_ideal_runoff_pv60_rest_ice_0.0Sv/pp/',
        '/Volumes/Antwater2/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv/pp/',
        '/Volumes/Antwater2/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.1Sv/pp/',
        '/Volumes/Antwater2/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.1Sv_tas_corrected/pp/',
        '/Volumes/Antwater2/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv/pp/']

VARS  = ['thetao',     'so',         'pot_rho_0']
ANOM = ['$\Delta$ T', '$\Delta$ S', '$\Delta$ $\sigma_0$']

period = ['(Years 41-50)', '(Years 91-100)', '(Years 141-150)']
date11  = "1641-07-02"
date20  = '1650-07-02'
date51  = "1691-07-02"
date60  = '1700-07-02'
date91  = '1741-07-02'
date100 = '1750-07-02'



def plot_zonal_mean(ii, varin, anomin, levels, cmap, SO=False):
    ''' Plot Global zonal-mean for any variable'''
    axissize=18
    labelsize=12
    if varin=='thetao':
       C = xr.open_dataset(PATH[0]+varin+'.nc')['temp']
       S = xr.open_dataset(PATH[ii+1]+varin+'.nc')['temp']
    elif varin=='so':
       C = xr.open_dataset(PATH[0]+varin+'.nc')['salt']
       S = xr.open_dataset(PATH[ii+1]+varin+'.nc')['salt']
    elif varin=='pot_rho_0':
       C = xr.open_dataset(PATH[0]+varin+'.nc')['pot_rho_0']
       S = xr.open_dataset(PATH[ii+1]+varin+'.nc')['pot_rho_0']
    else:
       print('variable not present!')
       exit()
    anom11 = (S.sel(time=slice(date11,date20)).cf.mean(dim=['longitude','time'])  - 
              C.sel(time=slice(date11,date20)).cf.mean(dim=['longitude','time'])).squeeze()
    anom51 = (S.sel(time=slice(date51,date60)).cf.mean(dim=['longitude','time'])  - 
              C.sel(time=slice(date51,date60)).cf.mean(dim=['longitude','time'])).squeeze()
    anom91 = (S.sel(time=slice(date91,date100)).cf.mean(dim=['longitude','time']) - 
              C.sel(time=slice(date91,date100)).cf.mean(dim=['longitude','time'])).squeeze()

    ax1 = fig.add_subplot(3,1,1)
    anom11.plot.contourf(ax=ax1, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom11.plot.contour(ax=ax1, colors='k', levels=levels,linewidths=0.5)
    ax2 = fig.add_subplot(3,1,2)
    anom51.plot.contourf(ax=ax2, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom51.plot.contour(ax=ax2, colors='k', levels=levels,linewidths=0.5)
    ax3 = fig.add_subplot(3,1,3)
    anom91.plot.contourf(ax=ax3, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom91.plot.contour(ax=ax3, colors='k', levels=levels,linewidths=0.5)

    ax1.set_title(period[0],  loc='right', size=labelsize)
    ax2.set_title(period[1],  loc='right', size=labelsize)
    ax3.set_title(period[2],  loc='right', size=labelsize)
    ax1.set_title(f'a) {FORCINGS[ii]}', loc='center', size=axissize)
    ax2.set_title(f'b) {FORCINGS[ii]}', loc='center', size=axissize)
    ax3.set_title(f'c) {FORCINGS[ii]}', loc='center', size=axissize)
    for ax in [ax1, ax2, ax3]:
        ax.set_facecolor([0.5, 0.5,0.5])
        ax.set_ylim([5500,0])
        if SO==True:
           ax.set_xlim([-80,-30])
        else:
           ax.set_xlim([-80,80])
        ax.set_ylabel('Depth [m]', fontsize=axissize)
    ax3.set_xlabel('Latitude [$^\circ$]', fontsize=axissize)
    ax1.set_xlabel('', fontsize=axissize)
    ax1.set_xticklabels([])
    ax2.set_xlabel('', fontsize=axissize)
    ax2.set_xticklabels([])

    norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
    ax = plt.gcf().add_axes((.25,.04,.5,.01))
    cb1 = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels,
                                      orientation='horizontal', extend='both',)

    for cb in [cb1]:
        cb.ax.set_xlabel(anomin, size=14)
        cb.ax.tick_params(labelsize=8)

def save_fig(var, forcing):
    plt.savefig(DIR_OUT / f'zonal_anomalies_Global_{var}_{forcing}.png', bbox_inches='tight',dpi=300)


# --- Plotting
plt.rc('figure', figsize=(8,10))
cmap = cmocean.cm.curl

# Zonal-mean Thetao anomalies
varin = VARS[0]
anomin = ANOM[0]
for ii in range(len(PATH)-1):
    if ii==0 or ii==1:
       vmin, vmax, ci = -2,2,0.2
    else:
       vmin, vmax, ci = -5,5,0.5
    levelsT = np.arange(vmin,vmax+ci,ci)
    
    fig = plt.figure(ii+1)
    plot_zonal_mean(ii,varin, anomin, levelsT, cmap, SO=False)
    if PRINTING:
       save_fig(varin, FORCINGS[ii])

# Zonal-mean Salt anomalies
varin = VARS[1]
anomin = ANOM[1]
levelsS = [-2.,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.0,
           0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.]
for ii in range(len(PATH)-1):
    fig = plt.figure(ii+4)
    plot_zonal_mean(ii,varin, anomin, levelsS, cmap)
    if PRINTING:
       save_fig(varin, FORCINGS[ii])

## Zonal-mean Density anomalies
varin = VARS[2]
anomin = ANOM[2]
for ii in range(len(PATH)-1):
    fig = plt.figure(ii+7)
    plot_zonal_mean(ii,varin, anomin, levelsS, cmap)
    if PRINTING:
       save_fig(varin, FORCINGS[ii])


plt.show()







