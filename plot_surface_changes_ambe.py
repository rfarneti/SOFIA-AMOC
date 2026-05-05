'''

plots surface anomalies at different decades

'''
import numpy as np
import xarray as xr
import cf_xarray as cfxr
import cftime
import datetime
import nc_time_axis 
#import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colorbar
#import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
from mpl_toolkits.axes_grid1 import make_axes_locatable
import my_style
from pathlib import Path

# choose among: sst, sss, mld
VARIABLE = 'sss'

FORCINGS  = [
             ('0.1Sv_ambe','0.1Sv-ambe',0),
             ('1.0Sv_ambe','1.0Sv-ambe',1)
            ]

PRINTING = True

DIR_OUT = Path('./Figures_tracers/')

nameplot = 'surface_anomalies'

BASE_PATH = Path('/Volumes/Antwater2/DATA/')
PREFIX = Path('FAFANTWATER_ideal_runoff_Boeira_rest_ice_')

period = ['(Years 11-20)','(Years 41-50)', '(Years 91-100)', '(Years 141-150)']
date11   = "1611-07-02"
date20   = '1620-07-02'
date41   = "1641-07-02"
date50   = '1650-07-02'
date91   = "1691-07-02"
date100  = '1700-07-02'
date141  = "1741-07-02"
date150  = '1750-07-02'

VAR  = ['sos', 'tos', 'mld']
ANOM = ['$\Delta$ SSS', '$\Delta$ SST', '$\Delta$ MLD']


labelsize=9
titlesize=16
axissize=16
textsize=12
clong = -65
linewidth = 0.15

if VARIABLE == 'sss':
   print('Processing', VAR[0])
   varin = VAR[0]
   anomin = ANOM[0]
   Sc   = xr.open_dataset(BASE_PATH / f'{PREFIX}0.0Sv' / 'pp' / f'{varin}.nc')['SSS']
   vmin,vmax,ci, cmap = -1.5,1.5,.25, cmocean.cm.curl
   levels = {}
   levels[0] = np.arange(vmin,vmax+ci,ci)
   levels[1] = np.arange(vmin*2.,vmax*2+ci*2,ci*2)
   levels[2] = np.arange(vmin,vmax+ci,ci)
   cbar_kwargs = {'orientation':'horizontal', 'label':'', 
           'ticks':[-3.0,-2.5,-2.0,-1.5,-1,-0.5,0.0,0.5,1,1.5,2.0,2.5,3.0]}
   #get also the barotropic stream function
   BSF = xr.open_dataset(BASE_PATH / f'{PREFIX}0.0Sv' / 'pp' / f'quasi_bsf.nc')['quasi_bsf']
   vminBSF,vmaxBSF,ciBSF = -80,180,10
   levelsBSF  = np.arange(vminBSF,vmaxBSF+ciBSF,ciBSF)
elif VARIABLE == 'sst':
   print('Processing', VAR[1])
   varin = VAR[1]
   anomin = ANOM[1]
   Sc   = xr.open_dataset(BASE_PATH / f'{PREFIX}0.0Sv' / 'pp' / f'{varin}.nc')['SST']
   vmin,vmax,ci,ci2,cmap = -2.,2.,.2,.4,cmocean.cm.curl
   levels = {}
   levels[0] = np.arange(vmin,vmax+ci,ci)
   levels[1] = np.arange(vmin,vmax+ci2,ci2)
   levels[2] = np.arange(vmin,vmax+ci,ci)
   cbar_kwargs = {'orientation':'horizontal', 'label':'', 
           'ticks':[-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0]}
elif VARIABLE == 'mld':
   print('Processing', VAR[2])
   varin = VAR[2]
   anomin = ANOM[2]
   Sc   = xr.open_dataset(BASE_PATH / f'{PREFIX}0.0Sv' / 'pp' / f'{varin}.nc')['mld']
   vmin,vmax,ci,ci2,cmap = -180.,180.,20,20,cmocean.cm.curl
   levels = {}
   levels[0]  = np.arange(vmin,vmax+ci,ci)
   levels[1] = np.arange(vmin,vmax+ci2,ci2)
   levels[2]  = np.arange(vmin,vmax+ci,ci)
   cbar_kwargs = {'orientation':'horizontal', 'label':'', 
           'ticks':[-160,-120,-80,-40,0.0,40,80,120,160]}

#
# MAIN LOOP THROUGH FORCINGS
# Compute Anomalies
for forcing, title, ii in FORCINGS:
    if VARIABLE == 'sss':
       Sa = xr.open_dataset(BASE_PATH / f'{PREFIX}{forcing}' / 'pp' / f'{VAR[0]}.nc')['SSS']
    elif VARIABLE == 'sst':
       Sa = xr.open_dataset(BASE_PATH / f'{PREFIX}{forcing}' / 'pp' / f'{VAR[1]}.nc')['SST']
    elif VARIABLE == 'mld':
       Sa = xr.open_dataset(BASE_PATH / f'{PREFIX}{forcing}' / 'pp' / f'{VAR[2]}.nc')['mld']

    S11  = (Sa.sel(time=slice(date11, date20)).mean(dim='time')   - 
            Sc.sel(time=slice(date11, date20)).mean(dim='time')).squeeze()    
    S41  = (Sa.sel(time=slice(date41, date50)).mean(dim='time')   - 
            Sc.sel(time=slice(date41, date50)).mean(dim='time')).squeeze()    
    S91  = (Sa.sel(time=slice(date91, date100)).mean(dim='time')  - 
            Sc.sel(time=slice(date91, date100)).mean(dim='time')).squeeze()    
    S141 = (Sa.sel(time=slice(date141, date150)).mean(dim='time') - 
            Sc.sel(time=slice(date141, date150)).mean(dim='time')).squeeze()    

    fig = plt.figure(ii+1, figsize=(8,6))
    ax1 = fig.add_subplot(2,2,1, projection=ccrs.Robinson(central_longitude=clong))
    ax2 = fig.add_subplot(2,2,2, projection=ccrs.Robinson(central_longitude=clong))
    ax3 = fig.add_subplot(2,2,3, projection=ccrs.Robinson(central_longitude=clong))
    ax4 = fig.add_subplot(2,2,4, projection=ccrs.Robinson(central_longitude=clong))

    S11.plot.contourf(ax=ax1, levels=levels[ii],extend='both', 
                      cmap=cmap, transform=ccrs.PlateCarree(), add_colorbar=False)
    
    S41.plot.contourf(ax=ax2, levels=levels[ii],extend='both', 
                      cmap=cmap, transform=ccrs.PlateCarree(), add_colorbar=False)
    
    S91.plot.contourf(ax=ax3, levels=levels[ii],extend='both', 
                      cmap=cmap, transform=ccrs.PlateCarree(), add_colorbar=False)

    S141.plot.contourf(ax=ax4, levels=levels[ii],extend='both', 
                       cmap=cmap, transform=ccrs.PlateCarree(), add_colorbar=False)

    if VARIABLE == 'sss':
       # Add the climatological Barotropic StreamFunction from the Control
       BSF.plot.contour(ax=ax1, colors='k', levels=levelsBSF,linewidths=0.2, transform=ccrs.PlateCarree())
       BSF.plot.contour(ax=ax2, colors='k', levels=levelsBSF,linewidths=0.2, transform=ccrs.PlateCarree())
       BSF.plot.contour(ax=ax3, colors='k', levels=levelsBSF,linewidths=0.2, transform=ccrs.PlateCarree())
       BSF.plot.contour(ax=ax4, colors='k', levels=levelsBSF,linewidths=0.2, transform=ccrs.PlateCarree())

    ax1.set_title(period[0], loc='right', size=labelsize)
    ax2.set_title(period[1], loc='right', size=labelsize)
    ax3.set_title(period[2], loc='right', size=labelsize)
    ax4.set_title(period[3], loc='right', size=labelsize)
    ax1.set_title(f'a) {title}', loc='center', size=axissize)
    ax2.set_title(f'b) {title}', loc='center', size=axissize)
    ax3.set_title(f'c) {title}', loc='center', size=axissize)
    ax4.set_title(f'd) {title}', loc='center', size=axissize)

    for ax in [ax1, ax2, ax3, ax4]:
        #ax.set_title(forcing, loc='center', size=titlesize)
        ax.set_global()
        ax.coastlines()
        #ax.add_feature(cfeature.LAND)

    norm = matplotlib.colors.BoundaryNorm(boundaries=levels[ii], ncolors=cmap.N)
    ax = plt.gcf().add_axes((.20,.07,.6,.01))
    cb = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels[ii],
                                      extend='both', **cbar_kwargs)
    cb.ax.set_xlabel(anomin, size=textsize)

    plt.tight_layout()
    if PRINTING:
       plt.savefig(DIR_OUT / f'surface_anomalies_{VARIABLE}_{forcing}.png', bbox_inches='tight',dpi=300)



plt.show()








