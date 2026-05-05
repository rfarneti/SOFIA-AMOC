'''

 Plotting Runoff, total Runoff, ideal Runoff  

'''
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar
import matplotlib.path as mpath
from pylab import *
import xarray as xr
import cftime
import datetime
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean as cm
from pathlib import Path
import my_style_latex

PRINTING = True
DIR_OUT = Path('./Figures_sfc_fluxes/')
nameplot = 'ideal_runoff'


BASE_PATH = Path('/Volumes/Antwater2/DATA/')
control     = 'FAFANTWATER_ideal_runoff_pv60_rest_ice_0.0Sv/pp/river.nc'            # RUNOFF from CORE
runoff      = 'FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv/pp/runoff.nc'         # total runoff
irunoff     = 'FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv/pp/ideal_runoff.nc'      # ideal_runoff antwater
irunoff60   = 'FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv_60S/pp/ideal_runoff.nc'  # ideal_runoff antwater-60S
irunoffambe = 'FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv_ambe/pp/ideal_runoff.nc'  # ideal_runoff antwater-ambe

ds_c     = xr.open_dataset(BASE_PATH / f"{control}",   decode_times = False)
ds_r     = xr.open_dataset(BASE_PATH / f"{runoff}",    decode_times = False)
ds_ir    = xr.open_dataset(BASE_PATH / f"{irunoff}",   decode_times = False)
ds_ir60  = xr.open_dataset(BASE_PATH / f"{irunoff60}", decode_times = False)
ds_irambe  = xr.open_dataset(BASE_PATH / f"{irunoffambe}", decode_times = False)

LABEL = ['a) CORE-I Runoff', 'a) {\it antwater}', 'b) {\it antwater-60S}', 'c) {\it antwater-ambe}']

titlesize   = 18


def polar_plot(ii):
    """This function returns prepared axes for the polar plot"""
    #fig = plt.figure(figsize=(10,8))
    #ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax = fig.add_subplot(1,3,ii, projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
    ax.set_facecolor('lightgrey')
    ax.coastlines(); 
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.5, \
                ylocs=range(-90,-50,10), xlocs=None, \
                color='gray', alpha=0.5, linestyle='--', zorder=10)
    theta = np.linspace(0, 2*np.pi, 100)
    map_circle = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + [0.5, 0.5])
    ax.set_boundary(map_circle, transform=ax.transAxes)
    return ax #,fig


fig = plt.figure(figsize=[10,5])

vmin,vmax,ci,cmap = 0,3.8,0.1,cm.cm.ice_r
cmap.set_extremes(under='white')
levels = np.arange(vmin,vmax,ci)
scale_factor = 1e4

#ax = polar_plot(1)
#var = ds_c.river.mean(dim='time') * scale_factor 
#var.plot.contourf(ax=ax, levels=levels,extend='both', cmap=cmap, transform=ccrs.PlateCarree(), add_colorbar=False)
#ax.set_title(LABEL[0], loc='center', size=titlesize)

ax = polar_plot(1)
var = ds_ir.ideal_runoff * scale_factor
var.plot.contourf(ax=ax, levels=levels,extend='both', cmap=cmap, transform=ccrs.PlateCarree(), add_colorbar=False)
ax.set_title(LABEL[1], loc='center', size=titlesize)

ax = polar_plot(2)
var = ds_ir60.ideal_runoff * scale_factor
var.plot.contourf(ax=ax, levels=levels,extend='both', cmap=cmap, transform=ccrs.PlateCarree(), add_colorbar=False)
ax.set_title(LABEL[2], loc='center', size=titlesize)

ax = polar_plot(3)
var = ds_irambe.ideal_runoff * scale_factor
var.plot.contourf(ax=ax, levels=levels,extend='both', cmap=cmap, transform=ccrs.PlateCarree(), add_colorbar=False)
ax.set_title(LABEL[3], loc='center', size=titlesize)

norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
ax = plt.gcf().add_axes((.25,.05,.5,.01))
cb = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels,
                                      orientation='horizontal', extend='both')
cb.ax.set_title('10$^{-4}$ [kg~m$^{-2}$~s$^{-1}$]')


plt.tight_layout()
if PRINTING:
   plt.savefig(DIR_OUT / f'{nameplot}.png', bbox_inches='tight',dpi=300)


plt.show()

