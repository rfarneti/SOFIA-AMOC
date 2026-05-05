'''

 Plot zonal-mean anomalies in and S and pot_rho_0
 using a Basin mask (not a sector Ocean all the way to Antarctica)

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

#---------------------------------------------------------------------------------

printing = False
dirout = './Figures_boeira/'
nameplot = ['zonal_anomalies_Atlantic','zonal_anomalies_IndoPacific']


path = [
       #'/Volumes/Antwater/DATA/FAFANTWATER_ideal_runoff_pv60_rest_ice_0.0Sv/pp/',
        '/Volumes/Antwater/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv/pp/',
        '/Volumes/Antwater/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.1Sv/pp/',
        '/Volumes/Antwater/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv/pp/']

var = ['thetao','so', 'pot_rho_0']
field = ['$\Delta$ T', '$\Delta$ S', '$\Delta$ $\sigma_0$']

label  = ['a) 0.1 Sv','b) 1.0 Sv']
label2 = ['c) 0.1 Sv','d) 1.0 Sv']
label3 = ['e) 0.1 Sv','f) 1.0 Sv']
label4 = ['g) 0.1 Sv','h) 1.0 Sv']

period = ['(Years 41-50)', '(Years 91-100)', '(Years 141-150)']
date41  = "1641-07-02"
date50  = '1650-07-02'
date91  = "1691-07-02"
date100 = '1700-07-02'
date141 = '1741-07-02'
date150 = '1750-07-02'


# Control variables, both Restoring and Boeira
ctlTb   = xr.open_dataset(path[0]+var[0]+'.nc')['temp']
ctlTb41 = ctlTb.sel(time=slice(date41,date50)).mean(dim='time').compute()
ctlTb91 = ctlTb.sel(time=slice(date91,date100)).mean(dim='time').compute()
ctlTb141 = ctlTb.sel(time=slice(date141,date150)).mean(dim='time').compute()

ctlSb = xr.open_dataset(path[0]+var[1]+'.nc')['salt']
ctlSb41  = ctlSb.sel(time=slice(date41,date50)).mean(dim='time').compute()
ctlSb91  = ctlSb.sel(time=slice(date91,date100)).mean(dim='time').compute()
ctlSb141  = ctlSb.sel(time=slice(date141,date150)).mean(dim='time').compute()

ctlDb = xr.open_dataset(path[0]+var[2]+'.nc')['pot_rho_0']
ctlDb41  = ctlDb.sel(time=slice(date41,date50)).mean(dim='time').compute()
ctlDb91  = ctlDb.sel(time=slice(date91,date100)).mean(dim='time').compute()
ctlDb141  = ctlDb.sel(time=slice(date141,date150)).mean(dim='time').compute()



axissize=12
labelsize=8




# Now do the zonal-mean in the Atlantic basin only
#=================================================
def basin_masks(ocean):
    """ Returns mask for Atlantic basin or IndoPacific basin"""
    land_mask = xr.open_dataset('/Volumes/Antwater/DATA/grids/regionmask_v6.nc')
    basin = land_mask.variables['tmask'][:]
    mask = np.zeros(basin.shape)
    print('Ocean =', ocean)
    if ocean=='Atlantic':
       #mask[(basin==2) | (basin==4) | (basin==6) | (basin==7) | (basin==8) | (basin==9) ] = 1
       mask[(basin==2) | (basin==4) ] = 1
    elif ocean=='IndoPacific':
       mask[(basin==3) | (basin==5) ] = 1
    elif ocean=='Southern':
       so_mask[(basin==1) ] = 1
    else:
        print("mask not available!") 
        exit()  
    vmask = 1*mask#; vmask[:-1,:] = vmask[:-1,:] * vmask[1:,:]
    vmask_1 = np.ma.array( vmask, mask=vmask==0)
    return vmask_1


ocean = 'Atlantic'
#ocean = 'IndoPacific'
basin_mask = basin_masks(ocean)

#Salinity
fig = plt.figure(1, figsize=(8,10))
varin = var[1]
fieldin = field[1]

vmin, vmax, ci, cmap = -1.2,1.2,0.1, cmocean.cm.curl
levels = np.arange(vmin,vmax+ci,ci)
for ii in range(len(path)-1):
    S = xr.open_dataset(path[ii+1]+varin+'.nc')['salt']
    anom41 =  ((S.sel(time=slice(date41,date50))   * basin_mask).cf.mean(dim=['time','longitude']) - (ctlSb41  * basin_mask).cf.mean(dim='longitude')).squeeze()
    anom91 =  ((S.sel(time=slice(date91,date100))  * basin_mask).cf.mean(dim=['time','longitude']) - (ctlSb91  * basin_mask).cf.mean(dim='longitude')).squeeze()
    anom141 = ((S.sel(time=slice(date141,date150)) * basin_mask).cf.mean(dim=['time','longitude']) - (ctlSb141 * basin_mask).cf.mean(dim='longitude')).squeeze()

    ax1 = fig.add_subplot(3,2,ii+1)
    anom41.plot.contourf(ax=ax1, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom41.plot.contour(ax=ax1, colors='k', levels=levels,linewidths=0.5)
    ax2 = fig.add_subplot(3,2,ii+3)
    anom91.plot.contourf(ax=ax2, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom91.plot.contour(ax=ax2, colors='k', levels=levels,linewidths=0.5)
    ax3 = fig.add_subplot(3,2,ii+5)
    anom141.plot.contourf(ax=ax3, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom141.plot.contour(ax=ax3, colors='k', levels=levels,linewidths=0.5)

    ax1.set_title(period[0], loc='right', size=labelsize)
    ax1.set_title(label[ii], loc='left', size=axissize)
    ax2.set_title(period[1], loc='right', size=labelsize)
    ax2.set_title(label2[ii], loc='left', size=axissize)
    ax3.set_title(period[2], loc='right', size=labelsize)
    ax3.set_title(label3[ii], loc='left', size=axissize)

    for ax in [ax1, ax2, ax3]:
        ax.set_facecolor([0.5, 0.5,0.5])
        ax.set_ylim([5500,0])
        ax.set_xlim([-80,80])
        if ii==0:
           ax.set_ylabel('Depth (m)', fontsize=axissize)
        else:
           ax.set_ylabel('', fontsize=axissize)
           ax.set_yticklabels([])
        ax3.set_xlabel('Latitude', fontsize=axissize)
        ax1.set_xlabel('')
        ax1.set_xticklabels([])
        ax2.set_xlabel('')
        ax2.set_xticklabels([])

norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
ax = plt.gcf().add_axes((.3,.02,.4,.01))
cb1 = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels,
                                      orientation='horizontal', extend='both',)

for cb in [cb1]:
    cb.ax.set_title(fieldin)
    cb.ax.tick_params(labelsize=8)

if printing:
   plt.savefig('{dirout}{nameplot}_{var}.png'.format(dirout=dirout, nameplot=nameplot[0], var=var[1]),
            format='png', transparent = True, bbox_inches='tight',dpi=300)



#Density
fig = plt.figure(2, figsize=(8,10))
varin = var[2]
fieldin = field[2]

for ii in range(len(path)-1):
    D = xr.open_dataset(path[ii+1]+varin+'.nc')['pot_rho_0']
    anom41 =  ((D.sel(time=slice(date41,date50))   * basin_mask).cf.mean(dim=['time','longitude']) - (ctlDb41  * basin_mask).cf.mean(dim='longitude')).squeeze()
    anom91 =  ((D.sel(time=slice(date91,date100))  * basin_mask).cf.mean(dim=['time','longitude']) - (ctlDb91  * basin_mask).cf.mean(dim='longitude')).squeeze()
    anom141 = ((D.sel(time=slice(date141,date150)) * basin_mask).cf.mean(dim=['time','longitude']) - (ctlDb141 * basin_mask).cf.mean(dim='longitude')).squeeze()

    ax1 = fig.add_subplot(3,2,ii+1)
    anom41.plot.contourf(ax=ax1, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom41.plot.contour(ax=ax1, colors='k', levels=levels,linewidths=0.5)
    ax2 = fig.add_subplot(3,2,ii+3)
    anom91.plot.contourf(ax=ax2, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom91.plot.contour(ax=ax2, colors='k', levels=levels,linewidths=0.5)
    ax3 = fig.add_subplot(3,2,ii+5)
    anom141.plot.contourf(ax=ax3, levels=levels,extend='both', cmap=cmap, add_colorbar=False)
    anom141.plot.contour(ax=ax3, colors='k', levels=levels,linewidths=0.5)

    ax1.set_title(period[0], loc='right', size=labelsize)
    ax1.set_title(label[ii], loc='left', size=axissize)
    ax2.set_title(period[1], loc='right', size=labelsize)
    ax2.set_title(label2[ii], loc='left', size=axissize)
    ax3.set_title(period[2], loc='right', size=labelsize)
    ax3.set_title(label3[ii], loc='left', size=axissize)

    for ax in [ax1, ax2, ax3]:
        ax.set_facecolor([0.5, 0.5,0.5])
        ax.set_ylim([5500,0])
        ax.set_xlim([-80,80])
        if ii==0:
           ax.set_ylabel('Depth (m)', fontsize=axissize)
        else:
           ax.set_ylabel('', fontsize=axissize)
           ax.set_yticklabels([])
        ax3.set_xlabel('Latitude', fontsize=axissize)
        ax1.set_xlabel('')
        ax1.set_xticklabels([])
        ax2.set_xlabel('')
        ax2.set_xticklabels([])

norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
ax = plt.gcf().add_axes((.3,.02,.4,.01))
cb1 = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels,
                                      orientation='horizontal', extend='both',)

for cb in [cb1]:
    cb.ax.set_title(fieldin)
    cb.ax.tick_params(labelsize=8)

if printing:
   plt.savefig('{dirout}{nameplot}_{var}.png'.format(dirout=dirout, nameplot=nameplot[0], var=var[2]),
            format='png', transparent = True, bbox_inches='tight',dpi=300)

plt.show()







