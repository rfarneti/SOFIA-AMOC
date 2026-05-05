'''
   Calculate the global Meridional Overturning Circulation 
   - in density space - using output from MOM5
   One can apply a mask for the Atlantic or IndoPacific basin
   Then 

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
import cmocean
import cartopy.feature as cfeature
from ocean_basins import basin_mask_ht

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif"
})


# Choose the basin: Global, Atlantic, or IndoPacific
basin = 'Atlantic'

expt = ['0.1Sv','1.0Sv','1.0Sv_60S']
#expt = ['0.1Sv']
forcing = ['(0.1Sv)','(1.0Sv)','(1.0Sv_60S)']
path   = '/Volumes/Antwater/DATA/'
dirout = '/Volumes/Antwater/DATA/MOC'

printing = True
dirout = './Figures_transports/'
nameplot = ['AMOC_hovmoller_rho']


# create dictionaries for all experiments
psi_args = {
#    "0.0Sv": {
#        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv",
#        "file1": "ty_trans_rho.nc",
#        "file2": "ty_trans_rho_gm.nc",
#        "start_time": "1601-07-02",
#        "end_time"  : "1650-07-02",
#        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 200,"grid_xt_ocean": -1,}
#    },
    "0.1Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.1Sv",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1601-07-02",
        "end_time"  : "1750-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 200,"grid_xt_ocean": -1,}
    },
    "1.0Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1601-07-02",
        "end_time"  : "1750-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 1080,"grid_xt_ocean": -1,}
    },
    "1.0Sv_60S": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv_60S",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1601-07-02",
        "end_time"  : "1750-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 1080,"grid_xt_ocean": -1,}
    },
}


def compute_psi_rho(expt, mask=None):
    """Load ty_trans_rho - and sum zonally. Also, if there is a ty_trans_rho_gm variable saved, 
    assume that GM is switched on and load that as well. 
    Most MOM simulations save transport with units of kg/s - convert to Sv."""
    dirin      = psi_args[expt]["dirin"]
    file1      = psi_args[expt]["file1"]
    file2      = psi_args[expt]["file2"]

    psi   = xr.open_dataset(path+dirin+'/pp/'+file1, decode_times=False)['ty_trans_rho']
    psiGM = 0 * psi.copy(deep=True)
    psiGM = xr.open_dataset(path+dirin+'/pp/'+file2, decode_times=False)['ty_trans_rho_gm']
    #convert to volume transport
    #rho0 = 1025
    #psi = psi / rho0
    #psiGM = psiGM / rho0
    if mask is not None:
       psi = psi * mask
       psiGM = psiGM * mask
    psi = psi.cf.sum("longitude")
    psiGM = psiGM.cf.sum("longitude")
    """This is equivalent to the Ferret computation
       psi[i=@sum,k=@rsum]-psi[i=@sum,k=@sum]
       but I'm not sure why here I don't need
       to remove the vertical sum""" 
    psi = psi.cf.cumsum("vertical")  ##- psi.cf.sum("vertical")
    psi = psi + psiGM
    psi.load()
    return psi


# Apply a mask if needed
if basin != 'Global':
   print('Basin is', basin)
   mask = basin_mask_ht(basin, SO=True)

   dirin = psi_args[expt[0]]["dirin"]
   file1 = psi_args[expt[0]]["file1"]
   tyrho = xr.open_dataset(path+dirin+'/pp/'+file1)['ty_trans_rho']
   mask.coords['xt_ocean'] = tyrho.grid_xt_ocean.values
   mask.coords['yt_ocean'] = tyrho.grid_yu_ocean.values
   mask = mask.rename({'xt_ocean': 'grid_xt_ocean', 'yt_ocean': 'grid_yu_ocean'})
   del tyrho
   print('I am appying a mask for', basin)
   print('mask shape', mask.shape)
else:
   print('Basin is', basin)
   print('I am not applying a mask')



number = ['a)','b)','c)']
titlesize = 12
axissize = 16

vmin, vmax, ci, cmap = -25,25,2.5, cmocean.cm.delta
levels = np.arange(vmin,vmax+ci,ci)
yticks = np.array([1030, 1033, 1034, 1035, 1035.5,1036,1036.5, 1036.8,1037])
scfac = 5.5  ## A power to set the stretching

#fig,ax = plt.subplots(nrows = 1, ncols = 1, figsize = (8,5))
fig = plt.figure(1, figsize=(8,10))
for ii in range(len(expt)):
    ax = fig.add_subplot(3,1,ii+1)
    if basin != 'Global':
       psi   = compute_psi_rho(expt[ii], mask=mask)
    else:
       psi   = compute_psi_rho(expt[ii])
    psi['time'] = psi.time/365-1600
    psi_max = psi.sel(grid_yu_ocean=slice(40,60)).max('grid_yu_ocean').squeeze()
    #print('psi_max = ',psi_max)
    #psi_max = psi_max.where(psi_max >=1.0)
    ax.contourf(psi_max.time, (psi_max.potrho-1028)**scfac, psi_max.transpose(), 
              cmap=cmap, levels=levels)
    ax.contour(psi_max.time, (psi_max.potrho-1028)**scfac, psi_max.transpose(), 
              levels=[0.0,], colors='k', linewidths=0.5)
    #print('dens_max = ',dens_max)
    #ax.plot(psi_max.time, dens_max, color='w', linestyle='-')
    #ax.axvline(x=100, color = 'k', linestyle='dashed', linewidth = 0.5)
    ax.set_title(forcing[ii], loc='right', size=titlesize)
    ax.set_title('{number}'.format(number=number[ii]), loc='left', size=axissize)
    ax.set_yticks((yticks-1028)**scfac)
    ax.set_yticklabels(yticks, fontsize=8)
    ax.set_ylim([0.5**scfac, 9.2**scfac])
    ax.invert_yaxis()
    ax.set_xlim([0,150])
    ax.set_ylabel(r'$\sigma_2$ [kg m$^{-3}$]', fontsize=axissize)
    if ii==2:
       ax.set_xlabel('Time [years]', fontsize=axissize)
    else:
       ax.set_xlabel('')
       ax.set_xticklabels([])

norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
ax1 = plt.gcf().add_axes((.25,.04,.5,.01))
cb1 = matplotlib.colorbar.ColorbarBase(ax=ax1, cmap=cmap, norm=norm, boundaries=levels,
                                    orientation='horizontal', extend='both',)
for cb in [cb1]:
    cb.ax.set_xlabel('[Sv]')
    cb.ax.tick_params(labelsize=8)

    if printing:
       plt.savefig('{dirout}{nameplot}.png'.format(dirout=dirout,
                   nameplot=nameplot[0]),
                   format='png', transparent = True, bbox_inches='tight',dpi=300)

plt.show()



