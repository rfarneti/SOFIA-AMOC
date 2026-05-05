'''
   Calculate the global Meridional Overturning Circulation 
   - in density space - using output from MOM5

   One can apply a mask for the Atlantic or IndoPacific basin

'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
from scipy.interpolate import interp1d
#import cftime
#import datetime
#import nc_time_axis
#import netCDF4 as nc
import matplotlib.pyplot as plt
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature


expt = ['0.0Sv']

pathin   = '/Volumes/Antwater/DATA/'

printing = True
dirout = './Figures_Boeira/'
nameplot = ['Baro_Psi']


# create dictionaries for all experiments
psi_args = {
    "0.0Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv",
        "file1": "tx_trans_int_z.nc",
        "start_time": "1651-07-02",
        "end_time"  : "1690-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 200,"grid_xt_ocean": -1,}
    },
}


def compute_baro_psi(expt):
    """Compute the Barotropic Stream Function from the
       vertically integrated transport. 
       Adds the Drake Passage transport"""
    dirin      = psi_args[expt]["dirin"]
    file1      = psi_args[expt]["file1"]

    psi   = xr.open_dataset(pathin+dirin+'/pp/'+file1)['tx_trans_int_z']
    drake = psi.sel(xu_ocean=-70, yt_ocean=slice(-80,-50)).cf.sum("latitude") 
    
    mask = xr.where(np.isnan(psi), np.nan, 1)
    mask2 = xr.where(np.isnan(psi), 1, np.nan)
    #fig = plt.figure(1, figsize=(12, 8))
    #ax1 = plt.subplot(2,1,1) 
    #ax2 = plt.subplot(2,1,2) 
    #mask.isel(time=0).plot.contourf(ax=ax1)
    #mask2.isel(time=0).plot.contourf(ax=ax2)
    
    psiu = -1 * psi.cf.cumsum("latitude") + drake.cf.mean("time") 
    return psiu*mask


psiu = compute_baro_psi(expt[0])



labelsize = 10
titlesize = 12
axissize = 14

start_time = psi_args[expt[0]]["start_time"]
end_time   = psi_args[expt[0]]["end_time"]

levels = np.arange(-80, 180, 20) 
fig = plt.figure(1, figsize=(10, 5))
ax1 = fig.add_subplot(1,1,1, projection=ccrs.Robinson(central_longitude=-140))
psiu.sel(time=slice(start_time, end_time)).cf.mean("time").plot.contourf(ax=ax1, levels=levels,
         cmap=cm.cm.speed, extend='both', cbar_kwargs={'shrink': 0.7, 'label': '[Sv]'},
         transform=ccrs.PlateCarree())
psiu.sel(time=slice(start_time, end_time)).cf.mean("time").plot.contour(ax=ax1, levels=levels, 
         #linestyles=np.where(levels >= 0, "-", "--"), colors='k', linewidths=1.)
         colors='k', linewidths=1.,
         transform=ccrs.PlateCarree())

ax1.set_global()
ax1.coastlines()
#plt.ylim((-90, 90))
#plt.xlim([-280, 80])
#plt.ylabel('Latitude', fontsize=axissize)
#plt.xlabel('Longitude [$^\circ$]', fontsize=axissize)
#plt.title(r'{number}'.format(number=number[ii]), loc='left', size=axissize);
#plt.title(r'{basin} Residual Overturning Circulation'.format(basin=basin), loc='center', size=titlesize);
#plt.title('({expt})'.format(expt=expt[ii]), loc='right', size=axissize);

#plt.tight_layout()
if printing:
   plt.savefig('{dirout}{nameplot}_{expt}.png'.format(dirout=dirout,
               nameplot=nameplot[0], expt=expt[0]),
               format='png', transparent = True, bbox_inches='tight',dpi=300)



plt.show()

# Saves to a NetCDF file
dirin = psi_args[expt[0]]["dirin"]
path = pathin+dirin+'/pp/'+'quasi_bsf.nc'

ds = psiu.rename("quasi_bsf")
ds.sel(time=slice(start_time, end_time)).cf.mean("time").to_netcdf(path=path, 
       mode="w", format="NetCDF4", compute=True)
psiu.close()


