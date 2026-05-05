'''
   Calculate in situ density 
   and compares with model diagnostic
   
   - calculate rho

'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
#import glob
import momsofia as ms
#from scipy.interpolate import interp1d
#import cftime
#import datetime
#import nc_time_axis
#import netCDF4 as nc
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colorbar
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

expt = ['0.0Sv']

pathin   = '/Volumes/Antwater/DATA/'

printing = False
dirout = './Figures_Boeira/'

# Create dictionaries for all experiments
psi_args = {
    "0.0Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv",
        "u": "u.nc",
        "v": "v.nc",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
        "start_time": "1691-07-02",
        "end_time"  : "1700-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": -1,"grid_xt_ocean": -1,}
    },
}

start_time = psi_args[expt[0]]["start_time"]
end_time   = psi_args[expt[0]]["end_time"]

dirin   = psi_args[expt[0]]["dirin"]
u       = psi_args[expt[0]]["u"]
v       = psi_args[expt[0]]["v"]
t       = psi_args[expt[0]]["t"]
s       = psi_args[expt[0]]["s"]
p       = psi_args[expt[0]]["p"]
rho     = psi_args[expt[0]]["rho"]
potrho0 = psi_args[expt[0]]["potrho0"]
potrho2 = psi_args[expt[0]]["potrho2"]

# Get the variables
u      = xr.open_dataset(pathin+dirin+'/pp/'+u)['u']
v      = xr.open_dataset(pathin+dirin+'/pp/'+v)['v']
thetao = xr.open_dataset(pathin+dirin+'/pp/'+t)['temp']
so     = xr.open_dataset(pathin+dirin+'/pp/'+s)['salt']
press  = xr.open_dataset(pathin+dirin+'/pp/'+p)['press']
rho    = xr.open_dataset(pathin+dirin+'/pp/'+rho)['rho']
pot_rho_0 = xr.open_dataset(pathin+dirin+'/pp/'+potrho0)['pot_rho_0']
pot_rho_2 = xr.open_dataset(pathin+dirin+'/pp/'+potrho2)['pot_rho_2']

u      = u.sel(time=slice(start_time, end_time))
v      = v.sel(time=slice(start_time, end_time))
thetao = thetao.sel(time=slice(start_time, end_time))
so     = so.sel(time=slice(start_time, end_time))
press  = press.sel(time=slice(start_time, end_time))
ho    = rho.sel(time=slice(start_time, end_time))
pot_rho_0 = pot_rho_0.sel(time=slice(start_time, end_time))
pot_rho_2 = pot_rho_2.sel(time=slice(start_time, end_time))

# Get some static fields
dxu = xr.open_dataset(pathin+dirin+'/pp/'+'dxu.nc')['dxu']
dyt = xr.open_dataset(pathin+dirin+'/pp/'+'dyt.nc')['dyt']
dzt = xr.open_dataset(pathin+dirin+'/pp/'+'dht.nc')['dht']
area_t = xr.open_dataset(pathin+dirin+'/pp/'+'areacello_t.nc')['area_t']
area_u = xr.open_dataset(pathin+dirin+'/pp/'+'areacello_u.nc')['area_u']
lon = xr.open_dataset(pathin+dirin+'/history/'+'17010101.ocean_grid.nc')["geolon_c"]
lat = xr.open_dataset(pathin+dirin+'/history/'+'17010101.ocean_grid.nc')["geolat_c"]

# 1/ Compute in-situ density and compare with online diagnostic 
print('--- Computing in situ density')
rho_d = ms.derived.calc_rho(thetao, so, press, eos="Wright")
print('rho_d = ', rho_d.dims)

fig = plt.figure(6, figsize=(10,8))
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
levels = np.arange(1006, 1029, 1)
rho_d.isel(st_ocean=1).isel(time=1).plot.contourf(ax=ax1, levels=levels,
           cmap=cm.cm.thermal)
ax1.set_title('in situ density derived')

rho.isel(st_ocean=1).isel(time=1).plot.contourf(ax=ax2, levels=levels,
           cmap=cm.cm.thermal)
ax2.set_title('in situ density from model')
plt.tight_layout()

plt.show()


