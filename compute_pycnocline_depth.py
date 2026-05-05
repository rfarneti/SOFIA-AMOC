'''
 Compute and Plot the Pycnocline Depth in the Atlantic 
 from pot_rho_0 or pot_rho_2
 
'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
import momsofia as ms
#import pint_xarray
import cftime
import datetime
import nc_time_axis 
#import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colorbar
import cmocean
from ocean_basins import basin_mask_ht_full
from pathlib import Path
from scipy.ndimage import distance_transform_edt
import my_style

#SIGMA = 0.0
SIGMA = 2000.0

BASE_PATH  = Path('./DATA/')
DIR_OUT = Path('./Figures_tracers/')
PRINTING = True

OCEAN    = 'Atlantic'
FORCING  = [
    ("0.0Sv",     "CTL",        "black",      "solid", 1.0),  # Control
    #("0.01Sv",    "0.01Sv",     "salmon",     "solid", 0.5),
    #("0.02Sv",    "0.02Sv",     "tomato",     "solid", 0.5),
    #("0.05Sv",    "0.05Sv",     "orangered",  "solid", 0.5),
    #("0.07Sv",    "0.07Sv",     "brown",      "solid", 0.5),
    ("0.1Sv",     "0.1Sv",      "black",      "solid", 2.0),
    #("0.1Sv_60S", "0.1Sv(60S)", "black",      "dashed",2.0),
    #("0.1Sv_tas_corrected", "0.1Sv$^{*}$)", "black",      "dashdot",2.0),
    #("0.15Sv",    "0.15Sv",     "aquamarine", "solid", 0.5),
    ("0.2Sv",     "0.2Sv",      "cyan",       "solid", 0.5),
    ("0.3Sv",     "0.3Sv",      "lime",       "solid", 0.5),
    ("0.4Sv",     "0.4Sv",      "green",      "solid", 0.5),
    ("0.5Sv",     "0.5Sv",      "royalblue",  "solid", 0.5),
    ("0.6Sv",     "0.6Sv",      "blue",       "solid", 0.5),
    ("0.7Sv",     "0.7Sv",      "blueviolet", "solid", 0.5),
    ("0.8Sv",     "0.8Sv",      "violet",     "solid", 0.5),
    ("0.9Sv",     "0.9Sv",      "magenta",    "solid", 0.5),
    ("1.0Sv",     "1.0Sv",      "purple",     "solid", 2.0),
    ("1.0Sv_60S", "1.0Sv(60S)", "purple",     "dashed",2.0)
]


LAT_LIMS=(-60, 60)
D_LIMS=(0, 2500)

OFFSET_YEAR = 1600;
DAYS_IN_YEAR = 365


def mask_coastal_buffer(ds, degree_buffer=3.0):
    """ Masks out data within a specified degree distance from the coastline.
    Assumes a relatively regular lat/lon grid. """
    # Identify Land (where the data is NaN)
    land_mask = ds.isel(time=0).isnull() 
    
    # Compute Euclidean Distance to the nearest 'False' (Ocean) point.
    # We invert land_mask so land is True/1 and ocean is False/0.
    # distance_transform_edt calculates distance to the closest zero.
    pixel_dist = distance_transform_edt(land_mask == False)
    
    # Convert Degree threshold to Pixel threshold
    # Estimate resolution (degrees per pixel)
    dlat = np.abs(ds.yt_ocean.diff('yt_ocean').mean().values)
    pixel_threshold = degree_buffer / dlat
    
    # Create the mask: Keep only points where distance > threshold
    ocean_mask = xr.DataArray(pixel_dist > pixel_threshold, 
                             coords=land_mask.coords, 
                             dims=land_mask.dims)
    return ds.where(ocean_mask)

def calculate_pycnocline_depth(sigma):
    """Computes pycnocline depth (D) using Xarray integration.
    Handles 3D (lat, lon, depth) 
    or 4D (time, lat, lon, depth)"""
    sigma_max = sigma.max(dim='st_ocean') #local max
    delta_sigma = sigma - sigma_max
    #z = sigma['st_ocean']
    z = sigma['st_ocean'].where(sigma.notnull())    
    
    num = (z * delta_sigma).fillna(0).integrate(coord='st_ocean')
    den = delta_sigma.fillna(0).integrate(coord='st_ocean')
    D = num / den.where(den != 0)

    return D.where(np.isfinite(D))


# Generate Basin Mask for the chosen sector
basin_mask = ms.ocean_basins.basin_mask_ht_full(OCEAN, check=True)
local_mask = basin_mask.sel(yt_ocean=slice(*LAT_LIMS))

area_path = BASE_PATH / f'FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv' / 'pp'
area = xr.open_dataset(area_path / 'areacello_t.nc')['area_t']
local_area = area.sel(yt_ocean=slice(*LAT_LIMS)) * local_mask


fig, ax1 = plt.subplots(1, 1, figsize=(8,8))
for i, (forcing, label, color, style, width) in enumerate(FORCING):
    print(f"Processing {label}...")
    data_path = BASE_PATH / f'FAFANTWATER_ideal_runoff_Boeira_rest_ice_{forcing}' / 'pp'
    var_name = 'pot_rho_0' if SIGMA == 0.0 else 'pot_rho_2'
    ds = xr.open_dataset(data_path /f'{var_name}.nc', decode_times=False)
    sigma = (ds[var_name] - 1000.0).sel(yt_ocean=slice(*LAT_LIMS))
    local_sigma = sigma * local_mask

    if SIGMA == 0.0:
       D = calculate_pycnocline_depth(local_sigma)
    else:
       D = calculate_pycnocline_depth(local_sigma.sel(st_ocean=slice(*D_LIMS)))
    D = D * 2.0

    # Mask out coastal points
    D_ocean = mask_coastal_buffer(D, degree_buffer=3.0)
    ###D_ocean.isel(time=0).plot.contourf(ax=ax1)
    
    Depth = (D_ocean * local_area).cf.mean(dim=['longitude', 'latitude']) / local_area.cf.mean(['longitude','latitude'])
    Depth['time'] = Depth['time'] / DAYS_IN_YEAR - OFFSET_YEAR 

    Depth.plot(ax=ax1, color=color, linestyle=style, linewidth=width, label=label)
    ax1.set_ylim([800,1000]) if SIGMA == 0.0 else ax1.set_ylim([700,900]) 
    #ax1.set_yticks(np.arange(725,925,25)) 
    ax1.set_xlim([0,160])
    ax1.set_xticks(np.arange(0,180,20))
    ax1.tick_params(labelsize=12) 
    ax1.set_ylabel(rf'Depth [m]', fontsize=16)
    ax1.set_xlabel("Time [Years]",fontsize=16)
    ax1.legend(loc='upper left', ncol=4, fontsize=10, frameon=True)
    ax1.grid(True, alpha=0.3) 
    ax1.set_title(f'Pycnocline depth', loc='center', size=16);

    #--- save pycnocline depth
    Depth = Depth.rename("Depth")
    Depth.to_netcdf(BASE_PATH / f'pycnocline_depth_sigma{SIGMA}_{forcing}.nc')

if PRINTING:
   plt.savefig(DIR_OUT / f'pycnocline_depth_sigma{SIGMA}.png', bbox_inches='tight',dpi=300)

plt.show()

