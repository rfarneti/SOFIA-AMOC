'''

 Plot time series of volume-weighted zonal and vertical mean of pot_rho_0 in the Atlantic,
 for different freshwater forcings

 \Delta \sigma = \sigma_n - \sigma_b

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
from ocean_basins import basin_mask_ht_full
from pathlib import Path
import my_style_latex

# Choose basin, density (in-situ, sigma0, sigma2), strength of forcing
OCEAN   = 'Atlantic'
SIGMA   = 0.0
FORCING = ['0.1Sv','1.0Sv']#,'1.0Sv_60S']
LSTYLE  = ['-','--',':']
COLORS = {'north': 'blue', 'basin': 'red'}
SB_LIMS = (-65, 40)
SN_LIMS = (40, 65)

# Set path
BASE_PATH = Path('/Volumes/Antwater2/DATA/')
DIR_OUT = Path('./Figures_transports/')
PRINTING = True

# Constants for time conversion
DAYS_IN_YEAR = 365
OFFSET_YEAR = 1600

# Generate Basin Mask for the chosen sector
basin_mask = ms.ocean_basins.basin_mask_ht_full(OCEAN, check=False)


def sigma_zonal_mean(forcing, mask):
    """Compute volume-weighted zonal mean density anomaly over time.."""
    DATA_PATH = BASE_PATH / f'FAFANTWATER_ideal_runoff_Boeira_rest_ice_{forcing}'
    dens = xr.open_dataset(DATA_PATH / 'pp' / 'pot_rho_0.nc', decode_times = False)['pot_rho_0']
    area = xr.open_dataset(DATA_PATH / 'pp'/ 'areacello_t.nc')['area_t']
    dzt = xr.open_dataset(DATA_PATH / 'pp' / 'dht.nc', decode_times = False)['dht']
    volume = ms.derived.calc_volume(area, dzt) * mask
    dens_z = (dens * mask * volume).cf.mean(dim='longitude') / volume.cf.mean('longitude')
    time_years = dens_z.time / DAYS_IN_YEAR - OFFSET_YEAR
    dens_z = dens_z.assign_coords(time=time_years)
    return dens_z - 1000.0


fig, ax1 = plt.subplots(figsize=(10,6))
ax2 = ax1.twinx()
axissize  = 16
labelsize = 16

for forcing, ls in zip(FORCING, LSTYLE):
    densz = sigma_zonal_mean(forcing, basin_mask)
    rho_n = densz.sel(yt_ocean=(slice(*SN_LIMS))).cf.mean(['latitude', 'vertical'])
    rho_b = densz.sel(yt_ocean=(slice(*SB_LIMS))).cf.mean(['latitude', 'vertical'])
    delta_rho = rho_n - rho_b
    
    # save Delta rho
    delta_rho = delta_rho.rename("delta_rho")
    delta_rho.to_netcdf(BASE_PATH / f'delta_rho_sigma{SIGMA}_{forcing}.nc')

    anom_n = rho_n - rho_n.isel(time=0)
    anom_b = rho_b - rho_b.isel(time=0)
    rn, = ax1.plot(anom_b.time, anom_n,
             color=COLORS['north'], linestyle=ls, linewidth=2.0, label=f'{forcing}')
    rb, = ax2.plot(anom_b.time, anom_b,
             color=COLORS['basin'], linestyle=ls, linewidth=2.0, label=f'{forcing}')

for ax in [ax1, ax2]:
    ax.set_xlim([0,150])
    ax.set_ylim([-0.7,0.1])
    ax.set_xticks(np.arange(0,160,10), fontsize=axissize)
    ax.set_yticks(np.arange(-0.7,0.1,0.10), fontsize=axissize)
    ax.tick_params(labelsize=axissize)
    ax.grid(color='grey', linestyle='solid', linewidth=0.2)
ax1.set_ylabel(r'$\delta \langle \sigma_n \rangle$ [kg m$^{-3}]$', fontsize=labelsize, color=COLORS['north'])
ax2.set_ylabel(r'$\delta \langle \sigma_b \rangle$ [kg m$^{-3}]$', fontsize=labelsize, color=COLORS['basin'])
ax1.set_xlabel('Time [Years]', fontsize=labelsize)
ax1.spines['left'].set_color(COLORS['north'])
ax2.spines['right'].set_color(COLORS['basin'])
ax1.tick_params(axis='y', colors=COLORS['north'])
ax2.tick_params(axis='y', colors=COLORS['basin'])
ax1.set_title(r"a)", loc='left', fontsize=20)
ax1.axhline(y=-0.51,  xmin=0.07, xmax=0.12, color='blue', linestyle='-', linewidth=2)
ax1.axhline(y=-0.52, xmin=0.07, xmax=0.12, color='red',  linestyle='-', linewidth=2)
ax1.axhline(y=-0.57,  xmin=0.07, xmax=0.12, color='blue', linestyle='--', linewidth=2)
ax1.axhline(y=-0.58, xmin=0.07, xmax=0.12, color='red',  linestyle='--', linewidth=2)
ax1.text(20,-0.525, '(0.1 Sv)', fontsize=12)
ax1.text(20,-0.585, '(1.0 Sv)', fontsize=12)
if PRINTING:
   plt.savefig(DIR_OUT / f'sigman_sigmab_boxes_time_series_sigma{SIGMA}.png',
            bbox_inches='tight',dpi=300)


plt.show()


