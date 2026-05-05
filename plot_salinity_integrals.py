import xarray as xr
import os
import cftime
#import xesmf as xe
import numpy as np
import matplotlib.pyplot as plt
import my_style
from pathlib import Path

BASE_PATH = Path('/Volumes/Antwater2/DATA/')
DIR_OUT = Path('/Volumes/Antwater2/Figures_tracers/')
PREFIX = Path('FAFANTWATER_ideal_runoff_Boeira_rest_ice_')

nameplot   = 'salinity_integral'


area = xr.load_dataset(BASE_PATH / f'{PREFIX}0.0Sv' / 'pp' / 'areacello_t.nc')['area_t']
dsc  = xr.load_dataset(BASE_PATH / f'{PREFIX}0.0Sv' / 'pp' / 'so.nc')
ds01 = xr.load_dataset(BASE_PATH / f'{PREFIX}0.1Sv' / 'pp' / 'so.nc')
ds1  = xr.load_dataset(BASE_PATH / f'{PREFIX}1.0Sv' / 'pp' / 'so.nc')
#dhtc = xr.load_dataset('/Volumes/Antwater/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv/pp/dht.nc')
#dht1 = xr.load_dataset('/Volumes/Antwater/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.1Sv/pp/dht.nc')

NZ = 49
NT = 149
print('area=', area.shape)
print('st_ocean.diff=', dsc.st_ocean.diff(dim='st_ocean').shape)
print('salt=', dsc.salt.isel(st_ocean=slice(0,NZ)).shape)

# Test the vertical levels
fig = plt.figure(1)
ax1 = fig.add_subplot(1,1,1)
ax1.plot(dsc.st_ocean.diff(dim='st_ocean'), color='b', marker='s', label='st_ocean.diff')
#ax1.plot(dhtc.dht.isel(xt_ocean=190, yt_ocean=80).isel(time=0), color='r', marker='o', label='dht')


#--- Get the volume-integrated Salinities
fig = plt.figure(2)
ax1 = fig.add_subplot(1,1,1)

diff01 = []
diff1  = []
for i in np.arange(0,NT,1):    
    print(f'Computing Salinity integral for time {i}')
    salc  =  dsc.salt.isel(st_ocean=slice(0,NZ)).isel(time=i)
    sal01 = ds01.salt.isel(st_ocean=slice(0,NZ)).isel(time=i)
    sal1  =  ds1.salt.isel(st_ocean=slice(0,NZ)).isel(time=i)
    #dzc   = dhtc.dht.isel(st_ocean=slice(0,NZ)).isel(time=i)
    #dza   = dht1.dht.isel(st_ocean=slice(0,NZ)).isel(time=i)
    
    soc  = ( salc  * dsc.st_ocean.diff(dim='st_ocean')  * area ).sum(dim=['yt_ocean', 'xt_ocean', 'st_ocean'])
    so01 = ( sal01 * ds01.st_ocean.diff(dim='st_ocean') * area ).sum(dim=['yt_ocean', 'xt_ocean', 'st_ocean'])
    so1  = ( sal1  *  ds1.st_ocean.diff(dim='st_ocean') * area ).sum(dim=['yt_ocean', 'xt_ocean', 'st_ocean'])
    #soc2 = ( salc * dzc                               * area ).sum(dim=['yt_ocean', 'xt_ocean', 'st_ocean'])
    #soa2 = ( sala * dza                               * area ).sum(dim=['yt_ocean', 'xt_ocean', 'st_ocean'])
    
    diff01.append((so01 - soc).values)
    diff1.append((so1 - soc).values)
  
    #print('soc=', soc.values) 
    ax1.plot(i, soc,  color='k',    marker='o', label='control')
    ax1.plot(i, so01, color='blue', marker='o', label='0.1Sv')
    ax1.plot(i, so1,  color='red',  marker='o', label='1.0Sv')


#--- Calculate the expected freshening
# Constants
flow_rate_sverdrup_01 = 0.1    # Flow rate in Sv [m^3/s]
flow_rate_sverdrup_1  = 1.0    # Flow rate in Sv [m^3/s]
seconds_per_month  = 2592000   # Number of seconds in a month
seconds_per_year   = 31104000  # Number of seconds in a year
density_freshwater = 1000      # Density of freshwater in kg/m³
total_years        = 150       # Total number of years

# Compute the volume of freshwater added per year in cubic meters
volume_per_year_01 = flow_rate_sverdrup_01 * seconds_per_year  # in m³
volume_per_year_1  = flow_rate_sverdrup_1  * seconds_per_year  # in m³
# Compute the mass   of freshwater added per year in kilograms
mass_per_year_01 = volume_per_year_01 * density_freshwater  # in kg
mass_per_year_1  = volume_per_year_1  * density_freshwater  # in kg
# Compute the cumulative amount of freshwater added each year
cumulative_fw_added_01 = [mass_per_year_01 * (year + 1) for year in range(total_years)]
cumulative_fw_added_1  = [mass_per_year_1  * (year + 1) for year in range(total_years)]

expected_01 = np.array(cumulative_fw_added_01)* -1 * 30000
expected_1  = np.array(cumulative_fw_added_1) * -1 * 30000


fig = plt.figure(3)
ax = fig.add_subplot(1,1,1)
ax.plot(expected_01[0:NT:1], color='black', linestyle='solid', label='0.1Sv')
ax.plot(diff01,              color='black', linestyle='dashed',label='MOM5 - 0.1Sv')
ax.plot(expected_1[0:NT:1],  color='gray',  linestyle='solid', label='1.0Sv')
ax.plot(diff1,               color='gray',  linestyle='dashed',label='MOM5 - 1.0Sv')
ax.set_xlabel("Time [Years]",fontsize=14)
ax.legend(loc='lower left', ncol=1, fontsize=12)
ax.set_title('volume integrated Salinity anomaly', loc='center', size=14)

plt.savefig(DIR_OUT / f'salinity_integral.png', bbox_inches='tight',dpi=300)

plt.show()
