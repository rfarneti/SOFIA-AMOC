'''
   Calculate indeces for the Meridional Overturning Circulation 
   - in density space - using output from MOM5

   One can apply a mask for the Atlantic or IndoPacific basin

'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
#import cftime
#import datetime
#import nc_time_axis
from pylab import *
import matplotlib.pyplot as plt
import momsofia as ms
from pathlib import Path
import my_style


BASIN = 'Atlantic'

BASE_PATH = Path('/Volumes/Antwater2/DATA/')
PREFIX    = "FAFANTWATER_ideal_runoff_"

DIR_OUT = Path('./Figures_transports/')

PRINTING = True
          

EXPERIMENTS = [
    ("0.0Sv",     "CTL",        "black",      "solid", 1.0),  # Control
    ("0.01Sv",    "0.01Sv",     "salmon",     "solid", 0.5),
    ("0.02Sv",    "0.02Sv",     "tomato",     "solid", 0.5),
    ("0.05Sv",    "0.05Sv",     "orangered",  "solid", 0.5),
    ("0.07Sv",    "0.07Sv",     "brown",      "solid", 0.5),
    ("0.1Sv",     "0.1Sv",      "black",      "solid", 2.0),
    ("0.1Sv_60S", "0.1Sv(60S)", "black",      "dashed",2.0),
    ("0.15Sv",    "0.15Sv",     "aquamarine", "solid", 0.5),
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


# --- Generate and later apply a mask, if needed
if BASIN != 'Global':
   print(f'---> Setting up mask for {BASIN}')
   mask = ms.ocean_basins.basin_mask_ht_full(BASIN, check=False)
   #mask = ms.ocean_basins.basin_mask_ht(BASIN, SO=True, check=False)
   ds = xr.open_dataset(BASE_PATH / f"{PREFIX}pv60_rest_ice_0.0Sv" / 'pp' /'ty_trans_rho.nc')
   mask = mask.rename({'xt_ocean': 'grid_xt_ocean', 'yt_ocean': 'grid_yu_ocean'})
   mask = mask.assign_coords(grid_xt_ocean=ds.grid_xt_ocean, grid_yu_ocean=ds.grid_yu_ocean)
else:
   print(f'---> Basin is {BASIN} and I am not applying a mask')


LATITUDES = [-32, 30, 40, 50]
D_RANGE = (1036.0, 1037.2) # NADW density range

OFFSET_YEAR = 1600
DAYS_IN_YEAR = 365

start_time,end_time,dt = 0, 160, 10
ymin, ymax, dy = 4, 22, 2.
titlesize, axissize, legendsize, max_psi  = 16, 14, 10, 30

figs = []
axes = []
for i, lat in enumerate(LATITUDES):
    fig, ax = plt.subplots(figsize=(10, 6), num=i+1)
    figs.append(fig)
    axes.append(ax)

# --- Main Processing and Plotting
for exp_id, label, color, style, width in EXPERIMENTS:    
    print(f"Processing Experiment: {exp_id}") 
 
    suffix = "pv60_rest_ice_0.0Sv/pp" if exp_id == "0.0Sv" else f"Boeira_rest_ice_{exp_id}/pp"
    file_path = BASE_PATH / f"{PREFIX}{suffix}"
    
    ds_psi = xr.open_dataset(file_path / 'ty_trans_rho.nc', chunks={'time': 10}, decode_times=False)
    ds_gm = xr.open_dataset(file_path / 'ty_trans_rho_gm.nc', chunks={'time': 10}, decode_times=False)
    psi = ds_psi['ty_trans_rho']
    psi_gm = ds_gm['ty_trans_rho_gm']
    
    amoc = ms.MOC.calc_psi_rho(psi, psi_gm, mask=mask if BASIN != 'Global' else None)
    time_years = amoc.time / DAYS_IN_YEAR - OFFSET_YEAR
    amoc = amoc.assign_coords(time=time_years)

    for i, lat in enumerate(LATITUDES):
        print(f"Plotting Latitude: {lat}") 
        ax = axes[i]
        fig = figs[i]
        plt.figure(fig.number)

        amoc_lat = amoc.sel(grid_yu_ocean=float(lat), method='nearest')
        amoc_index = amoc_lat.sel(potrho=slice(*D_RANGE)).max(dim="potrho")
        amoc_index.plot(ax=ax, color=color, linestyle=style, linewidth=width, label = label)

for i, lat in enumerate(LATITUDES):
       ax = axes[i]
       fig = figs[i]
       plt.figure(fig.number)

       ax.set_ylim([ymin,ymax])
       ax.set_yticks(np.arange(ymin,ymax+dy,dy), fontsize=axissize)
       ax.set_xlim([start_time,end_time])
       ax.set_xticks(np.arange(start_time,end_time+dt,dt), fontsize=axissize)
       ax.tick_params(axis='both', labelsize=axissize)
       latt = int(lat)
       hem = 'N' if latt >= 0 else 'S'
       ax.set_title(rf"AMOC at {abs(latt)}$^\circ${hem}", loc='center', size=titlesize)
       ax.set_ylabel('[Sv]', fontsize=titlesize)
       ax.set_xlabel("Time [Years]",fontsize=titlesize)
       ax.legend(loc='lower left', ncol=4, fontsize=legendsize, frameon=True)
       ax.grid(True, alpha=0.3)

       if PRINTING:
          plt.savefig(DIR_OUT / f'AMOC_timeseries_{abs(latt)}{hem}_rho.png',bbox_inches='tight',dpi=300)

#plt.show()



