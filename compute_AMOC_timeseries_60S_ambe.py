'''
 Plots AMOC time series at different latitudes:
       30N, 40N, 50N and 32S
'''
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import xarray as xr
import cftime
import datetime
import netCDF4 as nc
from pathlib import Path
import my_style_latex



BASE_PATH   = Path('/Volumes/Antwater2/DATA/MOC/')
DIR_OUT   = Path('./Figures_transports/')
PRINTING = True


EXPERIMENTS = [
    ("0.1Sv",     "0.1Sv",       "black",    "solid",  4.0),
    ("0.1Sv_60S", "0.1Sv-60S",   "black",    "dashed", 2.0),
    ("0.1Sv_ambe", "0.1Sv-ambe", "black",    "solid",  2.0),
    ("1.0Sv",     "1.0Sv",       "purple",   "solid",  4.0),
    ("1.0Sv_60S", "1.0Sv-60S",   "purple",   "dashed", 2.0),
    ("1.0Sv_ambe","1.0Sv-ambe",  "purple",   "solid",  2.0),
    ("0.0Sv",     "CTL",         "black",    "solid",  1.0)  # Control
]
EXPERIMENTS_SENS = [
    ("0.1Sv",     "0.1Sv",       "black",    "solid",  4.0),
    ("0.1Sv_60S", "0.1Sv-60S",   "black",    "dashed", 2.0),
    ("0.1Sv_ambe","0.1Sv-ambe",  "black",    "solid",  2.0),
    ("1.0Sv",     "1.0Sv",       "purple",   "solid",  4.0),
    ("1.0Sv_60S", "1.0Sv-60S",   "purple",   "dashed", 2.0),
    ("1.0Sv_ambe","1.0Sv-ambe",  "purple",   "solid",  2.0)
]  

#LATITUDE = ['-32','30','50']
LATITUDE = ['30']

OFFSET_YEAR = 1600; 
DAYS_IN_YEAR = 365


def compute_amoc_index(exp_id, lat, z_range=(0, 4000)):
    """Computes the vertical maximum of the AMOC at a specific latitude over all time."""
    """Loads a single experiment's AMOC index at a given latitude."""
    prefix = "MOC_FAFANTWATER_ideal_runoff_"
    suffix = "pv60_rest_ice_0.0Sv.nc" if exp_id == "0.0Sv" else f"Boeira_rest_ice_{exp_id}.nc"
    file_path = BASE_PATH / f"{prefix}{suffix}"
    
    with xr.open_dataset(file_path, decode_times = False) as ds:
         amoc_slice = ds['AMOC'].sel(YU_OCEAN=lat, method='nearest').sel(ST_OCEAN=slice(*z_range))
         amoc_index = amoc_slice.max(dim="ST_OCEAN")

         time_years = amoc_index.TIME / DAYS_IN_YEAR - OFFSET_YEAR
         return amoc_index.assign_coords(TIME=time_years)


# --- Main Processing and Plotting
rc('figure', figsize=(10,6))
start_time,end_time,dt = 0, 160, 10
ymin, ymax, dy = 4, 22, 2.
titlesize, axissize, legendsize  = 18, 14, 10
for i, lat in enumerate(LATITUDE):
    fig = plt.figure(i+1)
    ax = fig.add_subplot(1,1,1)

    for exp_id, label, color, style, width in EXPERIMENTS:
        amoc_index = compute_amoc_index(exp_id, lat)
        if exp_id == "0.0Sv":
           amoc_index.plot(ax=ax, color=color, linestyle=style, linewidth=width, label = None)
        else:
           amoc_index.plot(ax=ax, color=color, linestyle=style, linewidth=width, label = label)
        if lat == '30' or lat == '50' or lat == '40':
           amoc_index.to_netcdf(BASE_PATH / f'AMOC_index_{lat}N_{exp_id}.nc')
        else:
           amoc_index.to_netcdf(BASE_PATH / rf'AMOC_index_{abs(int(lat))}S_{exp_id}.nc')

    ax.set_ylim([ymin,ymax])
    ax.set_yticks(np.arange(ymin,ymax+dy,dy), fontsize=axissize)
    ax.set_xlim([start_time,end_time])
    ax.set_xticks(np.arange(start_time,end_time+dt,dt), fontsize=axissize)
    ax.tick_params(labelsize=axissize)
    ax.set_title("")
    latt = int(lat)
    hem = 'N' if latt >= 0 else 'S'  
    ax.set_title('a)', fontweight='bold', loc='left', size=titlesize)
    ax.set_title(rf"$\Psi_z$({abs(latt)}$^\circ${hem})", fontweight='bold', loc='center', size=titlesize)
    ax.set_ylabel('[Sv]', fontsize=titlesize)
    ax.set_xlabel("Time [Years]",fontsize=titlesize)
    ax.legend(loc='lower left', ncol=2, fontsize=legendsize, frameon=True)
    ax.grid(True, alpha=0.3)
    
    if PRINTING:
       plt.savefig(DIR_OUT / f'AMOC_timeseries_{abs(latt)}{hem}_60S_ambe.png', bbox_inches='tight',dpi=300)

# --- Main Processing and Plotting
rc('figure', figsize=(10,6))
start_time,end_time,dt = 0, 160, 10
ymin, ymax, dy = -12, 4, 2.
titlesize, axissize, legendsize  = 18, 14, 10
for i, lat in enumerate(LATITUDE):
    fig = plt.figure(i+4)
    ax = fig.add_subplot(1,1,1)

    amoc_index_ctl = compute_amoc_index("0.0Sv", lat)
    for exp_id, label, color, style, width in EXPERIMENTS_SENS:
        amoc_index = compute_amoc_index(exp_id, lat)
        anomaly = amoc_index - amoc_index_ctl 
        anomaly.plot(ax=ax, color=color, linestyle=style, linewidth=width, label = label)
        if lat == '30' or lat == '50' or lat == '40':
           amoc_index.to_netcdf(BASE_PATH / f'AMOC_index_{lat}N_{exp_id}.nc')
        else:
           amoc_index.to_netcdf(BASE_PATH / rf'AMOC_index_{abs(int(lat))}S_{exp_id}.nc')
    ax.hlines(0,start_time, end_time, colors='black', linestyles='solid', linewidth=0.5)
    ax.set_ylim([ymin,ymax])
    ax.set_yticks(np.arange(ymin,ymax+dy,dy), fontsize=axissize)
    ax.set_xlim([start_time,end_time])
    ax.set_xticks(np.arange(start_time,end_time+dt,dt), fontsize=axissize)
    ax.tick_params(labelsize=axissize)
    ax.set_title("")
    latt = int(lat)
    hem = 'N' if latt >= 0 else 'S'
    ax.set_title('a)', fontweight='bold', loc='left', size=titlesize)
    ax.set_title(rf"$\Delta \Psi_z$({abs(latt)}$^\circ${hem})", fontweight='bold', 
                 loc='center', size=titlesize)
    ax.set_ylabel('[Sv]', fontsize=titlesize)
    ax.set_xlabel("Time [Years]",fontsize=titlesize)
    ax.legend(loc='lower left', ncol=2, fontsize=legendsize, frameon=True)
    ax.grid(True, alpha=0.3)
    # shade regions of strengthening and collapse
    x1 = np.array([20, 30, 30, 20])
    x2 = np.array([140, 150, 150, 140])
    y1 = np.array([-12, -12, 4, 4])
    ax.fill(x1,y1, facecolor ='whitesmoke', edgecolor ='none', linewidth = 1)
    ax.fill(x2,y1, facecolor ='whitesmoke', edgecolor ='none', linewidth = 1)


    if PRINTING:
       plt.savefig(DIR_OUT / f'AMOC_timeseries_{abs(latt)}{hem}_60S_ambe_anomalies.png', bbox_inches='tight',dpi=300)

plt.show()

