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
    ("0.0Sv",     "Control",        "black",      "solid", 1.0),  # Control
    ("0.01Sv",    "0.01Sv",     "salmon",     "solid", 0.5),
    ("0.02Sv",    "0.02Sv",     "tomato",     "solid", 0.5),
    ("0.05Sv",    "0.05Sv",     "orangered",  "solid", 0.5),
    ("0.07Sv",    "0.07Sv",     "brown",      "solid", 0.5),
    ("0.1Sv",     "0.1Sv",      "black",      "solid", 3.0),
#    ("0.1Sv_60S", "0.1Sv-60S", "black",      "dashed",2.0),
    ("0.15Sv",    "0.15Sv",     "aquamarine", "solid", 0.5),
    ("0.2Sv",     "0.2Sv",      "cyan",       "solid", 0.5),
    ("0.3Sv",     "0.3Sv",      "lime",       "solid", 0.5),
    ("0.4Sv",     "0.4Sv",      "green",      "solid", 0.5),
    ("0.5Sv",     "0.5Sv",      "royalblue",  "solid", 0.5),
    ("0.6Sv",     "0.6Sv",      "blue",       "solid", 0.5),
    ("0.7Sv",     "0.7Sv",      "blueviolet", "solid", 0.5),
    ("0.8Sv",     "0.8Sv",      "violet",     "solid", 0.5),
    ("0.9Sv",     "0.9Sv",      "magenta",    "solid", 0.5),
    ("1.0Sv",     "1.0Sv",      "purple",     "solid", 3.0)
 #   ("1.0Sv_60S", "1.0Sv-60S", "purple",     "dashed",2.0)
]

LATITUDES = [-32,30,50]

OFFSET_YEAR = 1600; 
DAYS_IN_YEAR = 365

TIME_PERIODS = [
    {"label": "(years 01-10)", "start": "1601", "end": "1605", "ls": "-", "color": 'black'},
    {"label": "(years 21-30)", "start": "1621", "end": "1630", "ls": "-", "color": 'olive'},
    {"label": "(years 141-150)", "start": "1741", "end": "1750", "ls": "-", "color": 'purple'}
]

FORCING_LEVELS = ['0.1Sv', '1.0Sv']


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

def get_amoc_profile(ds, start_yr, end_yr, lat):
    """Slices and averages the AMOC data for a specific time and latitude."""
    t_start, t_end = f"{start_yr}-07-02", f"{end_yr}-07-02"
    # Average over time and select nearest latitude
    return ds.sel(TIME=slice(t_start, t_end)).mean(dim='TIME').sel(YU_OCEAN=lat, method='nearest')


# --- Main Processing and Plotting
rc('figure', figsize=(12,6))
start_time,end_time,dt = 0, 160, 10
ymin, ymax, dy = 4, 22, 2.
titlesize, axissize, legendsize  = 18, 14, 10

for i, lat in enumerate(LATITUDES):
    fig, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [4, 1]})
    ax1.tick_params(left=False, labelleft=False, right=True, labelright=True)
    ax0.tick_params(left=True, labelleft=True, right=False, labelright=False)
 
    #First plot the time series 
    for exp_id, label, color, style, width in EXPERIMENTS:
        amoc_index = compute_amoc_index(exp_id, lat)
        amoc_index.plot(ax=ax0, color=color, linestyle=style, linewidth=width, label = label)
        if lat == '30' or lat == '50' or lat == '40':
           amoc_index.to_netcdf(BASE_PATH / f'AMOC_index_{lat}N_{exp_id}.nc')
        else:
           amoc_index.to_netcdf(BASE_PATH / rf'AMOC_index_{abs(int(lat))}S_{exp_id}.nc')

    ax0.set_ylim([ymin,ymax])
    ax0.set_yticks(np.arange(ymin,ymax+dy,dy), fontsize=axissize)
    ax0.set_xlim([start_time,end_time])
    ax0.set_xticks(np.arange(start_time,end_time+dt,dt), fontsize=axissize)
    ax0.tick_params(labelsize=axissize)
    ax0.set_title("")
    latt = int(lat)
    hem = 'N' if latt >= 0 else 'S'  
    ax0.set_title('a)', fontweight='bold', loc='left', size=titlesize)
    ax0.set_title(rf"$\Psi_z$({abs(latt)}$^\circ${hem})", fontweight='bold', loc='center', size=titlesize)
    ax0.set_ylabel('[Sv]', fontsize=titlesize)
    ax0.set_xlabel("Time [Years]",fontsize=titlesize)
    ax0.legend(loc='lower left', ncol=4, fontsize=legendsize, frameon=True)
    ax0.grid(True, alpha=0.3)
    
    # shade regions of strengthening and collapse
    x1 = np.array([20, 30, 30, 20])
    x2 = np.array([140, 150, 150, 140])
    y1 = np.array([4,   4, 22, 22])
    ax0.fill(x1,y1, facecolor ='whitesmoke', edgecolor ='none', linewidth = 1) 
    ax0.fill(x2,y1, facecolor ='whitesmoke', edgecolor ='none', linewidth = 1) 

    #Now plot the profiles
    #plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
    #plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = False
    xmin, xmax, dx = -4.0, 20.0, 4.0
    forcing = '1.0Sv' 
    file_path = BASE_PATH / f'MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_{forcing}.nc'
    with xr.open_dataset(file_path) as ds:
        amoc_data = ds['AMOC']
        #Plot each time period
        for period in TIME_PERIODS:
            profile = get_amoc_profile(amoc_data, period['start'], period['end'], lat)
            profile.plot(ax=ax1, y="ST_OCEAN", yincrease=False,
                             color=period['color'], linestyle=period['ls'],
                             label=period['label'])

            #Add Depth of Max
            m_val, m_dep = profile.max(), profile.idxmax(dim='ST_OCEAN')
            depth = m_dep.values 
            ax1.plot([xmin, m_val], [m_dep, m_dep], '--o', color=period['color'], markevery=[-1], markersize=2)
            depth = m_dep.values 
            ax1.text(xmin+1, m_dep-20, f"{depth:.0f}", size=10, color=period['color'])

            ax1.set_xlabel("[Sv]", fontsize=titlesize)
            ax1.set_title(f"{forcing}", fontweight='bold', loc='center', size=titlesize)
            ax1.set_title('b)', fontweight='bold', loc='left', size=titlesize)
            # Grid and limits
            ax1.set_xlabel("[Sv]", fontsize=titlesize)
            ax1.set_ylabel("Depth [m]", fontsize=titlesize)
            ax1.yaxis.set_label_position("right")
            ax1.set_xlim(xmin, xmax)
            ax1.set_xticks(np.arange(xmin, xmax+dx, dx))
            ax1.set_ylim(5400, 0)
            ax1.grid(True, alpha=0.3)
            ax1.legend(fontsize=10, loc='lower right')
 
    plt.tight_layout() 
    if PRINTING:
       plt.savefig(DIR_OUT / f'AMOC_timeseries_{abs(latt)}{hem}_with_profiles.png', bbox_inches='tight',dpi=300)

plt.show()

