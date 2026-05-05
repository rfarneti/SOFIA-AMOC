'''
Plots ACC time series for all SOFIA runs
'''
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import xarray as xr
from scipy import interpolate
import cftime
import datetime
import netCDF4 as nc
from pathlib import Path
import my_style_latex


PRINTING = True
BASE_PATH = Path('/Volumes/Antwater2/DATA/ACC/')
DIR_OUT = Path('./Figures_transports/')
nameplot   = ['ACC_timeseries','ACC_antwater_anomalies','ACC_antwater_scale']


EXPERIMENTS = [
    ("0.0Sv",     "Control",        "black",      "solid", 1.0, 0.5),  # Control
    ("0.01Sv",    "0.01Sv",     "salmon",     "solid", 0.5, 0.5),
    ("0.02Sv",    "0.02Sv",     "tomato",     "solid", 0.5, 0.5),
    ("0.05Sv",    "0.05Sv",     "orangered",  "solid", 0.5, 0.5),
    ("0.07Sv",    "0.07Sv",     "brown",      "solid", 0.5, 0.5),
    ("0.1Sv",     "0.1Sv",      "black",      "solid", 2.0, 2.0),
#    ("0.1Sv_60S", "0.1Sv-60S", "black",      "dashed",2.0, 0.5),
    ("0.1Sv_tas_corrected", "0.1Sv$^{*}$", "black",      "dashed",2.0, 0.5),
    ("0.15Sv",    "0.15Sv",     "aquamarine", "solid", 0.5, 0.5),
    ("0.2Sv",     "0.2Sv",      "cyan",       "solid", 0.5, 0.5),
    ("0.3Sv",     "0.3Sv",      "lime",       "solid", 0.5, 0.5),
    ("0.4Sv",     "0.4Sv",      "green",      "solid", 0.5, 0.5),
    ("0.5Sv",     "0.5Sv",      "royalblue",  "solid", 0.5, 0.5),
    ("0.6Sv",     "0.6Sv",      "blue",       "solid", 0.5, 0.5),
    ("0.7Sv",     "0.7Sv",      "blueviolet", "solid", 0.5, 0.5),
    ("0.8Sv",     "0.8Sv",      "violet",     "solid", 0.5, 0.5),
    ("0.9Sv",     "0.9Sv",      "magenta",    "solid", 0.5, 0.5),
    ("1.0Sv",     "1.0Sv",      "purple",     "solid", 2.0, 0.5)
 #   ("1.0Sv_60S", "1.0Sv-60S", "purple",     "dashed",2.0, 0.5)
]

OFFSET_YEAR = 1500;
DAYS_IN_YEAR = 365


# --- Plotting
rc('figure', figsize=(12,6))
fig, ax1 = plt.subplots()
y1, y2, dy = 85, 135, 5
t1, t2, dt1, dt2 = 0, 160, 10, 20
titlesize  = 18
axissize   = 16
legendsize = 10

prefix = "ACC_FAFANTWATER_ideal_runoff_"
for i, (exp_id, label, color, style, width, width2) in enumerate(EXPERIMENTS):
    suffix = "pv60_rest_ice_0.0Sv.nc" if exp_id=="0.0Sv" else f"Boeira_rest_ice_{exp_id}.nc"
    file_path = BASE_PATH / f"{prefix}{suffix}"
    ds = xr.open_dataset(file_path, decode_times = False) 
    time_years = ds.TIME / DAYS_IN_YEAR - OFFSET_YEAR
    ds.assign_coords(TIME=time_years)

    ax1.plot(ds.ACC, color=color, linestyle=style, linewidth=width, label = label)
ax1.set_ylim([y1,y2])
ax1.set_xlim([t1,t2])
ax1.set_yticks(np.arange(y1,y2+dy,dy))
ax1.tick_params(labelsize=axissize)
ax1.set_ylabel("[Sv]",fontsize=titlesize)
ax1.legend(loc='lower left', ncol=2, fontsize=8, frameon=True)
ax1.set_title('ACC mass transport', loc='center', fontsize=titlesize)

# shade regions of strengthening and collapse
x1 = np.array([20, 30, 30, 20])
x2 = np.array([140, 150, 150, 140])
y1 = np.array([85,   85, 135, 135])
ax1.fill(x1,y1, facecolor ='whitesmoke', edgecolor ='none', linewidth = 1)          
ax1.fill(x2,y1, facecolor ='whitesmoke', edgecolor ='none', linewidth = 1)

# --- A small inlet with percentage change
left, bottom, width, height = [0.65, 0.65, 0.25, 0.23]
ax2 = fig.add_axes([left, bottom, width, height], facecolor='w')
y1,y2,dy = -35, 0, 5

for i, (exp_id, label, color, style, width, width2) in enumerate(EXPERIMENTS):
    suffix = f"Boeira_rest_ice_{exp_id}.nc"
    file_path = BASE_PATH / f"{prefix}{suffix}"
    ds = xr.open_dataset(file_path, decode_times = False)
    time_years = ds.TIME / DAYS_IN_YEAR - OFFSET_YEAR
    ds.assign_coords(TIME=time_years)

    ax2.plot( (ds.ACC - ds.ACC[0])/ds.ACC[0] * 100,
        color=color, linestyle=style,  linewidth=width2, label = label)
ax2.set_ylim([y1,y2])
ax2.set_xlim([t1,t2])
ax2.set_yticks(np.arange(y1,y2+dy,dy))
ax2.tick_params(labelsize=10)
ax2.set_ylabel('$\%~\delta$ACC', fontsize=axissize)

## MMM from coupled SOFIA models: -14%
##  "ACCESS-ESM1-5","AWI-ESM-1-REcoM","CESM2","CanESM5","FOCI",  "GFDL-CM4",
##  -9.4502         -6.4747           -3.1800 -10.4841  -10.6215 -32.6640 
## "GFDL-ESM4", "GISS-E2.1-G","HadGEM3-GC31-LL","NorESM2-MM","MOM5"
##              -13.7461      -16.5294          -11.5571     -25.0080
##MMM = 14, yerr1 = 18.66, yerr2 = 10.82
acc_coupled_mom = (-9.5, -6.5, -3.2, -10.5, -10.6, -32.7, -16.5, -13.7, -16.5, -11.6, -25.)
acc_coupled     = (-9.5, -6.5, -3.2, -10.5, -10.6, -32.7, -16.5, -13.7, -16.5, -11.6)
mmm_c   = np.mean(acc_coupled)
mmm_c_m = np.mean(acc_coupled_mom)
mmmedian_c_m = np.median(acc_coupled_mom)
print('MEAN', mmm_c_m)
print('MEDIAN', mmmedian_c_m)


yerr2 = np.min(acc_coupled_mom) - mmm_c_m
yerr1 = mmm_c_m - np.max(acc_coupled_mom)
spread = [yerr1, yerr2]
ax2.errorbar(100, mmm_c_m, yerr=[[yerr1],[yerr2]], fmt='o', mfc='w', mec='k', ms=4, 
             capsize=2, capthick=1, ecolor='k', elinewidth=1)

for ax in [ax1, ax2]:
    ax.set_xlim([t1,t2])
    ax.grid(True, alpha=0.3)
    if ax==ax1:
       ax.set_xticks(np.arange(t1,t2+dt1,dt1))
       ax.set_xlabel("Time [Years]",fontsize=titlesize)
    if ax==ax2:
       #ax.grid(True, alpha=0.3)
       ax.set_xticks(np.arange(t1,t2,dt2))
       ax.yaxis.set_label_position("right")
       ax.yaxis.tick_right()
       ax.xaxis.tick_top()


if PRINTING:
   plt.savefig(DIR_OUT / f'ACC_timeseries.png', bbox_inches='tight',dpi=600)

plt.show()
