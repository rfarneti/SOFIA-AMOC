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
    ("0.1Sv",     "0.1Sv",       "black",    "solid",  4.0, 0.5),
    ("0.1Sv_60S", "0.1Sv-60S",   "black",    "dashed", 2.0, 0.5),
    ("0.1Sv_ambe", "0.1Sv-ambe", "black",    "solid",  2.0, 0.5),
    ("1.0Sv",     "1.0Sv",       "purple",   "solid",  4.0, 0.5),
    ("1.0Sv_60S", "1.0Sv-60S",   "purple",   "dashed", 2.0, 0.5),
    ("1.0Sv_ambe", "1.0Sv-ambe", "purple",   "solid",  2.0, 0.5),
    ("0.0Sv",     "CTL",         "black",    "solid",  1.0, 0.5)  # Control
]
EXPERIMENTS_SENS = [
    ("0.1Sv",     "0.1Sv",       "black",   "solid",  4.0, 0.5),
    ("0.1Sv_60S", "0.1Sv-60S",   "black",   "dashed", 2.0, 0.5),
    ("0.1Sv_ambe", "0.1Sv-ambe", "black",   "solid",  2.0, 0.5),
    ("1.0Sv",     "1.0Sv",       "purple",  "solid",  4.0, 0.5),
    ("1.0Sv_60S", "1.0Sv-60S",   "purple",  "dashed", 2.0, 0.5),
    ("1.0Sv_ambe", "1.0Sv-ambe", "purple",  "solid",  2.0, 0.5)
]

OFFSET_YEAR = 1500;
DAYS_IN_YEAR = 365


def compute_acc_index(exp_id):
    prefix = "ACC_FAFANTWATER_ideal_runoff_"
    suffix = "pv60_rest_ice_0.0Sv.nc" if exp_id == "0.0Sv" else f"Boeira_rest_ice_{exp_id}.nc"
    file_path = BASE_PATH / f"{prefix}{suffix}"
    with xr.open_dataset(file_path, decode_times = False) as ds:
         acc_index = ds.ACC
         time_years = ds.TIME / DAYS_IN_YEAR - OFFSET_YEAR
         return acc_index.assign_coords(TIME=time_years)

# --- Plotting
rc('figure', figsize=(10,6))
fig, ax1 = plt.subplots()
y1, y2, dy = 85, 135, 5
t1, t2, dt = 0, 160, 10
titlesize, axissize, legendsize = 18, 14, 10

prefix = "ACC_FAFANTWATER_ideal_runoff_"
for i, (exp_id, label, color, style, width, width2) in enumerate(EXPERIMENTS):
    acc_index = compute_acc_index(exp_id)
    if exp_id == "0.0Sv":
       ax1.plot(acc_index, color=color, linestyle=style, linewidth=width, label = None)
    else:
       ax1.plot(acc_index, color=color, linestyle=style, linewidth=width, label = label)
ax1.set_ylim([y1,y2])
ax1.set_yticks(np.arange(y1,y2+dy,dy))
ax1.set_xlim([t1,t2])
ax1.set_xticks(np.arange(t1,t2+dt,dt), fontsize=axissize)
ax1.tick_params(labelsize=axissize)
ax1.set_title('b)', fontweight='bold', loc='left', size=titlesize)
ax1.set_title(rf"ACC mass transport", fontweight='bold', loc='center', size=titlesize)
ax1.set_ylabel("[Sv]",fontsize=titlesize)
ax1.set_xlabel("Time [Years]",fontsize=titlesize)
ax1.legend(loc='lower left', ncol=2, fontsize=legendsize, frameon=True)
ax1.grid(True, alpha=0.3)

if PRINTING:
   plt.savefig(DIR_OUT / f'ACC_timeseries_60S_ambe.png', bbox_inches='tight',dpi=600)

# --- Anomalies
rc('figure', figsize=(10,6))
fig, ax1 = plt.subplots()
y1, y2, dy = -45, 5, 5
t1, t2, dt = 0, 160, 10
titlesize, axissize, legendsize = 18, 14, 10

prefix = "ACC_FAFANTWATER_ideal_runoff_"

acc_index_ctl = compute_acc_index("0.0Sv")

for i, (exp_id, label, color, style, width, width2) in enumerate(EXPERIMENTS_SENS):
    acc_index = compute_acc_index(exp_id)
    ax1.plot(acc_index - acc_index_ctl, color=color, linestyle=style, linewidth=width, label = label)
ax1.set_ylim([y1,y2])
ax1.set_yticks(np.arange(y1,y2+dy,dy))
ax1.set_xlim([t1,t2])
ax1.set_xticks(np.arange(t1,t2+dt,dt), fontsize=axissize)
ax1.tick_params(labelsize=axissize)
ax1.set_title('b)', fontweight='bold', loc='left', size=titlesize)
ax1.set_title(rf"$\Delta$ ACC mass transport", fontweight='bold', loc='center', size=titlesize)
ax1.set_ylabel("[Sv]",fontsize=titlesize)
ax1.set_xlabel("Time [Years]",fontsize=titlesize)
ax1.legend(loc='lower left', ncol=2, fontsize=legendsize, frameon=True)
ax1.grid(True, alpha=0.3)
ax1.hlines(0,t1, t2, colors='black', linestyles='solid', linewidth=0.5)
# shade regions of strengthening and collapse
x1 = np.array([20, 30, 30, 20])
x2 = np.array([140, 150, 150, 140])
y1 = np.array([-45, -45, 5, 5])
ax1.fill(x1,y1, facecolor ='whitesmoke', edgecolor ='none', linewidth = 1)
ax1.fill(x2,y1, facecolor ='whitesmoke', edgecolor ='none', linewidth = 1)


if PRINTING:
   plt.savefig(DIR_OUT / f'ACC_timeseries_60S_ambe_anomalies.png', bbox_inches='tight',dpi=600)

plt.show()
