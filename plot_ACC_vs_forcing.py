'''
Plots ACC time series for all SOFIA runs
'''
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import xarray as xr
#from scipy import interpolate
from scipy.stats import linregress
import cftime
import datetime
import netCDF4 as nc
from pathlib import Path
import my_style_latex


PRINTING = True
BASE_PATH = Path('/Volumes/Antwater2/DATA/ACC/')
DIR_OUT = Path('./Figures_transports/')


prefix = "ACC_FAFANTWATER_ideal_runoff_"
EXPERIMENTS = [
    #("0.0Sv",     "CTL",        "black",      "solid", 1.0, 0.5),  # Control
    ("0.01Sv",    "0.01",     "salmon",     "solid", 0.5, 0.5),
    ("0.02Sv",    "0.02",     "tomato",     "solid", 0.5, 0.5),
    ("0.05Sv",    "0.05",     "orangered",  "solid", 0.5, 0.5),
    ("0.07Sv",    "0.07",     "brown",      "solid", 0.5, 0.5),
    ("0.1Sv",     "0.1",      "black",      "solid", 2.0, 2.0),
    ("0.15Sv",    "0.15",     "aquamarine", "solid", 0.5, 0.5),
    ("0.2Sv",     "0.2",      "cyan",       "solid", 0.5, 0.5),
    ("0.3Sv",     "0.3",      "lime",       "solid", 0.5, 0.5),
    ("0.4Sv",     "0.4",      "green",      "solid", 0.5, 0.5),
    ("0.5Sv",     "0.5",      "royalblue",  "solid", 0.5, 0.5),
    ("0.6Sv",     "0.6",      "blue",       "solid", 0.5, 0.5),
    ("0.7Sv",     "0.7",      "blueviolet", "solid", 0.5, 0.5),
    ("0.8Sv",     "0.8",      "violet",     "solid", 0.5, 0.5),
    ("0.9Sv",     "0.9",      "magenta",    "solid", 0.5, 0.5),
    ("1.0Sv",     "1.0",      "purple",     "solid", 2.0, 0.5)
]



# --- Plotting
y1, y2, dy = 100, 125, 5
titlesize  = 18
axissize   = 14

s_yr = '1621'
e_yr = '1640'
t_start, t_end = f"{s_yr}-07-02", f"{e_yr}-07-02"

fig, ax = plt.subplots(1, 1, figsize=(6,5))
for i, (exp_id, forcing, color, style, width, width2) in enumerate(EXPERIMENTS):
    suffix = f"Boeira_rest_ice_{exp_id}.nc"
    file_path = BASE_PATH / f"{prefix}{suffix}"
    ds = xr.open_dataset(file_path, decode_times = True) 
    ds = ds.sel(TIME=slice(t_start, t_end)).mean('TIME')
    ax.scatter(float(forcing), ds.ACC, s=80, facecolors=color, marker='o', edgecolors='white', alpha=1.0)
ax.grid(color='grey', linestyle='solid', linewidth=0.2)
ax.tick_params(axis='both', rotation=0, labelsize=axissize)
ax.set_xticks(np.arange(0.0,1.1,0.1))
ax.set_xlabel(r"$\mathcal{F}_w$ [Sv]",fontsize=titlesize)
ax.set_ylim([y1,y2])
ax.set_yticks(np.arange(y1,y2,dy))
ax.set_ylabel("[Sv]", fontsize=titlesize)
ax.set_title(f"ACC transport vs. Freshwater forcing", loc='center', size=16)

if PRINTING:
   plt.savefig(DIR_OUT / f'ACC_vs_forcing_scatterplot.png', bbox_inches='tight',dpi=600)

plt.show()
