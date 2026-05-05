'''
Plots SST and SSS time series

'''
import numpy as np
import xarray as xr
import cf_xarray as cfxr
import cftime
import datetime
import nc_time_axis
#import netCDF4 as nc
import matplotlib.pyplot as plt
import cmocean
from pathlib import Path
import my_style_latex

PRINTING = True
DIR_OUT = Path('./Figures_tracers/')
BASE_PATH = Path('/Volumes/Antwater2/DATA/')
PREFIX = 'FAFANTWATER_ideal_runoff_Boeira_rest_ice_'


offset = 1600
year = 365

EXPERIMENTS = [
    ("0.0Sv",     "Control",        "black",      "solid", 1.0),  # Control
    ("0.1Sv",     "0.1Sv",      "black",      "solid", 2.0),
    ("0.1Sv_tas_corrected",     "0.1Sv$^{*}$",      "black",      "dashed", 2.0),
    #("0.2Sv",     "0.2Sv",      "cyan",       "solid", 0.5),
    #("0.3Sv",     "0.3Sv",      "lime",       "solid", 0.5),
    #("0.4Sv",     "0.4Sv",      "green",      "solid", 0.5),
    #("0.5Sv",     "0.5Sv",      "royalblue",  "solid", 0.5),
    #("0.6Sv",     "0.6Sv",      "blue",       "solid", 0.5),
    #("0.7Sv",     "0.7Sv",      "blueviolet", "solid", 0.5),
    #("0.8Sv",     "0.8Sv",      "violet",     "solid", 0.5),
    #("0.9Sv",     "0.9Sv",      "magenta",    "solid", 0.5),
    #("1.0Sv",     "1.0Sv",      "purple",     "solid", 0.5),
]

YLIM = (-60,-45)

fontaxis=14
fonttitle=18
fontlegend=12
plt.rc('figure', figsize=(8,10))



#SST
fig1 = plt.figure()
ax1 = fig1.add_subplot(2,1,1)
for exp_id, label, color, style, width in EXPERIMENTS:
    ds_t   = xr.open_dataset(BASE_PATH / f"{PREFIX}{exp_id}" / 'pp' / 'tos.nc', decode_times = False)
    tos    = ds_t.sel(yt=[*YLIM], method='nearest').cf.mean(dim=['latitude','longitude']).compute()
    ax1.plot(tos.time/year-offset, tos.SST - tos.SST[0], color=color, linestyle=style, linewidth=width, label=label)
ax1.grid(alpha=0.3)
ax1.set_ylabel("[ $^{\circ}$C ]", fontsize=fonttitle)
v = [0,150,-0.6,0.1]
plt.yticks(np.arange(-0.6,0.1,0.1), fontsize=fontaxis)
plt.xticks(np.arange(0,160,10), fontsize=fontaxis)
ax1.axis(v)
#ax1.set_xlabel("Time [Years]",fontsize=fonttitle)
ax1.legend(loc='lower right', ncol=1, fontsize=fontlegend)
ax1.set_title(r'a) $\delta$ SST', fontsize=fonttitle, loc='left')

#SSS
ax2 = fig1.add_subplot(2,1,2)
for exp_id, label, color, style, width in EXPERIMENTS:
    ds_s   = xr.open_dataset(BASE_PATH / f"{PREFIX}{exp_id}" / 'pp' / 'sos.nc', decode_times = False)
    sos    = ds_s.sel(yt=[*YLIM], method='nearest').cf.mean(dim=['latitude','longitude']).compute()
    ax2.plot(sos.time/year-offset, sos.SSS - sos.SSS[0], color=color, linestyle=style, linewidth=width, label=label)
ax2.grid(alpha=0.3)
ax2.set_ylabel("[ g~Kg$^{-1}$ ] ", fontsize=fonttitle)
plt.xticks(np.arange(0,160,20), fontsize=fontaxis)
v = [0,150,-0.5,0.1]
plt.yticks(np.arange(-0.5,0.2,0.1), fontsize=fontaxis)
plt.xticks(np.arange(0,160,10), fontsize=fontaxis)
ax2.axis(v)
ax2.set_xlabel("Time [Years]",fontsize=fonttitle)
ax2.set_title(r'b) $\delta$ SSS', fontsize=fonttitle, loc='left')

if PRINTING:
   plt.savefig(DIR_OUT / f'SST_SSS_timeseries.png', bbox_inches='tight',dpi=300)


plt.show()
