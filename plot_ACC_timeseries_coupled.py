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
DIR_OUT = Path('./Figures_transports/')


BASE_PATH = Path('/Volumes/Antwater2/DATA/ACC/')

MODELS  = [
        ("ACCESS-ESM12-5",   'blue',   'solid','1.'), 
        ("AWI-ESM-1-REcoM",  'salmon', 'solid','1.'),
        ("CESM2",            'red', 'solid','1.'),
        ("CanESM5",          'brown',  'solid','1.'),
        ("FOCI",             'black',  'solid','1.'),
        ("GFDL-CM4",         'blueviolet','solid','2.'),
        ("GFDL-ESM4",        'indigo',    'solid','2.'),
        ("GISS-E2-1-G",      'aquamarine','solid','1.'),
        ("HadGEM3-GC3.1-LL", 'cyan',      'solid','1.'),
        ("NorESM2-MM",       'lime',      'solid','1.'),
        ("MOM5",             'purple',    'solid','2.')
] 

# --- Load the data set
ds = xr.open_dataset(BASE_PATH / f'ACC_ant_ctl.nc', decode_times = False)
ds_r = ds.rolling(year=4, center=True).mean()

decade = ds.ACC_ant[:,91:100].mean('year') - ds.ACC_ctl[:,91:100].mean('year')
print('anomalies = ',decade)


titlesize  = 18
axissize   = 16
legendsize = 10

#offset = 1600; year = 365
#ds['TIME'] = ds.TIME/year - offset

# --- Plotting
rc('figure', figsize=(12,6))
fig, ax1 = plt.subplots()
y1, y2, dy = -50, 10, 10
t1, t2, dt1, dt2 = 0, 100, 10, 20

for i, (exp_id, color, lstyle, width) in enumerate(MODELS[:-1]):
    ax1.plot(ds.year, ds.ACC_ant[i,:] - ds.ACC_ctl[i,:], 
             color=color, linestyle=lstyle, linewidth=0.2)
    ax1.plot(ds_r.year, ds_r.ACC_ant[i,:] - ds_r.ACC_ctl[i,:], 
             color=color, linestyle=lstyle, linewidth=width, label = exp_id)
# MOM5
ax1.plot(ds.year, ds.ACC_ant[10,:] - ds.ACC_ctl[10,:], 
         color='purple', linestyle='dashed', linewidth=2.0, label = 'MOM5')

ax1.set_xlim([t1,t2])
ax1.set_xticks(np.arange(t1,t2+dt1,dt1))

ax1.set_ylim([y1,y2])
ax1.set_yticks(np.arange(y1,y2+dy,dy))
ax1.set_ylabel("[Sv]",fontsize=titlesize)
ax1.set_xlabel("Time [Years]",fontsize=titlesize)
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
plt.tick_params(axis='y', which='both', direction='in', 
                left=True, right=True, labelleft=False, labelright=True)
ax1.tick_params(labelsize=axissize)
ax1.legend(loc='lower left', ncol=3, fontsize=legendsize, frameon=True)
ax1.grid(True, alpha=0.3)

if PRINTING:
   plt.savefig(DIR_OUT / f'ACC_anomalies_coupled.png', bbox_inches='tight',dpi=600)

plt.show()
