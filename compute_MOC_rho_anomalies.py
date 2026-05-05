'''
   Calculate the global Meridional Overturning Circulation 
   - in density space - using output from MOM5

   One can apply a mask for the Atlantic or IndoPacific basin

   indeces for the MOC are computed as:
   AMOC = max{ \Psi(\sigma)[45N] }
   AABW = max{ \Psi(\sigma)[50S:40S ; 37.25:36.8] }
'''
import numpy as np
import xarray as xr
import cf_xarray as cfxr
from scipy.interpolate import interp1d
#import cftime
#import datetime
#import nc_time_axis
#import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colorbar
import cmocean as cm
import cartopy.feature as cfeature
from ocean_basins import basin_mask_ht
from pathlib import Path
import momsofia as ms
import my_style_latex


#BASIN = 'IndoPacific'
BASIN = 'Atlantic'
#BASIN = 'Global'

EXPERIMENTS = [
         ("0.0Sv",     "a)", "d)"),
         ("0.1Sv",     "a)", "c)"),
         ("1.0Sv",     "b)", "d)")
]

BASE_PATH  = Path('/Volumes/Antwater2/DATA/')

PRINTING = True
DIR_OUT = Path('./Figures_transports/')


# create dictionaries for all experiments
psi_args = {
    "0.0Sv": {
        "dirin": "0.0Sv",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1691-07-02",
        "end_time"  : "1700-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 200,"grid_xt_ocean": -1,}
    },
    "0.1Sv": {
        "dirin": "0.1Sv",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1741-07-02",
        "end_time"  : "1750-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 200,"grid_xt_ocean": -1,}
    },
    "1.0Sv": {
        "dirin": "1.0Sv",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1741-07-02",
        "end_time"  : "1750-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 1080,"grid_xt_ocean": -1,}
    },
    "1.0Sv_60S": {
        "dirin": "1.0Sv_60S",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1741-07-02",
        "end_time"  : "1750-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 1080,"grid_xt_ocean": -1,}
    },
}

# --- Generate and later apply a mask, if needed
if BASIN != 'Global':
   print(f'---> Setting up mask for {BASIN}')
   dir_ctl = BASE_PATH / f'FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv' / 'pp'
   mask = ms.ocean_basins.basin_mask_ht(BASIN, SO=True, check=False)
   ds = xr.open_dataset(dir_ctl / 'ty_trans_rho.nc')
   mask = mask.rename({'xt_ocean': 'grid_xt_ocean', 'yt_ocean': 'grid_yu_ocean'})
   mask = mask.assign_coords(grid_xt_ocean=ds.grid_xt_ocean, grid_yu_ocean=ds.grid_yu_ocean)
else:
   print(f'---> Basin is {BASIN} and I am not applying a mask')


def compute_psi_rho(exp, mask=None):
    """Load ty_trans_rho - and sum zonally. Also, if there is a ty_trans_rho_gm variable saved, 
    assume that GM is switched on and load that as well. 
    Most MOM simulations save transport with units of kg/s - convert to Sv."""
    dirin  = BASE_PATH / f'FAFANTWATER_ideal_runoff_Boeira_rest_ice_{exp}' / 'pp'
    file1  = psi_args[exp]["file1"]
    file2  = psi_args[exp]["file2"]
    s_time = psi_args[exp]["start_time"]
    e_time = psi_args[exp]["end_time"]

    ds_psi    = xr.open_dataset(dirin / 'ty_trans_rho.nc')
    ds_psi_gm = xr.open_dataset(dirin / 'ty_trans_rho_gm.nc')
    psi = ds_psi['ty_trans_rho'].sel(time=slice(s_time, e_time))
    psi_gm = ds_psi_gm['ty_trans_rho_gm'].sel(time=slice(s_time, e_time))

    if mask is not None:
       psi = psi * mask
       psi_gm = psi_gm * mask
    psi = psi.cf.sum("longitude")
    psi_gm = psi_gm.cf.sum("longitude")
    """This is equivalent to the Ferret computation
       psi[i=@sum,k=@rsum]-psi[i=@sum,k=@sum]
       but I'm not sure why here I don't need
       to remove the vertical sum""" 
    psi = psi.cf.cumsum("vertical")  ##- psi.cf.sum("vertical")
    psi = psi.cf.mean("time") + psi_gm.mean("time")

    psi.load()

    return psi


def levels_and_colorbarticks(max_value,delta):
    """ Return the levels and the colorbarticks for the streamfunction plot.
    It may seem complicated but the truth is we just want to avoid the 0 contour
    so that the plot looks soothing to the eye"""
    levels =  np.hstack((np.arange(-max_value, 0, delta), np.flip(-np.arange(-max_value, 0, delta))))
    cbarticks = np.hstack((np.flip(-np.arange(delta*2, max_value+(delta*2), delta*2)), np.arange(delta*2, max_value+(delta*2), delta*2)))    
    return levels, cbarticks


labelsize = 14
titlesize = 16
max_psi = 20 # Sv
delta_psi = 2
cmap=cm.cm.delta

# --- Plotting the MOC for the Control
fig = plt.figure(1, figsize=(8,8))

psi_c = compute_psi_rho("0.0Sv", mask=mask if BASIN !='Global' else None)

# --- Plotting changes in MOC
for ii, (exp_id, number_atl, number_pac) in enumerate(EXPERIMENTS[1:]):
    ax = fig.add_subplot(2,1,ii+1)
    psi_s = compute_psi_rho(exp_id, mask=mask if BASIN !='Global' else None)
    diff = psi_s - psi_c

    levels, cbarticks = levels_and_colorbarticks(max_psi,delta_psi)

    diff.plot.contourf(ax=ax, levels=levels, cmap=cmap, extend='both',add_colorbar=False)
                  #cbar_kwargs={'pad':0.02, 'shrink': 0.9, 'label': '[Sv]', 'ticks': cbarticks2})
    psi_c.plot.contour(ax=ax, levels=levels, colors='k', linewidths=0.20)
    
    ax.invert_yaxis()
    ax.set_ylim((1037.2, 1036))
    ax.set_yticklabels([36.0,36.2,36.4,36.6,36.8,37.0,37.2])
    ax.set_xlim([-80, 80])
    ax.set_xticks(np.arange(-80,81,20))
    ax.set_xticklabels(('80S','60S','40S','20S','0','20N','40N','60N','80N'))
    ax.set_ylabel(r'$\sigma_2$ [kg m$^{-3}$]', fontsize=labelsize)
    if ii==1:
       ax.set_xlabel('Latitude [$^\circ$]', fontsize=labelsize) 
    else:
       ax.set_xlabel('')
       ax.set_xticklabels([])
    if BASIN == 'Atlantic': 
       ax.set_title(f'{number_atl}', loc='left', size=titlesize);
    else:
       ax.set_title(f'{number_pac}', loc='left', size=titlesize);
    ax.set_title(f'{BASIN}', loc='center', size=titlesize);
    ax.set_title(f'({exp_id})', loc='right', size=labelsize);

norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
ax = plt.gcf().add_axes((.2,.02,.6,.01))
cb1 = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels,
                                      orientation='horizontal', extend='both',)
cb1.ax.set_title('[Sv]', size=10)
cb1.set_ticks(cbarticks)
#plt.tight_layout()

if PRINTING:
   plt.savefig(DIR_OUT / f'MOC_rho_anomalies_{BASIN}.png',bbox_inches='tight',dpi=300)

plt.show()



