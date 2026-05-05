'''
   Calculate the global Meridional Overturning Circulation 
   - in density space - using output from MOM5

   One can apply a mask for the Atlantic or IndoPacific basin

   indeces for the MOC are computed as:
   AMOC = max{ \Psi(\sigma)[30N:65N ; 35.00:37.00] }
   AABW = max{ \Psi(\sigma)[60S:35S ; 36.80:37.25] }

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
import my_style_latex


# Choose the basin: Global, Atlantic, or IndoPacific
#BASIN = 'IndoPacific'
BASIN = 'Atlantic'
#BASIN = 'Global'

EXPERIMENTS = [
         ("0.0Sv",     "a)", "d)"),
         ("0.1Sv",     "b)", "e)"),
         ("1.0Sv",     "c)", "f)")
]  

BASE_PATH   = '/Volumes/Antwater2/DATA/'
PRINTING = True
DIR_OUT = Path('./Figures_transports/')

d_range_amoc=(1035.0, 1037.0)
d_range_pmoc=(1035.0, 1037.0)
d_range_aabw=(1036.8, 1037.25)
lat_range_amoc = (30, 65)
lat_range_pmoc = (-32, 0)
lat_range_aabw = (-60, -35)

# create dictionaries for all experiments
psi_args = {
    "0.0Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1691-07-02",
        "end_time"  : "1700-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 200,"grid_xt_ocean": -1,}
    },
    "0.1Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.1Sv",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1741-07-02",
        "end_time"  : "1750-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 200,"grid_xt_ocean": -1,}
    },
    "1.0Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv",
        "file1": "ty_trans_rho.nc",
        "file2": "ty_trans_rho_gm.nc",
        "start_time": "1741-07-02",
        "end_time"  : "1750-07-02",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": 1080,"grid_xt_ocean": -1,}
    }
}


def compute_psi_rho(expt, mask=None):
    """Load ty_trans_rho - and sum zonally. Also, if there is a ty_trans_rho_gm variable saved, 
    assume that GM is switched on and load that as well. 
    Most MOM simulations save transport with units of kg/s - convert to Sv."""
    dirin      = psi_args[expt]["dirin"]
    file1      = psi_args[expt]["file1"]
    file2      = psi_args[expt]["file2"]
    start_time = psi_args[expt]["start_time"]
    end_time   = psi_args[expt]["end_time"]

    psi   = xr.open_dataset(BASE_PATH + dirin+'/pp/'+file1)['ty_trans_rho']
    psiGM = 0 * psi.copy(deep=True)
    psiGM = xr.open_dataset(BASE_PATH + dirin+'/pp/'+file2)['ty_trans_rho_gm']
    
    psi   = psi.sel(time=slice(start_time, end_time))
    psiGM = psiGM.sel(time=slice(start_time, end_time))
    #convert to volume transport
    #rho0 = 1025
    #psi = psi / rho0
    #psiGM = psiGM / rho0
    if mask is not None:
       psi = psi * mask
       psiGM = psiGM * mask
    psi = psi.cf.sum("longitude")
    psiGM = psiGM.cf.sum("longitude")
    """This is equivalent to the Ferret computation
       psi[i=@sum,k=@rsum]-psi[i=@sum,k=@sum]
       but I'm not sure why here I don't need
       to remove the vertical sum""" 
    psi = psi.cf.cumsum("vertical")  ##- psi.cf.sum("vertical")
    psi = psi.cf.mean("time") + psiGM.mean("time")
    psi.load()

    return psi



# Apply a mask if needed
if BASIN != 'Global':
   print('Basin is', BASIN)
   mask = basin_mask_ht(BASIN, SO=True)

   dirin = psi_args["0.0Sv"]["dirin"]
   file1 = psi_args["0.0Sv"]["file1"]
   tyrho = xr.open_dataset(BASE_PATH+dirin+'/pp/'+file1)['ty_trans_rho']
   mask.coords['xt_ocean'] = tyrho.grid_xt_ocean.values
   mask.coords['yt_ocean'] = tyrho.grid_yu_ocean.values
   mask = mask.rename({'xt_ocean': 'grid_xt_ocean', 'yt_ocean': 'grid_yu_ocean'})
   del tyrho
   print('I am appying a mask for', BASIN)
   print('mask shape', mask.shape)
else:
   print('Basin is', BASIN)
   print('I am not applying a mask')


def levels_and_colorbarticks(max_value, delta):
    """ Return the levels and the colorbarticks for the streamfunction plot.
    It may seem complicated but the truth is we just want to avoid the 0 contour
    so that the plot looks soothing to the eye"""
    levels =  np.hstack((np.arange(-max_value, 0, delta), np.flip(-np.arange(-max_value, 0, delta))))
    cbarticks = np.hstack((np.flip(-np.arange(delta*2, max_value+(delta*2), delta*2)), np.arange(delta*2, max_value+(delta*2), delta*2)))

    return levels, cbarticks


labelsize = 14
titlesize = 16
max_psi = 28 # Sv
delta_psi = 2
cmap=cm.cm.delta

fig = plt.figure(1, figsize=(8,10))
for ii, (exp_id, number_a, number_ip) in enumerate(EXPERIMENTS):
    if BASIN != 'Global':
       psi   = compute_psi_rho(exp_id, mask=mask)
    else:
       psi   = compute_psi_rho(exp_id)

    #Print values for AMOC and AABW
    amoc = psi.sel(grid_yu_ocean=slice(*lat_range_amoc)).sel(potrho=slice(*d_range_amoc))
    amoc_max = amoc.max(dim=["potrho","grid_yu_ocean"])
    aabw = psi.sel(grid_yu_ocean=slice(*lat_range_aabw)).sel(potrho=slice(*d_range_aabw))
    aabw_max = aabw.min(dim=["potrho","grid_yu_ocean"])
    pmoc = psi.sel(grid_yu_ocean=slice(*lat_range_pmoc)).sel(potrho=slice(*d_range_pmoc))
    pmoc_max = pmoc.max(dim=["potrho","grid_yu_ocean"])
    print('##===================================#')
    print(f'AMOC max for {exp_id} =', amoc_max)
    print(f'PMOC max for {exp_id} =', pmoc_max)
    print(f'AABW max for {exp_id} =', aabw_max)
    print('##===================================#')

    ax = fig.add_subplot(3,1,ii+1)

    levels, cbarticks = levels_and_colorbarticks(max_psi,delta_psi)
    #levels = np.arange(-max_psi,max_psi+delta_psi,delta_psi)

    psi.plot.contourf(levels=levels,
                      cmap=cmap, extend='both', add_colorbar=False)
                      #cbar_kwargs={'pad':0.02, 'shrink': 0.9, 'label': '[Sv]', 'ticks': cbarticks})
    psi.plot.contour(levels=levels, colors='k', linewidths=0.20)

    ax.invert_yaxis()
    ax.set_ylim((1037.2, 1032.0))
    ax.set_yticklabels([32.0,33.0,34.0,35.0,36.0,37.0])
    ax.set_xlim([-80, 80])
    ax.set_xticks(np.arange(-80,81,20))
    ax.set_xticklabels(('80S','60S','40S','20S','0','20N','40N','60N','80N'))
    ax.set_ylabel(r'$\sigma_2$ [kg m$^{-3}$]', fontsize=labelsize)
    if ii == 2: 
       ax.set_xlabel('Latitude [$^\circ$]', fontsize=labelsize)
    else:
       ax.set_xlabel('')
       ax.set_xticklabels([])
    if BASIN == 'Global' or BASIN == 'Atlantic':
       ax.set_title(f'{number_a}', loc='left', size=titlesize);
    else:
       ax.set_title(f'{number_ip}', loc='left', size=titlesize);
    ax.set_title(f'({exp_id})', loc='right', size=labelsize);
    ax.set_title(f'{BASIN}', loc='center', size=titlesize);

norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
ax = plt.gcf().add_axes((.2,.02,.6,.01))
cb1 = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels,
                                      orientation='horizontal', extend='both',)
cb1.ax.set_title('[Sv]', size=10)
cb1.set_ticks(cbarticks)
#cb1.ax.tick_params(labelsize=8)
#plt.tight_layout()

if PRINTING:
   plt.savefig(DIR_OUT / f'MOC_rho_{BASIN}.png',bbox_inches='tight',dpi=300)

plt.show()



