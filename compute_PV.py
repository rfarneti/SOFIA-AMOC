'''
   Compute ocean potential vorticity 

   PV = (f + zeta) / rho * d(rho)/d(z) 
      = (f + zeta) * N2 / g
      =  f * N2 / g

   - calculate f
   - calculate zeta (Relative Vorticity)
   - calculate alpha
   - calculate beta
   - calculate N2

'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
#import glob
import momsofia as ms
#from scipy.interpolate import interp1d
#import cftime
#import datetime
#import nc_time_axis
#import netCDF4 as nc
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colorbar
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from pathlib import Path
import my_style

BASE_PATH = Path('/Volumes/Antwater2/DATA/')
DIR_OUT = Path('./Figures_N2_PV/')
PRINTING = True
# check all calculations by plotting
check = False


EXPT = ['0.0Sv','0.1Sv','1.0Sv','1.0Sv_60S']
BASIN = 'Atlantic'
DECADE = '41'

Y_LIMS = (-75, -30)       #meridional limits
D_LIMS = (100, 2500)     #vertical limits
X_LIMS_PAC = (-180, -90) #zonal limits for Pacific
X_LIMS_ATL = (-50, 0)    #zonal limts for Atlantic 


# Create dictionaries for all experiments
psi_args = {
    "0.0Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv",
        "u": "u.nc",
        "v": "v.nc",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": -1,"grid_xt_ocean": -1,}
    },
    "0.1Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.1Sv",
        "u": "u.nc",
        "v": "v.nc",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": -1,"grid_xt_ocean": -1,}
    },
    "1.0Sv": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv",
        "u": "u.nc",
        "v": "v.nc",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": -1,"grid_xt_ocean": -1,}
    },
    "1.0Sv_60S": {
        "dirin": "FAFANTWATER_ideal_runoff_Boeira_rest_ice_1.0Sv_60S",
        "u": "u.nc",
        "v": "v.nc",
        "t": "thetao.nc",
        "s": "so.nc",
        "p": "press.nc",
        "rho": "rho.nc",
        "potrho0": "pot_rho_0.nc",
        "potrho2": "pot_rho_2.nc",
        "chunks": {"time": 1,"potrho": -1,"grid_yu_ocean": -1,"grid_xt_ocean": -1,}
    },
}


if DECADE == '41':
   years = '(Years 41-50)'
   decade = 'year_41'
   year1 = '1641' ; year2   = '1650'
elif DECADE == '91':
   years = '(Years 91-100)'
   decade = 'year_91'
   year1 = '1691' ; year2   = '1700'
elif DECADE == '141':
   years = '(Years 141-150)'
   decade = 'year_141'
   year1 = '1741' ; year2   = '1750'
else:
   print('Wrong selection of decade')
   exit()
start_time = year1 + '-07-02'
end_time   = year2 + '-07-02' 

# Start main loop through experiments
# ===================================
for i in range(len(EXPT)):
    print('---------- Evaluating expt:', EXPT[i])
    expt = EXPT[i]

    #start_time = psi_args[expt[i]]["start_time"]
    #end_time   = psi_args[expt[i]]["end_time"]

    dirin   = psi_args[expt]["dirin"]
    u       = psi_args[expt]["u"]
    v       = psi_args[expt]["v"]
    t       = psi_args[expt]["t"]
    s       = psi_args[expt]["s"]
    p       = psi_args[expt]["p"]
    rho     = psi_args[expt]["rho"]
    potrho0 = psi_args[expt]["potrho0"]
    potrho2 = psi_args[expt]["potrho2"]

    # Get the variables
    u      = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / f'{u}')['u']
    v      = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / f'{v}')['v']
    thetao = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / f'{t}')['temp']
    so     = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / f'{s}')['salt']
    press  = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / f'{p}')['press']
    rho    = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / f'{rho}')['rho']
    pot_rho_0 = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / f'{potrho0}')['pot_rho_0']
    pot_rho_2 = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / f'{potrho2}')['pot_rho_2']

    u      = u.sel(time=slice(start_time, end_time))
    v      = v.sel(time=slice(start_time, end_time))
    thetao = thetao.sel(time=slice(start_time, end_time))
    so     = so.sel(time=slice(start_time, end_time))
    press  = press.sel(time=slice(start_time, end_time))
    ho    = rho.sel(time=slice(start_time, end_time))
    pot_rho_0 = pot_rho_0.sel(time=slice(start_time, end_time))
    pot_rho_2 = pot_rho_2.sel(time=slice(start_time, end_time))

    # Get some static fields
    dirinctl = psi_args[EXPT[0]]["dirin"]
    dxu = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / 'dxu.nc')['dxu']
    dyt = xr.open_dataset(BASE_PATH / f'{dirin}' / 'pp' / 'dyt.nc')['dyt']
    lat = xr.open_dataset(BASE_PATH / f'{dirinctl}' / 'history' / '17010101.ocean_grid.nc')["geolat_c"]
    #dzt = xr.open_dataset(pathin+dirin+'/pp/'+'dht.nc')['dht']
    #area_t = xr.open_dataset(pathin+dirin+'/pp/'+'areacello_t.nc')['area_t']
    #area_u = xr.open_dataset(pathin+dirin+'/pp/'+'areacello_u.nc')['area_u']
    #lon = xr.open_dataset(pathin+psi_args[expt[0]]["dirin"]+'/history/'+'17010101.ocean_grid.nc')["geolon_c"]



    # 1/ Compute Coriolis parameter 
    print('--- Computing Coriolis parameter')
    coriolis = ms.derived.calc_coriolis(lat)
    if check is True:
       ms.plotting.plot_coriolis(coriolis)

    # 2/ Compute Relative Vorticity 
    print('--- Computing Relative Vorticity')
    # Merge all necessary data first
    ds_zeta = xr.merge([u, v, dxu, dyt], compat="override")
    zeta = ms.derived.calc_rel_vort_MOM5(ds_zeta)
    if check is True:
       ms.plotting.plot_relvort(zeta)

    # 3/ Compute N^2 (alpha and beta)
    print('--- Computing alpha, beta, and N^2')
    n2, alpha, beta = ms.derived.calc_n2(thetao, so, eos="Wright", gravity=-9.8,
                      patm=101325.0, zcoord="st_ocean", interfaces=None,
                      adjust_negative=True)
    if check is True:
       ms.plotting.plot_alpha_beta(alpha, beta)

       ms.plotting.plot_n2(n2.sel(yt_ocean=slice(-80,0)).sel(st_ocean=slice(0,2200)).mean('xt_ocean').mean('time'), expt)
       plt.savefig('{dirout}{nameplot}_{expt}.png'.format(dirout=dirout,
                   nameplot='N2', expt=expt), bbox_inches='tight',dpi=300)

    # 4/ Compute PV
    print('--- Computing PV')
    PV = ms.derived.calc_pv(zeta, coriolis, n2, gravity=9.8,           
         coord_dict={"xcenter": "xu_ocean","ycenter": "yu_ocean",
                     "xcorner": "xt_ocean","ycorner": "yt_ocean"},
         symmetric=False, units="m", interp_f=True, full=False)



    if expt != '0.0Sv':
       print('--- Computing sigma_2 also for the Control')
       dirinctl    = psi_args[EXPT[0]]["dirin"]
       potrho_ctl = psi_args[EXPT[0]]["potrho2"]
       pot_rho_ctl = xr.open_dataset(BASE_PATH / dirinctl / 'pp' / potrho_ctl)['pot_rho_2']
       pot_rho_ctl = pot_rho_ctl.sel(time=slice(start_time, end_time))
       if BASIN == 'Pacific':
          sigma_ctl = pot_rho_ctl.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS)).sel(xt_ocean=slice(*X_LIMS_PAC))
       elif BASIN == 'Atlantic':
          sigma_ctl = pot_rho_ctl.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS)).sel(xt_ocean=slice(*X_LIMS_ATL))
       elif BASIN == 'Global':
          sigma_ctl = pot_rho_ctl.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS))
    else:
       sigma_ctl = 1000. #dummy value   

    # Plot the PV in the Global, Pacific or Atlantic sectors
    # ======================================================
    if BASIN == 'Pacific':
       PV =           PV.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS)).sel(xt_ocean=slice(*X_LIMS_PAC))
       sigma = pot_rho_2.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS)).sel(xt_ocean=slice(*X_LIMS_PAC))
    elif BASIN == 'Atlantic':
       PV =           PV.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS)).sel(xt_ocean=slice(*X_LIMS_ATL))
       sigma = pot_rho_2.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS)).sel(xt_ocean=slice(*X_LIMS_ATL))
    elif BASIN == 'Global':
       PV =           PV.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS))
       sigma = pot_rho_2.sel(yt_ocean=slice(*Y_LIMS)).sel(st_ocean=slice(*D_LIMS))
    else:
       print('basin Sector not present!')
       exit()
    
    factor = 1e11
    PV *= -1.
    
    fig = plt.figure(1, figsize=(8,4)) 
    ms.plotting.plot_zonalmean_PV(PV * factor, sigma-1000., sigma_ctl-1000., BASIN , expt, years, fig)
    if PRINTING:
       plt.savefig(DIR_OUT / f'PV_{BASIN}_{expt}_{decade}.png', bbox_inches='tight',dpi=300)
   

    # Plot the meridional gradient of PV (PVG) 
    PVG = PV.differentiate('yt_ocean', edge_order=2)
    fig = plt.figure(2, figsize=(8,4))
    ms.plotting.plot_zonalmean_PVG(PVG * factor, sigma-1000., sigma_ctl-1000, BASIN , expt, years, fig)
    if PRINTING:
       plt.savefig(DIR_OUT / f'PVG_{BASIN}_{expt}_{decade}.png', bbox_inches='tight',dpi=300)
   


    plt.show()




# Saves to a NetCDF file
#dirin = psi_args[expt[0]]["dirin"]
#path = pathin+dirin+'/pp/'+'quasi_bsf.nc'

#ds = psiu.rename("quasi_bsf")
#ds.sel(time=slice(start_time, end_time)).cf.mean("time").to_netcdf(path=path, 
#       mode="w", format="NetCDF4", compute=True)
#psiu.close()


