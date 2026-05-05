'''

 Plot Atlantic zonal-mean pot_rho_0 for different forcings

 One can choose 
 the basin: Atlantic, Pacific, Indian
 the field: so, thetao, pot_rho_0
 the forcing: 0.1Sv, ..., 1.0Sv

'''

import numpy as np
import xarray as xr
import cf_xarray as cfxr
import momsofia as ms
#import pint_xarray
import cftime
import datetime
import nc_time_axis 
#import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colorbar
import cmocean
from ocean_basins import basin_mask_ht_full
from pathlib import Path
import my_style_latex
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{xcolor}'

# Set path, plot name
BASE_PATH  = Path('./DATA/')
DIR_OUT = Path('./Figures_tracers/')
PRINTING = True

# Choose basin, density (in-situ, sigma0, sigma2), strength of forcing
OCEAN    = 'Atlantic'
#FORCING  = ['0.0Sv','0.1Sv','1.0Sv']
FORCING  = ['1.0Sv']
SIGMA = 2000.

# limits for \sigma_{basin} and \sigma_{north}
SB_LIMS = (-65, 40) 
SN_LIMS = (40, 65) 

# Set the decades
PERIOD = ['(years 01-05)', '(years 21-30)', '(years 141-150)']
YEARS  = ["1601-07-02", '1605-07-02',
          '1621-07-02', '1630-07-02',
          '1741-07-02', '1750-07-02'
         ]

LSTYLE = ['-','--', ':']

if SIGMA == 0.0:
   vmin, vmax, ci, cmap = 22,27.8,0.2, cmocean.cm.thermal_r
   levels = np.arange(vmin,vmax,ci)
   levels_c = [26.0,27.0,27.2,27.4,27.5,27.7]
   NAMEPLOT  = 'zonal_mean_Atlantic_pot_rho0'
else:
   vmin, vmax, ci, cmap = 33.0,37.4,0.2, cmocean.cm.thermal_r
   levels = np.arange(vmin,vmax,ci)
   levels_c = [35.0,35.6,36.0,36.4,36.6,36.8]
   NAMEPLOT  = 'zonal_mean_Atlantic_pot_rho2'

axissize  = 18
labelsize = 10

print('*****************************************')
print('Forcing =', FORCING)
print('The potential density is sigma', SIGMA)
print('Basin and North Box limits =', SB_LIMS, SN_LIMS)
print('*****************************************')

# Generate Basin Mask for the chosen sector
basin_mask = ms.ocean_basins.basin_mask_ht_full(OCEAN, check=False)

def sigma_zonal_mean(start, end, area, dzt, sigma, mask):
    """Compute time-mean volume-weighted zonal mean density."""
    volume = ms.derived.calc_volume(area, dzt)
    volume = volume.sel(time=slice(start,end)) * mask
    sigma  = sigma.sel(time=slice(start,end))  * mask
    sigma_z = (sigma * volume).cf.mean(dim='longitude') / volume.cf.mean('longitude')
    return sigma_z.mean('time') - 1000.0


# --- Plot the zonal-mean density and boxes profiles for each Forcing---
for i, forcing in enumerate(FORCING):
   print('Processing Forcing =', forcing)

   # Load the data sets
   data_path = BASE_PATH / f'FAFANTWATER_ideal_runoff_Boeira_rest_ice_{forcing}'
   area = xr.open_dataset(data_path / 'pp' /'areacello_t.nc')['area_t']
   dzt = xr.open_dataset(data_path / 'pp' / 'dht.nc')['dht']
   if SIGMA == 0.0:
      sigma = xr.open_dataset(data_path / 'pp' /'pot_rho_0.nc')['pot_rho_0']
   else:
      sigma = xr.open_dataset(data_path / 'pp' /'pot_rho_2.nc')['pot_rho_2']
   
   sigmaz_init = sigma_zonal_mean(YEARS[0], YEARS[1], area, dzt, sigma, basin_mask)
   sigmaz = sigma_zonal_mean(YEARS[4], YEARS[5], area, dzt, sigma, basin_mask)

   fig = plt.figure(i+1, figsize=(12,6))
   ax1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)
   ax2 = plt.subplot2grid((1, 5), (0, 4), colspan=2)
   sigmaz.plot.contourf(ax=ax1, 
              levels=levels, cmap=cmap, add_colorbar=False)
   sigmaz.plot.contour(ax=ax1, 
              colors='w', levels=levels_c,linewidths=0.5)
   sigmaz_init.plot.contour(ax=ax1,
              colors='w', levels=levels_c,linewidths=0.5, linestyles='dashed')
   ax1.axvline(SN_LIMS[0], linewidth=1, color='yellow', linestyle="dotted")
   ax1.axvline(SN_LIMS[1], linewidth=1, color='yellow', linestyle="dotted")
   ax1.axvline(SB_LIMS[0], linewidth=1, color='yellow', linestyle="dotted")
   ax1.axvline(SB_LIMS[1], linewidth=1, color='yellow', linestyle="dotted")
   #ax1.text(-15,5450,r'$\sigma_b$', size=14, color='y')
   #ax1.text(47,5450, r'$\sigma_n$', size=14, color='y')
   ax1.text(-17, 5000, r'$\sigma_b$',
         {'color': 'red', 'fontsize': 14, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="white", ec="black", pad=0.2)})
   ax1.text(52, 5000, r'$\sigma_n$',
         {'color': 'blue', 'fontsize': 14, 'ha': 'center', 'va': 'center',
          'bbox': dict(boxstyle="round", fc="white", ec="black", pad=0.2)})

   #ax1.set_title(f'{forcing}', loc='center', size=axissize)
   ax1.set_facecolor([0.5, 0.5,0.5])
   ax1.set_ylim([5500,0])
   ax1.set_xlim([-80,80])
   ax1.set_xticks(np.arange(-80,81,20))
   ax1.set_xticklabels(('80S','60S','40S','20S','0','20N','40N','60N','80N'))
   #ax1.set_xticks(np.arange(-60,80,20))
   ax1.set_ylabel('Depth [m]', fontsize=axissize)
   ax1.set_title(f'a) {forcing}', loc='left', size=axissize)
   ax1.set_xlabel('Latitude [$^{\circ}$]', fontsize=axissize)
   #ax1.set_xticklabels([-60,-40,-20,0,20,40,60])

   norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
   ax = plt.gcf().add_axes((.125,.020,.25,.01))
   cb1 = matplotlib.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=norm, boundaries=levels,
                                      orientation='horizontal', extend='neither',)
   if SIGMA==0.0:
      cb1.ax.set_title('$\sigma_0$ [Kg m$^{-3}$]', size=12)
   else:
      cb1.ax.set_title('$\sigma_2$ [Kg m$^{-3}$]', size=12)
   cb1.ax.tick_params(labelsize=8)


   # Area-mean vertical profile of sigma_n(z) and sigma_b(z)
   for j in range(len(PERIOD)):
       sigmaz = sigma_zonal_mean(YEARS[j+j*1], YEARS[j+j*1+1], area, dzt, sigma, basin_mask)
       if SIGMA==0.0:
          smin,smax = 23.5,28.0
       else:    
          smin,smax = 32.0,37.0
       rho_n = sigmaz.sel(yt_ocean=(slice(*SN_LIMS))).cf.mean('latitude') 
       rho_b = sigmaz.sel(yt_ocean=(slice(*SB_LIMS))).cf.mean('latitude')
       rho_n.plot(ax=ax2, y='st_ocean', yincrease=False, 
                       ylim=[3000,0], xlim=[smin,smax], color='blue', 
                       linewidth=1.0, linestyle=LSTYLE[j], label=PERIOD[j])
       rho_b.plot(ax=ax2, y='st_ocean', yincrease=False,
                       ylim=[3000,0], xlim=[smin,smax], color='red', 
                       linewidth=1.0, linestyle=LSTYLE[j])
       ##delta_rho = rho_n - rho_b
       ##delta_rho.plot.line(ax=ax2, y='st_ocean', yincrease=False, 
       ##                ylim=[5500,0], xlim=[-1.0,3.0], color=color[jj], 
       ##                linewidth=2.0, label=period[jj])
       ##ax2.axvline(delta_rho.cf.mean('vertical'), linewidth=1.0, 
       ##            color=color[jj], linestyle='--')
       ax2.set_xticks(np.arange(smin,smax+0.5,0.5))
       ax2.tick_params(labelsize=labelsize)
       ax2.tick_params("x", rotation=45)
       ax2.set_title(r"b) $\textcolor{red}{\sigma_b(z)}$ and $\textcolor{blue}{\sigma_n(z)}$", loc='left', size=16)
       #ax2.set_title(r"b) $\color{red}{\sigma_b}$ and $\color{blue}{\sigma_n}$", loc='left', size=axissize)
       ax2.legend(loc='lower left', ncol=1, fontsize=9, frameon=True)
       ax2.set_ylabel('Depth [m]', fontsize=axissize)
       ax2.tick_params(labelleft=False,labelright=True)
       ax2.yaxis.set_label_position('right')
       ax2.grid() 

   if PRINTING:
      plt.savefig(DIR_OUT / f'{NAMEPLOT}_{forcing}.png',
                  bbox_inches='tight',dpi=300)

plt.show()

