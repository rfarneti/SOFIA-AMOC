
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
import matplotlib.colors as mcolors
import scipy.stats as stats
import my_style_latex


DEPTH_PATH = Path('./DATA/')
DIR_OUT = Path('./Figures_transports/')
PRINTING = True


FORCING = ['1.0Sv','1.0Sv_60S']

VARIABLE = [('0.0',    'sigma0.0'), 
            ('2000.0', 'sigma2000.0')
]
label = r'$D \sim \Delta$b$^{-1/2}$'
#cmap = cmocean.cm.thermal_r
cmap = cmocean.cm.haline_r

g = 9.8
rho_0 = 1035

for forcing in FORCING:
    for density, sigma in VARIABLE: 
        ds_depth    = xr.open_dataset(DEPTH_PATH /f'pycnocline_depth_{sigma}_{forcing}.nc')
        ds_deltarho = xr.open_dataset(DEPTH_PATH /f'delta_rho_{sigma}_{forcing}.nc')
    
        Depth = ds_depth['Depth']
        time = ds_depth['time']
        Delta_rho = ds_deltarho['delta_rho']
        Delta_b = (g / rho_0) * Delta_rho
    
        #get a scaling: D \sim (\Delta b)**(-1/2) 
        C = Depth[0] * Delta_b[0]**(1/2)
        if sigma == 'sigma0.0':
           Delta_b_ext = np.arange(0.0060,0.010,0.0001)
        else:
           Delta_b_ext = np.arange(0.0065,0.009,0.0001)
        scaling = C.values * (1/Delta_b_ext)**(1/2)
        
 
        fig, ax = plt.subplots(figsize=(6,5))
        sc = ax.scatter(Delta_b * 10**2, Depth, c=time, marker='o', 
             cmap=cmap, s=40, edgecolors='0.1')
        ax.plot(Delta_b_ext * 10**2, scaling, color='black', linestyle='--', linewidth=1.0, label=label)

        bounds = np.linspace(0, 150.0, 16)
        norm = mcolors.BoundaryNorm(bounds, cmap.N)
        cbar = plt.colorbar(sc, ticks=bounds)
        cbar.set_label('Time [Years]', fontsize=14)

        if density == '2000.0':
           ax.set_ylim(800,825) if forcing == '0.1Sv' else ax.set_ylim(750,900)
           ax.set_yticks(np.arange(800,830,5)) if forcing == '0.1Sv' else ax.set_yticks(np.arange(750,925,25))
           ax.set_xlim(0.7,0.9) if forcing == '0.1Sv' else ax.set_xlim(0.6,1.1)
        elif density == '0.0':
           ax.set_ylim(850,950) if forcing == '0.1Sv' else ax.set_ylim(750,950)
           ax.set_yticks(np.arange(850,1000,25)) if forcing == '0.1Sv' else ax.set_yticks(np.arange(750,1000,25))
           ax.set_xlim(0.64,0.74) if forcing == '0.1Sv' else ax.set_xlim(0.4,1.0)
        
        plt.xlabel(r'$\Delta$b [10$^{-2}$ m s$^{-2}$]', fontsize=18)
        plt.ylabel(rf'D [m]', fontsize=18)
        plt.title(rf'b)', loc='left', fontsize=16)
        plt.title(rf'Pycnocline Depth vs $\Delta$b ({forcing})', loc='center', fontsize=16)
        #plt.title(f'({forcing})', loc='right', fontsize=14)
        plt.grid(alpha=0.3)
        plt.legend(loc='lower left', fontsize=14)
        ax.tick_params(labelsize=14)
 
        if PRINTING:
           plt.savefig(DIR_OUT / f'Delta_rho_vs_pycnocline_depth_{sigma}_{forcing}.png', 
                       bbox_inches='tight',dpi=300)

plt.show()


#Perform a linear regression on the log-transformed data
#amoc_arr = np.array(amoc)
#depth_arr = np.array(depth)

#log_psi = np.log(amoc_arr)
#log_d = np.log(depth_arr)
#slope, intercept, r_value, p_value, std_err = stats.linregress(log_psi, log_d)

#psi_range = np.linspace(min(amoc_arr), max(amoc_arr), 100)
#fit_line = np.exp(intercept) * (psi_range**slope)
#plt.plot(fit_line, psi_range, 'k--', label=f'Fit: $D \propto \Psi^{{{slope:.2f}}}$')

#print(f"Calculated scaling exponent: {slope:.2f}")
#print(f"Correlation (R^2): {r_value**2:.3f}")
