'''
Compute the streamfunction from Thermal Wind balance
for the SOFIA sensitivity experiments
'''
import numpy as np
import xarray as xr
import cf_xarray as cfxr
import momsofia as ms
import matplotlib.pyplot as plt
import cmocean
from pathlib import Path
from ocean_basins import basin_mask_ht_full
import my_style
# --- Set basin, density (in-situ, sigma0, sigma2), strength of forcing, and Paths ---

# Experiment Settings
OCEAN = 'Atlantic'
FORCING = '0.0Sv'  # Options: '0.0Sv', '0.1Sv', '1.0Sv'
DENSITY_TYPE = 'sigma0'  # Options: 'insitu', 'sigma0', 'sigma2'
PRINTING = True

# Constants
C = 1.0 # constant of integration for the TWB
G = 9.8 # gravitational acceleration 
RHO0 = 1025.0 # reference density
F0 = 1.2e-4 # Coriolis parameter

# Limits for \sigma_{basin} and \sigma_{north} 
#SB_LIMS = (-60, 40) #also good?
SB_LIMS = (-65, 40)
SN_LIMS = (40, 65)

DECADES = ['(years 01-05)', '(years 36-40)', '(years 146-150)']
START_YEARS = ['1601', '1636', '1746']
END_YEARS   = ['1605', '1640', '1750']

DENS_NAME_MAP = {
    'insitu': 'rho',
    'sigma0': 'pot_rho_0',
    'sigma2': 'pot_rho_2'
}

BASE_PATH = Path('/Volumes/Antwater2/DATA/')
DATA_PATH = BASE_PATH / f'FAFANTWATER_ideal_runoff_Boeira_rest_ice_{FORCING}'
MOC_PATH = BASE_PATH / 'MOC'
DIR_OUT = Path('./Figures_transports/')
DIR_OUT.mkdir(exist_ok=True)

print('*****************************************')
print('Forcing =', FORCING)
print('Density for TWB reconstruction =', DENSITY_TYPE)
print('Basin and North Box limits =', SB_LIMS, SN_LIMS)
print('*****************************************')


def get_dataset_paths(forcing, density_type):
    """Maps forcing and density types to specific filenames."""
    moc_file = MOC_PATH / f'MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_{forcing}.nc'
    
    density_map = {
        'insitu': 'rho.nc',
        'sigma0': 'pot_rho_0.nc',
        'sigma2': 'pot_rho_2.nc'
    }
    dens_file = DATA_PATH / 'pp' / density_map[density_type]
    return moc_file, dens_file

def dens_zonal_mean(dens, start, end, area, dzt, mask):
    """Compute volume-weighted zonal mean."""
    volume = ms.derived.calc_volume(area, dzt).sel(time=slice(start, end)) * mask
    dens_sliced = dens.sel(time=slice(start, end)) * mask
    return (dens_sliced * volume).cf.mean('longitude') / volume.cf.mean('longitude')

def buoyancy(rho):
    return -g/rho0 * rho


# --- Data Loading ---
moc_path, dens_path = get_dataset_paths(FORCING, DENSITY_TYPE)

amoc_ds = xr.open_dataset(moc_path)['AMOC']
target_dens = DENS_NAME_MAP[DENSITY_TYPE]
dens_ds = xr.open_dataset(dens_path)[target_dens]
area = xr.open_dataset(DATA_PATH / 'pp' / 'areacello_t.nc')['area_t']
dzt  = xr.open_dataset(DATA_PATH / 'pp' / 'dht.nc')['dht']
basin_mask = ms.ocean_basins.basin_mask_ht_full(OCEAN, check=False)

# --- Main Processing & Plotting ---
fig, axes = plt.subplots(1, 3, figsize=(10, 6), sharey=True)
xmin, xmax, dx = -4.0, 20.0, 2.0

for i, (label, s_yr, e_yr) in enumerate(zip(DECADES, START_YEARS, END_YEARS)):
    ax = axes[i]
    t_start, t_end = f"{s_yr}-07-02", f"{e_yr}-07-02"

    # Process AMOC
    amoc_m = amoc_ds.sel(TIME=slice(t_start, t_end)).cf.mean('time')
    amoc50 = amoc_m.sel(YU_OCEAN=50, method='nearest')
    #amoc40 = amoc_m.sel(YU_OCEAN=40, method='nearest')
    amoc30 = amoc_m.sel(YU_OCEAN=30, method='nearest')

    # Process Density
    dz_mean = dens_zonal_mean(dens_ds, t_start, t_end, area, dzt, basin_mask).cf.mean('time')
    rho_n = dz_mean.sel(yt_ocean=slice(*SN_LIMS)).cf.mean('latitude')
    rho_b = dz_mean.sel(yt_ocean=slice(*SB_LIMS)).cf.mean('latitude')
    # Handle NaNs if necessary
    rho_n = np.where(np.isnan(rho_n), 1027.82, rho_n) #not really needed
    rhob = np.asarray(rho_b)
    rhon = np.asarray(rho_n)
    # Thermal-wind based solution for overturning
    rhob = rhob[::-1]
    rhon = rhon[::-1] 
    z = -rho_b.st_ocean[::-1].values
    # Create column model instance and solve
    twb = ms.psi_thermwind.Psi_Thermwind(z=z, f=F0, rho0=RHO0, bb=rhob, bn=rhon)
    twb.solve()

    # Convert TWB to Xarray
    da_twb = xr.DataArray(twb.Psi, coords={'depth': -twb.z}, dims=['depth'],
             attrs=dict(description="TWB streamfunction",units="Sv"),
             )

    # --- Plotting
    #ax.plot(twb.Psi * C, -twb.z, color='k', linewidth=3.0, label='$\hat{\Psi}$')
    da_twb.plot(ax=ax, y="depth",    color='black',linestyle='-', linewidth=2, label=r'$\hat{\Psi}$')
    amoc50.plot(ax=ax, y="ST_OCEAN", color='grey', linestyle='-', linewidth=1,  label=r'$\Psi$(50$^\circ$N)')
    amoc30.plot(ax=ax, y="ST_OCEAN", color='grey', linestyle='--', linewidth=1, label=r'$\Psi$(30$^\circ$N)')

    # Add Depth of Max
    for data, d_name, col in [(da_twb, 'depth', 'black'), (amoc50, 'ST_OCEAN', 'grey'), (amoc30, 'ST_OCEAN', 'grey')]:
        m_val, m_dep = data.max(), data.idxmax(dim=d_name)
        ax.plot([xmin, m_val], [m_dep, m_dep], '--o', color=col, markevery=[-1], markersize=3)

    ax.set_title("")
    ax.set_title(f"{FORCING}", loc='left', fontsize=10)
    ax.set_title(f"{label}",   loc='right', fontsize=10)

    if ax==axes[0]:
       ax.set_ylabel('Depth [m]', fontsize=12)
    else:
       ax.set_ylabel('')
    #   ax.set_yticklabels([])

for ax in axes:
    ax.grid(alpha=0.3)
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.arange(xmin,xmax+dx,dx))
    ax.set_ylim(5400, 0)
    ax.set_xlabel('[Sv]', fontsize=12)
    ax.tick_params(labelbottom=False, labeltop=True)
    ax.xaxis.set_label_position('top')

axes[-1].legend(loc='lower right', fontsize=9)
#plt.tight_layout()

if PRINTING:
    plt.savefig(DIR_OUT / f'psi_profile_twb_{FORCING}.png', dpi=300, bbox_inches='tight')

plt.show()
