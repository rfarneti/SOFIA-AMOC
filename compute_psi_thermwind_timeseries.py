'''
Compute the streamfunction from Thermal Wind balance
for the SOFIA sensitivity experiments
'''
import numpy as np
import xarray as xr
import momsofia as ms
from ocean_basins import basin_mask_ht_full
import matplotlib.pyplot as plt
from pathlib import Path
import my_style_latex

PRINTING = True

# --- Experiment Settings
FORCINGS = ['0.0Sv','0.1Sv','1.0Sv']#,'1.0Sv_60S']
DENSITY_TYPE = 'sigma0'  # Options: 'insitu', 'sigma0', 'sigma2'
OCEAN = 'Atlantic'

# Constants and parameters for TWB
C = 1.0 # constant of integration for the TWB
G = 9.8 # gravitational acceleration 
RHO0 = 1025.0 # reference density
F0 = 1.2e-4 # Coriolis parameter

# Limits for \sigma_{basin} and \sigma_{north} 
SB_LIMS = (-65, 40)
SN_LIMS = (40, 65)

# Constants for time conversion
DAYS_IN_YEAR = 365
OFFSET_YEAR = 1600

# Densities
DENS_NAME_MAP = {
    'insitu': 'rho',
    'sigma0': 'pot_rho_0',
    'sigma2': 'pot_rho_2'
}


def get_dataset_paths(forcing, density_type):
    """Maps forcing and density types to specific filenames."""
    moc_file = MOC_PATH / f'MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_{forcing}.nc'

    density_map = {
        'insitu': 'rho.nc',
        'sigma0': 'pot_rho_0.nc',
        'sigma2': 'pot_rho_2.nc'
    }
    dens_file = DATA_PATH / density_map[density_type]
    return moc_file, dens_file

def dens_zonal_mean(dens, area, dzt, mask):
    """Compute volume-weighted zonal mean."""
    volume = ms.derived.calc_volume(area, dzt) * mask
    dens_sliced = dens * mask
    return (dens_sliced * volume).cf.mean('longitude') / volume.cf.mean('longitude')

def compute_amoc_timeseries(moc_file, lat=26, z_range=(0, 4000)):
    """Computes the vertical maximum of the AMOC at a specific latitude over all time."""
    with xr.open_dataset(moc_file, decode_times=False) as ds:
        amoc_slice = ds['AMOC'].sel(YU_OCEAN=lat, method='nearest').sel(ST_OCEAN=slice(*z_range))
        amoc_index = amoc_slice.max(dim='ST_OCEAN')
        amoc_depth = amoc_slice.idxmax(dim='ST_OCEAN')

        time_years = amoc_index.TIME / DAYS_IN_YEAR - OFFSET_YEAR
        amoc_index = amoc_index.assign_coords(TIME=time_years)
        
        return amoc_index




BASE_PATH = Path('/Volumes/Antwater2/DATA/')
MOC_PATH  = BASE_PATH / 'MOC'
DIR_OUT   = Path('./Figures_transports/')
DIR_OUT.mkdir(exist_ok=True)


for forcing in FORCINGS:
    DATA_PATH = BASE_PATH / f'FAFANTWATER_ideal_runoff_Boeira_rest_ice_{forcing}' / 'pp'
    print(f"--- Computing AMOC Index time series for {forcing}")

    moc_path, dens_path = get_dataset_paths(forcing, DENSITY_TYPE)

    # --- Processing the AMOC from the model ---
    amoc_index_50 = compute_amoc_timeseries(moc_path, lat=50)
    amoc_index_30 = compute_amoc_timeseries(moc_path, lat=30)

    # --- Processing the TWB reconstruction ---
    print(f"--- Computing the TWB AMOC for {forcing}")
    target_dens = DENS_NAME_MAP[DENSITY_TYPE]
    dens_ds = xr.open_dataset(dens_path, decode_times=False)[target_dens]
    area = xr.open_dataset(DATA_PATH / 'areacello_t.nc')['area_t']
    dzt  = xr.open_dataset(DATA_PATH / 'dht.nc', decode_times=False)['dht']
    basin_mask = ms.ocean_basins.basin_mask_ht_full(OCEAN, check=False)

    dz = dens_zonal_mean(dens_ds, area, dzt, basin_mask)
    rho_n = dz.sel(yt_ocean=slice(*SN_LIMS)).cf.mean('latitude')
    rho_b = dz.sel(yt_ocean=slice(*SB_LIMS)).cf.mean('latitude')

    n_time = len(rho_b.time)
    n_depth = len(rho_b.st_ocean)
    z = -rho_b.st_ocean[::-1].values
    time_years = rho_b.time.values / DAYS_IN_YEAR - OFFSET_YEAR

    amoc_profiles_twb = np.zeros((n_depth, n_time))
    amoc_index_twb = np.zeros(n_time)

    print(f"--- Starting TWB reconstruction for {n_time} time steps")
    for tt in range(n_time):
        """Thermal-wind based solution for overturning
           Create column model instance and solve """
        rhon = rho_n.isel(time=tt).values[::-1]
        rhob = rho_b.isel(time=tt).values[::-1]
        twb = ms.psi_thermwind.Psi_Thermwind(z=z, f=F0, rho0=RHO0, bb=rhob, bn=rhon)
        twb.solve()

        #amoc_profiles_twb[:, tt] = twb.Psi
        #amoc_index_twb[tt] = np.max(twb.Psi[z < -4000])
        amoc_index_twb[tt] = np.max(twb.Psi)

    time_years = rho_b.time.values / DAYS_IN_YEAR - OFFSET_YEAR
    # Convert to xarray
    amoc_index_twb = xr.DataArray(
        amoc_index_twb, 
        coords={'time': time_years},
        dims=['time'],
        name='AMOC_index_TWB'
    )


    # --- Plotting
    axissize, legendsize = 16, 16
    xmin, xmax, dx = 0, n_time, 10
    ymin, ymax, dy = 4, 22, 2.0

    fig, ax = plt.subplots(figsize=(10, 6))
    amoc_index_50.plot(ax=ax, label=r'$\Psi_z$(50$^{\circ}$N)',  linewidth= 1.0, linestyle='-', color='grey')
    amoc_index_30.plot(ax=ax, label=r'$\Psi_z$(30$^{\circ}$N)',  linewidth= 1.0, linestyle='-', color='black')
    amoc_index_twb.plot(ax=ax, label=r'$\hat{\Psi}$',            linewidth= 3.0, linestyle='--', color='black')
    ax.set_ylim([ymin,ymax])
    ax.set_yticks(np.arange(ymin,ymax+dy,dy), fontsize=axissize)
    ax.set_xlim([0,n_time])
    ax.set_xticks(np.arange(xmin,xmax+dx,dx), fontsize=axissize)
    ax.tick_params(labelsize=axissize)
    ax.set_title(r"b)", loc='left', fontsize=18)
    ax.set_title(f"({forcing})", loc='right', fontsize=axissize)
    ax.set_title(r"$\Psi_z(t)$ and TWB reconstruction $\hat{\Psi}(t)$", loc='center', fontsize=axissize)
    ax.set_ylabel("[Sv]", fontsize=axissize)
    ax.set_xlabel("Time [Years]", fontsize=axissize)
    ax.legend(loc='lower left', fontsize=legendsize, frameon=True)
    ax.grid(True, alpha=0.3)

    if PRINTING:
        plt.savefig(DIR_OUT / f'AMOC_timeseries_twb_{forcing}.png', dpi=300, bbox_inches='tight')


#plt.tight_layout()
plt.show()

#amoc_ts.rolling(TIME=10, center=True).mean().plot(ax=ax, color='black', linewidth=2, label='10-year Rolling Mean')
