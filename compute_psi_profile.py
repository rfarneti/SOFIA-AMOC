'''
Compute the streamfunction 
at different latitudes as a function of depth
from the SOFIA sensitivity experiments

'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import my_style
# --- Configuration

FORCING_LEVELS = ['0.1Sv', '1.0Sv']
LATITUDES = [-32, 30, 50]
TIME_PERIODS = [
    {"label": "(years 01-05)", "start": "1601", "end": "1605", "ls": "-", "color": 'black'},
    {"label": "(years 36-40)", "start": "1636", "end": "1640", "ls": "-", "color": 'olive'},
    {"label": "(years 146-150)", "start": "1746", "end": "1750", "ls": "-", "color": 'purple'}
]

BASE_PATH = Path('/Volumes/Antwater2/DATA/MOC/')
DIR_OUT = Path('./Figures_transports/')
PRINTING = True

def get_amoc_profile(ds, start_yr, end_yr, lat):
    """Slices and averages the AMOC data for a specific time and latitude."""
    t_start, t_end = f"{start_yr}-07-02", f"{end_yr}-07-02"
    # Average over time and select nearest latitude
    return ds.sel(TIME=slice(t_start, t_end)).mean(dim='TIME').sel(YU_OCEAN=lat, method='nearest')



# --- Plotting ---
fig, axes = plt.subplots(len(FORCING_LEVELS), len(LATITUDES), 
                         figsize=(8,10), sharex=True, sharey=True, 
                         constrained_layout=False)

xmin, xmax, dx = -6.0, 18.0, 3.0

for r, forcing in enumerate(FORCING_LEVELS):
    file_path = BASE_PATH / f'MOC_FAFANTWATER_ideal_runoff_Boeira_rest_ice_{forcing}.nc'
    with xr.open_dataset(file_path) as ds:
        amoc_data = ds['AMOC']
        
        for c, lat in enumerate(LATITUDES):
            ax = axes[r, c]
            
            # Plot each time period
            for period in TIME_PERIODS:
                profile = get_amoc_profile(amoc_data, period['start'], period['end'], lat)
                profile.plot(ax=ax, y="ST_OCEAN", yincrease=False, 
                             color=period['color'], linestyle=period['ls'], 
                             label=period['label'])
            
            ax.set_title("")
            if r == 0:
               ax.set_xlabel("")
            else:
               ax.set_xlabel("[Sv]", fontsize=12)
            
            lat_label = f"{abs(lat)}$^\circ$N" if lat >= 0 else f"{abs(lat)}$^\circ$S"
            ax.set_title(lat_label, fontweight='bold', loc='right')

            if c == 0:
                ax.set_title(f"({forcing})", fontweight='bold', loc='left')
                ax.set_ylabel("Depth [m]", fontsize=12)
            else:
                ax.set_ylabel("")

            # Grid and limits
            ax.set_xlim(xmin, xmax)
            ax.set_xticks(np.arange(xmin, xmax+dx, dx))
            ax.set_ylim(5400, 0)
            ax.grid(True, alpha=0.3)

axes[0, 0].legend(fontsize=10, loc='lower right')

plt.tight_layout()
if PRINTING:
   DIR_OUT.mkdir(exist_ok=True)
   plt.savefig(DIR_OUT / f'psi_profiles_{forcing}_{lat}.png', bbox_inches='tight', dpi=300)

plt.show()
