'''
 Plots AMOC anomalies at different latitudes (30N, 40N, 50N and 32S)
 and at different decades
'''
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import xarray as xr
import cftime
import datetime
import netCDF4 as nc
from pathlib import Path
import my_style_latex

PRINTING = True
DIR_OUT   = Path('./Figures_transports/')
BASE_PATH  = Path('/Volumes/Antwater2/DATA/MOC/')

LATITUDES = [-32, 30, 50]

EXPERIMENTS = [
    ("0.01Sv",    0.01, "salmon"),
    ("0.02Sv",    0.02, "tomato"),
    ("0.05Sv",    0.05, "orangered"),
    ("0.07Sv",    0.07, "brown"),
    ("0.1Sv",     0.1,  "black"),
    ("0.1Sv_60S", 0.1,  "black"),
    ("0.1Sv_ambe", 0.1, "black"),
    ("0.15Sv",    0.15, "aquamarine"),
    ("0.2Sv",     0.2,  "cyan"),
    ("0.3Sv",     0.3,  "lime"),
    ("0.4Sv",     0.4,  "green"),
    ("0.5Sv",     0.5,  "royalblue"),
    ("0.6Sv",     0.6,  "blue"),
    ("0.7Sv",     0.7,  "blueviolet"),
    ("0.8Sv",     0.8,  "violet"),
    ("0.9Sv",     0.9,  "magenta"),
    ("1.0Sv",     1.0,  "purple"),
    ("1.0Sv_60S", 1.0,  "purple"),
    ("1.0Sv_ambe", 1.0, "purple")
]

START_YEARS = ['1621', '1741']
END_YEARS   = ['1630', '1750']
DECADES = ['(years 21-30)', '(years 141-150)']
LETTERS = ['a)', 'b)']

def compute_amoc_index(exp_id, lat, z_range=(0, 4000)):
    """Computes the vertical maximum of the AMOC at a specific latitude over all time."""
    """Loads a single experiment's AMOC index at a given latitude."""
    prefix = "MOC_FAFANTWATER_ideal_runoff_"
    suffix = "pv60_rest_ice_0.0Sv.nc" if exp_id == "0.0Sv" else f"Boeira_rest_ice_{exp_id}.nc"
    file_path = BASE_PATH / f"{prefix}{suffix}"

    with xr.open_dataset(file_path, decode_times = True) as ds:
         amoc_index = ds['AMOC'].sel(YU_OCEAN=lat, method='nearest').sel(ST_OCEAN=slice(*z_range)).max(dim="ST_OCEAN")
         return amoc_index.compute()

titlesize = 18
axissize  = 14
for i, lat in enumerate(LATITUDES):
    fig, axes = plt.subplots(1, 2, figsize=(10,5), sharex=True, constrained_layout=True)
    
    amoc_index_c = compute_amoc_index("0.0Sv", lat)

    for exp_id, forcing, color in EXPERIMENTS:
        amoc_index = compute_amoc_index(exp_id, lat)

        marker = 'D' if ("_60S" in exp_id or "_ambe" in exp_id) else 'o'
        if "_60S" in exp_id or "_ambe" in exp_id:
           size = 40
        else:
           size = 80
        if "_ambe" in exp_id:
           facecolor = 'none'
           edgecolor = 'black'
        else:
           facecolor = color
           edgecolor = 'white'

        for i, (s_yr, e_yr, decade, letter) in enumerate(zip(START_YEARS, END_YEARS, DECADES, LETTERS)):
            ax = axes[i]
            t_start, t_end = f"{s_yr}-07-02", f"{e_yr}-07-02" 
            ctl = amoc_index_c.sel(TIME=slice(t_start, t_end)).mean('TIME')   
            exp = amoc_index.sel(TIME=slice(t_start, t_end)).mean('TIME')   
            anom = exp - ctl 
            ax.scatter(forcing, anom, s=size, facecolors=facecolor, marker=marker, edgecolors=edgecolor, alpha=1.0)
            latt = int(lat)
            hem = 'N' if latt >= 0 else 'S' 
            ax.grid(color='grey', linestyle='solid', linewidth=0.2)
            ax.tick_params(axis='both', rotation=0, labelsize=axissize)
            #ax.set_xlim([0,1.1])
            #ax.set_xticklabels(np.arange(0.0,1.1,0.1))
            ax.set_xticks(np.arange(0.0,1.1,0.1))
            ax.set_title(f"{decade}", loc='right', size=titlesize)
            ax.set_title(f"{letter}", loc='left', size=titlesize)
            ax.set_xlabel(r"$\mathcal{F}_w$ [Sv]",fontsize=titlesize)

            if ax==axes[0]:
                ax.axvline(x=0.6, color="grey", linestyle=(0, (5, 5)), linewidth=0.5)
                ax.set_ylim([0.0,1.75])
                ax.set_yticks(np.arange(0.0,2.2,0.20))
                ax.set_ylabel(rf"$\Delta \Psi_z$ ({abs(latt)}$^\circ${hem}) [Sv]", fontsize=titlesize)
            if ax==axes[1]:
                ax.axvline(x=0.6,  color="grey", linestyle=(0, (5, 5)), linewidth=0.5)   
                ax.axvline(x=0.15, color="grey", linestyle=(0, (5, 5)), linewidth=0.5)   
                ax.set_ylim([-12,4.0])
                ax.set_yticks(np.arange(-12.,6.,2.))
                ax.set_ylabel('')

    plt.tight_layout()
    if PRINTING:
       plt.savefig(DIR_OUT / f'AMOC_vs_forcing_scatterplot_{abs(latt)}{hem}.png',
                   bbox_inches='tight',dpi=600)

plt.show()

