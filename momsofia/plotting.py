""" plotting.py - plotting utilities for momsofia """

import warnings

import numpy as np
import xarray as xr
import xgcm
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colorbar
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

__all__ = [
    "plot_coriolis",
    "plot_relvort",
    "plot_alpha_beta",
    "plot_n2",
    "plot_zonalmean_pv",
    "plot_zonalmean_PVG",
    "plot_amoc_scatter",
]

def plot_coriolis(coriolis):
    fig = plt.figure(10, figsize=(8, 5))
    ax1 = fig.add_subplot(1,1,1)
    f = coriolis * 1e4
    f.isel(xu_ocean=1).plot.line(ax=ax1)
    plt.title('Coriolis parameter\n'
              + r'$f = 2 \Omega sin(\phi)$ [10${-4}$ s$^{-1}$]')

def plot_relvort(zeta):
    fig = plt.figure(11, figsize=(10, 5))
    ax1 = fig.add_subplot(1,1,1)
    maxvalue = 1.
    levels = np.linspace(-maxvalue, maxvalue, 24)
    relvort = zeta * 1e6
    relvort.isel(st_ocean=2).cf.mean("time").plot.contourf(ax=ax1, levels=levels,
              cmap="RdBu_r", extend="both")
    plt.title('relative vorticity at surface\n'
        + r'$\zeta = \partial_x v - \partial_y u$ [10$^{-6}$ s$^{-1}$]')
    plt.xlim(-280, 80)
    plt.ylim(-90, 90);

def plot_alpha_beta(alpha, beta):
    ''' Plot alpha and beta at surface 
        and their zonal-mean time-averaged
    '''
    labelsize = 12
    titlesize = 10
    fig = plt.figure(12, figsize=(10,8))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    maxalpha = 4.0
    levels = np.arange(0.0, maxalpha, 0.2)
    thermal = alpha * 1e4
    thermal.isel(st_ocean=1).cf.mean("time").plot.contourf(ax=ax1, levels=levels,
            cmap=cm.cm.thermal)
    ax1.set_title('Thermal expansion at ocean surface [10$^{-4}$ C$^{-1}]$',
                 size=titlesize)

    thermal.mean('xt_ocean').cf.mean("time").plot.contourf(ax=ax3, levels=levels,
            cmap=cm.cm.thermal)
    ax3.set_title('Zonal mean Thermal expansion [10$^{-4}$ C$^{-1}]$', size=titlesize)
    ax3.invert_yaxis()
 
    maxbeta = 8.0
    levels = np.arange(7, maxbeta, 0.05)
    haline = beta * 1e4
    haline.isel(st_ocean=1).cf.mean("time").plot.contourf(ax=ax2, levels=levels,
             cmap=cm.cm.haline)
    ax2.set_title('Haline contraction at ocean surface [10$^{-4}$ PSU$^{-1}]$',
        size=titlesize)
    
    haline.mean('xt_ocean').cf.mean("time").plot.contourf(ax=ax4, levels=levels,
             cmap=cm.cm.haline)
    ax4.set_title('Zonal mean Haline contraction [10$^{-4}$ PSU$^{-1}]$',
                  size=titlesize)
    ax4.invert_yaxis()
    
    for ax in [ax1, ax2, ax3, ax4]:
        if ax==ax3 or ax==ax4:
           ax.set_ylabel(r'Depth [m]', fontsize=labelsize)
           ax.set_xlabel(r'Latitude [$^\circ$]', fontsize=labelsize)
        else:
           ax.set_ylabel(r'Latitude [$^\circ$]', fontsize=labelsize)
           ax.set_xlabel(r'Longitude [$^\circ$]', fontsize=labelsize)

    plt.tight_layout()

def plot_n2(n2, expt):
    labelsize = 14
    titlesize = 14
    cmap = cm.cm.deep
    fig = plt.figure(13, figsize=(8,5))
    ax1 = fig.add_subplot(1,1,1)
    maxn2 = 3.4
    levels = np.arange(0.0, maxn2, 0.20)
    n2 *= 1e4
    n2.plot.contourf(ax=ax1, levels=levels,
            cmap=cmap, add_colorbar=False)
    n2.plot.contour(ax=ax1, colors='k',
           levels=[0.2, 0.6, 1.0], linewidths=0.5)
    plt.gca().invert_yaxis()
    ax1.set_title(r'N$^{2}$', loc='center', size=titlesize)
    ax1.set_title('({expt})'.format(expt=expt), loc='right', size=titlesize);
    ax1.set_ylabel(r'Depth [m]', fontsize=labelsize)
    ax1.set_xlabel(r'Latitude [$^\circ$]', fontsize=labelsize)
    ax1.set_facecolor([0.5, 0.5,0.5])
    norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
    ax_cb = plt.gcf().add_axes((.92,.2,.02,.6))
    cb = matplotlib.colorbar.ColorbarBase(ax=ax_cb, cmap=cmap, norm=norm, boundaries=levels,
         extend='max', orientation='vertical', label=r'10$^{-4} [$s$^{-2}$]')

def plot_zonalmean_PV(pv, s2, s2ctl, basin, expt, years, fig):
    labelsize = 14
    titlesize = 14
    titlesize2 = 12
    ax = fig.add_subplot(1,1,1)
    
    cmap = cm.cm.deep
    if expt == '1.0Sv' or expt == '1.0Sv_60S':
       levels = (list(np.arange(0.0, 16.0, 4.0)) +
                 list(np.arange(16.0, 20.0, 2.0)) +
                 list(np.arange(20.0, 35.0, 5.0)))
    else:
       levels = (list(np.arange(0.0, 6.0, 3.0)) +
                 list(np.arange(6.0, 10.0, 1.0)) +
                 list(np.arange(10.0, 22.0, 5.0)))

    sigma_levels = [36.4,36.5,36.6,36.7,36.8,36.9]
 
    pv.mean('xt_ocean').mean('time').plot.contourf(ax=ax,
                        levels=levels, cmap=cmap, add_colorbar=False)
    #pv.mean('xt_ocean').mean('time').plot.contour(ax=ax,
    #                    colors='w', levels=levels, linewidths=0.5)
    cs=s2.mean('xt_ocean').mean('time').plot.contour(ax=ax,
          colors='purple', levels=sigma_levels, 
          linewidths=1.5, linestyles='-')
    ax.clabel(cs, inline=False, fmt='%1.1f', rightside_up=False, fontsize=8)
    if expt != '0.0Sv':
       print('--- Plotting also isopycnals for the Control')
       cs2=s2ctl.mean('xt_ocean').mean('time').plot.contour(ax=ax,
          colors='black', levels=sigma_levels, 
          linewidths=0.5, linestyles='-')
       #ax.clabel(cs2, inline=False, fmt='%1.1f', rightside_up=False, fontsize=8)
    plt.gca().invert_yaxis() 
    plt.ylabel(r'Depth [m]', fontsize=labelsize)
    plt.xlabel('Latitude [$^\circ$]', fontsize=labelsize)
    norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=cmap.N)
    ax_cb = plt.gcf().add_axes((.92,.2,.01,.6))
    cb = matplotlib.colorbar.ColorbarBase(ax=ax_cb, cmap=cmap, norm=norm, boundaries=levels,
                             orientation='vertical', label=r'PV [m$^{-1}$ s$^{-1}$] $\times$10$^{-11}$')
    ax.set_title('{expt}'.format(expt=expt), loc='center', size=titlesize); 
    ax.set_title('{years}'.format(years=years), loc='right', size=titlesize2); 
    ax.set_title('{basin}'.format(basin=basin), loc='left', size=titlesize); 
    

def plot_zonalmean_PVG(pvg, s2, s2ctl, basin, expt, years, fig):
    labelsize = 14
    titlesize = 14
    titlesize2 = 12
    ax = fig.add_subplot(1,1,1)
   
    cmap = cm.cm.curl
    levels = np.linspace(-4., 4., 21)
    sigma_levels = [36.4,36.5,36.6,36.7,36.8,36.9]
    pvg.mean('xt_ocean').mean('time').plot.contourf(ax=ax,
                        levels=levels, cmap=cmap, add_colorbar=False)
    cs=s2.mean('xt_ocean').mean('time').plot.contour(ax=ax,
          colors='purple', levels=sigma_levels, 
          linewidths=1.5, linestyles='-')
    ax.clabel(cs, inline=False, fmt='%1.1f', rightside_up=True, fontsize=8)
    if expt != '0.0Sv':
       cs2=s2ctl.mean('xt_ocean').mean('time').plot.contour(ax=ax,
          colors='black', levels=sigma_levels, 
          linewidths=0.5, linestyles='-')
       #ax.clabel(cs2, inline=False, fmt='%1.1f', rightside_up=True, fontsize=8)
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth [m]', fontsize=labelsize)
    plt.xlabel('Latitude [$^\circ$]', fontsize=labelsize)
    norm = matplotlib.colors.BoundaryNorm(boundaries=levels, ncolors=256, extend='both')
    ax_cb = plt.gcf().add_axes((.92,.2,.01,.6))
    cb = matplotlib.colorbar.ColorbarBase(ax=ax_cb, cmap=cmap, norm=norm, boundaries=levels,
         orientation='vertical', label=r'$\partial_y$PV [m$^{-2}$ s$^{-1}$] $\times$10$^{-11}$')
    ax.set_title('{expt}'.format(expt=expt), loc='center', size=titlesize)
    ax.set_title('{years}'.format(years=years), loc='right', size=titlesize2); 


def plot_amoc_scatter(amoc26N_c, amoc32S_c, amoc26N, amoc32S, exp, forcing, color):
    year_t1 = '1621' ; year_t2  = '1640'
    year_t3 = '1731' ; year_t4 = '1750'
    start_time1 = year_t1 + '-07-02'
    end_time1   = year_t2 + '-07-02'
    start_time2 = year_t3 + '-07-02'
    end_time2   = year_t4 + '-07-02'

    diff1 = amoc26N.sel(TIME=slice(start_time1, end_time1)).mean('TIME').squeeze() - amoc26N_c.sel(TIME=slice(start_time1, end_time1)).mean('TIME').squeeze()
    diff2 = amoc26N.sel(TIME=slice(start_time2, end_time2)).mean('TIME').squeeze() - amoc26N_c.sel(TIME=slice(start_time2, end_time2)).mean('TIME').squeeze()
    diff3 = amoc32S.sel(TIME=slice(start_time1, end_time1)).mean('TIME').squeeze() - amoc32S_c.sel(TIME=slice(start_time1, end_time1)).mean('TIME').squeeze()
    diff4 = amoc32S.sel(TIME=slice(start_time2, end_time2)).mean('TIME').squeeze() - amoc32S_c.sel(TIME=slice(start_time2, end_time2)).mean('TIME').squeeze()

    if exp == '1.0Sv_60S' or exp == '0.1Sv_60S':
       m = 'D'
       a = 1.0
       s = 30
    else:
       m = 'o'
       a = 1.0
       s = 50
    ax1.scatter(forcing, diff1, s=s, c=color, marker=m, alpha=a)
    ax2.scatter(forcing, diff2, s=s, c=color, marker=m, alpha=a)
    ax3.scatter(forcing, diff3, s=s, c=color, marker=m, alpha=a)
    ax4.scatter(forcing, diff4, s=s, c=color, marker=m, alpha=a)

    for ax in [ax1,ax3]:
        ax.set_ylim([0.0,2.0])
        ax.set_yticks(np.arange(0.0,2.25,0.25), fontsize=yaxissize)
    for ax in [ax2,ax4]:
        ax.set_ylim([-12,3.0])
        ax.set_yticks(np.arange(-12.,6.,2.), fontsize=yaxissize)
    for ax in [ax1]:
        ax.set_ylabel(r'$\Psi$ (26N) [Sv]', fontsize=titlesize)
    for ax in [ax3]:
        ax.set_ylabel(r'$\Psi$ (32S) [Sv]', fontsize=titlesize)
    for ax in [ax1,ax2,ax3,ax4]:
        ax.grid(color='grey', linestyle='solid', linewidth=0.2)
        ax.tick_params(labelsize=xaxissize)
        ax.set_xticks(np.arange(0.0,1.1,0.1), fontsize=xaxissize)
    for ax in [ax3,ax4]:
        ax.set_xlabel(r"Freshwater forcing [Sv]",fontsize=titlesize)
    for ax in [ax1]:
        ax.set_title('(years 21-40)', loc='center', size=titlesize)
        ax.set_title('a)', loc='left', size=titlesize)
    for ax in [ax2]:
        ax.set_title('b)', loc='left', size=titlesize)
        ax.set_title('(years 131-150)', loc='center', size=titlesize)
    for ax in [ax3]:
        ax.set_title('c)', loc='left', size=titlesize)
    for ax in [ax4]:
        ax.set_title('d)', loc='left', size=titlesize)

#fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),

