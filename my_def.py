import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import netCDF4 as nc
#from netCDF4 import Dataset as nc
import cartopy.crs as ccrs
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.dates as mdates
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib.dates as mdates
import os
import cartopy
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import array as arr
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
from datetime import date
import math

params = {'xtick.labelsize': 'small',
          'ytick.labelsize': 'small',
          'text.usetex': False, 'font.size': 20,
          'font.family': 'serif', 'font.weight': 'normal'}
plt.rcParams.update(params)

# import data
def import_data_3D_data(file_name,var):
        file = nc.Dataset(file_name)
        x    = file.variables[var][:]
        lat  = file.variables['lat'][:]
        lon  = file.variables['lon'][:]
        return ma.mean(x,axis=0),lat,lon

def import_data_3D_data2(file_name,var):
        file = nc.Dataset(file_name)
        xx    = file.variables[var][:]
        lat  = file.variables['lat'][:]
        lon  = file.variables['lon'][:]
        f    = lambda x: ((x+180) % 360) - 180
        lon  = f(lon)
        ind  = np.argsort(lon)
        lon  = lon[ind]
        xx   = xx[:,:,ind]
        return ma.mean(xx,axis=0),lat,lon

def import_data_3D_month(file_name,var):
        file = nc.Dataset(file_name)
        xx   = file.variables[var][:]
        lat  = file.variables['lat'][:]
        lon  = file.variables['lon'][:]
        f    = lambda x: ((x+180) % 360) - 180
        lon  = f(lon)
        ind  = np.argsort(lon)
        lon  = lon[ind]
        xx   = xx[:,:,ind]
        return xx,lat,lon

def import_data_3D(file_name,var):
        file = nc.Dataset(file_name)
        x    = file.variables[var][:]
        lat  = file.variables['lat'][:]
        lon  = file.variables['lon'][:]
        return x,lat,lon

def import_data_3D_depth(file_name,var):
    file = nc.Dataset(file_name)
    xx   = file.variables[var][:]
    lat  = file.variables['lat'][:]
    lon  = file.variables['lon'][:]
    lev  = file.variables['lev'][:]
    f    = lambda x: ((x+180) % 360) - 180
    lon  = f(lon)
    ind  = np.argsort(lon)
    lon  = lon[ind]
    xx   = xx[:,:,:,ind]
    return xx,lev,lat,lon

def import_data_3D_depth2(file_name,var):
    file = nc.Dataset(file_name)
    xx   = file.variables[var][:]
    lat  = file.variables['lat'][:]
    lon  = file.variables['lon'][:]
    lev  = file.variables['depth'][:]
    f    = lambda x: ((x+180) % 360) - 180
    lon  = f(lon)
    ind  = np.argsort(lon)
    lon  = lon[ind]
    xx   = xx[:,:,:,ind]
    return xx,lev,lat,lon

def import_data_3D_depth3(file_name,var):
    file = nc.Dataset(file_name)
    xx   = file.variables[var][:]
    lat  = file.variables['lat'][:]
    lon  = file.variables['lon'][:]
    lev  = file.variables['olevel'][:]
    f    = lambda x: ((x+180) % 360) - 180
    lon  = f(lon)
    ind  = np.argsort(lon)
    lon  = lon[ind]
    xx   = xx[:,:,:,ind]
    return xx,lev,lat,lon

def mask_closest_to_coast_2D(var,n):
    mask = ma.ones(var.shape)
    for i in range(len(var[:,0])):
        z=ma.argmax(var[i,:]) #ma.nonzero(var[i,:])[0][-n:]
        mask[i,:][z] = 0
    return mask

def mask_closest_to_coast_3D(var,n):
        mask = ma.ones(var.shape)
        for i in range(len(var[:,0,0])):
                for j in range(len(var[0,:,0])):
                        z=ma.nonzero(var[i,j,:])[0][-n:]
                        #print(mask[i,j,:][z])
                        mask[i,j,:][z]=0
        return mask


def boxed_data(region,var,lat,lon):
        lat1 = region[0]
        lat2 = region[1]
        lon1 = region[2]
        lon2 = region[3]
        filt_lat = ((lat>=lat1)&(lat<=lat2))
        filt_lon = ((lon>=lon1)&(lon<=lon2))
        var = var[:,filt_lat,:]
        var = var[:,:,filt_lon]
        return var

def boxed_data2(region,var,lat,lon):
    lat1 = region[0]
    lat2 = region[1]
    lon1 = region[2]
    lon2 = region[3]
    f    = lambda x: ((x+180) % 360) - 180
    lon  = f(lon)
    ind  = np.argsort(lon)
    lon  = lon[ind]
    var  = var[:,:,ind]
    filt_lat = ((lat>=lat1)&(lat<=lat2))
    filt_lon = ((lon>=lon1)&(lon<=lon2))
    var = var[:,filt_lat,:]
    var = var[:,:,filt_lon]
    return var

def boxed_data_lon_depth(region,var,lev,lon):
    lev1 = region[0]
    lev2 = region[1]
    lon1 = region[2]
    lon2 = region[3]
    f    = lambda x: ((x+180) % 360) - 180
    lon  = f(lon)
    ind  = np.argsort(lon)
    lon  = lon[ind]
    var  = var[:,ind]
    filt_lev = ((lev>=lev1)&(lev<=lev2))
    filt_lon = ((lon>=lon1)&(lon<=lon2))
    var = var[filt_lev,:]
    var = var[:,filt_lon]
    return var

def boxed_data2d(region,var,lat,lon):
    lat1 = region[0]
    lat2 = region[1]
    lon1 = region[2]
    lon2 = region[3]
    f    = lambda x: ((x+180) % 360) - 180
    lon  = f(lon)
    ind  = np.argsort(lon)
    lon  = lon[ind]
    var  = var[:,ind]
    filt_lat = ((lat>=lat1)&(lat<=lat2))
    filt_lon = ((lon>=lon1)&(lon<=lon2))
    var = var[filt_lat,:]
    var = var[:,filt_lon]
    return var

def boxed_data_with_lat(region,var,lat,lon):
    lat1 = region[0]
    lat2 = region[1]
    lon1 = region[2]
    lon2 = region[3]
    filt_lat = ((lat>=lat1)&(lat<=lat2))
    filt_lon = ((lon>=lon1)&(lon<=lon2))
    var = var[:,filt_lat,:]
    var = var[:,:,filt_lon]
    l   = lat[filt_lat]
    return var,l

def boxed_data_with_lat_lon(region,var,lat,lon):
    lat1 = region[0]
    lat2 = region[1]
    lon1 = region[2]
    lon2 = region[3]
    filt_lat = ((lat>=lat1)&(lat<=lat2))
    filt_lon = ((lon>=lon1)&(lon<=lon2))
    var = var[:,filt_lat,:]
    var = var[:,:,filt_lon]
    lt   = lat[filt_lat]
    ln   = lon[filt_lon]
    return var,lt,ln

def boxed_data_with_lat_lon_depth(region,var,lat,lon,lev):
    lat1 = region[0]
    lat2 = region[1]
    lon1 = region[2]
    lon2 = region[3]
    filt_lat = ((lat>=lat1)&(lat<=lat2))
    filt_lon = ((lon>=lon1)&(lon<=lon2))
    filt_lev = (lev<=300)
    var = var[:,filt_lat,:]
    var = var[:,:,filt_lon]
    var = var[filt_lev,:,:]
    lt  = lat[filt_lat]
    ln  = lon[filt_lon]
    lv  = lev[filt_lev]
    return var,lv,lt,ln

def filtering_coordinate(region,lon):
    lat1 = region[0]
    lat2 = region[1]
    lon1 = region[2]
    lon2 = region[3]
    f    = lambda x: ((x+180) % 360) - 180
    lon  = f(lon)
    ind  = np.argsort(lon)
    lon  = lon[ind]

    if(ma.max(lon)<=90):
        filt_lat = ((lon>=lat1)&(lon<=lat2))
        lon  = lon[filt_lat]
        return lon
    else:
        filt_lon = ((lon>=lon1)&(lon<=lon2))
        lon = lon[filt_lon]
        return lon
    return lon



# plot with hatching
def plot_bias_hatched(lon,lat,var,z1,vmin,vmax,ci,cmap,fontsize):
    clevs=np.arange(vmin,vmax,ci)
    aspect = plt.figaspect(0.42)
    scale = '110m'
    fig=plt.figure(figsize=aspect)
    plt.subplots_adjust(top=1,
    bottom=0.02,
    left=0.081,
    right=0.974,
    hspace=0.21,
    wspace=0.225)
    ax=plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_global()
    land = cfeature.NaturalEarthFeature('physical', 'land', scale,edgecolor='face',
                                         facecolor=cfeature.COLORS['land'])
    ax.add_feature(land, facecolor='grey')
    ax.set_xticks([0, 60, 120, 240, 300, 360], crs=ccrs.PlateCarree())
    ax.set_yticks([-60, -30, 0, 30, 60],       crs=ccrs.PlateCarree())
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    lon_formatter = LongitudeFormatter(zero_direction_label=False)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())
    #ax.grid()

    f = lambda x: ((x+180) % 360) - 180
    lon = f(lon)
    ind = np.argsort(lon)
    lon = lon[ind]
    var = var[:, ind]
    z1 = z1[:, ind]

    mask = ma.ones(var.shape)
    cbar1 = ax.contourf(lon,lat,mask,cmap=cm.gray)
    cbar1 = ax.contourf(lon,lat,var,clevs,extend='both',            cmap=cmap,transform=ccrs.PlateCarree())
    cbar1 = ax.contourf(lon,lat,z1,1,hatches=['', '....'],  alpha=0,cmap=cmap,transform=ccrs.PlateCarree())
    cbar = ax.contour(lon,lat,var,clevs,colors='k',linewidths=0.5,linestyles='solid',transform=ccrs.PlateCarree())
    #ax.clim(-6,6)
    return fig,ax


# plot bias with a projection
def plot_bias_hatched2(lon,lat,var,z1,vmin,vmax,ci,cmap,fontsize):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1,1,1, projection=ccrs.EqualEarth(central_longitude=0))
    clevs=np.arange(vmin,vmax,ci)
    #norm = colors.BoundaryNorm(boundaries=clevs, ncolors=cmap.N)
    
    f = lambda x: ((x+180) % 360) - 180
    lon = f(lon)
    ind = np.argsort(lon)
    lon = lon[ind]
    var = var[:, ind]
    z1 = z1[:, ind]
   
    mask = ma.ones(var.shape)
    cs  = plt.contourf(lon,lat,mask,cmap=cm.gray)
    cs  = plt.contourf(lon,lat,var,clevs,extend='both',cmap=cmap,transform=ccrs.PlateCarree())
    cs1 = plt.contourf(lon,lat,z1, 1,hatches=['', '....'],  alpha=0,cmap=cmap,transform=ccrs.PlateCarree())
    cs2 = plt.contour(lon,lat,var,clevs,colors='k',linewidths=0.5,linestyles='solid',transform=ccrs.PlateCarree())
    ax.coastlines()
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='110m',edgecolor='face',
                                         facecolor=cfeature.COLORS['land'])
    #ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='black')
    ax.add_feature(land, facecolor='grey')
    #ax.set_extent([-180, 180, -60, 60], crs=ccrs.PlateCarree())
    #plt.title(title, fontsize=fonttitle)
    #cbar = plt.colorbar(cs, ticks=clevs,orientation='vertical', extend='both', pad=0.05)
    #cbar.ax.set_title(title_cb)
    #cbar.ax.tick_params(labelsize=8)    
    return cs,ax

# plot bias with a projection
def plot_bias_hatched_equator(lon,lat,var,z1,vmin,vmax,ci,cmap,fontsize):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_global()
    clevs=np.arange(vmin,vmax,ci)
    #norm = colors.BoundaryNorm(boundaries=clevs, ncolors=cmap.N)
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='110m',edgecolor='face',
                                         facecolor=cfeature.COLORS['land'])
    f = lambda x: ((x+180) % 360) - 180
    lon = f(lon)
    ind = np.argsort(lon)
    lon = lon[ind]
    var = var[:, ind]
    z1 = z1[:, ind]
    mask = ma.ones(var.shape)
    cs  = plt.contourf(lon,lat,mask,cmap=cm.gray)
    cs  = plt.contourf(lon,lat,var,clevs,extend='both',cmap=cmap,transform=ccrs.PlateCarree())
    cs1 = plt.contourf(lon,lat,z1, 1,hatches=['', '....'],  alpha=0,cmap=cmap,transform=ccrs.PlateCarree())
    cs2 = plt.contour(lon,lat,var,clevs,colors='k',linewidths=0.5,linestyles='solid',transform=ccrs.PlateCarree())
    #ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='black')
    #ax.add_feature(land, facecolor='grey')
    ax.set_xticks([0,10,20,30,40,50,60,90,120,150,180,210,240,270,290,300,310,320,330,340,350,360],
                   crs=ccrs.PlateCarree())
    ax.set_yticks([-60,-40,-30,-20,-10,0,10,20,30, 60],       crs=ccrs.PlateCarree())
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    lon_formatter = LongitudeFormatter(zero_direction_label=False)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_extent([-60, 20, -35, 20], crs=ccrs.PlateCarree())    
    return cs,ax


# plot diff in bias
def plot_bias(lon,lat,var,vmin,vmax,ci,cmap,fontsize):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1,1,1, projection=ccrs.Robinson(central_longitude=0))
    ax.set_global()
    ax.coastlines()
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    clevs=np.arange(vmin,vmax,ci)
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='110m',edgecolor='face',
                                         facecolor=cfeature.COLORS['land'])
    ax.add_feature(land, facecolor='grey')
    #norm = colors.BoundaryNorm(boundaries=clevs, ncolors=cmap.N)
    f = lambda x: ((x+180) % 360) - 180
    lon = f(lon)
    ind = np.argsort(lon)
    lon = lon[ind]
    var = var[:, ind]
    #mask = ma.ones(var.shape)
    #cs  = plt.contourf(lon,lat,mask,cmap=cm.gray)
    cs  = plt.contourf(lon,lat,var,clevs,extend='both',cmap=cmap,transform=ccrs.PlateCarree())
    cs2 = plt.contour(lon,lat,var,clevs,colors='k',linewidths=0.5,linestyles='solid',transform=ccrs.PlateCarree())
    #ax.set_xticks([0, 60, 120, 240, 300, 360], crs=ccrs.PlateCarree())
    #ax.set_yticks([-60, -30, 0, 30, 60],       crs=ccrs.PlateCarree())
    #plt.yticks(fontsize=fontsize)
    #plt.xticks(fontsize=fontsize)
    #lon_formatter = LongitudeFormatter(zero_direction_label=False)
    #lat_formatter = LatitudeFormatter()
    #ax.xaxis.set_major_formatter(lon_formatter)
    #ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_extent([-180, 180, -65, 60], crs=ccrs.PlateCarree())
    return cs,ax

def plot_bias_seta(lon,lat,var,vmin,vmax,ci,cmap,fontsize):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_global()
    #ax.coastlines()
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    clevs=np.arange(vmin,vmax,ci)
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='110m',edgecolor='face',
                                         facecolor=cfeature.COLORS['land'])
    ###ax.add_feature(land, facecolor='grey')
    #norm = colors.BoundaryNorm(boundaries=clevs, ncolors=cmap.N)
    f = lambda x: ((x+180) % 360) - 180
    lon = f(lon)
    ind = np.argsort(lon)
    lon = lon[ind]
    var = var[:, ind]
    mask = ma.ones(var.shape)
    cs1 = plt.contourf(lon,lat,mask,cmap=cm.gray)
    cs  = plt.contourf(lon,lat,var,clevs,extend='both',cmap=cmap,transform=ccrs.PlateCarree())
    cs2 = plt.contour(lon,lat,var,clevs,colors='k',
          linewidths=0.5,linestyles='solid',transform=ccrs.PlateCarree())
    ax.set_xticks([0,10,20,30,40,50,60,90,120,150,180,210,240,270,290,300,310,320,330,340,350,360],
                   crs=ccrs.PlateCarree())
    ax.set_yticks([-60,-40,-30,-20,-10,0,10,20,30, 60],       crs=ccrs.PlateCarree())
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    lon_formatter = LongitudeFormatter(zero_direction_label=False)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_extent([-60, 20, -35, 10], crs=ccrs.PlateCarree())
    return cs,ax

def plot_sst_seta(lon,lat,var,vmin,vmax,ci,cmap,fontsize):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_global()
    #ax.coastlines()
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    clevs=np.arange(vmin,vmax,ci)
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='110m',edgecolor='face',
                                         facecolor=cfeature.COLORS['land'])
    ##ax.add_feature(land, facecolor='grey')
    #norm = colors.BoundaryNorm(boundaries=clevs, ncolors=cmap.N)
    f = lambda x: ((x+180) % 360) - 180
    lon = f(lon)
    ind = np.argsort(lon)
    lon = lon[ind]
    var = var[:, ind]
    mask = ma.ones(var.shape)
    cs1 = plt.contourf(lon,lat,mask,cmap=cm.gray)
    cs  = plt.contourf(lon,lat,var,clevs,extend='both',cmap=cmap,transform=ccrs.PlateCarree())
    #cs2 = plt.contour(lon,lat,var,clevs,colors='k',
    #      linewidths=0.5,linestyles='solid',transform=ccrs.PlateCarree())
    ax.set_xticks([0,10,20,30,40,50,60,90,120,150,180,210,240,270,290,300,310,320,330,340,350,360],
                   crs=ccrs.PlateCarree())
    ax.set_yticks([-60,-40,-30,-20,-10,0,10,20,30, 60],       crs=ccrs.PlateCarree())
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    lon_formatter = LongitudeFormatter(zero_direction_label=False)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_extent([-50, 20, -35, 10], crs=ccrs.PlateCarree())
    return cs,ax


# color bar
def color_bar(levels,cmap):
    lower = np.min(levels)
    gap = abs(levels[1]-levels[0])
    upper = np.max(levels) + gap
    N = int((upper-lower)/gap)
    deltac = (upper-lower)/(2*(N-1))
    norm = colors.Normalize()
    cmap = plt.get_cmap(cmap, N+1)
    mapper = cm.ScalarMappable(norm=colors.Normalize(), cmap=cmap)
    mapper.set_array(np.linspace(lower-deltac,upper+deltac,N+1))
    #mapper.set_array([lower-deltac,upper+deltac])
    ticks = np.arange(lower,(upper+(gap/2)),gap)
    clb = plt.colorbar(mapper,orientation='horizontal',ticks=ticks,pad=0.10,shrink=0.70,aspect=50)
    return clb

def color_bar_v(cs, levels, cmap, **kwargs):
    lower = np.min(levels)
    gap = abs(levels[1]-levels[0])
    upper = np.max(levels) + gap
    N = int((upper-lower)/gap)
    deltac = (upper-lower)/(2*(N-1))
    norm = colors.Normalize()
    cmap = plt.get_cmap(cmap, N+1)
    mapper = cm.ScalarMappable(norm=colors.Normalize(), cmap=cmap)
    mapper.set_array(np.linspace(lower-deltac,upper+deltac,N+1))
    #mapper.set_array([lower-deltac,upper+deltac])
    ticks = np.arange(lower,(upper+(gap/2)),gap)
    #clb = plt.colorbar(mapper,orientation='vertical',ticks=ticks, cax=cbaxes, shrink=0.90)
    clb = plt.colorbar( mapper,orientation='vertical',ticks=ticks, pad=0.05, shrink=0.90, aspect=50)
    return clb


