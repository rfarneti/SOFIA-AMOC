""" ocean_basins.py - module for generating oean basins

"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cf_xarray as cfxr

__all__=[
    "basin_tmask",
    "basin_mask_ht_full",
    "basin_mask_ht",
]


def basin_tmask(ocean):
    """ Returns mask for Atlantic basin or IndoPacific basin
        appropriate for zonal integrations
    """
    land_mask = xr.open_dataset('/Volumes/Antwater2/DATA/grids/regionmask_v6.nc')
    basin = land_mask.variables['tmask'][:]
    bmask = np.zeros(basin.shape)
    print('Generating a mask for the ', ocean)
    if ocean=='Atlantic':
       #bmask[(basin==2) | (basin==4) | (basin==6) | (basin==7) | (basin==8) | (basin==9) ] = 1
       bmask[(basin==2)] = 1
    elif ocean=='AtlanticArctic':
       #bmask[(basin==2) | (basin==4) | (basin==6) | (basin==7) | (basin==8) | (basin==9) ] = 1
       bmask[(basin==2) | (basin==4) ] = 1
    elif ocean=='IndoPacific':
       bmask[(basin==3) | (basin==5) ] = 1
    elif ocean=='Southern':
       bmask[(basin==1) ] = 1
    else:
       print("mask not available!")
       exit()
    vmask = 1*bmask; vmask[:-1,:] = vmask[:-1,:] * vmask[1:,:]
    vmask_1 = np.ma.array( vmask, mask=vmask==0)
    
    return vmask_1



def basin_mask_ht_full(basin):
    """Returns mask for Atlantic, Pacific or Indian Ocean sector,
       including the corresponding portion of Southern Ocean.
       Useful for zonal-means. Not suitable for zonal integrations!
       Modified from 
       github.com/COSIMA/cosima-recipes/blob/main/Recipes/Atlantic_IndoPacific_Basin_Overturning_Circulation.ipynb"""

    path = '/Volumes/Antwater2/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv/pp/'
    ht = xr.open_dataset(path+'ht.nc')['ht']
    land_mask = ~ht.isnull()

    #fig=plt.figure(1, (10, 5))
    #ax = plt.subplot()
    #land_mask.plot.contour(levels=[0.5], colors='k')
    #ax.plot([-280, 80], [-34,-34],      'r', linewidth = 1) #southern ocean
    #ax.plot([-69, -69], [-80, 9],       'b', linewidth = 3) #south america 
    #ax.plot([-83.7, -83.7], [9, 15.5],  'b', linewidth = 3) #panama
    #ax.plot([-93.3, -93.3], [15.5, 17], 'b', linewidth = 3) #panama2
    #ax.plot([-99, -99], [17, 90],       'b', linewidth = 3) #north america
    #ax.plot([22, 22], [-80,31],       'b', linewidth = 3)   #africa
    #ax.plot([79, 79], [31, 90],       'b', linewidth = 3)   #asia
    #ax.plot([-217, -217], [-80, -5],  'k', linewidth = 3)   #australia
    #ax.plot([-256, -217], [-5, -5],   'k', linewidth = 3)   #indonesia
    #ax.plot([-256, -256], [-5, 90],   'k', linewidth = 3)   #asia 
    #ax.set_xlim([-280, 80])
    #ax.set_ylim([-80, 90]);
    #ax.grid()
    
    print('Ocean basin = ',basin)

    ## create masks out of the above chunks
    if basin=='Pacific':
       pac_map1 = (land_mask.where(land_mask.yt_ocean < -5).where(land_mask.yt_ocean > -90).where(land_mask.xt_ocean > -217).where(land_mask.xt_ocean < -68)).fillna(0)
       pac_map2 = (land_mask.where(land_mask.yt_ocean > -5).where(land_mask.yt_ocean < 9).where(land_mask.xt_ocean > -256).where(land_mask.xt_ocean < -68)).fillna(0)
       pac_map3 = (land_mask.where(land_mask.yt_ocean > 9).where(land_mask.yt_ocean < 15.5).where(land_mask.xt_ocean > -256).where(land_mask.xt_ocean < -83.7)).fillna(0)
       pac_map4 = (land_mask.where(land_mask.yt_ocean > 15.5).where(land_mask.yt_ocean < 17).where(land_mask.xt_ocean > -256).where(land_mask.xt_ocean < -93.3)).fillna(0)
       pac_map5 = (land_mask.where(land_mask.yt_ocean > 17).where(land_mask.yt_ocean < 90).where(land_mask.xt_ocean > -256).where(land_mask.xt_ocean < -99)).fillna(0)
       pac_sector_map = pac_map1 + pac_map2 + pac_map3 + pac_map4 + pac_map5
       mask = pac_sector_map.where(pac_sector_map>0)
    elif basin=='Atlantic':
       atl_map1 = (land_mask.where(land_mask.yt_ocean < 9).where(land_mask.yt_ocean > -90).where(land_mask.xt_ocean > -69).where(land_mask.xt_ocean < 22)).fillna(0)
       atl_map2 = (land_mask.where(land_mask.yt_ocean > 9).where(land_mask.yt_ocean < 15.5).where(land_mask.xt_ocean > -83.7).where(land_mask.xt_ocean < 22)).fillna(0)
       atl_map3 = (land_mask.where(land_mask.yt_ocean > 15.5).where(land_mask.yt_ocean < 17).where(land_mask.xt_ocean > -93.3).where(land_mask.xt_ocean < 22)).fillna(0)
       atl_map4 = (land_mask.where(land_mask.yt_ocean > 17).where(land_mask.yt_ocean < 31).where(land_mask.xt_ocean > -99).where(land_mask.xt_ocean < 22)).fillna(0)
       atl_map5 = (land_mask.where(land_mask.yt_ocean > 31).where(land_mask.yt_ocean < 90).where(land_mask.xt_ocean > -99).where(land_mask.xt_ocean < 79)).fillna(0)
       atl_sector_map = atl_map1 + atl_map2 + atl_map3 + atl_map4 + atl_map5
       mask = atl_sector_map.where(atl_sector_map>0)
    elif basin=='Indian':
       ind_map1 = (land_mask.where(land_mask.yt_ocean < 31).where(land_mask.yt_ocean > -90).where(land_mask.xt_ocean > 22).where(land_mask.xt_ocean < 79)).fillna(0)
       ind_map2 = (land_mask.where(land_mask.yt_ocean > -90).where(land_mask.yt_ocean < -5).where(land_mask.xt_ocean > -280).where(land_mask.xt_ocean < -217)).fillna(0)
       ind_map3 = (land_mask.where(land_mask.yt_ocean > -5).where(land_mask.yt_ocean < 31).where(land_mask.xt_ocean > -280).where(land_mask.xt_ocean < -256)).fillna(0)
       ind_sector_map = ind_map1 + ind_map2 + ind_map3
       mask = ind_sector_map.where(ind_sector_map>0)
    else:
       print("mask not available!")
       exit()

    #fig=plt.figure(2, (10, 5))
    #ax = plt.subplot()
    #levels=np.arange(0,1.5,0.5)
    #mask.plot.contourf(ax=ax, levels=levels)
    #ax.set_xlim([-280, 80])
    #ax.set_ylim([-80, 90])
    #ax.set_title('Mask');

    return mask

def basin_mask_ht(basin, SO=False):
    """Returns mask for Atlantic or IndoPacific Ocean,
       Modified from 
       github.com/COSIMA/cosima-recipes/blob/main/Recipes/Atlantic_IndoPacific_Basin_Overturning_Circulation.ipynb"""

    path = '/Volumes/Antwater2/DATA/FAFANTWATER_ideal_runoff_Boeira_rest_ice_0.0Sv/pp/'
    ht = xr.open_dataset(path+'ht.nc')['ht']
    land_mask = ~ht.isnull()
   
    print('Ocean basin = ',basin)

    ## create masks out of the above chunks
    if basin=='IndoPacific':
       south_map = (land_mask.where(land_mask.yt_ocean < -34)).fillna(0)

       pac_map1 = (land_mask.where(land_mask.yt_ocean < -5).where(land_mask.yt_ocean > -34).where(land_mask.xt_ocean > -217).where(land_mask.xt_ocean < -68)).fillna(0)
       pac_map2 = (land_mask.where(land_mask.yt_ocean > -5).where(land_mask.yt_ocean < 9).where(land_mask.xt_ocean > -256).where(land_mask.xt_ocean < -68)).fillna(0)
       pac_map3 = (land_mask.where(land_mask.yt_ocean > 9).where(land_mask.yt_ocean < 15.5).where(land_mask.xt_ocean > -256).where(land_mask.xt_ocean < -83.7)).fillna(0)
       pac_map4 = (land_mask.where(land_mask.yt_ocean > 15.5).where(land_mask.yt_ocean < 17).where(land_mask.xt_ocean > -256).where(land_mask.xt_ocean < -93.3)).fillna(0)
       pac_map5 = (land_mask.where(land_mask.yt_ocean > 17).where(land_mask.yt_ocean < 90).where(land_mask.xt_ocean > -256).where(land_mask.xt_ocean < -99)).fillna(0)
       pac_map = pac_map1 + pac_map2 + pac_map3 + pac_map4 + pac_map5
       
       ind_map1 = (land_mask.where(land_mask.yt_ocean < 31).where(land_mask.yt_ocean > -34).where(land_mask.xt_ocean > 22).where(land_mask.xt_ocean < 79)).fillna(0)
       ind_map2 = (land_mask.where(land_mask.yt_ocean > -34).where(land_mask.yt_ocean < -5).where(land_mask.xt_ocean > -280).where(land_mask.xt_ocean < -217)).fillna(0)
       ind_map3 = (land_mask.where(land_mask.yt_ocean > -5).where(land_mask.yt_ocean < 31).where(land_mask.xt_ocean > -280).where(land_mask.xt_ocean < -256)).fillna(0)
       ind_map = ind_map1 + ind_map2 + ind_map3 
       if SO==True:
          indopac_map = pac_map + ind_map + south_map
       else:
          indopac_map = pac_map + ind_map
       mask = indopac_map.where(indopac_map>0)
    elif basin=='Atlantic':
       south_map = (land_mask.where(land_mask.yt_ocean < -34)).fillna(0)
       atl_map1 = (land_mask.where(land_mask.yt_ocean < 9).where(land_mask.yt_ocean > -34).where(land_mask.xt_ocean > -69).where(land_mask.xt_ocean < 22)).fillna(0)
       atl_map2 = (land_mask.where(land_mask.yt_ocean > 9).where(land_mask.yt_ocean < 15.5).where(land_mask.xt_ocean > -83.7).where(land_mask.xt_ocean < 22)).fillna(0)
       atl_map3 = (land_mask.where(land_mask.yt_ocean > 15.5).where(land_mask.yt_ocean < 17).where(land_mask.xt_ocean > -93.3).where(land_mask.xt_ocean < 22)).fillna(0)
       atl_map4 = (land_mask.where(land_mask.yt_ocean > 17).where(land_mask.yt_ocean < 31).where(land_mask.xt_ocean > -99).where(land_mask.xt_ocean < 22)).fillna(0)
       atl_map5 = (land_mask.where(land_mask.yt_ocean > 31).where(land_mask.yt_ocean < 90).where(land_mask.xt_ocean > -99).where(land_mask.xt_ocean < 79)).fillna(0)
       if SO==True:
          atl_sector_map = atl_map1 + atl_map2 + atl_map3 + atl_map4 + atl_map5 + south_map
       else:
          atl_sector_map = atl_map1 + atl_map2 + atl_map3 + atl_map4 + atl_map5 
       mask = atl_sector_map.where(atl_sector_map>0)
    else:
       print("mask not available!")
       exit()

    fig=plt.figure(11, (10, 5))
    ax = plt.subplot()
    levels=np.arange(0,1.5,0.5)
    mask.plot.contourf(ax=ax, levels=levels)
    ax.set_xlim([-280, 80])
    ax.set_ylim([-80, 90])
    ax.set_title('Mask');

    return mask

