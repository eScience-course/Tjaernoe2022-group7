import xarray as xr
import intake
import cftime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import dask
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import cartopy as cy



ds_4xco2_10 = xr.open_mfdataset('../../Data/all_NorESM2-MM_abrupt-4xCO2_r1i1p1f1_000101-001012.nc', use_cftime=True)

def convert(ds, list_of_files, time):
    # converts a given dataset, ds, to have the same coordinates as the 10 first years from aburbt 4xco2 datatset 
    
    levels = ds_4xco2_10.lev.values
    lat = ds_4xco2_10.lat.values
    lon = ds_4xco2_10.lon.values


    ds_all = xr.Dataset({'time': time, 'lev': levels, 'lat': lat,'lon': lon})
    for path in list_of_files:

        ds = xr.open_dataset(path)
        for var in ds.data_vars:
                a = str(var)
        try:
                ds_all[a] = (['time', 'lat', 'lon'], ds[a].values)
        except:
                print(path)
                ds_all[a] = (['time', 'lev', 'lat', 'lon'], ds[a].values)

    ds_all = ds_all.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
    return ds_all


def fix_units(ds):
    # fixing units, and giving then new attributes
    
    ds['pr'] = ds['pr']*86400 # mm/year
    ds['pr'] = ds['pr'].assign_attrs(units='mm/year')
    ds['ts'] = ds['ts']-273.15 # C
    ds['ts'] = ds['ts'].assign_attrs(units='C$^\circ$')
    ds['ccn'] = ds['ccn']/100 # cm-3
    ds['ccn'] = ds['ccn'].assign_attrs(units='cm$^{-3}$')
    return ds



def annual(ds):
    # resample the data into yearly averages
    
    return ds.resample(time="Y").mean()


def weighted_mean(ds):
    # taking the weighted mean over lon and lat 
    
    weights = np.cos(np.deg2rad(ds.lat))
    weights.name = 'weights'
    air_weighted = ds.weighted(weights)
    weighted_mean = air_weighted.mean(('lon', 'lat'))
    return weighted_mean


def plot_map(ds_var, title, reverse=False, levels=6, extent=[-180,180,90,50], cmap='RdBu', vmin=-0.12, vmax=0.12):
    # Plotting a variables from a dataset on a map of the Arctic 
    
    if reverse:
            cmap=cmap+'_r'
    else:
            cmap=cmap
    
    sat_proj = ccrs.NorthPolarStereo()
    fig, ax = plt.subplots(figsize=(7,7),subplot_kw={'projection':sat_proj})
    ds_var.plot(facecolor="gray",
        ax = ax,
        cbar_kwargs={
            'orientation':'vertical',
            'shrink':.8
            },                 
        transform=ccrs.PlateCarree(),
        #levels=levels,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap
        )

    
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.gridlines(draw_labels=True)
    ax.coastlines()
    ax.set_title(title, fontweight='bold', fontsize=15)
    fig.tight_layout()




def fancy(ax, fontsize):
    # thickning the axes spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
        ax.spines[axis].set_color('k')
        
    # set the fontsize for all your ticks
    #fontsize = 15
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        
    # properties of the ticks
    ax.tick_params(direction='out', length=8, width=2, pad=10, bottom=True, top=False, left=True, right=False, color='k')
    
    # add a grid to the plot
    ax.grid(True, alpha=0.5)
    
    # mask top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    
    
