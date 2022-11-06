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

def trend(ds, variable, season):
    
    # create a dataset with same lon and lat as the precipitation dataset
    ds_trend = xr.Dataset({'lat': ds.lat,'lon': ds.lon})
    # add ds_trend_DJF to ds_NorESM2_MM_precip_season_trend
    ds_trend['mk_trend'] = (['lat', 'lon'], np.zeros((ds.lat.shape[0], ds.lon.shape[0]))*np.nan)
    ds_trend['mk_intercept'] = (['lat', 'lon'], np.zeros((ds.lat.shape[0], ds.lon.shape[0]))*np.nan)
    ds_trend['mk_p_val'] = (['lat', 'lon'], np.zeros((ds.lat.shape[0], ds.lon.shape[0]))*np.nan)  
    
    for ilat in range(ds[variable].shape[1]):
        for ilon in range(ds[variable].shape[2]):
            ts = ds[variable][:, ilat, ilon]
            results = mk.original_test(ts)
            ds_trend['mk_trend'][ilat, ilon] = results[7]
            ds_trend['mk_intercept'][ilat, ilon] = results[8]
            ds_trend['mk_p_val'][ilat, ilon] = results[2]
    return ds_trend

