import numpy as np
import xarray as xr
import pandas as pd


# Function for integrating aerosol concentrations
def pnsd_integration(ds, size_thresh, diameter_var = 'D'):
    '''
    - D is the mean diameter of the size bins
    - time is the time dimension
    - ds is the dataset containing the size distribution
    - pnsd is the variable in which the particle number size distribution is stored (dN/dlogDp)
    - size_thresh is the lower limit of the size bins to be integrated in nanometers
    '''
    Dp = ds.D.values
    logDp = np.log10(Dp)
    interval = np.array([logDp[i]-logDp[i-1] for i in range(1,np.size(Dp))])/2
    centers = logDp[:-1]+interval
    centers_bis = np.append(logDp[0]-interval[0], centers)
    centers_bis = np.append(centers_bis, logDp[-1]+interval[-1])
    bound_bin = 10**(centers_bis)#*10**(-9)
    dlogDp = np.array([np.log10(bound_bin[i+1])-np.log10(bound_bin[i]) for i in range(0, bound_bin.shape[0]-1)])
    pnsd_nolog = np.zeros(ds['pnsd'].values.shape)
    for i, Dp_i in enumerate(Dp):
        pnsd_nolog[:,i] = ds['pnsd'].sel(D=Dp_i).values * dlogDp[i]
    ds['pnsd_unlog'] = (['time', 'D'], pnsd_nolog)
    ds['N'+str(np.round(size_thresh, 2))] = ds['pnsd_unlog'].sel(D=slice(size_thresh, 1000)).sum(dim='D')


# Function for collocating the model data
def collocate(filelist):
    for ifile, file in enumerate(filelist):
        model = xr.open_mfdataset([s3.open('s3://'+file)])[['PRECC']]

        # change lon from 0-360 to -180-180 in model
        model['lon'] = ((model.lon + 180) % 360) - 180
        model = model.sortby('lon')

        print('file loaded: ' , ifile)

        for itime_mod, time_mod in enumerate(model.time):

            # find index of the time in the obs dataset
            idx = np.where(ds.time == (pd.to_datetime(time_mod.values).round('H') + timedelta(hours=1)))[0][0]

            # create a list of indexes to fill with corresponding values from the obs (diagonal)
            coord_list = [[idx+time_traj, time_traj] for time_traj in range(96)]

            for icoord, coord in enumerate(coord_list):
                # if list_lat[icoord] not nan and list_lon[icoord] not nan:
                lat = lat_array[coord[0], coord[1]]
                lon = lon_array[coord[0], coord[1]]
                if not np.isnan(lat) and not np.isnan(lon):
                    slice = model.sel(time=time_mod, lat=lat, lon=lon, method='nearest')
                    ds['Convective_precip'][coord[0], coord[1]] = slice['PRECC'].values
                    ds['Large_precip'][coord[0], coord[1]] = slice['PRECL'].values
                    ds['Total_CF'][coord[0], coord[1]] = slice['CLDTOT'].values
                else:
                    pass


# Function for normalizing the data
def normalize(df):
    df_norm = (df-df.min())/ (df.max() - df.min())
    return df_norm