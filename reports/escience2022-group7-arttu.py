import numpy as np
import pickle

# ##############################################################################
# Open txt files written by Appendix A-C codes
# ##############################################################################

def open_data(file_name):
    with open(file_name, 'rb') as handle:
        data = pickle.load(handle)
    return data


# ##############################################################################
# ERA5 precipitation calculation
# ##############################################################################

def calculate_total_precipitation_from_era5(sta_time_inter,sta_traj_inter,percipitation_data,zeppelin_data):
    zeppelin_vec = zeppelin_data["time"].values
    percipitation_inds = list()

    for idx, time_point in enumerate(sta_time_inter):
        #print(time_point)
        time_lim_min = time_point - np.timedelta64(30, "m")
        time_lim_max = time_point + np.timedelta64(30, "m")
        #print(time_lim_min)
        #print(time_lim_max)
    
        time_mask_lim_min = zeppelin_vec >= time_lim_min
        time_mask_lim_max = zeppelin_vec <= time_lim_max
    
        time_mask = time_mask_lim_min & time_mask_lim_max
    
        time_ind = np.where(time_mask == True)[0][0]
        #print(time_ind)
        percipitation_inds.append(time_ind)

    percipitation_data = list( percipitation_data[i] for i in percipitation_inds)

    # Calculating accumulated percipitation
    accumulated_percip_era5 = list()

    for i, trajectory_percipitation in enumerate(percipitation_data):
        if sta_traj_inter[i] == 0:
            traj_length = 0
        else:
            traj_length = np.arange(sta_traj_inter[i])
    
        percip_vec = np.sum(trajectory_percipitation[traj_length])
        #print(percip_vec)
        accumulated_percip_era5.append(percip_vec)

    accumulated_percip_era5 = np.array(accumulated_percip_era5)
    return accumulated_percip_era5


# ##############################################################################
# Zepplin station data precipitation calculation
# ##############################################################################

def calculate_accumulated_precipitation_zepplin_data(sta_time_inter,sta_traj_inter,zeppelin_data):
    accumulated_percip = list()

    for idx, trajectories in enumerate(sta_time_inter):
        #print(sta_traj_inter[idx])
        if sta_traj_inter[idx] == 0:
            traj_length = 0
        else:
            traj_length = np.arange(sta_traj_inter[idx])
    
        percip_vec = np.sum(zeppelin_data["Rainfall"].sel(time = trajectories,time_traj = traj_length).values)
        #print(percip_vec)
        accumulated_percip.append(percip_vec)

    accumulated_percip = np.array(accumulated_percip)
    return accumulated_percip


# ##############################################################################
# Distance calcultion functions
# ##############################################################################

def calc_dist_two_points(lat1,lon1,lat2,lon2):
    R = 6371*1e3
    p1 = lat1*np.pi/180
    p2 = lat2*np.pi/180
    dp = (lat2-lat1)*np.pi/180
    dL = (lon2-lon1)*np.pi/180
    a = np.sin(dp/2)*np.sin(dp/2) + np.cos(p1)*np.cos(p2)*np.sin(dL/2)*np.sin(dL/2)
    c = 2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
    return R*c
    
def calc_traj_distance(sta_time_inter,zeppelin_data):
    trajectory_length = list()

    for idx, trajectories in enumerate(sta_time_inter):
        trajectory_calculated_length = 0
        for j in range(95):
            lat1 = zeppelin_data["x"].sel(time=trajectories,time_traj=j).values
            lon1 = zeppelin_data["y"].sel(time=trajectories,time_traj=j).values
        
            lat2 = zeppelin_data["x"].sel(time=trajectories,time_traj=j+1).values
            lon2 = zeppelin_data["y"].sel(time=trajectories,time_traj=j+1).values

            distance_between_points = calc_dist_two_points(lat1,lon1,lat2,lon2)

            trajectory_calculated_length =+ distance_between_points
        
        trajectory_length.append(trajectory_calculated_length)

    trajectory_length = np.array(trajectory_length)
    return trajectory_length




