"""
Author: Thomas Sangrey
Date:   09/10/2023
Name: Observations_util.py
Summary:
    A set of observational procedures and functions for use by flake.py. These are used
    to curate and analyze observational and or meteorological data in cunjunction with
    FLake simulations run in fortran.
"""

import f90nml
import ctypes as ct
import numpy as np
from math import exp as exp
from math import pi as pi
from math import sin as sin
from matplotlib import pyplot as plt
from matplotlib import rc as rc
import datetime as date
from os import getcwd as cwd
from os import listdir
from os.path import isfile, join
import netCDF4
import pytz
import error_handling as error
from scipy.signal import lombscargle

T_w_freeze = 273.15    # Kelvin

def e_vp_from_DP(T_a, dp):
    """
    Input: 1) <float(T_a)> is air temperature at 2 meters (K)
    2) <float(dp)> is the dew point (K)
    Output: <float(e_a)> returns vapor pressure (Pa or Kg m-1 s-2)

    Converts dew point at given temp T_a to vapor pressure
    See  Bolton, David (1980). "The Computation of Equivalent Potential 
    Temperature" (PDF). Monthly Weather Review. 108 (7): 1046–1053. 
    Also see https://en.wikipedia.org/wiki/Vapour_pressure_of_water
    for August_Roche_Magnus Equation. Magnus equation returns kiloPascal
    so we need a prefix of 1000 to convert to Pa for use in FLake.
    """
    #=========== Empirical formula (Bolton, 1980) requires Celsius
    e_vp = 1000*0.6112*exp(17.67*(dp-T_w_freeze)/((dp-T_w_freeze)+243.5))
    return e_vp

def plt_rc(n, q, p):
    """
    Input: 1) <int(n)> (number of plots), 2) <int(q)> (number of rows),
    3) <int(p)> (number of columns)
    
    Output: (<int(row)>, <int(column)>) of nth plot in n-plot matrix

    Generates the plt row/column position (r, c) for nth plot
    in a qXp rowXcolumn plot matrix
    """
    r = int(n/p)   # modulo function
    c = n%p - 1    # div function
    if c == -1:
    	c= p-1
    	r=r-1
    return (r,c)

def plot_lk_progs_old(progs):
    """
    Input: <dict(progs)> is a dict containing prognastic variables (1-d arrays) 
    parameterized by time dimension
        i.e, progs = dict({"time": time, "T_mnw_out_arr": T_mnw_out_arr)
    Output: plot object- Makes plot matrix of some prognostic variables from
    a FLake calculation
    """
    time = progs.pop('time')
    progs_keys = [key for key in progs.keys()]
    progs_vals = [val for val in progs.values()]
    q = int(len(prog_keys)**(0.5))              # sets up the number of matrix rows and columns
    if len(prog_keys)**(0.5)-q > 10**(-6):      # adds one column if len(progs) is not perfect square
        p=q+1
    else:
        q=p
    fig, ax = plt.subplots(q,p)
    n = 1
    for prog_key, prog_val in zip(progs_keys, progs_vals):
        prog_title = prog_key + "vs time"
        ax[plt_rc(n,q,p)].set_title(prog_title)
        label = f'prog_key'
        ax[plt_rc(n,q,p)].plot(time, prog_val,label = label)
        ax[plt_rc(n,q,p)].set_xlabel(f'time (days)')
        ax[plt_rc(n,q,p)].set_ylabel(f'{prog_key}')
        ax[plt_rc(n,q,p)].legend(bbox_to_anchor=(1, 1.0), loc='upper left', \
            prop={'size': 10})
        n += 1
    #plt.tight_layout()
    return plt


def get_nml_meteo(nml):
    """
    Input:  <dict> ported from an nml file. Contains lake features of individual lakes
            from a limnology text file including local meteorological data vs time. 
    Output: 6 <float(arrays)> of 1-D containing lake data along time axis.
            
            1) time (sec); 2) T_a_in_arr: air temp(K); 3) I_atm_in_arr: short-wave 
            flux to ground (W/m^2); 4) Q_a_in_arr: vapor pressure (Pa); 5) U_a_in_arr: 
            windspeed (m/s); 6) CL_tot_arr: fractional cloud cover (dimensionless)

    Unpacks 5 variables from a single meteorological text file (defined within nml file 
    for individual lake and has extension .csv or .dat). .csv has one line header. .dat 
    has no header with or without header, but with columns: time_str, I_atm_in_arr, 
    T_a_in_arr, Q_a_in_arr,U_a_in_arr, CL_tot_arr for time (d), short wave solar flux 
    (W/m^2), air temp (K), vapor pressure (Pa), wind speed (m/s), and fractional
    cloud cover (dimensionless, 0-1).
    """
    meteo_file = nml['meteo']['meteofile']     # grabs 
    ext = meteo_file.split(".")[-1]
    #=========== determine whether the file is .dat or .csv ===============
    #=========== .csv file has a single header row ========================
    if(ext=="dat"):
        return np.loadtxt(meteo_file, unpack=True, usecols = [0,1,2,3,4,5])
    elif(ext=="csv"):
    #=========== Convert time strings to datetime objects and then to seconds. ===
    #=========== Then recompute time array relative to initial time ==============
        time_str, I_atm_in_arr, T_a_in_arr, Q_a_in_arr,U_a_in_arr, CL_tot_arr =\
        np.loadtxt(nml['meteo']['meteofile'], unpack=True, usecols = [0,1,2,3,4,5], \
        	delimiter=",", dtype=str, skiprows = 1)          # load the meteo data into arrays
        
        time_0_str = time_str[0]
        dt = date.datetime.strptime(time_0_str, '%Y-%m-%d %H:%M')
        time_0 = dt.timestamp()
        time = [0 for i in range(len(time_str))]             # initiaze/allocate time array
        i = 0
        for t in time_str:
            t = date.datetime.strptime(t, '%Y-%m-%d %H:%M')  #time string is converted to days
            time[i] = float((t.timestamp() - time_0)/86400)    # Convert time to days
            i += 1
    else:
    	pass
    I_atm_in_arr = [float(I) for I in I_atm_in_arr]
    T_a_in_arr   = [float(t) for t in T_a_in_arr]
    Q_a_in_arr   = [float(q) for q in Q_a_in_arr]
    U_a_in_arr   = [float(u) for u in U_a_in_arr]
    CL_tot_arr   = [float(cl) for cl in CL_tot_arr]
    return time, I_atm_in_arr, T_a_in_arr, Q_a_in_arr,U_a_in_arr, CL_tot_arr

def get_daily_Lake_temp(daily_lake_file_fp, pos, r, plot_area = False):
    """
    Input: 
    1) <str(rel_path)>: relative pathway (without leading '/') of directory containing Lake temperature 
    files. Each file contains daily global lake surface temp. Files are of type netCDF4
    (.nc). 
    2) pos = <[float(lat), float(lon)]> (degrees)- the center postion of a geographic region 
    (small) to plot or average daily temperature from a single file (day). 
    3) float(r) (degrees): The +/- boundaries of region described in 2), 
    i.e., max_lat = pos[0]+r, min_lon - pos[1] - r
    4) <bool(plot_area)> is flag to determine whether to make plots or no

    Output: None  (produces plots)

    Retrieves satellite lake_surface_water_temperature data from all files in ./rel_path at 
    pos=[lat +/- r,long +/- r] in lat and long. 

    Each input file is a day of global satellite data and is in netCDF format.
    Flake/LAKES/USGS/<files> used by TS are presently from CEDA (WGS84 projection) for years
    1995, 1996, 1997, 1998, 1999, 2000. 
    
    Carrea, L.; Embury, O.; Merchant, C.J. (2015): GloboLakes: high-resolution 
    global limnology dataset v1. Centre for Environmental Data Analysis
    lake_surface_water_temperature. 

    Problem: Lots of daily data gaps. Likely due to clouds.
    """
   
    print(f'file for get_daily_lake_temp() is : {daily_lake_file_fp}\n')
    ds = netCDF4.Dataset(file)   # form the dataset
    lat_low = pos[0]-r           # form the box boundaries for framing a lake
    lat_hi  = pos[0]+r           # r box length and width is 2*r
    lon_low = pos[1]-r
    lon_hi  = pos[1]+r
    lats = np.array(ds["lat"][:])   #All lat and lon
    lons = np.array(ds["lon"][:])
    lat_bool = np.logical_and(lats > lat_low, lats < lat_hi)  #mask the data according to box
    lon_bool = np.logical_and(lons > lon_low, lons < lon_hi)
    lon_bool_lst = [i for i, b in enumerate(lon_bool) if b==True] #form list[indices] of data within box
    lat_bool_lst = [i for i, b in enumerate(lat_bool) if b==True]
    #??????------How do I use booleans or logical operators directly here? ----------
    T_sfc_arr = np.array(ds["lake_surface_water_temperature"][0,\
    lat_bool_lst,lon_bool_lst])         # LSWT from within box
    T_sfc_arr[T_sfc_arr<-32000]=300     # set filled data to 300K for good plot contrast
    if plot_area:
        plt.imshow(T_sfc_arr)           # heatmap of LSWT
        plt.show()
    else:
        pass

def plt_rc(n, q, p):
    """
    Generates the plt (row, column) position (r, c) for nth plot
    in a q X p row-major plot matrix
    """

    if n > q*p:
        raise Exception('The plot number, n, exceeds available plots.')
    else:
        r = int(n/p)
        c = n%p - 1
        if c == -1:
            c = p-1
            r = r-1
        if p == 1:
            return r
        else:
            return (r, c)

def plot_lake_matrix(data):
    time = data.pop('time')          # pop time and title to leave variables
    title = data.pop('title')        # Not used yet
    data_keys = [key for key in data.keys()]
    data_vals = [val for val in data.values()]
    q = 3
    p = 3
    fig, ax = plt.subplots(q,p)
    n = 1
    for data_key, data_val in zip(data_keys, data_vals):
        data_title = data_key + "vs time"
        ax[plt_rc(n,q,p)].set_title(data_title)
        label = f'data_key'
        ax[plt_rc(n,q,p)].plot(time, data_val,label = label)
        ax[plt_rc(n,q,p)].set_xlabel(f'time (days)')
        ax[plt_rc(n,q,p)].set_ylabel(f'{data_key}')
        ax[plt_rc(n,q,p)].legend(bbox_to_anchor=(1, 1.0), loc='upper left', prop={'size': 10})
        n += 1
    #plt.tight_layout()
    return plt

def Quality_mask_Lake_temp(ds_dict, t_idx, QL_min, unc_max, Q_dict):

    err = False
    
    QL_mask     = ds_dict["QL_mask_all"][t_idx,:,:]
    unc_mask    = ds_dict["unc_mask_all"][t_idx,:,:]
    lswt_mask   = ds_dict["lswt_mask_all"][t_idx,:,:]

    common_mask = np.logical_not(np.logical_and(np.logical_and(QL_mask, lswt_mask), \
                  unc_mask))
    
    QL_full_box      = ds_dict["QL_all"][t_idx,:,:]
    unc_lswt_full_box   = ds_dict["unc_all"][t_idx,:,:]
    lswt_full_box        = ds_dict["lswt_all"][t_idx,:,:]
    """
    ===============================================================================
    !!!!! Possible BUG fix here. (see "except" case below).Make a common mask for all
    three variables They should always be the same mask,but for some reason, 
    quality_level is not masked for some masked values of lswt and lswt_uncertainty.
    This happens for several pixels in some files.

    Old Code:
    QL          = ds.variables["quality_level"][t_idx,:,:].data[QL_mask] 
    unc_lswt    = ds.variables["lswt_uncertainty"][t_idx,:,:].data[unc_mask]
    lswt        = ds.variables["lake_surface_water_temperature"][t_idx,:,:].data[lswt_mask]
    ===============================================================================
    """  
    QL          = QL_full_box[common_mask] 
    unc_lswt    = unc_lswt_full_box[common_mask]
    lswt        = lswt_full_box[common_mask]
      
    Q_mask      = np.logical_and(QL >= QL_min, unc_lswt < unc_max)
    #=========== Change sign of fillvalue for masked uncertainty values. This way, 
    #=========== an upper threshold (unc_max) will unmask a simply connected region
    #=========== we can call a "region of interest" 
    unc_lswt_full_box[np.logical_not(common_mask)] = +32768
    Q_mask_full_box  = np.logical_and(QL_full_box >= QL_min, unc_lswt_full_box < unc_max)
    try: 
        Q_lswt  = lswt[Q_mask]
    except:
        """
        =============================================================================
        !!!!BUG case here. Unclear how the number of lswt unfilled values can be less
        than QL or unc_lswt unfilled values. I would think these would always be the same.
        But here, I have to catch the case when there are too few lswt values compared
        to QL or unc_lswt values which, ostensibly, depend on unfilled lswt values
        =============================================================================
        """
        err = error.WARN_lswt_too_small(QL_mask, unc_mask, lswt_mask, QL, unc_lswt, \
                     lswt, Q_mask, t_idx)
        Q_lswt    = np.array([])
                
    if (Q_lswt.size > 0) and not err:
        mn_lswt   = Q_lswt.mean()
        num_lswt  = Q_lswt.size
    else:
        Q_dict["not_warned"] = error.WARN_Q_mask_no_pixels(Q_dict["not_warned"])
        Q_lswt    = np.array([])
        mn_lswt   = float('nan')
        num_lswt  = 0.0

    Q_dict["lswt_full_box"]     = lswt_full_box
    Q_dict["Q_mask_full_box"]   = Q_mask_full_box
    Q_dict["lswt_mask"]         = lswt_mask
    Q_dict["Q_lswt"]            = Q_lswt
    Q_dict["Q_mask"]            = Q_mask
    Q_dict["mn_lswt"]           = mn_lswt
    Q_dict["num_lswt"]          = num_lswt
    
    return Q_dict


def plot_lake_vars(ds_dict, t_idx, lswt_props, T_sfc_bool, time_array, k):
    
    lats = lswt_props["lswt_all"]["lats"]
    lons = lswt_props["lswt_all"]["lons"]
    q = 3; p = 1;

    fig, ax = plt.subplots(q,p)
    font1 = {'weight' : 'medium',
            'size'   : 9}
    font2 = {'weight' : 'bold',
            'size'   : 12}
    
    fourcorners = [lons[0], lons[-1], lats[0], lats[-1]]
    plt_kwargs = {'interpolation': 'none', 'extent': fourcorners, 'aspect': 'auto'}
    for n, key in enumerate(lswt_props.keys()):
        data = ds_dict[key][t_idx,:,:]
        data[T_sfc_bool] = lswt_props[key]["fill_val"]
        data = data/lswt_props[key]["norm_factor"]
        ax[plt_rc(n,q,p)].set_title(lswt_props[key]["abbrev"], **font2)
        vmin = lswt_props[key]["cbar_range"]["vmin"]
        vmax = lswt_props[key]["cbar_range"]["vmax"]
        an_im = ax[plt_rc(n,q,p)].imshow(data, **plt_kwargs, \
                    vmin = vmin, vmax = vmax, cmap = plt.get_cmap('jet'))
        cb = fig.colorbar(an_im,ax = ax[plt_rc(n,q,p)])
        cb.set_label(lswt_props[key]["z_label"], **font1)
        ax[plt_rc(n,q,p)].set_xlabel(f'lon', **font1)
        ax[plt_rc(n,q,p)].set_ylabel(f'lat', **font1)
        plt.xticks(**font1)
        plt.yticks(**font1)
        plt.subplots_adjust(wspace=None, hspace=0.5)
        x = (fourcorners[0] + 2*fourcorners[1])/4
        y = (fourcorners[2] + 2*fourcorners[3])/4
        date_str = f'{date.datetime.fromtimestamp(time_array[k]*86400)}'
        ax[plt_rc(n,q,p)].text(x, y, date_str, fontdict=font1, color = "white")
        plt.show(block=False)
    return (fig, ax)    

def get_ROI(Q_dict, lswt_props):
    """
    Input: Q_dict = <dict([[float or bool,],])> 
    Q_dict = {Q_dict["lswt_full_box"]     = lswt_full_box      ,   
              Q_dict["Q_mask_full_box"]   = Q_mask_full_box    ,
              Q_dict["lswt_mask"]         = lswt_mask          ,     
              Q_dict["Q_lswt"]            = Q_lswt             ,
              Q_dict["Q_mask"]            = Q_mask             ,
              Q_dict["mn_lswt"]           = mn_lswt            ,
              Q_dict["num_lswt"]          = num_lswt}          ,
    Output: ROI = <plt>
    Use Q_dict values to make a region of interest (ROI) plot.  
    """
    font1 = {'weight' : 'medium',
            'size'   : 9}
    font2 = {'weight' : 'bold',
            'size'   : 12}
    lats = lswt_props["lswt_all"]["lats"]
    lons = lswt_props["lswt_all"]["lons"]
    Q_dict["lswt_full_box"][Q_dict["lswt_mask"]]                       = 0
    Q_dict["lswt_full_box"][np.logical_not(Q_dict["lswt_mask"])]       = 1
    Q_dict["lswt_full_box"][Q_dict["Q_mask_full_box"]] = 2
    fig, ax = plt.subplots(1,1)
    ax.set_title("Region of Interest", **font2)
    ROI = plt.imshow(Q_dict["lswt_full_box"])
    ax.set_xlabel(f'lon', **font1)
    ax.set_ylabel(f'lat', **font1)
    plt.xticks(**font1)
    plt.yticks(**font1)
    plt.show(block = False)
    return ROI 

def get_Lake_temp_vs_time(lake_file_fp, cavu_thresh, QL_min, unc_max, \
                          startday=None, endday=None, plot_area = False):
    """
    Input: 
    1) <str(relative pathway)>: Without leading '/' where a single lake_file lives. 
    2) <str(lake_file)>: A single lake temperature file containing daily lake surface
       temp for several decades (1995-06-01 - 2016-12-31) Files must be of type netCDF4 
       (.nc). Time is given in seconds since Jan 1st, 1981. Python gives time referenced
       from posix standard time of January 1st, 1970 at 00-00-00 UTC.
    3) <float> cavu_thresh is a visability cutoff: Only keep T_obs for images with at 
       least cavu_thresh fraction of visible pixels since some satellite imagery is sparse
    4) <str(startday)> give as a string readable by Datetime. Format "YYYY-mm-dd-HH-MM-SS".
    5) <str(endday)> give as a string readable by Datetime. Format "YYYY-mm-dd-HH-MM-SS".
    6) <bool(plot_area)> determines whether to plot the heat map of the temperature at
    a given instant.
    
    Output: to be added
    """
    #================================ Initialization ==================================
    print(f'file is {lake_file_fp}\n')
    ds = netCDF4.Dataset(lake_file_fp)   # form the dataset
    lats       = np.array(ds.variables["lat"][:])   #All lats in image box
    lons       = np.array(ds.variables["lon"][:])   #All lons in image box
    
    #=========== Initialize a set of properties to help with plotting
    lswt_props = {};
    lats       = np.array(ds.variables["lat"][:])   #All lats in image box
    lons       = np.array(ds.variables["lon"][:])   #All lons in image box
    lswt_props.update({'lswt_all': {'abbrev': 'LWST', \
            'fill_val': 273, 'norm_factor': 1, 'z_label': 'Temperature',   \
            "lats": lats, "lons": lons, "cbar_range": {"vmin": None, "vmax": None}}})
    lswt_props.update({'unc_all': {'abbrev': 'LWST_unc', \
             'fill_val': 0, 'norm_factor': 0.01, 'z_label': "%% uncertainty",\
             "cbar_range": {"vmin": 0, "vmax": 100}}})
    lswt_props.update({'QL_all': {'abbrev': 'QL', \
             'fill_val': 0, 'norm_factor': 1, 'z_label': 'Quality', \
             "cbar_range": {"vmin": 0, "vmax": 5}}}) 
    #=============================End Initialization ==================================
    
    #=========== Create an offset that reflects the difference in measured
    #=========== seconds of dataset (referenced to 1/1/1981)from the unix posix standard
    #=========== time (referenced to 1/1/1970 00:00:00<UTC>).
    timeref_dt = date.datetime.strptime("1981-01-01-00-00-00", '%Y-%m-%d-%H-%M-%S')
    timeref_dt = pytz.utc.localize(timeref_dt)
    offset     = timeref_dt.timestamp()
    time_array = ds.variables["time"][:].data
    #=========== Use start_sec and end_sec to set points for cropping the dataset 
    #=========== Make sure each are referenced to the posix standard time
    if startday is not None:
        start_dt   = date.datetime.strptime(startday, '%Y-%m-%d-%H-%M-%S') # referenced to 1/1/1970 00:00:00, local
        start_dt   = pytz.utc.localize(start_dt)
        start_sec  = start_dt.timestamp()               # referenced to 01-01-1970-00-00-00<UTC>
    else:
        start_sec = time_array[0] + offset
    if endday is not None:
        end_dt     = date.datetime.strptime(endday, '%Y-%m-%d-%H-%M-%S')  # referenced to 1/1/1970 00:00:00, local
        end_dt     = pytz.utc.localize(end_dt)
        end_sec    = end_dt.timestamp()    # referenced to 01-01-1970-00-00-00<UTC>
    else:
        end_sec = time_array[-1] + offset
    print(f'start_sec = {start_sec} and end_sec is {end_sec}\n')
    
    #=========== Reference time array to posix standard (UTC)
    time_array   = np.array([(t + offset) for t in time_array]) 
    t_idx_array  = np.array([t_idx for t_idx, t in enumerate(time_array)])
    
    #=========== Then crop time_array according to startsec and endsec
    crop_bool    = np.logical_and(time_array >= start_sec, time_array <= end_sec) 
    time_array   = time_array[crop_bool]
    t_idx_array  = t_idx_array[crop_bool]
    time_array   = np.array([t/86400 for t in time_array])
    
    #=========== cavu is number of filled pixels, later normalized to give
    #=========== estimate of image visibility cavu[i]= 1-(cloud_fraction[i]).
    #=========== In addition to quality_level and lswt_uncertainty, cavu is
    #=========== used to score and discard low visibility images. 
    cavu = np.zeros(time_array.size)           # Initialize visability
    T_obs= np.zeros(time_array.size)           # Initialize T_obs 
    Q_dict = {"not_warned": True}                      
    ds_dict={}
    
    QL_mask_all     = ds.variables["quality_level"][:,:,:].mask
    unc_mask_all    = ds.variables["lswt_uncertainty"][:,:,:].mask
    lswt_mask_all   = ds.variables["lake_surface_water_temperature"][:,:,:].mask
    QL_all          = ds.variables["quality_level"][:,:,:].data
    unc_all         = ds.variables["lswt_uncertainty"][:,:,:].data
    lswt_all        = ds.variables["lake_surface_water_temperature"][:,:,:].data    
    
    ds_dict.update({"lswt_all": lswt_all, "QL_all": QL_all, "unc_all": unc_all})
    ds_dict.update({"lswt_mask_all": lswt_mask_all, "QL_mask_all": QL_mask_all, \
                    "unc_mask_all": unc_mask_all})

    for k, t_idx in enumerate(t_idx_array):
        Q_dict     = Quality_mask_Lake_temp(ds_dict, t_idx, QL_min, unc_max, Q_dict)
        cavu[k]    = Q_dict["num_lswt"]
        T_obs[k]   = Q_dict["mn_lswt"]
        
        #=========== Plot the image (first d images)        
        
        if plot_area and 0<=k<7:
            T_sfc_bool   = Q_dict["lswt_mask"]
            lswt_plt     = plot_lake_vars(ds_dict, t_idx, lswt_props, T_sfc_bool, \
                           time_array, k)                                 
        else:
            pass
    
    #=========== Demonsrate the region of interestUse the maximum atmospheric clarity (cavu) which presents the most
    #=========== the most pixels of a lake image.
    cavu_max_idx = cavu.argmax()
    Q_dict     = Quality_mask_Lake_temp(ds_dict, cavu_max_idx, QL_min, unc_max, Q_dict)
    get_ROI(Q_dict, lswt_props)
    T_sfc_bool   = Q_dict["lswt_mask"]
    plot_lake_vars(ds_dict, cavu_max_idx, lswt_props, T_sfc_bool, \
                           time_array, k)
    #=========== cavu needs to be normalized to a fractional value of its largest member
    cavu_max = cavu.max()
    cavu = cavu/cavu_max     # normalizes the cavu array
    
    #=========== Discard (as 'nan') any T_obs based on visibility (cavu) less than
    #=========== cavu_thresh
    T_obs = [T if c > cavu_thresh else float('nan') for T, c in zip(T_obs, cavu)]
    #=========== Plot the region of interest that is defined by quality thresholds

    
    return time_array, T_obs


def PeriodoGram_lake_time(time, T_obs):
    #===========================================
    #=========== This is a sandbox area for now....
    #===========================================
    #=========== Make sure both inputs are numpy arrays
    T_obs = np.array(T_obs)
    time = np.array(time)
    nan_mask  = np.logical_not(np.isnan(T_obs))
    T_obs_orig = T_obs[nan_mask]
    time_orig = time[nan_mask]
    #============ Gets the pth batch of q equal sections of the full time series
    #============ and performs the periodogram. Might be useful to quantify
    #============ noise or otherwise anomolous signals related to quality factors
    fig, ax = plt.subplots(4,5)
    p=3; L= 20; r = 10
    for q in list(range(p+1,L)):
        time = time_orig[p*int(int(time_orig.size/q)):(p+1)*int(time_orig.size/q)]
        T_obs = T_obs_orig[p*int(T_obs_orig.size/q):(p+1)*int(T_obs_orig.size/q)]
        points = time.size*r
        lo_f = 1/(2*np.pi)/time.ptp()
        hi_f = points/(2*np.pi)/time.ptp()
        freq = np.linspace(lo_f, hi_f, points)
        periodogram = lombscargle(time, T_obs, freq)
        ax[plt_rc(q,4,5)].plot(freq, np.sqrt(4*periodogram/(points)/(2*np.pi)**2))
        ax[plt_rc(q,4,5)].set_xlabel('Frequency (1/day)')
    plt.show(block=False)


def select_lake(Tsfc_vs_time_DIR_fp, lake_lim_file_FP, **selection):
    """
    Input: 
    1) Tsfc_vs_time_rel_path: 
    2)lake_lim_file: <current directory>/<lakefile> where lake_file is relative path and filename
    of lake information. Currently: 
    "/LAKES/GLOBO/LAKES_limnology/globolakes-static_lake_centre_fv1.csv"
    3) keyword argument or dict of form: {"IDs": [int(ID)]},
    {"names": [str(lakename)]} or {"bounds_arrs": [ {"lat_lo": real(lat_lo),
    "lat_hi": real(lat_hi), "lon_lo": real(lon_lo), "lon_hi": real(lon_hi)}} ]
    Output: str(fullpath) gives the file name of the selected lake from 
    directory LAKES_vs_time.
    
    Forms a crude database (DB) of basic lake features from the limnology file and is 
    adjoined by the string values of T_sfc_vs_time files for one or more lakes from the
    directory LAKES/LAKES_vs_time. The DB can be built up 1) from an array of global lake IDs,
    2) from a series of lake names, or 3)from a box domain of lat/long that contains the lakes.

    Returns dict with lake_names as keys and "lake_file": <file>, "country": <country>, 
    {"centroid": {"lat": lat}, {"lon": lon}}, and 
    {"lake_box": {"minlat": minlat, ,"maxlat": ,maxlat, "maxlat": maxlat, ,"maxlon": ,maxlon} 
    as dict components of the output dict.
    """
    
    files = [f for f in listdir(Tsfc_vs_time_DIR_fp) if \
            isfile(join(Tsfc_vs_time_DIR_fp, f)) and f.split(".")[-1] == "nc"]
    print(f'The number of files is {len(files)}\n')
    #============  Pull lake info from Global Lakes and Wetland Database ====
    #============  (GLWD) database summary ==================================
    GLWDs, lakes, countries, lats, lons, minlats, maxlats, minlons, maxlons = \
    np.loadtxt(lake_lim_file_FP, unpack=True, usecols = [0,1,2,3,4,5,6,7,8], \
            delimiter=",", dtype=str, skiprows = 42,encoding="latin-1")
    GLWDs_lkup = {}                 # A lookup table to match GLWD to GLWD index
    for i in range(len(GLWDs)):
        GLWDs_lkup[str(GLWDs[i])] = i 
    lakes_dict = {}
    try:                       # grab GLWD IDs from **kwargs, if they exist.
        IDs = selection["IDs"] # Error exception if they don't.
        for f in files:
            for id_num, ID in enumerate(IDs):
                #====== extract GLWD from beginning of lake_file name
                GLWD = int((f.split("-")[0]).split("LAKE")[-1])
                if GLWD == ID:             # a lake file matches a requested ID         
                    indx = GLWDs_lkup[str(GLWD)]  # use lookup to map GLWD 
                                                  # to index of limnology variables
                    lakes_dict[lakes[indx]] = {"lake_file": f, \
                            "country": countries[indx],     # ======== Store values
                            "lake_centroid": {"lat": lats[indx], "lon": lons[indx]},\
                            "lake_box": {"minlat":minlats[indx],"maxlat":maxlats[indx],\
                            "minlon":maxlons[indx], "maxlon": maxlons[indx]}}
                    break    # break from ID[:] loop because ID<=>GLWD is one to one
                else:
                    pass
        return lakes_dict
    except Exception as error:   # if no **kwargs for IDs exist, try names
        print(f'Warning in \"IDs\" selection:\n {error}')
        try:
            names = selection["names"]   # try names
            #  lowercase and no whitespace for text comparison
            names_tmp = [name.replace(' ','').lower() for name in names]
            lakes_tmp = [lake.replace(' ','').lower() for lake in lakes]
            for name in names_tmp:
                for L, lake_name in enumerate(lakes_tmp):
                    if lake_name == name:              # lake_name matches name **kwars request
                        for f in files:                # now find matching lake_file
                            #====== extract GLWD from beginning of lake_file name
                            GLWD =int((f.split("-")[0]).split("LAKE")[-1])  
                            if int(GLWDs[L]) == GLWD:     # found lake match to file
                                lakes_dict[lakes[L]] = {"lake_file": f, \
                                "country": countries[L],    # Store values
                                "lake_centroid": {"lat": lats[L], "lon": lons[L]},\
                                "lake_box": {"minlat":minlats[L],"maxlat":maxlats[L],\
                                "minlon":minlons[L], "maxlon": maxlons[L]}}
                                break       # break from files[:] loop because file<=>GLWD is 
                                            # one to one mapping
                            else:
                                pass 
                    else:
                        pass
            return lakes_dict
        except Exception as error:
        #========= If no **kwargs for IDs ornames, try finding all lakes ========
        #========= inside a defined region ======================================
            print(f'Warning in \"names\" selection:\n {error}')
            try:
                bounds_arrs = selection["bounds_arrs"]    # defines a box in which to search for lakes
                for b, bounds in enumerate(bounds_arrs):
                    for f in files:
                        #====== extract GLWD from beginning of lake_file name
                        GLWD =int((f.split("-")[0]).split("LAKE")[-1])
                        indx = GLWDs_lkup[str(GLWD)]     # use lookup to map GLWD 
                                                          # to index of limnology variables
                        #======= define a few booleans to flag whether a lake's extent =====
                        #======= falls with boxed region ===================================
                        within_lats = bounds["lat_lo"] <= float(minlats[indx]) \
                                     and bounds["lat_hi"] >= float(maxlats[indx])
                        within_lons = bounds["lon_lo"] <= float(minlons[indx]) \
                                     and bounds["lon_hi"] >= float(maxlons[indx])
                        if within_lats and within_lons:
                            lakes_dict[lakes[indx]] = {"lake_file": f, \
                                "country": countries[indx],                 # append output struct
                                "lake_centroid": {"lat": lats[indx], "lon": lons[indx]},\
                                "lake_box": {"minlat":minlats[indx],"maxlat":maxlats[indx],\
                                "minlon":minlons[indx], "maxlon": maxlons[indx]}}
                                # no break statement here since we actually need to test 
                                # each file in directory
                        else:
                            pass
                return lakes_dict
            except Exception as error:       # Arriving here means no lakes were found
                print(f'Error in lake selection.:\n {error}')
                return None

def get_ERA5_met(path, pos):
    """
    Input: 
    1) <str(path)> location and filename of netCDF4 file containing ERA5
    data. File has variables:('longitude', 'latitude', 'expver', 'time', 'u10', 
    'v10', 'd2m', 't2m','lcc', 'msdwlwrf', 'msdwswrf', 'sp', 'tcc', 'tp')
    2) <[float(pos)]> 2 element array containing [lat, lon] for a single location
    Output:
    1) <dict(ERA_dict)> is a dict that contains information usefule to 
    the ERA5 dataset.
    2) (<[float]>,) 14 element tuple of 1-d arrays (lists) dimensioned by 
    time and containing meteorological variables ordered chronologically.
    Gets the meteorological variable information within the nearest 0.025 X 0.025 
    degree^2 region of a given [lat, lon] position (i.e. lake position). 

    See ERA5 information for this dataset: 
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

    netCDF.nc file has dimensions (N-time, M-expvers, P-lat, P-lon) where M=2 
    and refers to experimental version. ERA5 and ERA5T (near real-time) are two
    versions of reanalysis. With a 3 month delay, ERA5T covers the last three
    months of data before dataset's final date (also called "present date"). For
    example, if the final complete day of dataset is August 7th, the dataset will
    will have 68 days of ERA5T data and the rest will be ERA5. That means 30 days
    of June, 31 days of July and 7 days of August. The instantaneous (i.e. T_a, 
    wind)and accumulated variables (i.e. solar flux, precipitation) have slightly 
    different implementations of ERA5 and ERA5T (7 hours of padding for accumulated
    variables). See ERA5 documentation: 
    https://confluence.ecmwf.int/pages/viewpage.action?pageId=173385064

    In this particular implementation, I will drop the most recent three months
    to avoid the hassle of mixed expver ERA5/ERA5T data.
    """
        
    def smooth(an_array, run):
        """
        simple runnning average to reduce noise
        """
        smooth_arr = [sum(an_array[i:i+run])/run for i in range(len(an_array)-run)]
        return smooth_arr

    #=========== Initialize a dictionary to carry adjunct information useful for
    #=========== for the ERA5 dataset
    ERA5_dict = {}
    #=========== Get the meteo file as netCDF4 object
    met = netCDF4.Dataset(path)
    lats = met.variables["latitude"][:].data
    lons = met.variables["longitude"][:].data
    #=========== Find the "four corners" of lats/lons closest to target lat/lon
    #=========== Use the four corners to do a linear 2-d interpolation of meteo
    #=========== variables, if desired. Presently, for simplicity, I only grab 
    #=========== the closest corner and use its meteo information exclusively.
    lats_idx_lo = (np.where(lats < pos[0]))[0][0]    # lats is in descending order
    lats_idx_hi = (np.where(lats > pos[0]))[0][-1]   # lats is in descending order
    lons_idx_lo = (np.where(lons < pos[1]))[0][-1]   # lons is in ascending order
    lons_idx_hi = (np.where(lons > pos[1]))[0][0]    # lons is in ascending order
    
    if abs(pos[0]-lats[lats_idx_lo]) < abs(pos[0]-lats[lats_idx_hi]):
        lat_idx = lats_idx_lo
    else:
        lat_idx = lats_idx_hi 
    if abs(pos[0]-lons[lons_idx_lo]) < abs(pos[0]-lons[lons_idx_hi]):
        lon_idx = lons_idx_lo 
    else:
        lon_idx = lons_idx_hi
    ERA5_dict.update(lat_close = lats[lat_idx],lon_close = lons[lon_idx], \
              lats_idx_lo = lats_idx_lo, lats_idx_hi = lats_idx_hi, \
              lons_idx_lo = lons_idx_lo, lons_idx_hi = lons_idx_hi)
    # Boolean array with True elements when ERA5T (expver[1]) is used. Skip these.
    ERA5T_mask = met.variables["t2m"][:,1,lat_idx,lon_idx].mask 
    # Pick out ERA5 (expver[0]) variables at chosen lat/lon and numpy array is cast
    # to list object.
    #=========== Smooth the arrays by computing a running average
    run = 16
    #=========== Time in seconds from 1/1/1900 at 00Hr 00m 00s
    time = smooth(list(met.variables["time"][ERA5T_mask]),run)
    #=========== Velocity in Eastward direction at 10 meters (m s-1)  
    u10 = smooth(list(met.variables["u10"][ERA5T_mask,0,lat_idx,lon_idx].data),run)
    #=========== Temperature at 2 meters (K)
    t2m = smooth(list(met.variables["t2m"][ERA5T_mask,0,lat_idx,lon_idx].data),run)    
    #=========== Velocity in Northward direction at 10 meters (m s-1)
    v10 = smooth(list(met.variables["v10"][ERA5T_mask,0,lat_idx,lon_idx].data),run)
    #=========== Dew point at 2 meters (K)
    d2m = smooth(list(met.variables["d2m"][ERA5T_mask,0,lat_idx,lon_idx].data),run)
    #=========== Low cloud cover (dimensionless 0 to 1)
    lcc = smooth(list(met.variables["lcc"][ERA5T_mask,0,lat_idx,lon_idx].data),run)
    #=========== Mean downward longwave radiative flux (W m-2) averaged over 1 hour
    #=========== which is the processing period for ERA5 reanalysis
    msdwlwrf = smooth(list(met.variables["msdwlwrf"][ERA5T_mask,0,lat_idx,lon_idx].data),run)
    #=========== Mean downward shortwave radiative flux (W m-2) averaged over 1 hour
    #=========== which is the processing period for ERA5 reanalysis
    msdwswrf = smooth(list(met.variables["msdwswrf"][ERA5T_mask,0,lat_idx,lon_idx].data),run)
    #=========== Surface Pressure (Pa or Kg m-1 s-2)
    sp = smooth(list(met.variables["sp"][ERA5T_mask,0,lat_idx,lon_idx].data),run)
    #=========== Total Cloud Cover (fraction 0 to 1)
    tcc = smooth(list(met.variables["tcc"][ERA5T_mask,0,lat_idx,lon_idx].data),run)
    #===========  Total Precipitation (kg m-2 s-1) ==================
    tp = smooth(list(met.variables["tp"][ERA5T_mask,0,lat_idx,lon_idx].data),run) 
    #=========== Convert seconds from 1/1/1900 to days from start. This lines up
    #=========== time with Flake
    t0 = time[0]
    time = [(t-t0)/24 for t in time]
    #=========== time is in units of days. Keep t0 (in hours since 1/1/1900 00h 00m 00s) for absolute
    #=========== time, which may be needed later.
    ERA5_dict.update({"t0_(sec)": t0, "time_step_(day)": time[1]-time[0]})
    #=========== Dewpoint (K) must be converted to vapor pressure (Pa)
    vp2m = [e_vp_from_DP(t, dp) for t, dp in zip(t2m, d2m)]
    #print(f'Vapor Pressure is {vp2m}\n')
    return ERA5_dict, time, u10, v10, vp2m, t2m, lcc, msdwlwrf, msdwswrf, sp, tcc, tp

def UT_select_lake():
    #=========== Some lake positions
    pos =[60.781662, 32.4200]
    #pos = [48.78, -72.20]
    #pos = [45.3292014252899, -80.8325204268411]
    #=========== Path for lake temperature file based on sattelite observations
    Tsfc_vs_time_DIR_fp = join(cwd(), "LAKES/GLOBO/LAKES_vs_time")
    #=========== Lake limnology file: Lake centers in lat/lon
    lake_lim_file_fp = join(cwd(), "LAKES/GLOBO/LAKES_limnology/globolakes-static_lake_centre_fv1.csv")
    #=========== A direct selection of lakes according to lake ID
    IDs = [3]
    #=========== A direct selection of lakes according to lake names
    names = ["lake Erie","Aral Sea"]
    #=========== A direct selection of lakes according to box boundaries (could be many lakes)    
    bounds = [{"lat_lo":60, "lat_hi": 63, "lon_lo":0, "lon_hi":20}]
    
    #=========== Build up an array of dicts containing lake information

    #lakes_dict = select_lake(Tsfc_vs_time_DIR_fp, lake_lim_file_fp, bounds_arrs = bounds)
    #print(f'lakes_dict from bounds: {lakes_dict}\n') 

    #lakes_dict = select_lake(Tsfc_vs_time_DIR_fp, lake_lim_file_fp, IDs = IDs)
    #print(f'lakes_dict from ID selection: {lakes_dict}\n') 
    
    lakes_dict = select_lake(Tsfc_vs_time_DIR_fp, lake_lim_file_fp, names = names)
    print(f'lakes_dict from names: {lakes_dict}\n') 

def UT_get_ERA5_met():
    #=========== Some lake positions
    pos =[60.781662, 32.4200]
    #pos = [48.78, -72.20]
    #pos = [45.3292014252899, -80.8325204268411] 
    ERA5_fp = "/Users/thomassangrey/Documents/Ocean/Flake/METEO/ERA5_NorEur/ERA5.nc"
    get_ERA5_met(ERA5_fp, pos)

def UT_get_daily_Lake_temp():
    #=========== Some lake positions.
    pos =[60.781662, 32.4200]
    #pos = [48.78, -72.20]
    #pos = [45.3292014252899, -80.8325204268411] 
    #=========== A relative path to daily temperature files of globe
    #=========== director is relative to Flake working directory
    rel_path = "LAKES/GLOBO/LAKES_daily"
    #=========== Example file from 1998 is chosen here
    file = "19980104120000-GloboLakes-L3S-LSWT-v4.0-fv01.0.nc"
    daily_lake_file_fp = join(cwd(), rel_path, file)
    get_daily_Lake_temp(daily_lake_file_fp, pos, 10, plot_area = True)

def UT_get_Lake_temp_vs_time():
    """
    The ncdump -h information yields :
    time_coverage_start = "19950611T120000Z" ;
    time_coverage_end = "20161126T120000Z" ;
    but this appears incorrect since the time_array[0] and time_array[-1],
    referenced to 1970-01-01-00-00-00<UTC> do not reveal these dates. They reveal
    the following dates:
    start_day = "1995-06-06-12-00-00"
    end_day = "2016-12-30-12-00-00"
    """
    startday = "1999-01-03-00-00-00"
    endday = "2016-12-30-12-00-00"
    lake_file = "LAKE00000012-GloboLakes-L3S-LSWT-v4.0-fv01.0.nc"
    rel_path = "LAKES/GLOBO/LAKES_vs_time"
    lake_file_fp = join(cwd(),rel_path, lake_file)
    #=========== <cavu_thresh> (float) defines visibility of atmosphere. 1 is unlimited, 0 is full cover.
    cavu_thresh = 0.05
    #=========== <QL_min> (int) defines lowest tolerable quality level (>=). Range from 1 to 5 
    QL_min = 3
    #=========== <unc_max> (float) defines highest tolerable uncertainty level (<=). Range from 0.0 to 1.0 
    unc_max = 1
    
    time, T_obs = get_Lake_temp_vs_time(lake_file_fp, cavu_thresh, \
                  QL_min, unc_max, startday = startday, endday = endday, \
                  plot_area = True)
    PeriodoGram_lake_time(time, T_obs)
#    time, T_obs = get_Lake_temp_vs_time(rel_path, lake_file, cavu_thresh, \
#                  plot_area = False)
    fig, ax = plt.subplots(1,1)
    ax.set_title("T_obs")
    label = f'T_obs'
    ax.plot(time, T_obs,label = label, marker = "d")
    #ax.plot(np.array(list(range(len(time)))), time,label = label)
    ax.set_xlabel(f'time (days)')
    ax.set_ylabel(f'T_obs')
    ax.legend(bbox_to_anchor=(1, 1.0), loc='upper left', prop={'size': 10})
    plt.show(block=True)
