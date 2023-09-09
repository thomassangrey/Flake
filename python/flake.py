"""
Author: Thomas Sangrey
Date:   09/10/2023
Name: flake.py
Summary:
    Provides key driving components for running flake simulations and for preparing
    some meteorological input for the FLake simulator. Some flux variables like downware
    longwave atmospheric flux are caclculated here. The procedure assemble_inputs(...)
    prepares a ctypes connection between python and flake's fortran code.
"""
import f90nml
import ctypes as ct
import numpy as np
from math import exp as exp
from math import pi as pi
from math import sin as sin
from matplotlib import pyplot as plt
import datetime as date
from os import getcwd as cwd
from os import listdir
from os.path import isfile, join
import netCDF4
import Observations_util as Obs



T_w_freeze = 273.15    # Kelvin

def T_mean(Tm, h, D, CT, Tb):
    """
    Mironov 2008 identifies self-similar scaling relationship that constrains the
    variables listed here
    """
    T_mn = Tm - CT*(1 - h/D)*(Tm - Tb)
    return T_mn

def f(lat):
    """
    calculate f, the coriolis parameter
    """
    return 2*2*pi*sin(pi*lat/180)/86400


def e_vp_from_RH(T_a, RH_a):
    """
    From relative humidity (%). Returns Pa.
    Function returns water vapor pressure over water
    From Murray, F.W. 1966. ``On the computation of \
    Saturation Vapor Pressure'' J. Appl. Meteor. 6 p.204
    Q_a_rel is the relative humidity give as %
    T_a is atmospheric temperature given in Kelvin
    """
    a = 17.67
    b = 29.66
    e_s = 100*6.112*exp(a*(T_a-273.15)/(T_a-b) )
    e_a = Q_a*e_s/100
    return e_a

def e_vp_Q(T_a, Q_a, P_a):
	# From specific humidity (kg/Kg). Returns Pa.
    eps = 0.622           #molar mass ratio of vapor to dry air
    e_a = Q_a*P_a/(eps*(1-Q_a) + Q_a)
    return e_a

def q_spec(T_a, e_a, P_a):
	# From vapor pressure (Pa) to specific humidty (kg/Kg). Returns Pa.
    eps = 0.622           #molar mass ratio of vapor to dry air
    q_a = eps*e_a/(P_a + eps*e_a - e_a)
    return q_a

def make_q_arr(e_a_in_arr, T_a_in_arr, P_a_in_arr):
    """
    calculate specific humidity (g/kg) over water
    """
    q_arr = [q_spec(T_a,e_a, P_a) for T_a, e_a, P_a in zip(T_a_in_arr, \
    	       e_a_in_arr, P_a_in_arr)]
    return q_arr

def make_I_atm_lw_in_arr(e_a_in_arr, T_a_in_arr, CL_arr, P_a_in_arr):
    #=========== Use Mironov's SfcFlx_lwradatm function (modified as 
    #=========== f2py3_SfcFlx_lwradatm.f90) Use shared library to tunnel fortran
    #=========== code to python here.
    source = join(cwd(), "BIN/", 'FI_shared.so')
    fortlib = ct.CDLL(source) 
    lw_f = fortlib.f2py3_SfcFlx_lwradatm
    lw_f.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double)]
    Q_lw_atm_in_arr = [0 for i in range(len(e_a_in_arr))]
    for i in range(len(e_a_in_arr)):
        T_a = ct.c_double(T_a_in_arr[i])
        e_a = ct.c_double(e_a_in_arr[i])
        #+++++++++ DEBUG: Changed CL_low from default ct.c_double(-1)
        CL_low = ct.c_double(CL_arr[i])
        CL_tot = ct.c_double(CL_arr[i])
        Q_lw_atm = ct.c_double(Q_lw_atm_in_arr[i])
        lw_f(ct.byref(T_a), ct.byref(e_a),ct.byref(CL_low),\
    		               ct.byref(CL_tot),ct.byref(Q_lw_atm))
    #=========== I made the sign for Q_lw_atm explicitly negative so that downward
    #=========== flux is positive
        Q_lw_atm_in_arr[i] = -Q_lw_atm.value
    return Q_lw_atm_in_arr

def Flake_example_meteo(nml):
    time, I_atm_in_arr, T_a_in_arr, Q_a_in_arr,U_a_in_arr, \
    CL_tot_arr = Obs.get_nml_meteo(nml)   
    #=========== For now, let P_a_in be 1 bar
    P_a_in_arr = [1.013*10**5 for i in range(len(I_atm_in_arr))]
    T_a_in_arr = [T_a + T_w_freeze for T_a in T_a_in_arr]
    #=========== For e_a (vapor pressure) Convert mB to Pa by multiplying by 100
    e_a_in_arr = [Q_a_in*100 for Q_a_in in Q_a_in_arr]
    #=========== Get relative humidty in g/Kg
    q_a_in_arr = make_q_arr(e_a_in_arr, T_a_in_arr, P_a_in_arr)
    #=========== Compute downward long wave radiant flux from atmosphere
    #=========== Use a Flake subroutine to make this calculation
    Q_lw_atm_in_arr = make_I_atm_lw_in_arr(e_a_in_arr, T_a_in_arr, \
                   CL_tot_arr, P_a_in_arr)
    return time, I_atm_in_arr, T_a_in_arr, q_a_in_arr, U_a_in_arr, \
           CL_tot_arr, P_a_in_arr, Q_lw_atm_in_arr

def ERA5_meteo(nml, a_lake_dict):
    #=========== Get the first lake key (a_lake_dict is a single key-val pair)
    a_lake_key = list(a_lake_dict.keys())[0]
    print(f'The lake {a_lake_key} has information:\n{a_lake_dict[a_lake_key]}\n')
    #=========== Get the position of the lake's center
    lat = float(a_lake_dict[a_lake_key]["lake_centroid"]["lat"])
    lon = float(a_lake_dict[a_lake_key]["lake_centroid"]["lon"])
    pos = [lat, lon]
    print(f'Lake Coordinates = {pos}\n')
    #=========== meteo file's path name relative to current working directory
    meteo_file = nml['meteo']['meteofile']
    path = join(cwd(),meteo_file)
    #=========== Get a reanalysis variables to be used to force Flake
    #=========== Output variable names of get_ERA5_met:
    #=========== time, u10, v10, vp2m, t2m, lcc, msdwlwrf, msdwswrf, sp, tcc, tp
    ERA5_dict, time, u10, v10, vp2m, T_a_in_arr, lcc, Q_lw_atm_in_arr, I_atm_in_arr, \
    P_a_in_arr, CL_tot_arr, tp = Obs.get_ERA5_met(path, pos)
    U_a_in_arr   = [(u**2 + v**2)**(0.5) for u, v in zip(u10, v10)]
    q_a_in_arr   = make_q_arr(vp2m, T_a_in_arr, P_a_in_arr)
    #+++++++++++ DEBUG calculating lw radiation here instead of using the renalysis data
    #Q_lw_atm_in_arr = make_I_atm_lw_in_arr(vp2m, T_a_in_arr, \
    #               lcc, P_a_in_arr)
    ERA5_dict.update({"lat": lat, "lon": lon})
    return time, I_atm_in_arr, T_a_in_arr, q_a_in_arr, U_a_in_arr, \
           CL_tot_arr, P_a_in_arr, Q_lw_atm_in_arr, ERA5_dict

def assemble_inputs(nml_file, a_lake_dict = None):
    """
    Input: 1) <str> takes a file name of type nml to convey necessary lake
    configuration parameters. The nml file contains parameters as well as 
    file names for observational data and meteorological data. The specific
    type of meteo data available determines how the data should be processed
    for use by flake (via Flake_example_meteo(), ERA5_meteo(), or other function).
    For example, ERA5 provides wind velocity components but Flake_example (FE)
    provides wind speed. Also, FE provides vapor pressure in mB but no
    downward longwave radiation information, whereas ERA5 provides dew point
    information and longwave radiation data.
    2) optional dict(a_lake_dict) has the example form: 
        {' lake Victoria': {'lake_file': 'LAKE00000003-GloboLakes-L3S-LSWT-v4.0-fv01.0.nc', 
        'country': 'Tanzania', 'lake_centroid': {'lat': '-0.8764', 'lon': '33.1431'}, 
        'lake_box': {'minlat': '-3.0403', 'maxlat': '0.4903', 'minlon': '34.8653', 
        'maxlon': '34.8653'}} }
    which is one key-value pair of possibly many obtained from a more comprehensive
    lake_dict. Returns 
    Output:<dict(progs)> a dict of prognostic variables of interest for
    post-processing and analysis:
    "time": time
    "T_mnw_out_arr": [mean water temp]                              (K)
    "T_wML_out_arr": [mixed layer water temp]                       (K)
    "T_bot_out_arr": [water temp at bottom]                         (K)
    "h_ML_out_arr":  [mixed layer depth]                            (m)
    "T_sfc_n_arr":   [water surface temp]                           (K)
    "h_ice_out_arr": [ice thickness]                                (m)
    "Q_lw_atm_in_arr":  [downward atmospheric longwave radiation flux] (W m-2 s-1)
    "T_a_arr":       [air temp]                                     (K) 
    "I_atm_in_arr":  [downward solar shortwave radiation flux]      (W m-2 s-1)
    
    Optica_params and albedos are assumed, for the purposes here, to 
    be internal to Flake. No attempt is made to provide access to those
    values or make use of them. Flake assigns default values for
    these parameters. The reason is, for now at least, that these
    variables are not necessary for this part of the project and because
    optica_params are stored as derived types within Flake. Derived
    types can't be tunneled to python so a different arrangement would
    be required to make use of these variables at the python level.
    """    
              
    #=========== Retrieve dict() from a lake nml file for configuration
    nml = f90nml.read(nml_file)   

    #=========== dict references From a lake nml file =============================
    height_u_in = nml['METEO']['z_wind_m']
    height_tq_in =  nml['METEO']['z_taqa_m']
    depth_w = nml['LAKE_PARAMS']['depth_w_lk']     # I shut off sediment: lflk_botsed_use = .FALSE.
    fetch = nml['LAKE_PARAMS']['fetch_lk']         # meters
    depth_bs = nml['Lake_PARAMS']['depth_bs_lk']   #recommended by Flake test runs
    par_Coriolis = f(nml['LAKE_PARAMS']['latitude_lk']) # s^-1
    T_wML_in = nml['SIMULATION_PARAMS']['T_wML_in'] + T_w_freeze 
    T_bot_in = nml['SIMULATION_PARAMS']['T_bot_in'] + T_w_freeze      
    h_ML_in = nml['SIMULATION_PARAMS']['h_ML_in']      
    del_time = nml['SIMULATION_PARAMS']['del_time_lk']       
    H_B1_in = nml['Lake_PARAMS']['depth_bs_lk']
    T_bs = nml['Lake_PARAMS']['T_bs_lk'] + T_w_freeze # Temperature of outer edge
                                                      # of thermally active bottom
                                                      # sediment layer
    
    #=========== Retrieve code for choosing meteo file
    code = nml_file.split(".")[0].split("_")[-1]
    #=========== Make a choice regarding what meteo file to use. Choices vary
    #=========== between simple test cases (Mueggelsee80-96_FE.nml), or
    #=========== a reanalysis file, or non-gridded weather station data (USGS)
    #=========== Use a code such as FE, ERA5, USGS to distinquish:
    #=========== i.e., "foofile_FE.nml", "foofile_ERA5.nml", "foofile_USGS.nml"
    
    if code == "FE":              # Flake example meteo file
        time, I_atm_in_arr, T_a_in_arr, q_a_in_arr, U_a_in_arr, \
        CL_tot_arr, P_a_in_arr, Q_lw_atm_in_arr = Flake_example_meteo(nml)
    elif code == "ERA5":          # ERA5 reanalysis meteo file
        if a_lake_dict is not None:
            time, I_atm_in_arr, T_a_in_arr, q_a_in_arr, U_a_in_arr, \
            CL_tot_arr, P_a_in_arr, Q_lw_atm_in_arr, ERA5_dict = ERA5_meteo(nml, a_lake_dict)
            # time step is set by ERA5 (in days) and coverted to seconds for Flake
            del_time = 86400*ERA5_dict["time_step_(day)"]          
            par_Coriolis = f(ERA5_dict["lat"])
        else:
            raise Exception('Failure to supply lake_dict info.') 
    elif code == "USGS":          # Future meteo implementation
        pass
    else:
        raise Exception('Failure to open necessary meteo file')    
    #=========== Some more lake parameters ========================================
    dMsnowdt_in = 0.0
    T_snow_in = T_w_freeze        # Recommended by Flake test runs (website)
    T_ice_in =  T_w_freeze        # Celsius - ice temperature 
    T_B1_in = T_bs + T_w_freeze   # Celsius - bottom of upper layer sediment temperature
    C_T_in = 0.5            # Temperature Profile dimensionless
    h_snow_in = 0           # meters - snow thickness
    h_ice_in = 0            # meters - ice thickness
    T_mnw_in = T_mean(T_wML_in, h_ML_in, depth_w, C_T_in, T_bot_in)  #Mirinov 2008
    T_sfc_p = T_wML_in             # Initial Surface Temp (updated each time step)
    step_number = 1         # Do loop iteration number in Flake interface
    #=========== Initialize these so that they can be passed to fortran ============
    #=========== without generating a "reference before assignment" error ==========
    T_snow_out   = 0
    T_ice_out    = 0
    T_mnw_out    = 0
    T_wML_out    = 0
    T_bot_out    = 0
    T_B1_out     = 0
    C_T_out      = 0
    h_snow_out   = 0
    h_ice_out    = 0
    h_ML_out     = 0
    H_B1_out     = 0
    T_sfc_n      = 300

    T_mnw_out_arr = [0 for i in range(len(I_atm_in_arr))]
    T_wML_out_arr = [0 for i in range(len(I_atm_in_arr))]
    T_bot_out_arr = [0 for i in range(len(I_atm_in_arr))]
    C_T_out_arr   = [0 for i in range(len(I_atm_in_arr))]
    h_ML_out_arr  = [0 for i in range(len(I_atm_in_arr))]
    T_sfc_n_arr   = [0 for i in range(len(I_atm_in_arr))]
    h_ice_out_arr = [0 for i in range(len(I_atm_in_arr))]
    T_ice_out_arr = [0 for i in range(len(I_atm_in_arr))]
    
    source = join(cwd(), "BIN/", 'FI_shared.so')
    fortlib = ct.CDLL(source) 
    FI_w = fortlib.FI_wrapper
    FI_w.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),\
                     ct.POINTER(ct.c_double), ct.POINTER(ct.c_long)]

    for i in range(len(I_atm_in_arr)):
        dMsnowdt_in  = ct.byref(ct.c_double(dMsnowdt_in))
        I_atm_in     = ct.byref(ct.c_double(I_atm_in_arr[i]))
        Q_lw_atm_in  = ct.byref(ct.c_double(Q_lw_atm_in_arr[i]))
        height_u_in  = ct.byref(ct.c_double(height_u_in))
        height_tq_in = ct.byref(ct.c_double(height_tq_in))
        U_a_in       = ct.byref(ct.c_double(U_a_in_arr[i]))
        T_a_in       = ct.byref(ct.c_double(T_a_in_arr[i]))
        q_a_in       = ct.byref(ct.c_double(q_a_in_arr[i]))
        P_a_in       = ct.byref(ct.c_double(P_a_in_arr[i]))                    
        depth_w      = ct.byref(ct.c_double(depth_w))
        fetch        = ct.byref(ct.c_double(fetch))
        depth_bs     = ct.byref(ct.c_double(depth_bs))
        T_bs         = ct.byref(ct.c_double(T_bs))
        par_Coriolis = ct.byref(ct.c_double(par_Coriolis))
        del_time     = ct.byref(ct.c_double(del_time))
        T_snow_in    = ct.byref(ct.c_double(T_snow_in))
        T_ice_in     = ct.byref(ct.c_double(T_ice_in))
        T_mnw_in     = ct.byref(ct.c_double(T_mnw_in))
        T_wML_in     = ct.byref(ct.c_double(T_wML_in))
        T_bot_in     = ct.byref(ct.c_double(T_bot_in))
        T_B1_in      = ct.byref(ct.c_double(T_B1_in))
        C_T_in       = ct.byref(ct.c_double(C_T_in))
        h_snow_in    = ct.byref(ct.c_double(h_snow_in))
        h_ice_in     = ct.byref(ct.c_double(h_ice_in))
        h_ML_in      = ct.byref(ct.c_double(h_ML_in))
        H_B1_in      = ct.byref(ct.c_double(H_B1_in))
        T_sfc_p      = ct.byref(ct.c_double(T_sfc_p))
        T_snow_out   = ct.byref(ct.c_double(T_snow_out))
        T_ice_out    = ct.byref(ct.c_double(T_ice_out))
        T_mnw_out    = ct.byref(ct.c_double(T_mnw_out))
        T_wML_out    = ct.byref(ct.c_double(T_wML_out))
        T_bot_out    = ct.byref(ct.c_double(T_bot_out))
        T_B1_out     = ct.byref(ct.c_double(T_B1_out))
        C_T_out      = ct.byref(ct.c_double(C_T_out))
        h_snow_out   = ct.byref(ct.c_double(h_snow_out))
        h_ice_out    = ct.byref(ct.c_double(h_ice_out))
        h_ML_out     = ct.byref(ct.c_double(h_ML_out))
        H_B1_out     = ct.byref(ct.c_double(H_B1_out))
        T_sfc_n      = ct.byref(ct.c_double(T_sfc_n))
        step_number  = ct.byref(ct.c_long(step_number))
        
        
        FI_w( dMsnowdt_in, I_atm_in, Q_lw_atm_in,                                 \
                height_u_in, height_tq_in,                                        \
                U_a_in, T_a_in, q_a_in, P_a_in,                                   \
                                                                                  \
                depth_w, fetch, depth_bs, T_bs, par_Coriolis, del_time,           \
                T_snow_in,  T_ice_in,  T_mnw_in,  T_wML_in,     \
                T_bot_in,  T_B1_in, C_T_in,  h_snow_in,  h_ice_in,  h_ML_in,      \
                H_B1_in, T_sfc_p,                                                 \
                                                                                  \
                T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out, T_B1_out, \
                C_T_out, h_snow_out, h_ice_out, h_ML_out, H_B1_out, T_sfc_n, step_number)

        #=========== All updates of prognostic variables (setting _p variable
        #=========== to value at end of last time step) are completed in Flak_driver.
        #=========== However, there is one exception: T_sfc_p = T_sfc_n is explicitly
        #=========== performed here.
        dMsnowdt_in  = dMsnowdt_in._obj.value
        height_u_in  = height_u_in._obj.value
        height_tq_in = height_tq_in._obj.value
        depth_w      = depth_w._obj.value
        fetch        = fetch._obj.value
        depth_bs     = depth_bs._obj.value
        T_bs         = T_bs._obj.value            
        par_Coriolis = par_Coriolis._obj.value
        del_time     = del_time._obj.value
        T_snow_in    = T_snow_out._obj.value      
        T_ice_in     = T_ice_out._obj.value       
        T_mnw_in     = T_mnw_out._obj.value       
        T_wML_in     = T_wML_out._obj.value       
        T_bot_in     = T_bot_out._obj.value       
        T_B1_in      = T_B1_out._obj.value        
        C_T_in       = C_T_out._obj.value         
        h_snow_in    = h_snow_out._obj.value      
        h_ice_in     = h_ice_out._obj.value       
        h_ML_in      = h_ML_out._obj.value        
        H_B1_in      = H_B1_out._obj.value        
        T_sfc_p      = T_sfc_n._obj.value         
        T_sfc_n      = T_sfc_p
        T_snow_out   = T_snow_out._obj.value
        T_ice_out    = T_ice_out._obj.value
        T_mnw_out    = T_mnw_out._obj.value
        T_wML_out    = T_wML_out._obj.value
        T_bot_out    = T_bot_out._obj.value
        T_B1_out     = T_B1_out._obj.value
        C_T_out      = C_T_out._obj.value
        h_snow_out   = h_snow_out._obj.value
        h_ice_out    = h_ice_out._obj.value
        h_ML_out     = h_ML_out._obj.value
        H_B1_out     = H_B1_out._obj.value
        step_number  = step_number._obj.value

        #=========== T_mnw_in = T_mean(T_wML_in, h_ML_in, depth_w, C_T_in, T_bot_in)  #Necessary ??? Mirinov 2008

        T_mnw_out_arr[i] = T_mnw_out
        T_wML_out_arr[i] = T_wML_out
        T_bot_out_arr[i] = T_bot_out
        C_T_out_arr[i]   = C_T_out_arr
        h_ML_out_arr[i]  = h_ML_out
        T_sfc_n_arr[i]   = T_sfc_n
        h_ice_out_arr[i] = h_ice_out
        T_ice_out_arr[i] = T_ice_out
    title = "PROGS"
    progs = dict({"title": title, "time": time, "T_mnw_out_arr": T_mnw_out_arr, \
                  "T_wML_out_arr": T_wML_out_arr, "T_bot_out_arr": T_bot_out_arr, \
                  "h_ML_out_arr": h_ML_out_arr, "T_sfc_n_arr": T_sfc_n_arr, \
                  "h_ice_out_arr": h_ice_out_arr, "T_ice_out_arr": T_ice_out_arr})
    title = "METEO"
    meteo = dict({"title": title, "time": time, "I_atm_in_arr": I_atm_in_arr, \
                 "T_a_in_arr": T_a_in_arr, "q_a_in_arr": q_a_in_arr, \
                 "U_a_in_arr": U_a_in_arr, "CL_tot_arr": CL_tot_arr, \
                 "P_a_in_arr": P_a_in_arr, "Q_lw_atm_in_arr": Q_lw_atm_in_arr })
    return progs, meteo


def UT_assemble_inputs_FE():
    #=========== Some lake nml files for individual lakes ----------- #
    #nml_file = 'Stechlin94-98.nml'
    nml_file = 'Flake_example_1_FE.nml'
    #nml_file = 'Mueggelsee_2017.nml'
    nml_file = join(cwd(), "Flake_nmls", nml_file)
    #=========== Assemble nml info for lake and develop prognostic variables vs
    #=========== time for lake (chosen from nml file). Plot the prognostic variables   
    progs, meteo = assemble_inputs(nml_file)       # assembled prognostic variables vs time for plotting
    plt1 = Obs.plot_lake_matrix(progs)
    plt1.show(block=False)
    plt2 = Obs.plot_lake_matrix(meteo)
    plt2.show(block=False)
    return plt2

def UT_assemble_inputs_ERA5():
    """
    Inputs: none
    Outputs: none
    A unit test of the main lake model driven by ERA5 meteorological data. 
    Visualize the metorological content and the sattelite-derived skin temperature.

    """
    #=========== Some lake nml files for individual lakes ----------- #
    nml_file = 'Flake_config_ERA5.nml'
    file_path = join(cwd(), "Flake_nmls", nml_file)
    #=========== Assemble nml info for lake and develop prognostic variables vs
    #=========== time for lake (chosen from nml file). Plot the prognostic variables   
    #=========== Path for lake temperature file based on sattelite observations
    Tsfc_vs_time_rel_path = "LAKES/GLOBO/LAKES_vs_time"
    #=========== Lake limnology file: Lake centers in lat/lon
    lake_lim_file = "LAKES/GLOBO/LAKES_limnology/globolakes-static_lake_centre_fv1.csv"
    lakes_dict = Obs.select_lake(Tsfc_vs_time_rel_path, lake_lim_file,\
        bounds_arrs = [{"lat_lo":52, "lat_hi": 63, "lon_lo":10, "lon_hi":20}])
    #print(f'lakes_dict from bounds: {lakes_dict}\n') 
    progs, meteo = assemble_inputs(file_path, lakes_dict )       # assembled prognostic variables vs time for plotting
    plt1 = Obs.plot_lake_matrix(progs)
    plt1.show(block=False)
    plt2 = Obs.plot_lake_matrix(meteo)
    plt2.show(block=False)
    return plt2

def main():
    #Obs.UT_select_lake()
    #Obs.UT_get_Lake_temp_vs_time()
    #plt_ERA5 = UT_assemble_inputs_ERA5()
    #plt_ERA5.show(block=True)
    plt_FE = UT_assemble_inputs_FE()
    plt_FE.show(block=True)
    

if __name__ == "__main__":
    main()
