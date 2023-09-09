! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------
!---Thomas Sangrey --- 7/24/2023
!---Creates a simple wrapper for flake_interface(...) that is absent any
!---derved types. Its purpose is for it and a python main() program to be 
!---converted to a shared library using the f2py3 command line tool. 
!---In this way, the flake_interface_wrapper can be imported as a module
!---inside a python program.
!---
!---Removed albedo_ice albedo_water, albedo_snow, and each opticpar_medium
!---variable from the flake_interface(...) INTENT(INOUT) variable list. 
!---These are now local variables within flake_interface. I did this
!---because derived types can't be tunneled to python and because 
!---these variables were set within flake_interface() anyway. 

SUBROUTINE FI_wrapper( dMsnowdt_in, I_atm_in, Q_atm_lw_in, &
                             height_u_in, height_tq_in,     &
                             U_a_in, T_a_in, q_a_in, P_a_in,                                    &
                             
                             depth_w, fetch, depth_bs, T_bs, par_Coriolis, del_time,            &
                             T_snow_in,  T_ice_in,  T_mnw_in,  T_wML_in,  T_bot_in,  T_B1_in,   &
                             C_T_in,  h_snow_in,  h_ice_in,  h_ML_in,  H_B1_in, T_sfc_p,        &
                             
                             T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out, T_B1_out,  & 
                             C_T_out, h_snow_out, h_ice_out, h_ML_out, H_B1_out, T_sfc_n, step_number )      &
                             bind(c, name='FI_wrapper')

USE iso_c_binding
USE data_parameters , ONLY : &
     ireals,                 & ! KIND-type parameter for real variables
     iintegers                 ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE


!==============================================================================
!
! Declarations

!  Input (procedure arguments)

INTEGER(c_long), INTENT(IN) ::   &
  step_number                           ! Controls the DO loop for days of simulation

REAL (c_double), INTENT(IN) ::   &
  dMsnowdt_in                       , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
  I_atm_in                          , & ! Solar radiation flux at the surface [W m^{-2}]
  Q_atm_lw_in                       , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  height_u_in                       , & ! Height above the lake surface where the wind speed is measured [m]
  height_tq_in                      , & ! Height where temperature and humidity are measured [m]
  U_a_in                            , & ! Wind speed at z=height_u_in [m s^{-1}]
  T_a_in                            , & ! Air temperature at z=height_tq_in [K]
  q_a_in                            , & ! Air specific humidity at z=height_tq_in
  P_a_in                                ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]

REAL (c_double), INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  fetch                             , & ! Typical wind fetch [m]
  depth_bs                          , & ! Depth of the thermally active layer of the bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of the bottom sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  del_time                              ! The model time step [s]

REAL (c_double), INTENT(IN)  :: &
  T_snow_in                        , & ! Temperature at the air-snow interface [K] 
  T_ice_in                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_in                         , & ! Mean temperature of the water column [K]
  T_wML_in                         , & ! Mixed-layer temperature [K]
  T_bot_in                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_in                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_in                           , & ! Shape factor (thermocline)
  h_snow_in                        , & ! Snow thickness [m]
  h_ice_in                         , & ! Ice thickness [m]
  h_ML_in                          , & ! Thickness of the mixed-layer [m]
  H_B1_in                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_p                              ! Surface temperature at the previous time step [K]  

!  Output (procedure arguments)

REAL (c_double), INTENT(OUT)  :: &
  T_snow_out                        , & ! Temperature at the air-snow interface [K] 
  T_ice_out                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_out                         , & ! Mean temperature of the water column [K]
  T_wML_out                         , & ! Mixed-layer temperature [K]
  T_bot_out                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_out                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_out                           , & ! Shape factor (thermocline)
  h_snow_out                        , & ! Snow thickness [m]
  h_ice_out                         , & ! Ice thickness [m]
  h_ML_out                          , & ! Thickness of the mixed-layer [m]
  H_B1_out                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_n                               ! Updated surface temperature [K]  

!  Local counter index
INTEGER (kind=iintegers):: &
  i_ctr

DO i_ctr = 1, step_number
    CALL flake_interface( dMsnowdt_in, I_atm_in, Q_atm_lw_in, height_u_in, height_tq_in,           &
                             U_a_in, T_a_in, q_a_in, P_a_in,                                   &

                             depth_w, fetch, depth_bs, T_bs, par_Coriolis, del_time,            &
                             T_snow_in,  T_ice_in,  T_mnw_in,  T_wML_in,  T_bot_in,  T_B1_in,   &
                             C_T_in,  h_snow_in,  h_ice_in,  h_ML_in,  H_B1_in, T_sfc_p,        &

                             T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out, T_B1_out,  & 
                             C_T_out, h_snow_out, h_ice_out, h_ML_out, H_B1_out, T_sfc_n )
END DO

END SUBROUTINE FI_wrapper