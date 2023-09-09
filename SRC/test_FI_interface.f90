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

SUBROUTINE FI_wrapper( dMsnowdt_in, T_sfc_n ) bind(c, name='FI_wrapper')

USE iso_c_binding


!==============================================================================

IMPLICIT NONE



!==============================================================================
!
! Declarations

!  Input (procedure arguments)

REAL (c_double), INTENT(INOUT) ::  dMsnowdt_in        ! The rate of snow accumulation [kg m^{-2} s^{-1}
  

!  Output (procedure arguments)

REAL (c_double), INTENT(OUT)  :: T_sfc_n           ! Updated surface temperature [K]  

dMsnowdt_in = 1.0
T_sfc_n = -1.0
CALL sum3( dMsnowdt_in, T_sfc_n )
print *, "T_sfc_n = ", T_sfc_n

END SUBROUTINE FI_wrapper