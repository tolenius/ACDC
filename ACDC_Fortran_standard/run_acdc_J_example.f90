program run_acdc_J_example

use acdc_simulation_setup, only : n_neutral_monomers                ! numbers of vapors and size bins
use get_acdc_J, only : acdc_plugin

	implicit none

	! All variables in SI units
    character(len=11), dimension(n_neutral_monomers) :: names_vapor ! vapor names
	real(kind(1.d0)) :: c_vapor(n_neutral_monomers)                 ! vapor concentrations (m^-3)
	real(kind(1.d0)) :: cs_ref			                            ! coagulation sink (s^-1)
	real(kind(1.d0)) :: temp			                            ! temperature (K)
	real(kind(1.d0)) :: ipr				                            ! ion production rate (m^-3 s^-1)
    real(kind(1.d0)) :: j_acdc			                            ! simulated formation rate (m^-3 s^-1)
	real(kind(1.d0)) :: diameter_acdc	                            ! mass diameter of the formed particles (m)
    
    character(len=11), dimension(n_neutral_monomers) :: str_vapor
    integer :: i
	
    
    ! Simulation conditions
    
    names_vapor(1)(:) = 'A'
	names_vapor(2)(:) = 'N'
	c_vapor = (/1.d7, 1.d9/)*1.d6
	
	cs_ref = 1.d-3
    temp = 280.d0
	ipr = 3.d0*1.d6
    
!--------------------------------------------------------------------------------------------------
    
    ! Make some output neater
    do i = 1,size(names_vapor)
        str_vapor(i)(:) = '['//trim(names_vapor(i)(:))//'] (cm^-3)'
    end do
    
    
    ! Get the formation rate for the given system and conditions
    
    ! All in- and output in SI units
    ! The time variables are only dummies in this example as no time-dependent simulations are done
    call acdc_plugin(names_vapor,c_vapor,cs_ref,temp,ipr,60.d0,0.d0,j_acdc,diameter_acdc)
    
    ! Values in cm^-3 (easily readable)
    write(*,'(*(a20,1x))') str_vapor(:)(:), 'CS_ref (s^-1)', 'T (K)', 'IPR (cm^-3 s^-1)', 'J (cm^-3 s^-1)'
    write(*,'(*(es20.5e4,1x))') c_vapor*1.d-6, cs_ref, temp, ipr*1.d-6, j_acdc*1.d-6
	
	
end program run_acdc_J_example
