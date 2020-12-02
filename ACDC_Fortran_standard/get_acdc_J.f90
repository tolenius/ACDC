module get_acdc_J

implicit none

contains

subroutine acdc_plugin(names_vapor,c_vapor,cs_ref,temp,ipr,t_sim,t_tot,j_acdc,diameter_acdc,c_out,c_bins)

use acdc_simulation_setup, only : get_system_size, get_vapor_indices                ! parameters related to the simulation system size
use acdc_simulation_setup, only : group_size_bins, nbins
use acdc_simulation_setup, only : solve_ss, use_solver                              ! logicals for the solution approaches
use acdc_system, only : small_set_mode, n_neutral_monomers, n_charges, nclust_max   ! simulation mode and numbers of vapors and charging states
use driver_acdc_J, only : acdc_driver

	implicit none
	
	! Input: ambient conditions and simulation time
    character(len=11), dimension(n_neutral_monomers), intent(in) :: names_vapor ! vapor names
	real(kind(1.d0)), intent(inout) :: c_vapor(n_neutral_monomers)              ! vapor concentrations (1/m^3)
	real(kind(1.d0)), intent(in) :: cs_ref				! reference coagulation sink (1/s)
	real(kind(1.d0)), intent(in) :: temp				! temperature (K)
	real(kind(1.d0)), intent(in) :: ipr					! ion production rate (1/s/m^3)
    real(kind(1.d0)), intent(in) :: t_sim, t_tot		! simulation time and total accumulated time (s)
	
	! Output: formed particles/time and their acid content
	real(kind(1.d0)), intent(out) :: j_acdc				! simulated formation rate (1/s/m^3)
	real(kind(1.d0)), intent(out) :: diameter_acdc		! mass diameter of the formed particles (m)
    ! Optional output
    real(kind(1.d0)), intent(out), optional :: c_out    ! concentration of outgrown particles
	real(kind(1.d0)), intent(out), optional :: c_bins(nbins)   ! size bin concentrations

	real(kind(1.d0)), save, allocatable :: c(:)			! cluster concentrations
	integer, save :: neq_syst = 0, nclust_syst = 0, nout_syst(n_charges) = 0	! parameters for the simulation system size
	real(kind(1.d0)), save :: diameter_max_syst
	integer, save :: n1vapor(n_neutral_monomers) = 0, nmol_vapor(n_neutral_monomers) = 0    ! indices for vapor monomers
	integer, save :: grouped_indices(0:nbins,nclust_max)! indices of the clusters belonging to each bin
	integer, save :: nclust_per_bin(0:nbins)			! total number of clusters in each bin

	real(kind(1.d0)) :: j_out(n_charges)				! formation rate vector (for neu only, or for neu, neg and pos)
	real(kind(1.d0)) :: t_in							! input time for the driver (s)
	logical :: int_ok = .true.							! .false. if integration fails
	real(kind(1.d0)), save :: t_iter = 1.d-16			! iteration time step (s) for the Euler method
	integer, save :: ipar(4)							! parameters for re-calling the monomer settings and rate constants
	logical, save :: firstcall = .true.
	integer :: i

	
	if (firstcall) then
	
		firstcall = .false.
		
		ipar = 0
		
		! Determine the system size and vapor indices
		call get_system_size(neq_syst,nclust_syst,nout_syst,diameter_max_syst)
		call get_vapor_indices(names_vapor,n1vapor,nmol_vapor)
		
		! Initialize the concentrations
		allocate(c(neq_syst))
		c = 0.d0
        
        if (present(c_out)) then
            c_out = 0.d0
        end if
		
        if (present(c_bins)) then
            if (.not. small_set_mode) then
                call group_size_bins(grouped_indices,nclust_per_bin)
            end if
            c_bins = 0.d0
        end if
        
        ! Write out brief information on the simulation settings
        
        write(*,*)
		
		write(*,*) 'Solving a set of ', neq_syst, ' equations'
		
		if (neq_syst .le. 5000) then
			use_solver = .true.
			write(*,*) 'Using VODE to integrate the equations'
		else
			use_solver = .false.
			write(*,*) 'Large system: using Eulerian integration instead of VODE'
		end if
        
        if (solve_ss) then
            write(*,*) 'Applying the steady-state assumption'
        end if
        
        write(*,*)
		
	end if

	! Initialize the rate constants etc. at every call because of the varying ambient conditions
	ipar(1) = 0
	ipar(3) = 0
	
	if (solve_ss) then
		! Override the input time with a maximum time that the driver is allowed to try to get to the steady state
		t_in = 1.d8
	else
		t_in = t_sim
	end if
	
	! Set the vapor concentrations
    c(n1vapor) = c_vapor

	! Run ACDC to the steady state / for the given time
	call acdc_driver(neq_syst,nclust_syst,nout_syst,c,cs_ref,temp,ipr,t_in,t_tot,t_iter,ipar,int_ok,j_out)
	j_acdc = sum(j_out)
	!write(*,*) 'Neutral, neg. and pos. J: ',j_out*1.d-6,' cm^-3 s^-1'
	if (j_acdc .lt. 0.d0) then
		if (j_acdc .lt. -1.d-12 .and. int_ok) then
			write(*,*) 'Warning: negative J = ',j_acdc*1.d-6,' cm^-3 s^-1, something wrong?'
		end if
		j_acdc = 1.d-100
	end if
	
	! Update the vapor concentrations, as clusters may be assumed to act as a sink for vapors
    c_vapor = c(n1vapor)
	
	! The first outgrown size is ~equal to the largest simulated size
	diameter_acdc = diameter_max_syst
    
    if (present(c_out)) then
        !c_out = sum(c(nout_syst)) ! Arbitrary for solve_ss = .true.
        c_out = c_out + j_acdc*t_sim
    end if
	
    if (present(c_bins)) then
        if (.not. small_set_mode) then
            ! Concentration in bins
            do i = 1,nbins
                if (nclust_per_bin(i) .ne. 0) then
                    c_bins(i) = sum(c(grouped_indices(i,1:nclust_per_bin(i))))
                end if
            end do
        end if
    end if
	
end subroutine acdc_plugin

end module get_acdc_J
