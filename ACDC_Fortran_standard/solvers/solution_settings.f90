module solution_settings

implicit none

    ! Relative tolerance in the solver or when assessing the iteration solutions (see driver)
    real(kind(1.d0)), parameter :: rtol_solver = 1.d-5
	!real(kind(1.d0)), parameter :: rtol_solver = 1.d-3	! 1.d-3 is the default tolerance e.g. in Matlab
	
	! Settings for a solver
	real(kind(1.d0)), parameter :: atol_solver = 1.d-6	! absolute tolerance in the solver
	
    ! Settings for Eulerian iteration
    real(kind(1.d0)), parameter :: chmax = 1.d-2	    ! threshold relative change in the concentration for
														! determining if the iteration time step should be decreased
    real(kind(1.d0)), parameter :: chtol = 1.d-6        ! lowest conc. (m^-3) to consider when comparing the changes
	
	! Lowest negative concentration accepted
	real(kind(1.d0)), parameter :: negtol = -1.d-6		! (m^-3)
	
	! Criteria for monitoring the convergence to a steady state
	real(kind(1.d0)), parameter :: sstol = 1.d-5		! maximum relative change in the concentrations
	real(kind(1.d0)), parameter :: sstimech = 6.d2		! minimum time interval for checking the changes (s)
	real(kind(1.d0)), parameter :: sstimetot = 1.2d3	! minimum total simulation time required for a steady state (s)


end module solution_settings
