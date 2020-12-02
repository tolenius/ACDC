module driver_acdc_J
use acdc_system, only : small_set_mode, n_charges		! simulation mode and numbers of charging states
use acdc_system, only : n_variables, variable_temp		! information on parameters that vary based on input
use acdc_simulation_setup, only : solve_ss, j_flux, use_solver	! logicals for the solution approaches
use acdc_simulation_setup, only : sources_and_constants ! for determining fitted concentrations
use solution_settings

implicit none

contains

subroutine acdc_driver(neqn,nclust,nout_all,c,cs_ref,temperature,ipr,t_max,t_tot,t_iter,ipar,ok,j_out)

	implicit none
	
	! Input and output
	! Cluster distribution
	integer, intent(in) :: neqn, nclust, nout_all(n_charges)	! total length of the c vector, the number of actual clusters and indices for outgoing flux
	real(kind(1.d0)), intent(inout) :: c(neqn)			! initial concentrations -> final concentrations (1/m^3)
	! Ambient conditions
	real(kind(1.d0)), intent(in) :: cs_ref				! reference coagulation sink (1/s)
	real(kind(1.d0)), intent(in) :: temperature			! temperature (K)
	real(kind(1.d0)), intent(in) :: ipr					! ion production rate (1/s/m^3)
	! Simulation settings and outcome
	real(kind(1.d0)), intent(in) :: t_max, t_tot		! simulation time and total accumulated time (s)
	real(kind(1.d0)), intent(inout) :: t_iter       	! iteration time step for the Euler method (s)
	integer, intent(inout) :: ipar(4)					! parameters for re-calling the monomer settings and rate constants
	logical, intent(out) :: ok							! .false. if integration failed
	real(kind(1.d0)), intent(out) :: j_out(n_charges)	! simulated formation rates (neutral, neg and pos) (1/s/m^3)

	! Variables and parameters for the simulation
	real(kind(1.d0)) :: c00(neqn), ci(neqn), c_1(neqn), ci_tmp(neqn), f_out(neqn), f_out_mp(neqn), eps_max	! temporary c and dc/dt
	real(kind(1.d0)) :: c_out(n_charges)				! total concentration of outgrown particles at the beginning
	integer, parameter :: nt = 20, n_ss = 3			    ! max. numbers of timeranges for integration by the solver
	real(kind(1.d0)) :: t(nt)           				! start and end times
	real(kind(1.d0)) :: t00, t0, tfin, ti, tt(2), t_tmp
	real(kind(1.d0)) :: t_iter_tot			            ! the cumulative iterated time (s)
	real(kind(1.d0)) :: t_iter_tmp			            ! temporary iterated time for monitoring the convergence to steady state
	logical :: save_ci_tmp					        	! logical used in monitoring the convergence
	real(kind(1.d0)) :: parameters(n_variables)			! one or more of the following parameters (max. 4 param.):
														! temperature (1 param.), coagulation sink (1) and ion source rate (2; neg and pos)
	real(kind(1.d0)) :: parameters_def(4)
	integer :: i, j, k, n, nind, nt_tot
	
	! Parameters for the solver
	integer, parameter :: itol = 1						! same absolute tolerance for all components
	integer, parameter :: itask = 4						! stop at or beyond t=TOUT, but not beyond t_max, and return
	integer, parameter :: iopt = 1						! allow using optional input in the solver
	integer, parameter :: mf = 22						! full jacobian, computed numerically
	real(kind(1.d0)), parameter :: rtol = rtol_solver	! relative tolerance in the solver / error tolerance per step
														!(EPS) in the iteration
	real(kind(1.d0)), parameter :: atol = atol_solver	! absolute tolerance in the solver / iteration when assessing EPS
	integer :: lwork, liwork
	!integer(kind=selected_int_kind(10)) :: lwork, liwork
	integer :: istate, iwork(neqn+30)
	real(kind(1.d0)) :: work(22+9*neqn+2*neqn**2)
	
	real(kind(1.d0)) :: j_tot, j_by_charge(4), j_by_cluster(neqn), j_all(neqn,4) ! output of the formation rate subroutine
	external feval, jeval, formation					! subroutines for equations, jacobian and formation rate


	! Initialize
	t0 = t_tot
	t00 = t0
	call fix_fitted_c(neqn,nclust,c,c,ipar)						! Fix the concentrations of the fitted species
	c00 = c
	j_out = 0.d0

	if (variable_temp) then
		parameters_def = (/temperature, cs_ref, ipr, ipr/)
	else
		parameters_def = (/cs_ref, ipr, ipr, 0.d0/)
	end if
	do i=1,size(parameters)
		parameters(i) = parameters_def(i)						! Values for the varied parameters
	end do

	! Save the current concentrations of outgrown particles
	c_out = c(nout_all)

	! Use either the solver or the Eulerian iteration
	if (use_solver) then
    
        t = 0.d0

		if (solve_ss) then
            i = 4
            t(1:i) = (/1.d-8, 1.d-2, sstimetot-sstimech, sstimetot/)
            
            do while ((i .lt. nt) .and. (t(i) .lt. t_max))
                t_tmp = 10.d0*t(i)
                
                if (nt-i .le. n_ss+1) then
                    t_tmp = max(t_max-real(n_ss,kind=kind(1.d0))*sstimech, t_tmp)
                end if
                
                j = -1
                do while ((i .lt. nt) .and. (t(i) .lt. t_max) .and. (j .lt. n_ss))
                    i = i+1
                    j = j+1
                    t(i) = t_tmp + real(j,kind=kind(1.d0))*sstimech
                end do
			end do
        else
			if (t00 .lt. 1.d-4) then
                i = 2
				t(1:i) = (/min(1.d-8,1.d-1*t_max), t_max/)
			else
                i = 1
				t(i) = t_max
			end if
        end if
        
        t = t00+t
        nt_tot = i
        
		lwork = 22+9*neqn+2*neqn**2
		liwork = neqn+30
		work(5:10) = 0.d0						! use default parameter values in the solver
		work(1) = t00+t_max						! don't go beyond the end time
		iwork(5:10) = 0							! use default parameter values in the solver
		iwork(6) = 1000000						! allow more steps in the solver (probably not needed)
		istate = 1								! tells the solver that this is the first call
												! -> 2 in the solver after succesful integration		
	
		! break the simulation into several time intervals to improve convergence		
		do i=1,nt_tot
		
			ci = c												! save the previous concentrations
			ti = t0
			
			!tfin = min(max(t(i),t0*10.d0),t00+t_max)
			tfin = min(t(i),t00+t_max)
            
			ok = .true.
			call DVODE (feval,neqn,c,t0,tfin,itol,rtol,atol,&	! calling the solver
				& itask,istate,iopt,work,lwork,iwork,liwork,&
				& jeval,mf,parameters,ipar)
			if (istate .ne. 2) then								! checking everything went ok - if not, print error
				!write(*,*) 'ERROR: returned istate =', istate
				ok = .false.
			elseif (minval(c) .lt. negtol) then					! checking all concentrations > 0 - if not, print error
				!write(*,*) 'Negative concentrations: c_min = ', minval(c),' at t = ',t0
				ok = .false.
			end if
			if (.not. ok) then									! if there was a problem, try a shorter time interval
				!write(*,*) 'Trying to continue integration from the previous time point with a shorter time interval.'
				istate = 1
				c = ci
				t0 = ti
				tt(1) = t0+(tfin-t0)/1.d2
				tt(2) = tfin
				do j = 1, 2
					call DVODE (feval,neqn,c,t0,tt(j),itol,rtol,atol,&	! calling the solver
						& itask,istate,iopt,work,lwork,iwork,liwork,&
						& jeval,mf,parameters,ipar)
					if ((istate .ne. 2) .or. (minval(c) .lt. negtol)) then
						write(*,*) 'ERROR: returned istate =', istate
						write(*,*) 'Exiting the driver'
						t0 = t00
						c = c00
						j_out = 0.d0
						return									! if there is still a problem, give up
					else
						ok = .true.
					end if
				end do
			end if
			
			call fix_fitted_c(neqn,nclust,c,ci,ipar)			! Fix the concentrations of the fitted species
			
			! Check the convergence for the steady state case
			if (solve_ss) then
				if ((t0-t00 .ge. sstimetot) .and. (maxval(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol))) .le. sstol)) then
					!write(*,*) 'converged at t=',t0
					exit		! end the simulation when a steady state is reached
				end if
				if ((t0-t00.ge.t_max) .or. (i.eq.nt_tot)) then
                    write(*,*) "#################################"
					write(*,*) "### Steady state not reached! ###"
                    write(*,*) "#################################"
					!write(*,*) "monomer concentrations: ", c(n_monomers)*1.d-6
					t0 = t00
					c = c00
					j_out = 0.d0
					return
				end if
			end if
			
		end do
		
	else
	
		t_iter_tot = t00
        tt(1) = t_iter

        if (solve_ss) then
            t_iter_tmp = 0.d0
            save_ci_tmp = .true.
        end if

		do while (t_iter_tot .lt. t00+t_max)

			ci = c						! save the previous concentrations
			ti = t_iter_tot

			if (solve_ss .and. save_ci_tmp) then
				ci_tmp = ci
				save_ci_tmp = .false.
			end if
					
			tt(1) = min(tt(1),t00+t_max-t_iter_tot)
            ok = .true.
            call feval(neqn,ti,ci,f_out,parameters,ipar)
			
			! Adjust the time step according to the largest estimated change, if needed
			tt(2) = minval(chmax/abs(f_out(1:nclust)/max(ci(1:nclust),chtol)))
			tt(1) = min(tt(1),tt(2))
			
			do n = 1,neqn
				! 1-step Euler
				c_1(n) = ci(n) + f_out(n)*tt(1)
				! 2-step Euler at midpoint
				c(n) = ci(n) + f_out(n)*tt(1)/2.d0
			end do
			call feval(neqn,ti+tt(1)/2.d0,c,f_out_mp,parameters,ipar) ! f at midpoint
			do n = 1,neqn
				! 2-step Euler at ti + tt(1)
				c(n) = c(n) + f_out_mp(n)*tt(1)/2.d0
			end do
			
			! Estimate for the largest error per step relative to ci
			eps_max = maxval(abs((c_1(1:nclust)-c(1:nclust))/max(ci(1:nclust),atol)))
			! Estimate for the final solution based on the 1- and 2-step solutions and
			! the estimated error per step
			do n = 1,neqn
				c(n) = 2.d0*c(n)-c_1(n)
			end do
            t_iter_tot = t_iter_tot + tt(1)

			! See if the solution is acceptable
			! Checking for negative concentrations
            if (minval(c) .lt. negtol) then
				write(*,*) 'Negative concentrations at t = ',t_iter_tot
				!nind = minloc(c,1)
				!write(*,*) 'Lowest c at ',nind,ci(nind),c(nind)
				ok = .false.
			! Checking for estimated local errors wrt ci larger than accepted
			elseif (eps_max .gt. rtol) then
				!nind = maxloc(abs((c_1(1:nclust)-c(1:nclust))/max(ci(1:nclust),atol)),1)
				!write(*,*) 'Largest EPS at ',nind,ci(nind),c(nind)
				ok = .false.
			! Checking for changes in c larger than accepted (~ stability of the solution)
            elseif (maxval(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol))) .gt. chmax) then
                !write(*,*) 'Large changes in concentrations at t = ',t_iter_tot
				!nind = maxloc(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol)),1)
                !write(*,*) 'Largest rel. change at ',nind,ci(nind),c(nind)
                ok = .false.
			end if
            ! If there was a problem, try a shorter time interval
			if (.not. ok) then
				!write(*,*) 'Trying to continue integration from the previous time point with a shorter interval at t = ',t_iter_tot
                do j = 1,40
					c = ci
                    t_iter_tot = ti

					! Take the initial guess for the time step based on eps or simply halve the previous step, whatever is smaller
					tt(1) = tt(1)/2.d0
					if (eps_max .gt. 0.d0) then
						tt(1) = min(0.8d0*sqrt(rtol/eps_max)*tt(1)*2.d0,tt(1))
					end if

					do n = 1,neqn
						! 1-step Euler
						c_1(n) = ci(n) + f_out(n)*tt(1)
						! 2-step Euler at midpoint
						c(n) = ci(n) + f_out(n)*tt(1)/2.d0
					end do
					call feval(neqn,ti+tt(1)/2.d0,c,f_out_mp,parameters,ipar) ! f at midpoint
					do n = 1,neqn
						! 2-step Euler at ti + tt(1)
						c(n) = c(n) + f_out_mp(n)*tt(1)/2.d0
					end do

					eps_max = maxval(abs((c_1(1:nclust)-c(1:nclust))/max(ci(1:nclust),atol)))
					do n = 1,neqn
						c(n) = 2.d0*c(n)-c_1(n)
					end do
                    t_iter_tot = t_iter_tot + tt(1)
					
                    ! If there is still a problem, give up after a few tries, otherwise continue
					if ((.not. minval(c) .lt. negtol) .and.&
						&(.not. eps_max .gt. rtol) .and.&
						&(.not. maxval(abs((c(1:nclust)-ci(1:nclust))/max(ci(1:nclust),chtol)))&
						&.gt. chmax)) then
					    ok = .true.
                        exit
				    end if
                end do
                if (.not. ok) then
					write(*,*) 'Iteration still not ok at t = ',t_iter_tot
					write(*,*) 'Exiting the driver'
					t0 = t00
					c = c00
					j_out = 0.d0
					return
                end if
			end if
			
			call fix_fitted_c(neqn,nclust,c,ci,ipar)					! Fix the concentrations of the fitted species

			if (solve_ss) then
				t_iter_tmp = t_iter_tmp + tt(1)
				! Check the convergence
				if (t_iter_tmp .ge. sstimech) then
					if ((maxval(abs((c(1:nclust)-ci_tmp(1:nclust))/max(ci_tmp(1:nclust),chtol))) .le. sstol)&
						&.and. (t_iter_tot-t00 .ge. sstimetot)) then
                        ! If the concentrations haven't changed during the max. value given for t_iter_tmp,
                        ! consider them converged...
			            !write(*,*) 'Converged at t=',t_iter_tot
			            exit		! end the simulation when a steady state is reached
				    else
					    ! ... otherwise save the current concentrations and
					    ! check again after t_iter_tmp has reached the given value
					    t_iter_tmp = 0.d0
					    save_ci_tmp = .true.
		            end if
			    end if
			    if (t_iter_tot .ge. t00+t_max) then
                    write(*,*) "#########################"
                    write(*,*) "### Did not converge! ###"
                    write(*,*) "#########################"
				    !write(*,*) "monomer concentrations: ", c(n_monomers)*1.d-6, " cm^-3"
				    t0 = t00
					c = c00
					j_out = 0.d0
				    return
			    end if
			end if
			
			! Set the next time step based on eps
			if (eps_max .gt. 0.d0) then
				tt(1) = 0.8d0*sqrt(rtol/eps_max)*tt(1)
				!write(*,*) 'Changed the iteration time step to ',tt(1),' s based on the estimated error'
			end if
			
			t_iter = tt(1)				! save the time step for output

		end do
		
		t0 = t_iter_tot
	
	end if
	
	! Formation rate
	if (solve_ss .or. j_flux) then
		! Steady state or J directly from the equations: flux out of the system
		if (small_set_mode) then
			call formation(neqn,c,j_tot,j_by_charge,j_by_cluster,j_all,ipar,parameters)
			do i=1,n_charges
				j_out(i) = j_by_charge(i)
			end do
			j_out(1) = j_out(1)+j_by_charge(4)	! Recombination products
		else
			call formation(neqn,c,j_out,ipar,parameters)
		end if
	else
		! Fixed time step: clusters grown out of the system per integrated time
		j_out = (c(nout_all)-c_out)/(t0-t00)
	end	if
	
end subroutine acdc_driver

subroutine fix_fitted_c(neqn,nclust,c,ci,ipar)

	implicit none

	integer, intent(in) :: neqn, nclust					! sizes of c and fitted
	real(kind(1.d0)), intent(inout) :: c(neqn)			! solved concentrations after integration -> fixed concentrations (1/m^3)
	real(kind(1.d0)), intent(in) :: ci(neqn)			! initial concentrations before integration
	integer, intent(in) :: ipar(4)						! parameters for re-calling settings after they've been changed

	logical, save :: firstcall = .true.
	real(kind(1.d0)) :: source(nclust)
	logical :: isconst(nclust)
	integer, save, allocatable :: fitted(:,:)
	integer :: k, n
	real(kind(1.d0)) :: c_fit
	
	
	if (firstcall) then
		firstcall = .false.
		allocate(fitted(nclust,0:nclust))
		call sources_and_constants(nclust,source,isconst,fitted) ! get the fitting settings (variable "fitted")
	end if
	
	do n = 1, fitted(1,0)	! loop over the clusters whose concentrations are fitted (e.g. sum([1AnD])=const.)
		k = fitted(n+1,0)	! this concentration is fitted
		if (ipar(1) .eq. 0) then
			c_fit = ci(k)	! if the settings have changed, the sum([1AnD])=c_fit has been fed in as [1A] in the driver call;
							! feval will change ipar(1) to 1 when the integration begins...
		else
			c_fit = sum(ci((/k, fitted(n+1,2:fitted(n+1,1)+1)/))) ! ...after which c_fit is determined from the previous concentrations
		end if
		c(k) = c_fit-sum(c(fitted(n+1,2:fitted(n+1,1)+1)))
		if (c(k).lt.0.d0) then
			c(k) = 0.d0
			c(fitted(n+1,2:fitted(n+1,1)+1)) = c(fitted(n+1,2:fitted(n+1,1)+1))*&
				&c_fit/sum(c(fitted(n+1,2:fitted(n+1,1)+1)))
		end if
	end do
	
end subroutine fix_fitted_c


end module driver_acdc_J
