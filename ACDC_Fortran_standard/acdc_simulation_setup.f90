module acdc_simulation_setup
use acdc_system

implicit none

logical, parameter :: solve_ss = .true. 		! solve the steady state or run only for a given time
logical, parameter :: j_flux = .false.			! J calculated as a flux or as particles per time step; note that the former
                                                ! is automatically used for solve_ss = .true. regardless of the j_flux logical
logical, save :: use_solver						! use a solver or a simpler Euler method; determined in get_acdc_J based on the system size

integer, parameter :: nbins = 5					! number of size bins for size-classifying the clusters (if used)


contains

!--------------------------------------------------------------------------------------------------
! Parameters related to the size of the simulation system, including
! the numbers of clusters and equations, indices of outgoing fluxes, ...
!--------------------------------------------------------------------------------------------------

subroutine get_system_size(neq_syst,nclust_syst,nout_syst,diameter_max_syst)
!use restriction_criteria, only : diameter_max_loop

	implicit none
    
	integer :: nclust_syst, neq_syst			! numbers of clusters and equations
	integer :: nout_syst(n_charges)				! indices of outgoing fluxes
	real(kind(1.d0)) :: diameter_max_syst		! max. diameter in system
	integer, allocatable :: clust_comp(:,:)		! cluster composition in the loop mode
	
	if (small_set_mode) then
	
		neq_syst = neq
		nclust_syst = nclust
		
		nout_syst = nout_all
		
		diameter_max_syst = diameter_max*1.d-9
		
	else
	
		allocate(clust_comp(nclust_max,n_mol_types))
		if (n_mol_types .gt. 1) then
			call get_molecule_numbers(nclust_syst,clust_comp)
		else
			nclust_syst = nclust_max
		end if
		neq_syst = nclust_syst+(neq_max-nclust_max)
		
		nout_syst = nclust_syst+1				! loop mode: last, additional index
		
		!diameter_max_syst = diameter_max_loop
		
	end if
	
end subroutine get_system_size

!--------------------------------------------------------------------------------------------------
! Settings for the monomer concentrations
!--------------------------------------------------------------------------------------------------

subroutine sources_and_constants(neqn,source,isconst,fitted)

	implicit none
    
	real(kind(1.d0)) :: source(neqn)
	logical :: isconst(neqn)
	integer :: neqn, fitted(neqn,0:neqn), cluster_numbers(n_1A_clusters), n, i
    integer :: n_monomers(n_monomer_types)
	
	source = 0.d0								! Initialize all source terms to 0
	isconst = .false.							! Initialize all concentrations to vary freely according to the eqs.
	fitted = 0
    
    ! NOTE: If some concentrations are forced to be constant by excluding their equations (Perl option --no_eq),
    ! the below settings are NOT used


    ! Set some or all monomers to have constant or fitted concentrations
    
    ! If everything related to the constant or fitted concentrations is commented out,
    ! all concentrations will be explicitly determined from the equations
    
    if (solve_ss) then
    
        if (small_set_mode) then

            ! Use constant concentrations for electrically neutral monomers
            isconst(neutral_monomers) = .true.			! Neutral monomer concentrations are constant
                                                        ! (Ionic monomers can vary according to the ion production rate)

            ! ! Fit concentrations (here for 1A)
            ! isconst(n1A) = .false.
            ! fitted(1,0) = 1                                        ! 1 concentration is fitted
            ! call clusters_with_1_A(cluster_numbers)
            ! n = n_1A_clusters
            ! fitted(2,0:n) = (/n1A, n-1, cluster_numbers(2:n)/)		! Concentration of 1A is fitted using
                                                                    ! !  n-1 other concentrations: [1A]+... = const.
            ! ! Other monomer concentrations will remain as set above
        
        else
        
            ! Use constant concentrations for all monomers
            call monomer_indices(n_monomers)
            isconst(n_monomers) = .true.
        
        end if
    
    end if
	
end subroutine sources_and_constants

!--------------------------------------------------------------------------------------------------
! Indices of the vapor molecules in the cluster and molecule type arrays
!--------------------------------------------------------------------------------------------------

subroutine get_vapor_indices(names_vapor,n1vapor,nmol_vapor)

	implicit none
    
    character(len=11), dimension(n_neutral_monomers), intent(in) :: names_vapor
	integer, intent(out) :: n1vapor(n_neutral_monomers), nmol_vapor(n_neutral_monomers)
	character(len=11), dimension(n_monomer_types) :: names_monomers
    integer :: n_monomers(n_monomer_types)
    character(len=11), dimension(n_mol_types) :: names_mol_types
	integer :: i, j
    logical :: lfound
	
	! Find the cluster and molecule numbers for the included vapors
	
	n1vapor = 0
    nmol_vapor = 0
	
	call monomer_names(names_monomers)
	call monomer_indices(n_monomers)
    call molecule_names(names_mol_types)
	
	do i = 1,size(names_vapor)
    
        lfound = .false.
		do j = 1,size(names_monomers)
			if (trim('1'//names_vapor(i)(:)) .eq. trim(names_monomers(j)(:))) then
				n1vapor(i) = n_monomers(j)
                lfound = .true.
                !write(*,*) names_vapor(i)(:), ': cluster no. ', n1vapor(i)
				exit
			end if
		end do
        if (.not. lfound) then
            write(*,*) 'Monomer ', names_vapor(i)(:), ' not included in the ACDC set-up'
            stop
        end if
        
        lfound = .false.
		do j = 1,size(names_mol_types)
			if (trim(names_vapor(i)(:)) .eq. trim(names_mol_types(j)(:))) then
                nmol_vapor(i) = j
                lfound = .true.
                !write(*,*) names_vapor(i)(:), ': molecule no. ', nmol_vapor(i)
				exit
			end if
		end do
        if (.not. lfound) then
            write(*,*) 'Molecule ', names_vapor(i)(:), ' not included in the ACDC set-up'
            stop
        end if
        
	end do
	
end subroutine get_vapor_indices

!--------------------------------------------------------------------------------------------------
! Limits of the size bins for size-classifying the output cluster concentrations
!--------------------------------------------------------------------------------------------------

subroutine get_bin_limits(bin_limits)

	implicit none
    
	real(kind(1.d0)) :: bin_limits(nbins+1)     ! vector of the bin limits (m)

	! Bin limits according to mobility diameter
	bin_limits = (/1.05d0, 1.28d0, 1.73d0, 2.59d0, 4.27d0, 6.36d0/)*1.d-9 ! m

end subroutine get_bin_limits

!--------------------------------------------------------------------------------------------------
! Parameters related to size-classification, including the indices of clusters for each size bin,
! numbers of clusters in bins, ...
!--------------------------------------------------------------------------------------------------

subroutine group_size_bins(grouped_indices,nclust_per_bin)

	implicit none
    
	integer :: grouped_indices(0:nbins,nclust_max)			! indices of the clusters belonging to each bin
	integer :: nclust_per_bin(0:nbins)						! total number of clusters belonging to each bin
	real(kind(1.d0)) :: mi, ri, ivalue
	real(kind(1.d0)) :: bin_limits(nbins+1)					! vector of the bin limits (m)
	integer :: i, j, nclust_syst, indices(nclust_max,n_mol_types)
	logical :: lfound

	grouped_indices = 0
	nclust_per_bin =  0
	lfound = .false.
	
	if (.not. small_set_mode) then              ! The conditional is needed to compile also in the small set mode

		if(n_mol_types .gt. 1) then
			call get_molecule_numbers(nclust_syst,indices)
		else
			nclust_syst = nclust_max
			indices = reshape((/((i,i=1,nclust_syst),j=1,n_mol_types)/), shape(indices))
		end if

		call get_bin_limits(bin_limits)
		do i=1,nclust_syst
			if (sum(indices(i,:)) .eq. 1) then
				cycle                           ! don't add monomers to the size bins
			end if
			call get_masses_and_radii(mi,ri,indices(i,:))
			ivalue = 2.d0*ri + 0.3d0*1.d-9      ! assuming d_mob = d_mass + 0.3 nm
			if (ivalue .ge. bin_limits(1)) then
				lfound = .false.
				do j=1,nbins
					if ((ivalue .ge. bin_limits(j)) .and. (ivalue .lt. bin_limits(j+1))) then
						nclust_per_bin(j) = nclust_per_bin(j) + 1
						grouped_indices(j,nclust_per_bin(j)) = i
						lfound = .true.
						exit
					end if
				end do
				if (.not. lfound) then
					write(*,*) 'Could not group cluster ', i, ', molecules: ', indices(i,:)
					stop
				end if
			else
				nclust_per_bin(0) = nclust_per_bin(0) + 1
				grouped_indices(0,nclust_per_bin(0)) = i
				!write(*,*) 'Cluster below the lowest bin limit: ' , i, ', molecules: ', indices(i,:)
			end if
		end do
	
	end if

end subroutine group_size_bins


end module acdc_simulation_setup
