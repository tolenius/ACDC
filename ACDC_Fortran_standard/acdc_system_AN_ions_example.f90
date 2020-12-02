module acdc_system

implicit none

integer, parameter :: nclust = 54						! number of clusters, molecules and ions
integer, parameter :: neq = 63							! number of equations
integer, parameter :: nclust_max = 54, neq_max = 63		! same as nclust and neq in the small set mode
logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: n_variables = 4
logical, parameter :: variable_temp = .true.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %

real(kind(1.d0)), parameter :: cs_exponent_default = -1.60000000000000d+00			! parameters for the exponential loss
real(kind(1.d0)), parameter :: cs_coefficient_default = 2.60000000000000d-03

integer, parameter :: n_monomer_types = 4
integer, parameter :: n_charges = 3			! number of charging states
integer, parameter :: n1A = 1, n1N = 3, n1B = 17, n1N1P = 35, nneg = 53, npos = 54			! cluster indices for monomers and ions

integer, parameter :: nout_all(3) = (/60, 61, 62/), nout_neu = 60, nout_neg = 61, nout_pos = 62			! indices for outgoing fluxes

integer, parameter :: n_mol_types = 4
integer, parameter :: nmolA = 1, nmolB = 2, nmolN = 3, nmolP = 4			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 2				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 16			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 19			! negative
integer, parameter :: n_positives = 19			! positive

integer, parameter :: n_neutral_monomers = 2, neutral_monomers(2) = (/1, 3/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 1, negative_monomers(1) = (/17/)			! negative
integer, parameter :: n_positive_monomers = 1, positive_monomers(1) = (/35/)			! positive

integer, parameter :: n_neutral_clusters = 14			! number of neutral clusters
integer, parameter :: n_negative_clusters = 17			! negative
integer, parameter :: n_positive_clusters = 17			! positive

real(kind(1.d0)), parameter :: mass_max = 576.60
real(kind(1.d0)), parameter :: diameter_max = 1.07
real(kind(1.d0)), parameter :: mob_diameter_max = 1.37
integer, parameter :: ij_ind_max(4) = (/5, 1, 5, 1/)		! maximum molecular content
integer, parameter :: n_bound = 12		! number of clusters at the system boundary

integer, parameter :: nmols_out_neutral(1, 4) = reshape((/6, 0, 5, 0/),(/1, 4/))			! criteria for outgrowing neutrals
integer, parameter :: nmols_out_negative(1, 4) = reshape((/5, 1, 3, 0/),(/1, 4/))			! negatives
integer, parameter :: nmols_out_positive(1, 4) = reshape((/5, 0, 6, 1/),(/1, 4/))			! positives


contains

subroutine n_A_in_clusters(n_A)
	implicit none
	integer :: n_A(54)

	n_A = (/1, 2, 0, 1, 2, 3, 2, 3, 4, 3, &
		&4, 5, 4, 5, 4, 5, 1, 2, 3, 4, &
		&5, 1, 2, 3, 4, 5, 1, 2, 3, 4, &
		&5, 4, 5, 5, 0, 1, 2, 0, 1, 2, &
		&3, 0, 1, 2, 3, 4, 2, 3, 4, 3, &
		&4, 5, 0, 0/)

end subroutine n_A_in_clusters

subroutine clusters_with_1_A(cluster_numbers)
	implicit none
	integer :: cluster_numbers(2)

	cluster_numbers = (/1, 4/)

end subroutine clusters_with_1_A

subroutine arrays(neutrals, negatives, positives)
	implicit none
	integer :: neutrals(16), negatives(19), positives(19)

	neutrals = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)
	negatives = (/17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, &
		&32, 33, 34, 53/)
	positives = (/35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, &
		&50, 51, 52, 54/)

end subroutine arrays

subroutine cluster_arrays(neutral_clusters, negative_clusters, positive_clusters)
	implicit none
	integer :: neutral_clusters(14), negative_clusters(17), positive_clusters(17)

	neutral_clusters = (/2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)
	negative_clusters = (/18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, &
		&32, 33, 34/)
	positive_clusters = (/36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, &
		&50, 51, 52/)

end subroutine cluster_arrays

subroutine get_charging_state(charging_state)
	implicit none
	integer :: charging_state(54)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 3, 3, 3, 3, 3, 3, &
		&3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
		&3, 3, 2, 3/)

end subroutine get_charging_state

subroutine get_mass(mass)
	implicit none
	real(kind(1.d0)) :: mass(54)

	mass = (/98.08, 196.16, 17.04, 115.12, 213.20, 311.28, 230.24, 328.32, 426.40, 345.36, &
		&443.44, 541.52, 460.48, 558.56, 477.52, 575.60, 97.08, 195.16, 293.24, 391.32, &
		&489.40, 114.12, 212.20, 310.28, 408.36, 506.44, 131.16, 229.24, 327.32, 425.40, &
		&523.48, 442.44, 540.52, 557.56, 18.04, 116.12, 214.20, 35.08, 133.16, 231.24, &
		&329.32, 52.12, 150.20, 248.28, 346.36, 444.44, 265.32, 363.40, 461.48, 380.44, &
		&478.52, 576.60, 32.00, 19.02/)

end subroutine get_mass

subroutine get_diameter(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(54)

	 diameter = (/0.55, 0.70, 0.43, 0.63, 0.75, 0.84, 0.79, 0.87, 0.94, 0.91, &
		&0.97, 1.03, 1.00, 1.05, 1.02, 1.07, 0.55, 0.70, 0.80, 0.88, &
		&0.95, 0.63, 0.75, 0.84, 0.91, 0.97, 0.69, 0.79, 0.87, 0.94, &
		&1.00, 0.97, 1.03, 1.05, 0.43, 0.63, 0.75, 0.54, 0.69, 0.79, &
		&0.87, 0.62, 0.74, 0.83, 0.91, 0.97, 0.87, 0.94, 1.00, 0.96, &
		&1.02, 1.07, 0.45, 0.39/)	! dry value

end subroutine get_diameter

subroutine get_mob_diameter(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(54)

	 mob_diameter = (/0.85, 1.00, 0.73, 0.93, 1.05, 1.14, 1.09, 1.17, 1.24, 1.21, &
		&1.27, 1.33, 1.30, 1.35, 1.32, 1.37, 0.85, 1.00, 1.10, 1.18, &
		&1.25, 0.93, 1.05, 1.14, 1.21, 1.27, 0.99, 1.09, 1.17, 1.24, &
		&1.30, 1.27, 1.33, 1.35, 0.73, 0.93, 1.05, 0.84, 0.99, 1.09, &
		&1.17, 0.92, 1.04, 1.13, 1.21, 1.27, 1.17, 1.24, 1.30, 1.26, &
		&1.32, 1.37, 0.75, 0.69/)	! dry value

end subroutine get_mob_diameter

subroutine cluster_names(clust)
	implicit none
	character(len=11), dimension(63) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '1N'
	clust(4)(:) = '1A1N'
	clust(5)(:) = '2A1N'
	clust(6)(:) = '3A1N'
	clust(7)(:) = '2A2N'
	clust(8)(:) = '3A2N'
	clust(9)(:) = '4A2N'
	clust(10)(:) = '3A3N'
	clust(11)(:) = '4A3N'
	clust(12)(:) = '5A3N'
	clust(13)(:) = '4A4N'
	clust(14)(:) = '5A4N'
	clust(15)(:) = '4A5N'
	clust(16)(:) = '5A5N'
	clust(17)(:) = '1B'
	clust(18)(:) = '1A1B'
	clust(19)(:) = '2A1B'
	clust(20)(:) = '3A1B'
	clust(21)(:) = '4A1B'
	clust(22)(:) = '1B1N'
	clust(23)(:) = '1A1B1N'
	clust(24)(:) = '2A1B1N'
	clust(25)(:) = '3A1B1N'
	clust(26)(:) = '4A1B1N'
	clust(27)(:) = '1B2N'
	clust(28)(:) = '1A1B2N'
	clust(29)(:) = '2A1B2N'
	clust(30)(:) = '3A1B2N'
	clust(31)(:) = '4A1B2N'
	clust(32)(:) = '3A1B3N'
	clust(33)(:) = '4A1B3N'
	clust(34)(:) = '4A1B4N'
	clust(35)(:) = '1N1P'
	clust(36)(:) = '1A1N1P'
	clust(37)(:) = '2A1N1P'
	clust(38)(:) = '2N1P'
	clust(39)(:) = '1A2N1P'
	clust(40)(:) = '2A2N1P'
	clust(41)(:) = '3A2N1P'
	clust(42)(:) = '3N1P'
	clust(43)(:) = '1A3N1P'
	clust(44)(:) = '2A3N1P'
	clust(45)(:) = '3A3N1P'
	clust(46)(:) = '4A3N1P'
	clust(47)(:) = '2A4N1P'
	clust(48)(:) = '3A4N1P'
	clust(49)(:) = '4A4N1P'
	clust(50)(:) = '3A5N1P'
	clust(51)(:) = '4A5N1P'
	clust(52)(:) = '5A5N1P'
	clust(53)(:) = 'neg'
	clust(54)(:) = 'pos'
	clust(55)(:) = 'source'
	clust(56)(:) = 'coag'
	clust(57)(:) = 'wall'
	clust(58)(:) = 'dil'
	clust(59)(:) = 'rec'
	clust(60)(:) = 'out_neu'
	clust(61)(:) = 'out_neg'
	clust(62)(:) = 'out_pos'
	clust(63)(:) = 'bound'

end subroutine cluster_names

subroutine monomer_names(clust_mon)
	implicit none
	character(len=11), dimension(4) :: clust_mon

	clust_mon(1)(:) = '1A'
	clust_mon(2)(:) = '1N'
	clust_mon(3)(:) = '1B'
	clust_mon(4)(:) = '1N1P'

end subroutine monomer_names

subroutine molecule_names(labels)
	implicit none
	character(len=11), dimension(4) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'B'
	labels(3)(:) = 'N'
	labels(4)(:) = 'P'

end subroutine molecule_names

subroutine monomer_indices(n_monomers)
	implicit none
	integer :: n_monomers(4)

	n_monomers = (/1, 3, 17, 35/)

end subroutine monomer_indices

subroutine get_bound(bound_clusters,nmols_bound)
	implicit none
	integer :: bound_clusters(12), nmols_bound(12,4)

	nmols_bound(1,:) = (/5, 0, 3, 0/)
	nmols_bound(2,:) = (/5, 0, 4, 0/)
	nmols_bound(3,:) = (/4, 0, 5, 0/)
	nmols_bound(4,:) = (/5, 0, 5, 0/)
	nmols_bound(5,:) = (/4, 1, 0, 0/)
	nmols_bound(6,:) = (/4, 1, 1, 0/)
	nmols_bound(7,:) = (/4, 1, 2, 0/)
	nmols_bound(8,:) = (/4, 1, 3, 0/)
	nmols_bound(9,:) = (/4, 1, 4, 0/)
	nmols_bound(10,:) = (/3, 0, 5, 1/)
	nmols_bound(11,:) = (/4, 0, 5, 1/)
	nmols_bound(12,:) = (/5, 0, 5, 1/)

	bound_clusters = (/12, 14, 15, 16, 21, 26, 31, 33, 34, 50, 51, 52/)

end subroutine get_bound


end module acdc_system

