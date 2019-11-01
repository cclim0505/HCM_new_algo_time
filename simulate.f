        MODULE simulate
        USE initialise          
        USE hausdorff
        USE deriv
        USE optimization
        USE rotations
        CONTAINS

        SUBROUTINE single_simul
! start single computation session
        IMPLICIT NONE
        INTEGER         :: iter,jter
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: temp1_coord  
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: temp2_coord  
        REAL(KIND=DBL)                             :: haus_dist    ! chirality measure


!========================================================================
! Variables used sorting the lowest 1 percent value
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE    :: haus_array
        INTEGER,DIMENSION(:),ALLOCATABLE    :: loca_array
        INTEGER                             :: loca_array_size
        INTEGER,DIMENSION(1)                       :: loca
!========================================================================
! Variables used for two-stage optimization
        REAL(KIND=DBL),DIMENSION(3)     :: min_angle
!========================================================================

        CALL calc_haus_diameter
        CALL compute_haus_distances(mir_coord)
        CALL calc_haus_dist(haus_dist)


!========================================================================
! Test multiple rotated structures and perform LBFGS
!========================================================================
        CALL gen_possible_rotations(mir_coord)  

!       DO iter=1,struct_count
!         temp1_coord = rotated_struct(iter,:,:)
!         CALL optim_lbfgs(temp1_coord)
!       END DO

!========================================================================
! Calculate haus distances for all the possible rotations without
! minimization & saves all haus value into array for sorting
        IF (ALLOCATED(temp1_coord)) DEALLOCATE(temp1_coord)
        IF (ALLOCATED(temp2_coord)) DEALLOCATE(temp2_coord)
        ALLOCATE(temp1_coord(atoms,3))
        ALLOCATE(temp2_coord(atoms,3))
        IF (ALLOCATED(haus_array)) DEALLOCATE(haus_array)
        ALLOCATE (haus_array(struct_count))

        DO iter=1,struct_count
          temp1_coord = rotated_struct(iter,:,:)
          CALL compute_haus_distances(temp1_coord)
          CALL calc_haus_dist(haus_dist)
          haus_array(iter) = haus_dist
        END DO

        PRINT *, haus_array

!========================================================================
! Find lowest percentile of haus values
        loca_array_size = REAL(struct_count) * opt_percentile
        IF (loca_array_size < 1) loca_array_size = 1
        IF (loca_array_size > struct_count) loca_array_size=struct_count
        ALLOCATE (loca_array(loca_array_size))

        CALL get_lowest_percentile_haus(haus_array,loca_array)

        PRINT *, 'loca_array after subroutine'
        PRINT *, loca_array
!========================================================================
! Perform L-BFGS for only lowest percentile values
! Suggested algorithm
       DO iter=1,loca_array_size
         temp1_coord = rotated_struct(loca_array(iter),:,:)
         min_angle = 0.0D0
         DO jter=1,3
           CALL optim_lbfgs_angle(jter,temp1_coord,min_angle(jter))
           temp2_coord = temp1_coord
           CALL haus_rotate(jter,min_angle(jter),temp2_coord)
           CALL optim_lbfgs_xyz(temp2_coord)
         END DO
       END DO
!========================================================================
! Perform L-BFGS for only lowest percentile values
! testing some algorithm for angle deriv
!      DO iter=1,loca_array_size
!        temp1_coord = rotated_struct(loca_array(iter),:,:)
!        min_angle = 0.0D0
!        DO jter=1,3
!          CALL optim_lbfgs_angle(jter,temp1_coord,min_angle(jter))
!        END DO
!      END DO
!========================================================================
! Perform L-BFGS for only lowest percentile values
! testing some algorithm for xyz deriv
!      DO iter=1,loca_array_size
!        temp1_coord = rotated_struct(loca_array(iter),:,:)
!        CALL optim_lbfgs_xyz(temp1_coord)
!      END DO
!========================================================================
! L-BFGS for angles
!       DO iter=1,1
!         DO jter=1,struct_count
!           temp1_coord = rotated_struct(jter,:,:)
!           CALL optim_lbfgs_angle(iter,temp1_coord)
!         END DO
!       END DO
!========================================================================
! L-BFGS for displacements
!        DO iter=1,struct_count
!          temp1_coord = rotated_struct(iter,:,:)
!          CALL optim_lbfgs_xyz(temp1_coord)
!        END DO
!========================================================================

        CALL print_hcm

        END SUBROUTINE single_simul



        END MODULE simulate

!        PRINT *, ''
!        PRINT *, 'mir coord'
!        PRINT *, ''
!        DO iter=1,atoms
!          PRINT *, mir_coord(iter,1),mir_coord(iter,2),mir_coord(iter,3)
!        END DO
!        PRINT *, ''
