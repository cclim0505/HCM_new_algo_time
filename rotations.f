        MODULE rotations
        USE constants
        USE initialise,         ONLY:atoms,mir_coord,angle_interval
        USE hausdorff
        REAL(KIND=DBL),DIMENSION(:,:,:),ALLOCATABLE :: rotated_struct
        INTEGER         :: struct_count   ! the number of generated structures
        SAVE :: rotated_struct
        CONTAINS

        SUBROUTINE gen_possible_rotations(in_coord)
! generate possible rotations for an input structure
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)   :: in_coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: temp1_x
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: temp2_x
        REAL(KIND=DBL)  :: phi
        INTEGER                         :: iter,jter,kter,lter
        INTEGER                         :: counter=0

        IF (ALLOCATED(rotated_struct)) DEALLOCATE(rotated_struct)
        IF (ALLOCATED(temp1_x)) DEALLOCATE(temp1_x)
        IF (ALLOCATED(temp2_x)) DEALLOCATE(temp2_x)
        ALLOCATE(rotated_struct(angle_interval**3,atoms,3))
        ALLOCATE(temp1_x(atoms,3))
        ALLOCATE(temp2_x(3,atoms))

        temp1_x = in_coord
        phi = 2.0D0 * PI / DBLE(angle_interval)


        DO iter=1,angle_interval
          CALL haus_rotate(1,phi,temp1_x)
          DO jter=1,angle_interval
            CALL haus_rotate(2,phi,temp1_x)
            DO kter=1,angle_interval
              CALL haus_rotate(3,phi,temp1_x)
              counter = counter + 1
              rotated_struct(counter,:,:) = temp1_x
            END DO
          END DO
        END DO

        struct_count = counter
        CALL print_rotated_struct(counter)


        END SUBROUTINE gen_possible_rotations

        SUBROUTINE print_rotated_struct(no_of_structures)
        IMPLICIT NONE
        INTEGER         :: iter, jter
        INTEGER         :: no_of_structures
        OPEN(32,file='rotated_struct.dat',status='replace')
        DO iter=1,no_of_structures
          WRITE(32,*) iter
          DO jter=1,atoms
            WRITE(32,*)  rotated_struct(iter,jter,1), &
     &         rotated_struct(iter,jter,2),rotated_struct(iter,jter,3)
          END DO
        END DO
        CLOSE(32)
        END SUBROUTINE print_rotated_struct


        SUBROUTINE get_lowest_percentile_haus(haus_array,loca_array)
! Get the index of lowest percentile's haus values
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:),INTENT(IN)  :: haus_array
        INTEGER,DIMENSION(:),INTENT(OUT)   :: loca_array
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: temp_array

        INTEGER         :: temp_size
        INTEGER         :: iter
        INTEGER         :: loca_array_size
        INTEGER,DIMENSION(1)    :: loca

        temp_size = SIZE(haus_array)
        IF (ALLOCATED(temp_array)) DEALLOCATE(temp_array)
        ALLOCATE(temp_array(temp_size))
        temp_array = haus_array

        loca_array_size = SIZE(loca_array)
        DO iter=1,loca_array_size
          loca = MINLOC(temp_array)
          loca_array(iter) = loca(1)
          temp_array(loca) = HUGE(temp_array)
        END DO

        DEALLOCATE(temp_array)

        END SUBROUTINE get_lowest_percentile_haus

        END MODULE rotations
