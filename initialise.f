        MODULE initialise
        USE     constants               ,ONLY: DBL
!=========================================================================================
! VARIABLE DICTIONARY
!=========================================================================================
! Simulation parameters
        INTEGER         :: atoms
        REAL(KIND=DBL),DIMENSION(:,:), ALLOCATABLE      :: inp_coord   ! input structure coordinates, set Q
        REAL(KIND=DBL),DIMENSION(:,:), ALLOCATABLE      :: mir_coord   ! mirror structure coordinates, set Q_prime
        REAL(KIND=DBL),DIMENSION(:,:), ALLOCATABLE      :: mirror_2    ! to store temporary coordinates
        CHARACTER(20)       :: session_file='session_in.dat'
        CHARACTER(30)       :: input_struc_file
        CHARACTER(30)       :: mirror_struc_file
        CHARACTER(20)       :: output_file
        INTEGER             :: angle_interval
        LOGICAL             :: ismirror
        REAL                :: opt_percentile
!=========================================================================================
        CONTAINS
! List of subroutines:
! start_session
!       :start a new session, setting up structure coordinates
! read_session
!       :read input parameters from session file
! read_two_coord
!       :read in input and mirror image coordinates
! calc_centroid
!       :calculate centroid for set_coord_to_origin
! set_coord_to_origin
!       :set respective centroids of input and mirrorc coordinates to origin
! print_two coord
!       :output coordinates for checking
! mirror_operation
!       :produce mirror structure

        SUBROUTINE start_session
        IMPLICIT NONE
        CALL read_session
        CALL read_two_coord
        CALL set_coord_to_origin


        IF (ismirror) THEN
          IF (ALLOCATED(mirror_2)) DEALLOCATE(mirror_2)
          ALLOCATE(mirror_2(atoms,3))
          CALL mirror_operation(inp_coord,mirror_2)
          mir_coord = mirror_2
        END IF

        CALL print_two_coord

        END SUBROUTINE start_session

        SUBROUTINE read_session
! read session input file
        IMPLICIT NONE
        CHARACTER       :: dummy
        OPEN(20,file=session_file,status='old')
        READ(20,*) dummy, input_struc_file
        READ(20,*) dummy, ismirror
        READ(20,*) dummy, mirror_struc_file
        READ(20,*) dummy, output_file
        READ(20,*) dummy, angle_interval
        READ(20,*) dummy, opt_percentile
        CLOSE(20)
        END SUBROUTINE read_session
        
        SUBROUTINE read_two_coord
! read in input and mirror image coordinates
        IMPLICIT NONE
        INTEGER         :: iter
        CHARACTER       :: dummy
        OPEN(21,file=input_struc_file,status='old')
        OPEN(22,file=mirror_struc_file,status='old')
        READ(21,*) atoms
        READ(21,*) dummy
        READ(22,*) dummy
        READ(22,*) dummy
        IF (ALLOCATED(inp_coord)) DEALLOCATE (inp_coord)
        IF (ALLOCATED(mir_coord)) DEALLOCATE (mir_coord)
        ALLOCATE(inp_coord(atoms,3))
        ALLOCATE(mir_coord(atoms,3))
        DO iter=1,atoms
          READ(21,*) dummy,inp_coord(iter,1), &
     &       inp_coord(iter,2),inp_coord(iter,3)
          READ(22,*) dummy,mir_coord(iter,1), &
     &       mir_coord(iter,2),mir_coord(iter,3)
        END DO
        CLOSE(21)
        CLOSE(22)
        END SUBROUTINE read_two_coord

        SUBROUTINE print_two_coord
! output coordinates for checking
        IMPLICIT NONE
        INTEGER :: iter
        PRINT *, ''
        PRINT *, 'Input Coordinates in xyz'
        PRINT *, ''
        DO iter=1,atoms
          PRINT '(I6,2X,3(F13.6,2X))', iter, &
     &inp_coord(iter,1),inp_coord(iter,2),inp_coord(iter,3)
        END DO
        PRINT *, ''
        PRINT *, 'Mirror Coordinates in xyz'
        PRINT *, ''
        DO iter=1,atoms
          PRINT '(I6,2x,3(F13.6,2X))', iter, &
     &mir_coord(iter,1),mir_coord(iter,2),&
     &mir_coord(iter,3)
        END DO

        END SUBROUTINE print_two_coord

        SUBROUTINE calc_centroid(coord,centroid)
! calculate centroid for set_coord_to_origin
        IMPLICIT NONE
        INTEGER :: iter
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: coord
        REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)     :: centroid
        centroid (:) = 0.0
!        PRINT *, ''
!        PRINT *, 'Coord in calc_centroid'
!        PRINT *, ''
!        DO iter=1,atoms
!          PRINT *, iter,coord(iter,1), coord(iter,2), coord(iter,3)
!        END DO
        DO iter=1,atoms
          centroid(1) = centroid(1) + coord(iter,1)
          centroid(2) = centroid(2) + coord(iter,2)
          centroid(3) = centroid(3) + coord(iter,3)
        END DO

        centroid(:) = centroid(:) / REAL(atoms)

!        PRINT *, ''
!        PRINT *, 'centroid =', centroid(1), centroid(2), centroid(3)
!        PRINT *, ''
        END SUBROUTINE calc_centroid
 
        SUBROUTINE set_coord_to_origin
! set coordinates' centroid at origin
        IMPLICIT NONE
        INTEGER :: iter
        REAL(KIND=DBL),DIMENSION(3)     :: centroid
        CALL calc_centroid(inp_coord,centroid)
        DO iter=1,atoms
          inp_coord(iter,:) = inp_coord (iter,:) - centroid(:)
        END DO
        CALL calc_centroid(mir_coord,centroid)
        DO iter=1,atoms
          mir_coord(iter,:) = mir_coord (iter,:) - centroid(:)
        END DO
        END SUBROUTINE set_coord_to_origin

        SUBROUTINE mirror_operation(in_x,out_x)
! produce mirror structure
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)        :: in_x
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)::out_x
        IF (ALLOCATED(out_x)) DEALLOCATE (out_x)
        ALLOCATE(out_x(atoms,3))
        out_x = in_x
        out_x(:,1) = out_x(:,1) * (-1.0D0)
        END SUBROUTINE mirror_operation
       
        END MODULE initialise
