        MODULE deriv
        USE constants
        USE initialise,         ONLY: atoms
        USE hausdorff, &
     &    ONLY: compute_haus_distances,calc_haus_dist, &
     &          haus_rotate, haus_translate 

        REAL(KIND=DBL),PARAMETER:: disp_stride=0.0001D0  !translational displacement stride
        REAL(KIND=DBL),PARAMETER:: angle_stride=PI/(180.0D0*1000.0D0)  !rotational angle stride
        REAL(KIND=DBL)          :: x_deriv,y_deriv,z_deriv
        REAL(KIND=DBL)          :: phi_deriv,theta_deriv,omega_deriv


        CONTAINS
!===================================================================================
! List of subroutines:
!       haus_deriv
!       finite_diff
!       haus_increment
!===================================================================================
        SUBROUTINE haus_deriv(param,input_x,gradient)
! first derivates for displacements in x,y and z-axes and rotations
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(6),INTENT(IN)      :: param
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: input_x   ! coordinates
        REAL(KIND=DBL),DIMENSION(6),INTENT(OUT)     :: gradient   ! 1st  derivatives or gradients
        REAL(KIND=DBL),DIMENSION(3)                 :: translate_param
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE   :: coord
        INTEGER                 :: iter

        IF(ALLOCATED(coord)) DEALLOCATE(coord)
        ALLOCATE(coord(atoms,3))

        coord = input_x


! rotate along x-axis, perform finite difference calculation
        CALL haus_rotate(1,param(1),coord)
        CALL finite_diff(coord,1,phi_deriv)
! rotate along x-axis then y-axis, perform finite difference calculation
        CALL haus_rotate(2,param(2),coord)
        CALL finite_diff(coord,2,theta_deriv)
! rotate along x-axis then y-axis, then z-axis, perform finite difference calculation
        CALL haus_rotate(3,param(3),coord)
        CALL finite_diff(coord,3,omega_deriv)

! translate in 3D space, perform finite difference calculation
        DO iter=1,3
          translate_param(iter) = param(iter+3)
        END DO
        CALL haus_translate(translate_param,coord)
        CALL finite_diff(coord,4,x_deriv)
        CALL finite_diff(coord,5,y_deriv)
        CALL finite_diff(coord,6,z_deriv)
        
!        PRINT *, 'in haus_deriv'
!        PRINT *, x_deriv,y_deriv,z_deriv
!        PRINT *, phi_deriv,theta_deriv,omega_deriv

        gradient(1) = phi_deriv
        gradient(2) = theta_deriv
        gradient(3) = omega_deriv
        gradient(4) = x_deriv
        gradient(5) = y_deriv
        gradient(6) = z_deriv
        
        END SUBROUTINE haus_deriv
        

        SUBROUTINE haus_deriv_xyz(param,input_x,gradient)
! first derivates for displacements in x,y and z-axes 
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(IN)      :: param
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: input_x   ! coordinates
        REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)     :: gradient   ! 1st  derivatives or gradients
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE   :: coord
        INTEGER                 :: iter

        IF(ALLOCATED(coord)) DEALLOCATE(coord)
        ALLOCATE(coord(atoms,3))

        coord = input_x

! translate in 3D space, perform finite difference calculation
        CALL haus_translate(param,coord)
        CALL finite_diff(coord,4,x_deriv)
        CALL finite_diff(coord,5,y_deriv)
        CALL finite_diff(coord,6,z_deriv)
        
!        PRINT *, 'in haus_deriv'
!        PRINT *, x_deriv,y_deriv,z_deriv
!        PRINT *, phi_deriv,theta_deriv,omega_deriv

        gradient(1) = x_deriv
        gradient(2) = y_deriv
        gradient(3) = z_deriv
        
        END SUBROUTINE haus_deriv_xyz
        

        SUBROUTINE haus_deriv_angle(param,param_index,input_x,gradient)
! first derivates for a rotation along one chosen axis
        IMPLICIT NONE
        REAL(KIND=DBL),INTENT(IN)      :: param
        INTEGER,INTENT(IN)             :: param_index
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: input_x   ! coordinates
        REAL(KIND=DBL),INTENT(OUT)     :: gradient   ! 1st  derivatives or gradients
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE   :: coord
        INTEGER                 :: iter

        IF(ALLOCATED(coord)) DEALLOCATE(coord)
        ALLOCATE(coord(atoms,3))

        coord = input_x


        SELECT CASE(param_index)
! rotate along x-axis, perform finite difference calculation
        CASE(1)
          CALL haus_rotate(1,param,coord)
          CALL finite_diff(coord,1,phi_deriv)
          gradient = phi_deriv
! rotate along y-axis, perform finite difference calculation
        CASE(2)
          CALL haus_rotate(2,param,coord)
          CALL finite_diff(coord,2,theta_deriv)
          gradient = theta_deriv
! rotate along z-axis, perform finite difference calculation
        CASE(3)
          CALL haus_rotate(3,param,coord)
          CALL finite_diff(coord,3,omega_deriv)
          gradient = omega_deriv
        END SELECT
        
        END SUBROUTINE haus_deriv_angle



        SUBROUTINE finite_diff(coord,var_index,deriv)
! calculate first derivate using the symmeterized form finite difference
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN):: coord
        INTEGER       ,INTENT(IN)               :: var_index
        REAL(KIND=DBL),INTENT(OUT)              :: deriv
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE       :: temp_1
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE       :: temp_2
        REAL(KIND=DBL)  :: forward, backward,delta_h
        INTEGER         :: iter

        IF(ALLOCATED(temp_1)) DEALLOCATE(temp_1)
        IF(ALLOCATED(temp_2)) DEALLOCATE(temp_2)
        ALLOCATE(temp_1(atoms,3))
        ALLOCATE(temp_2(atoms,3))

        SELECT CASE (var_index)
        CASE(1:3)
          delta_h = angle_stride
        CASE(4:6)
          delta_h = disp_stride
        END SELECT

        CALL haus_increment(var_index,'forward',coord,temp_1)
        CALL haus_increment(var_index,'backward',coord,temp_2)
        CALL compute_haus_distances(temp_1)
        CALL calc_haus_dist(forward)
        CALL compute_haus_distances(temp_2)
        CALL calc_haus_dist(backward)

        deriv = forward - backward
        deriv = deriv / (2.0D0*delta_h)
        
        
!        PRINT *, 'diff, var=', var_index, 'finite diff=',deriv


!        PRINT *,''
!        PRINT *,'temp_1'
!        DO iter=1,atoms
!          PRINT *, temp_1(iter,1),temp_1(iter,2),temp_1(iter,3)
!        END DO
!        PRINT *,''
!        PRINT *,'temp_2'
!        DO iter=1,atoms
!          PRINT *, temp_2(iter,1),temp_2(iter,2),temp_2(iter,3)
!        END DO
!        PRINT *,''
        END SUBROUTINE finite_diff

        SUBROUTINE haus_increment(var,direction,in_coord,out_coord)
! Forward or backward increment of the coordinates, either by rotation
! or translation, for finite difference calculation
        IMPLICIT NONE
        INTEGER,INTENT(IN)  :: var  ! increment variable
                                    ! 1,2,3: phi,theta,omega increments
                                    ! 4,5,6: x,y,z increments
        CHARACTER(LEN=*),INTENT(IN)         :: direction ! Forward or backward increment
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)             ::in_coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)::out_coord
        LOGICAL                  :: isforward    ! True for forward increment, false for backward decrement
        REAL(KIND=DBL)                                :: phi
        REAL(KIND=DBL)                                :: trans_stride
        REAL(KIND=DBL)                                :: rotat_stride
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE     :: vec
        REAL(KIND=DBL),DIMENSION(3,3)                 :: rotate
        INTEGER                                       :: iter

        IF (TRIM(direction) == 'forward') THEN
          isforward = .TRUE.
        ELSE IF(TRIM(direction) =='backward') THEN
          isforward = .FALSE.
        ELSE 
          PRINT *, 'increment direction is messed up'
        END IF

        IF(ALLOCATED(out_coord)) DEALLOCATE(out_coord)
        ALLOCATE(out_coord(atoms,3))
        out_coord = in_coord


! logical test for forward or backward
        IF(isforward) THEN
          trans_stride = disp_stride     ! forward
          rotat_stride = angle_stride    ! anti-clockwise rotation
        ELSE
          trans_stride = -disp_stride    ! backward
          rotat_stride = -angle_stride   ! clockwise rotation
        END IF


! rotation operation
        phi = rotat_stride
        IF (var > 0 .AND. var < 4) THEN
          IF (ALLOCATED(vec)) DEALLOCATE(vec)
          ALLOCATE(vec(3,atoms))
          vec = TRANSPOSE(in_coord)
          SELECT CASE(var)
          CASE(1)    ! rotation about x-axis, phi
            rotate(1,1) = 1.0D0
            rotate(1,2) = 0.0D0
            rotate(1,3) = 0.0D0
            rotate(2,1) = 0.0D0
            rotate(2,2) = COS(phi) 
            rotate(2,3) = -SIN(phi)
            rotate(3,1) = 0.0D0
            rotate(3,2) = SIN(phi) 
            rotate(3,3) = COS(phi)
          CASE(2)   ! rotation about y-axis,theta
            rotate(1,1) = COS(phi)
            rotate(1,2) = 0.0D0
            rotate(1,3) = SIN(phi)
            rotate(2,1) = 0.0D0
            rotate(2,2) = 1.0D0
            rotate(2,3) = 0.0D0
            rotate(3,1) = -SIN(phi)
            rotate(3,2) = 0.0D0
            rotate(3,3) = COS(phi)
           CASE(3)  ! rotations about z-axis,omega
            rotate(1,1) = COS(phi)
            rotate(1,2) = -SIN(phi)
            rotate(1,3) = 0.0D0
            rotate(2,1) = SIN(phi)
            rotate(2,2) = COS(phi)
            rotate(2,3) = 0.0D0
            rotate(3,1) = 0.0D0
            rotate(3,2) = 0.0D0
            rotate(3,3) = 1.0D0
          END SELECT
          vec = MATMUL(rotate,vec)
          out_coord = TRANSPOSE(vec)
        END IF

! translation operation
        SELECT CASE(var)
        CASE(4)   ! displacement along x-axis
          out_coord(:,1) = in_coord(:,1) + trans_stride
        CASE(5)   ! displacement along y-axis
          out_coord(:,2) = in_coord(:,2) + trans_stride
        CASE(6)   ! displacement along z-axis
          out_coord(:,3) = in_coord(:,3) + trans_stride
        END SELECT

        END SUBROUTINE haus_increment


        END MODULE deriv
