        MODULE hausdorff        
        USE constants
        USE initialise
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: dist_btw_pts  ! distance between points in sets Q and Q_prime
! dist_btw_pts (1,:) denote distance between first fixed point q and set Q_prime
! dist_btw_pts (2,:) denote distance between second fixed point q and set Q_prime
! dist_btw_pts (:,1) denote distance between first fixed point q_prime and set Q
! dist_btw_pts (:,2) denote distance between second fixed point q_prime and set Q

        REAL(KIND=DBL)               :: hcm=HUGE(atoms)  ! Hausdorff chirality measure, the final answer
        REAL(KIND=DBL)               :: diameter     ! Largest distance between two points in set Q

        CONTAINS 
!===============================================================================================
! List of subroutines:
! calc_haus_diameter
!       :calculate the diameter of set Q, used only once and for all
! compute_haus_distances
!       :compute distances between points in sets Q and Q_prime
! calc_haus_dist
!       :calculate the hausdorff distances between sets Q and Q_prime
! smallest_haus_distance
!       :find and save the smallest hausdorff distance or HCM after each calc_haus_dist
! print_hcm
!       :output HCM value
!===============================================================================================



        SUBROUTINE calc_haus_diameter
! calculate the diameter of set Q, i.e. the largest distance between any two points of Q
! calculation is required only to perform once in a simulation
! and repeatedly used for normalization
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: dist
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: diff
        INTEGER                                   :: iter,jter

        IF (ALLOCATED(diff)) DEALLOCATE(diff)
        IF (ALLOCATED(dist)) DEALLOCATE(dist)
        ALLOCATE(diff(atoms,3))
        ALLOCATE(dist(atoms,atoms))
        DO iter=1,atoms
          DO jter=1,atoms
            diff(jter,:) = inp_coord(iter,:) - inp_coord(jter,:)
            dist(iter,jter) = NORM2(diff(jter,:))
          END DO
        END DO
        diameter = MAXVAL(dist)
!        PRINT *, ''
!        PRINT *, 'diameter'
!        PRINT *, ''
!        DO iter=1,atoms
!          PRINT *, iter,dist(iter,1),dist(iter,2),dist(iter,3)
!        END DO
!        PRINT *, dist
!        PRINT *, 'largest=', diameter
        END SUBROUTINE calc_haus_diameter

        SUBROUTINE compute_haus_distances(coord_set)
! compute distances between points in sets Q and Q_prime
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: coord_set
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: diff
        INTEGER                                   :: iter,jter

        IF (ALLOCATED(diff)) DEALLOCATE(diff)
        IF (ALLOCATED(dist_btw_pts)) DEALLOCATE(dist_btw_pts)
        ALLOCATE(diff(atoms,3))
        ALLOCATE(dist_btw_pts(atoms,atoms))
        DO iter=1,atoms 
          DO jter=1,atoms
            diff(jter,:) = inp_coord(iter,:) - coord_set(jter,:)
            dist_btw_pts(iter,jter) = NORM2(diff(jter,:))
          END DO
        END DO

!        PRINT *, ''
!        PRINT *, 'diff'
!        PRINT *, ' '
!        DO iter=1,atoms
!          PRINT '(I6,2X,3(F13.6,2X))', iter,    &
!     &     diff(iter,1), diff(iter,2), diff(iter,3)
!        END DO 

!        PRINT *, ''
!        PRINT *, 'distance between points'
!        PRINT *, 'atom one,two to others'
!        PRINT *, ''
!        DO iter=1,atoms 
!          PRINT *,dist_btw_pts(1,iter),dist_btw_pts(2,iter),&
!     &       dist_btw_pts(3,iter)
!        END DO

        END SUBROUTINE compute_haus_distances

        SUBROUTINE calc_haus_dist(haus_dist)
! Calculate hausdorff distance between sets Q and Q_prime
! only one output value, hausdorff distance
! Reference: "On Quantifying Chiraliy by A.B. Buda, T.Auf der Hyde & K. Mislow,
!            Angew. Chem. Int. Ed. Engl. 1992, 31, 989-1007. 
! Sets Q and Q_prime are variables inp_coord and mir_coord
! Points q and q_prime are contained in sets Q and Q_prime
! Delta(Q_prime,q) denote the shortest distance between a fixed point q and set Q_prime, delta(1,:)
! Delta(Q,q_prime) denote the shortest distance between a fixed point q_prime and set Q, delta(2,:)
        IMPLICIT NONE
        REAL(KIND=DBL),INTENT(OUT)  :: haus_dist
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  ::      delta
        REAL(KIND=DBL)  :: rho_1, rho_2
        INTEGER         :: iter
        IF(ALLOCATED(delta)) DEALLOCATE(delta)
        ALLOCATE(delta(2,atoms))
        DO iter=1,atoms
          delta(1,iter) = MINVAL(dist_btw_pts(iter,:))
          delta(2,iter) = MINVAL(dist_btw_pts(:,iter))
        END DO

!          PRINT *, ''
!        DO iter=1,atoms
!          PRINT *, 'shortest_1', delta(1,iter)
!        END DO
!          PRINT *, ''

!          PRINT *, ''
!        DO iter=1,atoms
!          PRINT *, 'shortest_2', delta(2,iter)
!        END DO
!          PRINT *, ''

        rho_1=MAXVAL(delta(1,:))
        rho_2=MAXVAL(delta(2,:))

        haus_dist = MAX(rho_1,rho_2) 
        haus_dist = haus_dist / diameter
        
        CALL smallest_haus_dist(haus_dist)
        
!        PRINT *,''
!        PRINT *,'Hausdorff distance', haus_dist
!        PRINT *,''
        OPEN(50,file='out_hausdorff.dat',access='append')
        WRITE(50,*) haus_dist
        CLOSE(50)
        END SUBROUTINE calc_haus_dist


        SUBROUTINE haus_translate(trans_param,coord_x)
! translate input coordinates along x,y,z-axes, and update input coordinates
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(IN)    :: trans_param
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)  :: coord_x
        INTEGER                                   :: iter
! Translation operation in x,y,z-axes
        DO iter = 1,3
          coord_x(:,iter) = coord_x(:,iter) + trans_param(iter)
        END DO
        
        END SUBROUTINE haus_translate


        SUBROUTINE haus_rotate(rotate_axis,phi,coord_x)
! rotate input coordinates either along x, y or z-axis, and update input coordinates
        IMPLICIT NONE
        INTEGER,INTENT(IN)                           :: rotate_axis
        REAL(KIND=DBL),INTENT(IN)                    :: phi
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)  :: coord_x
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: temp_x
        REAL(KIND=DBL),DIMENSION(3,3)             :: rotate_x,rotate_z
        REAL(KIND=DBL),DIMENSION(3,3)             :: rotate_y


        IF (ALLOCATED(temp_x)) DEALLOCATE(temp_x)
        ALLOCATE(temp_x(3,atoms))
        temp_x = TRANSPOSE(coord_x)

        SELECT CASE(rotate_axis) 
        CASE(1)
          rotate_x(1,1) = 1.0D0
          rotate_x(1,2) = 0.0D0
          rotate_x(1,3) = 0.0D0
          rotate_x(2,1) = 0.0D0
          rotate_x(2,2) = COS(phi)
          rotate_x(2,3) = -SIN(phi)
          rotate_x(3,1) = 0.0D0
          rotate_x(3,2) = SIN(phi)
          rotate_x(3,3) = COS(phi)
          temp_x = MATMUL(rotate_x,temp_x)
        CASE(2)
          rotate_y(1,1) = COS(phi)
          rotate_y(1,2) = 0.0D0
          rotate_y(1,3) = SIN(phi)
          rotate_y(2,1) = 0.0D0
          rotate_y(2,2) = 1.0D0
          rotate_y(2,3) = 0.0D0
          rotate_y(3,1) = -SIN(phi)
          rotate_y(3,2) = 0.0D0
          rotate_y(3,3) = COS(phi)
          temp_x = MATMUL(rotate_y,temp_x)
        CASE(3)
          rotate_z(1,1) = COS(phi)
          rotate_z(1,2) = -SIN(phi)
          rotate_z(1,3) = 0.0D0
          rotate_z(2,1) = SIN(phi)
          rotate_z(2,2) = COS(phi)
          rotate_z(2,3) = 0.0D0
          rotate_z(3,1) = 0.0D0
          rotate_z(3,2) = 0.0D0
          rotate_z(3,3) = 1.0D0
          temp_x = MATMUL(rotate_z,temp_x)
        END SELECT
        coord_x = TRANSPOSE(temp_x)

        END SUBROUTINE haus_rotate
        
        SUBROUTINE haus_transform(transform_param,&
     &               input_x,output_x,haus_val)
! Perform rotation and translation of coordinates, then calculate Hausdorff distance
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(6),INTENT(IN)    :: transform_param
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: input_x
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)::output_x
        REAL(KIND=DBL),INTENT(OUT)                :: haus_val
        REAL(KIND=DBL),DIMENSION(3)               :: translate_param
        INTEGER                                   :: iter
        IF (ALLOCATED(output_x)) DEALLOCATE(output_x)
        ALLOCATE(output_x(atoms,3))
! Copy transformation parameter to translational parameter
        DO iter=4,6
          translate_param(iter) = transform_param(iter)
        END DO
        output_x = input_x
! Perform rotations in sequence
        CALL haus_rotate(1,transform_param(1),output_x)
        CALL haus_rotate(2,transform_param(2),output_x)
        CALL haus_rotate(3,transform_param(3),output_x)
! Perform translation        
        CALL haus_translate(translate_param,output_x)
! Calculate Hausdorff distance
        CALL compute_haus_distances(output_x)
        CALL calc_haus_dist(haus_val)  


        END SUBROUTINE haus_transform

        SUBROUTINE haus_transform_angle(transform_param&
     &              ,transform_index,input_x,output_x,haus_val)
! Perform rotation along one particular axis, then calculate Hausdorff distance
        IMPLICIT NONE
        REAL(KIND=DBL),INTENT(IN)    :: transform_param
        INTEGER,INTENT(IN)           :: transform_index
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: input_x
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)::output_x
        REAL(KIND=DBL),INTENT(OUT)                :: haus_val
        IF (ALLOCATED(output_x)) DEALLOCATE(output_x)
        ALLOCATE(output_x(atoms,3))
        output_x = input_x
! Perform rotations by case
        SELECT CASE(transform_index)
        CASE(1)
          CALL haus_rotate(1,transform_param,output_x)
        CASE(2)
          CALL haus_rotate(2,transform_param,output_x)
        CASE(3)
          CALL haus_rotate(3,transform_param,output_x)
        END SELECT
! Calculate Hausdorff distance
        CALL compute_haus_distances(output_x)
        CALL calc_haus_dist(haus_val)  

        END SUBROUTINE haus_transform_angle

        SUBROUTINE haus_transform_xyz(transform_param,&
     &               input_x,output_x,haus_val)
! Perform translation of coordinates, then calculate Hausdorff distance
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(IN)    :: transform_param
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: input_x
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)::output_x
        REAL(KIND=DBL),INTENT(OUT)                :: haus_val
        IF (ALLOCATED(output_x)) DEALLOCATE(output_x)
        ALLOCATE(output_x(atoms,3))
        output_x = input_x
! Perform translation        
        CALL haus_translate(transform_param,output_x)
! Calculate Hausdorff distance
        CALL compute_haus_distances(output_x)
        CALL calc_haus_dist(haus_val)  


        END SUBROUTINE haus_transform_xyz



        SUBROUTINE smallest_haus_dist(haus_val)
! find and save the smallest hausdorff distance or HCM after each calc_haus_dist
        IMPLICIT NONE
        REAL(KIND=DBL),INTENT(IN)    :: haus_val
        IF (haus_val < hcm) THEN
          hcm = haus_val
          CALL print_hcm
        END IF
        END SUBROUTINE smallest_haus_dist



        SUBROUTINE print_hcm
! output HCM value
        IMPLICIT NONE
        PRINT *, 'hcm value =', hcm
        OPEN(40,file=TRIM(output_file),status='replace')
          WRITE(40,*) 'Atom size=', atoms
          WRITE(40,*) 'HCM value=', hcm
        CLOSE(40)

        END SUBROUTINE print_hcm


!==========================================================================================
! For all possible rotation scanning
! Testing only, not used
!==========================================================================================
        SUBROUTINE rotate_coord(input_x,output_x,haus_val)
        IMPLICIT NONE
!        REAl(KIND=DBL),DIMENSION(3),INTENT(IN)    :: rotate_param
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: input_x
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(OUT) :: output_x
        REAL(KIND=DBL),INTENT(OUT)                :: haus_val
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: temp_x
        REAL(KIND=DBL),DIMENSION(3,3)             :: rotate_x,rotate_z
        REAL(KIND=DBL),DIMENSION(3,3)             :: rotate_y
        INTEGER                                   :: iter,jter,kter
        INTEGER                                   :: lter
        INTEGER                                   :: intervals
        
        REAL(KIND=DBL)          :: rotate_stride, phi

        intervals = angle_interval
        rotate_stride = 2.0D0 * PI / REAL(intervals)
        phi = rotate_stride

        IF (ALLOCATED(temp_x)) DEALLOCATE(temp_x)
        ALLOCATE(temp_x(3,atoms))
        temp_x = TRANSPOSE(input_x)

             rotate_x(1,1) = 1.0D0
             rotate_x(1,2) = 0.0D0
             rotate_x(1,3) = 0.0D0
             rotate_x(2,1) = 0.0D0
             rotate_x(2,2) = COS(phi)
             rotate_x(2,3) = -SIN(phi)
             rotate_x(3,1) = 0.0D0
             rotate_x(3,2) = SIN(phi)
             rotate_x(3,3) = COS(phi)

             rotate_y(1,1) = COS(phi)
             rotate_y(1,2) = 0.0D0
             rotate_y(1,3) = SIN(phi)
             rotate_y(2,1) = 0.0D0
             rotate_y(2,2) = 1.0D0
             rotate_y(2,3) = 0.0D0
             rotate_y(3,1) = -SIN(phi)
             rotate_y(3,2) = 0.0D0
             rotate_y(3,3) = COS(phi)

             rotate_z(1,1) = COS(phi)
             rotate_z(1,2) = -SIN(phi)
             rotate_z(1,3) = 0.0D0
             rotate_z(2,1) = SIN(phi)
             rotate_z(2,2) = COS(phi)
             rotate_z(2,3) = 0.0D0
             rotate_z(3,1) = 0.0D0
             rotate_z(3,2) = 0.0D0
             rotate_z(3,3) = 1.0D0
!        temp_x = MATMUL(rotate,temp_x)
!        output_x = TRANSPOSE(temp_x)


!        CALL compute_haus_distances(temp_x)
!        CALL calc_haus_dist(haus_val)  
!        output_x = temp_x
!        PRINT *, 'Haus value='
!        PRINT *, haus_val


        OPEN(60,file='rotate.xyz',status='replace')
        outer: DO lter = 1,intervals


          temp_x = MATMUL(rotate_y,temp_x)
          output_x = TRANSPOSE(temp_x)
          WRITE(60,*) '15'
          WRITE(60,*) haus_val
          DO iter =1,atoms
            WRITE(60,'(A,2X,3(F13.6,2X))') 'Au', &
     &         output_x(iter,1),output_x(iter,2),output_x(iter,3)
          END DO

        longitude: DO kter =1, intervals
          temp_x = MATMUL(rotate_x,temp_x)
          output_x = TRANSPOSE(temp_x)
          WRITE(60,*) '15'
          WRITE(60,*) haus_val
          DO iter =1,atoms
            WRITE(60,'(A,2X,3(F13.6,2X))') 'Au', &
     &         output_x(iter,1),output_x(iter,2),output_x(iter,3)
          END DO


        latitude: DO jter = 1,intervals
           !phi = rotate_stride * REAL(iter)
          temp_x = MATMUL(rotate_z,temp_x)
          output_x = TRANSPOSE(temp_x)

          CALL compute_haus_distances(output_x)
          CALL calc_haus_dist(haus_val)  
          WRITE(60,*) '15'
          WRITE(60,*) haus_val
          DO iter =1,atoms
!        PRINT '(I6,2X,3(F13.6,2X))', iter, &
!     &   output_x(iter,1),output_x(iter,2),output_x(iter,3)

            WRITE(60,'(A,2X,3(F13.6,2X))') 'Au', &
     &         output_x(iter,1),output_x(iter,2),output_x(iter,3)
          END DO
        END DO latitude
        END DO longitude
        END DO outer

        CLOSE(60)



        
        END SUBROUTINE rotate_coord
        END MODULE hausdorff        
