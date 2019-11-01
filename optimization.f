!                                                                                      
!  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
!  or “3-clause license”)                                                              
!  Please read attached file License.txt                                               
!                                        
!
!                      DRIVER1 in Fortran 90
!     --------------------------------------------------------------
!
!        L-BFGS-B is a code for solving large nonlinear optimization
!             problems with simple bounds on the variables.
!
!        The code can also be used for unconstrained problems and is
!        as efficient for these problems as the earlier limited memory
!                          code L-BFGS.
!
!        This is the simplest driver in the package. It uses all the
!                    default settings of the code.
!
!
!     References:
!
!        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!        memory algorithm for bound constrained optimization'',
!        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!        Subroutines for Large Scale Bound Constrained Optimization''
!        Tech. Report, NAM-11, EECS Department, Northwestern University,
!        1994.
!
!
!          (Postscript files of these papers are available via anonymous
!           ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                              *  *  *
!
!         March 2011   (latest revision)
!         Optimization Center at Northwestern University
!         Instituto Tecnologico Autonomo de Mexico
!
!         Jorge Nocedal and Jose Luis Morales
!
!     --------------------------------------------------------------
!             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
!     --------------------------------------------------------------
!
!     n is an INTEGER variable that must be set by the user to the
!       number of variables.  It is not altered by the routine.
!
!     m is an INTEGER variable that must be set by the user to the
!       number of corrections used in the limited memory matrix.
!       It is not altered by the routine.  Values of m < 3  are
!       not recommended, and large values of m can result in excessive
!       computing time. The range  3 <= m <= 20 is recommended. 
!
!     x is a DOUBLE PRECISION array of length n.  On initial entry
!       it must be set by the user to the values of the initial
!       estimate of the solution vector.  Upon successful exit, it
!       contains the values of the variables at the best point
!       found (usually an approximate solution).
!
!     l is a DOUBLE PRECISION array of length n that must be set by
!       the user to the values of the lower bounds on the variables. If
!       the i-th variable has no lower bound, l(i) need not be defined.
!
!     u is a DOUBLE PRECISION array of length n that must be set by
!       the user to the values of the upper bounds on the variables. If
!       the i-th variable has no upper bound, u(i) need not be defined.
!
!     nbd is an INTEGER array of dimension n that must be set by the
!       user to the type of bounds imposed on the variables:
!       nbd(i)=0 if x(i) is unbounded,
!              1 if x(i) has only a lower bound,
!              2 if x(i) has both lower and upper bounds, 
!              3 if x(i) has only an upper bound.
!
!     f is a DOUBLE PRECISION variable.  If the routine setulb returns
!       with task(1:2)= 'FG', then f must be set by the user to
!       contain the value of the function at the point x.
!
!     g is a DOUBLE PRECISION array of length n.  If the routine setulb
!       returns with taskb(1:2)= 'FG', then g must be set by the user to
!       contain the components of the gradient at the point x.
!
!     factr is a DOUBLE PRECISION variable that must be set by the user.
!       It is a tolerance in the termination test for the algorithm.
!       The iteration will stop when
!
!        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
!
!       where epsmch is the machine precision which is automatically
!       generated by the code. Typical values for factr on a computer
!       with 15 digits of accuracy in double precision are:
!       factr=1.d+12 for low accuracy;
!             1.d+7  for moderate accuracy; 
!             1.d+1  for extremely high accuracy.
!       The user can suppress this termination test by setting factr=0.
!
!     pgtol is a double precision variable.
!       On entry pgtol >= 0 is specified by the user.  The iteration
!         will stop when
!
!                 max{|proj g_i | i = 1, ..., n} <= pgtol
!
!         where pg_i is the ith component of the projected gradient.
!       The user can suppress this termination test by setting pgtol=0.
!
!     wa is a DOUBLE PRECISION  array of length 
!       (2mmax + 5)nmax + 11mmax^2 + 8mmax used as workspace.
!       This array must not be altered by the user.
!
!     iwa is an INTEGER  array of length 3nmax used as
!       workspace. This array must not be altered by the user.
!
!     task is a CHARACTER string of length 60.
!       On first entry, it must be set to 'START'.
!       On a return with task(1:2)='FG', the user must evaluate the
!         function f and gradient g at the returned value of x.
!       On a return with task(1:5)='NEW_X', an iteration of the
!         algorithm has concluded, and f and g contain f(x) and g(x)
!         respectively.  The user can decide whether to continue or stop
!         the iteration. 
!       When
!         task(1:4)='CONV', the termination test in L-BFGS-B has been 
!           satisfied;
!         task(1:4)='ABNO', the routine has terminated abnormally
!           without being able to satisfy the termination conditions,
!           x contains the best approximation found,
!           f and g contain f(x) and g(x) respectively;
!         task(1:5)='ERROR', the routine has detected an error in the
!           input parameters;
!       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
!         contains additional information that the user can print.
!       This array should not be altered unless the user wants to
!          stop the run for some reason.  See driver2 or driver3
!          for a detailed explanation on how to stop the run 
!          by assigning task(1:4)='STOP' in the driver.
!
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     csave  is a CHARACTER working array of length 60.
!
!     lsave is a LOGICAL working array of dimension 4.
!       On exit with task = 'NEW_X', the following information is
!         available:
!       lsave(1) = .true.  the initial x did not satisfy the bounds;
!       lsave(2) = .true.  the problem contains bounds;
!       lsave(3) = .true.  each variable has upper and lower bounds.
!
!     isave is an INTEGER working array of dimension 44.
!       On exit with task = 'NEW_X', it contains information that
!       the user may want to access:
!         isave(30) = the current iteration number;
!         isave(34) = the total number of function and gradient
!                         evaluations;
!         isave(36) = the number of function value or gradient
!                                  evaluations in the current iteration;
!         isave(38) = the number of free variables in the current
!                         iteration;
!         isave(39) = the number of active constraints at the current
!                         iteration;
!
!         see the subroutine setulb.f for a description of other 
!         information contained in isave
!
!     dsave is a DOUBLE PRECISION working array of dimension 29.
!       On exit with task = 'NEW_X', it contains information that
!         the user may want to access:
!         dsave(2)  = the value of f at the previous iteration;
!         dsave(5)  = the machine precision epsmch generated by the code;
!         dsave(13) = the infinity norm of the projected gradient;
!
!         see the subroutine setulb.f for a description of other 
!         information contained in dsave
!
!     --------------------------------------------------------------
!           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
!     --------------------------------------------------------------
!    
      MODULE optimization
      USE constants
      USE initialise
      USE hausdorff,    ONLY:haus_transform,haus_transform_xyz,&
     &                    haus_transform_angle
      USE deriv,        ONLY:haus_deriv,haus_deriv_xyz,&
     &                    haus_deriv_angle

      CONTAINS
      SUBROUTINE optim_lbfgs(in_coord)

!
!     This simple driver demonstrates how to call the L-BFGS-B code to
!       solve a sample problem (the extended Rosenbrock function 
!       subject to bounds on the variables). The dimension n of this
!       problem is variable.

      IMPLICIT NONE
!
!     Declare variables and parameters needed by the code.
!       Note thar we wish to have output at every iteration. 
!          iprint=1
!
!       We also specify the tolerances in the stopping criteria.
!          factr  = 1.0d+7, pgtol  = 1.0d-5
!
!       A description of all these variables is given at the beginning
!       of  the driver 
!
!======================================================================
!     Variables essential to LBFGS subroutine
!======================================================================
      integer                :: n 
      integer,  parameter    :: m = 20, iprint = 1
      real(KIND=DBL), parameter    :: factr  = 1.0d+1, pgtol  = 1.0d-12
!
      character(len=60)      :: task, csave
      logical                :: lsave(4)
      integer                :: isave(44)
      real(KIND=DBL)               :: f
      real(KIND=DBL)               :: dsave(29)
      integer,  allocatable  :: nbd(:), iwa(:)
      real(KIND=DBL), allocatable  :: x(:), l(:), u(:), g(:), wa(:)

!======================================================================
!     Declare a few additional variables for this problem
!======================================================================
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: in_coord
      INTEGER                      :: iter
      REAL(KIND=DBL)               :: haus_val
      REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE  :: grad
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: out_coord
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: temp_coord
!======================================================================
!     Determine the number of variables, n
!======================================================================
      n=6   ! 3 translations + 3 rotations
!======================================================================
!     Allocate dynamic arrays
!======================================================================
      ALLOCATE ( grad(n) )
      ALLOCATE ( nbd(n), x(n), l(n), u(n), g(n) )
      ALLOCATE ( iwa(3*n) )
      ALLOCATE ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
      IF (ALLOCATED(out_coord)) DEALLOCATE(out_coord)
      IF (ALLOCATED(temp_coord)) DEALLOCATE(temp_coord)
      ALLOCATE(out_coord(atoms,3))
      ALLOCATE(temp_coord(atoms,3))
!======================================================================
!     Set bounds for variables, x is unbounded, set to 0 
!======================================================================
      DO iter=1, n
         nbd(iter) = 0
      END DO
!======================================================================
!     Define the starting Hausdorff transformation values
!     all set to 0, i.e. no transformation done initially
!     x(1) = x-axis translation
!     x(2) = y-axis translation
!     x(3) = z-axis translation
!     x(4) = rotation about x-axis 
!     x(5) = rotation about y-axis 
!     x(6) = rotation about z-axis 
      x(1) = 0.0D0
      x(2) = 0.0D0
      x(3) = 0.0D0
      x(4) = 0.0D0
      x(5) = 0.0D0
      x(6) = 0.0D0
!======================================================================
!     Optimiztion loop starts here
!======================================================================
!     We start the iteration by initializing task.
 
      task = 'START'

 
      do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
     &          task.eq.'START') 
         
!     This is the call to the L-BFGS-B code.
        
        
         
         call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
     &                  wa, iwa, task, iprint, &
     &                  csave, lsave, isave, dsave )
         
         if (task(1:2) .eq. 'FG') then
!======================================================================
!     Calculate Hausdorff distance
!======================================================================
            temp_coord = in_coord        !use temp_coord, because haus_transform updates input coordinates
            CALL haus_transform(x,temp_coord,out_coord,haus_val)
            f = haus_val
!======================================================================
!     Call first deriv
!======================================================================
            temp_coord = in_coord        !use temp_coord, because haus_deriv updates input coordinates
            CALL haus_deriv(x,temp_coord,grad)
            g = grad
         end if
      end do

!======================================================================
!     End of do while loop
!======================================================================
      END SUBROUTINE optim_lbfgs

      SUBROUTINE optim_lbfgs_xyz(in_coord)

!
!     This simple driver demonstrates how to call the L-BFGS-B code to
!       solve a sample problem (the extended Rosenbrock function 
!       subject to bounds on the variables). The dimension n of this
!       problem is variable.

      IMPLICIT NONE
!
!     Declare variables and parameters needed by the code.
!       Note thar we wish to have output at every iteration. 
!          iprint=1
!
!       We also specify the tolerances in the stopping criteria.
!          factr  = 1.0d+7, pgtol  = 1.0d-5
!
!       A description of all these variables is given at the beginning
!       of  the driver 
!
!======================================================================
!     Variables essential to LBFGS subroutine
!======================================================================
      integer                :: n 
      integer,  parameter    :: m = 20, iprint = 1
      real(KIND=DBL), parameter    :: factr  = 1.0d+1, pgtol  = 1.0d-12
!
      character(len=60)      :: task, csave
      logical                :: lsave(4)
      integer                :: isave(44)
      real(KIND=DBL)               :: f
      real(KIND=DBL)               :: dsave(29)
      integer,  allocatable  :: nbd(:), iwa(:)
      real(KIND=DBL), allocatable  :: x(:), l(:), u(:), g(:), wa(:)

!======================================================================
!     Declare a few additional variables for this problem
!======================================================================
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: in_coord
      INTEGER                      :: iter
      REAL(KIND=DBL)               :: haus_val
      REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE  :: grad
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: out_coord
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: temp_coord
!======================================================================
!     Determine the number of variables, n
!======================================================================
      n=3   ! 3 translations 
!======================================================================
!     Allocate dynamic arrays
!======================================================================
      ALLOCATE ( grad(n) )
      ALLOCATE ( nbd(n), x(n), l(n), u(n), g(n) )
      ALLOCATE ( iwa(3*n) )
      ALLOCATE ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
      IF (ALLOCATED(out_coord)) DEALLOCATE(out_coord)
      IF (ALLOCATED(temp_coord)) DEALLOCATE(temp_coord)
      ALLOCATE(out_coord(atoms,3))
      ALLOCATE(temp_coord(atoms,3))
!======================================================================
!     Set bounds for variables, x is unbounded, set to 0 
!======================================================================
      DO iter=1, n
         nbd(iter) = 0
      END DO
!======================================================================
!     Define the starting Hausdorff transformation values
!     all set to 0, i.e. no transformation done initially
!     x(1) = x-axis translation
!     x(2) = y-axis translation
!     x(3) = z-axis translation
      x(1) = 0.0D0
      x(2) = 0.0D0
      x(3) = 0.0D0
!======================================================================
!     Optimiztion loop starts here
!======================================================================
!     We start the iteration by initializing task.
 
      task = 'START'

 
      do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
     &          task.eq.'START') 
         
!     This is the call to the L-BFGS-B code.
        
        
         
         call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
     &                  wa, iwa, task, iprint, &
     &                  csave, lsave, isave, dsave )
         
         if (task(1:2) .eq. 'FG') then
!======================================================================
!     Calculate Hausdorff distance
!======================================================================
            temp_coord = in_coord        !use temp_coord, because haus_transform updates input coordinates
            CALL haus_transform_xyz(x,temp_coord,out_coord,haus_val)
            f = haus_val
!======================================================================
!     Call first deriv
!======================================================================
            temp_coord = in_coord        !use temp_coord, because haus_deriv updates input coordinates
            CALL haus_deriv_xyz(x,temp_coord,grad)
            g = grad
         end if
      end do

!======================================================================
!     End of do while loop
!======================================================================
      END SUBROUTINE optim_lbfgs_xyz



      SUBROUTINE optim_lbfgs_angle(angle_index,in_coord,min_angle)

!
!     This simple driver demonstrates how to call the L-BFGS-B code to
!       solve a sample problem (the extended Rosenbrock function 
!       subject to bounds on the variables). The dimension n of this
!       problem is variable.

      IMPLICIT NONE
!
!     Declare variables and parameters needed by the code.
!       Note thar we wish to have output at every iteration. 
!          iprint=1
!
!       We also specify the tolerances in the stopping criteria.
!          factr  = 1.0d+7, pgtol  = 1.0d-5
!
!       A description of all these variables is given at the beginning
!       of  the driver 
!
!======================================================================
!     Variables essential to LBFGS subroutine
!======================================================================
      integer                :: n 
      integer,  parameter    :: m = 20, iprint = 1
      real(KIND=DBL), parameter    :: factr  = 1.0d+1, pgtol  = 1.0d-12
!
      character(len=60)      :: task, csave
      logical                :: lsave(4)
      integer                :: isave(44)
      real(KIND=DBL)               :: f
      real(KIND=DBL)               :: dsave(29)
      integer,  allocatable  :: nbd(:), iwa(:)
      real(KIND=DBL), allocatable  :: x(:), l(:), u(:), g(:), wa(:)

!======================================================================
!     Declare a few additional variables for this problem
!======================================================================
      INTEGER,INTENT(IN)        :: angle_index
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: in_coord
      REAL(KIND=DBL),INTENT(OUT)   :: min_angle
      INTEGER                      :: iter
      REAL(KIND=DBL)               :: haus_val
      REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE  :: grad
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: out_coord
      REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: temp_coord
!convert from array to scalar for functions haus_transform_angle
! & haus_deriv_angle
      REAL(KIND=DBL)               :: x_scalar  = 0.0D0
      REAL(KIND=DBL)               :: g_scalar  = 0.0D0
!======================================================================
!     Determine the number of variables, n
!======================================================================
      n=1   ! 1 rotation at a time
!======================================================================
!     Allocate dynamic arrays
!======================================================================
      ALLOCATE ( grad(n) )
      ALLOCATE ( nbd(n), x(n), l(n), u(n), g(n) )
      ALLOCATE ( iwa(3*n) )
      ALLOCATE ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
      IF (ALLOCATED(out_coord)) DEALLOCATE(out_coord)
      IF (ALLOCATED(temp_coord)) DEALLOCATE(temp_coord)
      ALLOCATE(out_coord(atoms,3))
      ALLOCATE(temp_coord(atoms,3))
!======================================================================
!     Set bounds for variables, x is unbounded, set to 0 
!======================================================================
      DO iter=1, n
         nbd(iter) = 0
      END DO
!======================================================================
!     Define the starting Hausdorff transformation values
!     all set to 0, i.e. no transformation done initially
!     x(1) = x-axis translation
!     x(2) = y-axis translation
!     x(3) = z-axis translation
      x(1) = 0.0D0
!======================================================================
!     Optimiztion loop starts here
!======================================================================
!     We start the iteration by initializing task.
 
      task = 'START'

 
      do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
     &          task.eq.'START') 
         
!     This is the call to the L-BFGS-B code.
        
        
         
         call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
     &                  wa, iwa, task, iprint, &
     &                  csave, lsave, isave, dsave )
         
         if (task(1:2) .eq. 'FG') then
!======================================================================
!from x array to x scalar for the next two subroutine calls
!======================================================================
            x_scalar = x(1) 
!======================================================================
!     Calculate Hausdorff distance
!======================================================================
            temp_coord = in_coord        !use temp_coord, because haus_transform updates input coordinates
            CALL haus_transform_angle(x_scalar,angle_index,&
     &         temp_coord,out_coord,haus_val)
            f = haus_val
!======================================================================
!     Call first deriv
!======================================================================
            temp_coord = in_coord        !use temp_coord, because haus_deriv updates input coordinates
            CALL haus_deriv_angle(x_scalar,angle_index,&
     &          temp_coord,g_scalar)
            g = g_scalar
         end if
      end do

!======================================================================
!     End of do while loop
!======================================================================
      min_angle = x(1)
      END SUBROUTINE optim_lbfgs_angle


      END MODULE optimization

