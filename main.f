!
! Driver programme for chirality measure and Haussforf distance
! Update: 1 Nov 2019
!
        PROGRAM main

        USE initialise
        USE simulate

        IMPLICIT NONE
        INTEGER :: c1,c2,c_rate
        REAL    :: t1,t2,t3,t4

        CALL SYSTEM_CLOCK(c1,c_rate)
        PRINT *, '**************************************'
        PRINT *, 'Driver programme for chirality measure'
        PRINT *, '**************************************'
        PRINT *, ''

        CALL start_session       
        CALL single_simul

        CALL SYSTEM_CLOCK(c2,c_rate)
        PRINT *, "c1", c1
        PRINT *, "c2", c2
        PRINT *, "c_rate", c_rate
        PRINT *, "Time used, SYSTEM_CLOCK, in seconds:", c2-c1
        PRINT *, "Time used, SYSTEM_CLOCK, in seconds:", &
     &              REAL(c2-c1)/REAL(c_rate)
        END PROGRAM main
