
        SUBROUTINE read_one_coord(struc_file)
! read xyz coordinates file
        IMPLICIT NONE
        CHARACTER(20),INTENT(IN)   :: struc_file
        INTEGER         :: iter
        CHARACTER       :: dummy
        OPEN(21,file=struc_file,status='old')
        READ(21,*) atoms
        READ(21,*) dummy
        IF (ALLOCATED(coord)) DEALLOCATE(coord)
        ALLOCATE(coord(atoms,3))
        DO iter=1,atoms
          READ(21,*) dummy, coord(iter,1),coord(iter,2),coord(iter,3)
        END DO
        CLOSE(21)
        END SUBROUTINE read_one_coord

        
        SUBROUTINE print_coord
! print xyz coordinate files in columns
        IMPLICIT NONE
        INTEGER  :: iter
        PRINT *, ''
        PRINT *, 'Coordinates in xyz'
        PRINT *, ''
        DO iter =1,atoms 
          PRINT *, iter, coord(iter,1), coord(iter,2), coord(iter,3)
        END DO
        PRINT *, ''
        END SUBROUTINE print_coord


