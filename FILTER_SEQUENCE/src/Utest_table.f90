! Chiara Tardioli (tardioli@mail.dm.unipi.it)
! --------------------------------------------------------------------- 
! Simple program to write on file one evolution and its inerpolation tecnique
! using the Ephemeris Interpolation Table
PROGRAM Utest_table
  USE fund_const
  USE option_file
  USE read_interpol_ephem
  IMPLICIT NONE
  ! ------------------------------------------------------------------
  ! HINT : elements are converted in EQU, evolution time in MJD
  ! ------------------------------------------------------------------
  INTEGER :: iunda           !direct access iunit
  INTEGER :: idx             !index of the orbit to plot
  REAL(KIND=dkind) :: tevol(2)

  ! ------------------------------------------------------------------
  !        Main program
  ! ------------------------------------------------------------------

! OPTIONS
  CALL read_options
  write(*,*)'************************************'
  write(*,*)'options'
!  write(*,*)'rhs=',rhs
!  write(*,*)'tmax=',tmax
  write(*,*)'initial epoch=',teph0
  write(*,*)'************************************'

! READ ORBIT INDEX FROM TERMINAL
  write(*,*)'Index orbit to plot : '
  read(*,*) idx
!  write(*,*)'************************************'

! EPHEMERIDES INTERPOLATION TABLE
  CALL read_eph_table(iunda,nevol)
  write(*,*)'no. of time nodes =',nevol
!  write(*,*)'************************************'

  READ(iunda,REC=1)tevol(1)
  READ(iunda,REC=2)tevol(2)
  write(*,*)'Time step in the integration: ',tevol(2)-tevol(1)
  write(*,*)'Evolution time(days)        : ',(tevol(2)-tevol(1)) * nevol

! READ time step
  write(*,*)'Set a time step (default is ',dt,') : '
  read(*,*) dt
!  write(*,*)'************************************'

! Set number of evolution accordinto the time step
  CALL set_hevol(iunda) !,nevol,dt,hevol)

! Write on file ephemeris interpolation table of object at pos = (idx-1)*nevol
  CALL check_eph_table(iunda,(idx-1)*nevol)


END PROGRAM Utest_table
