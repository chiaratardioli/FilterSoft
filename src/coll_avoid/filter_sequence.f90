!
! Chiara Tardioli (tardioli@mail.dm.unipi.it)
! Version: August 29, 2012
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *           T H R E E   F I L T E R    S E Q U E N C E          * 
!  *                                                               *    
!  *               Collision warnings for input orbits             *
!  *               in a ephemeris interpolation table              *    
!  *                                                               *
!  *****************************************************************    
!
PROGRAM filter_sequence
  USE fund_const
  USE option_file
  USE read_interpol_ephem
  USE filter1
  USE filter2
  USE filter3
  USE filter_output
  IMPLICIT NONE
! ------------------------------------------------------------------
! HINT : elements are converted in EQU, evolution time in MJD
! ------------------------------------------------------------------
  INTEGER :: iunda           !direct access iunit
  INTEGER :: orb1,orb2,norb  !no. of orbits
! for hoots_1
  LOGICAL :: apoperi
  INTEGER :: ntot,sum1,sum2,sum3 !results
! auxiliar variables
!  INTEGER :: j !loop index
  LOGICAL :: check_diff=.true.
  REAL(KIND=dkind) :: tstep
  INTEGER :: k,kmax,cross1,cross_curr
  INTEGER,DIMENSION(norbx,norbx) :: hoots_save
! ------------------------------------------------------------------
! OPTIONS
  CALL read_options
  write(*,*)'************************************'
  write(*,*)'options'
  write(*,*)'rhs=',rhs
!  write(*,*)'tmax=',tmax
  write(*,*)'initial epoch=',teph0
!  write(*,*)'dt=',dt
!  write(*,*)'d_thres=',d_thres
  write(*,*)'************************************'

! READ ORBIT INDEX FROM TERMINAL
  CALL read_orbit_index(orb1,orb2,norb)
  write(*,*)'first and last orbit index : ',orb1,orb2
  write(*,*)'no. of orbits in input=',norb
!  write(*,*)'************************************'

! EPHEMERIDES INTERPOLATION TABLE
  CALL read_eph_table(iunda,nevol)
  write(*,*)'no. of time nodes =',nevol
  write(*,*)'************************************'

! ================================================
! START THREE FILTER SEQUENCE
! ================================================
! initilal
  ntot=norb*(norb-1)/2
!  write(*,*)'************************************'
  write(*,*)'no. of couples =',ntot
  write(*,*)'************************************'

! ================================================
!          FILTER I (geocentric distance)
! output NxN matrix : hoots (declare in module filter_output)
! ================================================
! setting (change dt if different from the option file)
  CALL set_hevol(iunda,nevol,dt,hevol)

  CALL check_eph_table(iunda,0)
!  stop

  apoperi=.FALSE. !compute the apogee-perigee filter and compare
  CALL filter_one(iunda,orb1,orb2,apoperi)

! results
  CALL write_results(1)
! ================================================

! ================================================
!        FILTER II (orbit distance)
! output NxN matrix : dcros (declare in module filter_output)
! ================================================
! settings (change dt if different from the option file)
! set dt=0.025d0 : same step as for the ephemeris evolution

! save hoots matrix and dt
  hoots_save = hoots
  tstep = dt

! initialise dcros matrix
  dcros = 0.d0
  cross1 = sum(dcros)

! Extract maximum no. of iterations
  kmax = 1   !CEILING(1.d0/(2.d0*pig) / dt)
!  write(*,*)'kmax=',kmax

  DO k=1,kmax
     dt = 1.d0/(2.d0*pig) / k
     write(*,*)'it',k,'dt = ',dt

     CALL set_hevol(iunda,nevol,dt,hevol)
     CALL filter_two(iunda,orb1,orb2)

!  write(*,*)'evolution step =',dt
!  write(*,*)'no. of evolutions =',hevol

! Tmp matrix
     hoots = hoots - dcros
     if(sum(hoots).eq.0.d0)then
        GOTO 5
     endif

!!!!!!! THIS NOT TRUE BUT FAST !!!!!!!
     cross_curr = sum(dcros)
     write(*,*)cross_curr,cross1,cross_curr-cross1
     if(cross_curr-cross1.eq.0)then
        GOTO 5 ! no more cross are detected
     elseif(cross_curr-cross1.lt.0)then
        write(*,*)'error filter II: no. of crossings',cross_curr,cross1
     endif
     cross1 = cross_curr
     CALL write_results(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ENDDO
!  write(*,*)'filter II : all couples pass!'

5 CONTINUE
! Restore hoots matrix
  hoots = hoots_save

! results
  CALL write_results(2)
! ================================================

! ================================================
!      DIFFERENCE BETWEEN FILTER I AND II
! ================================================
  IF(check_diff)THEN
     CALL diff_filter_1_2(orb1,orb2)
  ENDIF
! ================================================

! ================================================
!        FILTER III (time filter)
! output NxN matrix : dtime (declare in module filter_output)
! ================================================
! settinsg (change dt if different from the option file)
  CALL set_hevol(iunda,nevol,dt,hevol)
  
  CALL filter_three(iunda,orb1,orb2)
  
! results
  CALL write_results(3)
! ================================================

! ================================================
! END THREE FILTER SEQUENCE
! ================================================
  CALL filclo(iunda,' ')

END PROGRAM filter_sequence


! ************************************************************
!        READ ORBIT INDEX FROM TERMINAL
! ************************************************************
SUBROUTINE read_orbit_index(orb1,orb2,norb)
  USE fund_const
  USE option_file,only:norbx
  IMPLICIT NONE
! -------------------- interface --------------------
  INTEGER,INTENT(OUT) :: orb1,orb2,norb
! ------------------ end interface ------------------
  write(*,*)'Insert first and last orbit index (from 1 to',norbx,')'
  read(*,*)orb1,orb2
  IF(orb1.gt.norbx.OR.orb1.lt.0.OR.orb2.gt.norbx.OR. & 
       & orb2.lt.0.OR.orb1.eq.orb2)THEN
     write(*,*)'number not in range'
     STOP
  ELSE
     IF(orb1.lt.orb2)THEN
        norb=orb2-orb1+1
     ELSE
        write(*,*)orb1,orb2,norb
        norb=orb1-orb2
        orb2=orb1
        orb1=orb2-norb
        norb=norb+1
     ENDIF
  ENDIF
END SUBROUTINE read_orbit_index


! ************************************************************
!        READ ORBIT EVOLUTIONS
! ************************************************************
SUBROUTINE set_hevol(iunda,nvl,dt,hvl)
  USE fund_const
  IMPLICIT NONE
! -------------------- interface --------------------
  INTEGER,INTENT(IN) :: iunda,nvl
  REAL(KIND=dkind),INTENT(IN) :: dt
  INTEGER,INTENT(OUT) :: hvl
! ------------------ end interface ------------------
  REAL(KIND=dkind) :: tevol
! ---------------------------------------------------
  READ(iunda,REC=nvl)tevol
  hvl=FLOOR(tevol/dt)
END SUBROUTINE set_hevol
