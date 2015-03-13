!
! Chiara Tardioli (tardioli@mail.dm.unipi.it)
! Version: August 29, 2012
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *   E P H E M E R I S   I N T E R P O L A T I O N   T A B L E   *
!  *                                                               *    
!  *        Generation of an orbit evolution for all orbits        *
!  *                 in a initial EQU ortbit catalog               *
!  *                    (no orbit uncetainty)                      *
!  *                                                               *
!  *        READ ORBIT EVOLUTIONS COMPUTED BY ORDER PROGRAM        *
!  *                                                               *    
!  *****************************************************************    
! Reference: Charles Hubaux (2012)
!
PROGRAM ephem_table
  USE fund_const
  IMPLICIT NONE
! ------------------------------------------------------------------
! HINT : elements are converted in EQU, evolution time in MJD
! ------------------------------------------------------------------
! for direct access
  REAL(KIND=dkind) :: tevol,elem(6)
  CHARACTER*30 :: elefil
  INTEGER :: le,iunda,length
! for read_elkep
  CHARACTER*20 :: filename
  INTEGER :: n1,nevol !no. of evolution for the orbital elements
  INTEGER :: norb,j !orbit number and index
  INTEGER :: pos !position to start a new evolution in binary file
! ------------------------------------------------------------------
  rhs=2
!  write(*,*)'rhs=',rhs
!  teph0=55000
  write(*,*)'initial epoch=',55000

3 write(*,*)'Insert no. of orbits:'
  read(*,*)norb
!  write(*,*)'no. orbits=',norb
  IF(norb.le.0.d0)THEN
     write(*,*)'error : no. of orbits <= 0',norb
     GOTO 3
  ENDIF

! ephemerides interpolaton table
  INQUIRE(IOLENGTH=length)tevol,elem(1:6)
  CALL filnam('.','all_orbit_evolutions','a',elefil,le)
  CALL fildiropn(iunda,elefil(1:le),'UNKNOWN',length)

! Read orbit evolutions and write them in a binary file
  CALL read_elkep(1,0,iunda,n1) !elements are EQU in output
  nevol=n1
  pos=nevol
  DO j=2,norb
     CALL read_elkep(j,pos,iunda,n1) !elements are EQU in output
     IF(n1.ne.nevol)THEN
        write(*,*)'error : different evolution span',n1,nevol
        STOP
     ENDIF
     pos=pos+nevol
  ENDDO

! check direct access
!  READ(iunda,REC=1)tevol,elem
!  write(*,*)tevol,elem

  CALL filclo(iunda,' ')

END PROGRAM ephem_table

! ************************************************************
!    READ ORBIT EVOLUTIONS COMPUTED BY ORDER (charles)
! ************************************************************
SUBROUTINE read_elkep(j,pos,iunda,nvl)
  USE fund_const
  IMPLICIT NONE
  INTEGER,PARAMETER :: tmax=3000
! -------------------- interface --------------------
  INTEGER,INTENT(IN) :: j !file index
  INTEGER,INTENT(IN) :: pos,iunda !for direct access
!  REAL(KIND=dkind),INTENT(OUT) :: tvl(tmax),evolel(6,tmax)
  INTEGER,INTENT(OUT) :: nvl
! ------------------ end interface ------------------
  CHARACTER*60 :: elefil !file names
  INTEGER :: le
  CHARACTER*10 :: jchar !integer j trasformed into a character variable
  REAL(KIND=dkind) :: tvl,evolel(6)
  INTEGER :: iunorb,i
! ---------------------------------------------------

! Open file
  write(jchar,'(I10)')j
  jchar=adjustl(jchar)
  CALL filnam('./datafile','results_kepl_'//trim(jchar),'dat',elefil,le)
!  write(*,*)elefil(1:le)
  CALL filopn(iunorb,elefil(1:le),'OLD')

! Read file
  DO i=1,tmax
     READ(iunorb,*,END=998) tvl,evolel(1:6)
     CALL convert_elem(tvl,evolel) !elements are EQU in output
! Write one line evolution
     WRITE(iunda,REC=pos+i)tvl,evolel(1:6)
  ENDDO
  write(*,*)'read_elkep : do loop too short!',i
  STOP
998 CONTINUE
  CALL filclo(iunorb,' ')
  nvl=i-1

ENDSUBROUTINE read_elkep
! ============================================================

! ********************************************************
!              ORBITAL ELEMENTS UNIT CONVERSION
!
! First, convert to KEP : a[km],angles[rad]
! Then, convert to EQU  : a,h=e*cos(om),k=e*sin(om),p=tan(i/2)*cos(Om)
!                         q=tan(i/2)*sin(Om),l=mean anomaly
! Time : conversion from [days/pig] to [days]
! ********************************************************
SUBROUTINE convert_elem(tvl,elem)
  USE fund_const
  USE orbit_elements
!  USE filter1,only:tmax,teph0
  IMPLICIT NONE
  REAL(KIND=dkind),PARAMETER :: ys=3.6d3*24.d0
  REAL(KIND=dkind),PARAMETER :: UD=42164.1697748545d0
  INTEGER,PARAMETER :: tmax=3000
  REAL(KIND=dkind),PARAMETER :: teph0=55000
! -------------------- interface --------------------
  REAL(KIND=dkind),INTENT(INOUT) :: tvl,elem(6)
! ------------------ end interface ------------------
  TYPE(orbit_elem) :: ele
  INTEGER :: i,fail_flag
! ===================================================
  tvl=tvl*dpig        !from days/2pi to days
  elem(1)=elem(1)*UD      !to km

  ele=undefined_orbit_elem
  ele%t=teph0+tvl !in MJD
  ele%coo='KEP'
  ele%coord(1:6)=elem(1:6)
  CALL coo_cha(ele,'EQU',ele,fail_flag)
  IF(fail_flag.ge.5) THEN
     WRITE(*,*)'convert_elem, coo_cha: fail_flag=',fail_flag
     STOP
  ENDIF
  elem(1:6)=ele%coord(1:6)

ENDSUBROUTINE convert_elem
! ================================================

