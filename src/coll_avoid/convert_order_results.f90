!
! Chiara Tardioli (c.tardioli@strath.ac.uk)
! Version: August 31, 2014
! --------------------------------------------------------------------- 
!
PROGRAM convert_order_results
  USE fund_const
  IMPLICIT NONE
! ------------------------------------------------------------------
! HINT : elements are converted in EQU, evolution time in MJD
! ------------------------------------------------------------------
  INTEGER,PARAMETER :: tmax=5000
! ------------------- interface ----------------------------
  REAL(KIND=dkind) :: tevol,elem(6)
! ----------------------------------------------------------
  CHARACTER*30 :: elefil
  INTEGER :: le,iun,iunout,i
! output elems
  REAL(KIND=dkind),DIMENSION(6) :: elkep,elequ,elcar ![UD,days]
! ------------------------------------------------------------------
  rhs=2
!  write(*,*)'rhs=',rhs
!  write(*,*)'initial epoch=','boh'

! Output file
  CALL filnam('.','results_in','fla',elefil,le)
  CALL filopn(iun,elefil(1:le),'OLD')

! Output file
  CALL filnam('.','results_out','fla',elefil,le)
  CALL filopn(iunout,elefil(1:le),'UNKNOWN')

  write(*,*)'-----------------------------------'
  write(*,*)'  convert_order_results'
  write(*,*)'INPUTS : t[days/2pi], KEP[m,rad]'
  write(*,*)'-----------------------------------'

! Read and write on file
  DO i=1,tmax
     READ(iun,*,END=999) tevol,elem(1:6) ![m,rad,days/2pi]
     CALL convert_from_kep(tevol,elem,elkep,elequ,elcar) ![km,rad,days]
     write(iunout,*)tevol,elkep,elequ,elcar
  ENDDO
  write(*,*)'convert_order_results : do loop too short!',i
  STOP
999 CONTINUE

  CALL filclo(iun,' ')
  CALL filclo(iunout,' ')

END PROGRAM convert_order_results

! ********************************************************
!              ORBITAL ELEMENTS UNIT CONVERSION
!
! Conversion in KEP : a[km],angles[rad]
! Conversion in EQU : a[km],h=e*cos(om),k=e*sin(om),p=tan(i/2)*cos(Om)
!                         q=tan(i/2)*sin(Om),l=mean anomaly[rad]
! Conversion in CAR : x,y,x[km],dx,dy,dz[km/days]
!
! Time : conversion from [days/pig] to [days]
! ********************************************************
SUBROUTINE convert_from_kep(tvl,elem,elkep,elequ,elcar)
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE
  REAL(KIND=dkind),PARAMETER :: ys=3.6d3*24.d0
  REAL(KIND=dkind),PARAMETER :: UD=42164.1697748545d0
  REAL(KIND=dkind),PARAMETER :: teph0=0.d0
! -------------------- interface --------------------
  REAL(KIND=dkind),INTENT(INOUT) :: tvl
  REAL(KIND=dkind),INTENT(IN) :: elem(6)
  REAL(KIND=dkind),INTENT(OUT) :: elkep(6),elequ(6),elcar(6)
! ------------------ end interface ------------------
  TYPE(orbit_elem) :: ele
  INTEGER :: i,fail_flag
! ===================================================
  tvl=tvl*2.d0*pig        !from days/2pi to days

  elkep = elem         !KEP [m,rad]
  elkep(1)=elem(1)*UD  !KEP [km,rad]

  ele=undefined_orbit_elem
  ele%t=teph0+tvl !in MJD
  ele%coo='KEP'
  ele%coord(1:6)=elkep(1:6)
  CALL coo_cha(ele,'EQU',ele,fail_flag)
  IF(fail_flag.ge.5) THEN
     WRITE(*,*)'convert_elem, coo_cha: fail_flag=',fail_flag
     STOP
  ENDIF
  elequ(1:6)=ele%coord(1:6) ![km,rad]

  CALL coo_cha(ele,'CAR',ele,fail_flag)
  IF(fail_flag.ge.5) THEN
     WRITE(*,*)'convert_elem, coo_cha: fail_flag=',fail_flag
     STOP
  ENDIF
  elcar(1:6)=ele%coord(1:6) ![km,km/days]

ENDSUBROUTINE convert_from_kep
! ================================================

