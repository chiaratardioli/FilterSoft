! Chiara Tardioli (tardioli@mail.dm.unipi.it)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                           !!
!!              EPHEMERIS INTERPOLATION TABLE                !!
!!                                                           !!
!!           Chiara Tardioli, Aug 2012, Namur                !!
!!                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ============================================================
! Contains
! read_eph_table : read ephemeris interpolation table
! read_elkep     : read no. of evolution steps (nevol)
! read_evol      : read first orbit evolution from ephemeris table
! interpolate_elem : compute an orbit by linear interpolation
!
! out of date
! convert_elem   : to convert in EQU
! ============================================================
MODULE read_interpol_ephem
  USE fund_const
  USE orbit_elements
  USE option_file
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_eph_table,read_evol,interpolate_elem,check_eph_table

CONTAINS

! ************************************************************
!      READ EPHEMERIS INTERPOLATION TABLE
! ************************************************************
  SUBROUTINE read_eph_table(iunda,nvl)
! -------------------- interface --------------------
    INTEGER,INTENT(OUT) :: iunda !direct access iunit
    INTEGER,INTENT(OUT) :: nvl !no. of evolutions
! ------------------ end interface ------------------
! for direct access
    INTEGER :: length
    REAL(KIND=dkind) :: tevol,evolelem(6)
! for read_elkep
    CHARACTER*20 :: filename
    INTEGER :: n1 !no. of evolution for the orbital elements
!  INTEGER :: k !input data (1=nimstep,2=charles' order)
! ---------------------------------------------------
! ephemerides interpolaton table
    INQUIRE(IOLENGTH=length)tevol,evolelem(1:6)
    CALL fildiropn(iunda,'all_orbit_evolutions.a','OLD',length)

! check direct access
!  READ(iunda,REC=1)tevol,evolelem
!  write(*,*)tevol,evolelem

! Evolution span : I take this info from one particular evolution
    filename='results_kepl_1.dat'
    CALL read_elkep(filename,n1)
    nvl=n1
  END SUBROUTINE read_eph_table


! ************************************************************
!        READ NUMBER OF EVOLUTION STEPS
! ************************************************************
  SUBROUTINE read_elkep(filename,nvl)
!    USE fund_const
!  USE option_file,only:tmax
!  IMPLICIT NONE
! -------------------- interface --------------------
    CHARACTER*20,INTENT(IN) :: filename
    INTEGER,INTENT(OUT) :: nvl
! ------------------ end interface ------------------
    REAL(KIND=dkind) :: tvl,evolel(6)
    INTEGER :: iunorb,i
! ---------------------------------------------------
    CALL filopn(iunorb,filename,'OLD')
    DO i=1,tmax
       READ(iunorb,*,END=998) tvl,evolel(1:6)
! check
!     CALL convert_elem(2,tvl,evolel) !elements are EQU in output
!     WRITE(*,*)tvl,evolel(1:6)
    ENDDO
    write(*,*)'read_elkep : do loop too short!',i
    STOP
998 CONTINUE
    CALL filclo(iunorb,' ')
    nvl=i-1
  ENDSUBROUTINE read_elkep
! ============================================================

! ********************************************************
! READ FIRST ORBIT EVOLUTION FROM THE EPHEMERIS TABLE
! ********************************************************
  SUBROUTINE read_evol(iunda,pos,elem)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: iunda,pos
    TYPE(orbit_elem) :: elem
! ------------------- end interface ------------------
    REAL(KIND=dkind) :: t
! ---------------------------------------------------------
    elem=undefined_orbit_elem
    elem%coo='EQU'
    READ(iunda,REC=pos+1)t,elem%coord(1:6)
    elem%t=teph0+t !MJD
  END SUBROUTINE read_evol
! ============================================================

! **************************************************
!  COMPUTE ORBITAL ELEMENTS BY LINEAR INTERPOLATION
! **************************************************
  SUBROUTINE interpolate_elem(iunda,pos,time,elout)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: iunda,pos !for direct access
    REAL(KIND=dkind),INTENT(IN) :: time !interpolation time (in MJD)
    REAL(KIND=dkind),DIMENSION(6),INTENT(OUT) :: elout
! ------------------- end interface ------------------
    REAL(KIND=dkind) :: tevol(2),princ
    REAL(KIND=dkind) :: t1,t2 !extrema of the interpolation interval
    REAL(KIND=dkind),DIMENSION(6) :: el1,el2 !orbits at t1,t2
    REAL(KIND=dkind) :: lambda !coeff of convex combination
    REAL(KIND=dkind) :: intstep ! orbit integration time step (in sec)
    INTEGER :: nstep ! no. of time steps to reach t1
! check
    LOGICAL :: err=.false.,check=.true.
    INTEGER :: i
! =========================================================
    READ(iunda,REC=pos+1)tevol(1)
    READ(iunda,REC=pos+2)tevol(2)
! intstep is not a parameter, allowing to change its value
    intstep=tevol(2)-tevol(1) !in days
! *** assuming that time >= teph0 ***
    IF(time.lt.teph0)THEN
       write(*,*)'interpolate_elem : negative time'
    ENDIF
    nstep= FLOOR((time-(teph0+tevol(1)))/intstep)

!  *** t0+tevol(nstep+1) <= time <= t0+tevol(nstep+2) ***
    READ(iunda,REC=pos+nstep+1)t1,el1(1:6)
    READ(iunda,REC=pos+nstep+2)t2,el2(1:6)
    lambda=(time-teph0-t1)/(intstep) !interpolating parameter,t1->t1+teph0

!  IF(lambda.lt.-1.d-10.OR.lambda.gt.1.d0+1.d-10)THEN
    IF(lambda.lt.-1.d-10.OR.lambda.gt.1.d0)THEN
       WRITE(*,*) 'error! lambda = ',lambda
       WRITE(*,*) 'time = ',time,'nstep = ',nstep
       WRITE(*,*) intstep
       err=.true.
    ENDIF

    IF(el2(6).lt.el1(6)) THEN
       el2(6) = el1(6) + dpig
    ENDIF
    elout = (1.d0-lambda) * el1 + lambda * el2
    elout(6) = princ(elout(6))

    IF(check.AND.err)THEN
       write(*,*)'POS',pos+nstep+1,pos+nstep+2
       write(*,*)el1
       write(*,*)el2
       write(*,*)'intstep',intstep
       write(*,*)'nstep',nstep
       write(*,*)'t1=',t1,'t2=',t2
       write(*,*)'time=',time-teph0,'epoch=',time
       write(*,*)'lambda',lambda
       DO i=1,6
          write(*,*)i,el1(i),el2(i),elout(i),elout(i)-el1(i),elout(i)-el2(i)
       ENDDO
       STOP
    ENDIF

  END SUBROUTINE interpolate_elem
! ============================================================


! ********************************************************
!              ORBITAL ELEMENTS UNIT CONVERSION
! k=1 NIMSTEP (Nicolas Delsate's program), k=2 ORDER (Charles Hubaux's program)
!
! First, convert to KEP : a[km],angles[rad],time[h]
! Then, convert to EQU  : a,h=e*cos(om),k=e*sin(om),p=tan(i/2)*cos(Om)
!                         q=tan(i/2)*sin(Om),l=mean anomaly
! ********************************************************
  SUBROUTINE convert_elem(k,tvl,elem)
!    USE fund_conast
!  USE orbit_elements
!  USE option_file,only:tmax,teph0
!  IMPLICIT NONE
    REAL(KIND=dkind),PARAMETER :: ys=3.6d3*24.d0
    REAL(KIND=dkind),PARAMETER :: UD=42164.1697748545d0
! -------------------- interface --------------------
    INTEGER,INTENT(IN) :: k
    REAL(KIND=dkind),INTENT(INOUT) :: tvl,elem(6)
! ------------------ end interface ------------------
    TYPE(orbit_elem) :: ele
    INTEGER :: i,fail_flag
! ===================================================
    IF(k.eq.1.OR.k.eq.3)THEN
       tvl=tvl/ys              !from sec to days
       elem(1)=elem(1)*1.d-3   !from meters to km
    ELSEIF(k.eq.2)THEN
       tvl=tvl*2.d0*pig        !from days/2pi to days
       elem(1)=elem(1)*UD      !to km
    ELSE
       write(*,*)'convert_elem : k too big !!!',k
    ENDIF
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
! ===================================================

! ================================================
! CHECK EPHEMERIS INTERPOLATION TABLE
! ================================================
  subroutine check_eph_table(iunda,pos)
    USE option_file,ONLY: tmax,dt,hevol,nevol
! -------------------- interface --------------------
    INTEGER,INTENT(IN) :: iunda,pos
! ------------------ end interface ------------------
    TYPE(orbit_elem) :: equ
    REAL(KIND=dkind),DIMENSION(tmax) :: tevol
    INTEGER :: k
! ===================================================

    write(*,*)dt,hevol,nevol

    READ(iunda,REC=pos+1)tevol(1),equ%coord
    equ%t=teph0+tevol(1)

    READ(iunda,REC=pos+2)tevol(2),equ%coord
    write(*,*)tevol(2)-tevol(1)

    do k=2,hevol
       tevol(k)=tevol(k-1)+dt
       equ%t=teph0+tevol(k) !in MJD
       CALL interpolate_elem(iunda,pos,equ%t,equ%coord(1:6))
       write(1,*)tevol(k),equ%coord(1:6)
    enddo

    do k=1,nevol
       READ(iunda,REC=pos+k)tevol(k),equ%coord
       write(2,*)tevol(k),equ%coord
    enddo

  END subroutine check_eph_table

END MODULE read_interpol_ephem
