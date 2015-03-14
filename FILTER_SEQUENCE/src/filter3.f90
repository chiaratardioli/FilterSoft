! Chiara Tardioli (tardioli@mail.dm.unipi.it)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                           !!
!!                        FILTER III                         !!
!!                                                           !!
!!         Chiara Tardioli, August 2012, Namur               !!
!!                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ============================================================
! Contains
! filter_three : main subroutine, call car_distance,
!                contains the filter condition
! --> A close encounter is excluded if 'd>D', D threshold  <--
! car_distance : compute the Keplerian distance d
!                at each node of the evolution
! write_equ_3  : write dist, EQU elements on screen
!
! OUTPUT
! dtime  : dtime(i,j)=1 if i and j pass filter III, otherwise it is 0
!           (declare in module filter_output)
! ============================================================
MODULE filter3
  USE fund_const
  USE orbit_elements
  USE critical_points
  USE option_file
  USE read_interpol_ephem
  USE filter_output,only:dcros,dtime
  IMPLICIT NONE
  PRIVATE
! ---------------------------------------------------------
  PUBLIC :: filter_three
! ============================================================
CONTAINS
! ********************************************************
! FILTER III : the time filter
! ********************************************************
  SUBROUTINE filter_three(iunda,orb1,orb2)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: iunda     !direct acces iunit
    INTEGER,INTENT(IN) :: orb1,orb2 !orbit indexes for doloop
! ------------------- end interface -------------------
    TYPE(orbit_elem) :: elem1   !EQU orbital elements 1
    TYPE(orbit_elem) :: elem2   !EQU orbital elements 2
    REAL(KIND=dkind),DIMENSION(tmax) :: tevol   !evolution vector
    REAL(KIND=dkind) :: dist    !Keplerian distance
    INTEGER :: i,j,k,pos
! =====================================================
! initialization
    dtime=0.d0
! ----------------------------------------------------------
! time vector
    tevol(1)=0.d0
    DO k=2,hevol
       tevol(k)=tevol(k-1)+dt
    ENDDO
! ----------------------------------------------------------
! Loop on the orbits
    DO i=orb1,orb2-1
!       write(*,*)i 
       IF(SUM(dcros(i,:)).eq.0.d0)THEN
!          write(*,*)'no crossings for',i
          GO TO 7 !next i
       ENDIF
       IF(mod(i,100).eq.0.d0)THEN
          write(*,*)i
       ENDIF
       CALL read_evol(iunda,(i-1)*nevol,elem1)
       DO j=i+1,orb2
          IF(dcros(i,j).eq.0.d0)THEN
!             write(*,*)'no cross:',i,j
             GO TO 5 !next j
          ENDIF
!          write(*,*)i,j
          CALL read_evol(iunda,(j-1)*nevol,elem2)
! ----------------------------------------------------------
! Loop on the evolution
          DO k=1,hevol
! ----------------------------------------------------------
! set the current time
             elem1%t=teph0+tevol(k) !in MJD
             elem2%t=elem1%t
! ----------------------------------------------------------
! interpolated orbit 1
             pos=(i-1)*nevol
             CALL interpolate_elem(iunda,pos,elem1%t,elem1%coord(1:6))
! interpolated orbit 2
             pos=(j-1)*nevol
             CALL interpolate_elem(iunda,pos,elem2%t,elem2%coord(1:6))
! ----------------------------------------------------------
! Keplerian distance funcion
             CALL car_distance(elem1,elem2,dist)
! ----------------------------------------------------------
! check : write results on file
!    CALL write_equ_3(elem1,elem2,dist)
! -----------------------------------------------
! filter condition
             IF(dist.lt.d_thres)THEN
                dtime(i,j)=1.d0
                GOTO 5 !next j
             ENDIF
! -----------------------------------------------
          ENDDO
5      CONTINUE
       ENDDO
7     CONTINUE
    ENDDO
  END SUBROUTINE filter_three

! ********************************************************
!      KEPLERIAN DISTANCE FUNCTION
! ********************************************************
  SUBROUTINE car_distance(elem1,elem2,dist)
! --------------------- interface ---------------------
    TYPE(orbit_elem),INTENT(IN) :: elem1   !EQU orbital elements 1
    TYPE(orbit_elem),INTENT(IN) :: elem2   !EQU orbital elements 2
    REAL(KIND=dkind),INTENT(OUT) :: dist   !Keplerian distance
! ------------------- end interface -------------------
    TYPE(orbit_elem) :: elcar1  !CAR orbital elements 1
    TYPE(orbit_elem) :: elcar2  !CAR orbital elements 2
    REAL(KIND=dkind),DIMENSION(3) :: car1,car2 !Cartesian positions
    REAL(KIND=dkind) :: dist2   !square Keplerian distance
    INTEGER :: fail_flag
! =====================================================
! CAR orbit 1
    CALL coo_cha(elem1,'CAR',elcar1,fail_flag)
    IF(fail_flag.ge.5)THEN
       WRITE(*,*)'time_filter, coo_cha: fail_flag=',fail_flag
       STOP
    ENDIF
    car1(1:3)=elcar1%coord(1:3)
! ----------------------------------------------------------
! CAR orbit 2
    CALL coo_cha(elem2,'CAR',elcar2,fail_flag)
    IF(fail_flag.ge.5)THEN
       WRITE(*,*)'time_filter, coo_cha: fail_flag=',fail_flag
       STOP
    ENDIF
    car2(1:3)=elcar2%coord(1:3)
! ----------------------------------------------------------
! Keplerian distance funcion
    dist2=(car1(1)-car2(1))**2+(car1(2)-car2(2))**2+(car1(3)-car2(3))**2
    dist=sqrt(dist2)
! ----------------------------------------------------------
! check : write results on file
!    CALL write_equ_3(elem1,elem2,dist)
! -----------------------------------------------
! filter condition
!    IF(dist.lt.d_thres)THEN
!       dtime(i,j)=1.d0
!!!!! new exit condition !!!!!!!
!    ENDIF
! -----------------------------------------------
  END SUBROUTINE car_distance


! ********************************************************
!     WRITE RESULTS OF FILTER III ON SCREEN
! ********************************************************
  SUBROUTINE write_equ_3(el1,el2,dist)
! --------------------- interface ---------------------
    TYPE(orbit_elem),INTENT(IN) :: el1  !EQU orbital elements 1
    TYPE(orbit_elem),INTENT(IN) :: el2  !EQU orbital elements 2
    REAL(KIND=dkind),INTENT(IN) :: dist !Keplerian distance
! ------------------- end interface -------------------
    write(0,100)el1%t,dist,el1%coord(1:5),el1%coord(1:5)
100 FORMAT(f15.5,1x,f15.8,2(1x,f10.3,3(1x,e15.8),1x,f15.8))
  END SUBROUTINE write_equ_3

END MODULE filter3
