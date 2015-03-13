! Chiara Tardioli (tardioli@mail.dm.unipi.it)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                           !!
!!                      FILTER II                            !!
!!                                                           !!
!!         Chiara Tardioli, May-Jun 2012, Namur              !!
!!                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ============================================================
! Contains
! filter_two : main subroutine, call orbit_distance
! orbit_distance   : compute the distance with sign (dmintil)
!                    at each node of the evolution
! --> A cross is excluded if dmintil does not change sign <--
! exchange_minimum : exchange local minima if needed
! exchange_value   : given values a,b exchange a with b
! write_elqu       : write the evolution of an EQU orbit on file
! write_equ_3      : write dist, EQU elements on screen
!
! OUTPUT
! dcros  : dcros(i,j)=1 if i and j pass filter II, otherwise it is 0
!           (declare in module filter_output)
! -----------------------------------------------------------
! Reference: Gronchi & Tommei (2007)
! ============================================================
MODULE filter2
  USE fund_const
  USE orbit_elements
  USE critical_points
  USE option_file
  USE read_interpol_ephem
  USE filter_output,only:hoots,dcros
  IMPLICIT NONE
  PRIVATE
! ---------------------------------------------------------
  PUBLIC :: filter_two
! ============================================================
CONTAINS
! ********************************************************
!                   FILTER II
! ********************************************************
  SUBROUTINE filter_two(iunda,orb1,orb2)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: iunda     !direct acces iunit
    INTEGER,INTENT(IN) :: orb1,orb2 !orbit indexes for doloop
! ------------------- end interface ------------------
    TYPE(orbit_elem),DIMENSION(norbx,tmax) :: elem !EQU orbital elements 1
    INTEGER :: pos
    LOGICAL :: check=.false.
! ============================================================
    IF(check)THEN
       write(*,*)norbx,tmax,nevol,hevol
       DO pos=1,nevol
          CALL read_evol(iunda,pos,elem(1,pos))
       ENDDO
       CALL write_elqu(0,elem(1,:),nevol)
    ENDIF
! ----------------------------------------------------------
! main call
    CALL local_orbit_distances(iunda,orb1,orb2)
! ----------------------------------------------------------
  END SUBROUTINE filter_two

! ********************************************************
!          LOCAL ORBIT DISTANCES WITH SIGN
! ********************************************************
  SUBROUTINE local_orbit_distances(iunda,orb1,orb2)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: iunda     !direct acces iunit
    INTEGER,INTENT(IN) :: orb1,orb2 !orbit indexes for doloop
! ------------------- end interface ------------------
    REAL(KIND=dkind),DIMENSION(tmax) :: tevol
! for dmintil_rms
    TYPE(orbit_elem) :: elem1 !EQU orbital elements 1
    TYPE(orbit_elem) :: elem2 !EQU orbital elements 2
    INTEGER :: nummin
    REAL(KIND=dkind),DIMENSION(nminx) :: dmintil
    REAL(KIND=dkind),DIMENSION(3,nminx) :: c1min,c2min
! results
    REAL(KIND=dkind),DIMENSION(2) :: dmint   !local orbit distances
    REAL(KIND=dkind),DIMENSION(3,2) :: c1,c2 !minimum points
! auxilia variables
    LOGICAL :: succ
    INTEGER :: i,j,k,n,pos
    REAL(KIND=dkind) :: aux
! ============================================================
! initialization
!  dcros=0.d0
! ----------------------------------------------------------
! time vector
    tevol(1)=0.d0
    DO k=2,hevol
       tevol(k)=tevol(k-1)+dt
    ENDDO
! ----------------------------------------------------------
! Loop on the orbits
    DO i=orb1,orb2-1
       IF(SUM(hoots(i,:)).eq.0.d0)THEN
!          write(*,*)'no crossings for',i
          GO TO 7 !next i
       ENDIF
       IF(mod(i,100).eq.0.d0)THEN
          write(*,*)i
       ENDIF
       CALL read_evol(iunda,(i-1)*nevol,elem1)
       DO j=i+1,orb2
!          IF(hoots(i,j).eq.0.d0)THEN
!             write(*,*)'no cross:',i,j
!             GO TO 5 !next j
!          ENDIF
          IF(hoots(i,j)-dcros(i,j).ne.1.d0)THEN
!             write(*,*)'no cross or cross already detected:',i,j
             GO TO 5 !next j
          ENDIF
!          write(*,*)i,j
          CALL read_evol(iunda,(j-1)*nevol,elem2)
! ----------------------------------------------------------
! first step and initializations
          CALL dmintil_rms(elem1,elem2,nummin,dmintil,c1min,c2min)
! ----------------------------------------------------------
! check L write results on screen
!             CALL write_equ_dist(elem1,elem2,dmint,c1,c2)
! ----------------------------------------------------------
! update variables
          dmint(1)=dmintil(1)
          dmint(2)=dmintil(2)
          c1(:,1:2)=c1min(:,1:2)
          c2(:,1:2)=c2min(:,1:2)
! ----------------------------------------------------------
! Loop on the evolution, dt is the time step, hevol is the number of steps
          DO k=2,hevol
             elem1%t=teph0+tevol(k) !in MJD
             elem2%t=elem1%t
! ----------------------------------------------------------
! linar interpolation to get the orbital elements at time dt*k 
             pos=(i-1)*nevol
             CALL interpolate_elem(iunda,pos,elem1%t,elem1%coord(1:6))
             pos=(j-1)*nevol
             CALL interpolate_elem(iunda,pos,elem2%t,elem2%coord(1:6))
! ----------------------------------------------------------
! MAIN PART
! compute the minima distances and the minima points
             CALL dmintil_rms(elem1,elem2,nummin,dmintil,c1min,c2min)
!             write(*,*)i,j,nummin,dmintil(1:2)
! check: too many minima for dmintil_rms
             IF(nummin.ge.5)THEN
                dcros(i,j)=1.d0
                write(*,*)i,j,nummin
                GOTO 5 !next j
             ENDIF
! check for minima exchange
! ATTENTION: the following 'if' do not works when minimal points collide
             CALL exchange_minimum(c1(:,1:2),c2(:,1:2),c1min(:,1:2), &
                  & c2min(:,1:2),succ)
             IF(succ)THEN
!                write(*,*)'******************************'
!                write(*,*)'exchange of minima'
!                write(*,*)'before :',dmintil(1),dmintil(2)
                CALL exchange_value(dmintil(1),dmintil(2))
!                write(*,*)'after :',dmintil(1),dmintil(2)
!                write(*,*)'******************************'
             ENDIF
! ----------------------------------------------------------
! check : write results on screen
!             CALL write_equ_dist(elem1,elem2,dmintil,c1min,c2min)
! ----------------------------------------------------------
! filter condition
             IF(dmint(1)*dmintil(1).le.0.d0)THEN
!                write(*,*)'cross 1',i,j,dmint(1),dmintil(1)
                dcros(i,j)=1.d0
                GOTO 5  !next j
             ELSEIF(dmint(2)*dmintil(2).le.0.d0)THEN
!                write(*,*)'cross 2',i,j,dmint(2),dmintil(2)
                dcros(i,j)=1.d0
                GOTO 5  !next j
             ENDIF
! ----------------------------------------------------------
! update variables
             dmint(1)=dmintil(1)
             dmint(2)=dmintil(2)
             c1(:,1:2)=c1min(:,1:2)
             c2(:,1:2)=c2min(:,1:2)
! ----------------------------------------------------------
          ENDDO
5      CONTINUE
       ENDDO
7     CONTINUE
    ENDDO
  END SUBROUTINE local_orbit_distances


! ********************************************************
!        CONTROL OF EXCHANGING OF LOCAL MINIMA 
! ********************************************************
  SUBROUTINE exchange_minimum(p1,p2,q1,q2,succ)
! --------------------- interface ---------------------
    REAL(KIND=dkind),DIMENSION(3,2),INTENT(INOUT) :: p1,p2,q1,q2
    LOGICAL,INTENT(OUT) :: succ
! ------------------- end interface -------------------
    LOGICAL :: check=.false.
    INTEGER :: i
! -----------------------------------------------------
    succ=.false.
!    write(*,*)p1(:,1:2),q1(:,1:2)!,p2(:,1),q2(:,1)

! ******* ANOTHER EXCHANGE CONDITION *************
!    IF((vsize(p1(:,1)-q1(:,1)).gt.vsize(p1(:,2)-q1(:,1)).and. &
!         & vsize(p2(:,1)-q2(:,1)).gt.vsize(p2(:,2)-q2(:,1))).or. &
!         & (vsize(p1(:,2)-q1(:,2)).gt.vsize(p1(:,1)-q1(:,2)).and. &
!         & vsize(p2(:,2)-q2(:,2)).gt.vsize(p2(:,1)-q2(:,2))))THEN
!       write(*,*)'minima exchange'
!    ENDIF
! ***********************************************

    IF(abs(q1(1,1)-p1(1,2)).lt.abs(q1(1,1)-p1(1,1)).AND. &
         & abs(q1(2,1)-p1(2,2)).lt.abs(q1(2,1)-p1(2,1)).AND. &
         & abs(q1(3,1)-p1(3,2)).lt.abs(q1(3,1)-p1(3,1)))THEN

       IF(check)THEN
          IF(abs(q2(1,1)-p2(1,2)).lt.abs(q2(1,1)-p2(1,1)).AND. &
               & abs(q2(2,1)-p2(2,2)).lt.abs(q2(2,1)-p2(2,1)).AND. &
               & abs(q2(3,1)-p2(3,2)).lt.abs(q2(3,1)-p2(3,1)))THEN
             write(*,*)'c2min(:,1) and c2min(:,2) exchange'
             write(*,*)'p2=',p2
             write(*,*)'q2=',q2
          ELSE
             write(*,*)'error exchange minima'
             write(*,*)'p2=',p2
             write(*,*)'q2=',q2
             STOP
          ENDIF
       ENDIF

!       write(*,*)'c1min(:,1) and c1min(:,2) exchange'
!       write(*,*)'q1=',q1
!       write(*,*)'q2=',q2
       DO i=1,3
          CALL exchange_value(q1(i,1),q1(i,2))
          CALL exchange_value(q2(i,1),q2(i,2))
       ENDDO
!       write(*,*)'q1=',q1
!       write(*,*)'q2=',q2
       succ=.true.
!       write(*,*)succ
    ENDIF
  END SUBROUTINE exchange_minimum

! ********************************************************
! EXCHANGE TWO REAL VALUES
! ********************************************************
  SUBROUTINE exchange_value(a,b)
! --------------------- interface ---------------------
    REAL(KIND=dkind),INTENT(INOUT) :: a,b
! ------------------- end interface -------------------
    REAL(KIND=dkind) :: aux
! -----------------------------------------------------
    aux=a
    a=b
    b=aux
  END SUBROUTINE exchange_value

! ********************************************************
!        WRITE EQU ORBIT EVOLUTION ON FILE
! ********************************************************
  SUBROUTINE write_elqu(iunit,elqu,n)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: iunit
    TYPE(orbit_elem),DIMENSION(tmax),INTENT(IN) :: elqu
    INTEGER,INTENT(IN) :: n
! ------------------- end interface -------------------
    INTEGER :: i
! -----------------------------------------------------
    write(*,*)iunit,elqu(1)
    DO i=1,n
       write(iunit,100)elqu(i)%t,elqu(i)%coord(1:5)
    ENDDO
100 FORMAT(f15.5,2x,f10.3,3(2x,e15.8),2x,f15.8)
  END SUBROUTINE write_elqu

! ********************************************************
!     WRITE RESULTS OF FILTER II ON SCREEN
! ********************************************************
  SUBROUTINE write_equ_dist(el1,el2,dmint,c1,c2)
! --------------------- interface ---------------------
    TYPE(orbit_elem),INTENT(IN) :: el1  !EQU orbital elements 1
    TYPE(orbit_elem),INTENT(IN) :: el2  !EQU orbital elements 2
    REAL(KIND=dkind),DIMENSION(nminx),INTENT(IN) :: dmint   !orbit distances
    REAL(KIND=dkind),DIMENSION(3,nminx),INTENT(IN) :: c1,c2 !minimum points
! ------------------- end interface -------------------
! =====================================================
    write(0,101)el1%t,dmint(1:2),el1%coord(1:5),el1%coord(1:5), & 
         & c1(:,1),c2(:,1),c1(:,2),c2(:,2)
101 FORMAT(f15.5,2(1x,f15.8),2(1x,f10.3,3(1x,e15.8),1x,f15.8),4(1x,f15.8))
  END SUBROUTINE write_equ_dist


END MODULE filter2
