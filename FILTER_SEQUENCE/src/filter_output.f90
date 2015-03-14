!
! Chiara Tardioli, Aug 2012, Namur
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                           !!
!!                                                           !!
!!           M O D U L E   F O R   O U T P U T S             !!
!!                                                           !!
!!                                                           !!
!!                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ============================================================
! Contains
! hoots : output of filter I
! dcros : output of filter II
! dtime : output of filter III
!
! write_radii     : write on file radial distance evolution
! write_results   : write on screen the outup sum for each filter
! diff_filter_1_2 : difference in the output sums of filter I and II
! ============================================================
MODULE filter_output
  USE fund_const
  USE option_file
  IMPLICIT NONE
  PRIVATE
  INTEGER,DIMENSION(norbx,norbx) :: hoots ! hoots(i,j)=1 if i and j pass
                                          ! filter 1, otherwise it is 0
  INTEGER,DIMENSION(norbx,norbx) :: dcros ! dcros(i,j)=1 if i and j pass
                                          ! filter 2, otherwise it is 0
  INTEGER,DIMENSION(norbx,norbx) :: dtime ! dtime(i,j)=1 if i and j pass
                                          ! filter 3, otherwise it is 0

  PUBLIC :: hoots,dcros,dtime
  PUBLIC :: write_results,diff_filter_1_2,write_radii

CONTAINS
! ********************************************************
!         WRITE ON SCREEN RESULTS OF THE 3 FILTERS
! ********************************************************
  SUBROUTINE write_results(fil)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: fil
! ------------------- end interface ------------------
    INTEGER,DIMENSION(norbx,norbx) :: m
    INTEGER :: sum1,sum2,sum3
! ============================================================
    IF(fil.eq.1)THEN
       sum1=SUM(hoots,hoots.gt.0)
       write(*,*)'Filetr I (geo dist)'
       write(*,*)'evolution step =',dt
       write(*,*)'no. of evolutions =',hevol
       write(*,*)'d_thres=',d_thres
       write(*,*)'cross =',sum1
       write(*,*)'****************************************'
    ELSEIF(fil.eq.2)THEN
       sum2=SUM(dcros)
       write(*,*)'Filetr II (orbit distance)'
       write(*,*)'evolution step =',dt
       write(*,*)'no. of evolutions =',hevol
       write(*,*)'cross    =',sum2
       write(*,*)'****************************************'
    ELSEIF(fil.eq.3)THEN
       sum3=SUM(dtime)
       write(*,*)'Filetr III (time filter)'
       write(*,*)'evolution step =',dt
       write(*,*)'no. of evolutions =',hevol
       write(*,*)'d_thres=',d_thres
       write(*,*)'cross    =',sum3
       write(*,*)'****************************************'
    ELSE
       write(*,*)'error write_results: fil =',fil
    ENDIF

  END SUBROUTINE write_results
! ============================================================

! ********************************************************
!       DIFFERENCE IN THE OUTPUTS OF FILTER I AND II
! ********************************************************
  SUBROUTINE diff_filter_1_2(orb1,orb2)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: orb1,orb2
! ------------------- end interface ------------------
    INTEGER,DIMENSION(norbx,norbx) :: m
    INTEGER :: couples,sumf1,sumf2,sumf12
    INTEGER :: i,j,count
! ============================================================
    write(*,*)'Filter 1 VS Filter 2'
    m=hoots-dcros
    couples=(orb2-orb1+1)*(orb2-orb1)/2
    sumf1=SUM(m,m.eq.1.d0)
    sumf2=SUM(m,m.eq.-1.d0)
    m=hoots+dcros
    sumf12=SUM(m,m.eq.2.d0)/2
    write(*,*)'no fI si fII :',sumf2
    write(*,*)'si fI no fII :',sumf1
    write(*,*)'si fI si fII :',sumf12
    write(*,*)'no fI no fII :',couples-sumf1-sumf2-sumf12
    write(*,*)'****************************************'

!    count=0
!    DO i=orb1,orb2-1
!       DO j=i+1,orb2
!          IF(m(i,j).eq.0.d0)THEN
!             count=count+1
!             write(*,*)i,j
!          ENDIF
!       ENDDO
!    ENDDO
!    write(*,*)'count =',count,orb1,orb2
  END SUBROUTINE diff_filter_1_2
! ============================================================

! ********************************************************
!        WRITE GEOCENTRIC DISTANCE ON FILE
! ********************************************************
  SUBROUTINE write_radii(iun,tevol,vecr,rmin,rmax)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: iun
    REAL(KIND=dkind),DIMENSION(tmax),INTENT(IN) :: tevol,vecr
    REAL(KIND=dkind),INTENT(IN) :: rmin,rmax
! ------------------- end interface -------------------
    INTEGER :: k,m
! -----------------------------------------------------
    write(*,*)'write on iunit',iun
    write(iun,*)rmin,rmax
    DO k=1,hevol
       write(iun,100)tevol(k),vecr(k)
    ENDDO
100 FORMAT(f10.3,1x,f15.5)
  END SUBROUTINE write_radii


END MODULE filter_output
