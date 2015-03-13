! Chiara Tardioli (tardioli@mail.dm.unipi.it)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                           !!
!!                       FILTER I                            !!
!!                                                           !!
!!           Chiara Tardioli, May 2012, Namur                !!
!!                                                           !!
!! Hoots  (1984) : apogee-perigee computed at time knots     !!
!! Milani (2005) : max-min of the geocentric distance        !!
!!                 computed by interpolation                 !!
!!                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Revised on Jun 13, 2014 CT
! ============================================================
! Contains
! filter_one : main subroutine, call hoot_1
! hoots_1    : contains the filter condition
!              compare r_max and r_min between a pair of orbits
!              by Hoots' filter1: q=min(r1_min,r2_min), Q=max(r1_max,r2_max)
!-> if |q-Q|>D (D>0 threshold) then no orbit crossing may occur <-
! apogee_perigee : SET min(r) and max(r), it will be useful for orbit uncertainty
! maxmin_geodist : cointains the main part; compute the actual max(r) and min(r) 
!                  by a regula-falsi and interpolation
! geodist    : compute the geocentric distance r(t) (suggestions on uncertainty)
! minmax     : contains spline interpolation; compute max and min of r(t)
!              on the biggest monotonal interval
! regula_falsi : regula_falsi algorithm
! evalspline   : evaluate f(x)
!
! HOOTS ET AL. FILTER
! maxmin_apoperi : compute the maximum apogee and minimum perigee
! apoperi        : compute apogee and perigee for fixed EQU elements
!
! OUTPUT 
! hoots   : hoots(i,j)=1 if i and j pass filter I, otherwise it is 0
!           (declare in module filter_output)
! ------------------------------------------------------------
! Reference: Hoots and Roerich (1984)
! ============================================================
MODULE filter1
  USE fund_const
  USE orbit_elements
  USE option_file
  USE read_interpol_ephem
  USE filter_output,only:hoots,write_radii
  IMPLICIT NONE
  PRIVATE
! ---------------------------------------------------------
  PUBLIC :: filter_one
! ============================================================
CONTAINS
! ********************************************************
! FILTER I : the geocentric distance
! output NxN matrix : hoots
! ********************************************************
  SUBROUTINE filter_one(iunda,orb1,orb2,check_apo_peri)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: iunda
    INTEGER,INTENT(IN) :: orb1,orb2 !first and last orbit index for doloop
    LOGICAL,INTENT(IN) :: check_apo_peri
! ------------------- end interface ------------------
    INTEGER :: pos   !index for reading binary file
    REAL(KIND=dkind),DIMENSION(norbx) :: r_min,r_max
    INTEGER :: i,j
    INTEGER,DIMENSION(norbx,norbx) :: m !auxiliar matrix
    INTEGER :: s !for results
! =======================================================
! Compute min,max of geocentric distance r(t)
! HINT :  r(t) change with the anomalies
! initializations
    r_min(1:orb2)=0.d0
    r_max(1:orb2)=0.d0
! loop on orbits
    DO j=orb1,orb2
       pos=(j-1)*nevol
       CALL maxmin_geodist(iunda,pos,r_min(j),r_max(j))
!       write(*,*)j,'***', r_min(j),r_max(j),'***'
    ENDDO
! ----------------------------------------------------
! main call
    CALL hoots_1(orb1,orb2,r_min(1:orb2),r_max(1:orb2))
! ----------------------------------------------------
! ----------------------------------------------------
! check : coparison with Hoot et al. filter
    IF(check_apo_peri)THEN
! save Filter I output
       m=hoots
! Compute min,max of perigee, apogee
       r_min(1:orb2)=0.d0
       r_max(1:orb2)=0.d0
       DO j=orb1,orb2
          pos=(j-1)*nevol
          CALL maxmin_apoperi(iunda,pos,r_min(j),r_max(j))
!     write(*,*)j,'***', r_min(j),r_max(j),'***'
       ENDDO
       CALL hoots_1(orb1,orb2,r_min(1:orb2),r_max(1:orb2))
! results
       s=SUM(hoots,hoots.gt.0)
       write(*,*)'********* filetr 1 (apo/peri) **********'
       write(*,*)'cross =',s
       write(*,*)'****************************************'
! ----------------------------------------------------
       write(*,*) 'Diff < 2'
       DO i=orb1,orb2-1
          DO j=i+1,orb2
             IF(m(i,j)-hoots(i,j).eq.1.AND.abs(i-j).lt.2)THEN
                write(*,*)i,j,m(i,j),hoots(i,j)
             ELSEIF(m(i,j)-hoots(i,j).eq.-1)THEN
                write(*,*)'-1',i,j
             ENDIF
          ENDDO
       ENDDO
! ----------------------------------------------------
! restore Filter I output
       hoots=m
! ----------------------------------------------------
    ENDIF
! end check
! ---------------------------------------------------
  END SUBROUTINE filter_one


! ********************************************************
! main part of FILTER I, output : hoots matrix
! ********************************************************
  SUBROUTINE hoots_1(orb1,orb2,rmin,rmax)
! --------------------- interface ---------------------
    INTEGER,INTENT(IN) :: orb1,orb2         !first and last orbit index
    REAL(KIND=dkind),DIMENSION(norbx),INTENT(IN) :: rmin,rmax
!    REAL(KIND=dkind),DIMENSION(norbx),INTENT(IN),OPTIONAL :: sigmin,sigmax
! ------------------- end interface ------------------
    REAL(KIND=dkind) :: apo1,peri1,apo2,peri2,qmax,Qmin
    INTEGER :: i,j,ntot
! ============================================================
! initialization
    hoots=0
! ----------------------------------------------------
! Loop on orbits
    DO i=orb1,orb2-1
! orbit 1 : max(r),min(r) -> apo1,peri1 , without uncertainty
       CALL apogee_perigee(rmin(i),rmax(i),apo1,peri1) 
       DO j=i+1,orb2
! orbit 2 : max(r),min(r) -> apo2,peri2 , without uncertainty
          CALL apogee_perigee(rmin(j),rmax(j),apo2,peri2)
! ----------------------------------------------------
! filter condition
          qmax=max(peri1,peri2)
          Qmin=min(apo1,apo2)
          IF(qmax-Qmin.le.d_thres)THEN
!             write(*,*)'orbit crossing : q-Q=',qmax-Qmin,'<',d_thres,'km'
             hoots(i,j)=1
          ELSE
!             write(*,*)i,j,': q=',qmax,'Q=',Qmin
!             write(*,*)'no orbit crossing :',qmax-Qmin,'>',d_thres,'km'
          ENDIF
! ----------------------------------------------------
       ENDDO
    ENDDO
  END SUBROUTINE hoots_1


! ********************************************************
! SET min(r) and max(r) in variables apo,peri
! it will be useful for orbit uncertainty
! ********************************************************
  SUBROUTINE apogee_perigee(minr,maxr,apo,peri) !,minsig,maxsig)
! --------------------- interface ---------------------
    REAL(KIND=dkind),INTENT(IN) :: minr,maxr
    REAL(KIND=dkind),INTENT(OUT) :: apo,peri
!    REAL(KIND=dkind),INTENT(IN),OPTIONAL :: minsig,maxsig
! ------------------- end interface ------------------
! ============================================================
!    IF(PRESENT(minsig))THEN
!       IF(PRESENT(maxsig))THEN
!          apo =maxr+maxsig
!          peri=minr-minsig
!          RETURN
!       ELSE
!          write(*,*)'uncertainty of r_max not present: sigmax'
!          STOP
!       ENDIF
!    ENDIF
!    IF(PRESENT(maxsig))THEN
!       write(*,*)'uncertainty of r_min not present: sigmin'
!       STOP
!    ENDIF
!    write(*,*)'no uncertainties for r_min,r_max'
    apo = maxr
    peri= minr
  END SUBROUTINE apogee_perigee


! ************************************************
! COMPUTATATION OF ACTUAL min(r) AND min(r), r geocentric distance
! throught the zero of its derivative dr(t)
! ************************************************
  SUBROUTINE maxmin_geodist(iunda,pos,rmin,rmax)
! ---------------- interface ------------------------
    INTEGER,INTENT(IN) :: iunda !direct access iunit
    INTEGER,INTENT(IN) :: pos   !integer for position of next 
                                !orbital element evolution
    REAL(KIND=dkind),INTENT(OUT) :: rmin,rmax  !,sigmin,sigmax
! ---------------- end interface -------------------
    REAL(KIND=dkind),DIMENSION(tmax) :: tevol !days from teph0
    TYPE(orbit_elem) :: elequ
    REAL(KIND=dkind) :: r,dr
    REAL(KIND=dkind),DIMENSION(tmax) :: vecr,vecdr
    REAL(KIND=dkind),DIMENSION(tmax) :: maxr,minr
    REAL(KIND=dkind) :: tcrit
    INTEGER :: nmin,nmax,ini
    CHARACTER*3 :: tipo !'min' or 'max'
    INTEGER :: kaux(1),kmin,kmax
    INTEGER :: k !loop index
! =======================================================
! first step
    elequ=undefined_orbit_elem
    elequ%coo='EQU'
    READ(iunda,REC=pos+1)tevol(1),elequ%coord(1:6) !check:pos=k*nevol=>t=0.d0
    elequ%t=teph0+tevol(1) !in MJD
    CALL geodist(elequ,r,dr)
!    write(*,*)'          1   r=',r,'dr=',dr
!       write(*,*)pos/400+1
!    write(pos/400+51,100)elequ%t-teph0,elequ%coord(1:5),r
!100 FORMAT(f15.5,2x,f10.3,3(2x,e15.8),2x,f15.8,2x,f15.5)
! ---------------------------------------------------
! initializations
    ini  = 1  !initial index for the spline curve of dr
    nmin = 1 !counting no. of minima
    nmax = 1 !counting no. of maxima
    maxr = r
    minr = r 
    vecdr(1) = dr
    vecr(1)  = r !for output
! ---------------------------------------------------
! Loop on the number of evolutions
    DO k=2,hevol
! ---------------------------------------------------
! set evolution time
       tevol(k)=tevol(k-1)+dt
       elequ%t=teph0+tevol(k) !in MJD
!       write(*,*)elequ%t
! ---------------------------------------------------
! EQU orbit interpolated at time tevol(k)
       CALL interpolate_elem(iunda,pos,elequ%t,elequ%coord(1:6))
! ---------------------------------------------------
! compute geocentric distance
       CALL geodist(elequ,r,dr)
!       write(*,*)k,'  r=',r,'dr=',dr
!       write(pos/400+51,100)elequ%t-teph0,elequ%coord(1:5),r
       vecdr(k)=dr
       vecr(k)=r
! ---------------------------------------------------
! check for STATIONARY POINTS dr(t)=0 in the interval [ t(ini),t(k) )
       IF(vecdr(k)*vecdr(k-1).le.0.d0)THEN
!          write(*,*)h,'  t1=',tevol(k-1),'t2=',tevol(k)
!          write(*,*)'dr1=',vecdr(k-1),'dr2=',vecdr(k)
          CALL minmax(k-ini+1,tevol(ini:k),vecdr(ini:k),tcrit,tipo)
!          write(*,*)'minmax :',tcrit,tipo
! compute the actual geocentric distance at time tcrit
          CALL interpolate_elem(iunda,pos,teph0+tcrit,elequ%coord(1:6))
          CALL geodist(elequ,r,dr)
!          write(*,*)'at tcrit : r=',r,'dr=',dr
! update
          IF(tipo.eq.'min')THEN
             nmin=nmin+1
             minr(nmin)=r
!             write(*,*)'min',nmin,r
          ELSEIF(tipo.eq.'max')THEN
             nmax=nmax+1
             maxr(nmax)=r
!             write(*,*)'max',nmax,r
          ELSE
             write(*,*)'minmax: saddle point'
          ENDIF
          ini=k
       ENDIF
! ---------------------------------------------------
    ENDDO
! Values of r(t) at the boundary
    nmin = nmin+1
    minr(nmin) = r
    nmax = nmax+1
    maxr(nmax) = r
!    write(*,*)k,'  r=',r,'dr=',dr
! End loop on the evolution
! ---------------------------------------------------
! Minimizzo i "pericentri" e massimizzo gli "apocentri"
    kaux=MINLOC(minr(1:nmin))
    kmin=kaux(1)
    rmin=minr(kmin)

    kaux=MAXLOC(maxr(1:nmax))
    kmax=kaux(1)
    rmax=maxr(kmax)
! ---------------------------------------------------
! write partial results for orb1=iun1,orb2=iun2 (see option_file)
!    IF(iun1.eq.INT(pos/nevol)+1)THEN
!       CALL write_radii(iun1,tevol,vecr,rmin,rmax)
!    ELSEIF(iun2.eq.INT(pos/nevol)+1)THEN
!       CALL write_radii(iun2,tevol,vecr,rmin,rmax)
!    ENDIF
! ---------------------------------------------------    
  ENDSUBROUTINE maxmin_geodist

! ---------------------------------------------------------
! **********************************************
! GEOCENTRIC DISTANCE (WITHOUT UNCERTAINTY)
! **********************************************
!  SUBROUTINE geodist(ele,r,dr,unckep,uncr)
  SUBROUTINE geodist(ele,r,dr)
! --------------------- interface ---------------------
    TYPE(orbit_elem),INTENT(IN) :: ele
    REAL(KIND=dkind),INTENT(OUT) :: r,dr
!    REAL(KIND=dkind),DIMENSION(5,5),INTENT(IN),OPTIONAL :: unckep
!    REAL(KIND=dkind),INTENT(OUT),OPTIONAL :: uncr
! ------------------- end interface ------------------
    TYPE(orbit_elem) :: elcar
!    REAL(KIND=dkind),DIMENSION(6,6) :: dcar_ele
!    REAL(KIND=dkind),DIMENSION(3,1) :: dr_dvecr
!    REAL(KIND=dkind),DIMENSION(3,3) :: uncvecr
!    REAL(KIND=dkind),DIMENSION(1,1) :: uncr_aux
    REAL(KIND=dkind) :: vsize
    INTEGER :: fail_flag
! for check_der
    REAL(KIND=dkind) :: incr,r_incr,vecr_incr(3)
! for check
    LOGICAL :: check_der=.false.
! =========================================================
!  write(*,*)'input geodist:'
!  write(*,*)ele,unckep
!  IF(PRESENT(unckep).AND.PRESENT(uncr))THEN
!     write(*,*)'unckep presente'
!  ELSE
!     write(*,*)'unckep non presente'
!  ENDIF

!    IF(PRESENT(unckep).AND.PRESENT(uncr))THEN
!     write(*,*)'qui unckep',unckep
!       CALL coo_cha(ele,'CAR',elcar,fail_flag,dcar_ele)
!       r=vsize(elcar%coord(1:3))
!       dr_dvecr(1:3,1)=elcar%coord(1:3)/r
!     write(*,*)dr_dvecr
!       uncvecr=MATMUL(TRANSPOSE(dcar_ele(1:5,1:3)),MATMUL(unckep(1:5,1:5), &
!            & dcar_ele(1:5,1:3)))
!       uncr_aux=MATMUL(TRANSPOSE(dr_dvecr),MATMUL(uncvecr,dr_dvecr))
!       uncr=uncr_aux(1,1)
!     write(*,*)'unc=',uncr
!    ELSE
! -----------------------------------------------------
    CALL coo_cha(ele,'CAR',elcar,fail_flag)
    IF(fail_flag.gt.5)THEN
       write(*,*)'coo_cha : fail_flag',fail_flag
       STOP
    ENDIF
    r=vsize(elcar%coord(1:3))
!    ENDIF
!  write(*,*)'car',elcar%coord(4:6)
! -----------------------------------------------------
!  dr=(elcar%coord(4)+elcar%coord(5)+elcar%coord(6))/r
    dr=DOT_PRODUCT(elcar%coord(1:3),elcar%coord(4:6))/r
!  write(*,*)'r=',r,'dr=',dr
! -----------------------------------------------------
! CHECK r-derivative
    IF(check_der)THEN
       incr=1.d-6
       vecr_incr=elcar%coord(1:3)
       vecr_incr(1)=vecr_incr(1)+incr
       vecr_incr(2)=vecr_incr(2)+incr
       vecr_incr(3)=vecr_incr(3)+incr
       r_incr=vsize(vecr_incr(1:3))
       write(*,*)'***********************************'
       write(*,*)'incremental ratio', (r_incr-r)/incr
       write(*,*)'derivative       ', &
            & (elcar%coord(1)+elcar%coord(2)+elcar%coord(3))/r
       write(*,*)'***********************************'
       STOP
    ENDIF
! -----------------------------------------------------
  END SUBROUTINE geodist


! ********************************************************
! APPROXIMATION OF max(r) and min(r) BY splines
! call regula_falsi
! ********************************************************
  SUBROUTINE minmax(n,te,dre,tcrit,tipo)
! ---------------- interface ----------------
    INTEGER,INTENT(IN) :: n
    REAL(KIND=dkind),DIMENSION(tmax),INTENT(IN) :: te
    REAL(KIND=dkind),DIMENSION(tmax),INTENT(IN) :: dre
! --------------- end interface -------------
    REAL(KIND=dkind),INTENT(OUT) :: tcrit
    CHARACTER*3 :: tipo !'min' or 'max'
! -------------------------------------------
    INTEGER :: i,j,h
! for spline
    REAL(KIND=dkind),DIMENSION(n) :: x,y,b,c,d,dmin_rms,t,dr
! for regula falsi
    INTEGER :: succ !0=ok, 1=zero in un nodo, 2=errore
! ============================================================
    t=te(1:n)
    dr=dre(1:n)
    CALL spline(n,t,dr,b,c,d)
    i=n-1
    CALL regula_falsi(t(i),dr(i),b(i),c(i),d(i),t(i),t(i+1),tcrit,succ)
    IF(tcrit.lt.t(i).OR.tcrit.gt.t(i+1))THEN
       write(*,*)'ERROR!! check_cross : tcrit not in [t(i),t(i+1)]'
       write(*,*)tcrit,t(i),t(i+1)
       STOP
    ENDIF
    IF((succ.eq.1.AND.dr(i).eq.0d0).OR.succ.eq.2)THEN
       write(*,*)'zero at first knot'
!    ELSEIF(succ.eq.2)THEN
!       write(*,*)'double crossing at knots, take the first'
    ELSEIF(succ.gt.2)THEN
       write(*,*)'check_cross : error in regula falsi'
       STOP
    ENDIF
    IF(dr(i+1).gt.0.d0)THEN
       tipo='min'
    ELSE
       tipo='max'
    ENDIF
5   CONTINUE
!    write(*,*)'*****************************'
!    write(*,*)'minmax'
!    write(*,*)n,t(n-1:n)
!    write(*,*)dr(n-1:n)
!    write(*,*)tipo,tcrit
!    write(*,*)'*****************************'
  END SUBROUTINE minmax


! ******************************************************************
! REGULA FALSI ALGORITHM (order 3)
! check zero for the function:
! s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
! ******************************************************************
  SUBROUTINE regula_falsi(xi,yi,bi,ci,di,xa,xb,x0,succ)
    INTEGER,PARAMETER :: itmax=10000
    REAL(KIND=dkind),PARAMETER :: eps=1.d-15
! ----------------- interface ------------------------
    REAL(KIND=dkind),INTENT(IN) :: xi,yi,bi,ci,di,xa,xb
    REAL(KIND=dkind),INTENT(OUT) :: x0
    INTEGER,INTENT(OUT) :: succ !0=ok, 1=zero in un nodo, 5=errore
! --------------- end interface ----------------------
    REAL(KIND=dkind) :: x,c,x1,fx,fc,fx1
    INTEGER :: k
! =========================================================
    succ=0
    x=xa
    c=xb
!  write(*,*)'x0=',x,'c0=',c

    DO k=1,itmax
!     fx = yi + bi*(x-xi) + ci*(x-xi)**2 + di*(x-xi)**3
       CALL evalspline(bi,ci,di,xi,yi,x,fx)
!     fc = yi + bi*(c-xi) + ci*(c-xi)**2 + di*(c-xi)**3
       CALL evalspline(bi,ci,di,xi,yi,c,fc)
       IF(fx*fc.ge.0d0)THEN
          IF(fx.eq.0.d0.and.fc.ne.0.d0)THEN
             succ=1
             GOTO 5
          ELSEIF(fc.eq.0.d0.and.fx.ne.0.d0)THEN
             succ=1
             x=c
             fx=fc
             GOTO 5
          ELSEIF(fx.eq.0.d0.and.fc.eq.0.d0)THEN
             succ=2
             write(*,*)'two zero???? x=',x,'c=',c
             GOTO 5
          ELSE
             succ=5
             write(*,*)'ERROR!!! f(x)*f(c)>0 :'&
                  &'f(',x,')=',fx,'f(',c,')=',fc,'it=',k
             x=-1000.d0
             GOTO 5
          ENDIF
       ENDIF
       IF(ABS(fx).lt.eps)THEN
!          write(*,*)'=',x-c,fx
          GOTO 5
       ENDIF
       x1=x-(fx*(x-c)/(fx-fc))
       
!     fx1= yi + bi*(x1-xi) + ci*(x1-xi)**2 + di*(x1-xi)**3 
       CALL evalspline(bi,ci,di,xi,yi,x1,fx1)
       IF(fx1*fc.gt.0.d0)THEN
          c=x
       ENDIF
       x=x1
    ENDDO
!  write(*,*)'regula_falsi : do loop too shrto!'
5   CONTINUE
!  write(*,100)'result fegula falsi (0=ok,1=zero in 1 knot,2=double zero
!      & 5=error) : ',succ
100 FORMAT(a67,1x,i2)
    x0=x
!  write(*,*)'f(',x0,')=',fx,'it=',k
  END SUBROUTINE regula_falsi


! ********************************************************
!   EVALUATE f(x)
! ********************************************************
  SUBROUTINE evalspline(bi,ci,di,xi,yi,x,fx)
! --------------- interface ------------
    REAL(KIND=dkind),INTENT(IN) :: bi,ci,di,xi,yi,x 
    REAL(KIND=dkind),INTENT(OUT) :: fx
! ----------- end interface ------------
    fx= yi + bi*(x-xi) + ci*(x-xi)**2 + di*(x-xi)**3
  END SUBROUTINE evalspline


! ********************************************************
! COMPUTATION of max apogee and min perigee at each evolution node
! (comparison with Hoot et al. filter)
! ********************************************************
  SUBROUTINE maxmin_apoperi(iunda,pos,peri,apo)
! ---------------- interface ------------------------
    INTEGER,INTENT(IN) :: iunda !direct access iunit
    INTEGER,INTENT(IN) :: pos !evolution position index
    REAL(KIND=dkind),INTENT(OUT) :: peri,apo
! ---------------- end interface -------------------
    REAL(KIND=dkind),DIMENSION(tmax) :: tevol !days from teph0
    REAL(KIND=dkind),DIMENSION(3) :: elequ
    REAL(KIND=dkind),DIMENSION(tmax) :: rmin,rmax
    INTEGER :: h,kaux(1)
! --------------------------------------
    DO h=1,nevol
       READ(iunda,REC=pos+h)tevol(h),elequ(1:3)
       CALL apoperi(elequ,rmax(h),rmin(h))
!       write(pos/400+11,*)tevol(h),rmin(h),rmax(h)
    ENDDO
! Minimizzo i "pericentri" e massimizzo gli "apocentri"
    kaux=MINLOC(rmin(1:nevol))
    peri=rmin(kaux(1))
    kaux=MAXLOC(rmax(1:nevol))
    apo=rmax(kaux(1))
!    write(*,*)'peri=',peri,'apo=',apo
  END SUBROUTINE maxmin_apoperi


! ********************************************************
! COMPUTE apogee and perigee given EQU elements
! ********************************************************
  SUBROUTINE apoperi(elqu,apo,peri)
! --------------- interface ------------
    REAL(KIND=dkind),DIMENSION(3),INTENT(IN) :: elqu
    REAL(KIND=dkind),INTENT(OUT) :: apo,peri
! ----------- end interface ------------
    REAL(KIND=dkind) :: aa,ea
! --------------------------------------
    aa=elqu(1)
    ea=sqrt(elqu(2)**2+elqu(3)**2)
    peri=aa*(1.d0-ea)
    apo=aa*(1.d0+ea)
  END SUBROUTINE apoperi

END MODULE filter1
