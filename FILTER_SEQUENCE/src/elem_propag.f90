! Chiara Tardioli (tardioli@mail.dm.unipi.it)
!
! ****************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                         !!
!!  Secular perturbations of the first order of the orbit  !! 
!!   of a satellite/debris about the Earth togheter with   !!
!!          their derivatives w.r.t. input elements        !!
!!                                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Chiara Tardioli, July 2011
! revised on March 2012, C.T.
!
! Reference: A.E. Roy, 'Orbital Motion'
! ****************************************************
SUBROUTINE elem_propag(deltat,elkep,unc_in,elkout,unc_out)
  USE fund_const
  USE orbit_elements
  USE planet_masses, ONLY: gmearth
  USE spher_harm, ONLY: harmco
  IMPLICIT NONE
! ----------------------- interface ----------------------
  REAL(kind=dkind),INTENT(IN) :: deltat !propagation time step
  TYPE(orbit_elem),INTENT(IN) :: elkep  !input KEP elems
 ! REAL(kind=dkind),DIMENSION(6),INTENT(IN) :: elkep  !input KEP elems
  REAL(kind=dkind),DIMENSION(5,5),INTENT(IN) :: unc_in !input unc matrix
! --------------------------------------------------------
! output KEP elems, obtained from Roy, 'Orbital Motion', pag.285
  TYPE(orbit_elem),INTENT(OUT) :: elkout
!  REAL(kind=dkind),DIMENSION(6),INTENT(OUT) :: elkout
  REAL(kind=dkind),DIMENSION(5,5),INTENT(OUT) :: unc_out !output inc matrix
! -------------------- end  interface --------------------
! constant term in the second harmonic of Earth's potential
  REAL(kind=dkind) :: j2
! elkep%coord: 1-a, 2-ecc, 3-inc, 4-Omega.nod, 5-omega, 6-Mean anomaly
  REAL(kind=dkind),DIMENSION(6) :: elem !auxiliar KEP elems
  REAL(KIND=dkind) :: a0AU
  REAL(kind=dkind) :: pe !parameter of the conic
  REAL(kind=dkind) :: n0 !input mean motion
  REAL(kind=dkind) :: nbar !propag mean motion
!  REAL(kind=dkind) :: deltat
  REAL(kind=dkind) :: ci,si,si2,beta,beta2 !function of orbit elems
  REAL(kind=dkind) :: fc,nt !terms in averaged solutions
  REAL(kind=dkind) :: da,de,dni !variables for derivatives
! Derivatives of propagated elements w.r.t. input elements
  REAL(kind=dkind),DIMENSION(5,5) :: derE
! check with incremental ratio
  LOGICAL check_der
  REAL(kind=dkind),DIMENSION(6) :: elem_inc
  REAL(kind=dkind) :: incr(5)
  INTEGER :: i ! loop index
! ====================================================================

!  check_der=.true.
  check_der=.false.

  elkout=elkep

  IF(deltat.eq.0.d0)THEN
!     elkout=elkep
     unc_out=unc_in
     RETURN
     write(*,*)'elem_propag : you shoud not be here'
     STOP
  ENDIF

  elkout%t=elkep%t+deltat

  j2 = -harmco(5)

! ***** parameter *******
!  write(*,*)'J2 =',j2,' equatorial radius = ',reau
!  write(*,*)'Earth mass =',gmearth
!  write(*,*)'aukm=',aukm
! ***********************

  elem=elkep%coord(1:6)
  a0AU = elem(1)/aukm  !elem(1) in AU  
  beta2 = 1-elem(2)**2
  pe = elem(1)*beta2
  ci = cos(elem(3))
  si = sin(elem(3))
  si2 = si**2
  beta = sqrt(beta2)
  n0 = sqrt(gmearth)/a0AU**1.5d0
  fc = 1.5d0*j2*(reau*aukm)**2/pe**2
  nbar = n0*(1.d0+fc*(1.d0-1.5d0*si2)*beta)
  nt = nbar*deltat

! a,ecc,inc are costants
! Omega.nod, omega and Mean anomaly change linearly with time 
  elem(4) = elem(4) - fc*ci*nt !averaged Omega.nod
  elem(5) = elem(5) + fc*(2.d0-2.5d0*si2)*nt !averaged omega
  elem(6) = elem(6) + nt !averaged mean anomaly

!  write(*,*)'elem_propag:', -fc*ci*nbar,fc*(2.d0-2.5d0*si2)*nbar
!  write(*,*)'elem_propag: elem out',elkout%coord(1:6) 

! Propagated elements at step time deltat
  elkout%coord(1:6) = elem(1:6)
!  write(*,*)'elkep',elkep%coord-elkout%coord
! ---------------------------------------------------------------

! Elements derivatives w.r.t. input elements
  derE = 0.d0 !5x5 matrix

! (derivatives of nbar/p^2 wrt a,e)*p^2  
  da = (-11.d0/2.d0*nbar+2.d0*n0)/elem(1) 
!  de = elem(2)/beta2*(7.d0*nbar-3.d0*n0)
  de = elem(1)*elem(2)/pe*(7.d0*nbar-3.d0*n0)
! derivative of nbar wrt I
  dni = -3.d0*n0*fc*si*ci*beta 

  derE(1,1) = 1.d0 !da/da0
  derE(2,2) = 1.d0 !de/de0
  derE(3,3) = 1.d0 !di/di0

  derE(4,1) =-fc*da*ci*deltat            !dOmnod/da0
  derE(4,2) =-fc*de*ci*deltat            !dOmnod/de0
  derE(4,3) =-fc*(dni*ci-nbar*si)*deltat !dOmnod/di0
  derE(4,4) = 1.d0                       !dOmnod/dOmnod0

  derE(5,1) = fc*da*(2.d0-2.5d0*si2)*deltat                    !dom/da0
  derE(5,2) = fc*de*(2.d0-2.5d0*si2)*deltat                    !dom/de0
  derE(5,3) = fc*(dni*(2.d0-2.5d0*si2)-5.d0*nbar*si*ci)*deltat !dom/di0
  derE(5,5) = 1.d0                                             !dom/dom0
     
  IF(check_der)THEN
     i=2
     CALL check_der_inc(i,elkep%coord(1:6),deltat,elem_inc,incr)
     WRITE(*,*)'derivatives of Omnod'
     WRITE(*,*)derE(4,i)
     write(*,*)'derivatives of Omnod with incremental ratio'
     write(*,*)i,(elem_inc(4)-elem(4))/incr(i)
     write(*,*)'--------------------------------------------'
     WRITE(*,*)'derivatives of omega'
     WRITE(*,*)derE(5,i)
     write(*,*)'derivatives of omega with incremental ratio'
     write(*,*)i,(elem_inc(5)-elem(5))/incr(i)
  ENDIF

! compute the normal matrix Gamma of the propagated elements at epoch t1
  unc_out(1:5,1:5) = MATMUL(MATMUL(derE,unc_in(1:5,1:5)),TRANSPOSE(derE))

END SUBROUTINE elem_propag


SUBROUTINE check_der_inc(i,elkepin,dt,elem,incr)
  USE fund_const
  USE orbit_elements
  USE planet_masses, ONLY: gmearth
  USE spher_harm, ONLY: harmco
  IMPLICIT NONE
! -------- interface ---------------------
  INTEGER,INTENT(IN) :: i
  REAL(Kind=dkind),INTENT(IN) :: elkepin(6),dt
  REAL(kind=dkind),INTENT(OUT) :: elem(6),incr(5)
! -------- end interface ----------------
  REAL(KIND=dkind) :: eps
  REAL(KIND=dkind) :: a0AU,j2
  REAL(kind=dkind) :: pe
  REAL(kind=dkind) :: n0
  REAL(kind=dkind) :: nbar
  REAL(kind=dkind) :: deltat
  REAL(kind=dkind) :: ci,si,si2,beta,beta2
  REAL(kind=dkind) :: fc,nt
  REAL(kind=dkind) :: da,de,dni
! ========================================

  j2 = -harmco(5)
  elem = elkepin

  incr=0.d0
  eps = 1.d-4
  
  incr(i) = elkepin(i)*eps
  elem(i) = elkepin(i)+ incr(i)
!  write(*,*) 'incr of ',i,':',incr(i)
  a0AU = elem(1)/aukm  !elem(1) in AU  
  beta2 = 1-elem(2)**2
  pe = elem(1)*beta2
  ci = cos(elem(3))
  si = sin(elem(3))
  si2 = si**2
  beta = sqrt(beta2)
  n0 = sqrt(gmearth)/a0AU**1.5d0
  fc = 1.5d0*j2*(reau*aukm)**2/pe**2
  nbar = n0*(1.d0+fc*(1.d0-1.5d0*si2)*beta)
  nt = nbar*dt
  elem(4) = elem(4) - fc*ci*nt
  elem(5) = elem(5) + fc*(2.d0-2.5d0*si2)*nt
  elem(6) = elem(6) + nt 
!  write(*,*)'elem',elem

END SUBROUTINE check_der_inc
