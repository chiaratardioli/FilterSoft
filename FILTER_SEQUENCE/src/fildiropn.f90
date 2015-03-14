!
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 5, 1997                                             
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                      F I L D I R O P N                        * 
!  *                                                               *    
!  *               Unit allocation and file opening                *    
!  *                                                               *
!  *                   D I R E C T    A C C E S S                  *
!  *****************************************************************    
!                                                                       
! INPUT:    NAME      -  File name to be opened                         
!           STATUS    -  Open status                                    
!           LENGTH    -  Length of a record

! OUTPUT:   IUN       -  Allocated unit                                 
!         
! revised on March 2012, C.T.
SUBROUTINE fildiropn(iun,name,status,length)
  USE file_oper
  IMPLICIT NONE 
  INTEGER iun,ll,ls,i,i0,iunold,length 
  CHARACTER*(*) name,status 
  CHARACTER*100 name_rescue
  CHARACTER*40 name1
  LOGICAL opnd,exis, nmd 
                                                                        
  INTEGER lench,io_stat_no,ntry 
  EXTERNAL lench 
  LOGICAL first
  DATA first /.true./
  SAVE first,i0
  IF(first)THEN
     i0=iunf1-1 
     first=.false.
  ENDIF
  IF(iicfil.NE.36) THEN 
     DO  i=iunf1,iunf2 
        allunt(i)=.false.
     ENDDO 
     iicfil=36 
  END IF
  ntry=0
7 INQUIRE(FILE=name,OPENED=opnd,EXIST=exis,NUMBER=iunold) 
  IF(opnd) THEN 
     ll=lench(name) 
     WRITE(*,102) name(1:ll), opnd,exis,iunold,allunt(iunold)
102 FORMAT(' **** filopn: INQUIRE error ****'/                  &
     &       ' **** FILE: ',A,' opnd=',L1,' exis=',L1,' ****'/        &
     &       ' **** UNIT: ',i3,' allunt; ',i3)
     WRITE(*,*)'NUMBER=',iunold
     INQUIRE(UNIT=iunold,OPENED=opnd,NAMED=nmd,NAME=name1,EXIST=exis,ERR=99)
     WRITE(*,*)' UNIT= ',iunold,' NAMED= ',nmd,' NAME= ',name1
!     CLOSE(iunold,ERR=4,IOSTAT=io_stat_no)
!     GOTO 7
     GOTO 5
!!          STOP '**** filopn: abnormal end ****' 
  ELSE
     GOTO 5                       
  END IF
!4 WRITE(*,*)' **** filopn: CLOSE error, IOSTAT: ',io_stat_no
5 CONTINUE                                                          
  DO  iun=i0+1,iunf2 
     IF(allunt(iun)) CYCLE 
     OPEN(iun,FILE=name,ACCESS='DIRECT',RECL=length,STATUS=status, &
          & ERR=3,IOSTAT=io_stat_no)
!  OPEN(11,FILE='initial_propag.db',ACCESS='DIRECT',RECL=length, &
!       & FORM='UNFORMATTED',STATUS='REPLACE') 
     filna_fo(iun)=name 
     allunt(iun)=.true. 
     i0=iun
     ll=lench(name)
!     WRITE(*,*)'filopn: ', iun, name(1:ll)
     RETURN 
  ENDDO
  DO  iun=iunf1,i0-1 
     IF(allunt(iun)) CYCLE 
     OPEN(iun,FILE=name,ACCESS='DIRECT',RECL=length,STATUS=status, &
          & ERR=3,IOSTAT=io_stat_no) 
     filna_fo(iun)=name 
     allunt(iun)=.true.
     i0=iun 
     ll=lench(name)
!     WRITE(*,*)'filopn: ', iun, name(1:ll)
     RETURN 
  ENDDO
  WRITE(*,*)'filopn: all units are already allocated'
  STOP
                                                                        
3 CONTINUE 
  ll=lench(name) 
  ls=lench(status) 
  WRITE(*,101) name(1:ll),status(1:ls),iun,exis,opnd,io_stat_no 
101 FORMAT(' **** filopn: cannot OPEN file "',A,'" (status=',A,       &
     &    '   as ',I4,' EXISTS=',L1,' OPENED=',L1, ' IOSTAT=',i8,') ****') 
  STOP '**** filopn: abnormal end ****' 
! try desperate move
!   name_rescue='rescue/'//name
!   CALL rmsp(name_rescue,ll)
!   OPEN(iun,FILE=name_rescue(1:ll),STATUS=status,ERR=33,IOSTAT=io_stat_no) 
!   WRITE(*,*)'filopn: rescue attempt ', iun, name_rescue(1:ll)
!   filna_fo(iun)=name 
!   allunt(iun)=.true.
!   i0=iun 
!    RETURN
!33  ls=lench(status) 
!    WRITE(*,103) name_rescue(1:ll),status(1:ls),iun,io_stat_no
!103 FORMAT(' **** filopn: cannot OPEN rescue file "',A,'" (status=',A,       &
!      &    '   as ',I4,' IOSTAT=',i8,') ****') 
!    ntry=ntry+1
!    IF(ntry.gt.10)  STOP '**** filopn: abnormal end ****'
!    CALL waste_time(100) 
!    GOTO 5
99 continue
  WRITE(*,*)' error in second INQUIRE, file=',name(1:ll)
  STOP '**** filopn: abnormal end ****' 
END SUBROUTINE fildiropn
