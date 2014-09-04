PROGRAM pmain1
!************************************
! main program for testing the free thermodynamic system
!************************************
!
  use cmon1oc
!
  implicit none
!
  character*12 linkdate
  TYPE(gtp_equilibrium_data), pointer :: ceq
! these will be used later for dimensioning things and efaults
  integer intvar(10)
  double precision dblvar(10)
!
! the next line overwritten with current linkdate by linkocdate
  linkdate='01-01-2012'
! intvar(1) must not be negative
  intvar(1)=20
  call init_gtp(intvar,dblvar)
  if(gx%bmperr.ne.0) then
     stop 'Error initiating GTP data structures'
  endif
!  
  call oc_command_monitor(linkdate)
!
! This is the data structure for the default equilibrium
  ceq=>eqlista(1)
!
  write(*,*)'A bientot'
!
  end
