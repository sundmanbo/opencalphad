PROGRAM pmain1
!************************************
! main program for testing the free thermodynamic system
!************************************
!
  use oc_cmon1
!
  implicit none
!
  character*12 linkdate
  TYPE(gtp_equilibrium_data), pointer :: ceq
!
! the next line overwritten with current linkdate by linkocdate
  linkdate='01-01-2012'
  call init_gtp
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
