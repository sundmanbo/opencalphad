PROGRAM pmain1
!************************************
! main program for the free Open Calphad software
!************************************
!
  use cmon1oc
!
! For parallel processing
!$  use omp_lib
!
  implicit none
!
  character linkdate*12,version*8
  TYPE(gtp_equilibrium_data), pointer :: ceq
! these will be used later for dimensioning things and efaults
  integer intvar(10)
  double precision dblvar(10)
!
! the next line overwritten with current linkdate by linkocdate
  linkdate='01-01-2015'
! this is the overall version identifier
  version='  4.015 '
! intvar and dblvar will eventually be used for allocations and defaults
  intvar(1)=30
  call init_gtp(intvar,dblvar)
  if(gx%bmperr.ne.0) then
     stop 'Error initiating GTP data structures'
  endif
!  
  call oc_command_monitor(version,linkdate)
!
! This is the data structure for the default equilibrium
! additional code can be added below for some particular application.
  ceq=>eqlista(1)
!
  write(*,*)'A bientot'
  call deallocate_gtp(intvar,dblvar)
!
  end
