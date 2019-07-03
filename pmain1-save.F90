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
  integer i,narg,intvar(10)
  double precision dblvar(10)
  character arginline(12)*64,arg*64
!
! the next line overwritten with current linkdate by linkocdate
  linkdate='2015-01-01'
! this is the overall version identifier
  version='  5.043 '
! intvar and dblvar will eventually be used for allocations and defaults
  intvar(1)=30
  call init_gtp(intvar,dblvar)
  if(gx%bmperr.ne.0) then
     stop 'Error initiating GTP data structures'
  endif
! extract arguments from the line of invocation
  narg=iargc()
  if(narg.gt.12) then
     write(*,*)'OC accepts max 12 inline arguments'
     narg=12
!  else
!     write(*,*)'Inline arguments: ',narg
  endif
  do i=1,narg
     call getarg(i,arginline(i))
!     write(*,*)trim(arginline(i))
  enddo
!  
  call oc_command_monitor(version,linkdate,narg,arginline)
!
! we come back here with the "back" command in the user i/f
! The data structure for the default equilibrium is in eqlista(1)
  ceq=>eqlista(1)
! additional code can be added below for some particular application.
!
  write(*,*)'A bientot'
  call deallocate_gtp(intvar,dblvar)
!
end PROGRAM pmain1
