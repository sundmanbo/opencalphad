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
! these will be used later for dimensioning things and efaul
  integer i,narg,intvar(10)
  double precision dblvar(10)
  character arginline(4)*256,arg*64,date*16
!
! save the data of linking the program
!  call date_and_time(date)
!  write(*,*)'Stored linking date: ',date
! This line replaced by linkocdate to the date when compiling this program.
  linkdate=date(1:4)//'-'//date(5:6)//'-'//date(7:8)
! for example: linkdate='2019-11-27'
! this is the overall version identifier
  version='  6.030 '
! intvar and dblvar will eventually be used for allocations
  intvar(1)=30
  call init_gtp(intvar,dblvar)
  if(gx%bmperr.ne.0) then
     stop 'Error initiating GTP data structures'
  endif
! extract arguments from the line of invocation
! at present just a macro file name
   narg=iargc()
  if(narg.gt.4) then
     write(*,*)'OC accepts max 4 inline arguments'
     narg=4
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
! The data structure for the default equilibrium is in eqlista
  ceq=>eqlista(1)
! additional code can be added below for some particular app
!
  write(*,*)'A bientot'
  call deallocate_gtp(intvar,dblvar)
!
end PROGRAM pmain1
