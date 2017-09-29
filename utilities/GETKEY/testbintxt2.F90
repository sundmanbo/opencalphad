program onechar
  !
  use bintxt
  !
  implicit none
  integer ll
  character cline*80
!
  type(history) :: hist
  !
  ll=0
  hist%hpos=0
!
5 continue
  write(*,10,advance='no')
10 format('type something:')
  call bintxt2(5,cline,hist)
  ll=ll+1
  write(*,11)ll,hist%hpos,trim(cline)
11 format(/'line: ',2i3,': ',a)
  if(cline(1:1).eq.'q' .or. ll.gt.10) stop
  goto 5
end program onechar

