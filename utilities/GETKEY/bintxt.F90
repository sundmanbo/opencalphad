!
! The c program getkey-OK.c must be compiled with
! cc -c -DBSD getkey-OK.c
!
! The interface must be compiled with
! gfortran -c M_getkey-OK.F90
!
! This module must be compiled with
! gfortran -c bintxt.F90
!
! The test program can be compiled with
! gfortran -o tin testbintxt2.F90 bintxt.o getkey-OK.o
!
module bintxt
  
  use M_getkey

  implicit none

  type history
     character*128 hline(20)
     integer hpos
  end type history
  
contains
  
!\begin verbatim
  subroutine bintxt2(lin,cline,hist)
! subroutine to read a line with history and editing on LINUX a la emacs
!
    character cline*(*)
    integer lin
    type(history) :: hist
!\end verbatim
!--------------------  
! CONTROL CHARACTERS FROM KEYBOARD
! DEL delete curret character
    integer, parameter :: ctrla=1        ! CTRLA move cursor to first position
    integer, parameter :: backspace2=2   ! CTRLB move cursor one step left
    integer, parameter :: ctrlc=3        ! CTRLC terminate program
    integer, parameter :: ctrld=4        ! CTRLD delete char at cursor
    integer, parameter :: ctrle=5        ! CTRLE move cursor to last position
    integer, parameter :: forward=6      ! CTRLF move cursor one step right
    integer, parameter :: HELP=8         ! CTRLH give coordinates and update
    integer, parameter :: TAB=9          ! CTRLI end of input
    integer, parameter :: ctrlk=11       ! CTRLK delete to end of line
    integer, parameter :: return=13      ! CTRLM end of input
    integer, parameter :: DEL=127        ! DEL delete char left of cursor
    integer, parameter :: mode=17        ! CTRLQ toggle insert/replace
! on MAC same as UP DOWN FORWARD suck
    integer, parameter :: backspace=27   ! CTRL[ previous in history
!--------------------  
! UP previous history line (if any)
! DOWN and LF next history line (if any)
    integer, parameter :: CTRLP=16       ! CTRLP previous in history
!    integer, parameter :: UP=27         ! uparrow previous in history
    integer, parameter :: LF=10          ! CTRLJ next in history
!--------------------
! backspace on a MAC screen
    integer, parameter :: tbackspace=8
!-----------
!
    integer ip,lastp,kiud,jj,kou,size,hlast
    character line*128,ch1*1
    logical insert
!
    kiud=5; kou=6
    size=1
!
!      write(*,*)'Reading input using getkey',lin,kiud
    if(lin.ne.kiud) then
! reading macro from file
       read(lin,10)cline
10     format(a)     
       goto 1000
    endif
!
! input trom terminal with editing, always set intert TRUE at start
!
    insert=.TRUE.
    ip=0
    lastp=0
    line=' '
    hlast=hist%hpos+1
! read one character at a time without echo
100 continue
! read one character at a time without echo and allow editing and history
    ch1=getkey()
!    write(*,*)'got from getkey: ',ichar(ch1)
! handle control character
    if(ichar(ch1).ge.32 .and. ichar(ch1).lt.127) then
! printable character, write on screen and store inline     
       if(ip.eq.lastp .or. .not.insert) then
          write(kou,10,advance='no')ch1
          ip=ip+1
          lastp=lastp+1
          line(ip:ip)=ch1
       else
          line(ip+2:)=line(ip+1:)
          line(ip+1:ip+1)=ch1
          lastp=lastp+1
          write(kou,10,advance='no')tbackspace
          write(kou,10,advance='no')line(ip:lastp)
          do jj=ip,lastp-2
             write(kou,10,advance='no')tbackspace
          enddo
          ip=ip+1
       endif
       goto 100
    endif
!=======================  
!    write(*,*)'control character: ',ichar(ch1)
    select case(ichar(ch1))
    case default
! ignore
!       write(*,*)'Ignoring ',ichar(ch1)
       goto 100
!............. OK
    case(ctrla)
! move cursor to first character
       do jj=1,ip
          write(kou,10,advance='no')tbackspace
       enddo
       ip=0
!............. OK
    case(backspace,backspace2) ! ctrlb leftarrow (also up/down/right arrow)
! move cursor one step back
       if(ip.ge.1) then
          write(kou,10,advance='no')tbackspace
          ip=ip-1
       endif
!............. OK
    case(ctrle)
! move cursor after last character
       do jj=ip+1,lastp
          write(kou,10,advance='no')line(jj:jj)
       enddo
       ip=lastp
!............. OK
    case(ctrld)
! delete character at cursor (ctrld)
       if(ip.eq.lastp) goto 100
       jj=ip+1
! remove the character at position jj and write the whole line from jj to end
       line(jj:)=line(jj+1:)
       write(kou,10,advance='no')line(jj:lastp)
       lastp=lastp-1
! NOTE ip and lastp can be zero here
       do jj=lastp,ip,-1
          write(kou,10,advance='no')tbackspace
       enddo
!............. OK
    case(del)
! delete character to the left of cursor (del)
       if(ip.ge.1) then
          write(kou,10,advance='no')tbackspace
! remove the character at position jj and write the whole line from jj to end
          jj=ip
          line(jj:)=line(jj+1:)
          write(kou,10,advance='no')line(jj:lastp)
          ip=ip-1
       endif
       lastp=lastp-1
! NOTE ip and lastp can be zero here
! othewise we should backspace lastp-ip positions
       do jj=lastp,ip,-1
          write(kou,10,advance='no')tbackspace
       enddo
!............. OK
    case(ctrlk)
! delete all characters from cursor to end of line
       if(ip.lt.lastp) then
          line(ip+1:)=' '
          if(ip.gt.0) write(kou,10,advance='no')tbackspace
          write(kou,10,advance='no')line(ip:lastp)
          do jj=ip,lastp-1
             write(kou,10,advance='no')tbackspace
          enddo
          lastp=ip
       endif
!.............
    case(help)
       write(kou,77, advance='no')ip,lastp,line(1:lastp+1)
77     format(/'Current local values are: ',2i4/a)
!       write(kou,10,advance='no')'xyz'
       do jj=lastp,ip,-1
          write(kou,10,advance='no')tbackspace
       enddo
!.............
    case(mode) ! crtlQ
! toggle insert/overwrite mode
       if(insert) then
          insert=.FALSE.
       else
          insert=.TRUE.
       endif
!.............
    case(return,tab)
! save line (if not empty) finish editing and return current line
       cline=line
       if(len_trim(line).eq.0 .or. &
            line(1:ip+1).eq.hist%hline(hist%hpos)(1:ip+1)) then
          continue
!            write(*,*)'Not saving same or empty line'
       else
          if(hist%hpos.eq.20) then
             do jj=2,20
                hist%hline(jj-1)=hist%hline(jj)
             enddo
          else
             hist%hpos=hist%hpos+1
          endif
          hist%hline(hist%hpos)=line
       endif
       goto 1000
!............. OK
    case(forward)
! move cursor one step right
       if(ip.le.lastp) then
          write(kou,10,advance='no')line(ip:ip)
          ip=ip+1
       endif
!.............
    case(ctrlp)
! copy previous history line to current
! first remove anything on the line (not the question ...)
       if(hlast.gt.1) then
          do jj=1,ip
             write(kou,10,advance='no')tbackspace
          enddo
          line=' '
          write(kou,10,advance='no')line(1:ip)
          do jj=1,ip
             write(kou,10,advance='no')tbackspace
          enddo
          hlast=hlast-1
          line=hist%hline(hlast)
          ip=len_trim(line)
          lastp=ip
          write(kou,10,advance='no')line(1:ip)
       endif
!.............CTRLJ
    case(lf)
! copy next history line to current
       if(hlast.lt.hist%hpos) then
          do jj=1,ip
             write(kou,10,advance='no')tbackspace
          enddo
          line=' '
          write(kou,10,advance='no')line(1:ip)
          do jj=1,ip
             write(kou,10,advance='no')tbackspace
          enddo
          hlast=hlast+1
          line=hist%hline(hlast)
          ip=len_trim(line)
          lastp=ip
          write(kou,10,advance='no')line(1:ip)
       endif
    end select
!-----------------
    goto 100
!=================  
1000 continue
!      call system('stty echo ')
    return
  end subroutine bintxt2
end module bintxt
  
