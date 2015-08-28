!
! gtp3Z included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>      20. Subrotines used by applications
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine makeoptvname(name,indx)
    implicit none
    character name*(*)
    integer indx
!\end{verbatim}
    if(indx.lt.99) then
       if(indx.le.9) then 
          name(1:2)='A0'
          name(3:3)=char(indx+ichar('0'))
       else
          name(1:1)='A'
          name(2:2)=char(indx/10+ichar('0'))
          name(3:3)=char(mod(indx,10)+ichar('0'))
       endif
    else
       name='A99'
    endif
1000 continue
    return
  end subroutine makeoptvname

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!-\begin{verbatim}
  subroutine dummy(apphead)
!
    implicit none
    TYPE(gtp_applicationhead), pointer :: apphead
!-\end{verbatim}
!   
1000 continue
    return
  end subroutine dummy

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

