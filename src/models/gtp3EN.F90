!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!
! These are dummy subroutines replacing those needed for reading
! encrypted databases.
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine readencrypt
!\begin{verbatim}
 subroutine readencrypt(line,nr)
! This is a dummy routine because the real one is not part of open source
! line contans "ENCRYPTED filename (and maybe more)
! nr is missing functions reading here
   character line*(*)
   integer nr
!\end{verbatim}
   nr=0
   write(*,*)'3Z This OC version cannot read encrypted databases'
   gx%bmperr=4399
   return
 end subroutine readencrypt

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine notallowlisting
!\begin{verbatim}
 logical function notallowlisting(privil)
! check if user is allowed to list data
   double precision privil
!\end{verbatim}
   logical ok
! false means listing allowed
   ok=.FALSE.
   notallowlisting=ok
   return
 end function notallowlisting

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

