!=======================================================================--------
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!=======================================================================--------
! These routines are available for general use. I ask that you send me
! interesting alterations that are available for public use; and that you
! include a note indicating the original author --  John S. Urban
! Last updated May 5th, 2009
!=======================================================================--------
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!=======================================================================--------
! make Fortran/C interface for C routine getkey(3C)
module M_getkey
   use iso_c_binding
   implicit none
   public
      interface
         function getkey() bind(c, name='getkeyC')
            use iso_c_binding
            implicit none
            character(kind=c_char) :: getkey
         end function getkey
      end interface
end module M_getkey
!=======================================================================--------
!()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
!=======================================================================--------
!-------------------------------------------------------------------------------

program test_getkey
   use M_getkey
   character :: A
   integer :: icount
   icount=0
   write(*,*)'begin striking keys to demonstrate interactive raw I/O mode'
   write(*,*)'q to quit; up to 40 tries allowed'
   do
      A=getkey()
      icount=icount+1
      write(*,*)icount,' f03:key=',A,'->',ichar(A)
      !flush(6)
      if(A.eq.'q')stop
      if(icount.gt.40)stop
   enddo
end program test_getkey

