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

