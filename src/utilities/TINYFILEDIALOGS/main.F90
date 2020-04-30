program test3
  use iso_c_binding
  implicit none

! A C function that returns a string need a pointer to the array of single char 
  type (c_ptr) :: C_String_ptr
! This is the Fortran equivalent to a string of single char
  character (len=1, kind=c_char), dimension(:), pointer :: filchar => null()

! Interface to a C routine which opens a window for browsing a file to open
  interface
     function tinyopen(typ) bind(c, name="tinyopen")
       use iso_c_binding
       implicit none
       integer(c_int), value :: typ
       type (C_Ptr) :: tinyopen
     end function tinyopen
  end interface

! variables used in Fortran
  character (len=256) :: filename
  integer typ,jj

! ------------------------------------
! specify a name of a TDB file (typ=1)
  typ=1
  C_String_ptr = tinyopen(typ)
! convert C pointer to Fortran pointer
  call c_f_pointer(C_String_ptr,filchar,[256])
  filename=' '
  if(.not.associated(filchar)) then
! if no characters give error message
     write(*,*)'No file name'
  else
! convert the array of single characters to a Fortran character
     jj=1
     do while(filchar(jj).ne.c_null_char)
        filename(jj:jj)=filchar(jj)
        jj=jj+1
     enddo
  endif
  write(*,*)'File name is: ',trim(filename)
! ------------------------------------
! specify a name of an unformatted file (typ=2)
  typ=2
  C_String_ptr = tinyopen(typ);
! convert C pointer to Fortran pointer
  call c_f_pointer(C_String_ptr,filchar,[256])
  filename=' '
  if(.not.associated(filchar)) then
! if no characters give error message
     write(*,*)'No file name'
  else
! convert the array of single characters to a Fortran character
     jj=1
     do while(filchar(jj).ne.c_null_char)
        filename(jj:jj)=filchar(jj)
        jj=jj+1
     enddo
  endif
  write(*,*)'File name is: ',trim(filename)
  write(*,90)
90 format(/'All well that ends well'/)
end program test3

