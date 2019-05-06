module ftinyopen
  use iso_c_binding
  implicit none

! A C function that returns a string need a pointer to the array of single char 
  type (c_ptr) :: C_String_ptr
! This is the Fortran equivalent to a string of single char
  character (len=1, kind=c_char), dimension(:), pointer :: filchar => null()

!\begin{verbatim}
! Interface to a C routine which opens a window for browsing a file to open
  interface
     function tinyopen(typ) bind(c, name="tinyopen")
       use iso_c_binding
       implicit none
       integer(c_int), value :: typ
       type (C_Ptr) :: tinyopen
     end function tinyopen
  end interface
!\end{verbatim}

contains

!\begin{verbatim}
  subroutine getfilename(typ,filename)
! DUMMY version of Fortran routine to call a C routine to browse for a file name
    implicit none
    integer typ
    character filename*(*)
    write(*,*)'getfilename dummy routine'
    return
  end subroutine getfilename

end module ftinyopen


