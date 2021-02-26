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
! Fortran routine to call a C routine to browse for a file name
! typ if default extension:
! 1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=PLT ,6=PDB, 7=DAT, 8=LOG
    character (len=256) :: filename
    integer typ
!\end{verbatim}
    integer jj
! the current directory can be found by 
! character directory*128
! call getcwd(directory)
! specify a name of a file type:
!    write(*,*)'ftinyopen ',typ
! 1=TDB, 2=OCU, 3=OCM, 4=OCD, 5=plt, 6=PDB, 7=DAT, 8=LOG
    C_String_ptr = tinyopen(typ)
! convert C pointer to Fortran pointer
    call c_f_pointer(C_String_ptr,filchar,[256])
    filename=' '
    if(associated(filchar)) then
! convert the array of single characters to a Fortran character
       jj=1
       do while(filchar(jj).ne.c_null_char)
          filename(jj:jj)=filchar(jj)
          jj=jj+1
       enddo
    endif
!    write(*,*)'ftinyopen getfilename: ',trim(filename),typ
1000 continue
    return
  end subroutine getfilename

end module ftinyopen


