module cstr
contains
function c_to_f_string(s) result(str)
  use iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  integer i, nchars
  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str)
  str = transfer(s(1:nchars), str)
end function c_to_f_string
end module cstr


MODULE c_f_string_ptr
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: C_F_STRING
CONTAINS
  FUNCTION C_F_STRING(c_str) RESULT(f_str)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_CHAR
    TYPE(C_PTR), INTENT(IN) :: c_str
    CHARACTER(:,KIND=C_CHAR), POINTER :: f_str
    CHARACTER(KIND=C_CHAR), POINTER :: arr(:)
    INTERFACE
      ! Steal std C library function rather than writing our own.
      FUNCTION strlen(s) BIND(C, NAME='strlen')
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_SIZE_T
        IMPLICIT NONE
        !----
        TYPE(C_PTR), INTENT(IN), VALUE :: s
        INTEGER(C_SIZE_T) :: strlen
      END FUNCTION strlen
    END INTERFACE
    !****
    print *,"Fstring"
    print *,"Length strlen", strlen(c_str)
    CALL C_F_POINTER(c_str, arr, [strlen(c_str)])
   
    CALL get_scalar_pointer(SIZE(arr), arr, f_str)
    
  END FUNCTION C_F_STRING
  SUBROUTINE get_scalar_pointer(scalar_len, scalar, ptr)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
    INTEGER, INTENT(IN) :: scalar_len
    CHARACTER(KIND=C_CHAR,LEN=scalar_len), INTENT(IN), TARGET :: scalar(1)
    CHARACTER(:,KIND=C_CHAR), INTENT(OUT), POINTER :: ptr
    !***
    ptr => scalar(1)
  END SUBROUTINE get_scalar_pointer
END MODULE c_f_string_ptr

module test
        use iso_c_binding
        use cstr
        use liboctq
        use general_thermodynamic_package
        implicit none

        TYPE, bind(c) :: c_gtp_equilibrium_data 
! this contains all data specific to an equilibrium like conditions,
! status, constitution and calculated values of all phases etc
! Several equilibria may be calculated simultaneously in parallell threads
! so each equilibrium must be independent 
! NOTE: the error code must be local to each equilibria!!!!
! During step and map these records with results are saved
! values of T and P, conditions etc.
! Values here are normally set by external conditions or calculated from model
! local list of components, phase_varres with amounts and constitution
! lists of element, species, phases and thermodynamic parameters are global
! tpval(1) is T, tpval(2) is P, rgas is R, rtn is R*T
! status: not used yet?
! multiuse: used for various things like direction in start equilibria
! eqno: sequential number assigned when created
! next: index of next equilibrium in a sequence during step/map calculation.
! eqname: name of equilibrium
! tpval: value of T and P
! rtn: value of R*T
     integer(c_int) :: status,multiuse,eqno,next
     character(c_char) :: eqname*24
     real(c_double) :: tpval(2),rtn
! svfunres: the values of state variable functions valid for this equilibrium
     type(c_ptr) :: svfunres
! the experiments are used in assessments and stored like conditions 
! lastcondition: link to condition list
! lastexperiment: link to experiment list
     TYPE(c_ptr) :: lastcondition,lastexperiment
! components and conversion matrix from components to elements
! complist: array with components
! compstoi: stoichiometric matrix of compoents relative to elements
! invcompstoi: inverted stoichiometric matrix
     TYPE(c_ptr) :: complist
     real(c_double) :: compstoi
     real(c_double) :: invcompstoi
! one record for each phase+composition set that can be calculated
! phase_varres: here all calculated data for the phase is stored
     TYPE(c_ptr) :: phase_varres
! index to the tpfun_parres array is the same as in the global array tpres 
! eq_tpres: here local calculated values of TP functions are stored
     TYPE(c_ptr) :: eq_tpres
! current values of chemical potentials stored in component record but
! duplicated here for easy acces by application software
     real(c_double) :: cmuval
! xconc: convergence criteria for constituent fractions and other things
     real(c_double) :: xconv
! delta-G value for merging gridpoints in grid minimizer
! smaller value creates problem for test step3.BMM, MC and austenite merged
     real(c_double) :: gmindif=-5.0D-2
! maxiter: maximum number of iterations allowed
     integer(c_int) :: maxiter
! this is to save a copy of the last calculated system matrix, needed
! to calculate dot derivatives, initiate to zero
     integer(c_int) :: sysmatdim=0,nfixmu=0,nfixph=0
     integer(c_int) :: fixmu
     integer(c_int) :: fixph
     real(c_double) :: savesysmat
  END TYPE c_gtp_equilibrium_data
contains
  subroutine examine_gtp_equilibrium_data(c_ceq) bind(c, name='examine_gtp_equilibrium_data')
     type(c_ptr), intent(in), value :: c_ceq
     type(gtp_equilibrium_data), pointer :: ceq
     integer :: i,j
     call c_f_pointer(c_ceq, ceq)
     write(*,10) ceq%status, ceq%multiuse, ceq%eqno
10   format(/'gtp_equilibrium_data: status, multiuse, eqno, next'/, 3i4)
     write(*,20) ceq%eqname
20   format(/'Name of equilibrium'/,a)
     write(*,30) ceq%tpval, ceq%rtn
30   format(/'Value of T and P'/, 2f8.3, /'R*T'/, f8.4)
     do i = 1, size(ceq%compstoi,1)
        write(*,*) (ceq%compstoi(i,j), j=1,size(ceq%compstoi,2))
     end do
     write(*,*) ceq%cmuval
     write(*,*) ceq%xconv
     write(*,*) ceq%gmindif
     write(*,*) ceq%maxiter
     write(*,*) ceq%sysmatdim, ceq%nfixmu, ceq%nfixph
     write(*,*) ceq%fixmu, ceq%fixph, ceq%savesysmat
  end subroutine examine_gtp_equilibrium_data
!\begin{verbatim}
  subroutine c_tqini(n, c_ceq) bind(c, name='c_tqini')
    integer(c_int), intent(in) :: n
    type(c_ptr), intent(out) :: c_ceq
!\end{verbatim}  
    type(gtp_equilibrium_data), pointer :: ceq
    integer :: i1,i2
    print *, "c_ceq pointer", i1
    print *, "ceq pointer", i2
    
    call tqini(n, ceq)
    i2 = transfer(c_ceq, i2)
    print *, "c_ceq pointer before transfer", i2
    c_ceq = c_loc(ceq)
    i2 = transfer(c_ceq, i2)
    print *, "c_ceq pointer after transfer", i2
  end subroutine c_tqini

!\begin{verbatim}
  subroutine c_tqrfil(filename,c_ceq) bind(c, name='c_tqrfil')
    character(kind=c_char,len=1), intent(in) :: filename(*)
    character(len=:), allocatable :: fstring
    type(gtp_equilibrium_data), pointer :: ceq
    type(c_ptr), intent(inout) :: c_ceq
      integer :: i2
      i2 = transfer(c_ceq,i2)
!\end{verbatim}
! convert type(c_ptr) to fptr
    print *,"c_ceq pointer", i2
    call c_f_pointer(c_ceq, ceq)
    fstring = c_to_f_string(filename)
    print *,"Filename:", fstring
    print *,"ceq%status", ceq%status
    call tqrfil(fstring, ceq)
    c_ceq = c_loc(ceq)
  end subroutine c_tqrfil
  
!\begin{verbatim}
  subroutine c_tqrpfil(filename,nel,c_selel,c_ceq) bind(c)
    character(kind=c_char), intent(in) :: filename
    integer(c_int), intent(in) :: nel
    type(c_ptr), intent(in), dimension(nel), target :: c_selel
    type(c_ptr), intent(inout) :: c_ceq  
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=:), allocatable :: fstring
    character, pointer :: selel(:)
    integer :: i
    fstring = c_to_f_string(filename)
    call c_f_pointer(c_ceq, ceq)
    ! convert the c type selel strings to f-selel strings
    ! note: additional character is for C terminated '\0'
    do i = 1, nel
        call c_f_pointer(c_selel(i), selel, [3])
    end do
    call tqrpfil(fstring, nel,selel,ceq)
    c_ceq = c_loc(ceq)
  end subroutine c_tqrpfil
!\begin{verbatim}
  subroutine c_tqgcom(n,components,c_ceq) bind(c, name='c_tqgcom')
! get system components
    integer(c_int), intent(inout) :: n
    !character(kind=c_char, len=24), dimension(24), intent(out) :: c_components
    type(c_ptr), intent(inout) :: c_ceq  
!\end{verbatim}
    integer, target :: nc
    character(kind=c_char, len=1), dimension(maxel) :: components
    type(gtp_equilibrium_data), pointer :: ceq  
    integer :: i2
    i2 = transfer(c_ceq,i2)
    print *,"c_ceq pointer from tqgcom", i2
    call c_f_pointer(c_ceq, ceq)
    print *,"status ", ceq%status
    print *, "Number of components", nc
    !do i = 0,10
    !  c_components = c_to_f_string(components(i))
    !end do
    call tqgcom(nc, components, ceq)
    ! convert the F components strings to C 
    !do i = 1, nc
       
    !end do
    write(*,*) 'Return from tqgcom',nc
    c_ceq = c_loc(ceq)
    write(*,*) 'c_loc sucessed'
    n = nc
    write(*,*) 'c_ceq location'
    write(*,*) 'End of tqgcom'
    end subroutine c_tqgcom  
!\begin{verbatim}
  subroutine c_tqgnp(n, c_ceq) bind(c, name='c_tqgnp')
    integer(c_int), intent(inout) :: n
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    call tqgnp(n, ceq)
    c_ceq = c_loc(ceq)
  end subroutine c_tqgnp
!\begin{verbatim}
  subroutine c_tqgpn(n,phasename, c_ceq) bind(c, name='c_tqgpn')
! get name of phase n,
! NOTE: n is phase number, not extended phase index
    integer(c_int), intent(in) :: n
    character(len=1), intent(inout) :: phasename(24)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fstring
    integer :: i
    call c_f_pointer(c_ceq, ceq)
    !fstring = c_to_f_string(phasename)
    call tqgpn(n, fstring, ceq)
    print *, "Name of phase:", fstring
    do i=1,len(fstring)
        phasename(i)(1:1) = fstring(i:i)
    end do
    c_ceq = c_loc(ceq)
  end subroutine c_tqgpn 

end module test
