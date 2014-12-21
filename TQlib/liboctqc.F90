! The name of this library
module liboctq_interface
  use iso_c_binding
  implicit none
  use general_thermodynamic_package
  
  TYPE, bind(c) :: gtp_equilibrium_data 
! this contains all data specific to an equilibrium like conditions,
! status, constitution and calculated values of all phases etc
! Several equilibria may be calculated simultaneously in parallell threads
! so each equilibrium must be independent 
! NOTE: the error code must be local to each equilibria!!!!
! During step and map thses records with results are saved
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
     double precision(c_double) :: tpval(2),rtn
! svfunres: the values of state variable functions valid for this equilibrium
     double precision(c_double), dimension(:), allocatable :: svfunres
! the experiments are used in assessments and stored like conditions 
! lastcondition: link to condition list
! lastexperiment: link to experiment list
     TYPE(gtp_condition), pointer :: lastcondition,lastexperiment
! components and conversion matrix from components to elements
! complist: array with components
! compstoi: stoichiometric matrix of compoents relative to elements
! invcompstoi: inverted stoichiometric matrix
     TYPE(gtp_components), dimension(:), allocatable :: complist
     double precision(c_double), dimension(:,:), allocatable :: compstoi
     double precision(c_double), dimension(:,:), allocatable :: invcompstoi
! one record for each phase+composition set that can be calculated
! phase_varres: here all calculated data for the phase is stored
     TYPE(gtp_phase_varres), dimension(:), allocatable :: phase_varres
! index to the tpfun_parres array is the same as in the global array tpres 
! eq_tpres: here local calculated values of TP functions are stored
     TYPE(tpfun_parres), dimension(:), pointer :: eq_tpres
! current values of chemical potentials stored in component record but
! duplicated here for easy acces by application software
     double precision(c_double), dimension(:), allocatable :: cmuval
! xconc: convergence criteria for constituent fractions and other things
     double precision(c_double) :: xconv
! delta-G value for merging gridpoints in grid minimizer
! smaller value creates problem for test step3.BMM, MC and austenite merged
     double precision(c_double) :: gmindif=-5.0D-2
! maxiter: maximum number of iterations allowed
     integer(c_int) :: maxiter
! this is to save a copy of the last calculated system matrix, needed
! to calculate dot derivatives, initiate to zero
     integer(c_int) :: sysmatdim=0,nfixmu=0,nfixph=0
     integer(c_int), allocatable :: fixmu(:)
     integer(c_int), allocatable :: fixph(:,:)
     double precision(c_double), allocatable :: savesysmat(:,:)
  END TYPE gtp_equilibrium_data

!
! access to main OC library for equilibrium calculations
! using C language
!  use liboceq
!
  implicit none
!
!
contains
!
!\begin{verbatim}
  subroutine c_tqini(n, ceq) bind(c, name='c_tqini')
    integer(c_int), intent(in) :: n
!    type(gtp_equilibrium_data), pointer, intent(out) :: ceq 
    type(c_ptr), pointer, intent(out) :: ceq 
!\end{verbatim}  
  call tqini(n, ceq)
  end subroutine c_tqini

!\begin{verbatim}
  subroutine c_tqrfil(filename,ceq) bind(c)
    character(kind=c_char), intent(in) :: filename
    type(gtp_equilibrium_data), pointer, intent(in) :: ceq
!\end{verbatim}
    call tgrfil(filename, ceq)
  end subroutine c_tqrfil
!\begin{verbatim}
  subroutine c_tqrpfil(filename,nel,selel,ceq) bind(c)
    character(kind=c_char), intent(in) :: filename
    integer(c_int), intent(in) :: nel
    character(c_char), intent(in) :: selel
    type(gtp_equilibrium_data), pointer, intent(in) :: ceq
!\end{verbatim}
    call tqrpfil(filename, nel,selel,ceq)
  end subroutine c_tqrpfil
 
!\begin{verbatim}
  subroutine c_tqgcom(n,components,ceq) bind(c)
! get system components
    integer(c_int), intent(out) :: n
    character(kind=c_char, len=24), dimension(*), intent(out) :: &
    components
  type(gtp_equilibrium_data), pointer, intent(in) :: ceq
!\end{verbatim}
  end subroutine c_tqgcom

!\begin{verbatim}
  subroutine c_tqgnp(n,ceq) bind(c)
! get name of phase n,
! NOTE: n is phase number not extended phase index
    integer(c_int), intent(out) :: n 
    type(qtp_equilibrium_data), intent(in) :: ceq
!\end{verbatim}
  end subroutine c_tqgpn
end module liboctq_interface


