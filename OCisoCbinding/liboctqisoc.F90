! 
!
! Part of iso-C binding for OC TQlib from Teslos
! modified by Matthias Stratmann, Christophe Sigli,
! and Bo Sundman
!
! Update proposed by Romain Le Tellier and Clément Introïni
!
MODULE cstr
!
! convert characters from Fortran to C and vice versa
contains
  function c_to_f_string(s) result(str)
    use iso_c_binding
    implicit none
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
  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

  subroutine c_to_f_str(s,sty)
    use iso_c_binding
    implicit none
    character(kind=c_char,len=1), intent(in) :: s(*)
	character(len=24), intent(out) :: sty
    character(len=:), allocatable :: str
	
    integer i, nchars
    i = 1
    do
       if (s(i) == c_null_char) exit
       i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(len=nchars) :: str)
    sty = transfer(s(1:nchars), str)
	deallocate (str)
  end subroutine c_to_f_str

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

  subroutine f_to_c_string(fstring, cstr)
    use iso_c_binding
    implicit none
    character(len=24) :: fstring
    character(kind=c_char, len=1), intent(out) :: cstr(*)
    integer i
    do i = 1, len(fstring)
       cstr(i) = fstring(i:i)
       cstr(i+1) = c_null_char
    end do
  end subroutine f_to_c_string
  
END MODULE cstr

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
!
! module liboctqisoc
!
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

MODULE liboctqisoc
! 
! OCTQlib with iso-C binding
!
  use iso_c_binding
  use cstr
  use liboctq
!  use general_thermodynamic_package
  implicit none


  integer(c_int), bind(c) :: c_niter=-1

  integer(c_int), bind(c) :: c_nel=-1
  integer(c_int), bind(c) ::c_maxc=40, c_maxp=500
  type(c_ptr), bind(c), dimension(maxc) :: c_cnam
  character(len=25), dimension(maxc), target :: cnames
  real(c_double), bind(c), dimension(maxc) :: c_mass
  
  integer(c_int), bind(c) :: c_ntup
   
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
     character(c_char) :: eqname(24)
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
     real(c_double) :: gmindif
! maxiter: maximum number of iterations allowed
     integer(c_int) :: maxiter     
!CCI
! conv_iter: number of iterations reached after the equilibrium calculation
     integer(c_int) :: conv_iter
!CCI
! this is to save a copy of the last calculated system matrix, needed
! to calculate dot derivatives, initiate to zero
     integer(c_int) :: sysmatdim=0,nfixmu=0,nfixph=0
     integer(c_int) :: fixmu
     integer(c_int) :: fixph
     real(c_double) :: savesysmat
  END TYPE c_gtp_equilibrium_data

contains

! functions
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
  integer function c_noofcs(iph) bind(c, name='c_noofcs')
    integer(c_int), value :: iph
    c_noofcs = noofcs(iph)
    return 
  end function c_noofcs
  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
  integer function c_noconst(iph,ics,c_ceq) bind(c, name='c_noconst')
    integer(c_int), intent(in), value :: iph
    integer(c_int), intent(in), value :: ics
    type(c_ptr), intent(inout)  :: c_ceq
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    c_noconst = noconst(iph,ics,ceq)
    nullify(ceq)
    return 
  end function c_noconst

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

  subroutine examine_gtp_equilibrium_data(c_ceq) &
       bind(c, name='examine_gtp_equilibrium_data')
    type(c_ptr), intent(in), value :: c_ceq
    type(gtp_equilibrium_data), pointer :: ceq
    integer :: i,j
    call c_f_pointer(c_ceq, ceq)
    write(*,10) ceq%status, ceq%multiuse, ceq%eqno
10  format(/'gtp_equilibrium_data: status, multiuse, eqno, next'/, 3i4)
    write(*,20) ceq%eqname
20  format(/'Name of equilibrium'/,a)
    write(*,30) ceq%tpval, ceq%rtn
30  format(/'Value of T and P'/, 2f8.3, /'R*T'/, f8.4)
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

!CCI
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
!\begin{verbatim}
! Get the stoichiometric factor of an element in a species given by name
! species_name: name of the species (input character)
! iel: index of the element (input intger)
! el_name: name of the element (output character)
! stoi: value of the stoichiometric coefficient (output real)
!
  subroutine get_stoichiometric_coef(species_name, iel, el_name, c_stoi) bind(c, name='c_get_stoichiometry')
    character(kind=c_char), intent(in) :: species_name
    integer(c_int), intent(in), value :: iel
    character(kind=c_char), intent(inout) :: el_name(24)
    real(c_double), intent(inout) ::  c_stoi
!\end{verbatim}
    integer :: loksp
    character(len=:), allocatable :: eq_species_name
    character :: fstring*24

    ! Get the index of the species by its name
    eq_species_name = c_to_f_string(species_name)
    call find_species_record_exact(eq_species_name,loksp)
    if(gx%bmperr.ne.0) goto 1000

    ! Get the stoichiometric coefficient
    call get_stoichiometry(loksp, iel, fstring, c_stoi)
    call f_to_c_string(fstring, el_name)

1000 continue
    deallocate(eq_species_name)
    return
  end subroutine get_stoichiometric_coef
!CCI

!CCI
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
!\begin{verbatim}

! Change the stoichiometric factor of a species given by name
! species_name: name of the species (input character)
! new_stoi: new value of the stoichiometric coefficient (input real)
!
  subroutine change_stoichiometric(species_name,new_stoi) bind(c, name='c_change_stoichiometric')
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
    character(kind=c_char), intent(in) :: species_name
    real(c_double), intent(in), value :: new_stoi
!\end{verbatim}
    integer :: loksp
    character(len=:), allocatable :: eq_species_name

    eq_species_name = c_to_f_string(species_name)

    ! Get the index of the species by its name
    call find_species_record(eq_species_name,loksp)
    if(gx%bmperr.ne.0) goto 1000
    ! Change the stoichiometric factor (new_stoi) of the loksp-th species
    call set_new_stoichiometry(loksp,new_stoi)

1000 continue
    deallocate(eq_species_name)
    return
  end subroutine change_stoichiometric
!CCI

!\begin{verbatim}
  subroutine getelem()
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
!
! When an external unformatted file is read without reading any data base
! it is necessary to get the name of each component for setting new conditions
! The number of elements is then updated ut also the number of phase tuples
! for possible using before doing another equilibrium calculation
    integer iz
    character elname*2,name*24,refs*24
    double precision a1,a2,a3

! by default, c_nel is initialized to -1
    if(c_nel.lt.0) then
        nel=noel()
        c_nel=noel()
        do iz=1,nel
            call get_element_data(iz,elname,name,refs,a1,a2,a3)
            cnam(iz)=elname
            cnames(iz)=trim(elname) // c_null_char
            c_cnam(iz) = c_loc(cnames(iz))
         enddo
        ntup=nooftup()
        c_ntup=nooftup()
    endif
!
  end subroutine getelem
!\end{verbatim}

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqini(n, c_ceq) bind(c, name='c_tqini')
    integer(c_int), intent(in) :: n
    type(c_ptr), intent(out) :: c_ceq
!\end{verbatim}  
    type(gtp_equilibrium_data), pointer :: ceq
    integer :: i1,i2
   
    call tqini(n, ceq)
    c_ceq = c_loc(ceq)
    
    nullify(ceq)
    return 
	
  end subroutine c_tqini

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
!\begin{verbatim}
  subroutine c_tqvalpfu(phtupx, c_molepfu, c_masspfu, c_napfu, c_ncpfu, c_ceq) bind(c, name='c_tqvalpfu')
! get the number of moles, the mass of components and the number of components/atoms per Formula Units of phase n,
! NOTE: n is phase number, not extended phase index
    integer(c_int), intent(in), value :: phtupx   ! IN: index in phase tuple array
    real(c_double), intent(inout) :: c_molepfu    ! INOUT: moles per FU
    real(c_double), intent(inout) :: c_masspfu    ! INOUT: mass of components per FU
    real(c_double), intent(inout) :: c_ncpfu      ! INOUT: number of components per FU
    real(c_double), intent(inout) :: c_napfu      ! INOUT: number of atoms per FU
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fstring
    double precision :: napfu
    character name*(24)

    call c_f_pointer(c_ceq, ceq)

    c_ncpfu=ceq%phase_varres(phasetuple(phtupx)%lokvares)%abnorm(1)
    c_masspfu=ceq%phase_varres(phasetuple(phtupx)%lokvares)%abnorm(2)
    c_napfu=ceq%phase_varres(phasetuple(phtupx)%lokvares)%abnorm(3)
    c_molepfu=ceq%phase_varres(phasetuple(phtupx)%lokvares)%amfu

    c_ceq = c_loc(ceq)
    nullify(ceq)
  end subroutine c_tqvalpfu
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!  

!\begin{verbatim}
  subroutine c_tqrfil(filename,c_ceq) bind(c, name='c_tqrfil')
    character(kind=c_char,len=1), intent(in) :: filename(*)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
	type(gtp_equilibrium_data), pointer :: ceq
	character(len=:), allocatable :: fstring
    integer :: i,j,l
    character(kind=c_char, len=1),dimension(24), target :: f_pointers
! convert type(c_ptr) to fptr
    call c_f_pointer(c_ceq, ceq)
    fstring = c_to_f_string(filename)
!CCI : turn off warnings from reading the TDB file
    !call readtdbsilent
!CCI
    call tqrfil(fstring, ceq)
! after tqrfil ntup variable is defined
    c_ntup = ntup
    c_nel = nel
    do i = 1, nel
       cnames(i) = trim(cnam(i)) // c_null_char
       c_cnam(i) = c_loc(cnames(i))
	   c_mass(i)=cmass(i)
	   write(*,*) cmass(i)
    end do
    c_ceq = c_loc(ceq)
    deallocate(fstring)
    nullify(ceq)
  end subroutine c_tqrfil

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
  
!\begin{verbatim}
  subroutine c_tqrpfil(filename,nel,c_selel,c_ceq) bind(c, name='c_tqrpfil')
!change   
    character(kind=c_char), intent(in) :: filename
    integer(c_int), intent(in), value :: nel
    type(c_ptr), intent(in), dimension(nel), target :: c_selel
    type(c_ptr), intent(inout) :: c_ceq  
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=:), allocatable :: fstring
    character, pointer :: selel(:)
    integer :: i
    character elem(nel)*2
    fstring = c_to_f_string(filename)
    call c_f_pointer(c_ceq, ceq)
! convert the c type selel strings to f-selel strings
! note: additional character is for C terminated '\0'
    do i = 1, nel
       call c_f_pointer(c_selel(i), selel, [3])
       elem(i) = c_to_f_string(selel)
    end do
!CCI : turn off warnings from reading the TDB file
    call readtdbsilent
!CCI
    call tqrpfil(fstring, nel, elem, ceq)
! after tqrpfil ntup variable is defined
    c_ntup = ntup
    c_nel = nel
    do i = 1, nel
       cnames(i) = trim(cnam(i)) // c_null_char
       c_cnam(i) = c_loc(cnames(i))
	   c_mass(i)=cmass(i)
    end do
    c_ceq = c_loc(ceq)
    deallocate (fstring)
    nullify(ceq)
  end subroutine c_tqrpfil

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgcom(n,components,c_ceq) bind(c, name='c_tqgcom')
! get system components
    integer(c_int), intent(inout) :: n
    !character(kind=c_char, len=24), dimension(24), intent(out) :: c_components
    type(c_ptr), intent(inout) :: c_ceq  
!\end{verbatim}
    integer, target :: nc
    character(len=24) :: fcomponents(maxel)
    character(kind=c_char, len=1), dimension(maxel*24) :: components
    type(gtp_equilibrium_data), pointer :: ceq  
    integer :: i,j,l
    call c_f_pointer(c_ceq, ceq)
    call tqgcom(nc, fcomponents, ceq)
! convert the F components strings to C 
    l = len(fcomponents(1))
    do i = 1, nc
       do j = 1, l
          components((i-1)*l+j)(1:1) = fcomponents(i)(j:j)
       end do
! null termination
       components(i*l) = c_null_char 
    end do
    c_ceq = c_loc(ceq)
    n = nc
    
    nullify(ceq)
    return 
  end subroutine c_tqgcom

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgnp(n, c_ceq) bind(c, name='c_tqgnp')
    integer(c_int), intent(inout) :: n
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    call tqgnp(n, ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqgnp

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgpn(n,phasename, c_ceq) bind(c, name='c_tqgpn')
! get name of phase n,
! NOTE: n is phase number, not extended phase index
    integer(c_int), intent(in), value :: n
    character(kind=c_char, len=1), intent(inout) :: phasename(24)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fstring
    integer :: i
    call c_f_pointer(c_ceq, ceq)
! fstring = c_to_f_string(phasename)
    call tqgpn(n, fstring, ceq)
! copy the f-string to c-string and end with '\0'
    call f_to_c_string(fstring, phasename)
!    do i=1,len(trim(fstring))
!       phasename(i)(1:1) = fstring(i:i)
!       phasename(i+1)(1:1) = c_null_char
!    end do
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_tqgpn 

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgpi(n,phasename,c_ceq) bind(c, name='c_tqgpi')
! get index of phase phasename
    integer(c_int), intent(out) :: n
    character(c_char), intent(in) :: phasename(24)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fstring
    call c_f_pointer(c_ceq, ceq)
    fstring = c_to_f_string(phasename)
    call tqgpi(n, fstring, ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_tqgpi

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgpi2(iph,ics,phasename,c_ceq) bind(c, name='c_tqgpi2')
! get index of phase phasename
    integer(c_int), intent(out) :: iph
    integer(c_int), intent(out) :: ics
    character(c_char), intent(in) :: phasename(24)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fstring
    call c_f_pointer(c_ceq, ceq)
    fstring = c_to_f_string(phasename)
    call tqgpi2(iph, ics, fstring, ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_tqgpi2

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgpcn2(n, c, csname) bind(c, name='c_tqgpcn2')
! get name of constituent c in phase n
    integer(c_int), intent(in), value :: n  ! phase number
    integer(c_int), intent(in), value :: c  ! extended constituent index: 
!                                             10*species_number + sublattice
    character(kind=c_char, len=1), intent(inout) :: csname(24)
!\end{verbatim}
    character(len=24) :: fstring
    integer :: i
    call tqgpcn2(n,c,fstring)
    call f_to_c_string(fstring, csname)
! copy the f-string to c-string and end with '\0'
!   do i=1,len(trim(fstring))
!      csname(i)(1:1) = fstring(i:i)
!      csname(i+1)(1:1) = c_null_char
!   end do
  end subroutine c_tqgpcn2

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgpci(n,c, constituentname, c_ceq) bind(c, name='c_tqgpci')
! get index of constituent with name in phase n
    integer(c_int), intent(in) :: n 
    integer(c_int), intent(out) :: c ! exit: extended constituent index:
!                                      10*species_number+sublattice
    character(c_char), intent(in) :: constituentname(24)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fstring
    fstring = c_to_f_string(constituentname)
    call c_f_pointer(c_ceq, ceq)
    call tqgpci(n, c, fstring, ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqgpci

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgpcs(c,nspel,ielno,stoi,smass,qsp) bind(c, name='c_tqgpcs')
!get stoichiometry of constituent c in phase n
    integer(c_int), intent(in), value :: c ! in: extended constituent index:
!                                     10*species_number + sublattice
    integer(c_int), intent(out) :: nspel 
    
    integer(c_int), intent(out) :: ielno(*)
    real(c_double), intent(out) :: stoi(*) ! exit: stoichiometry of elements
    real(c_double), intent(out) :: smass     ! exit: total mass
    real(c_double), intent(out) :: qsp   
!\end{verbatim}
    call tqgpcs(c,nspel,ielno,stoi,smass,qsp)
  end subroutine c_tqgpcs

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgccf(n1,n2,elnames,stoi,mass,c_ceq)
! get stoichiometry of component n1
! n2 is number of elements ( dimension of elements and stoi )
    integer(c_int), intent(in) :: n1  ! in: component number
    integer(c_int), intent(out) :: n2 ! exit: number of elements in component
    character(c_char), intent(out) :: elnames(2) ! exit: element symbols
    real(c_double), intent(out) :: stoi(*) ! exit: element stoichiometry
    real(c_double), intent(out) :: mass    ! exit: component mass
!                                           (sum of element mass)
    type(c_ptr), intent(inout) :: c_ceq  
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    call tqgccf(n1,n2,elnames,stoi, mass, ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqgccf

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgnpc(n,c,c_ceq) bind(c, name='c_tqgnpc')
! get number of constituents of phase n
    integer(c_int), intent(in) :: n ! in: phase number 
    integer(c_int), intent(out) :: c ! exit: number of constituents
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq,ceq)
    call tqgnpc(n,c,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqgnpc

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqphsts(phtupx,newstat,val,c_ceq) bind(c, name='c_tqphsts')
! set status of phase tuple: SUSPEND, DORMANT, ENTERED, FIX
    integer(c_int), intent(in), value :: phtupx
    integer(c_int), intent(in), value :: newstat
    real(c_double), intent(in) :: val
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq,ceq)
    call tqphsts(phtupx,newstat,val,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_tqphsts

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqphsts2(phnames,newstat,val,c_ceq) bind(c, name='c_tqphsts2')
! set status of phase tuple: SUSPEND, DORMANT, ENTERED, FIX
    character(c_char), intent(in) :: phnames
    integer(c_int), intent(in), value :: newstat
    real(c_double), intent(in) :: val
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=:), allocatable :: fphnames
    call c_f_pointer(c_ceq,ceq)
    fphnames = c_to_f_string(phnames)
    call tqphsts2(fphnames,newstat,val,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    deallocate(fphnames)
    return
  end subroutine c_tqphsts2

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqsetc(statvar, n1, n2, mvalue, cnum, c_ceq) &
       bind(c, name='c_tqsetc')
! set condition
! stavar is state variable as text
! n1 and n2 are auxilliary indices
! value is the value of the condition
! cnum is returned as an index of the condition.
! to remove a condition the value sould be equial to RNONE ????
! when a phase indesx is needed it should be 10*nph + ics
! SEE TQGETV for doucumentation of stavar etc.
!>>>> to be modified to use phase tuplets
    integer(c_int), intent(in),value :: n1 !in: 0 or extended phase index:
!                                       10*phase_number+comp.set
                                     ! or component set
    integer(c_int), intent(in),value :: n2 !
    integer(c_int), intent(out) :: cnum !exit: 
!                                        sequential number of this condition
    character(c_char), intent(in) :: statvar !in: character
!                                             with state variable symbol
    real(c_double), intent(in),value :: mvalue  !in: value of condition
   
    type(c_ptr), intent(in) :: c_ceq ! in: current equilibrium
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=:), allocatable :: fstatvar
!
    call c_f_pointer(c_ceq, ceq)
    fstatvar = c_to_f_string(statvar)
    call tqsetc(fstatvar, n1, n2, mvalue, cnum, ceq)
    nullify(ceq)
    deallocate(fstatvar)
  end subroutine c_tqsetc
!\end{verbatim}

!\begin{verbatim}
  subroutine c_tqcalc(c_ceq,mode) bind(c,name='c_tqcalc')
! calculate equilibrium with different methods
! mode=0 means calculate without grid minimizer
! mode=1 means start values using global gridminimization
! mode=2 means calculate carefully (default)
    integer(c_int), intent(in),value :: mode
    integer n
    logical confirm

    double precision, allocatable, dimension(:) ::xknown,aphl,cmu
    integer, allocatable, dimension(:) ::iphl,icsl,nyphl
    integer nv
    double precision totam

    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq,ceq)

! mode=0 means calculate without global minimizer
    if(mode.eq.0) then
       call calceq2(mode,ceq)
! mode=1 means start values using global gridminimization
    elseif(mode.eq.1) then
       call calceq2(mode,ceq)
! mode=2 means calculate carefully (default)
    else
    ! first parameter 0 means bosses_method, 1 means carefully
       n=1
       !
       allocate(xknown(noel()))
       xknown(:)=0.0
       allocate(aphl(noel()))
       aphl(:)=0.0
       allocate(cmu(noel()))
       cmu(:)=0.0
       allocate(iphl(noel()))
       iphl(:)=0
       allocate(icsl(noel()))
       icsl(:)=0
       allocate(nyphl(noel()))
       nyphl(:)=0
       !
       call extract_massbalcond(ceq%tpval,xknown,totam,ceq)
       if(gx%bmperr.eq.0) then
            call global_gridmin(1,ceq%tpval,xknown,nv,iphl,icsl,&
            aphl,nyphl,cmu,ceq)
          if(gx%bmperr.eq.0) then
             call calculate_carefully(n,ceq)
          endif
       endif
       !
       deallocate(xknown)
       deallocate(aphl)
       deallocate(cmu)
       deallocate(iphl)
       deallocate(icsl)
       deallocate(nyphl)
    endif

    if(gx%bmperr.ne.0) goto 1000

    ntup=nooftup()
    c_ntup=nooftup()
    c_niter=ceq%conv_iter

1000 continue
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqcalc


!\begin{verbatim}
  subroutine c_tqce(mtarget,n1,n2,mvalue,c_ceq) bind(c,name='c_tqce')
! calculate equilibrium with possible target
! Target can be empty or a state variable with indicies n1 and n2
! value is the calculated value of target
    integer(c_int), intent(in),value :: n1
    integer(c_int), intent(in),value :: n2
    type(c_ptr), intent(inout) :: c_ceq
    character(c_char), intent(inout) :: mtarget  
    real(c_double), intent(inout) :: mvalue
!\end{verbatim}
	type(gtp_equilibrium_data), pointer :: ceq
    character(len=:), allocatable :: fstring
    call c_f_pointer(c_ceq,ceq)
    fstring = c_to_f_string(mtarget)
    call tqce(fstring,n1,n2,mvalue,ceq)
    if(gx%bmperr.ne.0) goto 1000
    c_ntup=ntup
!CCI    
    c_niter=ceq%conv_iter
!CCI
1000 continue
    c_ceq = c_loc(ceq)
    deallocate(fstring)
    nullify(ceq)
    return 
  end subroutine c_tqce


!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\
!\begin{verbatim}
  subroutine c_tqdceq(ceqname) bind(c,name='c_tqdceq')
    character(kind=c_char), intent(in) :: ceqname
! delete equilibrium with name
!\end{verbatim}
    character(len=:), allocatable :: name
    integer n1
    name = c_to_f_string(ceqname)
    call tqdceq(name)
    deallocate(name)
  end subroutine c_tqdceq


!\begin{verbatim}
  subroutine c_tqfree() bind(c, name='c_tqfree')
!\end{verbatim}
    integer intv(10)
    double precision dblv(10)
    call deallocate_gtp(intv,dblv)
  end subroutine c_tqfree



!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgetv(statvar,n1,n2,n3,values,c_ceq) bind(c,name='c_tqgetv')
! get equilibrium results using state variables
! stavar is the state variable IN CAPITAL LETTERS with indices n1 and n2 
! n3 at the call is the dimension of values, changed to number of values
! value is the calculated value, it can be an array with n3 values.
    implicit none
    integer(c_int), intent(in),value ::  n1,n2
    integer(c_int), intent(inout) :: n3
    character(c_char), intent(in) :: statvar
    real(c_double), intent(inout) :: values(*)
    type(c_ptr), intent(inout) :: c_ceq  !IN: current equilibrium
!========================================================
! >>>> implement use of phase tuples 
! stavar must be a symbol listed below
! IMPORTANT: some terms explained after the table
! Symbol  index1,index2                     Meaning (unit)
!.... potentials
! T     0,0                                             Temperature (K)
! P     0,0                                             Pressure (Pa)
! MU    component,0 or phase-tuple*1,constituent*2  Chemical potential (J)
! AC    component,0 or phase-tuple,constituent      Activity = EXP(MU/RT)
! LNAC  component,0 or phase-tuple,constituent      LN(activity) = MU/RT
!...... extensive variables
! U     0,0 or phase-tuple,0       Internal energy (J) whole system or phase
! UM    0,0 or phase-tuple,0       same per mole components
! UW    0,0 or phase-tuple,0       same per kg
! UV    0,0 or phase-tuple,0       same per m3
! UF    phase-tuple,0              same per formula unit of phase
! S*3   0,0 or phase-tuple,0       Entropy (J/K) 
! V     0,0 or phase-tuple,0       Volume (m3)
! H     0,0 or phase-tuple,0       Enthalpy (J)
! A     0,0 or phase-tuple,0       Helmholtz energy (J)
! G     0,0 or phase-tuple,0       Gibbs energy (J)
! ..... some extra state variables
! NP    phase-tuple,0              Moles of phase
! BP    phase-tuple,0              Mass of moles (kg)
! Q     phase-tuple,0              Internal stability/RT (dimensionless)
! DG    phase-tuple,0              Driving force/RT (dimensionless)
!....... amounts of components
! N     0,0 or component,0 or phase-tuple,component   Moles of component
! X     component,0 or phase-tuple,component          Mole fraction of component
! B     0,0 or component,0 or phase-tuple,component   Mass of component
! W     component,0 or phase-tuple,component          Mass fraction of component
! Y     phase-tuple,constituent*1                     Constituent fraction
!........ some parameter identifiers
! TC    phase-tuple,0              Magnetic ordering temperature
! BMAG  phase-tuple,0              Aver. Bohr magneton number
! MQ&   phase-tuple,constituent    Mobility
! THET  phase-tuple,0              Debye temperature
! LNX   phase-tuple,0              Lattice parameter
! EC11  phase-tuple,0              Elastic constant C11
! EC12  phase-tuple,0              Elastic constant C12
! EC44  phase-tuple,0              Elastic constant C44
!........ NOTES:
! *1 The phase-tuple is   is structure with 2 integers: phase and comp.set
! *2 The constituent index is 10*species_number + sublattice_number
! *3 S, V, H, A, G, NP, BP, N, B and DG can have suffixes M, W, V, F also
!--------------------------------------------------------------------
! special addition for TQ interface: d2G/dyidyj
! D2G + extended phase index
!------------------------------------
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fstring
    integer :: n
    integer :: i

    call c_f_pointer(c_ceq, ceq)
!    call list_conditions(6,ceq)
!    call list_phase_results(1,1,0,6,ceq)
!    write(*,*)'Phase and error code: ',1,gx%bmperr
!    call list_phase_results(2,1,0,6,ceq)
!    write(*,*)'Phase and error code: ',2,gx%bmperr
!    write(*,*)

    call c_to_f_str(statvar,fstring)

    call tqgetv(fstring, n1, n2, n3, values, ceq)
! debug ...
!   write(*,55)fstring(1:len_trim(fstring)),n1,n2,n3,(values(i),i=1,n3)
!55  format(/'From c_tqgetv: ',a,': ',3i3,6(1pe12.4))
!    write(*,*)
! end debug
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
	
  end subroutine c_tqgetv

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgnsubl(n1,nsub,c_ceq) bind(c,name='c_tqgnsubl')
! This subroutine returns the number of sublattices (1 if no sublattices)
! of phase identified by its phase tuple index
    implicit none
    integer(c_int), intent(in), value :: n1
    integer(c_int), intent(out) :: nsub
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    call get_sublattice_number(phasetuple(n1)%ixphase,nsub,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqgnsubl

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!


!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgsubstruc(n1,nsub,nkl,nsites,c_ceq) bind(c,name='c_tqgsubstruc')
! This subroutine returns structures of each sublattice
! of phase identified by its phase tuple and composition set indexes
! (number of constituents in each sublattice and number of sites)
    implicit none
    integer(c_int), intent(in), value :: n1,nsub
    integer(c_int), intent(out), dimension(nsub) :: nkl
    real(c_double), intent(out), dimension(nsub) :: nsites
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)

    call get_sublattice_structure(phasetuple(n1)%ixphase,phasetuple(n1)%compset,nsub,nkl,nsites,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqgsubstruc

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgconsdata(n1,icons,yarr,charge,csname,ncel,c_ceq) bind(c,name='c_tqgconsdata')
! This subroutine returns mole fraction, charge and name of constituent
! of phase identified by its phase tuple and composition set indexes
! (index of constituents and number of sites)
    implicit none
    integer(c_int), intent(in), value :: n1,icons
    real(c_double), intent(inout) :: yarr
    integer(c_int), intent(inout) :: charge, ncel
    character(kind=c_char, len=1), intent(inout) :: csname(24)
    type(c_ptr), intent(inout) :: c_ceq
    character(len=24) :: fstring
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    call get_constituent_data(phasetuple(n1)%ixphase,phasetuple(n1)%compset,icons,yarr,charge,&
                              fstring,ncel,ceq)
    call f_to_c_string(fstring, csname)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqgconsdata

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqgphc1(n1,nsub,cinsub,spix,yfrac,sites,extra,c_ceq)&
 bind(c,name='c_tqgphc1')
! tq_get_phase_constitution
! This subroutine returns the sublattices and constitution of a phase
! n1 is phase tuple index
! nsub is the number of sublattices (1 if no sublattices)
! cinsub is an array with the number of const\EDtuents in each sublattice
! spix is an array with the species index of the constituents in all sublattices
! sites is an array of the site ratios for all sublattices.  
! yfrac is the constituent fractions in same order as in spix
! extra is an array with some extra values: 
!    extra(1) is the number of moles of components per formula unit
!    extra(2) is the net charge of the phase
    implicit none
    !integer n1,nsub,cinsub(*),spix(*)
    integer(c_int), intent(in), value :: n1
    integer(c_int), intent(out) :: nsub
    integer(c_int), intent(out) :: cinsub(*)
    integer(c_int), intent(in) :: spix(*)
    !double precision sites(*),yfrac(*),extra(*)
    real(c_double), intent(in) :: sites(*)
    real(c_double), intent(in) :: yfrac(*)
    real(c_double), intent(in) :: extra(*)
    !type(gtp_equilibrium_data), pointer :: ceq
    type(c_ptr), intent(inout) :: c_ceq  
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    !call tqgphc1(n1,nsub2,cinsub2,spix2,yfrac2,sites2,extra2,ceq)
    call tqgphc1(n1,nsub,cinsub,spix,yfrac,sites,extra,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqgphc1

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqsphc1(n1,yfra,extra,c_ceq) bind(c,name='c_tqsphc1')
! tq_set_phase_constitution
! To set the constitution of a phase
! n1 is phase tuple index
! yfra is an array with the constituent fractions in all sublattices
! in the same order as obtained by tqgphc1
! extra is an array with returned values with the same meaning as in tqgphc1
! NOTE The constituents fractions are normallized to sum to unity for each
!      sublattice and extra is calculated by tqsphc1
! T and P must be set as conditions.
    implicit none
    integer(c_int), intent(in), value :: n1
    real(c_double), intent(in) ::yfra(*)
    real(c_double), intent(out) :: extra(*)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    call set_constitution(phasetuple(n1)%ixphase,phasetuple(n1)%compset,&
         yfra,extra,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqsphc1

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqcph1(n1,n2,n3,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,c_ceq) &
       bind(c,name='c_tqcph1')
! tq_calculate_phase_properties
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! WARNIG: this is not a subroutine to calculate chemical potentials
! those can only be made by an equilibrium calculation.
! The values returned are partial derivatives of G for the phase at the
! current T, P and phase constitution.  The phase constitution has been
! obtained by a previous equilibrium calculation or 
! set by the subroutine tqsphc
! It corresponds to the "calculate phase" command.
!
! NOTE that values are per formula unit divided by RT, 
! divide also by extra(1) in subroutine tqsphc1 to get them per mole component
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! calculate G and some or all derivatives for a phase at current composition
! n1 is the phase tuple index
! n2 is 0 if only G and derivatves wrt T and P, 1 also first drivatives wrt 
!    compositions, 2 if also 2nd derivatives
! n3 is returned as number of constituents (dimension of returned arrays)
! gtp is an array with G, G.T, G:P, G.T.T, G.T.P and G.P.P
! dgdy is an array with G.Yi
! d2gdydt is an array with G.T.Yi
! d2gdydp is an array with G.P.Yi
! d2gdy2 is an array with the upper triangle of the symmetrix matrix G.Yi.Yj 
! reurned in the order:  1,1; 1,2; 1,3; ...           
!                             2,2; 2,3; ...
!                                  3,3; ...
! for indexing one can use the integer function ixsym(i1,i2)
    implicit none
    integer(c_int), intent(in), value :: n1
    integer(c_int), intent(in), value :: n2
    integer(c_int), intent(out) :: n3
    real(c_double), intent(out) :: gtp(6)
    real(c_double), intent(out) :: dgdy(*)
    real(c_double), intent(out) :: d2gdydt(*)
    real(c_double), intent(out) :: d2gdydp(*)
    real(c_double), intent(out) :: d2gdy2(*)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    call tqcph1(n1,n2,n3,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqcph1

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_tqcph3(n1,n2,g,c_ceq) bind(c,name='c_tqcph3')
    implicit none
    integer(c_int), intent(in), value :: n1
    integer(c_int), intent(in), value :: n2
    real(c_double), intent(out) :: g(*)
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    call tqcph3(n1,n2,g,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_tqcph3
  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

 !\begin{verbatim}  
  subroutine c_reset_conditions(cline,c_ceq) bind(c, name='c_reset_conditions')
    implicit none
    character(c_char), intent(in) :: cline(24) 
    type(c_ptr), intent(inout) :: c_ceq ! in: current equilibrium
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fstring
    fstring = c_to_f_string(cline)
    call c_f_pointer(c_ceq, ceq)
    
    call reset_conditions(fstring,ceq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return 
  end subroutine c_reset_conditions
  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_Change_Status_Phase(phasename,nystat,myval,c_ceq)&
       bind(c, name='c_Change_Status_Phase') 
!change the status Fixed or Entered of a phase 
!PHFIXED=2
!PHENTERED=0
    implicit none
    character(c_char), intent(in) :: phasename(24)
    integer(c_int), intent(in), value :: nystat
    real(c_double), intent(in),value :: myval
    type(c_ptr), intent(inout) :: c_ceq ! in: current equilibrium
!\end{verbatim}	
    type(gtp_equilibrium_data), pointer :: ceq 
    character(len=24) :: fstring
    call c_f_pointer(c_ceq, ceq)
    call c_to_f_str(phasename,fstring)
    call change_many_phase_status(fstring,nystat,myval,ceq)
!    call Change_Status_Phase(fstring,nystat,myval,ceq)
    c_ceq = c_loc(ceq)
    
1000 continue
    nullify(ceq)	
    return
  end subroutine c_Change_Status_Phase
  

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine c_Set_Reference_State(iel,c_phase,tpref,c_ceq)&
       bind(c, name='c_Set_Reference_State')
! set component reference state
    integer(c_int), intent(in), value :: iel
    character(c_char), intent(in) :: c_phase(24)
    real(c_double), intent(in) :: tpref(2)
    type(c_ptr), intent(inout) :: c_ceq
    character(len=24) :: phase
    type(gtp_equilibrium_data), pointer :: ceq
    integer phtupx
!\end{verbatim}

    call c_f_pointer(c_ceq, ceq)
    phase = c_to_f_string(c_phase)
    call find_phasetuple_by_name(phase,phtupx)
    if(gx%bmperr.ne.0) goto 1000
    call set_reference_state(iel,phtupx,tpref,ceq)
1000 continue
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_Set_Reference_State


!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}  
  subroutine c_List_Conditions(c_ceq)&
       bind(c, name='c_List_Conditions') 
!change the status Fixed or Entered of a phase 
!PHFIXED=2
!PHENTERED=0
    implicit none
	
    type(c_ptr), intent(inout) :: c_ceq ! in: current equilibrium
!\end{verbatim}	
    type(gtp_equilibrium_data), pointer :: ceq 
    call c_f_pointer(c_ceq, ceq)
    call list_conditions(6,ceq)
1000 continue
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_List_Conditions
  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim} 
  subroutine c_checktdb(tdbfile)&
       bind(c, name='c_checktdb') 
    character(kind=c_char), intent(in) :: tdbfile
!\end{verbatim}
    integer:: nel,i
    character selel(maxel)*2
    character(len=:), allocatable :: fstring
    character(len=:), allocatable :: ext
    ext='.tdb'
    fstring = c_to_f_string(tdbfile)
    call checkdb2(fstring,ext,nel,selel)
    c_nel = nel
    do i = 1, nel
       cnames(i) = trim(selel(i)) // c_null_char
       c_cnam(i) = c_loc(cnames(i))
    end do
    deallocate(fstring)
    return
  end subroutine c_checktdb
  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!
  
!\begin{verbatim} 
  subroutine c_newEquilibrium(ceqname,ieq) bind(c, name='c_newEquilibrium') 
    character(kind=c_char), intent(in) :: ceqname
    integer(c_int), intent(out):: ieq
!\end{verbatim}
    character(len=:), allocatable :: fstring
    fstring = c_to_f_string(ceqname)
    call enter_equilibrium(fstring,ieq)
    deallocate(fstring)
  end subroutine c_newEquilibrium


!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim} 
  subroutine c_tqcceq(ceqname,n1,c_newceq,c_ceq) &
       bind(c, name='c_tqcceq')
    character(kind=c_char), intent(in) :: ceqname
    integer(c_int), intent(out) :: n1
    type(c_ptr), intent(inout) :: c_newceq
    type(c_ptr), intent(in) :: c_ceq
!\end{verbatim}
    character(len=:), allocatable :: name

    type(gtp_equilibrium_data), pointer :: newceq,ceq
    call c_f_pointer(c_ceq, ceq)
    call c_f_pointer(c_newceq, newceq)
    name = c_to_f_string(ceqname)
    
    call tqcceq(name,n1,newceq,ceq)
    c_newceq=c_loc(newceq)
    deallocate(name)
    nullify(ceq)
  end subroutine c_tqcceq


!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim} 
  subroutine c_tqselceq(ceqname,c_ceq) &
       bind(c, name='c_tqselceq')
    character(kind=c_char), intent(in) :: ceqname
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    character(len=:), allocatable :: name
    type(gtp_equilibrium_data), pointer :: ceq
    call c_f_pointer(c_ceq, ceq)
    name = c_to_f_string(ceqname)
    call tqselceq(name,ceq)
    c_ceq=c_loc(ceq)
    deallocate(name)
    nullify(ceq)
    return
  end subroutine c_tqselceq

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine c_tqgdmat(phtupx,tpval,xknown,cpot,tyst,nend,mugrad,mobval,consnames,n1,c_ceq) &
       bind(c, name='c_tqgdmat') 
    
    integer(kind=c_int), intent(in), value :: phtupx
    real(kind=c_double), intent(in) :: tpval(2)
    real(kind=c_double), intent(in) :: xknown(*)
    real(kind=c_double), intent(out) :: cpot(*)
    integer(kind=c_int), intent(in), value :: tyst
    integer(kind=c_int), intent(out) :: nend
    real(kind=c_double), intent(out) :: mugrad(*)
    real(kind=c_double), intent(out) :: mobval(*)
    character(kind=c_char, len=1), intent(out), dimension(maxconst*24) :: consnames
    integer(kind=c_int), intent(out) :: n1
    type(c_ptr), intent(inout) :: c_ceq 
!\end{verbatim}
    logical btyst
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=24) :: fconsnames(maxconst)
    integer :: i,j,l
    call c_f_pointer(c_ceq, ceq)
    if (tyst.eq.1) then
       btyst=.TRUE.
       else 
       btyst=.FALSE.
    endif 
    call tqgdmat(phtupx,tpval,xknown,cpot,btyst,nend,mugrad,mobval,fconsnames,n1,ceq)
! convert the F fconsnames strings to C 
    l = len(fconsnames(1))
    do i = 1, n1
       do j = 1, l
          consnames((i-1)*l+j)(1:1) = fconsnames(i)(j:j)
       end do
! null termination
       consnames(i*l) = c_null_char 
    end do
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_tqgdmat
  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim} 
  subroutine c_copy_equilibrium(c_neweq,ceqname,c_ceq) &
       bind(c, name='c_copy_equilibrium') 
    type(c_ptr), intent(inout) :: c_neweq  
    character(kind=c_char), intent(in) :: ceqname
    type(c_ptr), intent(in) :: c_ceq  
!\end{verbatim}
    character(len=:), allocatable :: fstring
    type(gtp_equilibrium_data), pointer :: ceq
    type(gtp_equilibrium_data), pointer :: neweq
    call c_f_pointer(c_ceq, ceq)
    fstring = c_to_f_string(ceqname)
    call copy_equilibrium(neweq,fstring,ceq)
    c_neweq=c_loc(neweq)
    deallocate(fstring)
    nullify(ceq)
    return
  end subroutine c_copy_equilibrium

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_selecteq(ieq,c_ceq) bind(c, name='c_selecteq')
    integer(c_int), intent(in),value :: ieq
    type(c_ptr), intent(out) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
!call c_f_pointer(c_ceq, ceq)
!call selecteq(ieq,ceq)
    ceq=>eqlista(ieq)
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_selecteq

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_enter_svf(c_tpfun,c_ceq) bind(c, name='c_enter_svf')
! enter a state variable function like CP=H.T;
    character(kind=c_char), intent(in) :: c_tpfun
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=:), allocatable :: tpfun
    integer ip
    tpfun = c_to_f_string(c_tpfun)
    ip=1
    call c_f_pointer(c_ceq, ceq)
    call enter_svfun(tpfun,ip,ceq)
    !call evaluate_all_svfun_old(-1,ceq) ! mandatory ?
    c_ceq = c_loc(ceq)
    deallocate(tpfun)
    return
  end subroutine c_enter_svf

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_get_value_svf(c_tpfun,c_svfvalue,c_ceq) bind(c, name='c_get_value_svf')
! evaluate all state variable funtions
! actual_arg are names of phases, components or species as @Pi, @Ci and @Si
! (NOT IMPLEMENTED YET see minimizer/matsmin.F90 for more details)
! if mode=1 always evaluate, if mode=0 several options
    character(kind=c_char), intent(in) :: c_tpfun
    real(c_double), intent(inout) :: c_svfvalue
    type(c_ptr), intent(inout) :: c_ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: ceq
    character(len=:), allocatable :: tpfun
    character actual_arg(2)*16
    double precision value
    integer ip,mode
    tpfun = c_to_f_string(c_tpfun)
    call c_f_pointer(c_ceq, ceq)
    call capson(tpfun)
    call find_svfun(tpfun,ip)
    mode=1 ! always is evaluated (see minimizer/matsmin.F90 for more details)
    actual_arg = ' '
    c_svfvalue=meq_evaluate_svfun(ip,actual_arg,mode,ceq)
    c_ceq = c_loc(ceq)
    deallocate(tpfun)
    return
  end subroutine c_get_value_svf

  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}  
  subroutine c_set_grid_density(ngrid) bind(c, name='c_set_grid_density') 
!\end{verbatim}
    integer(c_int), intent(in), value :: ngrid
    if(ngrid.eq.0) then
       ! this set GSOGRID, small grid and clears GSXGRID
       globaldata%status=ibset(globaldata%status,GSOGRID)
       globaldata%status=ibclr(globaldata%status,GSXGRID)
       globaldata%status=ibclr(globaldata%status,GSYGRID)
 !      write(*,*)'Sparse grid set'
    elseif(ngrid.eq.1) then
       ! DEFAULT, all gridbits are cleared
       globaldata%status=ibclr(globaldata%status,GSXGRID)
       globaldata%status=ibclr(globaldata%status,GSOGRID)
       globaldata%status=ibclr(globaldata%status,GSYGRID)
 !      write(*,*)'Normal grid set'
    elseif(ngrid.eq.2) then
       ! set GSXGRID (and clear GSOGRID and GSYGRID)
       globaldata%status=ibclr(globaldata%status,GSOGRID)
       globaldata%status=ibset(globaldata%status,GSXGRID)
       globaldata%status=ibclr(globaldata%status,GSYGRID)
 !      write(*,*)'Dense grid set'
    elseif(ngrid.eq.3) then
       ! set GSYGRID (and clear GSXGRID and GSOGRID)
       globaldata%status=ibclr(globaldata%status,GSOGRID)
       globaldata%status=ibclr(globaldata%status,GSXGRID)
       globaldata%status=ibset(globaldata%status,GSYGRID)
 !      write(*,*)'Very dense grid set'
    else
       write(*,*)'Only level 0, 1, 2 and implemented'
    endif
    
    return
  end subroutine c_set_grid_density
  
!\begin{verbatim}  
  subroutine c_set_status_globaldata() bind(c, name='c_set_status_globaldata') 
!\end{verbatim}
!globaldata%status=ibclr(globaldata%status,GSADV)
!globaldata%status=ibclr(globaldata%status,GSNOPAR)
!globaldata%status=ibclr(globaldata%status,GSXGRID)
    globaldata%status=ibclr(globaldata%status,GSNOACS)
    
    return
  end subroutine c_set_status_globaldata

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  integer function c_errors_number() bind(c, name='c_errors_number')
!\end{verbatim}
    c_errors_number=0
    if(gx%bmperr.ne.0) then
       c_errors_number=gx%bmperr
    endif
    return
  end function c_errors_number

!\begin{verbatim} 	
  subroutine c_reset_errors_number() bind(c, name='c_reset_errors_number')
!\end{verbatim}	
    gx%bmperr=0
    return
  end subroutine c_reset_errors_number
  
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
  subroutine c_new_gtp() bind(c, name='c_new_gtp')
!\end{verbatim}
    call new_gtp
  end subroutine c_new_gtp

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
! Save OC environment in the c_filename at c_specification format
! c_specification : UNFORMATTED, DIRECT, TDB, MACRO or LaTeX
  subroutine c_gtpsave(c_filename,c_specification) bind(c, name='c_gtpsave')
    character(kind=c_char), intent(in) :: c_filename,c_specification
!\end{verbatim}
    character(len=:), allocatable :: filename,specification
    filename = c_to_f_string(c_filename)
    specification = c_to_f_string(c_specification)
    call gtpsaveu(filename,specification)
    deallocate(filename,specification)
  end subroutine c_gtpsave
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
! Read OC environment from the c_filename at c_specification format
! c_specification : UNFORMATTED, DIRECT, TDB, MACRO or LaTeX
  subroutine c_gtpread(c_filename,c_specification) bind(c, name='c_gtpread')
    character(kind=c_char), intent(in) :: c_filename,c_specification
!\end{verbatim}
    character(len=:), allocatable :: filename,specification
    filename = c_to_f_string(c_filename)
    specification = c_to_f_string(c_specification)
    call gtpread(filename,specification)
    deallocate(filename,specification)
!CCI
    call getelem()
!CCI
  end subroutine c_gtpread



!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine c_tqcheckphstab(is_stable,c_phtupx,c_ceq) bind(c, name='c_tqcheckphstab')
! check if a phase if stable
    implicit none
    logical(c_bool), intent(inout) :: is_stable
    integer(c_int), intent(in), value :: c_phtupx
    type(c_ptr), intent(inout) :: c_ceq
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer lokvares

    call c_f_pointer(c_ceq, ceq)
    lokvares=phasetuple(c_phtupx)%lokvares
    if(ceq%phase_varres(lokvares)%phstate.ge.phentstab) then
        is_stable = .TRUE.
    endif
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_tqcheckphstab


!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine c_tqlr(c_mode,c_ceq) bind(c, name='c_tqlr')
! list the equilibrium results like in OC
    implicit none
    type(c_ptr), intent(inout) :: c_ceq
    integer(c_int), intent(in), value :: c_mode
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer phtupx,iph,ics,lokvares,mode
    logical once

    call c_f_pointer(c_ceq, ceq)
    write(6,10)
10  format(/20('*')/'Start debug output from TQLR: ')
    call list_conditions(6,ceq)
    call list_global_results(6,ceq)
    call list_components_result(6,1,ceq)
    call list_all_elements(6)
    call list_all_species(6)
    call list_all_phases(6,ceq)
    once=.TRUE.
    mode=max(0,c_mode)

    do phtupx=1,nooftup()
       lokvares=phasetuple(phtupx)%lokvares
       if(ceq%phase_varres(lokvares)%phstate.ge.phentstab) then
          iph=phasetuple(phtupx)%ixphase
          ics=phasetuple(phtupx)%compset
          call list_phase_results(iph,ics,mode,lut,once,ceq)
       endif
    enddo
    write(6,20)
20  format('End debug output from TQLR'/20('*')/)
1000 continue
    c_ceq = c_loc(ceq)
    nullify(ceq)
    return
  end subroutine c_tqlr


!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine c_tqtgsw(i) bind(c, name='c_tqtgsw')
  integer(c_int), intent(in),value :: i
!\end{verbatim}
  call tqtgsw(i)  
  end subroutine c_tqtgsw

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\
	
end module liboctqisoc
