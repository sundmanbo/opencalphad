!
! Minimal TQ interface.
!
! To compile and link this with an application one must first compile
! and form a library with of the most OC subroutines (lib\liboceq.a)
! and copy this and the corresponding "liboceqplus.mod" file
! from this compilation to the folder with this library
!
! NOTE that for the identification of phase and composition sets this
! TQ interface use a Fortran TYPE called gtp_phasetuple containing two
! integers, "phase" with the phase number and "compset" with the
! comp.set The number of phase tuples is initially equal to the number
! of phases and have the same index.  This represent comp.set 1 of the
! phases as each phase has just one composition set.  A phase may have
! several comp.sets created by calculations or by commands and these will
! have phase tuple index higher than the number of phases and their index
! is in the order of which they were created.
! This may cause some problems if composition sets are deleted because that
! will change the phase tuple index for those with higher index.  So do not
! delete comp.sets or at least be very careful when deleting comp.sets
!
! 210328 BOS Tested
! 191101 BOS Updates some routines and added two dummy modules for C routines
! 181030 BOS Updates some routines
! 150520 BOS added a few subroutines for single phase data and calculations
! 141210 BOS changed to use phase tuples
! 140128 BOS added D2G and phase specific V and G
! 140128 BOS added possibility to calculate without invoking grid minimizer
! 140125 BOS Changed name to liboctq
! 140123 BOS Added ouput of MQ G, V and normalized
!------------------------------------------------------------
! subroutines and functions
! tqini    ok initiate
! tqrfil   ok read a database file
! tqrpfil  ok read specified elements from database file
! -------------------------
! tqgcom   ok get number of system components and their names
! tqgnp    ok get number of phase tuples (phases and comp. sets)
! tqgpn    ok get name of phase tuple
! tqgpi    ok get phase tuple index of phase using its name
! tqgpcn   -  get name of constituent of a phase using index
! tqgpci   -  get index of constituent of a phase using name
! tqgpcs   -  get stoichiometry of species as system components 
! tqgccf   -  get stoichiometry of system component as elements
! tqgnpc   -  get number of constituents in phase
! -------------------------
! tqcref  -  set reference state for component
! tqphsts  ok set status of phase tuple
! tqsetc   ok set condition
! tqce     ok calculate equilibrium
! tqgetv   ok get equilibrium results as state variable values
! -------------------------
! tqgphc1  ok get phase constitution
! tqsphc1  ok set phase constitution
! tqcph1   ok calculate phase properties and return arrays
! tqcph2   ok calculate phase properties and return index
! tqdceq   ok delete equilibrium record
! tqcceq   ok copy current equilibrium to a new one
! tqselceq ok select new current equilibrium
! tqlr     ok list results 
! tqlc     ok list conditions
!
!------------------------------------------------------------
!
! The name of this library
module liboctq
!
! access to main OC library for equilibrium calculations and models
  use liboceqplus
!
  implicit none
!
  integer, parameter :: maxc=maxel,maxp=maxph
!
! This is for storage and use of components
  integer nel
  character, dimension(maxc) :: cnam*24
! Number of phase tuples
  integer ntup
! use the array PHASETUPLE available from OC
! save phase constitution to speed up calculation by interpolation
  double precision, allocatable, dimension(:,:) :: ysave
!
contains
!
!\begin{verbatim}
  subroutine tqini(n,ceq)
! initiate workspace
    implicit none
    integer n ! Not nused, could be used for some initial allocation
    type(gtp_equilibrium_data), pointer :: ceq ! EXIT: current equilibrium
!\end{verbatim}
! these should be provide linits and defaults
    integer intv(10)
    double precision dblv(10)
    intv(1)=-1
! This call initiates the OC package
!@CC
    if (allocated(eqlista)) then
       call new_gtp
    endif
    call init_gtp(intv,dblv)
!@CC
    ceq=>firsteq
    write(*,*)'tqini created: ',ceq%eqname
1000 continue
    return
  end subroutine tqini

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqrfil(filename,ceq)
! read all elements from a TDB file
    implicit none
    character*(*) filename  ! IN: database filename
    character ellista(10)*2  ! dummy
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim} %+
    integer iz
    character elname*2,name*24,refs*24
    double precision a1,a2,a3
! second argument 0 means ellista is ignored, all element read
    call readtdb(filename,0,ellista)
!    ceq=>firsteq
    nel=noel()
    do iz=1,nel
! store the element name in the cname array
       call get_element_data(iz,elname,name,refs,a1,a2,a3)
       cnam(iz)=elname
    enddo
! store phase tuples
    ntup=nooftup()
1000 continue
    return
  end subroutine tqrfil

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim} %-
  subroutine tqrpfil(filename,nsel,selel,ceq)
! read TDB file with selection of elements
    implicit none
    character*(*) filename  ! IN: database filename
    integer nsel
    character selel(*)*2  ! IN: elements to be read from the database
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    integer iz
    character elname*2,name*24,refs*24
    double precision a1,a2,a3
!
    call readtdb(filename,nsel,selel)
    if(gx%bmperr.ne.0) goto 1000
! is this really necessary??
!    ceq=>firsteq
    nel=noel()
    do iz=1,nel
! store element name in module array components
       call get_element_data(iz,elname,name,refs,a1,a2,a3)
       cnam(iz)=elname
    enddo
! store phase tuples and indices
    ntup=nooftup()
1000 continue
    return
  end subroutine tqrpfil
 
!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgcom(n,compnames,ceq)
! get system component names. At present the elements
    implicit none
    integer n                               ! EXIT: number of components
    character*24, dimension(*) :: compnames ! EXIT: names of components
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
    integer iz
    character elname*24,refs*24
    double precision a1,a2,a3
    do iz=1,nel
       compnames(iz)=' '
       call get_element_data(iz,compnames(iz),elname,refs,a1,a2,a3)
! store name in module array components also (already done when reading TDB)
       cnam(iz)=compnames(iz)
    enddo
    n=nel
1000 continue
    return
  end subroutine tqgcom

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgnp(n,ceq)
! get total number of phase tuples (phases and composition sets)
! A second composition set of a phase is normally placed after all other
! phases with one composition set
    implicit none
    integer n    !EXIT: n is number of phases
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
! NOTE the number composition sets may change at a calculation or if new
! composition sets are added or deleted explicitly
! This changes the number of phase tuples!
    ntup=nooftup()
    n=ntup
1000 continue
    return
  end subroutine tqgnp

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpn(phtupx,phasename,ceq)
! get name of phase tuple with index phtupx (ceq redundant)
    implicit none
    integer phtupx                  ! IN: index in phase tuple array
    character phasename*(*)      !EXIT: phase name, max 24+8 for pre/suffix
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    call get_phasetup_name(phtupx,phasename)
1000 continue
    return
  end subroutine tqgpn

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpi(phtupx,phasename,ceq)
! get index of phase phasename (including comp.set (ceq redundant)
    implicit none
    integer phtupx           !EXIT: phase tuple index
    character phasename*(*) !IN: phase name
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    call find_phasetuple_by_name(phasename,phtupx)
1000 continue
    return
  end subroutine tqgpi

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpcn2(n,c,constituentname,ceq)
! get name of consitutent with index c in phasetuple n
! NOTE An identical routine with different constituent index is tqgpcn
    implicit none
    integer n !IN: phase number (not phase tuple)
    integer c !IN: constituent index sequentially over all sublattices
    character constituentname*(24) !EXIT: costituent name
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    double precision mass
    call get_constituent_name(n,c,constituentname,mass)
!    write(*,*)'tqgpcn not implemented yet'
!    gx%bmperr=8888
1000 continue
    return
  end subroutine tqgpcn2

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpci(n,c,constituentname,ceq)
! get index of constituent with name in phase n
    implicit none
    integer n !IN: phase index
    integer c !IN: sequantial constituent index over all sublattices
    character constituentname*(*)
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
    write(*,*)'tqgpci not implemented yet'
    gx%bmperr=8888
1000 continue
    return
  end subroutine tqgpci

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpcs(n,c,stoi,mass,ceq)
! get stoichiometry of constituent c in phase n 
!? missing argument number of elements????
    implicit none
    integer n !IN: phase number
    integer c !IN: sequantial constituent index over all sublattices
    double precision stoi(*) !EXIT: stoichiometry of elements 
    double precision mass    !EXIT: total mass
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
    write(*,*)'tqgpcs not implemented yet'
    gx%bmperr=8888
1000 continue
    return
  end subroutine tqgpcs

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgccf(n1,n2,elnames,stoi,mass,ceq)
! get stoichiometry of component n1
! n2 is number of elements (dimension of elnames and stoi)
    implicit none
    integer n1 !IN: component number
    integer n2 !EXIT: number of elements in component
    character elnames(*)*(2) ! EXIT: element symbols
    double precision stoi(*) ! EXIT: element stoichiometry
    double precision mass    ! EXIT: component mass (sum of element mass)
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
    write(*,*)'tqgccf not implemented yet'
    gx%bmperr=8888
1000 continue
    return
  end subroutine tqgccf

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgnpc(n,c,ceq)
! get number of constituents of phase n
    implicit none
    integer n !IN: Phase number
    integer c !EXIT: number of constituents
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
    write(*,*)'tqgnpc not implemented yet'
    gx%bmperr=8888
1000 continue
    return
  end subroutine tqgnpc

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqcref(cix,phase,tpref,ceq)
! set component reference state
    integer cix
    character phase*(*)
    double precision tpref(*)
    type(gtp_equilibrium_data), pointer :: ceq  ! IN: current equilibrium
!\end{verbatim}
    integer phtupx
    call find_phasetuple_by_name(phase,phtupx)
    if(gx%bmperr.ne.0) goto 1000
    call set_reference_state(cix,phtupx,tpref,ceq)
1000 continue
    return
  end subroutine tqcref

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqphsts(phtupx,newstat,val,ceq)
! set status of phase tuple: SUSPEND, DORMANT, ENTERED, FIX
    integer phtupx,newstat
    double precision val
    type(gtp_equilibrium_data), pointer :: ceq  ! IN: current equilibrium
!\end{verbatim}
    integer n
    if(phtupx.le.0) then
! if tup<0 change status of all phases
       do n=1,ntup
          call change_phtup_status(n,newstat,val,ceq)
          if(gx%bmperr.ne.0) goto 1000
       enddo
    elseif(phtupx.le.ntup) then
       call change_phtup_status(phtupx,newstat,val,ceq)
    else
       write(*,*)'Illegal phase tuple index',phtupx
       gx%bmperr=8888
    endif
1000 continue
    return
  end subroutine tqphsts

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqsetc(stavar,n1,n2,value,cnum,ceq)
! set condition
! stavar is state variable as text
! n1 and n2 are auxilliary indices
! value is the value of the condition
! cnum is returned as an index of the condition.
! to remove a condition the value sould be equial to RNONE ????
! phase index is phase tuple index (include composition set)
! see TQGETV for doucumentation of stavar etc.
    implicit none
    integer n1             ! IN: 0 or phase tuple index or component number
    integer n2             ! IN: 0 or component number
    integer cnum           ! EXIT: sequential number of this condition
    character stavar*(*)   ! IN: character with state variable symbol
    double precision value ! IN: value of condition
    type(gtp_equilibrium_data), pointer :: ceq  ! IN: current equilibrium
!\end{verbatim}
    integer ip,ip2
    character cline*60,selvar*4,cval*24
!
!    write(*,11)'In tqsetc ',stavar(1:len_trim(stavar)),n1,n2,value
11  format(a,a,2i5,1pe14.6)
    cline=' '
! extract a value after an =
    ip=index(stavar,'=')
    if(ip.gt.0) then
       selvar=stavar(1:ip-1)
       cval=stavar(ip:)
!@CC
       ip2=index(stavar,'(')
       if(ip2.gt.0) then
          ip = ip2
          selvar=stavar(1:ip-1)
          cval=stavar(ip:)
       endif
!@CC
!       write(*,*)'Value after = :',cval
    else
       selvar=stavar
       cval=' '
    endif
    call capson(selvar)
    select case(selvar)
    case default
       write(*,*)'Condition wrong, not implemented or illegal: ',stavar
       gx%bmperr=8888; goto 1000
! Potentials T and P
    case('T   ','P   ')
       if(ip.gt.0) then
          cline=' '//stavar
       else
          write(cline,110)selvar(1:1),value
110       format(' ',a,'=',E15.8)
       endif
! Total amount or amount of a component in moles
    case('N   ')
       if(ip.gt.0) then
          cline=' '//stavar
       else
          if(n1.gt.0) then
!          call get_component_name(n1,name,ceq)
!          if(gx%bmperr.ne.0) goto 1000
             write(cline,112)selvar(1:1),cnam(n1)(1:len_trim(cnam(n1))),value
112       format(' ',a,'(',a,')=',E15.8)
!          write(*,*)'Setting condition: ',cline(1:len_trim(cline))
          else
             write(cline,110)selvar(1:1),value
          endif
       endif
! Overall fraction of a component 
    case('X   ','W   ')
! ?? fraction of phase component not implemented, n1 must be component number
!       call get_component_name(n1,cnam,ceq)
!       if(gx%bmperr.ne.0) goto 1000
       if(ip.gt.0) then
          cline=' '//stavar
       else
          write(cline,120)selvar(1:1),cnam(n1)(1:len_trim(cnam(n1))),value
120       format(1x,a,'(',a,')=',1pE15.8)
       endif
    case('H  ','V  ')
! enthalpy or volume of system
       if(cval(1:1).eq.'=') then
          cline=' '//stavar
       else
          write(cline,130)selvar(1:1),value
130       format(1x,a,'=',1pE15.8)
       endif
! case ....
! ?? MORE CONDITIONS WILL BE ADDED ...
    end select
!    write(*,*)'tqsetc condition: ',trim(cline)
    ip=1
    call set_condition(cline,ip,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error setting condition: ',cline(1:len_trim(cline)),ip
    endif
1000 continue
    return
  end subroutine tqsetc

!@CC
!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine toggle_dense_grid()
    if(btest(globaldata%status,GSXGRID)) then
       globaldata%status=ibclr(globaldata%status,GSXGRID)
       write(*,3110)'reset'
3110   format('Dense grid ',a)
    else
       globaldata%status=ibset(globaldata%status,GSXGRID)
       write(*,3110)'dense grid set'
    endif
    return
  end subroutine toggle_dense_grid
!@CC

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqce(target,n1,n2,value,ceq)
! calculate quilibrium with possible target
! Target can be empty or a state variable with indices n1 and n2
! value is the calculated value of target
    implicit none
    integer n1,n2,mode
    character target*(*)
    double precision value
    logical confirm
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
    integer nyfas,j1,j2
! mode=1 means start values using global gridminimization
    if(n1.lt.0) then
! this means calculate without grid minimuzer
       mode=0
       confirm=.FALSE.
! calcqeq3 is silent, no listing of phase changes etc.
       call calceq3(mode,confirm,ceq)
    else
       mode=1
       call calceq2(mode,ceq)
       if(gx%bmperr.eq.4204) then
! if the error code is "too many iterations" try without grid minimizer
! it converges in many cases
!          write(*,2048)gx%bmperr
2048      format('Error ',i5,', cleaning up and trying harder')
          gx%bmperr=0
          call calceq2(0,ceq)
       endif
    endif
    if(gx%bmperr.ne.0) goto 1000
! there may be new composition sets, update ntup
!    write(*,*)'Number of phase tuples: ',ntup
    nyfas=nooftup()
!    write(*,*)'Number of phase tuples: ',ntup,nyfas
    if(nyfas.ne.ntup) then
!       write(*,*)'Number of phase tuples changed: ',nyfas,ntup
       ntup=nyfas
!       if(allocated(ysave)) deallocate(ysave)
!       allocate(ysave(nyfas,maxconst))
    endif
! copy the constitution to a local save array
!    if(.not.allocated(ysave)) then
!       allocate(ysave(nyfas,maxconst))
!    endif
    if(allocated(ysave)) deallocate(ysave)
    allocate(ysave(nyfas,maxconst))
! the intention of saving constitution is to make it possible to interpolate
! the calculation of G if the constitution is changed very little
   do j1=1,nyfas
       do j2=1,size(ceq%phase_varres(phasetuple(j1)%lokvares)%yfr)
          ysave(j1,j2)=ceq%phase_varres(phasetuple(j1)%lokvares)%yfr(j2)
       enddo
    enddo
1000 continue
    return
  end subroutine tqce

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgetv(stavar,n1,n2,n3,values,ceq)
! get equilibrium results using state variables
! stavar is the state variable IN CAPITAL LETTERS with indices n1 and n2 
! n1 can be a phase tuple index, n2 a component index
! n3 at the call is the dimension of the array values, 
! changed to number of values on exit
! value is an array with the calculated value(s), n3 set to number of values.
    implicit none
    integer n1,n2,n3
    character stavar*(*)
    double precision values(*)
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!========================================================
! stavar must be a symbol listed below
! IMPORTANT: some terms explained after the table
! Symbol  index1,index2                     Meaning (unit)
!.... potentials
! T     0,0                                             Temperature (K)
! P     0,0                                             Pressure (Pa)
! MU    component,0 or ext.phase.index*1,constituent*2  Chemical potential (J)
! AC    component,0 or ext.phase.index,constituent      Activity = EXP(MU/RT)
! LNAC  component,0 or ext.phase.index,constituent      LN(activity) = MU/RT
!...... extensive variables
! U     0,0 or ext.phase.index,0   Internal energy (J) whole system or phase
! UM    0,0 or ext.phase.index,0       same per mole components
! UW    0,0 or ext.phase.index,0       same per kg
! UV    0,0 or ext.phase.index,0       same per m3
! UF    ext.phase.index,0              same per formula unit of phase
! S*3   0,0 or ext.phase.index,0   Entropy (J/K) 
! V     0,0 or ext.phase.index,0   Volume (m3)
! H     0,0 or ext.phase.index,0   Enthalpy (J)
! A     0,0 or ext.phase.index,0   Helmholtz energy (J)
! G     0,0 or ext.phase.index,0   Gibbs energy (J)
! ..... some extra state variables
! NP    ext.phase.index,0          Moles of phase
! BP    ext.phase.index,0          Mass of moles (kg)
! Q     ext.phase.index,0          Internal stability/RT (dimensionless)
! DG    ext.phase.index,0          Driving force/RT (dimensionless)
!....... amounts of components
! N     0,0 or component,0 or ext.phase.index,component    Moles of component
! X     component,0 or ext.phase.index,component   Mole fraction of component
! B     0,0 or component,0 or ext.phase.index,component     Mass of component
! W     component,0 or ext.phase.index,component   Mass fraction of component
! Y     ext.phase.index,constituent*1                    Constituent fraction
!........ some parameter identifiers
! TC    ext.phase.index,0                Magnetic ordering temperature
! BMAG  ext.phase.index,0                Aver. Bohr magneton number
! MQ&   ext.phase.index,constituent    Mobility
! THET  ext.phase.index,0                Debye temperature
! LNX   ext.phase.index,0                Lattice parameter
! EC11  ext.phase.index,0                Elastic constant C11
! EC12  ext.phase.index,0                Elastic constant C12
! EC44  ext.phase.index,0                Elastic constant C44
!........ NOTES:
! *1 The phase index is the phase tuple index (extra composition sets at end)
! *2 The constituent index is 10*species_number + sublattice_number
! *3 S, V, H, A, G, NP, BP, N, B and DG can have suffixes M, W, V, F also
!--------------------------------------------------------------------
! special addition for TQ interface: d2G/dyidyj
! D2G + phase tuple
!--------------------------------------------------------------------
!\end{verbatim}
    integer ics,mjj,nph,ki,kj,lp,lokph,lokcs
    character statevar*60,encoded*2048,name*24,selvar*4,norm*4
! mjj should be the dimension of the array values ...
    mjj=n3
    selvar=stavar
    call capson(selvar)
! for state variables like MQ&FE remove the part from & before the select
!    write(*,11)'In tqgetv: ',selvar,n1,n2,n3
11  format(a,a,3i5)
    norm=' '
    lp=index(selvar,'&')
    if(lp.gt.0) then
       selvar(lp:)=' '
    else
! check if variable is normallized, only M (per mole) allowed
       ki=len_trim(selvar)
       if(ki.ge.2) then
          if(selvar(ki:ki).eq.'M') then
             norm='M'
             selvar(ki:)=' '
             ki=ki-1
          endif
       endif
    endif
!=======================================================================
    kj=index(selvar,'(')
    if(kj.gt.0) then
       selvar=selvar(1:kj-1)
    endif
!    write(*,*)'tqgetv 0: ',kj,selvar,'>',stavar,'<'
    select case(selvar)
    case default
       write(*,*)'Unknown state variable: ',stavar(1:20),'>:<',selvar
       gx%bmperr=8888; goto 1000
!--------------------------------------------------------------------
! T or P
    case('T  ','P  ')
       call get_state_var_value(selvar,values(1),encoded,ceq)
!--------------------------------------------------------------------
! chemical potential for a component
    case('MU  ','MUS ')
       if(n1.lt.-1 .or. n1.eq.0) then
          write(*,*)'tqgetv 17: component number must be positive'
          gx%bmperr=8888; goto 1000
       elseif(n1 .eq.-1) then
! this means all components
          statevar=trim(selvar)//'(*)'
          call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
       elseif(n1.le.noel()) then
          statevar=trim(selvar)//'('//trim(cnam(n1))//') '
!       write(*,*)'tqgetv 4: ',statevar(1:len_trim(statevar))
! we must use index value(1) as the subroutine expect a single variable
          call get_state_var_value(statevar,values(1),encoded,ceq)
       else
          write(*,*)'No such component'
       endif
!--------------------------------------------------------------------
!@CC
! Amount of moles /mass of components in a phase
    case('NP  ', 'BP  ')
       if(n1.lt.0) then
! all phases
          statevar=stavar(1:2)//'(*)'
!@CC
! this returns all composition sets for all phases
          call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
! this output gives the amounts for all compsets of a phase sequentially
! but here we want them in phase tuple order
! the second argument is the number of values for each phase, here is 1 but
! it can be for example compositions, then it should be number of components
          call sortinphtup(n3,1,values)
       else
! NP for just one phase
          call get_phasetup_name(n1,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar='NP('//trim(name)//') '
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       endif
!--------------------------------------------------------------------
! Mole or mass fractions
!@CC
    case('N   ','B    ','X   ','W   ')
!@CC
!       write(*,*)'in tqgetv n,x,w: ',n1,n2,n3
       if(n2.eq.0) then
          if(n1.lt.0) then
! moles, mole or mass fraction of all components for all phases
             statevar=stavar(1:1)//'(*) '
!             write(*,*)'tqgetv 3: ',mjj,statevar(1:len_trim(statevar))
             call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
          elseif(n1.eq.0) then
! mole fraction for the state variable written as X(FE)
! n1 and n2 not used, just check for wildcard
!             write(*,*)'tqgetv 20: ',stavar(1:len_trim(stavar))
             if(index(stavar,'*').gt.0) then
                call get_many_svar(stavar,values,mjj,n3,encoded,ceq)
             else
                call get_state_var_value(stavar,values(1),encoded,ceq)
             endif
          else
! mole fraction of a single component, no phase specification
             n3=1
             ics=1
!             call get_component_name(n1,name,ceq)
!             if(gx%bmperr.ne.0) goto 1000
             statevar=stavar(1:1)//'('//trim(cnam(n1))//')'
!             write(*,*)'tqgetv 4: ',statevar(1:len_trim(statevar))
             call get_state_var_value(statevar,values(1),encoded,ceq)
          endif
       elseif(n1.lt.0) then
!........................................................
! for all phases one or several components
          if(n2.lt.0) then
! this means all components all phases, for example x(*,*)
             statevar=stavar(1:1)//'(*,*) '
!             write(*,*)'tqgetv 5: ',mjj,statevar(1:len_trim(statevar))
             call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
! this output gives the composition for all compsets of a phase sequentially
! but we want them in phase tuple order
! The second argument is the number of values for each phase, noel()
! in this case
             ics=noel()
             call sortinphtup(n3,ics,values)
          else
! a single component in all phases. n2 must not be zero
!             call get_component_name(n2,name,ceq)
!             if(gx%bmperr.ne.0) goto 1000
             if(n2.le.0 .or. n2.ge.noel()) then
                write(*,*)'No such component'
                goto 1000
             endif
! state variable like w(*,cr), the Cr content in all (stable) phases
             statevar=stavar(1:1)//'(*,'//cnam(n2)(1:len_trim(cnam(n2)))//')'
!             write(*,*)'tqgetv 6: ',mjj,statevar(1:len_trim(statevar))
             call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
! this output gives the composition for all compsets of a phase sequentially
! but we want them in phase tuple order
! The second argument is the number of values for each phase, in this case 1
!             ics=noel()
! THIS MUST BE CHECKED !!!
             call sortinphtup(n3,1,values)
          endif
       elseif(n2.lt.0) then
! this means all components in one phase
          call get_phasetup_name(n1,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=stavar(1:1)//'('//trim(name)//',*) '
          call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
       else
! one component (n2) of one phase (n1)
          call get_phasetup_name(n1,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=stavar(1:1)//'('//trim(name)//','
          call get_component_name(n2,name,ceq)
          if(gx%bmperr.ne.0) goto 1000
          statevar(len_trim(statevar)+1:)=trim(name)//') '
!          write(*,*)'tqgetv 8: ',statevar
          call get_state_var_value(statevar,values(1),encoded,ceq)
       endif
!--------------------------------------------------------------------
! volume
    case('V   ')
       if(norm(1:1).ne.' ') then
          statevar='V'//norm
          ki=2
       else
          statevar='V '
          ki=1
       endif
       if(n1.gt.0) then
! Volume for a specific phase
          call get_phasetup_name(n1,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=statevar(1:ki)//'('//trim(name)//') '
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       else
! Total volume
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       endif
!--------------------------------------------------------------------
! Enthalpy
    case('H   ')
! phase specifier not allowed
       if(norm(1:1).ne.' ') then
          statevar='H'//norm
          ki=2
       else
          statevar='H '
          ki=1
       endif
!       write(*,*)'tqgetv 1: ',n1,ki
       if(n1.gt.0) then
! Gibbs energy for a specific phase
          call get_phasetup_name(n1,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=statevar(1:ki)//'('//trim(name)//') '
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       else
! Total enthalpy
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       endif
!--------------------------------------------------------------------
! Gibbs energy
    case('G   ')
! phase specifier not allowed
       if(norm(1:1).ne.' ') then
          statevar='G'//norm
          ki=2
       else
          statevar='G '
          ki=1
       endif
!       write(*,*)'tqgetv 1: ',n1,ki
       if(n1.gt.0) then
! Gibbs energy for a specific phase
          call get_phasetup_name(n1,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=statevar(1:ki)//'('//trim(name)//') '
!          write(*,*)'tqgetv 3: ',statevar
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       else
! Total Gibbs energy 
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       endif
!--------------------------------------------------------------------
! Driving force relative stable equilibrium
    case('DG  ')
! Always normalized per mole
       if(norm(3:3).ne.' ') then
          statevar='DG'//norm
          ki=3
       else
          statevar='DG '
          ki=2
       endif
       write(*,*)'tqgetv DGM: ',n1,ki
       if(n1.gt.0) then
! The driving force for a specific phase
          call get_phasetup_name(n1,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=statevar(1:ki)//'M('//trim(name)//') '
!          write(*,*)'tqgetv 3: ',statevar
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       else
! For all phases
          n3=0
          if(nooftup().gt.mjj) then
             write(*,*)'TQGETV error, array too small for DGM',mjj,nooftup()
             gx%bmperr=8888
             goto 1000
          endif
          statevar='DGM(#) '
          write(*,*)'tqgetv 3: ',statevar
          call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
          write(*,'(a,10(1pe12.4))')'TQGETV: ',(values(ki),ki=1,n3)
          write(*,*)'gx%bmperr: ',gx%bmperr
       endif
!--------------------------------------------------------------------
! Mobilities
    case('MQ   ')
       call get_phasetup_name(n1,name)
       if(gx%bmperr.ne.0) goto 1000
       statevar=stavar(1:len_trim(stavar))//'('//trim(name)//')'
!       write(*,*)'statevar: ',statevar
       call get_state_var_value(statevar,values(1),encoded,ceq)
!--------------------------------------------------------------------
! Second derivatives of the Gibbs energy of a phase
    case('D2G   ')
       lokcs=phasetuple(n1)%lokvares
! this gives wrong value!! ??
       n3=size(ceq%phase_varres(lokcs)%yfr)
!       write(*,*)'D2G 3: ',n3
       kj=(n3*(n3+1))/2
       if(kj.gt.mjj) then
          write(*,*)'TQGETV error, array too small for D2G',mjj,kj
          gx%bmperr=8888
          goto 1000
       endif
!       write(*,*)'D2G 3: ',kj
       do ki=1,kj
          values(ki)=ceq%phase_varres(lokcs)%d2gval(ki,1)
       enddo
    end select
!===========================================================================
1000 continue
    return
  end subroutine tqgetv

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

  subroutine tqgetg(lokres,n1,n2,values,ceq)
! the partial derivative of the Gibbs energy ....??
    implicit none
    integer n1,n2,lokres
    double precision values(*)
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!    
    double precision napfu, rgast
    integer count
    integer jl,size
    TYPE(gtp_phase_varres), pointer :: parres
!
    count = 1
!
    napfu=ceq%phase_varres(lokres)%abnorm(1)
    rgast=globaldata%rgas*ceq%tpval(1)
    parres=>ceq%phase_varres(lokres)
!  
!    write(*,100)(rgast*parres%gval(jl,1),jl=1,4)
!    write(*,200)parres%gval(1,1)/parres%abnorm(1),parres%abnorm(1)
100 format('G/N, dG/dT:',4(1PE16.8))
200 format('G/N/RT, N:',2(1PE16.8))
!   G_m^\alpha = G_M^\alpha/N^\alpha, \frac{\partial G_m^\alpha}{\partial T},
! \frac{\partial G_m^\alpha}{\partial P},
! \frac{\partial^2 G_m^\alpha}{\partial T^2}
    values(count:count+3) = rgast*parres%gval(1:4,1)/napfu
    count = count + 4
    if (n1>0) then
!      1/N^\alpha * \frac{\partial G_M^\alpha}{\partial y_i}
       values(count:count+n1-1) = rgast*parres%dgval(1,1:n1,1)/napfu
       count = count + n1
       if (n2>0) then
!         1/N^\alpha * \frac{\partial^2 G_M^\alpha}{\partial y_i\partial y_j}
          values(count:count+n2-1) = rgast*parres%d2gval(1:n2,1)/napfu
       endif
    endif
  end subroutine tqgetg

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\
  
  subroutine tqgdmat(phtupx,tpval,xknown,cpot,tyst,nend,mugrad,mobval,&
       consnames,n1,ceq)
! equilibrates the constituent fractions of a phase for mole fractions xknown
! and calculates the Darken matrix and unreduced diffusivities
! phtup is phase tuple
! tpval is T and P
! ceq is a datastructure with all relevant thermodynamic data
! cpot are the (calculated) chemical potentials
! tyst is TRUE means no outut
! nend is the number of values returned in mugrad
! mugrad are the derivatives of the chemical potentials wrt mole fractions??
! mobval are the mobilities
    implicit none
    integer phtupx                  ! IN: index in phase tuple array
    integer nend
    logical tyst
    double precision tpval(*),xknown(*),cpot(*),mugrad(*),mobval(*)
    character*24, dimension(*) :: consnames 
    integer n1
    TYPE(gtp_phasetuple), pointer :: phtup
    TYPE(gtp_equilibrium_data), pointer :: ceq

    integer iph, ics, ll
    double precision mass
    character*24 spname
             
    phtup=>phasetuple(phtupx)    
    call equilph1d(phtup,tpval,xknown,cpot,tyst,nend,mugrad,mobval,ceq)
    
    iph=phasetuple(phtupx)%ixphase
    ics=1   
    n1 = noconst(iph,ics,firsteq)
    do ll=1,n1
       call get_constituent_name(iph,ll,consnames(ll),mass)
    enddo

  end subroutine tqgdmat
!@CC

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgphc1(n1,nsub,cinsub,spix,yfrac,sites,extra,ceq)
! tq_get_phase_constitution
! This subroutine returns the sublattices and constitution of a phase
! n1 is phase tuple index
! nsub is the number of sublattices (1 if no sublattices)
! cinsub is an array with the number of constítuents in each sublattice
! spix is an array with the species index of the constituents in all sublattices
! sites is an array of the site ratios for all sublattices.  
! yfrac is the constituent fractions in same order as in spix
! extra is an array with some extra values: 
!    extra(1) is the number of moles of components per formula unit
!    extra(2) is the net charge of the phase
    implicit none
    integer n1,nsub,cinsub(*),spix(*)
    double precision sites(*),yfrac(*),extra(*)
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    call get_phase_data(phasetuple(n1)%ixphase,phasetuple(n1)%compset,&
         nsub,cinsub,spix,yfrac,sites,extra,ceq)
1000 continue
    return
  end subroutine tqgphc1
  
!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqsphc1(n1,yfra,extra,ceq)
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
    integer n1
    double precision yfra(*),extra(*)
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    call set_constitution(phasetuple(n1)%ixphase,phasetuple(n1)%compset,&
         yfra,extra,ceq)
1000 continue
    return
  end subroutine tqsphc1

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqcph1(n1,n2,n3,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,ceq)
! tq_calculate_phase_properties
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! WARNING: this is not a subroutine to calculate chemical potentials
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
    integer n1,n2,n3
    double precision gtp(6),dgdy(*),d2gdydt(*),d2gdydp(*),d2gdy2(*)
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
    integer ij,lokres,nofc
!    write(*,*)'tqcph1 1: ',ceq%eqname
!    write(*,*)'tqcph1 2',phasetuple(n1)%ixphase,phasetuple(n1)%compset
!----------------------------------------------------------------------
! THIS IS NO EQUILIBRIUM, JUST G AND DERIVATIVES FOR CURRENT T, P AND Y
    call calcg(phasetuple(n1)%ixphase,phasetuple(n1)%compset,n2,lokres,ceq)
!----------------------------------------------------------------------
!    write(*,*)'tqcph1 3A',lokres,gx%bmperr
! The inital size here can be 1000
!    n3=size(ceq%phase_varres(lokres)%yfr)
! the actual number of constituents is better to take from this call
    n3=noconst(phasetuple(n1)%ixphase,1,ceq)
!    write(*,*)'tqcph1 3C',n3
! gval last index is the property, other properties can also be extracted
! t.ex. mobilites 
! The application program can also access these data directly ...
    if(gx%bmperr.eq.0) then
       do ij=1,6
          gtp(ij)=ceq%phase_varres(lokres)%gval(ij,1)
       enddo
       do ij=1,n3
          dgdy(ij)=ceq%phase_varres(lokres)%dgval(1,ij,1)
          d2gdydt(ij)=ceq%phase_varres(lokres)%dgval(2,ij,1)
          d2gdydp(ij)=ceq%phase_varres(lokres)%dgval(3,ij,1)
       enddo
! size of upper triangle of symetrix matrix
       nofc=n3*(n3+1)/2
       do ij=1,nofc
          d2gdy2(ij)=ceq%phase_varres(lokres)%d2gval(ij,1)
       enddo
    else
       gtp=zero
       do ij=1,nofc
          dgdy(ij)=zero
          d2gdydt(ij)=zero
          d2gdydp(ij)=zero
       enddo
       nofc=nofc*(nofc+1)/2
       do ij=1,nofc
          d2gdy2(ij)=zero
       enddo
    endif
1000 continue
    return
  end subroutine tqcph1

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim} %-
  subroutine tqcph2(n1,n2,n3,n4,ceq)
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
! n2 is type of calculation (0, 1 or 2)
! n3 is returned as number of constituents
! n4 is index to ceq%phase_varres(lokres)% with all results
! for indexing one can use the integer function ixsym(i1,i2)
    implicit none
    integer n1,n2,n3,n4
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer ij,lokres,nofc
!    write(*,*)'tqcph1 1: ',ceq%eqname
!    write(*,*)'tqcph1 2',phasetuple(n1)%ixphase,phasetuple(n1)%compset
!----------------------------------------------------------------------
! THIS IS NO EQUILIBRIUM, JUST G AND DERIVATIVES FOR CURRENT T, P AND Y
    call calcg(phasetuple(n1)%ixphase,phasetuple(n1)%compset,n2,lokres,ceq)
!----------------------------------------------------------------------
!    write(*,*)'tqcph1 3A',lokres,gx%bmperr
! this should work but gave segmentation fault, find this a more cumbersum way
    n3=size(ceq%phase_varres(lokres)%yfr)
    n4=lokres
! Uer can access results like
! ceq%phase_varres(n4)%gval(1..6,1..prop)
! prop=1 is G, other can be t.ex. Curie T, mobilites etc
! ceq%phase_varres(lokres)%dgval(1,ij,1) are dG/dy(ij)
! ceq%phase_varres(lokres)%dgval(2,ij,1) are d2G/dy(ij)dT
! ceq%phase_varres(lokres)%dgval(3,ij,1) are d2G/dy(ij)dP
! ceq%phase_varres(lokres)%d2gval(ij,1) are d2G/dy(i)dy(j)
! arranged as a single dimenion array indexed by ixsym(i,j)
!
! NEVER CHANGE THE CONSTITUTION DIRECTLY, using n4, ALWAYS CALL tqsph1(...)
!
1000 continue
    return
  end subroutine tqcph2

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqdceq(name)
! delete equilibrium with name
    implicit none
    character name*24
!    integer n1
    type(gtp_equilibrium_data), pointer :: newceq,ceq
!\end{verbatim}
    integer n1
    call findeq(name,n1)
    if(gx%bmperr.ne.0) goto 1000
! do not allow delete equilibrium 1
    if(n1.eq.1) then
       write(*,*)'No allowed to delete default equilibrium'
       gx%bmperr=4333
       goto 1000
    endif
!    ceq=>eqlista(n1)
    call delete_equilibria(name,ceq)
1000 continue
    return
  end subroutine tqdceq

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqcceq(name,n1,newceq,ceq)
! copy_current_equilibrium to newceq
! creates a new equilibrium record with name with values same as ceq
! n1 is returned as index
    implicit none
    character name*24
    integer n1
    type(gtp_equilibrium_data), pointer :: newceq,ceq
!\end{verbatim}
    !call enter_equilibrium(name,n1)
    !if(gx%bmperr.ne.0) goto 1000
    !newceq=>eqlista(n1)
    call copy_equilibrium(newceq,name,ceq)
1000 continue
    return
  end subroutine tqcceq

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqcneq(name,n1,newceq)
! creates a new equilibrium record, same but simpler call than tqcceq
! n1 is returned as index in eqlista
    implicit none
    character*(*), intent(in) :: name
    integer, intent(out) :: n1
    type(gtp_equilibrium_data), pointer, intent(out) :: newceq
!\end{verbatim}
    call enter_equilibrium(name,n1)
    if(gx%bmperr.ne.0) goto 1000
    newceq=>eqlista(n1)
1000 continue
    return
  end subroutine tqcneq

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqselceq(name,ceq)
! select current equilibrium to be that with name.
! Note that equilibria can be deleted and change number but not name
    implicit none
    character name
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer n1
    call findeq(name,n1)
    if(gx%bmperr.ne.0) goto 1000
    call selecteq(n1,ceq)
1000 continue
    return
  end subroutine tqselceq

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqlr(lut,ceq)
! list the equilibrium results like in OC
    implicit none
    integer lut
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer phtupx,iph,ics,lokvares,mode
    logical once
    write(lut,10)
10  format(/20('*')/'Start debug output from TQLR: ')
    call list_conditions(lut,ceq)
    call list_global_results(lut,ceq)
    call list_components_result(lut,1,ceq)
    once=.TRUE.
    mode=0
    do phtupx=1,nooftup()
       lokvares=phasetuple(phtupx)%lokvares
       if(ceq%phase_varres(lokvares)%phstate.ge.phentstab) then
          iph=phasetuple(phtupx)%ixphase
          ics=phasetuple(phtupx)%compset
          call list_phase_results(iph,ics,mode,lut,once,ceq)
       endif
    enddo
    write(lut,20)
20  format('End debug output from TQLR'/20('*')/)
1000 continue
    return
  end subroutine tqlr

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqlc(lut,ceq)
! list conditions like in OC
    implicit none
    integer lut
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    write(lut,10)
10  format(/'Debug output from TQLC: ')
    call list_conditions(lut,ceq)
1000 continue
    return
  end subroutine tqlc

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqquiet(yes)
! if argument TRUE spurious output should be suppressed
    implicit none
    logical yes
!\end{verbatim}
    if(yes) then
       globaldata%status=ibclr(globaldata%status,GSVERBOSE)
       globaldata%status=ibset(globaldata%status,GSSILENT)
    else
       globaldata%status=ibset(globaldata%status,GSVERBOSE)
    endif
    return
  end subroutine tqquiet

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqchange_globalbit(bit,onoff)
! set a global bit
    implicit none
    integer bit,onoff
!\end{verbatim}
! list here taken from models/gtp3.F90, only some allowed!!
! BEWHEARE, the meaning of bits may have changed !!! check with gtp3.F90
!  4 NOMERGE: no merge of gridmin result, 
!  5 NODATA: not any data, 
!  6 NOPHASE: no phase in system, 
!  7 NOACS: no automatic creation of composition set for any phase
!  8 NOREMCS: do not remove any redundant unstable composition sets
!  9 NOSAVE: data changed after last save command
! 10 VERBOSE: maximum of listing
! 11 SETVERB: permanent setting of verbose
! 12 SILENT: as little output as possible
! 13 NOAFTEREQ: no manipulations of results after equilibrium calculation
! 14 XGRID: extra dense grid for all phases
! 15 NOPAR: do not run in parallel
! 16 NOSMGLOB do not test global equilibrium at node points
! 17 NOTELCOMP the elements are not the components
! 18 TGRID use grid minimizer to test if global after calculating equilibrium
! 19 OGRID use old grid generator
! 20 NORECALC do not recalculate equilibria even if global test after fails
! 21 OLDMAP use old map algorithm
! 22 NOAUTOSP do not generate automatic start points for mapping
! 23 GSYGRID extra dense grid
! 24 GSVIRTUAL (CCI) enables calculations with a virtual element
    if((bit.ge.7 .and. bit.le.16) .or. (bit.ge.18 .and. bit.le.23)) then
       if(onoff.gt.0) then
! set bit
          globaldata%status=ibset(globaldata%status,bit)
       else
          globaldata%status=ibclr(globaldata%status,bit)
       endif
    else
       gx%bmperr=4326
    endif
    return
  end subroutine tqchange_globalbit

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqchange_phasebit(phtupx,bit,onoff)
! set a bit of phase
    implicit none
    integer phtupx,bit,onoff
!\end{verbatim}
! taken from models/gtp3.F90
!-Bits in PHASE record STATUS1 there are also bits in each phase_varres record!
! BEWHEARE, the meaning of bits may have changed !!! check with gtp3.F90
!  0 HID phase is hidden (not implemented)
!  1 IMHID phase is implictly hidden (not implemented)
!  2 ID phase is ideal, substitutional and no interaction
!  3 NOCV phase has no concentration variation (fix composition)
!  4 HASP phase has at least one parameter entered
!  5 FORD phase has 4 sublattice FCC ordering with parameter permutations
!  6 BORD phase has 4 sublattice BCC ordering with parameter permutations
!  7 SORD phase has TCP type ordering (like for sigma)
!  8 MFS phase has a disordered fraction set
!  9 GAS this is the gas phase (first in phase list) 
! 10 LIQ phase is liquid (can be several but listed directly after gas)
! 11 IONLIQ phase has ionic liquid model (I2SL)
! 12 AQ1 phase has aqueous model (not implemented)
! 13 STATE elemental liquid twostate (2-state) model parameter UNUSED?
! 14 QCE phase has quasichemical SRO configurational entropy (not implemented)
! 15 CVMCE phase has some CVM ordering entropy (not implemented)
! 16 EXCB phase need explicit charge balance (has ions)
! 17 XGRID use extra dense grid for this phase
! 18 FACTCE phase has FACT quasichemical SRO model (not implemented)
! 19 NOCS not allowed to create composition sets for this phase
! 20 HELM parameters are for a Helmholz energy model (not implemented),
! 21 PHNODGDY2 phase has model with no analytical 2nd derivatives
! 22 not implemented ELMA phase has elastic model A (not implemented)
! 23 EECLIQ the condensed phase (liquid) that should have highest entropy
! 24 PHSUBO special use testing models DO NOT USE
! 25 PALM interaction records numbered by PALMTREE NEEDED FOR PERMUTATIONS !!!
! 26 MULTI may be used with care
! 27 BMAV Xion magnetic model with average Bohr magneton number
! 28 UNIQUAC The UNIQUAC fluid model
! 29 DILCE phase has dilute configigurational entropy (not implemented)
! only bittar 3 left!
    integer lokph
    if(phtupx.le.0 .or. phtupx.gt.nooftup()) then
       gx%bmperr=4325
    elseif(bit.eq.17 .or. bit.eq.19) then
       lokph=phasetuple(phtupx)%lokph
       if(onoff.gt.0) then
          call set_phase_status_bit(lokph,bit)
       else
          call clear_phase_status_bit(lokph,bit)
       endif
    else
       gx%bmperr=4326
    endif
    return
  end subroutine tqchange_phasebit

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqset_gaddition(phtupx,gadd,ceq)
! set fix addition to Gibbs energy of a phase#compset
    implicit none
    integer phtupx
    double precision gadd
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! Provided by Christophe Sigli 2018?
    integer lokcs
    lokcs=phasetuple(phtupx)%lokvares
    if(.not.allocated(ceq%phase_varres(lokcs)%addg)) then
       allocate(ceq%phase_varres(lokcs)%addg(1))
    endif
    ceq%phase_varres(lokcs)%addg(1)=gadd
! set bit that this should be calculated
    ceq%phase_varres(lokcs)%status2=&
         ibset(ceq%phase_varres(lokcs)%status2,CSADDG)
    return
  end subroutine tqset_gaddition

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tq_add_const_energy(energy,phtupx,ceq)
! add a constant energy in J/mole
    double precision,intent(in) :: energy
    type(gtp_equilibrium_data), pointer :: ceq
    integer,intent(in) :: phtupx
!\end{verbatim}
! Provided by Jan Herrnring 2020.12.15
    integer :: lokcs
    lokcs=phasetuple(phtupx)%lokvares
    if(.not.allocated(ceq%phase_varres(lokcs)%addg)) then
       allocate(ceq%phase_varres(lokcs)%addg(1))
    endif
! add a constant term to G, value in J/FU
! Abnorm is the number of moles of the phase
    ceq%phase_varres(lokcs)%addg(1)=energy*ceq%phase_varres(lokcs)%abnorm(1)
! set bit that this should be calculated
    ceq%phase_varres(lokcs)%status2=&
         ibset(ceq%phase_varres(lokcs)%status2,CSADDG)
  end subroutine tq_add_const_energy

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

  
end MODULE LIBOCTQ

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\
! dummy modules
!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

module ftinyopen
  !
  ! This module replaces a C module for a popup window to open files
  ! used in the interactive OC.  If you want to use the original
  ! version for opening files please check the linkmake or Makefile
  !
contains

  subroutine getfilename(typ,sval)
    implicit none
    integer typ
    character sval*(*)
    sval=' '
    return
  end subroutine getfilename

end module ftinyopen

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\
! dummy module (only Linux)
!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

module M_getkey
  !
  ! This module replaces a C module fore single character input on Linux
  !
contains

  character function getkex()
    getkex=' '
    return
  end function getkex

end module M_getkey

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

