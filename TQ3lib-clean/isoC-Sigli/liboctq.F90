!
! Minimal TQ interface.
!
! To compile and link this with an application one must first compile
! and form a library with of the most OC subroutines (oclib.a)
!  and to copy this and the corresponding "mov" files from this compilation 
! to the folder with this library
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
! When not using Fortran 95 (or later) one can probably replace this
! with a 2-dimensional array with first index phase number and second
! the comp.set number.
!
! For constituents an EXTENDED CONSTITUENT INDEX is sometimes used 
! and equal to 10*species_number + sublattice
!
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
! tqgcom   ok get system component names
! tqgnp    ok get number of phase tuples (phases and comp. sets
! tqgpn    ok get name of phase tuple
! tqgpi    ok get phase tuple index of phase using its name
! tqgpcn   -  get name of constituent of a phase using index
! tqgpci   -  get index of constituent of a phase using name
! tqgpcs   -  get stoichiometry of species as system components 
! tqgccf   -  get stoichiometry of system component as elements
! tqgnpc   -  get number of constituents in phase
! -------------------------
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
  integer, parameter :: maxc=40,maxp=500
!
! This is for storage and use of components
  integer nel
  character, dimension(maxc) :: cnam*24
  double precision cmass(maxc)
! This is for storage and use of phase+composition tuples
  integer ntup
  type(gtp_phasetuple), dimension(maxp) :: phcs
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
    call init_gtp(intv,dblv)
    ceq=>firsteq
	if(allocated(firstash%eqlista)) deallocate(firstash%eqlista)

!    write(*,*)'tqini created: ',ceq%eqname
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
!\end{verbatim}
    integer iz
    character elname*2,name*24,refs*24
    double precision mass,h298,s298
! second argument 0 means ellista is ignored, all element read
    call readtdb(filename,0,ellista)
!    ceq=>firsteq
    nel=noel()
    do iz=1,nel
! store element name in module array components
       call get_element_data(iz,elname,name,refs,mass,h298,s298)
       cnam(iz)=elname
	   cmass(iz)=mass
	   write(*,*) iz,mass
    enddo
! store phase tuples and indices
    ntup=get_phtuplearray(phcs)
1000 continue
    return
  end subroutine tqrfil

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
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
     double precision mass,h298,s298
!
    call readtdb(filename,nsel,selel)
    if(gx%bmperr.ne.0) goto 1000
! is this really necessary??
!    ceq=>firsteq
    nel=noel()
    do iz=1,nel
! store element name in module array components
       call get_element_data(iz,elname,name,refs,mass,h298,s298)
       cnam(iz)=elname
	   cmass(iz)=mass
    enddo
! store phase tuples and indices
    ntup=get_phtuplearray(phcs)
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
1000 continue
    return
  end subroutine tqgcom

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgnp(n,ceq)
! get total number of phases and composition sets
! A second composition set of a phase is normally placed after all other
! phases with one composition set
    implicit none
    integer n    !EXIT: n is number of phases
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
! This call fills the module array phcs with phase and composition set indices
! NOTE the number composition sets may change at a calculation or if new
! composition sets are added or deleted explicitly
    ntup=get_phtuplearray(phcs)
    n=ntup
1000 continue
    return
  end subroutine tqgnp

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpn(phcsx,phasename,ceq)
! get name of phase+compset tuple with index phcsx 
    implicit none
    integer phcsx                  ! IN: index in phase tuple array
!    TYPE(gtp_phasetuple), pointer :: phcs  !IN: phase number and comp.set
    character phasename*(*)      !EXIT: phase name, max 24+8 for pre/suffix
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    call get_phase_name(phcs(phcsx)%phaseix,phcs(phcsx)%compset,phasename)
1000 continue
    return
  end subroutine tqgpn

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpi(phcsx,phasename,ceq)
! get index of phase phasename (including comp.set, ceq not needed ...
    implicit none
    integer phcsx           !EXIT: phase tuple index
    character phasename*(*) !IN: phase name
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    call find_phasetuple_by_name(phasename,phcsx)
1000 continue
    return
  end subroutine tqgpi

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpcn(n,c,constituentname,ceq)
! get name of consitutent c in phase n
! NOTE An identical routine with different constituent index is tqgpcn2
    implicit none
    integer n !IN: phase number
    integer c !IN: extended constituent index: 10*species_number+sublattice
    character constituentname*(24) !EXIT: costituent name
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    write(*,*)'tqgpcn not implemented yet'
    gx%bmperr=8888
1000 continue
    return
  end subroutine tqgpcn

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpcn2(n,c,constituentname,ceq)
! get name of consitutent c in phase n
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
    integer c !EXIT: extended constituent index: 10*species_number+sublattice
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
    integer c !IN: extended constituent index: 10*species_number+sublattice
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
  subroutine tqphsts(tup,newstat,val,ceq)
! set status of phase tuple, 
    integer tup,newstat
    double precision val
    type(gtp_equilibrium_data), pointer :: ceq  ! IN: current equilibrium
!\end{verbatim}
    integer n
    if(tup.le.0) then
       do n=1,ntup
          call change_phase_status(phcs(n)%phaseix,phcs(n)%compset,&
               newstat,val,ceq)
          if(gx%bmperr.ne.0) goto 1000
       enddo
    elseif(tup.le.ntup) then
       call change_phase_status(phcs(tup)%phaseix,phcs(tup)%compset,&
            newstat,val,ceq)
    else
       write(*,*)'Illegal phase tuple index'
       gx%bmperr=5001
    endif
1000 continue
    return
  end subroutine tqphsts

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\
  
!\begin{verbatim}
  subroutine tqsetc(stavar,n1,n2,myvalue,cnum,ceq)
! set condition
! stavar is state variable as text
! n1 and n2 are auxilliary indices
! value is the value of the condition
! cnum is returned as an index of the condition.
! to remove a condition the value sould be equial to RNONE ????
! when a phase indesx is needed it should be 10*nph + ics
! SEE TQGETV for doucumentation of stavar etc.
    implicit none
    integer n1             ! IN: 0 or phase tuple index or component number
    integer n2             ! IN: 0 or component number
    integer cnum           ! EXIT: sequential number of this condition
    character stavar*(*)   ! IN: character with state variable symbol
    double precision myvalue ! IN: value of condition
    type(gtp_equilibrium_data), pointer :: ceq  ! IN: current equilibrium
!\end{verbatim}
    integer ip
    character cline*60,selvar*4
!
!    write(*,11)'In tqsetc ',stavar(1:len_trim(stavar)),n1,n2,value
11  format(a,a,2i5,1pe14.6)
    cline=' '
    selvar=stavar
    call capson(selvar)
    select case(selvar)
    case default
       write(*,*)'Condition wrong, not implemented or illegal: ',stavar
       gx%bmperr=8888; goto 1000
! Potentials T and P
    case('T   ','P   ')
       write(cline,110)selvar(1:1),myvalue
110    format(' ',a,'=',E15.8)
! Total amount or amount of a component in moles
    case('N   ')
       if(n1.gt.0) then
!          call get_component_name(n1,name,ceq)
!          if(gx%bmperr.ne.0) goto 1000
          write(cline,112)selvar(1:1),cnam(n1)(1:len_trim(cnam(n1))),myvalue
112       format(' ',a,'(',a,')=',E15.8)
       else
          write(cline,110)selvar(1:1),myvalue
       endif
! Overall fraction of a component 
    case('X   ','W   ')
! ?? fraction of phase component not implemented, n1 must be component number
!       call get_component_name(n1,cnam,ceq)
!       if(gx%bmperr.ne.0) goto 1000
       write(cline,120)selvar(1:1),cnam(n1)(1:len_trim(cnam(n1))),myvalue
120    format(1x,a,'(',a,')=',1pE15.8)
! ?? MORE CONDITIONS WILL BE ADDED ...
    end select
!    write(*,*)'tqsetc condition: ',cline(1:len_trim(cline))
    ip=1
    call set_condition(cline,ip,ceq)
1000 continue
    return
  end subroutine tqsetc

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqce(target,n1,n2,myvalue,ceq)
! calculate quilibrium with possible target
! Target can be empty or a state variable with indices n1 and n2
! value is the calculated value of target
    implicit none
    integer n1,n2,mode
    character target*(*)
    double precision myvalue
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
! mode=1 means start values using global gridminimization
    mode=1
    if(n1.lt.0) then
! this means calculate without grid minimuzer
      ! write(*,*)'No grid minimizer'
       mode=0
    endif
	!call calceq2(mode,ceq)
     call calceq3(mode,.FALSE.,ceq) ! the same as calceq2 but silent
    if(gx%bmperr.ne.0) goto 1000
! there may be new composition sets, update tup and phcs
! this call updates both the number of tuples and the phcs array
    ntup=get_phtuplearray(phcs)
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
! myvalue is an array with the calculated myvalue(s), n3 set to number of values.
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
! *1 The ext.phase.index is   10*phase_number+comp.set_number
! *2 The constituent index is 10*species_number + sublattice_number
! *3 S, V, H, A, G, NP, BP, N, B and DG can have suffixes M, W, V, F also
!--------------------------------------------------------------------
! special addition for TQ interface: d2G/dyidyj
! D2G + extended phase index
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
! check if variable is normallized
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
       write(*,*)'Unknown state variable: ',stavar(1:20),'>:<',selvar,'>'
       gx%bmperr=8888; goto 1000
!--------------------------------------------------------------------
! chemical potential for a component
    case('MU  ')
       if(n1.le.0) then
          write(*,*)'tqgetv 17: component number must be positive'
          gx%bmperr=8888; goto 1000
       endif
!       call get_component_name(n1,name,ceq)
!       if(gx%bmperr.ne.0) goto 1000
       statevar=stavar(1:2)//'('//cnam(n1)(1:len_trim(cnam(n1)))//') '
!       write(*,*)'tqgetv 4: ',statevar(1:len_trim(statevar))
! we must use index myvalue(1) as the subroutine expect a single variable
       call get_state_var_value(statevar,values(1),encoded,ceq)
!--------------------------------------------------------------------
! Amount of moles of components in a phaase
    case('NP  ')
       if(n1.lt.0) then
! all phases
          statevar='NP(*)'
		 
!          write(*,*)'tqgetv 1: ',mjj,statevar(1:len_trim(statevar))
! hopefully this returns all composition sets for all phases ... YES!
          call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
		 
! this output gives the amounts for all compsets of a phase sequentially
! but here we want them in phase tuple order
! the second argument is the number of values for each phase, here is 1 but
! it can be for example compositions, then it should be number of components

        call sortinphtup(n3,1,values)
		    
       else
! NOTE in this case n1 is a phase tuple index
!          ics=mod(n1,10)
!          nph=n1/10
!          if(nph.eq.0 .or. ics.eq.0) then
!             write(*,*)'You must use extended phase index'
!             gx%bmperr=8887; goto 1000
!          endif
!          call get_phase_name(nph,ics,name)
          call get_phase_name(phcs(n1)%phaseix,phcs(n1)%compset,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar='NP('//name(1:len_trim(name))//') '
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       endif
!--------------------------------------------------------------------
! Mole or mass fractions
    case('X   ','W   ')
!       write(*,*)'tqgetv: ',n1,n2,n3
       if(n2.eq.0) then
          if(n1.lt.0) then
! mole ´fraction of all components, no phase specification
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
             statevar=stavar(1:1)//'('//cnam(n1)(1:len_trim(cnam(n1)))//')'
!             write(*,*)'tqgetv 4: ',statevar(1:len_trim(statevar))
             call get_state_var_value(statevar,values(1),encoded,ceq)
          endif
       elseif(n1.lt.0) then
!........................................................
! for all phases one or several components
          if(n2.lt.0) then
! this means all components all phases
             statevar=stavar(1:1)//'(*,*) '
!             write(*,*)'tqgetv 5: ',mjj,statevar(1:len_trim(statevar))
             call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
! this output gives the composition for all compsets of a phase sequentially
! but we want them in phase tuple order
! ??             call sortinphtup(n3,,values)
          else
! a single component in all phases. n2 must not be zero
!             call get_component_name(n2,name,ceq)
!             if(gx%bmperr.ne.0) goto 1000
             if(n2.le.0 .or. n2.ge.noel()) then
                write(*,*)'No such component'
                goto 1000
             endif
             statevar=stavar(1:1)//'(*,'//cnam(n2)(1:len_trim(cnam(n2)))//')'
!             write(*,*)'tqgetv 6: ',mjj,statevar(1:len_trim(statevar))
             call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
! this output gives the composition for all compsets of a phase sequentially
! but we want them in phase tuple order
             ics=noel()
             call sortinphtup(n3,ics,values)
          endif
       elseif(n2.lt.0) then
! this means all components in one phase
! NOTE in this case n1 is a phasetuple index
!          ics=mod(n1,10)
!          nph=n1/10
!          if(nph.eq.0 .or. ics.eq.0) then
!             write(*,*)'You must use extended phase index'
!             gx%bmperr=8887; goto 1000
!          endif
!          call get_phase_name(nph,ics,name)
!          write(*,*)'Phase: ',phcs(n1)%phaseix,phcs(n1)%compset
          call get_phase_name(phcs(n1)%phaseix,phcs(n1)%compset,name)
          if(gx%bmperr.ne.0) goto 1000
! added for composition sets
!          if(ics.gt.1) then
!             name=name//'#'//char(ichar('0')+ics)
!          endif
          statevar=stavar(1:1)//'('//name(1:len_trim(name))//',*) '
!          write(*,*)'tqgetv 7: ',mjj,statevar(1:len_trim(statevar))
          call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
       else
! one component (n2) of one phase (n1)
! NOTE in this case n1 is 10*phase number + composition set number
!          ics=mod(n1,10)
!          nph=n1/10
!          if(nph.eq.0 .or. ics.eq.0) then
!             write(*,*)'You must use extended phase index'
!             gx%bmperr=8887; goto 1000
!          endif
!          call get_phase_name(nph,ics,name)
          call get_phase_name(phcs(n1)%phaseix,phcs(n1)%compset,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=stavar(1:1)//'('//name(1:len_trim(name))//','
          call get_component_name(n2,name,ceq)
          if(gx%bmperr.ne.0) goto 1000
          statevar(len_trim(statevar)+1:)=name(1:len_trim(name))//') '
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
! NOTE in this case n1 is 10*phase number + composition set number
!          ics=mod(n1,10)
!          nph=n1/10
!          if(nph.eq.0 .or. ics.eq.0) then
!             write(*,*)'You must use extended phase index'
!             gx%bmperr=8887; goto 1000
!          endif
!          call get_phase_name(nph,ics,name)
          call get_phase_name(phcs(n1)%phaseix,phcs(n1)%compset,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=statevar(1:ki)//'('//name(1:len_trim(name))//') '
!          call get_state_var_value(statevar,values(ics),encoded,ceq)
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       else
! Total volume
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       endif
!--------------------------------------------------------------------
! TEMPERATURE
    case('T   ')
! phase specifier not allowed
       if(norm(1:1).ne.' ') then
          statevar='T'//norm
          ki=2
       else
          statevar='T '
          ki=1
       endif
	   call get_state_var_value(statevar,values(1),encoded,ceq)	
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
! NOTE in this case n1 is 10*phase number + composition set number
!          ics=mod(n1,10)
!          nph=n1/10
!          if(nph.eq.0 .or. ics.eq.0) then
!             write(*,*)'You must use extended phase index'
!             gx%bmperr=8887; goto 1000
!          endif
!          write(*,*)'tqgetv 2: ',nph,ics
!          call get_phase_name(nph,ics,name)
          call get_phase_name(phcs(n1)%phaseix,phcs(n1)%compset,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=statevar(1:ki)//'('//name(1:len_trim(name))//') '
!          write(*,*)'tqgetv 3: ',statevar
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       else
! Total Gibbs energy 
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
! H for a specific phase
! NOTE in this case n1 is 10*phase number + composition set number
!          ics=mod(n1,10)
!          nph=n1/10
!          if(nph.eq.0 .or. ics.eq.0) then
!             write(*,*)'You must use extended phase index'
!             gx%bmperr=8887; goto 1000
!          endif
!          write(*,*)'tqgetv 2: ',nph,ics
!          call get_phase_name(nph,ics,name)
			call get_phase_name(phcs(n1)%phaseix,phcs(n1)%compset,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar=statevar(1:ki)//'('//name(1:len_trim(name))//') '
!          write(*,*)'tqgetv 3: ',statevar
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       else
! Total Enthalpy 
          call get_state_var_value(statevar,values(1),encoded,ceq)
          n3=1
       endif
!--------------------------------------------------------------------	   
!--------------------------------------------------------------------
! Mobilities
    case('MQ   ')
!       ics=mod(n1,10)
!       nph=n1/10
!       if(nph.eq.0 .or. ics.eq.0) then
!          write(*,*)'You must use extended phase index: 10*phase+compset'
!          gx%bmperr=8887; goto 1000
!       endif
!       call get_phase_name(nph,ics,name)
       call get_phase_name(phcs(n1)%phaseix,phcs(n1)%compset,name)
       if(gx%bmperr.ne.0) goto 1000
       statevar=stavar(1:len_trim(stavar))//'('//name(1:len_trim(name))//')'
!       write(*,*)'statevar: ',statevar
       call get_state_var_value(statevar,values(1),encoded,ceq)
!--------------------------------------------------------------------
! Second derivatives of the Gibbs energy of a phase
    case('D2G   ')
!       ics=mod(n1,10)
!       nph=n1/10
!       if(nph.eq.0 .or. ics.eq.0) then
!          write(*,*)'You must use extended phase index: 10*phase+compset'
!          gx%bmperr=8887; goto 1000
!       endif
!       write(*,*)'D2G 1: ',nph,ics
!       call get_phase_compset(nph,ics,lokph,lokcs)
       call get_phase_compset(phcs(n1)%phaseix,phcs(n1)%compset,lokph,lokcs)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'D2G 2: ',lokph,lokcs
! this gives wrong myvalue!!
!       n3=ceq%phase_varres(lokcs)%ncc
       n3=size(ceq%phase_varres(lokcs)%yfr)
!       write(*,*)'D2G 3: ',n3
       kj=(n3*(n3+1))/2
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
    call get_phase_data(phcs(n1)%phaseix,phcs(n1)%compset,&
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
    call set_constitution(phcs(n1)%phaseix,phcs(n1)%compset,&
         yfra,extra,ceq)
1000 continue
    return
  end subroutine tqsphc1

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqcph1(n1,n2,n3,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,ceq)
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
    integer n1,n2,n3
    double precision gtp(6),dgdy(*),d2gdydt(*),d2gdydp(*),d2gdy2(*)
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer ij,lokres,nofc
!    write(*,*)'tqcph1 1: ',ceq%eqname
!    write(*,*)'tqcph1 2',phcs(n1)%phaseix,phcs(n1)%compset
!----------------------------------------------------------------------
! THIS IS NO EQUILIBRIUM, JUST G AND DERIVATIVES FOR CURRENT T, P AND Y
    call calcg(phcs(n1)%phaseix,phcs(n1)%compset,n2,lokres,ceq)
!----------------------------------------------------------------------
!    write(*,*)'tqcph1 3A',lokres,gx%bmperr
! this should work but gave segmentation fault, find this a more cumbersum way
    n3=size(ceq%phase_varres(lokres)%yfr)
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

!\begin{verbatim}
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
!    write(*,*)'tqcph1 2',phcs(n1)%phaseix,phcs(n1)%compset
!----------------------------------------------------------------------
! THIS IS NO EQUILIBRIUM, JUST G AND DERIVATIVES FOR CURRENT T, P AND Y
    call calcg(phcs(n1)%phaseix,phcs(n1)%compset,n2,lokres,ceq)
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
    ceq=>eqlista(n1)
    call delete_equilibrium(name,ceq)
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
    call enter_equilibrium(name,n1)
    if(gx%bmperr.ne.0) goto 1000
    newceq=>eqlista(n1)
    call copy_equilibrium(newceq,name,ceq)
1000 continue
    return
  end subroutine tqcceq

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
  subroutine reset_conditions(cline,ceq)
!reset any condition on temperature  
    implicit none
	character cline*24
    type(gtp_equilibrium_data), pointer :: ceq 
!\end{verbatim}	
 
    integer ip
  
	ip=0
!	write(*,*) cline
    call set_condition(cline,ip,ceq)
1000 continue	
	return
  end subroutine reset_conditions
 
!\begin{verbatim}
  subroutine Change_Status_Phase(myname,nystat,myval,ceq)
!PHFIXED=2
!PHENTERED=0
    implicit none
	character myname*24
	integer nystat
	double precision myval
	type(gtp_equilibrium_data), pointer :: ceq 
!\end{verbatim}	
	integer iph,ics
	
	
	call find_phase_by_name(myname,iph,ics)
	call change_phase_status(iph,ics,nystat,myval,ceq)

1000 continue	
	return
  end subroutine Change_Status_Phase
 

end MODULE LIBOCTQ
