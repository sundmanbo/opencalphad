!
! Minimal TQ interface.
!
! To compile and link this with an application one must first compile
! and form a library of most OC subroutines (oclib.a) and to copy this and
! the "mov" files from this compilation to the folder with this library
!
! NOTE that for identification of phases this TQ interface often use an
! EXTENDED PHASE INDEX created by  10*phase number + comp.set number
! This extended phase index is used in many of the subroutines.
! Also for constituents an EXTENDED CONSTITUENT INDEX is sometimes used 
! and equal to 10*species_number + sublattice
!
! 140128 BOS added D2G and phase specific V and G
! 140128 BOS added possibility to calculate without evoking grid minimizer
! 140125 BOS Changed name to liboctq
! 140123 BOS Added ouput of MQ G, V and normalized
!
! The name of this librqry
module liboctq
!
! access to main OC library for equilibrium calculations
  use liboceq
!
  implicit none
!
!
contains
!
!\begin{verbatim}
  subroutine tqini(n,ceq)
! initiate workspace
    implicit none
    integer n ! Not nused
    type(gtp_equilibrium_data), pointer :: ceq ! EXIT: current equilibrium
!\end{verbatim}
! these should eventually provide linits and defaults
    integer intv(10)
    double precision dblv(10)
    intv(1)=-1
    call init_gtp(intv,dblv)
1000 continue
    return
  end subroutine tqini

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqrfil(filename,ceq)
! read TDB file
    implicit none
    character*(*) filename  ! IN: database filename
    character ellista(10)*2  ! dummy
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
    integer :: i2
!\end{verbatim}
! second argument 0 means ellista is ignored, all element read
    call readtdb(filename,0,ellista)
    ceq=>firsteq
    i2 = transfer(firsteq,i2)
1000 continue
    return
  end subroutine tqrfil

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqrpfil(filename,nel,selel,ceq)
! read TDB file with selection of elements
    implicit none
    character*(*) filename  ! IN: database filename
    integer nel,i
    character selel(*)*2  ! IN: elements to be read from the database
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    call readtdb(filename,nel,selel)
    ceq=>firsteq
1000 continue
    return
  end subroutine tqrpfil
 
!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgcom(n,components,ceq)
! get system components
    implicit none
    integer n ! EXIT: number of components
    character*24, dimension(*) :: components ! EXIT: names of components
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
    integer iz
    character elname*24,refs*24
    double precision a1,a2,a3
    n=noel()
    do iz=1,n
       components(iz)=' '
       call get_element_data(iz,components(iz),elname,refs,a1,a2,a3)
    enddo
1000 continue
    return
  end subroutine tqgcom

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgnp(n,ceq)
! get number of phases.
! NOTE the number of composition sets is given in global array phcs
    implicit none
    integer n    !EXIT: n is number of phases
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    n=noph()
1000 continue
    return
  end subroutine tqgnp

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpn(n,phasename,ceq)
! get name of phase n, 
! NOTE: n is phase number, not extended phase index
    implicit none
    integer n  !IN: phase number (not extended phase index)
    character phasename*(24) !EXIT: phase name
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
!    integer ics
!    ics=1
    call get_phase_name(n,1,phasename)
1000 continue
    return
  end subroutine tqgpn

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpi(n,phasename,ceq)
! get index of phase phasename
    implicit none
    integer n  !EXIT: phase number
    character phasename*(*) !IN: phase name
    type(gtp_equilibrium_data), pointer :: ceq !IN: current equilibrium
!\end{verbatim}
    call find_phase_by_name(phasename,n,1)
1000 continue
    return
  end subroutine tqgpi

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgpcn(n,c,constituentname,ceq)
! get name of consitutent c in phase n
    implicit none
    integer n !IN: phase nubmer
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
  subroutine tqsetc(stavar,n1,n2,value,cnum,ceq)
! set condition
! stavar is state variable as text
! n1 and n2 are auxilliary indices
! value is the value of the condition
! cnum is returned as an index of the condition.
! to remove a condition the value sould be equial to RNONE ????
! when a phase indesx is needed it should be 10*nph + ics
! SEE TQGETV for doucumentation of stavar etc.
    implicit none
    integer n1 ! IN: 0 or extended phase index: 10*phase_number+comp.set
!                or component number
    integer n2 ! IN: 0 or component number
    integer cnum !EXIT: sequential number of this condition
    character stavar*(*) !IN: character with state variable symbol
    double precision value !IN: value of condition
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
    integer ip
    character cline*60,name*24,selvar*4
    cline=' '
    selvar=stavar
    call capson(selvar)
    select case(selvar)
    case default
       write(*,*)'Condition not implemented or illegal',stavar
       gx%bmperr=8888; goto 1000
    case('T   ','P   ','N   ')
       write(cline,110)selvar(1:1),value
110    format(' ',a,'=',E15.8)
    case('X   ','W   ')
       call get_component_name(n1,name,ceq)
       if(gx%bmperr.ne.0) goto 1000
       write(cline,120)selvar(1:1),name(1:len_trim(name)),value
120    format(1x,a,'(',a,')=',1pE15.8)
! MORE CONDITIONS WILL BE ADDED ...
    end select
!    write(*,*)'tqsetc: ',cline(1:len_trim(cline))
    ip=1
    call set_condition(cline,ip,ceq)
1000 continue
    return
  end subroutine tqsetc

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
    type(gtp_equilibrium_data), pointer :: ceq  !IN: current equilibrium
!\end{verbatim}
! mode=1 means start values using global gridminimization
    mode=1
    if(n1.lt.0) then
! this means calculate without grid minimuzer
       write(*,*)'No grid minimizer'
       mode=0
    endif
    call calceq2(mode,ceq)
1000 continue
    return
  end subroutine tqce

!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\

!\begin{verbatim}
  subroutine tqgetv(stavar,n1,n2,n3,values,ceq)
! get equilibrium results using state variables
! stavar is the state variable IN CAPITAL LETTERS with indices n1 and n2 
! n3 at the call is the dimension of values, changed to number of values
! value is the calculated value, it can be an array with n3 values.
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
    character statevar*60,encoded*60,name*24,selvar*4,norm*4
! mjj should be the dimension of the array values ...
    mjj=n3
    selvar=stavar
    call capson(selvar)
! for state variables like MQ&FE remove the part from & in the select
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
       write(*,*)'Unknown state variable: ',stavar(1:20),'>:<',selvar
       gx%bmperr=8888; goto 1000
!--------------------------------------------------------------------
! chemical potential for a component
    case('MU  ')
       if(n1.le.0) then
          write(*,*)'tqgetv 17: component number must be positive'
          gx%bmperr=8888; goto 1000
       endif
       call get_component_name(n1,name,ceq)
       if(gx%bmperr.ne.0) goto 1000
       statevar=stavar(1:2)//'('//name(1:len_trim(name))//') '
!       write(*,*)'tqgetv 4: ',statevar(1:len_trim(statevar))
! we must use index value(1) as the subroutine expect a single variable
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
       else
! NOTE in this case n1 is 10*phase number + composition set number
          ics=mod(n1,10)
          nph=n1/10
          if(nph.eq.0 .or. ics.eq.0) then
             write(*,*)'You must use extended phase index'
             gx%bmperr=8887; goto 1000
          endif
          call get_phase_name(nph,ics,name)
          if(gx%bmperr.ne.0) goto 1000
          statevar='NP('//name(1:len_trim(name))//') '
! why should one return the value in value(ics) ??
!          call get_state_var_value(statevar,values(ics),encoded,ceq)
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
! mole fraction for the actual state variable, n1 and n2 not used
! just check for wildcard
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
             call get_component_name(n1,name,ceq)
             if(gx%bmperr.ne.0) goto 1000
             statevar=stavar(1:1)//'('//name(1:len_trim(name))//') '
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
          else
! a single component in all phases. n2 cannot be zero
             call get_component_name(n2,name,ceq)
             if(gx%bmperr.ne.0) goto 1000
             statevar=stavar(1:1)//'(*,'//name(1:len_trim(name))//') '
!             write(*,*)'tqgetv 6: ',mjj,statevar(1:len_trim(statevar))
             call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
          endif
       elseif(n2.lt.0) then
! this means all components in one phase
! NOTE in this case n1 is 10*phase number + composition set number
          ics=mod(n1,10)
          nph=n1/10
          if(nph.eq.0 .or. ics.eq.0) then
             write(*,*)'You must use extended phase index'
             gx%bmperr=8887; goto 1000
          endif
          call get_phase_name(nph,ics,name)
          if(gx%bmperr.ne.0) goto 1000
! added for composition sets
          if(ics.gt.1) then
             name=name//'#'//char(ichar('0')+ics)
          endif
          statevar=stavar(1:1)//'('//name(1:len_trim(name))//',*) '
!          write(*,*)'tqgetv 7: ',mjj,statevar(1:len_trim(statevar))
          call get_many_svar(statevar,values,mjj,n3,encoded,ceq)
       else
! one component (n2) of one phase (n1)
! NOTE in this case n1 is 10*phase number + composition set number
          ics=mod(n1,10)
          nph=n1/10
          if(nph.eq.0 .or. ics.eq.0) then
             write(*,*)'You must use extended phase index'
             gx%bmperr=8887; goto 1000
          endif
          call get_phase_name(nph,ics,name)
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
          ics=mod(n1,10)
          nph=n1/10
          if(nph.eq.0 .or. ics.eq.0) then
             write(*,*)'You must use extended phase index'
             gx%bmperr=8887; goto 1000
          endif
          call get_phase_name(nph,ics,name)
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
          ics=mod(n1,10)
          nph=n1/10
          if(nph.eq.0 .or. ics.eq.0) then
             write(*,*)'You must use extended phase index'
             gx%bmperr=8887; goto 1000
          endif
!          write(*,*)'tqgetv 2: ',nph,ics
          call get_phase_name(nph,ics,name)
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
! Mobilities
    case('MQ   ')
       ics=mod(n1,10)
       nph=n1/10
       if(nph.eq.0 .or. ics.eq.0) then
          write(*,*)'You must use extended phase index: 10*phase+compset'
          gx%bmperr=8887; goto 1000
       endif
       call get_phase_name(nph,ics,name)
       if(gx%bmperr.ne.0) goto 1000
       statevar=stavar(1:len_trim(stavar))//'('//name(1:len_trim(name))//')'
!       write(*,*)'statevar: ',statevar
       call get_state_var_value(statevar,values(1),encoded,ceq)
!--------------------------------------------------------------------
! Second derivatives of the Gibbs energy of a phase
    case('D2G   ')
       ics=mod(n1,10)
       nph=n1/10
       if(nph.eq.0 .or. ics.eq.0) then
          write(*,*)'You must use extended phase index: 10*phase+compset'
          gx%bmperr=8887; goto 1000
       endif
!       write(*,*)'D2G 1: ',nph,ics
       call get_phase_compset(nph,ics,lokph,lokcs)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,*)'D2G 2: ',lokph,lokcs
! this gives wrong value!!
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

end MODULE LIBOCTQ

