! Hillert Minimizer as implemented by Sundman (HMS)
! Based on Mats Hillert paper in Physica 1981 and Bo Janssons thesis 1984
!
MODULE matsmin
!
! Copyright 2012-2013, Bo Sundman, France
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!---------------------------
!
! Done:
! 120121: single liquid phase binary Cr-Fe with mass balance conditions works !
!         bcc-C-Cr-Fe and sigma Cr-Fe works too.
! 120210: Single phase equilibria with 3 componets works for gas, intermetall,
!         interstial phases, miscibility gap in bcc.  Working on finding
!         the set of stable phases in 6 component steels
! 120220: Bug in gridmin fixed when restoring gridpoints. Global gidmin did
!         not find the Fcc-MC composition set. CrFe binary test case works
!         including phase changes.  Phase set changes in steel test case
!         does not work, neither having a closed system with several phases.
! 120222: Changed all massbalance conditions to be N(i)=value, not x(i)=value
!         together with N=1.  This gives simpler derivatives and now it
!         converges in most cases, even 6 component steel with 40 phases. WoW
! 120226: Modified extract_massbalance in PMOD so now one can calculate with
!         mixed massbalance conditions, N, X, B, W. Result presented with N=1
! 120301: Fixing minor problems like problem when 2 B2 phases want to be stable
!         in the Re-V system.  Non-convergece forced when correction in
!         constituent fractions of stable phaes increases (also in Re-V).
!         Created a new structure to contain global variables for the
!         equilibrium calculation. Not many at present but needed when more
!         conditions will be used.
! 120305: Charge balance added.  Rater slow convergence and bad precision in
!         chemical potential values.
!
! >> It is high time for a revision of the whole system and to discuss the
! >> data structure for the post processing and graphics.
!
!-----------------------------
!
! To be done before first release in January 2013:
! - adding variable T and P
! - adding fix phase condition
! - adding chemical potential condition
! - calculating partial derivatives (Cp, thermal expansion etc)
! - stability check (eigenvalues)
! - conditions representing normallized properties like x, VM etc.
! To be done later outside this module:
! - step
! - map
! - graphical postprocessing
! - assessments
! 
!------------------------------
!
  use general_thermodynamic_package
!
  implicit none
!  
!  double precision, dimension(2) :: tpvalx
  character*8, parameter :: hmsversion='HMS-1.00'
!
!\begin{verbatim}
  TYPE meq_setup
! one structure of this type is created when an equilibrium calculation
! is started and it holds all global data needed for handling the
! calculation of an equilibrium.  The phase specific data is in meq_data
! nv: initial guess of number of stable phases
! nphase: total number of phases and composition sets
! nstph: current number of stable phases
! dormlink: is start of list of phases temporarily set dormant
! typesofcond: types of conditions, =1 only massbal, =2 any conditions
     integer nv,nphase,nstph,dormlink,noofits
     integer nrel,typesofcond,maxsph,nfixmu,nfixph
! component numbers of fixed potentials
     integer, dimension(:), allocatable :: mufixel
     integer, dimension(:), allocatable :: mufixref
     double precision, dimension(:), allocatable :: mufixval
! fix phases
     integer, dimension(:,:), allocatable :: fixph
     double precision, dimension(:), allocatable :: fixpham
! iphl, icsl: phase and composition sets of intial guess of stable phases
! aphl: initial guess of amount of each stable phase
     integer iphl(maxel+2),icsl(maxel+2)
     double precision aphl(maxel+2)
! this is because I tried to scale the total amount of phases during iterations
     double precision antot
! stphl: current list of stable phases, value is index in phr array
     integer, dimension(maxel+2) :: stphl
! current values of chemical potentials
     double precision, dimension(:), allocatable :: curmu
! if variable T and P these are TRUE, otherwise FALSE
     logical tpindep(2)
     double precision deltat,deltap
! normallizing of amount changes
     double precision nnorm
! information about conditions should be stored here.  Note that conditions
! may change during STEP and MAP
  end TYPE meq_setup
!\end{verbatim}
!
!\begin{verbatim}
  TYPE meq_data
! parts of this structure should maybe be part of the gtp_equilibrium_data
! it contains phase specific results from various subroutines during
! equilibrium calculation
! iph: phase number
! ics: constituent set number
! idim: the dimension of phase matrix, 
! ncc: the number of constituents
! stable: is 1 for a stable phase
! xdone: set to 1 for stoichiometric phases after calculating xmol first time
! dormlink: used to link phases that temporarily been set dormant
     integer iph,ics,idim,stable,ncc,xdone,dormlink
! true if phase is fix (no variable amount)
     integer phasestatus
! inverted phase matrix
     double precision, dimension(:,:), allocatable :: invmat
! mole fractions of components and their sum
     double precision, dimension(:), allocatable :: xmol
     double precision :: sumxmol,sumwmol
! total moles and mass of components and their sum
!     double precision, dimension(:), allocatable :: nmol
!     double precision :: sumnmol,sumbmol
! Derivatives of moles of component wrt all constituent fractions of the phase
     double precision, dimension(:,:), allocatable :: dxmol
! link to phase_varres record
     TYPE(gtp_phase_varres), pointer :: curd
! value of amount and driving force at previous iteration
     double precision prevam, prevdg
! iteration when phase was added/removed
     integer itadd, itrem
! chargebal is 1 if external charge balance needed
     integer chargebal
     double precision charge !redundant
  end TYPE meq_data
!\end{verbatim}
!  
! IMPORTANT
! phase_varres(lokcs)%amfu is the number of formula units of the phase
! phase_varres(lokcs)%netcharge is the total charge  of the phase
! phase_varres(lokcs)%abnorm(1) is the number of real atoms per formula unit
! (may vary with composition like in (Fe,Cr,...)(Va,C,N,...) )
! phase_varres(lokcs)%abnorm(2) is the mass per formula unit
! NOTE: abnorm(1) and abnorm(2) is set by call to set_constitution)
! There has been a lot of wiggling around this and probably some errors
!
!
CONTAINS
  
!\begin{verbatim}
  subroutine calceq2(mode,ceq)
! calculates the equilibrium for the given set of conditions
! mode=0 means no global minimization
! ceq is a datastructure with all relevant thermodynamic data
    implicit none
    integer mode
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    TYPE(gtp_condition), pointer :: condition,lastcond
    TYPE(meq_setup), target :: meqrec
! conditions on T and P and mole fractions of components
    double precision, dimension(2) :: tpval
    double precision, dimension(maxel) :: xknown,vmu
    double precision xxx,antot,finish2,starting,cvalue
    logical alltid,gridtest
! for global minimization (change maybe to allocate dynamically)
    integer, dimension(maxph) :: nyphl
    double precision, dimension(maxconst) :: yarr
    integer starttid,endoftime,np,iph,ics,jph,lokph,lokcs,mode2,ierr
    integer mostcon,mph,nph,nvf
    integer, parameter :: mmu=5
    integer what,mjj,ij,kjj,cmix(10),cmode,ctype,mufixel(mmu),mufixref(mmu)
    integer fixph(2,maxel),oldorder(mmu),kst
    double precision fixpham(maxel)
    character statevar*40,encoded*60,name*24
!
    write(*,*)"Entering calceq2"
    if(gx%bmperr.ne.0) then
       write(*,*)'Error code set before calling grid minimizer',gx%bmperr
       goto 1000
    endif
    mode2=mode
    call cpu_time(starting)
    call system_clock(count=starttid)
! extract conditions
    call extract_massbalcond(tpval,xknown,antot,ceq)
     if(gx%bmperr.ne.0) then
! error 4144 means wrong number of conditions, other codes mean other things
!       if(gx%bmperr.eq.4151) goto 1000
       if(gx%bmperr.eq.4144) goto 1000
!       write(*,*)'Not only massbalance conditions: '
       ierr=gx%bmperr; gx%bmperr=0
       gridtest=.true.
       meqrec%typesofcond=2
    else
       meqrec%antot=antot
! no need for final grid minimizer as we will do one as start
       gridtest=.false.
       meqrec%typesofcond=1
    endif
! allocate curmu for chemical potentials
    meqrec%nrel=noel()
    allocate(meqrec%curmu(meqrec%nrel))
! set some initial values
    meqrec%maxsph=noel()+2
    meqrec%nfixph=0
    meqrec%nfixmu=0
    meqrec%tpindep=.TRUE.
    meqrec%nnorm=one
! now we calculate maxsph, nfixmu and maybe other things for later
    lastcond=>ceq%lastcondition
    condition=>lastcond
    cmix=0
    np=0
    mjj=0
! set default values
!    write(*,69)tpval,ceq%tpval
!69  format('T&P: ',4(1pe12.4))
    tpval(1)=ceq%tpval(1)
    tpval(2)=ceq%tpval(2)
    ceq%rtn=globaldata%rgas*tpval(1)
!---------------- loop
! loop through all conditions, end when the pointer condition is empty
! loop to investigate conditions, apply_condition:value in pmod25D.F90
70  continue
       cmode=-1
       condition=>condition%next
       mjj=mjj+1
       call apply_condition_value(condition,cmode,ctype,cvalue,cmix,ceq)
       if(gx%bmperr.ne.0) goto 1000
! cmix(1)=0 for inactive conditions
! cmix(1)=1 fix T, =2, fix P, =3 fix MU/AC/LNAC, =4 fix phase, =5 anything else
! if condition on T, P, potential or fix phase reduce maxsph
       select case(cmix(1))
       case default
          if(.not.associated(condition,lastcond)) goto 70
       case(1) ! fix T
          if(cvalue.le.zero) then
             write(*,*)'Condition on T must be larger than zero'
             gx%bmperr=7777; goto 1000
          endif
          meqrec%maxsph=meqrec%maxsph-1
          meqrec%tpindep(1)=.FALSE.
          ceq%tpval(1)=cvalue
       case(2) ! fix P
          if(cvalue.le.zero) then
             write(*,*)'Condition on P must be larger than zero'
             gx%bmperr=7777; goto 1000
          endif
          meqrec%maxsph=meqrec%maxsph-1
          meqrec%tpindep(2)=.FALSE.
          ceq%tpval(2)=cvalue
!-------------------------
       case(3) ! (MU,AC,LNAC) in cmix(2)=3,4,5
! The component is in cmix(3) and reference state in cmix(4)
! Handling of the reference state ignored at present
          np=np+1
          if(np.gt.mmu) then
             write(*,*)'Max conditions on potentials is ',mmu
             gx%bmperr=7777; goto 1000
          endif
          mufixel(np)=cmix(3)
          mufixref(np)=cmix(4)
! temporarily use yarr for something else
          if(cmix(2).eq.3) then
! Divide MU by RT
             yarr(np)=cvalue/ceq%rtn
          elseif(cmix(2).eq.4) then
! AC converted to chemical potential/RT
             if(cvalue.le.zero) then
                write(*,*)'Conditions on activity must be larger than zero'
                gx%bmperr=7777; goto 1000
             endif
             yarr(np)=LOG(cvalue)
          else
! LNAC is MU/RT which is the value used during minimization
             yarr(np)=cvalue
          endif
!          write(*,*)'Chemical potential condition: ',yarr(np)
          meqrec%maxsph=meqrec%maxsph-1
!             write(*,72)'MM, chemp: ',cmix(1),cmix(2),cmix(3),cvalue
!72           format(a,3i3,1pe12.4)
!-------------------------
       case(4) ! fix phase
          meqrec%nfixph=meqrec%nfixph+1
          fixph(1,meqrec%nfixph)=cmix(2)
          fixph(2,meqrec%nfixph)=cmix(3)
          fixpham(meqrec%nfixph)=cvalue
!          write(*,*)'Fix phase condition: ',cmix(2),cmix(3),cvalue
! debug output of fix phase composition
!          call calc_phase_mol(cmix(1),yarr,ceq)
!          write(*,83)'fix phasse: ',cmix(1),(yarr(mjj),mjj=1,noel())
!83        format(a,i2,10f7.4)
       end select
       if(.not.associated(condition,lastcond)) goto 70
!------------------------- end loop of conditions
!       write(*,*)'variable potentials, max variable phases: ',&
!               noel()-cmix(2),meqrec%maxphases
    meqrec%nfixmu=np
    if(np.gt.0) then 
       allocate(meqrec%mufixel(np))
       allocate(meqrec%mufixref(np))
       allocate(meqrec%mufixval(np))
       if(np.gt.1) then
! sort components in increasing order to simplify below
          call sortin(mufixel,np,oldorder)
          do mjj=1,np
             nvf=mufixel(mjj)
             meqrec%mufixel(mjj)=nvf
             meqrec%mufixref(mjj)=mufixref(oldorder(mjj))
             meqrec%mufixval(mjj)=yarr(oldorder(mjj))
! copy fixed chemical potential to curmu also
             meqrec%curmu(nvf)=yarr(oldorder(mjj))
             ceq%complist(nvf)%chempot(1)=meqrec%curmu(nvf)*ceq%rtn
!             write(*,*)'Fix lnac: ',nvf,meqrec%curmu(nvf)
          enddo
       else
          nvf=mufixel(1)
          meqrec%mufixel(1)=nvf
          meqrec%mufixref(1)=mufixref(1)
          meqrec%mufixval(1)=yarr(1)
! also copy fixed chemical potential to curmu
          meqrec%curmu(nvf)=yarr(1)
          ceq%complist(nvf)%chempot(1)=meqrec%curmu(nvf)*ceq%rtn
!          write(*,*)'Fix lnac: ',nvf,meqrec%curmu(nvf)
       endif
    endif
    if(meqrec%nfixph.gt.0) then
       allocate(meqrec%fixph(2,meqrec%nfixph))
       allocate(meqrec%fixpham(meqrec%nfixph))
       if(np.gt.1) then
! ?? sort phases in increasing order to simplify below
          write(*,*)'Cannot handle two fix phases ... '
          stop 'sorry'
       endif
       do mjj=1,meqrec%nfixph
          meqrec%fixph(1,mjj)=fixph(1,mjj)
          meqrec%fixph(2,mjj)=fixph(2,mjj)
          meqrec%fixpham(mjj)=fixpham(mjj)
       enddo
    endif
!----------------------------
!    call list_conditions(kou,ceq)
! skip global gridminimizer if bit set
    if(mode2.eq.0 .or. btest(globaldata%status,GSNOGLOB)) goto 110
! skip global gridminimizer if one component
    if(meqrec%nrel.eq.1) goto 110
! bad rstults for ge H2O1 gas case without global minimization
!    alltid=.TRUE.
!    do nph=1,noph()
!       do ics=1,noofcs(nph)
! skip global minimization if a phase+compset is fix
!          if(test_phase_status(nph,ics,xxx,ceq).eq.2) goto 110
!       enddo
! set alltid to FALSE if phase is non-ideal
!       if(.not.phase_bit(nph,PHID)) alltid=.FALSE.
!    enddo
! skip global gridminimizer if all phases ideal or fixsed composition
!    if(alltid) goto 110
! Try global gridminimization.  Returned values are:
! nv is number of stable phase, iphl, icsl list of stable  phases, aphl amounts
! nyphl(j) is number of constituent fractions in phase j, yarr are the 
! constituent fractions, vmu the chemical potentials
! THIS CALL MAY CREATE NEW COMPOSITION SETS unless GSNOACS set.
! loop through all phases and set amount=0 and CSSTABLE off
    ij=1
    call todo_before(ij,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(meqrec%typesofcond.eq.1) then
! with only massbalance condition make a global grid minimization
       call global_gridmin(1,tpval,xknown,meqrec%nv,meqrec%iphl,meqrec%icsl,&
            meqrec%aphl,nyphl,yarr,vmu,ceq)
!       write(*,*)'back from gridmin'
       if(gx%bmperr.ne.0) then
! if global fails reset error code and try a default start set of phases
          if(gx%bmperr.ge.4000 .and. gx%bmperr.le.4220) then
             write(kou,102)gx%bmperr,bmperrmess(gx%bmperr)
102          format('Tried but failed to use grid minimization, error: ',i5/a)
          endif
! no initial gridmin, make a gridtest at the end
          gridtest=.true.
          gx%bmperr=0; goto 110
       else
          write(*,*)'Done grid minimization, ',meqrec%nv,' stable phase(s):'
       endif
! copy vmu to curmu
       do mjj=1,meqrec%nrel
          meqrec%curmu(mjj)=vmu(mjj)
       enddo
       goto 200
    endif
!    write(*,103)(meqrec%iphl(j),meqrec%icsl(j),meqrec%aphl(j),j=1,meqrec%nv)
103 format('Phases: ',6(i3,i2,1pe12.4))
!--------------------
! no global gridmin or we come here if gridtest finds a new stable phase
110 continue
    meqrec%nv=0
! at least one phase must be set as stable
    mostcon=0
    mph=0
    jph=0
    do iph=1,noph()
       do ics=1,noofcs(iph)
! status 1=entered, 2=fix, 3=dormant, 4=suspended, 5=hidden
          kst=test_phase_status(iph,ics,xxx,ceq)
          if(kst.gt.3) goto 115
          if(gx%bmperr.ne.0) goto 1000
          call get_phase_compset(iph,ics,lokph,lokcs)
!          if(ceq%phase_varres(lokcs)%amount(1).gt.zero) then
          if(ceq%phase_varres(lokcs)%amfu.gt.zero) then
             meqrec%nv=meqrec%nv+1
             meqrec%iphl(meqrec%nv)=iph
             meqrec%icsl(meqrec%nv)=ics
             meqrec%aphl(meqrec%nv)=ceq%phase_varres(lokcs)%amfu
!             meqrec%aphl(meqrec%nv)=ceq%phase_varres(lokcs)%amount(1)
          elseif(kst.ne.2) then
! set default constitution of non-fix phases with zero amount,
! one may have a very bad constitution from a failed global minimization,
! ignore any errors
!             write(*,*)'phase and status: ',iph,ics,kst
             call set_default_constitution(iph,ics,0,ceq)
             if(gx%bmperr.ne.0) gx%bmperr=0
          endif
       enddo
       call get_phase_variance(iph,nvf)
       if(nvf.gt.mostcon) then
          mostcon=nvf
          jph=iph
       endif
115    continue
    enddo
    if(meqrec%nv.eq.0) then
! no phase with positive amount, set phase with most constituents as stable
       if(jph.gt.0) then
          call get_phase_compset(jph,1,lokph,lokcs)
          ceq%phase_varres(lokcs)%amfu=one
          meqrec%nv=1
          meqrec%iphl(1)=jph
          meqrec%icsl(1)=1
          meqrec%aphl(1)=one
          write(*,*)'No gridminimization, selecting phase ',jph,' as stable'
       else
!          write(*,*)'No phase to set stable'
          gx%bmperr=4200; goto 1000
       endif
    else
       write(*,*)'No gridminimization, using current phase set',meqrec%nv
    endif
! copy ceq%complist% to curmu
    do mjj=1,meqrec%nrel
       meqrec%curmu(mjj)=ceq%complist(mjj)%chempot(1)/ceq%rtn
    enddo
!
! we must make sure the fix phases are in the initial list of stable phases
! the order does not matter, the phases will be sorted later
    addfixph: do mjj=1,meqrec%nfixph
       jph=1
       do while (jph.le.meqrec%nv)
          if(meqrec%iphl(jph).eq.meqrec%fixph(1,mjj) .and. &
             meqrec%icsl(jph).eq.meqrec%fixph(2,mjj)) then
! found fix phase as already stable, just store the amount
             meqrec%aphl(jph)=meqrec%fixpham(mjj)
             cycle addfixph
          endif
          jph=jph+1
       enddo
! add this phase, check that not too many stable phases ...
       if(meqrec%nv.eq.meqrec%maxsph) then
          write(*,*)'Too many stable phases'
          gx%bmperr=9998; goto 1000
       endif
       write(*,*)'Adding fix phase to stable phase set',&
            meqrec%fixph(1,mjj),meqrec%fixph(2,mjj)
       meqrec%nv=meqrec%nv+1
       meqrec%iphl(meqrec%nv)=meqrec%fixph(1,mjj)
       meqrec%icsl(meqrec%nv)=meqrec%fixph(2,mjj)
       meqrec%aphl(meqrec%nv)=meqrec%fixpham(mjj)
    enddo addfixph
! debug output of fix phase composition
!    call calc_phase_mol(1,yarr,ceq)
!    write(*,83)'at 200: ',cmix(1),(yarr(mjj),mjj=1,noel())
! zero start of link to phases set temporarily dormant ....
! >>>> why does this not work!!!, has to be set to zero inside meq_phaseset
    meqrec%dormlink=0
!-------------------------------
200 continue
!
! this routine varies the set of phases and the phase constitutions
! until the stable set is found for the given set of conditions.
!    write(*,*)'calling meq_phaseset'
    call meq_phaseset(meqrec,ceq)
    if(gx%bmperr.ne.0) goto 1000
!------------------------------------------------------
!    gridtest=.true.
    gridtest=.false.
    gridcheck: if(gridtest .and. btest(ceq%status,EQGRIDTEST)) then
! if conditions are not massbalance apply gridtest as we now know the
! overall composition.  First extract the molefractions
       write(*,*)'Testing solution with new grid'
       do mjj=1,meqrec%nrel
          call get_component_name(mjj,name,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'error getting component name ',mjj,gx%bmperr
             goto 1000
          endif
          statevar='X('//name(1:len_trim(name))//') '
          call get_state_var_value(statevar,xknown(mjj),encoded,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error extracting x to test the solution',gx%bmperr
             goto 1000
          endif
          statevar='MU('//name(1:len_trim(name))//') '
          call get_state_var_value(statevar,vmu(mjj),encoded,ceq)
          if(gx%bmperr.ne.0) then
             write(*,*)'Error extracting mu to test the solution',gx%bmperr
             goto 1000
          endif
       enddo
!       write(*,667)'mu: ',mjj,(vmu(ij),ij=1,mjj)
667    format(a,i2,7(1pe11.3))
       tpval(1)=ceq%tpval(1)
       tpval(2)=ceq%tpval(2)
! what=-1 means this is a test, what is changed if new phase should be stable
! use same variables as earler call. New composition sets can be entered
       what=-1
       call global_gridmin(what,tpval,xknown,meqrec%nv,meqrec%iphl,&
            meqrec%icsl,meqrec%aphl,nyphl,yarr,vmu,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error calculating grid to test the solution',gx%bmperr
          goto 1000
       endif
       if(what.gt.0) then
! solution incorrect, a gridpoint is below the chem.pot surface
          iph=what/10
          ics=mod(what,10)
          write(*,*)'Recalculate as another stable phase ',iph,ics
          goto 110
       endif
    endif gridcheck
!--------------------------------------------------
! Here we have now an equilibrium calculated.  Do a cleanup of the structure
! for phases with several compsets the call below shifts the stable one
! to the lowest compset number unless the default constitution fits another
! For example to ensure fcc-carbonitrides always is the same compset.
    ij=1
    call todo_after_found_equilibrium(ij,ceq)
    if(gx%bmperr.ne.0) goto 1000
1000 continue
   call system_clock(count=endoftime)
   call cpu_time(finish2)
   if(gx%bmperr.ne.0) then
!      write(*,1005)gx%bmperr
1005  format('Failed calculation, error code: ',i5)
   else
      write(*,1010)finish2-starting,endoftime-starttid
1010  format('Equilibrium calculation ',1pe12.4,' s and ',i7,' clockcycles')
   endif
    return
  end subroutine calceq2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine meq_phaseset(meqrec,ceq)
! this subroutine can change the set of stable phase and their amounts
! and constitutions until equilibrium is found for the current conditions.
! nv is the number of stable phases initially, iphl, icsl and aphl are the
! phase number, compset number and amount.
    implicit none
    TYPE(meq_setup) :: meqrec
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! this is the initial set of stable phases
    character ch1*1
! should one use meqrec as pointer here???
    TYPE(meq_data), dimension(:), allocatable :: phr
    integer ok,iadd,iph,ics,irem,jj,jph,kk,lastchange,lokph,lokcs,minadd
    integer kph,minrem,mph,nip,ierr,nochange,zap
    double precision, parameter :: ylow=1.0D-3,addedphase_amount=1.0D-1
    double precision xxx
! number of iterations without adding or removing a phase
    minadd=meqrec%nrel
!    minrem=2
!    minrem=0
    minrem=4
! minum number if iterations between any change of stable phase set
!    nochange=2
    nochange=2
    lastchange=0
!
!    write(*,*)'entering meq_phaseset: '
    meqrec%dormlink=0
! nphase is set to total number of phases (phase+compset) to be calculated
! >>> parallellization ALERT, nphase may change when composition sets created
    call sumofphcs(meqrec%nphase,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'number of phases+compsets: ',meqrec%nphase,meqrec%nv
    allocate(phr(meqrec%nphase))
! set default T and P as not variable (change below if no condition on T or P)
!    meqrec%tpindep=.FALSE.
! order the inital set of stable phases in ascending order
! VERY CLUMSY SORTING
15  continue
    ok=0
!    write(*,16)meqrec%nv,meqrec%nphase,size(meqrec%iphl)
!16  format('sort: ',10i3)
    do iph=2,meqrec%nv
       if(meqrec%iphl(iph-1).gt.meqrec%iphl(iph)) then
          ok=1
          kk=meqrec%iphl(iph-1)
          meqrec%iphl(iph-1)=meqrec%iphl(iph)
          meqrec%iphl(iph)=kk
          kk=meqrec%icsl(iph-1)
          meqrec%icsl(iph-1)=meqrec%icsl(iph)
          meqrec%icsl(iph)=kk
          xxx=meqrec%aphl(iph-1)
          meqrec%aphl(iph-1)=meqrec%aphl(iph)
          meqrec%aphl(iph)=xxx
       endif
    enddo
    if(ok.ne.0) goto 15
17  continue
    ok=0
    do iph=2,meqrec%nv
       if(meqrec%iphl(iph-1).eq.meqrec%iphl(iph)) then
          if(meqrec%icsl(iph-1).gt.meqrec%icsl(iph)) then
             kk=meqrec%icsl(iph-1)
             meqrec%icsl(iph-1)=meqrec%icsl(iph)
             meqrec%icsl(iph)=kk
             xxx=meqrec%aphl(iph-1)
             meqrec%aphl(iph-1)=meqrec%aphl(iph)
             meqrec%aphl(iph)=xxx
             ok=1
          endif
       endif
    enddo
    if(ok.ne.0) goto 17
!    write(*,18)meqrec%nv,(meqrec%iphl(jj),meqrec%icsl(jj),&
!         meqrec%aphl(jj),jj=1,meqrec%nv)
!18  format('MEQ start:  ',i2,' phases: ',4(i3,i2,1PE10.2))
!    stop 'order'
    mph=0
    nip=1
    meqrec%nstph=0
    do iph=1,noph()
       do ics=1,noofcs(iph)
! ignore hidden and suspended phases (also ignored above in sumofphcs)
! entered, fixed and dormat has values 1, 2 and 3, suspended 4, hidden 5
          zap=test_phase_status(iph,ics,xxx,ceq)
          if(zap.le.3) then
             mph=mph+1
             phr(mph)%iph=iph
             phr(mph)%ics=ics
! compare with these the first time a phase wants to be added or removed
! if zero it means phase can be added/removed at iteration minadd/minrem
             phr(mph)%itadd=0
             phr(mph)%itrem=0
! initiate indicator for phases with fix composition, set to 1 later if so
             phr(mph)%xdone=0
! save status
             phr(mph)%phasestatus=zap
! set link to calculated values of G etc.
             call get_phase_compset(iph,ics,lokph,lokcs)
             phr(mph)%curd=>ceq%phase_varres(lokcs)
             if(iph.eq.meqrec%iphl(nip) .and. ics.eq.meqrec%icsl(nip)) then
! this phase is part of the initial stable set, increment nstph
                meqrec%nstph=meqrec%nstph+1
                meqrec%stphl(meqrec%nstph)=mph
                phr(mph)%stable=1
                phr(mph)%curd%amfu=meqrec%aphl(meqrec%nstph)
!                phr(mph)%curd%amount(1)=meqrec%aphl(meqrec%nstph)
! set "previous values"
                phr(mph)%prevam=meqrec%aphl(meqrec%nstph)
                phr(mph)%prevdg=zero
                nip=nip+1
             else
! unstable phase
                phr(mph)%stable=0
                phr(mph)%prevam=zero
                phr(mph)%prevdg=-one
                phr(mph)%curd%amfu=zero
!                phr(mph)%curd%amount(1)=zero
             endif
! mark that no data arrays allocated for this phase
             phr(mph)%idim=0
! initiate link to another phase temporarily set dormant zero
             phr(mph)%dormlink=0
          else
! make sure stable bit is cleared in phases not included in calculation
! maybe the whole status word should be zeroed?
             call get_phase_compset(iph,ics,lokph,lokcs)
             ceq%phase_varres(lokcs)%status2=&
                  ibclr(ceq%phase_varres(lokcs)%status2,CSSTABLE)
          endif
       enddo
    enddo
! problem phases suspended are restored!!
!    write(*,*)'at start, nonsuspenden phases: ',mph
    meqrec%noofits=0
! return here after change of set of stable phase
! irem nonzero if phase irem should be removed
! iadd nonzero if phase iadd should be added
! meqrec%noofits is iteration counter
! nphase is dimension of phr (redundant)
! nstph is number of stable phases
! meqrec%stphl is indices in phr of stable phases
200 continue
    iadd=0
    irem=0
!    if(meqrec%noofits.gt.0) then
!       write(*,217)meqrec%noofits,meqrec%nstph,&
!            (phr(meqrec%stphl(jj))%iph,&
!            phr(meqrec%stphl(jj))%ics,&
!            phr(meqrec%stphl(jj))%curd%amfu,&
!            jj=1,meqrec%nstph)
!217    format('meq_phaseset:',i3,',',i2,' phases:',4(i3,i2,1pe10.2))
!    endif
    call meq_sameset(irem,iadd,meqrec,phr,ceq)
!    write(*,*)'Back from sameset ',irem,iadd,meqrec%noofits
    if(gx%bmperr.ne.0) goto 1000
    if(irem.gt.0 .or. iadd.gt.0) then
       if(meqrec%noofits-lastchange.lt.nochange) then
!          write(*,221)' *** Phase set change not allowed: ',&
!               meqrec%noofits,lastchange,nochange,irem,iadd
221       format(a,10i3)
          goto 200
       endif
    endif
    remove: if(irem.gt.0) then
! remove a phase ---------------------------
!       write(*,223)'Phase to be removed: ',phr(irem)%iph,phr(irem)%ics,&
!            phr(irem)%curd%amfu,meqrec%noofits
       if(meqrec%nstph.eq.1) then
          write(*,*)'Attempt to remove the only stable phase!!!'
!          gx%bmperr=7777; goto 1000
          goto 200
       endif
       if(meqrec%noofits-phr(irem)%itadd.lt.minrem) then
! if phase was just added do not remove before minrem iterations
          write(*,*)'Too soon to remove phase',&
               phr(irem)%iph,meqrec%noofits,phr(irem)%itadd
          goto 200
       endif
! shift phases after irem down in meqrec%stphl
! irem is index to phr(), meqrec%stphl(jph) is index to phr
       meqrec%nstph=meqrec%nstph-1
       do iph=1,meqrec%nstph
          jj=meqrec%stphl(iph)
          if(jj.ge.irem) then
             meqrec%stphl(iph)=meqrec%stphl(iph+1)
          endif
       enddo
! we must zero the last phase !!
       meqrec%stphl(meqrec%nstph+1)=0
!
       phr(irem)%itrem=meqrec%noofits
       phr(irem)%prevam=zero
       phr(irem)%stable=0
!       phr(irem)%curd%amount(1)=zero
       phr(irem)%curd%amfu=zero
       if(phr(irem)%xdone.eq.1) &
            write(*,*)'**** Removing phase with fix composition: ',&
            irem,phr(irem)%iph
       irem=0
       lastchange=meqrec%noofits
! one can remove and add a phase at the same time !!!
       if(iadd.eq.0) then
          goto 200
       endif
    endif remove
!------------------------------------------- 
    add: if(iadd.gt.0) then
! add a phase.  This can be tricky
! NOTE it must be added so meqrec%stphl in ascending order
       write(*,223)'Phase to be added:   ',phr(iadd)%iph,phr(iadd)%ics,&
            phr(iadd)%curd%dgm,meqrec%noofits
223    format(a,2x,2i4,1pe15.4,i7)
       if(meqrec%noofits-phr(iadd)%itrem.lt.minadd) then
! if phase was just removed do not add it before minadd iterations
!          write(*,224)'Too soon to add phase: ',phr(iadd)%iph,&
!               phr(iadd)%ics,meqrec%noofits,phr(iadd)%itrem
224       format(a,2i4,5i5)
          goto 200
       endif
       if(meqrec%nstph.eq.meqrec%maxsph) then
! No more phases allowed, we must see if  some other phase may be removed
          write(*,*)'Attempt to set too many phases stable'
!          write(*,225)'Attempt to set too many phases stable',iadd,&
!               (phr(meqrec%stphl(j))%iph),j=1,meqrec%nstph)
!225       format(a,i3,5x,10i3)
          gx%bmperr=4201; goto 1000
       endif
! the phase must be added in sequential order of phase and composition set no
       findplace: do jph=1,meqrec%nstph
          jj=meqrec%stphl(jph)
          if(phr(iadd)%iph.gt.phr(jj)%iph) then
             cycle
          endif
          if(phr(iadd)%iph.lt.phr(jj)%iph) then
             exit
          endif
! if same phase number compare composition set numbers
          if(phr(iadd)%iph.eq.phr(jj)%iph) then
             if(phr(iadd)%ics.gt.phr(jj)%ics) then
                cycle
             else
                exit
             endif
          endif
       enddo findplace
! one should come here at exit, iadd should be inserted before 
! meqrec%stphl(jph), jph can be nstph+1 if added phase should be the last
! otherwise shift previous phases one step up.
       do kph=meqrec%nstph,jph,-1
          meqrec%stphl(kph+1)=meqrec%stphl(kph)
       enddo
! phase added at jph, (note jph may be equal to nstph+1)
       meqrec%stphl(jph)=iadd
       meqrec%nstph=meqrec%nstph+1
       phr(iadd)%itadd=meqrec%noofits
       phr(iadd)%curd%dgm=zero
       lastchange=meqrec%noofits
! maybe some more variables should be set?
!       phr(iadd)%curd%amount(1)=addedphase_amount
       phr(iadd)%curd%amfu=addedphase_amount
       phr(iadd)%stable=1
       if(phr(iadd)%xdone.eq.1) &
            write(*,*)'**** Adding phase with fix composition:   ',&
            iadd,phr(iadd)%iph
       iadd=0
       goto 200
    endif add
!---------------------------------------------------
! found stable phase set or error
1000 continue
    if(gx%bmperr.eq.0) then
! equilibrium calculation converged, one should add check on stability
!
! >> calculate eigenvalues of phase matrix to check stability, 
! >> a negative eigenvalue means inside spinodal
! >> Note charge problems for metastable phases, phase must be neutral ...
!
!------------------------------------------------------------
! clear bits: no equilibrium calculated/ inconsistent conditions and result/
! equilibrium calculation failed
       ceq%status=ibclr(ceq%status,EQNOEQCAL)
       ceq%status=ibclr(ceq%status,EQINCON)
       ceq%status=ibclr(ceq%status,EQFAIL)
! set stable bit in stable phases and clear it in all others
       kk=1
       do jj=1,mph
          if(jj.eq.meqrec%stphl(kk)) then
             phr(jj)%curd%status2=ibset(phr(jj)%curd%status2,CSSTABLE)
! the stable phase list should be ordered in increasing phase number
             kk=kk+1
!             write(*,*)'Stable phase: ',jj,phr(jj)%iph,phr(jj)%ics
          else
             phr(jj)%curd%status2=ibclr(phr(jj)%curd%status2,CSSTABLE)
!             write(*,*)'Unstable phase: ',jj,phr(jj)%iph,phr(jj)%ics
          endif
       enddo
    else
! set some failure bits
       ceq%status=ibset(ceq%status,EQINCON)
       ceq%status=ibset(ceq%status,EQFAIL)
    endif
! restore phases set dormant
    jj=meqrec%dormlink
1200 continue
    if(jj.ne.0) then
!       write(*,*)'Restore from dormant: ',jj,phr(jj)%iph,phr(jj)%ics
       phr(jj)%curd%status2=ibclr(phr(jj)%curd%status2,CSSUS)
       phr(jj)%curd%status2=ibclr(phr(jj)%curd%status2,CSFIXDORM)
       jj=phr(jj)%dormlink
       goto 1200
    endif
! try to find problem with listed chemical potential    
! chempot(2) should be value with user defined reference state, has
! to be implemented ...
    do jj=1,meqrec%nrel
       ceq%complist(jj)%chempot(2)=ceq%complist(jj)%chempot(1)
!       write(*,*)'Meq_phaseset chempot: ',ceq%complist(jj)%chempot(1)
    enddo
    deallocate(phr)
! >>>> here one can allow new composition set in parallelization
    return
  end subroutine meq_phaseset

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine meq_sameset(irem,iadd,meqrec,phr,ceq)
! iterate until phase set change, converged or error (incl too many its)
    implicit none
    integer irem,iadd
    TYPE(meq_setup) :: meqrec
    TYPE(meq_data), dimension(*), target :: phr
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer increase,ioff,im,jd,ik,jk,jj,jph,ie,ierr,jmaxy
    integer k1,k2,kk,kkz,level3,mph,negam,ncol,negamph,nj,nk,nl,nrow
    integer nz1,nz2
    TYPE(gtp_condition), pointer :: condition,lastcond
    TYPE(meq_data), pointer :: pmi
    character ch1*1
    double precision, dimension(maxel) :: ccm,xknown,xxmol,wmass,xfrac,pccm
    double precision, dimension(maxel) :: sccm
    double precision, dimension(5) :: qq
    double precision, dimension(2) :: tpvalz
    double precision, dimension(maxel+2) :: ycormax
    double precision phfrac(maxel),xpfrac(maxel,maxel)
    double precision, dimension(:,:), allocatable :: smat
    double precision, dimension(:), allocatable :: svar
! these arrays should maybe be allocated ....
    double precision, dimension(maxconst) :: ycorr,yarr,addont,addonp
    integer converged,jz,cmix(10),cmode,ctype
    double precision molescomp,amsum,antot,chargefact,chargerr,clhs,cvalue
    double precision dgm,summ,dgmmax,gsurf,pclhs,pfsum,phf,phscale,phs
    double precision prevmaxycorr,pv,sclhs,signerr,summol,tmass,tmol
    double precision xall,xmax,xxx,ycormax2,yprev,ys,ysmm,ysmt,yss,yst
    double precision, parameter :: ylow=1.0D-3,ymin=1.0D-12,ymingas=1.0D-30
    double precision mamu(maxel),mag,mat,map,totam,pham,zval,xval
    double precision zmat(2*maxel,2*maxel),zrhs(2*maxel),sum
    double precision, dimension(:), allocatable :: xcol,zcol
    double precision, dimension(:), allocatable :: cit
    double precision deltat,deltap,deltaam
    integer je,ke,iz,iz0,iz1,stvix,stvnorm,tcol,pcol,sel,sph,scs
    integer notf,kph,kjj,dncol,iy,jy
    logical doubt,removeme
!    write(*,*)'entering meq_sameset',meqrec%nphase
    ycormax=one
!    doubt=.TRUE.
    doubt=.FALSE.
! force use of new way to calculate system matrix
!    meqrec%typesofcond=2
! dimension matrix for conditions, components+stable phases
    nz1=meqrec%nrel-meqrec%nfixmu+meqrec%nstph-meqrec%nfixph
    if(meqrec%tpindep(1)) nz1=nz1+1
    if(meqrec%tpindep(2)) nz1=nz1+1
!    write(*,11)meqrec%nrel,meqrec%nfixmu,meqrec%nstph,meqrec%nfixph,&
!         meqrec%tpindep,nz1
!11  format('System matrix dimensions: ',4i7,2l2,i5)
    nz2=nz1+1
    allocate(smat(nz1,nz2))
    allocate(svar(nz1))
! check if constituent fraction correction in stable phases increases
! for each iteration.  Needed for the Re-V case ....
    prevmaxycorr=zero
    increase=0
    level3=0
!-------------------------------------------------------------
! return here until converged or phase set change
100 continue
    meqrec%noofits=meqrec%noofits+1
101 format(a)
!    write(*,*)'Iteration: ',meqrec%noofits,' ----------------------------- '
!    write(*,199)(meqrec%stphl(jz),jz=1,meqrec%nstph)
199 format('Index: ',10i3)
!    write(*,102)(phr(meqrec%stphl(jz))%iph,phr(meqrec%stphl(jz))%ics,&
!         jz=1,meqrec%nstph)
102 format('stable phases: ',10(i3,i2))
!    write(*,103)(phr(meqrec%stphl(jz))%curd%amfu,jz=1,meqrec%nstph)
103 format('amounts: ',6(1pe12.4))
    if(meqrec%noofits.gt.ceq%maxiter) goto 1200
    converged=0
! loop for all phases and composition sets, loop over phr
    do mph=1,meqrec%nphase
       pmi=>phr(mph)
! this routine calculates the phase matrix and inverts it.
! it also calculates the amounts of moles of components in the phase
!       call meq_onephase(mph,pmi,ceq)
       call meq_onephase(meqrec%nrel,pmi,ceq)
       if(gx%bmperr.ne.0) then
          if(pmi%stable.eq.0) then
! if this happends for an unstable phase just continue but ensure it will
! not be stable (in a very crude way)
!             write(*,*)'Matrix inversion error for unstable phase',pmi%iph
             pmi%curd%gval(1,1)=one
             gx%bmperr=0
          else
! Inversion error for stable phase is fatal, error code already set
!             write(*,*)'Matrix inversion error for stable phase',pmi%iph
             goto 1000
          endif
       endif
!       write(*,*)'Back from meq_onephase'
!       do i=1,pmi%ncc
!          write(*,107)'pmi: ',(pmi%invmat(i,j),j=1,pmi%ncc)
107       format(a,6(1pe12.3))
!       enddo
! end of pmi% scope
    enddo
!=======================================================================
! step 2: calculation of system matrix
! Solve for chemical potentials and conditions using all stable phases
! The SYSTEM MATRIX (smat) has one row for each stable phase and
! one row for each component representing a condition
! (If a fix phase condition or chem.pot. condition slightly different??)
!    if(meqrec%typesofcond.ne.1) goto 300
!----------------------------------------
! for testing the new code always jump to 300
!    goto 300
!==============================================================
!
! removed old code
!
!==============================================================
! not all massbalance conditions, meqrec%typesofcond=2, generate system matrix
300 continue
!    write(*,301)'Calculating general system matrix',meqrec%nfixmu,&
!         meqrec%nfixph,meqrec%tpindep,meqrec%noofits
301 format(/a,2i2,2l2,i5)
!-------------------------------------------------------------------
! Formulating the system equation in general:
! Variables (one column per variable):
! - The chemical potentials of the components:     MEQREC%NREL
!   minus the number of fixed chemical potentials: -MEQREC%NFIXMU
! - The variation in T if not fixed                +1
! - The variation in P if not fixed                +1
! - The variation of the amounts of stable phases: MEQREC%NSTPH
!   minus those that have fixed amount:            -MEQREC%NFIXPH
!
! The variables will be ordered: MU, DeltaT, DeltaP, Delta Phase amounts
! this is important for the order of columns in the system matrix
!
! Equations (one row per equation):
! If T or P are variable extra columns and terms are needed
! - The expression for the Gibbs energy for each stable phase, M_A mu_A = G
!   if a fixed chemical potentials incorporated that incorporated
!   if T or P variable an extra term for these
! - The user defined conditions like:
!   - Amount of components, N(A)= or B(A)=
!   - The total amount of moles, N=, or mass, B=
!   - Overall mole fractions, x(A)=, or mass fractions, w(A)=
!   - Phase specific mole or mass fractions, x(FCC,C)= or w(LIQUID,B)=
!   - The volume V=; enthalpy H= etc., with phase spec and normallizing
!   - relations between state variables x(C14,Fe)-x(liq,Fe)= 0 etc.
!
! The equations will always have the G expressions first.  The other will 
! be random (or in order of the user entered them)
!
! There must be as many equations as there are variables and the construction
! of the equations can be rather complex.  
! At present only a limited set has been implemented.
!
!-------------------------------------------------------------------
! zero all values in system matrix, dimension (nz1)x(nz1)
    smat=zero
    tcol=0
    pcol=0
! step 2.1 the Gibbs energies for the stable phases (incl fixed)
    allstableph: do jph=1,meqrec%nstph
       jj=meqrec%stphl(jph)
!       write(*,*)'Stable phase: ',phr(jj)%iph,phr(jj)%ics
! column nz2 is the right hand side of the equation, to molar G
       smat(jph,nz2)=phr(jj)%curd%gval(1,1)
!       write(*,*)'Gm: ',jph,phr(jj)%curd%gval(1,1),ceq%tpval(1)
! one column with amount of component A for each variable chemical potential
! components with fixed chemical potential are automatically skipped
       ncol=1
       xxx=zero
       gloop: do je=1,meqrec%nrel
          do ie=1,meqrec%nfixmu
             if(meqrec%mufixel(ie).eq.je) then
! meqrec%mufixel(ie) is the component number with fix mu, better use curmu ?
! UNIFISHED: reference state must be handelled (may depend on T) ??
                xxx=smat(jph,nz2)
                smat(jph,nz2)=smat(jph,nz2)-&
                     phr(jj)%xmol(je)*meqrec%mufixval(ie)
!                write(*,312)'fix mu: ',jj,je,ie,xxx,smat(jph,nz2),&
!                     phr(jj)%xmol(je),meqrec%mufixval(ie)
312             format(a,3i3,6(1pe12.4))
                cycle gloop
             endif
          enddo
          smat(jph,ncol)=phr(jj)%xmol(je)
          ncol=ncol+1
       enddo gloop
       
       if(meqrec%tpindep(1)) then
! column for variable T, value is -dG/dT ??
          smat(jph,ncol)=-phr(jj)%curd%gval(2,1)
!          smat(jph,ncol)=+phr(jj)%curd%gval(2,1)
          tcol=ncol
!          write(*,353)'dG/dT: ',jph,tcol,phr(jj)%curd%gval(2,1)
353       format(a,2i3,1pe12.4)
          ncol=ncol+1
       endif
       if(meqrec%tpindep(2)) then
! column for variable P, value is +dG/dP ??
          smat(jph,ncol)=-phr(jj)%curd%gval(3,1)
!          smat(jph,ncol)=+phr(jj)%curd%gval(3,1)
          pcol=ncol
          write(*,353)'dG/dP: ',jph,pcol,phr(jj)%curd%gval(2,1)
          ncol=ncol+1
       endif
    enddo allstableph
! we have generated meqrec%nstph rows with ncol columns and rhs in column nz2
! The columns for delta_phase-amounts should be zero
    dncol=ncol-1
!    do iz=1,dncol
!       write(*,228)'smat 1: ',(smat(iz,jz),jz=1,nz2)
!    enddo
228    format(a,6(1pe12.4))
    nrow=meqrec%nstph
!-------------------------------------------------------------------
! step 2.2 equations due to user conditions
    lastcond=>ceq%lastcondition
    condition=>lastcond
    jz=0
350 continue
    cmode=0
    cmix=0
    condition=>condition%next
    jz=jz+1
    call apply_condition_value(condition,cmode,ctype,cvalue,cmix,ceq)
    if(gx%bmperr.ne.0) goto 1000
! Only cmix(1)=5 is interesting here
    if(cmix(1).ne.5) then
! loop if not the last condition
!       write(*,*)'Taking next condition: ',cmix(1)
       if(.not.associated(condition,lastcond)) goto 350
       goto 380
    endif
! do something with the condition ... it can be N=1, x(A)=.1, VM(GAS)=1e-6 etc.
! THE MASTER VERSION OF THIS TABLE in PMOD25C.F90
! symb cmix(2) indices                   irrelevant Property
! U       10   (phase#set)                    6     Internal energy (J)
! UM      11    "                             6     per mole components
! UW      12    "                             6     per kg
! UV      13    "                             6     per m3
! UF      14    "                             6     per formula unit
! S       2x    "                             7     entropy
! V       3x    "                             8     volume
! H       4x    "                             9     enthalpy
! A       5x    "                            10     Helmholtz energy
! G       6x    "                            11     Gibbs energy
! NP      7x    "                            12     moles of phase
! BP      8x    "                            13     mass of moles
! DG      9x    "                            15     Driving force
! Q       19x   "                            14     Internal stability
! N       11x  (component/phase#set,component) 16  moles of components
! X       111   "                            17     mole fraction of components
! B       12x   "                            18     mass of components
! W       122   "                            19     mass fraction of components
! Y       13    phase#set,constituent#subl   20     constituent fraction
!----- model variables <<<< these now treated differently
    stvix=cmix(2)/10
    stvnorm=mod(cmix(2),10)
    select case(stvix)
    case default
       write(*,*)'not a condition:',stvix,stvnorm,cmix(1),cmix(2),cmix(3)
       gx%bmperr=9988; goto 1000
    case(1:6) 
! stvix=1..6: U, S, V, H, A, G conditions not implemented
       write(*,*)'Not implemented yet: ',stvix,stvnorm
       gx%bmperr=9988; goto 1000
!------------------------------------------------------------------
    case(7) ! NP
! Amount of phase in moles, use fix phase instead
       write(*,*)'Not implemented yet, use set phase fix: ',stvix,stvnorm
       gx%bmperr=9988; goto 1000
       nrow=nrow+1
       if(nrow.gt.nz1) stop 'too many equations 7A'
!------------------------------------------------------------------
    case(8) ! BP
! Amount of phase in mass, use fix phase instead
       write(*,*)'Not implemented yet, use set phase fix: ',stvix,stvnorm
       gx%bmperr=9988; goto 1000
       nrow=nrow+1
       if(nrow.gt.nz1) stop 'too many equations 8A'
!------------------------------------------------------------------
! 9 and 10 (DG and Q) not allowed as conditions
!------------------------------------------------------------------
    case(11) ! N or X with or without indices and normalization
       if(stvnorm.eq.0) then
          if(cmix(3).eq.0) then
! condition is N=fix
             sel=0; sph=0
          elseif(cmix(4).eq.0) then
! condition is N(A)=fix
             sel=cmix(3); sph=0
          else
! condition is N(phase#set,A)=fix;  how to handle if phase#set not stable?
             write(*,*)'Condition N(phase#set,A)=fix not allowed'
             gx%bmperr=9898; goto 1000
             sel=cmix(5); sph=cmix(3); scs=cmix(4)
          endif
! Formulate equation for total amount N:
! rhs:  N-N+\sum_alpha N^a + \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j
! \sum_B \sum_alpha N^a \sum_i \sum_j dM^_A/dy_i dM^a_B/dy_j*z^a_ij  *mu(B)
!        \sum_alpha N^a \sum_i d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j      *deltaT
!        \sum_alpha N^a \sum_i d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j      *deltaP
!        \sum_A M^a_A                                    *deltaN^a
          allocate(xcol(nz2))
          xcol=zero
          totam=zero
          nallph: do jph=1,meqrec%nstph
! sum over all stable phases
             jj=meqrec%stphl(jph)
             pmi=>phr(jj)
! moles formulat unit of phase
!             pham=pmi%curd%amount(1)
             pham=pmi%curd%amfu
!             write(*,*)'amfu: ',pham
             nallel: do ie=1,meqrec%nrel
! if sel=/=0 then skip all components except sel
                if(sel.gt.0 .and. ie.ne.sel) cycle nallel
! multiply terms with the inverse phase matrix
                call calc_dyterms1(meqrec%nrel,ie,meqrec%tpindep,&
                     mamu,mag,mat,map,pmi)
                if(gx%bmperr.ne.0) goto 1000
! the call above calculates (A is "ie", z_ij is the inverted phase matrix): 
! mamu_A(B=1..nrel) = \sum_i \sum_j dM^a_A/dy_i dM^a_B/dy_j z^a_ij
! mag_A             = \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j
! mat_A             = \sum_i \sum_j d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j
! map_A             = \sum_i \sum_j d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j
                ncol=1
! calculate a term for each column to be multiplied with chemical potential
! if the potential is fixed add the term to the rhs
                nloop1: do je=1,meqrec%nrel
                   do ke=1,meqrec%nfixmu
                      if(meqrec%mufixel(ke).eq.je) then
! components with fix chemical potential added to rhs, do not increment ncol!!!
                         xcol(nz2)=xcol(nz2)+pham*mamu(je)*meqrec%mufixval(ke)
                         cycle nloop1
                      endif
                   enddo
! mamu(B) = \sum_i \sum_j \sum_A dM^a_B/dy_i dM^a_A z^a_ij
                   xcol(ncol)=xcol(ncol)-pham*mamu(je)
                   ncol=ncol+1
                enddo nloop1
! If T or P are variable
                if(tcol.gt.0) then
                   xxx=xcol(tcol)
                   xcol(tcol)=xcol(tcol)+pham*mat
!                   write(*,363)'d2G/dTdy: ',nrow-1,ie,tcol,&
!                        xxx,xcol(tcol),pham,mat
                endif
                if(pcol.gt.0) then
                   xxx=xcol(pcol)
                   xcol(pcol)=xcol(pcol)+pham*map
!                   write(*,363)'d2G/dPdy: ',nrow-1,ie,pcol,&
!                        xxx,xcol(pcol),pham,mat
                endif
! last columns on lhs are amounts of element ie for all stable non-fix phases
                notf=1
                do kph=1,meqrec%nstph
                   kjj=meqrec%stphl(kph)
                   if(phr(kjj)%phasestatus.eq.1) then
                      if(sel.gt.0 .and. sel.eq.ie) then
                         xcol(dncol+notf)=phr(kjj)%xmol(ie)
                      else
                         xcol(dncol+notf)=xcol(dncol+notf)+phr(kjj)%xmol(ie)
                      endif
                      notf=notf+1
!                   elseif(phr(jj)%phasestatus.eq.2) then
! nothing to do in else link, stable fix phases have no column
!                      write(*,*)'Fixed phase 1: ',kjj
!                   else
! nothing to do in else link, stable fix phases have no column
                   endif
                enddo
! right hand side (rhs) contribution is
! - NP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij
                xcol(nz2)=xcol(nz2)-pham*mag
             enddo nallel
! sum of moles in phase will be multiplied with delta-phase_amount
             if(sel.gt.0) then
                totam=totam+pham*pmi%xmol(sel)
!                totam=totam+pham*pmi%nmol(sel)
!                write(*,363)'totam 1: ',jj,sel,0,pham,&
!                     pmi%xmol(sel),pmi%nmol(sel)
             else
!                write(*,363)'totam 2: ',jj,0,0,pham,pmi%sumxmol,pmi%sumnmol
!                totam=totam+pham*pmi%sumnmol
                totam=totam+pham*pmi%sumxmol
             endif
!             write(*,363)'363: ',jph,dncol,ncol,pham,totam,&
!                  (mamu(je),je=1,2),(xcol(k1),k1=1,nz2)
363          format(a,3i3,6(1pe12.4))
          enddo nallph
!
! in xcol are values summed over all phases and components
! copy summed columns to smat nrow
          nrow=nrow+1
          if(nrow.gt.nz1) then
             write(*,*)'too many equations 11A',nrow
             stop
          endif
          do ncol=1,nz2
             smat(nrow,ncol)=xcol(ncol)
          enddo
! add N^prescribed - N^current to rhs (right hand side)
          xxx=smat(nrow,nz2)
          smat(nrow,nz2)=smat(nrow,nz2)-cvalue+totam
!          write(*,363)'RHS: ',nrow,nz2,0,smat(nrow,nz2),xxx,cvalue,totam,&
!               cvalue-totam
          deallocate(xcol)
          if(abs(totam-cvalue).gt.ceq%xconv) then
!             write(*,266)'Unconverged condition N(A): ',sel,&
!                  cvalue,zval
!266       format(a,i3,2(1pe15.7))
             if(converged.lt.5) converged=5
          endif
!          if(sel.eq.0) then
!             write(*,363)'Condition N=fix   ',0,0,0,cvalue,totam
!          else
!             write(*,363)'Condition N(a)=fix',sel,0,0,cvalue,totam
!          endif
       elseif(stvnorm.gt.1) then
! only normallizing of N with respect to amount of moles (M) is allowed
          write(*,*)'N can only be normalled with M',stvix,stvnorm,cmix(2)
          gx%bmperr=9988; goto 1000
       else
!-------------------------------
! N=fix and N(A)=fix treated above as they have a "simple" summation, 
! Now handle all other cases like x(A)=fix, x(phase#set,A)=fix
! We must sum over all phases and constituents for the normallizing factor
! definition: X(A)=N(A)/N; 
! derivative: dX(A)=dN(A)/N - N(A)/N**2 *dN
! sum dN(A) and dN at the same time and multiply the sums with 1/N 
! and -N(A)/N**2 in the end.
          if(cmix(3).eq.0) then
             write(*,*)'Condition NM=fix is illegal'
             gx%bmperr=9898; goto 1000
          elseif(cmix(4).eq.0) then
! condition is x(A)=fix
             sel=cmix(3); sph=0
          else
! condition is x(phase#set,A)=fix
!             write(*,*)'Condition x(phase#set,A)=fix not yet allowed'
!             gx%bmperr=9898; goto 1000
             sel=cmix(5); sph=cmix(3); scs=cmix(4)
          endif
! two summations, zcol sums the term dN(A); xcol sums dN (as above)
          allocate(xcol(nz2))
          allocate(zcol(nz2))
          xcol=zero
          zcol=zero
          totam=zero
          zval=zero
          xval=zero
!          xallph: do jph=1,meqrec%nstph
!             jj=meqrec%stphl(jph)
! sum over all phases to handle conditions like x(phase#set,A)=fix
! as the phase#set may not be stable
          xallph: do jj=1,meqrec%nphase
             if(sph.eq.0) then
! skip this phase if not stable and condition not on a specific phase (sph)
                if(phr(jj)%stable.eq.0) cycle xallph
             else
! condition on specific phase, skip this phase if not the right one
                if(phr(jj)%iph.ne.sph .or. phr(jj)%ics.ne.scs) cycle xallph
             endif
             pmi=>phr(jj)
! moles formulat unit of phase
!             pham=pmi%curd%amount(1)
             pham=pmi%curd%amfu
             xallel: do ie=1,meqrec%nrel
! we cannot skip summation over all element as that is needed for normallizing
! calculate a term for each column to be multiplied with chemical potential
! we must sum xcol for all elemenets and add to zcol for element sel
! if sel=/=0 then we sum also zcol(sel) for all phases
                call calc_dyterms1(meqrec%nrel,ie,meqrec%tpindep,&
                     mamu,mag,mat,map,pmi)
                if(gx%bmperr.ne.0) goto 1000
                ncol=1
                nloop2: do je=1,meqrec%nrel
! Calculate one column for each component to be multiplied with chem.pot.
! components with fix chemical potential added to rhs, do not increment ncol!!!
                   do ke=1,meqrec%nfixmu
                      if(meqrec%mufixel(ke).eq.je) then
!                         xcol(nz2)=xcol(nz2)-pham*mamu(je)*meqrec%mufixval(ke)
                         xcol(nz2)=xcol(nz2)+pham*mamu(je)*meqrec%mufixval(ke)
                         if(sel.eq.je) then
                            zcol(nz2)=zcol(nz2)-&
                                 pham*mamu(je)*meqrec%mufixval(ke)
                         endif
                         cycle nloop2
                      endif
                   enddo
! mamu(B) = \sum_i \sum_j dM^a_B/dy_i dM^a_A z^a_ij
! sum over all elements for normallizing
                   xcol(ncol)=xcol(ncol)-pham*mamu(je)
                   if(sel.eq.ie) then
! if this is the specified element sum to zcol
                      zcol(ncol)=zcol(ncol)-pham*mamu(je)
                   endif
                   ncol=ncol+1
                enddo nloop2
! If T or P are variable
                if(tcol.gt.0) then
                   xcol(tcol)=xcol(tcol)+pham*mat
                   if(sel.eq.ie) then
                      zcol(tcol)=zcol(tcol)+pham*mat
                   endif
                endif
                if(pcol.gt.0) then
                   xcol(pcol)=xcol(pcol)+pham*map
                   if(sel.eq.ie) then
                      zcol(pcol)=zcol(pcol)+pham*map
                   endif
                endif
! last columns are amounts of element ie for all stable non-fix phase,
                notf=1
                do kph=1,meqrec%nstph
                   kjj=meqrec%stphl(kph)
                   if(phr(kjj)%phasestatus.eq.1) then
! phasestatus=1 for phases with variable amount, sum over all components
                      xcol(dncol+notf)=xcol(dncol+notf)+pham*phr(kjj)%xmol(ie)
!                      write(*,363)'xcol: ',dncol+notf,kjj,ie,&
!                           xcol(dncol+notf),phr(kjj)%xmol(ie),&
!                           phr(kjj)%curd%amfu
                      if(ie.eq.sel) then
!                         zcol(dncol+notf)=phr(kjj)%xmol(ie)
                         zcol(dncol+notf)=zcol(dncol+notf)+&
                              pham*phr(kjj)%xmol(ie)
                      endif
                      notf=notf+1
!                   elseif(phr(jj)%phasestatus.eq.2) then
!                      write(*,*)'Fixed phase 2: ',kjj
                   endif
                enddo
! right hand side (rhs) contribution is
! - NP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij 
                xcol(nz2)=xcol(nz2)-pham*mag
!                xcol(nz2)=xcol(nz2)+pham*mag
! the rhs sum the amount of component ie ... hm, different sign??
                if(sel.eq.ie) then
                   zcol(nz2)=zcol(nz2)-pham*mag
!                   zcol(nz2)=zcol(nz2)+pham*mag
                endif
             enddo xallel
!             totam=totam+pham*pmi%sumnmol
             totam=totam+pham*pmi%sumxmol
! if sph=/=0 next line must be changed
             zval=zval+pham*pmi%xmol(sel)
          enddo xallph
! in xcol is dN and in zcol dN(A) summed over all phases and components
! calculate the normallized values now
! xmat=dN(A)/N - N(A)*dN/N**2
          nrow=nrow+1
          if(nrow.gt.nz1) stop 'too many equations 11B'
!          write(*,228)'zcol: ',(zcol(jz),jz=1,nz2)
!          write(*,228)'xcol: ',(xcol(jz),jz=1,nz2)
!          write(*,228)'totam mm: ',totam,smat(nrow,nz2),cvalue,zval
! this was the trick!!
          zval=zval/totam
          do ncol=1,nz2
             smat(nrow,ncol)=(zcol(ncol)-xcol(ncol)*zval)/totam
          enddo
! add N^prescribed - N^current to rhs (right hand side)
          smat(nrow,nz2)=smat(nrow,nz2)-cvalue+zval
          deallocate(xcol)
          deallocate(zcol)
! check on convergence
          if(abs(zval-cvalue).gt.ceq%xconv) then
!             write(*,266)'Unconverged condition x(A): ',sel,&
!                  cvalue,zval
!266       format(a,i3,2(1pe15.7))
             if(converged.lt.5) converged=5
          endif
!          if(sph.eq.0) then
!             write(*,363)'Condition x(A)=fix',sel,0,0,cvalue,zval
!          else
!             write(*,363)'Condition x(phase#set,A)=fix',sph,sel,0,cvalue,zval
!          endif
       endif
!
!------------------------------------------------------------------
    case(12) ! B or W
! Amount of component in mass, can have indices and normallization
! code copied from the case(11) for N and X and modified
       if(stvnorm.eq.0) then
          if(cmix(3).eq.0) then
! condition is B=fix
             sel=0; sph=0
          elseif(cmix(4).eq.0) then
! condition is B(A)=fix
             sel=cmix(3); sph=0
          else
! condition is B(phase#set,A)=fix;  how to handle if phase#set not stable?
             write(*,*)'Condition B(phase#set,A)=fix not allowed'
             gx%bmperr=9898; goto 1000
             sel=cmix(5); sph=cmix(3); scs=cmix(4)
          endif
! Formulate equation for total amount B: each M_A multiplied with mass_A
! rhs:  B-B+\sum_alpha N^a + \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j 
! \sum_B \sum_alpha N^a \sum_i \sum_j dM^_A/dy_i dM^a_B/dy_j*z^a_ij  *mu(B)
!        \sum_alpha N^a \sum_i d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j      *deltaT
!        \sum_alpha N^a \sum_i d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j      *deltaP
!        \sum_A M^a_A                                    *deltaN^a
          allocate(xcol(nz2))
          xcol=zero
          totam=zero
          zval=zero
          ballph: do jph=1,meqrec%nstph
! sum over all stable phases
             jj=meqrec%stphl(jph)
             pmi=>phr(jj)
! amount of phase, amfu is moles formula units
! abnorm(1) is amount real atoms, abnorm(2) is mass
!             pham=pmi%curd%amount(1)
             pham=pmi%curd%amfu
             ballel: do ie=1,meqrec%nrel
! if sel=/=0 then skip all components except sel
                if(sel.gt.0 .and. ie.ne.sel) cycle
! multiply terms with the inverse phase matrix
                call calc_dyterms1(meqrec%nrel,ie,meqrec%tpindep,&
                     mamu,mag,mat,map,pmi)
                if(gx%bmperr.ne.0) goto 1000
! the call above calculates (A is "ie", z_ij is the inverted phase matrix): 
! mamu_A(B=1..nrel) = \sum_i \sum_j dM^a_A/dy_i dM^a_B/dy_j z^a_ij
! mag_A             = \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j
! mat_A             = \sum_i \sum_j d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j
! map_A             = \sum_i \sum_j d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j
                ncol=1
! calculate a term for each column to be multiplied with chemical potential
! if the potential is fixed add the term to the rhs
                bloop1: do je=1,meqrec%nrel
                   do ke=1,meqrec%nfixmu
                      if(meqrec%mufixel(ke).eq.je) then
! components with fix chemical potential added to rhs, do not increment ncol!!!
! NOTE: mamu includes summation of two components, multiply with two masses!!!
                         xcol(nz2)=xcol(nz2)+&
                              pham*mamu(je)*meqrec%mufixval(ke)*mass_of(ie,ceq)
                         cycle bloop1
                      endif
                   enddo
! mamu(B) = \sum_i \sum_j \sum_A dM^a_B/dy_i dM^a_A z^a_ij mass_A mass_B
                   xcol(ncol)=xcol(ncol)-pham*mamu(je)*mass_of(ie,ceq)
                   write(*,363)'xcola: ',ncol,ie,je,xcol(ncol),pham,mamu(je),&
                        mass_of(ie,ceq)
                   ncol=ncol+1
                enddo bloop1
! If T or P are variable
                if(tcol.gt.0) then
                   xxx=xcol(tcol)
                   xcol(tcol)=xcol(tcol)+pham*mat*mass_of(ie,ceq)
!                   write(*,363)'d2G/dTdy: ',nrow-1,ie,tcol,&
!                        xxx,xcol(tcol),pham,mat
                endif
                if(pcol.gt.0) then
                   xxx=xcol(pcol)
                   xcol(pcol)=xcol(pcol)+pham*map*mass_of(ie,ceq)
!                   write(*,363)'d2G/dPdy: ',nrow-1,ie,pcol,&
!                        xxx,xcol(pcol),pham,mat
                endif
! last columns are amounts of element ie for all stable non-fix phases
                notf=1
                do kph=1,meqrec%nstph
                   kjj=meqrec%stphl(kph)
                   if(phr(kjj)%phasestatus.eq.1) then
                      zval=zval+phr(kjj)%xmol(ie)*mass_of(ie,ceq)
                      if(sel.gt.0 .and. sel.eq.ie) then
                         xcol(dncol+notf)=&
                              phr(kjj)%xmol(ie)*mass_of(ie,ceq)
                      else
                         xcol(dncol+notf)=xcol(dncol+notf)+&
                              pham*phr(kjj)%xmol(ie)*mass_of(ie,ceq)
                      endif
                      write(*,363)'xcolb: ',dncol+notf,ie,kjj,&
                           xcol(dncol+notf),pham,phr(kjj)%xmol(ie),&
                           mass_of(ie,ceq)
                      notf=notf+1
!                   elseif(phr(jj)%phasestatus.eq.2) then
!                      write(*,*)'Fixed phase 3: ',kjj
                   endif
                enddo
! right hand side (rhs) contribution is
! - BP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij
                xcol(nz2)=xcol(nz2)-pham*mag*mass_of(ie,ceq)
             enddo ballel
! sum of mass in phase will be multiplied with delta-phase_amount
             if(sel.gt.0) then
                totam=totam+pham*pmi%xmol(sel)*mass_of(sel,ceq)
!                write(*,228)'totam 3: ',totam,pham,pmi%xmol(sel),&
!                     mass_of(sel,ceq)
             else
! sumwmol includes masses
                totam=totam+pham*pmi%sumwmol
!                totam=totam+pham*pmi%sumbmol
!                write(*,228)'totam 4: ',totam,pham,pmi%sumwmol,pmi%sumxmol,&
!                     pmi%sumbmol,zval
             endif
!             totam=zval
!             write(*,363)'363: ',jph,dncol,ncol,pham,totam,&
!                  (mamu(je),je=1,2),(xcol(k1),k1=1,nz2)
          enddo ballph
!
! in xcol are values summed over all phases and components
! copy summed columns to smat nrow
          nrow=nrow+1
          if(nrow.gt.nz1) then
             write(*,*)'too many equations 12A',nrow
             stop
          endif
          do ncol=1,nz2
             smat(nrow,ncol)=xcol(ncol)
          enddo
! add B^prescribed - B^current to rhs (right hand side)
          xxx=smat(nrow,nz2)
          smat(nrow,nz2)=smat(nrow,nz2)-cvalue+totam
!          write(*,363)'RHS: ',nrow,nz2,0,smat(nrow,nz2),xxx,cvalue,totam,&
!               cvalue-totam
          deallocate(xcol)
          if(abs(totam-cvalue).gt.ceq%xconv) then
!             write(*,266)'Unconverged condition N(A): ',sel,&
!                  cvalue,zval
!266       format(a,i3,2(1pe15.7))
             if(converged.lt.5) converged=5
          endif
          if(sel.eq.0) then
             write(*,363)'Condition B=fix',0,0,0,cvalue,totam
          else
             write(*,363)'Condition B(a)=fix',sel,0,0,cvalue,totam
          endif
       elseif(stvnorm.ne.2) then
! only normallizing of B with respect to mass (W) is allowed
          write(*,*)'Allowed normallizing with W only',stvix,stvnorm,cmix(2)
          gx%bmperr=9988; goto 1000
       else
!-------------------------------
! B=fix and B(A)=fix treated above as they have a "simple" summation, 
! Now handle all other cases like w(A)=fix, w(phase#set,A)=fix
! We must sum over all phases and constituents for the normallizing factor
! definition: W(A)=B(A)/B; 
! derivative: dW(A)=dB(A)/B - B(A)/N**2 *dB
! sum dB(A) and dB at the same time and multiply the sums with 1/B
! and -B(A)/B**2 in the end.
          if(cmix(3).eq.0) then
             write(*,*)'Condition BW=fix is illegal'
             gx%bmperr=9898; goto 1000
          elseif(cmix(4).eq.0) then
! condition is x(A)=fix
             sel=cmix(3); sph=0
          else
! condition is w(phase#set,A)=fix;  how to handle if phase#set not stable?
!             write(*,*)'Condition w(phase#set,A)=fix not yet allowed'
!             gx%bmperr=9898; goto 1000
             sel=cmix(5); sph=cmix(3); scs=cmix(4)
          endif
! two summations, zcol sums the term dN(A); xcol sums dN (as above)
          allocate(xcol(nz2))
          allocate(zcol(nz2))
          xcol=zero
          zcol=zero
          totam=zero
          zval=zero
          xval=zero
!          wallph: do jph=1,meqrec%nstph
!             jj=meqrec%stphl(jph)
! sum over all phases to handle conditions like x(phase#set,A)=fix
! as the phase#set may not be stable
          wallph: do jj=1,meqrec%nphase
             if(sph.eq.0) then
! skip this phase if not stable and condition not on a specific phase
                if(phr(jj)%stable.eq.0) cycle wallph
             elseif(sph.gt.0) then
                if(phr(jj)%iph.ne.sph .or. phr(jj)%ics.ne.scs) cycle wallph
             endif
             pmi=>phr(jj)
! mass of phase
!             pham=pmi%curd%amount(1)
             pham=pmi%curd%abnorm(2)
             wallel: do ie=1,meqrec%nrel
! calculate a term for each column to be multiplied with chemical potential
! we must sum xcol for all elemenets and add to zcol for element sel
! if sel=/=0 then we sum also zcol(sel) for all phases
                call calc_dyterms1(meqrec%nrel,ie,meqrec%tpindep,&
                     mamu,mag,mat,map,pmi)
                if(gx%bmperr.ne.0) goto 1000
                ncol=1
                bloop2: do je=1,meqrec%nrel
! Calculate one column for each component to be multiplied with chem.pot.
! components with fix chemical potential added to rhs, do not increment ncol!!!
                   do ke=1,meqrec%nfixmu
                      if(meqrec%mufixel(ke).eq.je) then
! different sign ???     xcol(nz2)=xcol(nz2)-pham*mamu(je)*meqrec%mufixval(ke)
                         xcol(nz2)=xcol(nz2)+&
                              pham*mamu(je)*meqrec%mufixval(ke)*mass_of(ie,ceq)
                         if(sel.eq.je) then
! different sign ???
                            zcol(nz2)=zcol(nz2)-&
                                 pham*mamu(je)*meqrec%mufixval(ke)*&
                                 mass_of(ie,ceq)
                         endif
                         cycle bloop2
                      endif
                   enddo
! mamu(B) = \sum_i \sum_j dM^a_B/dy_i dM^a_A z^a_ij
                   xcol(ncol)=xcol(ncol)-pham*mamu(je)*mass_of(ie,ceq)
                   if(sel.eq.ie) then
                      zcol(ncol)=zcol(ncol)-pham*mamu(je)*mass_of(ie,ceq)
                   endif
                   ncol=ncol+1
                enddo bloop2
! If T or P are variable
                if(tcol.gt.0) then
                   xcol(tcol)=xcol(tcol)+pham*mat*mass_of(ie,ceq)
                   if(sel.eq.ie) then
                      zcol(tcol)=zcol(tcol)+pham*mat*mass_of(ie,ceq)
                   endif
                endif
                if(pcol.gt.0) then
                   xcol(pcol)=xcol(pcol)+pham*map*mass_of(ie,ceq)
                   if(sel.eq.ie) then
                      zcol(pcol)=zcol(pcol)+pham*map*mass_of(ie,ceq)
                   endif
                endif
! last columns are amounts of element ie for all stable non-fix phase,
                notf=1
                do kph=1,meqrec%nstph
                   kjj=meqrec%stphl(kph)
                   if(phr(kjj)%phasestatus.eq.1) then
! phasestatus=1 for phases with variable amount, sum over all components
                      xcol(dncol+notf)=xcol(dncol+notf)+phr(kjj)%xmol(ie)*&
                           mass_of(ie,ceq)
                      if(ie.eq.sel) then
                         zcol(dncol+notf)=phr(kjj)%xmol(ie)*mass_of(ie,ceq)
                      endif
                      notf=notf+1
!                   elseif(phr(jj)%phasestatus.eq.2) then
!                      write(*,*)'Fixed phase 4: ',kjj
                   endif
                enddo
! right hand side (rhs) contribution is
! - NP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij * mass_ie
                xcol(nz2)=xcol(nz2)-pham*mag*mass_of(ie,ceq)
! the rhs sum the amount of component ie ... hm, different sign??
                if(sel.eq.ie) then
! different sign here and 2 lines above, why??
                   zcol(nz2)=zcol(nz2)-pham*mag*mass_of(ie,ceq)
                endif
             enddo wallel
             totam=totam+pham*pmi%sumwmol
! if sph=/=0 next line must be changed
             zval=zval+pham*pmi%xmol(sel)*mass_of(sel,ceq)
          enddo wallph
! in xcol is dB and in zcol dB(A) summed over all phases and components
! calculate the normallized values now
! xmat=dB(A)/B - B(A)*dB/B**2
          nrow=nrow+1
          if(nrow.gt.nz1) stop 'too many equations 12B'
!          write(*,228)'zcol: ',(zcol(jz),jz=1,nz2)
!          write(*,228)'xcol: ',(xcol(jz),jz=1,nz2)
!          write(*,228)'totam mm: ',totam,smat(nrow,nz2),cvalue,zval
! this was the missing statement!!
          zval=zval/totam
          do ncol=1,nz2
             smat(nrow,ncol)=(zcol(ncol)-xcol(ncol)*zval)/totam
          enddo
! add W^prescribed - W^current to rhs (right hand side)
          smat(nrow,nz2)=smat(nrow,nz2)-cvalue+zval
          deallocate(xcol)
          deallocate(zcol)
! check on convergence
          if(abs(zval-cvalue).gt.ceq%xconv) then
!             write(*,266)'Unconverged condition x(A): ',sel,&
!                  cvalue,zval
!266       format(a,i3,2(1pe15.7))
             if(converged.lt.5) converged=5
          endif
!          if(sph.eq.0) then
!             write(*,363)'Condition w(A)=fix',sel,0,0,cvalue,zval
!          else
!             write(*,363)'Condition w(phase#set,A)=fix',sph,sel,0,cvalue,zval
!          endif
       endif
!
!------------------------------------------------------------------
    case(13) ! Y
! Constituent fraction: phase#set, (subl.,) constituent index (over all subl)
       write(*,*)'Not implemented yet: ',stvix,stvnorm,cmix(2),cmix(3),cmix(4)
       gx%bmperr=9988; goto 1000
       nrow=nrow+1
       if(nrow.gt.nz1) stop 'too many equations 13A'
    end select
!
! loop if not the last condition
!    write(*,*)'Taking next condition',cmix(1)
    if(.not.associated(condition,lastcond)) goto 350
!=====================================================================
! debug output of system matrix
380 continue
!    do iz=1,nz1
!       write(*,228)'smat: ',(smat(iz,jz),jz=1,nz2)
!    enddo
    call lingld(nz1,nz2,smat,svar,nz1,ierr)
    if(ierr.ne.0) then
!       write(*,*)'Error solving system matrix',ierr
       gx%bmperr=4203; goto 1000
    endif
!    write(*,228)'svar: ',(svar(iz),iz=1,nz1)
!---------
! copy the chemical potentials, take care of fixed values ....
! new potentials are in svar(1..meqrec%nrel-meqrec%nfixmu)
    iz=1
    notf=1
    setmu: do ik=1,meqrec%nrel
       if(notf.le.meqrec%nfixmu) then
          if(ik.eq.meqrec%mufixel(notf)) then
! this potential is fixed, skip incrementing "iz", curmu(ik) is a condition
             notf=notf+1
             cycle setmu
          endif
       endif
       if(abs(svar(iz)-meqrec%curmu(ik)).gt.ceq%xconv) then
!          write(*,387)'Not converged: ',iz,ik,svar(iz),meqrec%curmu(ik)
387       format(a,2i3,6(1pe12.4))
          converged=7
       endif
       meqrec%curmu(ik)=svar(iz)
       iz=iz+1
    enddo setmu
    ioff=meqrec%nrel-meqrec%nfixmu+1
!------------
! update T and P if variable
    if(meqrec%tpindep(1)) then
       xxx=ceq%tpval(1)
! limit changes in T to +/-half its current value
!       write(*,*)'Tcol column check: ',tcol,ioff
       if(abs(svar(ioff)/ceq%tpval(1)).gt.0.5D0) then
          svar(ioff)=sign(0.5D0*ceq%tpval(1),svar(ioff))
       endif
       deltat=svar(ioff)
       ceq%tpval(1)=ceq%tpval(1)+deltat
!       write(*,388)'New T: ',tcol,ioff,ceq%tpval(1),xxx,deltat
388    format(a,2i3,3(1pe12.4))
       if(ceq%tpval(1).le.zero) then
          write(*,*)'Attempt to set temperature negative!!!'
          gx%bmperr=9996; goto 1000
       endif
       ioff=ioff+1
    endif
    if(meqrec%tpindep(2)) then
       xxx=ceq%tpval(2)
       ceq%tpval(2)=ceq%tpval(2)+svar(ioff)
! one should also limit too big changes in T
       if(ceq%tpval(2).le.zero) then
          write(*,*)'Pressure set negative!!!'
          gx%bmperr=9996; goto 1000
       endif
       ioff=ioff+1
    endif
!------------
! update phase amounts, take care of fixed phases ....
! the change in amounts are in svar(ioff+...)
    negamph=0
    negam=0
    irem=0
! dncol+1 should be the first Delta_phase-amount
!    write(*,*)'Column check dncol: ',dncol+1,ioff
    ioff=dncol+1
    phamount2: do jph=1,meqrec%nstph
       jj=meqrec%stphl(jph)
       phr(jj)%curd%damount=zero
!       kkz=test_phase_status(phr(jj)%iph,phr(jj)%ics,xxx,ceq)
       kkz=phr(jj)%phasestatus
!       write(*,379)'phase status: ',phr(jj)%iph,phr(jj)%ics,kkz
379    format(a,2i3,2i5)
       if(kkz.eq.1) then
! phase is entered so its amount can change
!          phs=phr(jj)%curd%amount(1)
          phs=phr(jj)%curd%amfu
! limit change in amount of phase
          deltaam=svar(ioff)
          if(abs(deltaam).gt.meqrec%nnorm) then
             deltaam=sign(meqrec%nnorm,deltaam)
          endif
          if(phr(jj)%curd%amfu-deltaam.le.zero) then
             if(meqrec%nstph.eq.1) then
! this is the only stable phase!  cannot have negative or zero amount!
                deltaam=phr(jj)%curd%amfu-1.0D-2
             endif
          endif
          phf=phr(jj)%curd%amfu-deltaam
!          write(*,363)'Stable phase: ',jj,phr(jj)%iph,phr(jj)%ics,&
!               phf,phs,phr(jj)%prevam
          phr(jj)%curd%damount=deltaam
          ioff=ioff+1
       elseif(kkz.eq.2) then
! phase is fix, there is no change in its amounts
!          write(*,*)'fixed phase: ',jj,phr(jj)%iph,kkz
!          phf=phr(jj)%curd%amount(1)
          phf=phr(jj)%curd%amfu
       else
! phase is dormant or suspended, must not be stable!!!!
          write(*,*)'A dormant or suspended phase stable!!!',kkz
          write(*,*)phr(jj)%iph,phr(jj)%ics
          gx%bmperr=7777; goto 1000
       endif
! make sure the driving force of stable phases to zero
!       write(*,363)'Stable phase: ',jj,phr(jj)%iph,phr(jj)%ics,phf,phs,&
!            phr(jj)%prevam
       phr(jj)%curd%dgm=zero
       if(phf.lt.zero) then
! phase has negative amount, NOT ALLOWED if it is the only stable phase 
          if(meqrec%nstph-meqrec%nfixph.eq.1) then
             write(*,*)'Trying to remove the only stable phase ',&
                  phr(jj)%curd%amfu
!                  phr(jj)%curd%amount(1)
!             phf=0.5D0*phr(jj)%curd%amount(1)
             phf=0.5D0*phr(jj)%curd%amfu
          else
!             write(*,363)'Phase with negative amount: ',jj,0,0,&
!                  phf,phs,phr(jj)%prevam
!             if(phf.lt.-1.0D-2) phf=1.0D-6
             if(phf.lt.-1.0D-2) phf=zero
!             if(phf.lt.-5.0D0*abs(phs)) then
!                write(*,*)'Force removal of ',jj,phf,phs
!                irem=jj
!                goto 1000
!             endif
             if(phr(jj)%prevam.lt.zero) then
! remove this phase right away if negative amount previous iteration also
                irem=jj
! jumping to 1000 here means constitutions not changed
!                goto 1000
             else
! mark this phase had negative amount this iteration
! PROBLEM removing one of two composition sets of the same phase,
! (miscibility gap), they may change which have negative amount each iteration
                phr(jj)%prevam=-one
                irem=jj
             endif
          endif
       else ! phase has positive amount, mark in prevam
          phr(jj)%prevam=one
       endif
! store the new phase fraction (moles formula units)
!       phr(jj)%curd%amount(1)=phf
       phr(jj)%curd%amfu=phf
    enddo phamount2
!-------------------------------------------------------
! After solving the system matrix and updating the chemical potentials,
! the phase amounts and possibly T and P we correct constitions of all phases
! - Now calculate correction of constituent fractions for all phases
! See BoJ thesis eq. 30 (also in metastable phases) (paper I)
! At the same time calculate the driving force for metastable phases
    ycorr=zero
    ycormax2=zero
! to handle charge balance correction of constituent fractions
    chargerr=zero
! chargerr fitted to fastest convergence using the ou test case
!    chargefact=1.0D-1 requires more than 100 iterations
!    chargefact=one requires more than 100 iterations
! this value requires about 40 iteration
    chargefact=5.0D-1
! kk is incremented for each stable phase
    kk=1
! iadd is set to the unstable phase with largest positive driving force
! dgmmax is the largest psoitive driving force
    iadd=0
    dgmmax=zero
    ysmm=zero
!-----------------------------------------------------
! Update the constitutions.  If rem>0 recalculate with this phase amount zero
    if(irem.gt.0) goto 1000
!--------------------------
    lap: do jj=1,meqrec%nphase
! The current chemical potentials are in curmu(i)
!       write(*,*)'Phase: ',phr(jj)%iph,phr(jj)%ics,phr(jj)%curd%amfu
       if(jj.eq.meqrec%stphl(kk)) then
! increment kk but do not make it larger than meqrec%nstph
          kk=min(kk+1,meqrec%nstph)
       else ! phase is not stable
! calculate driving force for unstable phases. First calculate the sum
! of the current phase composition and the calculated chemical potentials
          gsurf=zero; summ=zero
          do ie=1,meqrec%nrel
             gsurf=gsurf+phr(jj)%xmol(ie)*meqrec%curmu(ie)
             summ=summ+phr(jj)%xmol(ie)
          enddo
          gsurf=gsurf/summ
! calculte G_m plus any deltat and deltap terms
          dgm=phr(jj)%curd%gval(1,1)
!          write(*,33)'surf: ',jj,gsurf,summ,phr(jj)%curd%abnorm(1),&
!               phr(jj)%curd%qqsave(1),dgm
33        format(a,i3,6(1pe12.4))
          if(meqrec%tpindep(1)) then
             dgm=dgm+phr(jj)%curd%gval(2,1)*deltat
!             write(*,13)'Delta_T: ',jj,&
!                  dgm,phr(jj)%curd%gval(2,1),deltat
13           format(a,i3,6(1pe12.4))
          endif
          if(meqrec%tpindep(2)) then
             dgm=dgm+phr(jj)%curd%gval(3,1)*deltap
             write(*,13)'Delta_P: ',jj,&
                  dgm,phr(jj)%curd%gval(3,1),deltap
          endif
! scale dgm per mole atoms
          dgm=gsurf-dgm/phr(jj)%curd%abnorm(1)
!          write(*,11)
          if(dgm.gt.dgmmax) then
!             if(test_phase_status(phr(jj)%iph,phr(jj)%ics,xxx,ceq).eq.1) then
!             write(*,11)'dgm:  ',jj,dgm,dgmmax,phr(jj)%curd%qqsave(1),&
!                  phr(jj)%curd%abnorm(1)
11           format(a,i3,6(1pe12.4))
             if(phr(jj)%phasestatus.eq.1) then
! phase is entered, can have status changed
! if this is another constitution set of an already stable phase then check
! below if the constitution of this phase is very similar to the stable one
                iadd=jj
                dgmmax=dgm
             endif
          endif
          phr(jj)%prevdg=dgm
          phr(jj)%curd%dgm=dgm
       endif
! Update constituent fractions for all phases, stable or not
! if phr(jj)%xdone=1 then phase has no composition variation
!       write(*,*)'xdone: ',phr(jj)%xdone
       if(phr(jj)%xdone.eq.1) cycle
!----------------------------------------------------
       allocate(cit(phr(jj)%idim))
       cit=zero
       if(meqrec%tpindep(1)) then
! variable T, code copied from calc_dyterm, cit(nj) used below
!          write(*,44)'index 1: ',jj,phr(jj)%ncc,phr(jj)%idim,&
!               size(phr(jj)%invmat)
          do jy=1,phr(jj)%ncc
             sum=zero
             do iy=1,phr(jj)%ncc
!                write(*,44)'index 1A: ',jj,jy,iy,phr(jj)%ncc,sum
                ys=phr(jj)%invmat(iy,jy)
!                write(*,44)'index 1B: ',jj,jy,iy,0,sum
                ys=ys*phr(jj)%curd%dgval(1,iy,1)
!                write(*,44)'index 1C: ',jj,jy,iy,0,sum
                ys=ys*phr(jj)%curd%dgval(2,iy,1)
!                write(*,44)'index 1D: ',jj,jy,iy,0,sum
                sum=sum+phr(jj)%invmat(iy,jy)*phr(jj)%curd%dgval(2,iy,1)
             enddo
             cit(iy)=sum*deltat
!             write(*,44)'index 2: ',jj,jy,iy,0,sum
!44           format(a,4i3,6(1pe12.4))
          enddo
!! end copy
!          write(*,*)'Adding contribution from variable T to delta-y',&
!               phr(jj)%ncc
       endif
!
       moody: do nj=1,phr(jj)%ncc
          ys=zero
          do nk=1,phr(jj)%ncc
             pv=zero
             do nl=1,meqrec%nrel
! curmu(nl) is the chemical potential of element nl (divided by RT)
! phr(jj)%dxmol(nl,nk) is the derivative of component nl wrt constituent nk
                pv=pv+meqrec%curmu(nl)*phr(jj)%dxmol(nl,nk)
             enddo
             pv=pv-phr(jj)%curd%dgval(1,nk,1)
             ys=ys+phr(jj)%invmat(nj,nk)*pv
          enddo
          if(phr(jj)%chargebal.eq.1) then
! For charged phases add a term  phr(jj)%invmat(phr(jj)%idim,phr(jj)%idim)*Q
             ys=ys-chargefact*phr(jj)%invmat(nj,phr(jj)%idim)*phr(jj)%charge
             if(abs(phr(jj)%charge).gt.chargerr) then
                chargerr=abs(phr(jj)%charge)
                signerr=phr(jj)%charge
             endif
!             write(*,*)'Charge: ',jj,phr(jj)%charge
!             write(*,262)nj,phr(jj)%idim,ys,ysss,&
!                  phr(jj)%invmat(nj,phr(jj)%idim),phr(jj)%charge
          endif
          ycorr(nj)=ys+cit(nj)
!          write(*,51)'at 1A: ',jj,nj,phr(jj)%iph,phr(jj)%ics,phr(jj)%stable,&
!               ys,cit(nj),phr(jj)%curd%yfr(nj)
!51        format(a,5i4,6(1pe12.4))
          if(ycorr(nj).gt.ycormax2) then
             ycormax2=ycorr(nj)
          endif
! if the change in any constituent fraction larger than xconv continue iterate
          if(abs(ys).gt.ceq%xconv) then
!             write(*,*)'Convergence criteria, phase/const: ',jj,nk
             if(phr(jj)%stable.eq.0) then
! Phase is not stable
                if(abs(ys).gt.1.0D1*phr(jj)%curd%yfr(nj)) then
! for unstable phases the corrections must be smaller than ...????
                   if(converged.lt.3) then
                      converged=3
                      yss=ys
                      yst=phr(jj)%curd%yfr(nj)
                   endif
                elseif(abs(ys).gt.1.0D2*ceq%xconv) then
! maybe accept 100 times larger correction than for stable phases
!                   write(*,107)'metast ph ycorr: ',ys,phr(jj)%curd%yfr(nj)
                   if(converged.lt.2) then
                      converged=2
                      yss=ys
                      yst=phr(jj)%curd%yfr(nj)
                   endif
                else
                   if(converged.eq.0) then
                      converged=1
                      yss=ys
                      yst=phr(jj)%curd%yfr(nj)
                   endif
                endif
             elseif(converged.lt.4) then
! large correction in fraction of constituent fraction of stable phase
                converged=4
                yss=ys
                yst=phr(jj)%curd%yfr(nj)
             endif
          elseif(phr(jj)%stable.eq.1) then
! check to find good convergence criteria in Re-V test case
             if(abs(ycorr(nj)).gt.ysmm) then
                jmaxy=jj
                ysmm=abs(ycorr(nj))
                ysmt=phr(jj)%curd%yfr(nj)
             endif
          endif
! fetch constituent fractions directly from phase_varres record
          yprev=phr(jj)%curd%yfr(nj)
          yarr(nj)=phr(jj)%curd%yfr(nj)+ycorr(nj)
!          write(*,51)'at 3A: ',jj,nj,phr(jj)%iph,phr(jj)%ics,phr(jj)%stable,&
!               ys,cit(nj),phr(jj)%curd%yfr(nj)
          if(yarr(nj).lt.ymin) then
! this added to avoid too drastic jumps in small fractions
! The test case ccrfe1.BMM needs this
             if(yprev.gt.ylow) then
!                write(*,*)'Applying fraction change limitation'
                yarr(nj)=0.9*ylow
             elseif(phase_bit(phr(jj)%iph,PHGAS)) then
!             elseif(phase_bit(phr(jj)%iph,PHGAS,ceq)) then
! for gas phase one must allow smaller constituent fractions
                if(yarr(nj).lt.ymingas) then
                   yarr(nj)=ymingas
                endif
             else
                yarr(nj)=ymin
             endif
          endif
          if(yarr(nj).gt.one) then
             yarr(nj)=one
          endif
!          write(*,51)'at 5A: ',jj,nj,phr(jj)%iph,phr(jj)%ics,phr(jj)%stable,&
!               ys,cit(nj),phr(jj)%curd%yfr(nj),yarr(nj)
       enddo moody
!
       ycormax(jph)=ycormax2
!       write(*,263)'yarr3: ',jj,(phr(jj)%curd%yfr(jz),jz=1,phr(jj)%ncc)
!       write(*,*)'phase: ',phr(jj)%iph,phr(jj)%ics
!       write(*,107)'ycorr: ',(ycorr(jz),jz=1,phr(jj)%ncc)
!       write(*,107)'yarr: ',(yarr(i),i=1,phr(jj)%ncc)
! maybe one can set constitution directly but some internal arrays need
! updating for each new constitition so better use this subroutine
       phs=phr(jj)%curd%abnorm(1)
!       write(*,263)'yarr4: ',jj,(yarr(jz),jz=1,phr(jj)%ncc)
       call set_constitution(phr(jj)%iph,phr(jj)%ics,yarr,qq,ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,264)'abnorm: ',jj,phr(jj)%iph,phr(jj)%ics,&
!            phs,phr(jj)%curd%abnorm(1),yarr(1)
       deallocate(cit)
    enddo lap
! finished correction of all constituent fractions
!-------------------------------------------------------
!    do jph=1,meqrec%nstph
!       jj=meqrec%stphl(jph)
!       write(*,393)'Stable phase: ',phr(jj)%iph,phr(jj)%ics,&
!            phr(jj)%curd%amfu
!    enddo
393 format(a,2i4,6(1pe12.4))
! check if fraction corrections in stable phases increases
! it solved a problem in ReV when fractions initially changed very little
! but the change increased each iteration
    if(meqrec%noofits.gt.8) then
       increase=0
    elseif(abs(ysmm).gt.prevmaxycorr) then
! do this check only for the first 8 iterations
       increase=1
!       write(*,265)increase,ysmm,prevmaxycorr
!265    format('*** max stable phase ycorr: ',i3,2(1pe12.4))
    endif
    prevmaxycorr=abs(ysmm)
!-------------------------------------------------------
! check charge balance, must be 100 times better than fractions
! otherwise strange chemical potentials, why??
    if(chargerr.gt.1.0D-2*ceq%xconv) then
!       write(*,*)'Charge error: ',signerr
       if(converged.lt.6) converged=6
    endif
!-------------------------------------------------------
! - check if conditions fullfilled xknown no known ....
!    call calc_molmass(xxmol,wmass,tmol,tmass,ceq)
!    if(gx%bmperr.ne.0) goto 1000
!    do ik=1,meqrec%nrel
! check on convergence
!       if(abs(xxmol(ik)-xknown(ik)).gt.ceq%xconv) then
!          write(*,266)'Unconverged element fraction: ',&
!               ik,xxmol(ik),xknown(ik),ceq%xconv
!266       format(a,i3,3(1pe15.7))
!          if(converged.lt.5) converged=5
!       endif
!    enddo
!    write(*,107)'total mol/mass: ',tmol,tmass
! write Gibbs energy and potentials ....
!    gg=zero
!    do i=1,noel()
!       gg=gg+xxmol(i)*meqrec%curmu(i)
!    enddo
!    write(*,107)'gg: ',gg,(meqrec%curmu(i),i=1,noel())
!    write(*,107)'xx: ',tmol,(xxmol(i),i=1,noel())
!    write(*,107)'xpf: ',xall,(xfrac(i),i=1,noel())
!    write(*,*)'Convergence criteria: ',converged
    if(converged.eq.3) then
! force some iterations with large fraction variations in unstable phases
!       write(*,267)'End of iteration: ',meqrec%noofits,converged,&
!            increase,yss,yst
       level3=level3+1
    elseif(converged.eq.4) then
! this meas large fraction variations in stable phases
!       write(*,267)'End of iteration: ',meqrec%noofits,converged,&
!            increase,yss,yst
!267    format(a,3i4,2(1pe12.4))
       level3=0
    else
!       write(*,267)'End of iteration: ',meqrec%noofits,converged,increase
       level3=0
    endif
!----------------------------------------------
! continue iterate if phase change or not converged
    if(iadd.gt.0) then
! check if phase to be added is already stable as another composition set
! This check should maybe be above as maybe another phase want to be stable??
       if(same_composition(iadd,phr,meqrec,ceq,dgm)) iadd=0
    endif
    if(meqrec%noofits.gt.2 .and. (irem.gt.0 .or. iadd.gt.0)) then
! if a phase have negative amount remove it or if a phase has positive
! driving force add it
       write(*,363)'Try phase change remove/add: ',irem,iadd,0,phf,dgmmax
       goto 1100
    endif
! converged=1 or 2 means constituent fraction in metastable phase not converged
    if(converged.gt.3) goto 100
! converged 3 means a constituent fraction of a stable phase change a lot
    if(converged.eq.3 .and. level3.lt.4) goto 100
! always force 4 iterations, there is a minimum above forcing 9 iterations.
    if(meqrec%noofits.lt.4) goto 100
    if(increase.ne.0) then
! continue if corrections in constituent fractions in stable phases increases
       goto 100
    endif
!------------------------
! equilibrium calculation converged, do some common thing
    goto 800
!
!==============================================================
! equilibrium calculation converged, save chemical potentials (svar*RT)
800 continue
    ceq%rtn=globaldata%rgas*ceq%tpval(1)
    do ie=1,meqrec%nrel
!       ceq%complist(ie)%chempot(1)=svar(ie)*ceq%rtn
       ceq%complist(ie)%chempot(1)=meqrec%curmu(ie)*ceq%rtn
!       write(*,*)'Chempot: ',meqrec%curmu(ie),svar(ie)
    enddo
! list stable phases on exit
!    do jph=1,meqrec%nstph
!       jj=meqrec%stphl(jph)
!       write(*,393)'Stable phase Z: ',phr(jj)%iph,phr(jj)%ics,&
!            phr(jj)%curd%amfu
!    enddo
!----------------------
! save inverted phase matrix and more for future use when calculating H.T etc
! If already allocated then dealloc/alloc as number of constituents can change
    do jj=1,meqrec%nphase
       if(allocated(phr(jj)%curd%cinvy)) then
          deallocate(phr(jj)%curd%cinvy)
          deallocate(phr(jj)%curd%cxmol)
          deallocate(phr(jj)%curd%cdxmol)
       endif
! why is the dimension if invmat so different???
       ie=phr(jj)%idim
!       ie=int(sqrt(real(size(phr(jj)%invmat)))+0.1)
!       write(*,*)'Size: ',ie,phr(jj)%ncc
       allocate(phr(jj)%curd%cinvy(ie,ie))
       allocate(phr(jj)%curd%cxmol(meqrec%nrel))
       allocate(phr(jj)%curd%cdxmol(meqrec%nrel,phr(jj)%ncc))
       phr(jj)%curd%cinvy=phr(jj)%invmat
       phr(jj)%curd%cxmol=phr(jj)%xmol
       phr(jj)%curd%cdxmol=phr(jj)%dxmol
!----------------------
    enddo
1000 continue
   if(gx%bmperr.ne.0) then
      ceq%status=ibset(ceq%status,EQFAIL)
!      write(*,*)'minimization error: ',gx%bmperr
   elseif(irem.eq.0 .and. iadd.eq.0) then
      write(*,1001)meqrec%noofits
1001  format('Solution found after ',i5,' iterations.')
   endif
! jump here if phase change
1100 continue
   deallocate(smat)
   deallocate(svar)
   return
! too many iterations
1200 continue
   write(*,*)'Too many iterations: ',meqrec%noofits,ceq%maxiter
   gx%bmperr=4204
   goto 1000
 end subroutine meq_sameset

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} sumxmol=
  subroutine meq_onephase(nrel,pmi,ceq)
! this subroutine calculates new constituent fractions for a phase iph+ics
! with given T, P and chemical potentials for the components 
! For ionic liquids the sites on the sublattices varies with composition
! THIS IS A FIRST VERSION WITHOUT ANY TRICKS FOR SPEED
    implicit none
    integer nrel
    TYPE(meq_data), pointer :: pmi
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer mph,ik,iph,ics,jz,iz,jk,ierr,kk,kkk,ll,lokcs,ncc,loksp,ncl
    integer nd1,nd2,neq,nochange,nsl,nspel,nv,ncon
! needed for call to get_phase_data
    integer, dimension(maxsubl) ::  nkl
    integer, dimension(maxconst) :: knr
    double precision, dimension(5) :: qq
    double precision, dimension(maxsubl) :: sites
! needed for call to get_species_data
    integer, dimension(maxspel) :: ielno
    double precision, dimension(maxspel) :: stoi
! minimal y, charge
    double precision, parameter :: ymin=1.0D-12,ymingas=1.0D-30,qeps=1.0D-30
! derivative of moles of component wrt y_ks
    double precision, dimension(maxel) :: addmol
! for mass balance and charge
    double precision, dimension(maxconst) :: yarr(maxconst),dqsum
! phase matrix, its inverse is returned as part of pmi
    double precision, dimension(:,:), allocatable :: pmat
    double precision qsp,sumsit,ykvot,ysum,qsum,spmass
!
    iph=pmi%iph
    ics=pmi%ics
! extract phase structure
    call get_phase_data(iph,ics,nsl,nkl,knr,yarr,sites,qq,ceq)
!    write(*,*)'in meq_onephase: ',iph
    if(gx%bmperr.ne.0) then
       write(*,*)'get_phase_data error in meq_onephase',iph,gx%bmperr
       goto 1000
    endif
! make sure all fractions >ymin and sums in alll sublattices equail to unity
    nochange=0
    ncc=0
    do ll=1,nsl
       ysum=zero
       ncl=ncc
       do ik=1,nkl(ll)
          ncc=ncc+1
          if(yarr(ncc).lt.ymin) then
!             if(phase_bit(iph,PHGAS,ceq)) then
             if(phase_bit(iph,PHGAS)) then
                if(yarr(ncc).lt.ymingas) then
                   yarr(ncc)=ymingas
                   nochange=1
                endif
             else
                nochange=1
                yarr(ncc)=ymin
             endif
          endif
          ysum=ysum+yarr(ncc)
       enddo
       ykvot=one/ysum
       if(abs(ykvot-one).gt.ymingas) then
          nochange=1
          do ik=1,nkl(ll)
             yarr(ncl+ik)=yarr(ncl+ik)*ykvot
          enddo
       endif
    enddo
    if(nochange.ne.0) then
! if constitution changed save it. qq will be updated automatically
       call set_constitution(iph,ics,yarr,qq,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'never never error 17'
          goto 1000
       endif
    endif
    if(pmi%stable.eq.1) then
!       write(*,12)iph,ics,(yarr(i),i=1,ncc)
12     format('Const: ',2i3,6(1pe10.3))
    endif
!    write(*,17)'curd: ',pmi%curd%amfu,pmi%curd%netcharge,pmi%curd%abnorm
    if(phase_bit(iph,PHEXCB)) then
! If external charge balance phase matrix has one more line+column
       pmi%chargebal=1
       nd1=ncc+1
       pmi%charge=qq(2)
    else
       pmi%chargebal=0
       nd1=ncc
       pmi%charge=zero
    endif
!--------------------------
! sublattice rows, nd2=nd1+1 because I use Lukas matrix inverter
    nd1=nd1+nsl
    nd2=nd1+1
! Allocate phase matrix, one extra dimension if external charge balance
! last column of pmat is left hand side ?? (reminicent from Lukas program)
    allocate(pmat(nd1,nd2))
! return dimension of pmi%invmat
    if(pmi%idim.eq.0) then
       pmi%idim=nd1
       pmi%ncc=ncc
       allocate(pmi%invmat(nd1,nd1))
! meqrec is not available in this subroutine but meqrec%nrel passed in call
       allocate(pmi%xmol(nrel))
       allocate(pmi%dxmol(nrel,ncc))
!       allocate(pmi%nmol(nrel))
!       write(*,*)'Allocated phase matrix: ',nd2,noel(),ncc
    endif
! value of RT
    ceq%rtn=globaldata%rgas*ceq%tpval(1)
!--------------------------------------------------
! now treat different phase types
    call get_phase_variance(iph,nv)
    if(nv.eq.0) then
!----------------------------------------------- stoichiometric phase
! For stoichiometric phases calculate just G with T and P derivatives
! and driving force.  All pmi%dxmol=zero but one must also calculate 
! pmi%xmol and save it for all future iterations
! It must also be saved in curd%abnorm(1)
       if(pmi%xdone.eq.1) goto 90
       pmi%xmol=zero
       pmi%dxmol=zero
       pmi%sumxmol=zero
       pmi%sumwmol=zero
!       pmi%nmol=zero
!       pmi%sumnmol=zero
!       pmi%sumbmol=zero
       sumsit=zero
       do ll=1,nsl
          sumsit=sumsit+sites(ll)
       enddo
       kkk=0
       sublatt: do ll=1,nsl
          allconst: do ik=1,nkl(ll)
             kkk=kkk+1
             loksp=knr(kkk)
             call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp)
             addmol=zero
             do jz=1,nspel
                addmol(jz)=stoi(jz)
             enddo
             dqsum(kkk)=qsp
             qsum=qsum+qsp
             do jz=1,nspel
                if(ielno(jz).gt.0) then
! ignore vacancies, taken care of by using sumsit=qq(1) above
!                   pmi%dxmol(ielno(j),kkk)=sites(ll)*addmol(j)/sumsit
                   pmi%dxmol(ielno(jz),kkk)=zero
                   pmi%xmol(ielno(jz))=pmi%xmol(ielno(jz))+&
                        sites(ll)*addmol(jz)/sumsit
!                   pmi%nmol(ielno(jz))=pmi%nmol(ielno(jz))+&
!                        sites(ll)*addmol(jz)
                endif
             enddo
          enddo allconst
       enddo sublatt
! meqrec is not available in this subroutine
       do iz=1,nrel
          pmi%sumxmol=pmi%sumxmol+pmi%xmol(iz)
!          write(*,*)'sumwmol 1: ',pmi%xmol(iz),mass_of(iz,ceq)
          pmi%sumwmol=pmi%sumwmol+pmi%xmol(iz)*mass_of(iz,ceq)
!          pmi%sumnmol=pmi%sumnmol+pmi%nmol(iz)
!          pmi%sumbmol=pmi%sumbmol+pmi%nmol(iz)*mass_of(iz,ceq)
       enddo
!       write(*,88)pmi%iph,ll,(pmi%xmol(jz),jz=1,noel())
!       enddo sublatt
! some stoichiometric phases have wrong number of moles like Cr2VC2 ...
!       write(*,92)'onephase 1: ',pmi%iph,nsl,pmi%xdone,pmi%sumxmol,&
!            pmi%sumnmol,qq(1)
92     format(a,3i3,6(1pe12.4))
!       write(*,88)pmi%iph,pmi%ics,(pmi%xmol(jz),jz=1,noel())
88     format('x: ',2i2,6(1pe12.4))
!       write(*,*)'stoik: ',pmi%curd%abnorm(1),pmi%sumxmol
! phase_varres(lokcs)%abnorm already set by set_constitution
       pmi%xdone=1
!
90     continue
       call calcg(iph,ics,2,lokcs,ceq)
       if (gx%bmperr.ne.0) then
          write(*,*),'calcg error in meq_onephase ',iph,gx%bmperr
          goto 1000
       endif
! maybe some common ending
       goto 900
    endif
!--------------------------------------------- zero some arrays, ideal phase
    pmi%xmol=zero
    pmi%dxmol=zero
    pmi%sumxmol=zero
    pmi%sumwmol=zero
!    pmi%nmol=zero
!    pmi%sumnmol=zero
!    pmi%sumbmol=zero
    pmi%xdone=-1
!    if(phase_model(iph,ics,PHID,ceq)) then
!    if(phase_bit(iph,PHID,ceq)) then
    if(phase_bit(iph,PHID)) then
!--------------------------------------------- ideal phase (subst, no excess)
!       write(*,*)'Phase is ideal'
! special treatment of ideal phase (gas), sites assumed to be unity
! 1. Calculate M_i and dM_i/dy^s_k and the net charge charge Q and dQ/dy^s_k
       pmi%xmol=zero
       pmi%dxmol=zero
       qsum=zero
       dqsum=zero
       ncon=0
       do ik=1,nkl(1)
          loksp=knr(ik)
          call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp)
          addmol=zero
!          sumstoi=zero
          do jk=1,nspel
             addmol(jk)=stoi(jk)
!             sumstoi=sumstoi+stoi(jk)
          enddo
          dqsum(ik)=qsp
          qsum=qsum+qsp*yarr(ik)
!          addmol=addmol/sumstoi
          do jk=1,nspel
             pmi%dxmol(ielno(jk),ik)=addmol(jk)
             pmi%xmol(ielno(jk))=pmi%xmol(ielno(jk))+addmol(jk)*yarr(ik)
          enddo
          ncon=ncon+1
       enddo
! meqrec is not available in this subroutine
       do ik=1,nrel
          pmi%sumxmol=pmi%sumxmol+pmi%xmol(ik)
!          write(*,*)'sumwmol 2: ',pmi%xmol(ik),mass_of(ik,ceq)
          pmi%sumwmol=pmi%sumwmol+pmi%xmol(ik)*mass_of(ik,ceq)
       enddo
!       pmi%sumnmol=pmi%sumxmol
!       pmi%sumbmol=pmi%sumwmol
!       write(*,92)'onephase 2: ',pmi%iph,nsl,pmi%xdone,pmi%sumxmol,&
!            pmi%sumnmol,qq(1)
!       do ik=1,noel()
!          write(*,17)'dmx: ',(pmi%dxmol(ik,jk),jk=1,nkl(1))
!       enddo
17     format(a,6(1pe12.3))
! now calculate G and all 1st and 2nd derivatives
! This can be speeded up as all 2nd derivatives of constituents are RT/y
! The calculated values are used also in other parts of the code 
      call calcg(iph,ics,2,lokcs,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error calculating ideal gas'
          goto 1000
       endif
! calculate phase matrix elements
! temporarely ignore that the phase matrix is symmetric
! ceq%phase_varres(lokcs)%...
! gval(1:6,1) are G, G.T, G.P, G.T.T, G.T.P, G.P.P
! dgval(1,1:N,1) are first derivatives of G wrt constituent 1:N
! dgval(2,1:N,1) are second derivatives of G wrt constituent 1:N and T
! dgval(3,1:N,1) are second derivatives of G wrt constituent 1:N and P
! d2gval(ixsym(N*(N+1)/2),1) are 2nd derivatives of G wrt constituents N and M
! Last index is other properties than G like TC, BMAGN etc.
       pmat=zero
       do ik=1,nkl(1)
          do jk=ik,nkl(1)
             pmat(ik,jk)=ceq%phase_varres(lokcs)%d2gval(ixsym(ik,jk),1)
             if(jk.gt.ik) pmat(jk,ik)=pmat(ik,jk)
          enddo
       enddo
       neq=nkl(1)
! add one column and row for each sublattice (here only one)
       neq=neq+1
       do jk=1,neq-1
          pmat(jk,neq)=one
          pmat(neq,jk)=one
       enddo
       if(pmi%chargebal.eq.1) then
! if external charge balance add one column and one row
          neq=neq+1
          do jk=1,nkl(1)
! this is the row
             pmat(jk,neq)=dqsum(jk)
! this is the column
             pmat(neq,jk)=dqsum(jk)
          enddo
       endif
!       do jk=1,neq
!          write(*,17)'pmat: ',(pmat(ik,jk),ik=1,neq)
!       enddo
!       write(*,*)'Phase matrix: ',nd1,neq,pmi%chargebal
! invert the phase matrix (faster subroutine should be used)
       call mdinv(nd1,nd2,pmat,pmi%invmat,neq,ierr)
       if(ierr.eq.0) then
!          write(*,*)'Phase matrix singular',pmi%iph,pmi%ics,pmi%ncc,ierr
          gx%bmperr=4205; goto 1000
       endif
!       pmi%invmat=zero
! maybe some common ending
       goto 900
    endif
!---------------------------------------------- no analythical 2nd derivatives
! phases with models with no analythical second derivatives ....
!    if(phase_model(iph,ics,PHNODGDY2,ceq)) then
!    if(phase_bit(iph,PHNODGDY2,ceq)) then
    if(phase_bit(iph,PHNODGDY2)) then
!       write(*,*)'Models without 2nd derivatives not implemented'
       gx%bmperr=4206; goto 1000
    endif
!------------------------------------------------- ionic liquid phase
!    if(phase_model(iph,ics,PHIONLIQ,ceq)) then
!    if(phase_bit(iph,PHIONLIQ,ceq)) then
    if(phase_bit(iph,PHIONLIQ)) then
       write(*,*)'Ionic liquid model not implemented'
       gx%bmperr=4207; goto 1000
    endif
!------------------------------------------------- all other phase models
! For all other phases calculate G and all first and second derivatives
! for current composition
300 continue
! Calculate M_i and dM_i/dy^s_k and the net charge charge Q and dQ/dy^s_k
!   call get_phase_data(iph,ics,nsl,nkl,knr,yarr,sites,qq,ceq)
! how to normalize xmol?  use qq(1)!!, it handels vacancies .... ????
    sumsit=one
!    sumsit=pmi%curd%abnorm(1)
!    sumsit=qq(1)
!    sumsit=zero
!    do ll=1,nsl
!       sumsit=sumsit+sites(ll)
!    enddo
    pmi%xmol=zero
    pmi%dxmol=zero
    qsum=zero
    dqsum=zero
    ncon=0
    pmi%sumxmol=zero
    pmi%sumwmol=zero
!    pmi%nmol=zero
!    pmi%sumnmol=zero
!    pmi%sumbmol=zero
    subll: do ll=1,nsl
       constll: do ik=1,nkl(ll)
          ncon=ncon+1
          loksp=knr(ncon)
          call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp)
          addmol=zero
          do jk=1,nspel
             addmol(jk)=stoi(jk)
          enddo
          dqsum(ncon)=sites(ll)*qsp
          qsum=qsum+sites(ll)*qsp*yarr(ncon)
          do jk=1,nspel
!             write(*,*)'xmol: ',ncon,ik,jk,ielno(j)
             if(ielno(jk).gt.0) then
! ignore vacancies
! sumit is one above, can be removed here ....
                pmi%dxmol(ielno(jk),ncon)=sites(ll)*addmol(jk)/sumsit
                pmi%xmol(ielno(jk))=pmi%xmol(ielno(jk))+&
                     sites(ll)*addmol(jk)*yarr(ncon)/sumsit
!                pmi%nmol(ielno(jk))=pmi%nmol(ielno(jk))+&
!                     sites(ll)*addmol(jk)*yarr(ncon)
             endif
          enddo
       enddo constll
    enddo subll
! meqrec is not available in this subroutine
    do ik=1,nrel
       pmi%sumxmol=pmi%sumxmol+pmi%xmol(ik)
!       write(*,*)'sumwmol 3:',pmi%xmol(ik),mass_of(ik,ceq)
       pmi%sumwmol=pmi%sumwmol+pmi%xmol(ik)*mass_of(ik,ceq)
!       pmi%sumnmol=pmi%sumnmol+pmi%nmol(ik)
!       pmi%sumbmol=pmi%sumbmol+pmi%nmol(ik)*mass_of(ik,ceq)
    enddo
!    write(*,92)'onephase 3: ',pmi%iph,nsl,pmi%xdone,pmi%sumxmol,&
!         pmi%sumnmol,qq(1)
!    write(*,17)'Vacanies: ',qq
!       do i=1,noel()
!          write(*,17)'xm: ',pmi%xmol(i)
!          write(*,17)'dxm: ',(pmi%dxmol(i,j),j=1,ncon)
!       enddo
! now calculate G and all 1st and 2nd derivatives
! The calculated values are used also in other parts of the code 
    call calcg(iph,ics,2,lokcs,ceq)
    if(gx%bmperr.ne.0) then
       write(*,*)'Error calculating G'
       goto 1000
    endif
! calculate phase matrix elements, first second derivatives
    pmat=zero
    neq=ncon
    do ik=1,ncon
       do jk=ik,ncon
          pmat(ik,jk)=ceq%phase_varres(lokcs)%d2gval(ixsym(ik,jk),1)
! remove next line when using an inversion for symmetric matrix
          if(jk.gt.ik) pmat(jk,ik)=pmat(ik,jk)
       enddo
!       write(*,17)'row2A: ',(pmat(ik,jj),jj=1,nd1)
    enddo
! Then set the sublattice elements
    kk=0
    do ll=1,nsl
       do ik=1,nkl(ll)
! set the sublattice columns and rows
          kk=kk+1
          pmat(kk,neq+ll)=one
          pmat(neq+ll,kk)=one
       enddo
!       write(*,17)'row3: ',(pmat(ncon+ll,jj),jj=1,nd1)
    enddo
    neq=neq+nsl
    if(pmi%chargebal.eq.1) then
! if external charge balance add one column and one row (not yet tested)
       neq=neq+1
       do jk=1,ncon
! this is the row
          pmat(jk,neq)=dqsum(jk)
! this is the column
          pmat(neq,jk)=dqsum(jk)
       enddo
    endif
! debug output
!    write(*,*)'Phase matrix',nd1,neq,pmi%chargebal
!    do j=1,neq
!       write(*,17)'pmat: ',(pmat(i,j),i=1,neq)
!    enddo
! invert the phase matrix (faster subroutine should be used)
    call mdinv(nd1,nd2,pmat,pmi%invmat,neq,ierr)
    if(ierr.eq.0) then
!       write(*,*)'Phase matrix singular',pmi%iph,pmi%ics,pmi%ncc,ierr
!       do i=1,neq
!          write(*,17)'pmi:  ',(pmi%invmat(i,j),j=1,neq)
!       enddo
       gx%bmperr=4205; goto 1000
    endif
!    do i=1,neq
!       write(*,17)'pinv: ',(pmi%invmat(i,j),j=1,neq)
!    enddo
! maybe some common ending
    goto 900
!-------------------------------------------
900 continue
    goto 1000
!
1000 continue
    return
  end subroutine meq_onephase
! sumxmol=
 
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  logical function same_composition(jj,phr,meqrec,ceq,dgm)
! returns .TRUE. if phase phr(jj) has almost exactly the same composition
! as another composition set of the same phase that is stable
! dgm just for debug output
! =============================================================
! The composition of the phases are compared as ordered phases one can have
! the same constitution but distributed on different sets of sublattices ....
! ==============================================================
    implicit none
    integer jj
    double precision dgm
    TYPE(meq_data), dimension(*) :: phr
    TYPE(meq_setup) :: meqrec
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer jp,jy
! If the difference is larger than xdiff then the compositions are not the same
!    double precision, parameter :: xdiff=0.01D0
! FINETUNING: a large value of xdiff may mean you miss a miscibility gap
! a small value may create bad convergence
    double precision, parameter :: xdiff=0.05D0
    double precision, dimension(maxel) :: xmol1,xmol2,wmass
    double precision amount,totmol,totmass
! check if any phase with lower index than jj is the same phase
!    write(*,17)'testing if same_composition: ',phr(jj)%iph,phr(jj)%ics,dgm
17  format(a,i3,i2,5(1pe12.4))
    call calc_phase_molmass(phr(jj)%iph,phr(jj)%ics,xmol1,wmass,&
         totmol,totmass,amount,ceq)
    if(gx%bmperr.ne.0) goto 1000
    do jp=jj-1,1,-1
       if(phr(jp)%iph.eq.phr(jj)%iph) then
          if(phr(jp)%stable.eq.1) then
             call calc_phase_molmass(phr(jp)%iph,phr(jp)%ics,xmol2,wmass,&
                  totmol,totmass,amount,ceq)
             if(gx%bmperr.ne.0) goto 1000
!             write(*,118)phr(jj)%ics,dgm,(xmol1(j),j=1,noel())
!             write(*,118)phr(jp)%ics,zero,(xmol2(j),j=1,noel())
118          format('xtest: ',i2,1e12.4,8F9.6)
             do jy=1,meqrec%nrel
                if(abs(xmol1(jy)-xmol2(jy)).gt.xdiff) goto 110
             enddo
! we have found another stable composition set with same composition
             goto 300
          endif
       else
          exit
       endif
110    continue
    enddo
! check if any with higher index is same phase as jj
    do jp=jj+1,meqrec%nphase
       if(phr(jp)%iph.eq.phr(jj)%iph) then
          if(phr(jp)%stable.eq.1) then
             call calc_phase_molmass(phr(jp)%iph,phr(jp)%ics,xmol2,wmass,&
                  totmol,totmass,amount,ceq)
             if(gx%bmperr.ne.0) goto 1000
!             write(*,118)phr(jj)%ics,dgm,(xmol1(j),j=1,noel())
!             write(*,118)phr(jp)%ics,zero,(xmol2(j),j=1,noel())
             do jy=1,meqrec%nrel
                if(abs(xmol1(jy)-xmol2(jy)).gt.xdiff) goto 120
             enddo
! we have found another stable composition set with same composition
             goto 300
          endif
       else
          exit
       endif
120    continue
    enddo
! no composition set of same phase with same constitution
    same_composition=.FALSE.
    goto 1000
! we found a stable composition set with the same composition
300 continue
    same_composition=.TRUE.
    write(*,117)'Not added comp.set with same composition as stable phase: ',&
         phr(jj)%iph,phr(jj)%ics,phr(jp)%ics
117 format(a,i3,2i4)
! useless with two composition sets with same consitution, 
! try to reset this composition set to default constition
    call set_default_constitution(phr(jj)%iph,phr(jj)%ics,0,ceq)
    if(gx%bmperr.ne.0) goto 1000
! debug output
!    call get_phase_compset(phr(jj)%iph,phr(jj)%ics,lokph,lokcs)
!    if(gx%bmperr.ne.0) goto 1000
!    write(*,317)phr(jj)%iph,phr(jj)%ics,ceq%phase_varres(lokcs)%yfr
!317 format('New const: ',i3,i2/(11F7.4))
! The default constitution of this comp.set should be set if there is one
! Here just set this composition set dormant temporarily.
! It will be restored as entered at the end of the calculation.
!    phr(jj)%curd%status2=ibset(phr(jj)%curd%status2,CSSUS)
!    phr(jj)%curd%status2=ibset(phr(jj)%curd%status2,CSFIXDORM)
!    phr(jj)%dormlink=meqrec%dormlink
!    meqrec%dormlink=jj
    goto 1000
!
1000 continue
    return
  end function same_composition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine meq_deltaprop(iprop,pmi,jj,addon)
! Calculate the addition to a property in a condition or in a dot-derivative
! The property is not normallized (that will be handelled later) and
! if jj>0 that specifies the phase+comp.set.  If jj=0 the propery is
! summed over all stable phases.
    implicit none
    double precision, dimension(*) :: addon
    integer jj,iprop
    type(meq_data), pointer :: pmi
!\end{verbatim}
    integer ii
    double precision deltap
    write(*,*)'In meq_deltaprop: ',pmi%iph,pmi%ics,deltap
! multiply d2G/dPdy_i with deltaT, one value for each fraction
    do ii=1,pmi%ncc
       addon(ii)=pmi%curd%dgval(3,ii,1)*deltap
    enddo
1000 continue
    return
  end subroutine meq_deltaprop

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine calc_dyterms1(nrel,ia,tpindep,mamu,mag,mat,map,pmi)
! any change must also be made in
! subroutine calc_dyterms2(nrel,ia,tpindep,mamu,mag,mat,map,curd)
! calculate the terms in the deltay expression for amounts of component ia
! Equilibrium calculations can be speeded up if this routine calculates
! mamu as a symmeric array and mag, mat and map as an array for all 
! components.  But it was easier to write it like this initially to 
! understand what is going on.
!
! DM_A = \sum_B mu_B*MAMU(B) - MAG - MAT*dt - MAP*dp
!
! where MAMU=\sum_i dM_A/dy_i*\sum_j invmat(i,j)*dM_B/dy_j
!       c_iB=\sum_j invmat(i,j)*dM_B/dy_j etc etc
!
! it may not be very efficient but first get it right ....
! tpindep is TRUE if T or P are variable
    implicit none
    integer ia,nrel
    logical tpindep(2)
    double precision, dimension(*) :: mamu
    double precision mag,mat,map
! pmi is the phase data record for this phase
    type(meq_data), pointer :: pmi
!\end{verbatim} %+
! these are to be multiplied with mu(ib), nothing, deltaT, deltaP
    integer iy,jy,ib
    double precision sum,cig,cit,cip,cib
!
    do ib=1,nrel
       sum=zero
       do jy=1,pmi%ncc
          cib=zero
          do iy=1,pmi%ncc
             cib=cib+pmi%invmat(iy,jy)*pmi%dxmol(ib,iy)
          enddo
          sum=sum+cib*pmi%dxmol(ia,jy)
       enddo
       mamu(ib)=sum
    enddo
!-----------
    mag=zero
    mat=zero
    map=zero
    do jy=1,pmi%ncc
!       cig=zero
       cit=zero
       cip=zero
       do iy=1,pmi%ncc
!          write(*,11)jy,iy,mag,pmi%invmat(iy,jy),pmi%dxmol(ia,jy),&
!               pmi%curd%dgval(1,iy,1)
!11        format('cdy: ',2i3,6(1pe12.4))
          mag=mag+pmi%invmat(iy,jy)*pmi%dxmol(ia,jy)*pmi%curd%dgval(1,iy,1)
          if(tpindep(1)) cit=cit+pmi%invmat(iy,jy)*pmi%curd%dgval(2,iy,1)
          if(tpindep(2)) cip=cip+pmi%invmat(iy,jy)*pmi%curd%dgval(3,iy,1)
       enddo
!       mag=mag+pmi%dxmol(ia,jy)*cig
       mat=mat+pmi%dxmol(ia,jy)*cit
       map=map+pmi%dxmol(ia,jy)*cip
    enddo
1000 continue
    return
  end subroutine calc_dyterms1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
  subroutine calc_dyterms2(nrel,ia,tpindep,mamu,mag,mat,map,curd)
! any change must also be made in
! subroutine calc_dyterms1(nrel,ia,tpindep,mamu,mag,mat,map,pmi)
! Identical to calc_dyterms1 but data are taken from phase_varres record
!   call calc_dyterms2( .... ,phase_varres(lokcs))
! calculate the terms in the deltay expression for amounts of component ia
!
! DM_A = \sum_B mu_B*MAMU(B)-MAG-MAT*dt-MAP*dp
!
! where MAMU=\sum_i dM_A/dy_i*c_iB and 
!       c_iB=\sum_j invmat(i,j)*dM_B/dy_j etc etc
!
! it may not be very efficient but first get it right ....
! tpindep is TRUE if T or P are variable
    implicit none
    integer ia,nrel
    double precision, dimension(*) :: mamu
    double precision mag,mat,map
    logical tpindep(2)
    type(gtp_phase_varres), pointer :: curd
!\end{verbatim}
! pmi is the phase data record for this phase, replaced by curd ...
!    type(meq_data), pointer :: pmi
! these are to be multiplied with mu(ib), nothing, deltaT, deltaP
    integer iy,jy,ib,lokcs
    double precision sum,cig,cit,cip,cib
!
    do ib=1,nrel
       sum=zero
       do jy=1,curd%ncc
          cib=zero
          do iy=1,curd%ncc
             cib=cib+curd%cinvy(iy,jy)*curd%cdxmol(ib,iy)
          enddo
          sum=sum+cib*curd%cdxmol(ia,jy)
       enddo
       mamu(ib)=sum
    enddo
!-----------
    mag=zero
    mat=zero
    map=zero
    do jy=1,curd%ncc
       cig=zero
       cit=zero
       cip=zero
       do iy=1,curd%ncc
          cig=cig+curd%cinvy(iy,jy)*curd%dgval(1,iy,1)
          if(tpindep(1)) cit=cit+curd%cinvy(iy,jy)*curd%dgval(2,iy,1)
          if(tpindep(2)) cip=cip+curd%cinvy(iy,jy)*curd%dgval(3,iy,1)
       enddo
       mag=mag+curd%cdxmol(ia,jy)*cig
       mat=mat+curd%cdxmol(ia,jy)*cit
       map=map+curd%cdxmol(ia,jy)*cip
    enddo
1000 continue
    return
  end subroutine calc_dyterms2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

END MODULE matsmin
