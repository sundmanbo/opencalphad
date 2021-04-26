! Hillert's Minimizer as implemented by Sundman (HMS)
! Based on Mats Hillert paper in Physica 1981 and Bo Janssons thesis 1984
! Details of this implementation published in Computational Materials Science,
! vol 101, (2015) pp 127-137
!
! MODULE liboceq
!
MODULE liboceqplus
!
  use general_thermodynamic_package
  use minpack
!
! Copyright 2012-2021, Bo Sundman, France
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
! contact person: bo.sundman@gmail.com
!
!---------------------------
!
! To be implemented/improved
! - calculating dot derivatives (Cp, thermal expansion etc) PARTIALLY DONE
! - stability check (eigenvalues)
! - conditions for properties H, V, S etc. (partially done)
! - expressions as conditions (only for x(A) and N(A))
! - calculate gridminimizer after equilibrium as check DONE 
! - cleanup the use of chemical potentials. DONE
!
!
! For parallellization, also used in gtp3.F90
!$  use omp_lib
!
  implicit none
  character*8, parameter :: hmsversion='HMS-3.0'
!
!-------------------------------------------------------
! for single equilibrium
!
! BITS in meqrec status word
! MMQUIET means no output for the equilibrium calculation
! MMNOSTARTVAL means grid minimizer not called at start
  integer, parameter :: MMQUIET=0, MMNOSTARTVAL=1,MMSTEPINV=2
! NOTE in calceq7 status word is set to zero if more bits used because
! it seemed to have an arbitrary value and it created problems in macro map7
! I have now correceted the main reason (creating linehead records in SMP)
! but I kept this check
!
!\begin{verbatim}
  TYPE meq_phase
! parts of the data in this structure should be in the gtp_equilibrium_data
! it contains phase specific results from various subroutines during
! equilibrium calculation
! iph: phase number
! ics: composition set number
! idim: the dimension of phase matrix, 
! ncc: the number of constituents (same as idim??)
! stable: is 1 for a stable phase
! xdone: set to 1 for stoichiometric phases after calculating xmol first time
! dormlink: link to next phase that has temporarily been set dormant
! eec_check for equi-entropy check
! phtupix phase tuple index
     integer iph,ics,idim,stable,ncc,xdone,dormlink,eeccheck,phtupix
! value of phase status (-1,0=ent, 1=stable, 2=fix, -2=dorm, -3=sus, -4 hidden)
     integer phasestatus
! inverted phase matrix
     double precision, dimension(:,:), allocatable :: invmat
! mole fractions of components and their sum
     double precision, dimension(:), allocatable :: xmol
     double precision :: sumxmol,sumwmol
! Derivatives of moles of component wrt all constituent fractions of the phase
     double precision, dimension(:,:), allocatable :: dxmol
! link to phase_varres record
     TYPE(gtp_phase_varres), pointer :: curd
! value of amount and driving force at previous iteration
     double precision prevam, prevdg
! iteration when phase was added/removed
     integer itadd, itrem
! chargebal is 1 if external charge balance needed, ionliq<0 unless 
! ionic liquid when it is equal to nkl(1)=number of cations
     integer chargebal,ionliq,i2sly(2)
     double precision iliqcharge,yva
! end specific ionic liquids
  end TYPE meq_phase
!\end{verbatim}
!
!-------------------------------------------------------------------
!  
!\begin{verbatim}
  TYPE meq_setup
! one structure of this type is created when an equilibrium calculation
! is started and it holds all global data needed for handling the
! calculation of an equilibrium.  The phase specific data is in meq_phase
! nv: initial guess of number of stable phases
! nphase: total number of phases and composition sets
! nstph: current number of stable phases
! dormlink: is start of list of phases temporarily set dormant
! noofits current number of iterations
! status for various things
! nrel number of elements (components)
! typesofcond: types of conditions, =1 only massbal, =2 any conditions
! nfixmu number of fixed chemical potentials
! nfixph number of conditions representing fix phases
     integer nv,nphase,nstph,dormlink,noofits,status
     integer nrel,typesofcond,maxsph,nfixmu,nfixph
! component numbers of fixed potentials, reference and value 
     integer, dimension(:), allocatable :: mufixel
     integer, dimension(:), allocatable :: mufixref
! in this array the mu value as calculated from SER is stored
     double precision, dimension(:), allocatable :: mufixval
! in this array the mu value for user defined reference state is stored
     double precision, dimension(:), allocatable :: mufixvalref
! fix phases and amounts
     integer, dimension(:,:), allocatable :: fixph
     double precision, dimension(:), allocatable :: fixpham
! indices of axis conditions that has been inactivated
!     integer, dimension(:), allocatable :: inactiveaxis
! iphl, icsl: phase and composition sets of intial guess of stable phases
! aphl: initial guess of amount of each stable phase
     integer iphl(maxel+2),icsl(maxel+2)
     double precision aphl(maxel+2)
! stphl: current list of stable phases, value is index in phr array
     integer, dimension(maxel+2) :: stphl
! current values of chemical potentials stored in gtp_equilibrium_data
! if variable T and P these are TRUE, otherwise FALSE
     logical tpindep(2)
! these are the maximum allowed changes in T and P during iterations
     double precision tpmaxdelta(2)
! individual phase information
     type(meq_phase), dimension(:), allocatable :: phr
! this is used for EEC, pointer to liquid phr record and highest liquid entropy
     type(meq_phase), pointer :: pmiliq
     double precision seecliq
! information about conditions should be stored here.  Note that conditions
! may change during STEP and MAP
  end TYPE meq_setup
!\end{verbatim}
!
!------------------------------------------------------------------
!
! This is a connection to step/map
!\begin{verbatim}
  TYPE map_fixph
! provides information about phase sets for each line during mapping
     integer nfixph,nstabph,status
     type(gtp_phasetuple), dimension(:), allocatable :: fixph
     type(gtp_phasetuple), dimension(:), allocatable :: stableph
! most likely some of these variables are redundant stable_phr added 2020.03.05
     integer, dimension(:), allocatable :: stable_phr
     double precision, dimension(:), allocatable :: stablepham
! new 180814 to have nonzero fix phase amounts  ... not yet used
     double precision, dimension(:), allocatable :: fixphamap
  end TYPE map_fixph
!\end{verbatim}
! declared as mapfix in call to calceq7 and some other routines
!
! Added for debugging converge problems
  TYPE meqdebug
     integer mconverged,nvs,typ(10)
     integer :: flag=0
     double precision val(10),dif(10)
  end type meqdebug
  type(meqdebug) :: cerr
!
!\begin{verbatim}
! This is for returning the calculated value of an experimental property
! as we need an array to store the calculated values of the experimental  
! properties in order to calculate the Relative Standarad Deviation (RSD)
  double precision, allocatable, dimension(:) :: calcexp
! We cannot have EEC variabler here as it does bot work in parallel
! this is for EEC test
!  type(meq_phase), pointer :: pmiliq
! if several liquids check for largest S
!  type(meq_phase), pointer :: pmiliqsave
!  double precision eecliqentropy
! this is set TRUE when entering meq_onephase and false after one solid checked?
! it is now check for EEC
!  logical eecextrapol
! This is an (failed) attempt to limit Delta-T when having condition on y
  logical ycondTlimit
  double precision deltatycond
! TZERO and PARAEQUIL calculation need these (CANNOT BE USED IN PARALLEL)
  type(gtp_equilibrium_data), pointer :: tzceq
  type(gtp_condition), pointer :: tzcond
  type(gtp_state_variable), target :: musvr,xsvr
  integer tzph1,tzph2
! To prevent calculating a dot derivative at a given equilibrium
  integer :: special_circumstances=0
!\end{verbatim}
!
!--------------------------------------------------------------
!
! declared as part of phase_varres to be used in parallel
!  integer, dimension (:,:), allocatable :: phaseremoved
! debug output indicator
! mmdotder indicate dot derivative calculation, phase set may be different
! from the static memory
  integer :: mmdebug=0,mmdotder=0
! warning using B=value as condition
  logical bwarning
!--------------------------------------------------------------
!
! IMPORTANT
! phase_varres(lokcs)%amfu is the number of formula units of the phase
! phase_varres(lokcs)%netcharge is the current total charge  of the phase
! phase_varres(lokcs)%abnorm(1) is the number of real atoms per formula unit
! (may vary with composition like in (Fe,Cr,...)(Va,C,N,...) )
! phase_varres(lokcs)%abnorm(2) is the mass per formula unit
! NOTE: abnorm(1) and abnorm(2) are set by call to set_constitution
!
CONTAINS

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calceq2(mode,ceq)
!\begin{verbatim}
  subroutine calceq2(mode,ceq)
! calculates the equilibrium for the given set of conditions
! mode=0 means no global minimization
! ceq is a datastructure with all relevant thermodynamic data
    implicit none
    integer mode
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
    TYPE(meq_setup), allocatable, target :: meqrec1
    TYPE(meq_setup), pointer :: meqrec
    type(map_fixph), allocatable :: mapfix
!    type(map_fixph), pointer :: mapfix
    double precision starting,finish2,gtot
    integer starttid,endoftime,ij,addtuple,errall
    character name*16
!--------------------------------
    allocate(meqrec1,stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 1: ',errall
       gx%bmperr=4370; goto 1000
    endif
    meqrec=>meqrec1
    meqrec%status=0
    if(allocated(mapfix)) deallocate(mapfix)
    call cpu_time(starting)
    call system_clock(count=starttid)
! we may return here if gridcheck found a gridpoint below
100 continue
    call calceq7(mode,meqrec,mapfix,ceq)
    call system_clock(count=endoftime)
    call cpu_time(finish2)
1000 continue
    if(gx%bmperr.eq.0) then
! Gibbs energy using SER as reference state
       call get_state_var_value('GS ',gtot,name,ceq)
       if(gx%bmperr.ne.0) gx%bmperr=0
       if(.not.btest(globaldata%status,GSSILENT)) then
          if(ceq%eqno.ne.1) then
             write(*,1010)ceq%eqname(1:11),meqrec%noofits,&
                  finish2-starting,endoftime-starttid,gtot
          else
             write(*,1010)'Equilibrium',meqrec%noofits,&
                  finish2-starting,endoftime-starttid,gtot
          endif
1010   format(a,' result:',i4,' its, ',&
            1pe11.4,' s, ',i6,' cc, GS=',1pe15.7,' J/mol')
       endif
! Here we have now an equilibrium calculated.  Do a cleanup of the structure
! for phases with several compsets the call below shifts the stable one
! to the lowest compset number unless the default constitution fits another
! For example to ensure a fcc-carbonitrides is always the same compset.
       ij=1
! if meqrec%status indicate no initial startvalues set ij<0 to indicate test
! DO not test if mode=0
       if(mode.ne.0 .and. btest(meqrec%status,MMNOSTARTVAL)) ij=-ij
! OC went into a loop for a complex alloy calcumation here (once long ago ...)
!       write(*,*)'MM calling todo_after: 2',&
!            btest(meqrec%status,MMNOSTARTVAL),mode
       call todo_after_found_equilibrium(ij,addtuple,ceq)
       if(gx%bmperr.ne.0) then
          if(gx%bmperr.eq.4358) then
! gridpoint below current equilibrium found and set as stable (maybe new
! composition set).  Recalculate
             gx%bmperr=0
             write(*,*)'MM recalculating with this phase as stable 2: ',addtuple
             goto 100
          endif
       endif
!       write(*,*)'MM back in calceq2 after todo_after'
    endif
!CCI
! save the number of iterations needed to calculate the equilibrium
    ceq%conv_iter=meqrec%noofits
! maybe memory leak 2
!    write(*,*)'MM deallocate 2'
    deallocate(meqrec1)
!    write(*,*)'MM deallocated meqrec1'
    return
  end subroutine calceq2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calceq3
!\begin{verbatim} %-
  subroutine calceq3(mode,confirm,ceq)
! calculates the equilibrium for the given set of conditions
! mode=0 means no global minimization
! confirm is TRUE if output of CPU time
! ceq is a datastructure with all relevant thermodynamic data
    implicit none
    integer mode
    logical confirm
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    TYPE(meq_setup), allocatable, target :: meqrec1
    TYPE(meq_setup), pointer :: meqrec
    type(map_fixph), allocatable :: mapfix
!    type(map_fixph), pointer :: mapfix
    double precision starting,finish2
    integer starttid,endoftime,ij,addtuple,errall
!--------------------------------
    allocate(meqrec1,stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 2: ',errall
       gx%bmperr=4370; goto 1000
    endif
    meqrec=>meqrec1
    meqrec%status=0
    if(.not.confirm) meqrec%status=ibset(meqrec%status,MMQUIET)
    if(allocated(mapfix)) deallocate(mapfix)
!    nullify(mapfix)
    call cpu_time(starting)
    call system_clock(count=starttid)
! we may return here if gricheck found a new phase stable
100 continue
    call calceq7(mode,meqrec,mapfix,ceq)
    call system_clock(count=endoftime)
    call cpu_time(finish2)
1000 continue
    if(gx%bmperr.eq.0) then
! Here we have now an equilibrium calculated.  Do a cleanup of the structure
! for phases with several compsets the call below shifts the stable one
! to the lowest compset number unless the default constitution fits another
! For example to ensure a fcc-carbonitrides is always the same compset.
       ij=1
! if meqrec%status indicate no initial startvalues set ij<0 to indicate test
       if(mode.ne.0 .and. btest(meqrec%status,MMNOSTARTVAL)) ij=-ij
!       write(*,*)'MM Calling todo_after calceq3'
       call todo_after_found_equilibrium(ij,addtuple,ceq)
       if(gx%bmperr.eq.4358) then
! gridcheck after found a new phase stable!  recalculate
          gx%bmperr=0
!          write(*,*)'MM recalculate with new phase added as stable 3:',addtuple
          goto 100
       endif
       if(confirm) then
          write(*,1010)meqrec%noofits,finish2-starting,endoftime-starttid
1010      format('Equilibrium calculation ',i4,', its, ',&
               1pe12.4,' s and ',i7,' clockcycles')
       endif
    elseif(confirm) then
       write(*,1020)gx%bmperr
1020   format('Error return from equilibrium calculation ',i5)
    endif
! CCI save the number of iterations to calculate the equilibrium
    ceq%conv_iter=meqrec%noofits
! memory leak 2
!    write(*,*)'MM deallocate 3'
    deallocate(meqrec1)
!    write(*,*)'MM deallocated'
    return
  end subroutine calceq3

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calceq7
!\begin{verbatim}
  subroutine calceq7(mode,meqrec,mapfix,ceq)
! calculates the equilibrium for the given set of conditions
! mode=0 means no global minimization
! mode=-1 means used during step/map, no gridmin and do not deallocate phr
! ceq is a datastructure with all relevant thermodynamic data
! calling this routine instead of calceq2 makes it possible to extract
! additional information about the equilibrium from meqrec.
! Meqrec is also used for calculation of derivatives of state vatiables
    implicit none
    integer mode
    TYPE(meq_setup), pointer :: meqrec
    type(map_fixph), allocatable :: mapfix
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    TYPE(gtp_condition), pointer :: condition,lastcond
! conditions on T and P and mole fractions of components
    double precision, dimension(2) :: tpval
    double precision, dimension(maxel) :: xknown,vmu
! antot is total number of moles of atoms.  Needed to scale results from
! gridmin which assumes 1 mole of atoms
    double precision xxx,antot,cvalue,ccf(5)
    logical gridtest,formap
! for global minimization (change maybe to allocate dynamically)
    integer, dimension(maxph) :: nyphl
    double precision, dimension(maxconst) :: yarr
    integer np,iph,ics,jph,lokph,lokcs,mode2,errall
    integer mostcon,mph,nvf,mostconph(2,maxel),icc,jcc
! max number of potential conditions
    integer, parameter :: mmu=20
! dimension cmix(22) allows 5 terms: 2+4*5 
    integer mjj,ij,cmix(22),cmode,mufixel(mmu),mufixref(mmu),errout
    integer fixph(2,maxel),oldorder(mmu),kst,jj
! just for debugging
!    integer idum(1000)
    double precision fixpham(maxel),sumnp,props(5)
    logical ycond
    integer jq,ntup,saverr
!    character statevar*40
!
    ntup=nooftup()
!    write(*,*)'MM in calceq7',ntup
    ycond=.FALSE.
! this will be set to false when warning shown once for each calculation
    bwarning=.TRUE.
!    if(btest(meqrec%status,MMSTEPINV)) then
! this is the problem with map7? only bit 0 and 1 are used!!
!       write(*,'(a,z8)')'MM warning **** eqcalc7 meqrec%status: ',&
!            meqrec%status,' reset'
!       meqrec%status=0
!       meqrec%status may be set by STEPMAPPLOT to indicate 
!       write(*,*)'MM calceq7 meqrec%status: ',meqrec%status,meqrec%nstph
!    endif
    if(btest(globaldata%status,GSSILENT)) &
         meqrec%status=ibset(meqrec%status,MMQUIET)
    if(ocv()) write(*,*)"Entering calceq7",mode
    errout=0
! clear bit that start values has not been calculated
    meqrec%status=ibclr(meqrec%status,MMNOSTARTVAL)
    if(gx%bmperr.ne.0) then
       if(gx%bmperr.eq.4203 .or. gx%bmperr.eq.4204) then
! this means system matrix error and too many iterations respectivly
          write(kou,3)gx%bmperr
3         format('Error code ',i5,' reset before calling global minimizer')
          gx%bmperr=0
          errout=kou
       else
          write(kou,*)'Error code ',gx%bmperr,' prevents using global minimizer'
          goto 1000
       endif
    endif
    if(mode.ge.0) then
       mode2=mode
       formap=.FALSE.
    else
! formap .TRUE. means that phr will not be deallocated
! and that phr(jj)%phasestatus will be set from meqrec%fixph ....
       mode2=0
       formap=.TRUE.
    endif
! skip this if mode=-1, we may not have degrees of freedom equal to zero
! as the fix phase is not stored as condition ...
    if(mode.ge.0) then
!---------------------------
! extract conditions
       call extract_massbalcond(tpval,xknown,antot,ceq)
!       write(*,7)'MM xk: ',gx%bmperr,(xknown(mjj),mjj=1,noel())
7      format(a,i5,9(F8.4))
       if(gx%bmperr.ne.0) then
! error 4143 means no conditions, 4144 wrong number of conditions
          if(gx%bmperr.eq.4143 .or. gx%bmperr.eq.4144) then
!             write(*,*)'Degrees of freedom not zero',gx%bmperr
             goto 1000
          endif
! 4151 not only massbalance conditions
!       if(gx%bmperr.eq.4151) goto 1000
! these are other errors that makes it impossible to use gridminimizer
!          if(gx%bmperr.eq.4173 .or. &
!               gx%bmperr.eq.4174 .or. &
!               (gx%bmperr.ge.4176 .and. gx%bmperr.le.4185)) goto 1000
! if mode=0 we should not use grid minimizer
!          if(mode.ne.0 .or. .not.btest(meqrec%status,MMQUIET)) &
          if(mode.ne.0 .and. .not.btest(meqrec%status,MMQUIET)) &
               write(*,9)
9         format('Warning: global minimizer cannot be used for the current',&
               ' set of conditions')
          gx%bmperr=0
          gridtest=.true.
          meqrec%typesofcond=2
       else
!          meqrec%antot=antot
! no need for final grid minimizer as we will do one as start
          gridtest=.false.
          meqrec%typesofcond=1
       endif
       if(ocv()) write(*,*)'checked massbalance'
!------------------------------------
    endif
!    write(*,*)'In Calceq7 2'
    meqrec%nrel=noel()
! set some initial values
    meqrec%maxsph=noel()+2
    meqrec%nfixph=0
    meqrec%nfixmu=0
    meqrec%tpindep=.TRUE.
! limit change in T and P.  For P it should be a factor ...
    meqrec%tpmaxdelta(1)=2.0D2
    meqrec%tpmaxdelta(2)=1.0D2
! now we calculate maxsph, nfixmu and maybe other things for later
    lastcond=>ceq%lastcondition
    if(.not.associated(lastcond)) then
!       write(*,*)'No conditions'
       gx%bmperr=4143; goto 1000
    endif
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
! loop to investigate conditions, apply_condition:value in gtp3D.F90
70  continue
! comode=-1 means just check type of condition
! NOTE SPECIAL: condition on Y returns cmix(1)=6 to inhibit grid minimizer
       cmode=-1
       condition=>condition%next
       mjj=mjj+1
       if(ocv()) write(*,*)'check condition'
! a condtion can have several terms, ccf is coefficient for each term, 
! if just one term ccf (is assumed?) to be 1.0
       call apply_condition_value(condition,cmode,cvalue,cmix,ccf,ceq)
       if(gx%bmperr.ne.0) goto 1000
!       write(*,71)'MM apply 1: ',cmode,cvalue,cmix,ccf(1)
!71     format(a,i3,1pe14.4,10i4/12i4,1pe12.4)
!71     format(a,i3,1pe14.4,10i4/5(1pe12.4))
! cmix(1)=0 for inactive conditions
! cmix(1)=1 fix T, =2, fix P, =3 fix MU/AC/LNAC, =4 fix phase, =5 anything else
! if condition on T, P, potential or fix phase reduce maxsph
       select case(cmix(1))
       case default
          if(.not.associated(condition,lastcond)) goto 70
       case(1) ! fix T
          if(cvalue.le.1.0D-2) then
             write(*,*)'Condition on T must be larger than 0.01 K'
             gx%bmperr=4187; goto 1000
          endif
          meqrec%maxsph=meqrec%maxsph-1
          meqrec%tpindep(1)=.FALSE.
          ceq%tpval(1)=cvalue
       case(2) ! fix P
          if(cvalue.le.1.0D-2) then
             write(*,*)'Condition on P must be larger than 0.01 Pa'
             gx%bmperr=4187; goto 1000
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
             gx%bmperr=4189; goto 1000
          endif
          mufixel(np)=cmix(3)
          mufixref(np)=cmix(4)
! temporarily use yarr for something else
          if(cmix(2).eq.3) then
! Divide MU by RT
             yarr(np)=cvalue/ceq%rtn
          elseif(cmix(2).eq.4) then
! AC=exp(MU/RT) converted to chemical potential/RT
             if(cvalue.le.zero) then
                write(*,*)'Conditions on activity must be larger than zero'
                gx%bmperr=4191; goto 1000
             endif
             yarr(np)=LOG(cvalue)
          else
! LNAC=MU/RT which is the value used during minimization
             yarr(np)=cvalue
          endif
!          write(*,*)'Chemical potential condition: ',yarr(np)
          meqrec%maxsph=meqrec%maxsph-1
!             write(*,72)'MM, chemp: ',cmix(1),cmix(2),cmix(3),cvalue
!72           format(a,3i3,1pe12.4)
!-------------------------
       case(4) ! fix phase
! cmix(2) is phase index; cmix(2) is composition set
          meqrec%nfixph=meqrec%nfixph+1
          fixph(1,meqrec%nfixph)=cmix(2)
          fixph(2,meqrec%nfixph)=cmix(3)
          fixpham(meqrec%nfixph)=cvalue
!          write(*,*)'Fix phase condition: ',cmix(2),cmix(3),cvalue
! debug output of fix phase composition
!          call calc_phase_mol(cmix(1),yarr,ceq)
       case(5) ! mass balance condition
!          write(*,*)'MM cmix(1..4): ',cmix(1),cmix(2),cmix(3),cmix(4)
       case(6) ! Condition on Y, no grid minimizer
          ycond=.TRUE.
!          write(*,*)'MM condition on Y inhibit grid minimizer!'
       end select !-----------------------------------------------
       if(.not.associated(condition,lastcond)) goto 70
! end loop of conditions
!--------------------------------------------------------------
!       write(*,*)'variable potentials, max variable phases: ',&
!            noel()-cmix(2),meqrec%maxphases
    meqrec%nfixmu=np
    if(np.gt.0) then 
! number of fixed chemical potentials
       if(.not.allocated(meqrec%mufixel)) then
          allocate(meqrec%mufixel(np),stat=errall)
          allocate(meqrec%mufixref(np),stat=errall)
          allocate(meqrec%mufixval(np),stat=errall)
          allocate(meqrec%mufixvalref(np),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 3: ',errall
             gx%bmperr=4370; goto 1000
          endif
       else
! this can happen if activity condition and calculating without gridmin
          write(*,*)'Warning: meqrec has already mufixel allocated!'
       endif
       if(np.gt.1) then
! sort components with fix MU in increasing order to simplify below
          call sortin(mufixel,np,oldorder)
          do mjj=1,np
             nvf=mufixel(mjj)
             meqrec%mufixel(mjj)=nvf
             meqrec%mufixref(mjj)=mufixref(oldorder(mjj))
             meqrec%mufixval(mjj)=yarr(oldorder(mjj))
             meqrec%mufixvalref(mjj)=yarr(oldorder(mjj))
! copy fixed chemical potential (divided by RT) to ceq%cmuval also
             ceq%cmuval(nvf)=yarr(oldorder(mjj))
! in the component records multiply with RT
             ceq%complist(nvf)%chempot(1)=yarr(oldorder(mjj))*ceq%rtn
          enddo
       else
          nvf=mufixel(1)
          meqrec%mufixel(1)=nvf
          meqrec%mufixref(1)=mufixref(1)
          meqrec%mufixval(1)=yarr(1)
          meqrec%mufixvalref(1)=yarr(1)
! also copy fixed chemical potential to ceq%cmuval
          ceq%cmuval(nvf)=yarr(1)
          ceq%complist(nvf)%chempot(1)=ceq%cmuval(nvf)*ceq%rtn
       endif
    endif
    if(meqrec%nfixph.gt.0) then
! allocate 5 extra places for fix phase during mapping ...
       if(.not.allocated(meqrec%fixph)) then
!          write(*,*)'Allocate  meqrec%fixph'
          allocate(meqrec%fixph(2,meqrec%nfixph+5),stat=errall)
          allocate(meqrec%fixpham(meqrec%nfixph+5),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 4: ',errall
             gx%bmperr=4370; goto 1000
          endif
!          write(*,*)'Allocated meqrec%fixph'
       endif
       if(np.gt.1) then
! ?? sort phases in increasing order to simplify below
          write(*,*)'MM Cannot handle two fix phases ... '
          gx%bmperr=4192; goto 1000
       endif
       do mjj=1,meqrec%nfixph
          meqrec%fixph(1,mjj)=fixph(1,mjj)
          meqrec%fixph(2,mjj)=fixph(2,mjj)
          meqrec%fixpham(mjj)=fixpham(mjj)
       enddo
    else
! allocate 5 places for fix phase during mapping (one per axis)
       if(.not.allocated(meqrec%fixph)) then
          allocate(meqrec%fixph(2,5),stat=errall)
          allocate(meqrec%fixpham(5),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 5: ',errall
             gx%bmperr=4370; goto 1000
          endif
       endif
    endif
!----------------------------
!    call list_conditions(kou,ceq)
! skip if mode2=0 or global gridminimizer if bit set
!    write(*,*)'In Calceq7 4'
    if(mode2.eq.0 .or. btest(globaldata%status,GSNOGLOB)) then
! if errout set then grimin probably called to handel bad start point
      if(errout.eq.0) goto 110
!       write(*,*)'errout 2: ',errout
    endif
! skip global gridminimizer if only one component but make sure one phase
! has positive amount
    if(meqrec%nrel.eq.1) then
       goto 110
    endif
! skip global minimizer if ycond is true
    if(ycond) then
!       write(*,*)'MM condition on y(phase,const), no global minimizer'
       goto 110
    endif
!---------------------------------------------------------------
! Try global gridminimization.  Returned values are:
! nv is number of stable phase, iphl, icsl list of stable  phases, aphl amounts
! nyphl(j) is number of constituent fractions in phase j, yarr are the 
! constituent fractions, vmu the chemical potentials
! THIS CALL MAY CREATE NEW COMPOSITION SETS unless GSNOACS set.
! loop through all phases and set amount=0 and CSABLE off
    ij=1
    call todo_before(ij,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(meqrec%typesofcond.eq.1) then
! with only massbalance condition make a global grid minimization
!       call global_gridmin(1,tpval,xknown,meqrec%nv,&
!            meqrec%iphl,meqrec%icsl,meqrec%aphl,nyphl,yarr,vmu,idum,ceq)
!       write(*,*)'MM calling global gridmin'
       call global_gridmin(1,tpval,xknown,meqrec%nv,&
            meqrec%iphl,meqrec%icsl,meqrec%aphl,nyphl,vmu,ceq)
       if(ocv()) write(*,*)'MM back from gridmin'
!       write(*,*)'MM back from gridmin'
       if(gx%bmperr.ne.0) then
! if global fails reset error code and try a default start set of phases
!          if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
!             write(*,102)gx%bmperr,trim(bmperrmess(gx%bmperr))
!102          format('Error ',i5,': ',a/&
!                  'Minimizer tries using current or default start values')
!  write(kou,102)gx%bmperr,bmperrmess(gx%bmperr)
!             write(kou,102)bmperrmess(gx%bmperr)
!102          format(a/'Current constitution used as start values.')
!          else
!             write(kou,113)gx%bmperr 
!113          format('Cannot use grid minimazer, error: ',i5/&
!                  'Current constitution used as start values.')
!          endif
! no initial gridmin, make a gridtest at the end (not implemented ...)
!          else
!             write(*,*)'Grid minimizer cannot be used with these conditions'
!          endif
! set that grid minimizer is called after the equilibrium calculation
          gridtest=.true.
! problems using gridmin
! use current constitution or set default constitution (does not work well)
          gx%bmperr=0; goto 110
       endif
! multiply phase amounts with antot as global_grimin assumes 1 mole
       if(abs(antot-one).gt.1.0D-8) then
!          write(*,*)'From gridmin: ',meqrec%nv,antot
          do mph=1,meqrec%nv
             call get_phase_compset(meqrec%iphl(mph),meqrec%icsl(mph),&
                  lokph,lokcs)
             ceq%phase_varres(lokcs)%amfu=antot*ceq%phase_varres(lokcs)%amfu
          enddo
       endif
       if(ocv() .or. errout.gt.0) &
            write(*,103)(meqrec%iphl(mjj),meqrec%icsl(mjj),meqrec%aphl(mjj),&
            mjj=1,meqrec%nv)
103    format('Phases: ',12(i3,i2,F5.2))
       goto 200
    endif
!--------------------
! no global gridmin or we come here if gridtest finds a new stable phase
! UNFINISHED: A better start guess should be made!!!
!
110 continue
!    write(*,*)'starting without gridmin',errout
    meqrec%nv=0
! at least one phase must be stable
    mostcon=0
    mostconph=0
    mph=0
    jph=0
    sumnp=zero
    selph1: do iph=1,noph()
       selcs1: do ics=1,noofcs(iph)
          kst=test_phase_status(iph,ics,xxx,ceq)
          if(gx%bmperr.ne.0) goto 1000
! new: -4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fix
! skip loop selph1 for phases that are dormant or suspended
          if(kst.le.PHDORM) then
             if(ics.lt. noofcs(iph)) then
                cycle selcs1
             else
                cycle selph1
             endif
          endif
          call get_phase_compset(iph,ics,lokph,lokcs)
          if(ceq%phase_varres(lokcs)%amfu.gt.zero) then
             meqrec%nv=meqrec%nv+1
             meqrec%iphl(meqrec%nv)=iph
             meqrec%icsl(meqrec%nv)=ics
             meqrec%aphl(meqrec%nv)=ceq%phase_varres(lokcs)%amfu
             sumnp=sumnp+ceq%phase_varres(lokcs)%amfu
          endif
       enddo selcs1
! select the phases with most constituents
       call get_phase_variance(iph,nvf)
       if(mostcon.eq.0) then
          mostcon=mostcon+1
          mostconph(1,1)=nvf
          mostconph(2,1)=iph
       else
! very very clumsy
          do icc=1,mostcon
             if(nvf.le.mostconph(1,icc)) then
                if(icc.gt.1) then
! store this phase as a start phase if not in first position
! otherwise ignore it
                   if(mostcon.lt.noel()-meqrec%nfixmu) then
                      mostcon=mostcon+1
                      do jcc=icc+1,mostcon
                         mostconph(1,jcc)=mostconph(1,jcc-1)
                         mostconph(2,jcc)=mostconph(2,jcc-1)
                      enddo
                      mostconph(1,icc)=nvf
                      mostconph(2,icc)=iph
                   else
                      mostconph(1,icc-1)=nvf
                      mostconph(2,icc-1)=iph
                   endif
                endif
             endif
          enddo
       endif
    enddo selph1
    if(meqrec%nv.eq.0) then
! no phase with positive amount, set the noel()-meqrec%nfixmu-1 phases stable
! starting with those with highest number of constituents
       if(mostcon.eq.0) then
!          write(*,*)'MM no phase to set stable'
          gx%bmperr=4200; goto 1000
       endif
!       write(*,55)'Initial phases set stable: ',mostcon,&
!            (mostconph(1,icc),mostconph(2,icc),icc=1,mostcon)
!55     format(a,i3,10(2i3,2x))
       meqrec%nv=mostcon
!       write(*,56)(mostconph(1,icc),icc=1,mostcon)
!56     format('Setting start phases: ',20(i3))
       do icc=1,mostcon
          call get_phase_compset(mostconph(2,icc),1,lokph,lokcs)
!          ceq%phase_varres(lokcs)%amfu=one/mostcon
          ceq%phase_varres(lokcs)%amfu=one
          ceq%phase_varres(lokcs)%phstate=PHENTSTAB
          meqrec%iphl(icc)=mostconph(2,icc)
          meqrec%icsl(icc)=1
          meqrec%aphl(icc)=one
! this sets a default constitution 
          call set_default_constitution(mostconph(2,icc),1,ceq)
       enddo
    else
! hopefully set_constitution has been called ...
! normallize the sum of phase amounts assuming N=1 ... this did not help ...
!       if(sumnp.gt.one) then
!          sumnp=one/sumnp
!          do icc=1,meqrec%nv
!             meqrec%aphl(icc)=meqrec%aphl(icc)*sumnp
!          enddo
!       endif
!       write(*,57)(meqrec%iphl(icc),meqrec%icsl(icc),meqrec%aphl(icc),&
!            icc=1,meqrec%nv)
!57     format('Start phase set: ',10(i3,i2,F6.2))
       if(ocv()) write(*,*)'No global minimization, using current phase set',&
            meqrec%nv
    endif
! copy ceq%complist%chempot(1) to ceq%cmuval
    do mjj=1,meqrec%nrel
       if(abs(ceq%complist(mjj)%chempot(1)).ge.one) then
          ceq%cmuval(mjj)=ceq%complist(mjj)%chempot(1)/ceq%rtn
       else
          ceq%cmuval(mjj)=zero
       endif
    enddo
    if(ocv()) write(*,68)'MM cmuval: ',meqrec%nrel,&
         (ceq%cmuval(mjj),mjj=1,meqrec%nrel)
68  format(a,i3,6(1pe12.4))
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
! add this phase as stable, check that not too many stable phases ...
! meqrec%nv is the current number of stable phases
       if(meqrec%nv.eq.meqrec%maxsph) then
          write(*,*)'MM Too many stable phases'
          gx%bmperr=4193; goto 1000
       endif
!       write(*,*)'Adding fix phase to stable phase set',&
!            meqrec%fixph(1,mjj),meqrec%fixph(2,mjj)
       meqrec%nv=meqrec%nv+1
       meqrec%iphl(meqrec%nv)=meqrec%fixph(1,mjj)
       meqrec%icsl(meqrec%nv)=meqrec%fixph(2,mjj)
       meqrec%aphl(meqrec%nv)=meqrec%fixpham(mjj)
    enddo addfixph
!------------------------------- special for mapping and STEP
    mapfixdata: if(allocated(mapfix)) then
! for step only the status word is used to indicate an invarant node
!       if(mapfix%nfixph.eq.0) then
!          if(btest(mapfix,STEPINVARIANT)) then
!             exit mapfixdata
!          endif
!       endif
! the stable and fix phases copied from mapfix record.
       do ij=1,meqrec%nv
          meqrec%iphl(ij)=0
          meqrec%icsl(ij)=0
       enddo
       meqrec%nfixph=mapfix%nfixph
       meqrec%nv=0
       do ij=1,meqrec%nfixph
          meqrec%fixph(1,ij)=mapfix%fixph(ij)%ixphase
          meqrec%fixph(2,ij)=mapfix%fixph(ij)%compset
          meqrec%fixpham(ij)=zero
          if(allocated(mapfix%fixphamap)) then
! attempt 180814 to let fix phases have nonzero amount to improve mapping
             meqrec%fixpham(ij)=mapfix%fixphamap(ij)
             write(*,65)'MM fix mapphase: ',mapfix%fixph(ij)%ixphase,&
                  mapfix%fixph(ij)%compset,mapfix%fixphamap(ij)
65           format(a,2i5,1pe12.4)
!          else
!             write(*,65)'MM mapfix phase: ',mapfix%fixph(ij)%ixphase,&
!                  mapfix%fixph(ij)%compset
          endif
          meqrec%nv=meqrec%nv+1
          meqrec%iphl(meqrec%nv)=mapfix%fixph(ij)%ixphase
          meqrec%icsl(meqrec%nv)=mapfix%fixph(ij)%compset
! 180814 not sufficient to set aphl 
! because around line 1010 amfu is set to zero fix mapfix ... removed that!!
!          meqrec%aphl(meqrec%nv)=mapfix%fixpham(ij)
! I am not sure what value for mph  here
!          meqrec%phr(mph)%curd%amfu=zero
       enddo
       do ij=1,mapfix%nstabph
          meqrec%nv=meqrec%nv+1
          meqrec%iphl(meqrec%nv)=mapfix%stableph(ij)%ixphase
          meqrec%icsl(meqrec%nv)=mapfix%stableph(ij)%compset
          meqrec%aphl(meqrec%nv)=mapfix%stablepham(ij)
       enddo
!       write(*,64)'MM Stable mapphase: ',mapfix%nstabph,&
!            mapfix%stableph(1)%ixphase,mapfix%stableph(1)%compset,&
!            mapfix%stablepham(1)
64     format(a,i3,2i5,1pe12.4)
!    elseif(formap) then
! mapfixrecord not allocated for STEP calculations
! this dis not work for handling invariant nodes for STEP
!       write(*,*)'MM calceq7 formap MMSTEPINV:',btest(meqrec%status,MMSTEPINV)
!       if(btest(meqrec%status,MMSTEPINV)) then
! The line start at an invariant node for a STEP calculation, 
!          write(*,*)'MM invariant node with phases: ',meqrec%nstph
!          do jj=1,meqrec%nstph
!             jq=meqrec%stphl(jj)
!             write(*,*)'MM stable: ',jj,jq,meqrec%phr(jq)%curd%amfu
!          enddo
!       endif
    endif mapfixdata
!------------------------------- 
! zero start of link to phases set temporarily dormant ....
    meqrec%dormlink=0
!
!-------------------------------
! Now we calculate the equilibrium
200 continue
! allocate phaseremoved to avoid same phase stable again and again
!    write(*,*)'MM start interative minimizer',ceq%eqno
    if(allocated(ceq%phaseremoved)) deallocate(ceq%phaseremoved)
    ntup=nooftup()
    allocate(ceq%phaseremoved(2,ntup),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 6: ',errall
       gx%bmperr=4370; goto 1000
    endif
    ceq%phaseremoved=0
!
! this routine varies the set of phases and the phase constitutions
! until the stable set is found for the given set of conditions.
    if(ocv()) write(*,*)'MM calling meq_phaseset'
    call meq_phaseset(meqrec,formap,mapfix,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    gridtest=.false.
!------------------------------------------------------
!
! When we come here the equilibrium is calculated or calculation failed
! if failed or called from step/map (formap TRUE) just exit
    if(gx%bmperr.ne.0 .or. formap) goto 1000
!    write(*,*)'End of calceq7 ',gridtest
    if(gridtest) then
! gridtest value is set to .TRUE. if no gridmin done initially
       meqrec%status=ibset(meqrec%status,MMNOSTARTVAL)
    endif
!--------------------------------------------------
1000 continue
!    write(*,*)'MM back from meq_phaseset'
    if(gx%bmperr.ne.0) then
! test if total number of models > 10; that can create converge problems
       saverr=gx%bmperr; gx%bmperr=0
! This routine returns total G, S, V, N and B
       call sumprops(props,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Convergence error, check your conditions are reasonable'
       elseif(props(4).gt.1.0D1 .and. &
            .not.(saverr.eq.4210 .or. saverr.eq.4364)) then
          write(*,'(a,a,i5,1pe12.4)')'Convergence error, maybe reduce ',&
               'the size of your system!',saverr,props(4)
       endif
       gx%bmperr=saverr
    endif
! This error means T or P is less than 0.1
    if(gx%bmperr.eq.4187) write(*,*)'Exit calceq7 with error ',gx%bmperr
    return
  end subroutine calceq7

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine meq_phaseset
!\begin{verbatim}
  subroutine meq_phaseset(meqrec,formap,mapfix,ceq)
! this subroutine can change the set of stable phase and their amounts
! and constitutions until equilibrium is found for the current conditions.
    implicit none
    TYPE(meq_setup) :: meqrec
    type(map_fixph), allocatable :: mapfix
    TYPE(gtp_equilibrium_data), pointer :: ceq
    logical formap
!\end{verbatim}
! should one use meqrec as pointer here???
    integer ok,iadd,iph,ics,irem,jj,jph,kk,lastchange,lokph,lokcs,minadd
    integer kph,minrem,mph,nip,nochange,zap,toomanystable,jrem,krem,inmap
    double precision, parameter :: ylow=1.0D-3
    double precision, parameter :: addedphase_amount=1.0D-2
    double precision xxx,tpvalsave(2)
    integer iremsave,zz,tupadd,tuprem,samephase,phloopaddrem1,phloopaddrem2
! mapx is special for using meq_sameset for mapping
    integer phloopv,noremove,findtupix,saverr,mapx,errall
    character phnames*50
! prevent loop that a phase is added/removed more than 10 times
    integer, allocatable, dimension(:,:) :: addremloop
! replace always FALSE except when we must replace a phase as we have max stable
    logical replace,force
! number of iterations without adding or removing a phase
    replace=.FALSE.
    minadd=4
    minrem=4
    samephase=0
! minimum number of iterations between a change of the set of stable phases
! should not be smaller than minadd/minrem above ...
!    nochange=5
    nochange=4
! modification 180323/BoS
! allow removing a phase after 2 iterations
    noremove=2
    lastchange=0
!
    if(ocv()) write(*,*)'entering meq_phaseset: '
!    write(*,*)'MM entering meq_phaseset: '
    meqrec%dormlink=0
! nphase is set to total number of phases (phase+compset) to be calculated
! >>> parallellization ALERT, nphase may change when composition sets created
!    call sumofphcs(meqrec%nphase,ceq)
!    meqrec%nphase=totalphcs(ceq)
    meqrec%nphase=nonsusphcs(ceq)
    if(gx%bmperr.ne.0) goto 1000
! Nathalie had an error here "already allocated"
    if(allocated(meqrec%phr)) deallocate(meqrec%phr)
    allocate(meqrec%phr(meqrec%nphase),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 7: ',errall
       gx%bmperr=4370; goto 1000
    endif
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
!-----------------------------
    mph=0
    nip=1
!    krem=0
    meqrec%nstph=0
    do iph=1,noph()
       do ics=1,noofcs(iph)
! ignore hidden and suspended phases (also ignored above in sumofphcs)
! entered, fixed and dormat has values 1, 2 and 3, suspended 4, hidden 5
          zap=test_phase_status(iph,ics,xxx,ceq)
! new: -4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fixed
          if(zap.ge.PHDORM) then
             mph=mph+1
! this iph is the index in the phlista record
             meqrec%phr(mph)%iph=iph
             meqrec%phr(mph)%ics=ics
! compare with these the first time a phase wants to be added or removed
! if zero it means phase can be added/removed at iteration minadd/minrem
             meqrec%phr(mph)%itadd=0
             meqrec%phr(mph)%itrem=0
! initiate indicator for phases with fix composition, set to 1 later if so
             meqrec%phr(mph)%xdone=0
! save phasestatus, zap>-2 here so set all -1,0,1 set to 0
             if(abs(zap).le.1) zap=0
             meqrec%phr(mph)%phasestatus=zap
! set link to calculated values of G etc.
             call get_phase_compset(iph,ics,lokph,lokcs)
             meqrec%phr(mph)%curd=>ceq%phase_varres(lokcs)
! ALSO save phase tuple index
             findtupix=meqrec%phr(mph)%curd%phtupx
             meqrec%phr(mph)%phtupix=findtupix
!             write(*,'(a,4i6,5x,3i6)')'MM save tuple index: ',mph,iph,ics,&
!                  findtupix,phasetuple(findtupix)%ixphase,&
!                  phasetuple(findtupix)%compset,phasetuple(findtupix)%lokph
! set number of constituents, DO NOT USE size(...curd%size(yfr)!!!
             meqrec%phr(mph)%ncc=noconst(iph,ics,ceq)
             if(formap) then
! when mapping fix phases are used to replace axis conditions.  The
! fix phases are in the meqrec%fixph array
! They do not return PHFIXED for test_phase_status !!!
                do zz=1,meqrec%nfixph
                   if(iph.eq.meqrec%fixph(1,zz) .and. &
                        ics.eq.meqrec%fixph(2,zz)) then
                      meqrec%phr(mph)%phasestatus=PHFIXED
                      if(allocated(mapfix)) then
                         if(allocated(mapfix%fixphamap)) then
                            meqrec%phr(mph)%curd%amfu=mapfix%fixphamap(1)
                            write(*,*)'MM set fixamount: ',&
                                 mapfix%fixphamap(1)
                         endif
                      endif
                   endif
                enddo
! inmap=1 turns off converge control of T
                inmap=1
             else
! inmap=0 means not called from step/map routines
                inmap=0
             endif
             meqrec%phr(mph)%ionliq=-1
             meqrec%phr(mph)%i2sly=0
             if(test_phase_status_bit(iph,PHIONLIQ)) meqrec%phr(mph)%ionliq=1
! already done: set link to calculated values of G etc. 
!             call get_phase_compset(iph,ics,lokph,lokcs)
!             meqrec%phr(mph)%curd=>ceq%phase_varres(lokcs)
! causing trouble at line 3175 ???
             if(nip.le.meqrec%nv) then
                if(iph.eq.meqrec%iphl(nip) .and. ics.eq.meqrec%icsl(nip)) then
! this phase is part of the initial stable set, increment nstph
                   meqrec%nstph=meqrec%nstph+1
                   meqrec%stphl(meqrec%nstph)=mph
                   meqrec%phr(mph)%stable=1
                   if(meqrec%phr(mph)%phasestatus.eq.PHFIXED) then
! Rather confused here ...
! fixed phases as conditions have an amount in meqrec%fixpham
! fixed phases during mapping should have zero amount (maybe not ...)
!                   krem=krem+1
!                   write(*,*)'MM aphl for fix phase: ',krem,mph,&
!                        meqrec%fixpham(krem)
                      if(meqrec%phr(mph)%curd%phstate.ne.PHFIXED) then
! this is a phase set fix by mapping, set amount to zero unless mapfix%fixpham 
! but mapfix is not available in this routine ..
                         if(allocated(mapfix%fixphamap)) then
! 180814 tried to remove setting fix phase amount to zero
                            write(*,*)'MM nonzero mapfix amount !'
                            meqrec%phr(mph)%curd%amfu=mapfix%fixphamap(1)
                         else
                            meqrec%phr(mph)%curd%amfu=zero
                         endif
                      endif
                   else
! this is setting non-zero fixed amount of a phase as condition
! Trying to handle this in mapping ... but here it not the fix phase ...
                      if(allocated(mapfix)) then
                         if(allocated(mapfix%fixphamap)) &
                              write(*,*)'MM phase amount: ',&
                              meqrec%phr(mph)%iph,meqrec%aphl(meqrec%nstph)
                      endif
                      meqrec%phr(mph)%curd%amfu=meqrec%aphl(meqrec%nstph)
                   endif
! set "previous values"
                   meqrec%phr(mph)%prevam=meqrec%aphl(meqrec%nstph)
                   meqrec%phr(mph)%prevdg=zero
                   nip=nip+1
                else
! unstable phase
                   meqrec%phr(mph)%stable=0
                   meqrec%phr(mph)%prevam=zero
                   meqrec%phr(mph)%prevdg=-one
                   meqrec%phr(mph)%curd%amfu=zero
                endif
             else
!                write(*,312)'MM nip: ',nip,meqrec%nv
!312             format(a,5i4)
! unstable phase
                meqrec%phr(mph)%stable=0
                meqrec%phr(mph)%prevam=zero
                meqrec%phr(mph)%prevdg=-one
                meqrec%phr(mph)%curd%amfu=zero
             endif
! mark that no data arrays allocated for this phase
             meqrec%phr(mph)%idim=0
! initiate link to another phase temporarily set dormant zero
             meqrec%phr(mph)%dormlink=0
          else
! we are here for phases that are suspended, test_phase_status return -3
! make sure stable bit is cleared in phases not included in calculation
! maybe the whole status word should be zeroed?
             call get_phase_compset(iph,ics,lokph,lokcs)
             ceq%phase_varres(lokcs)%status2=&
                  ibclr(ceq%phase_varres(lokcs)%status2,CSABLE)
! check if suspended phase bits CSSUS set
!z             if(btest(ceq%phase_varres(lokcs)%status2,CSSUS)) then
!                write(*,*)'MM Suspended bit set',lokph,lokcs
!             else
! This should not be necessary but it fixes the problem using c n with
! suspended phases.  The CSSUS bit should no longer be used???
!                write(*,*)'MM warning, suspended bit NOT set',lokph,lokcs
!z                ceq%phase_varres(lokcs)%status2=&
!z                     ibset(ceq%phase_varres(lokcs)%status2,CSSUS)
!z             endif
          endif
       enddo
    enddo
! problem phases suspended are restored!!
!    write(*,*)'MM at start, nonsuspenden phases: ',mph
    meqrec%noofits=0
    toomanystable=0
    jrem=0
    krem=0
    iremsave=0
    phloopaddrem1=0
! code above executed only intially
!    write(*,*)'MM allocating addremloop',meqrec%nphase
    allocate(addremloop(meqrec%nphase,3),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 8: ',errall
       gx%bmperr=4370; goto 1000
    endif
    addremloop=0
!----------------------------------------------------------------
!
! meq_sameset calculate the equilibrium for a given set of stable phases
! if the phase set change we return to this routine to take some action and
! then call meq_sameset again
! irem nonzero if phase irem should be removed
! iadd nonzero if phase iadd should be added
! meqrec has the general information needed
! meqrec%phr is the array with phases
! ceq is the connection to the model package data
200 continue
!    iadd=-1 ! iadd =-1 turns on verbose in meq_sameset
    iadd=0
    irem=iremsave
! for debuging convergence
!    call list_stable_phases('MM call:',meqrec%noofits,iadd,irem,meqrec,ceq)
!    write(*,*)'MM calling meq_sameset ',meqrec%noofits
!    write(*,*)'MM calling list conditions'
!    call list_conditions(kou,ceq)
! meq_sameset varies amounts of stable phases and constitutions of all phases
! If there is a phase change (iadd or irem nonzeri) or error it exits 
! mapx is needed when using meq_sameset for mapping, irrelevant here
    mapx=0
    call meq_sameset(irem,iadd,mapx,meqrec,meqrec%phr,inmap,ceq)
    if(ocv()) write(*,*)'MM back from sameset ',irem,iadd,meqrec%noofits
!    write(*,*)'MM back from meq_sameset ',irem,iadd,meqrec%noofits
    if(gx%bmperr.ne.0) then
       if(gx%bmperr.eq.4364) then
!          write(*,*)'MM Two phases with same stoichiometry stable, to be fixed'
       endif
       goto 1000
    endif
!
    force=.false.
!    call list_stable_phases('MM back:',meqrec%noofits,iadd,irem,meqrec,ceq)
!    write(*,*)'MM line 1114:',irem,iadd
    if(irem.gt.0 .or. iadd.gt.0) then
       if(iremsave.gt.0 .and. iadd.eq.iremsave) then
! if iadd=iremsave>0 there was a equil matrix error when removing iremsave
          irem=0
          force=.true.
       elseif(meqrec%noofits-lastchange.lt.nochange) then
!          write(*,221)' *** Phase set change not allowed: ',&
!               meqrec%noofits,lastchange,nochange,irem,iadd
!221       format(a,10i4)
          goto 200
       endif
! keep record of adding and removing phases
       if(iadd.gt.0) then
          addremloop(iadd,1)=meqrec%noofits
          if(irem.eq.0) then
             addremloop(iadd,2)=addremloop(iadd,2)+1
!             write(*,'(a,4i5)')'MM adding:   ',addremloop(iadd,1),iadd,&
!                  addremloop(iadd,2),addremloop(iadd,3)
          endif
          if(addremloop(iadd,2).gt.5) then
             if(.not.btest(meqrec%status,MMQUIET)) &
                  write(*,'(a,2i4,"#",i1)')'MM Removing phase: ',iadd,&
                  meqrec%phr(iadd)%iph,meqrec%phr(iadd)%ics
             meqrec%phr(iadd)%phasestatus=PHDORM
             meqrec%phr(iadd)%curd%phstate=PHDORM
             meqrec%phr(iadd)%dormlink=meqrec%dormlink
             meqrec%dormlink=iadd
! iremsave keeps track of last removed phase, if equal to iadd set it to 0
             if(iremsave.eq.iadd) iremsave=0
             iadd=0
             goto 200
          endif
       else
          addremloop(irem,3)=addremloop(irem,3)+1
!             write(*,'(a,3i5)')'MM removing: ',addremloop(irem,1),irem,&
!                  addremloop(iadd,2),addremloop(iadd,3)
!          if(addremloop(irem,3).gt.5) then
!             write(*,'(a,3i5)')'MM Suspend ',addremloop(irem,1),irem,&
!                  meqrec%dormlink
!             meqrec%phr(irem)%phasestatus=PHDORM
!             meqrec%phr(irem)%curd%phstate=PHDORM
!             meqrec%phr(irem)%dormlink=meqrec%dormlink
!             meqrec%dormlink=irem
!             irem=0
!             goto 200
!          endif
       endif
       if(iadd.gt.0) then
! check if phase to be added is already stable as another composition set
! This check should maybe be above as maybe another phase want to be stable??
! The last argument is not used
          if(same_composition(iadd,meqrec%phr,meqrec,ceq,zero)) then
!             write(*,*)'MM ignoring the same phase twice: ',iadd
             goto 200
          endif
! do not add phases with net charge
          if(meqrec%phr(iadd)%curd%netcharge.gt.1.0D-2) then
             if(iadd.ne.samephase) then
                write(*,218)'MM ignoring phase with net charge: ',iadd
!                meqrec%phr(iadd)%curd%phtupx,meqrec%phr(iadd)%curd%netcharge
218             format(a,2i5,1pe14.6)
                samephase=iadd
             endif
             goto 200
          elseif(phloopaddrem1.gt.4) then
! reset this phase to a default constitution
             if(.not.btest(meqrec%status,MMQUIET)) &
                  write(*,*)'MM phloopaddrem: ',phloopaddrem2
             iadd=phloopaddrem2
             phloopv=phasetuple(iadd)%lokph
!             if(ceq%phlista(phloopv)%tnooffr-ceq%phlista(phloopv)%noofsubl &
!                  .gt. 0) then
! reset troublesome phase constitution if it can vary
                call set_default_constitution(phasetuple(iadd)%ixphase,&
                     phasetuple(iadd)%compset,ceq)
!             else
! set phase dormant ... Hm I do not understand meqrec%phr any longer ...
!                phloopv=phasetuple(iadd)%lokvares
!                ceq%phase_varres(phloopv)%phstate=PHDORM
!             endif
             iadd=0
             phloopaddrem1=0
             phloopaddrem2=0
             goto 200
!          elseif(meqrec%phr(iadd)%curd%netcharge.gt.1.0D-8) then
!             write(*,231)'MM adding phase with net charge: ',iadd,&
!                  meqrec%phr(iadd)%curd%phtupx,meqrec%phr(iadd)%curd%netcharge
!231          format(a,2i5,1pe14.6)
          endif
       endif
       tupadd=0
       tuprem=0
       xxx=0.0D0
!       if(iadd.gt.0) tupadd=meqrec%phr(iadd)%curd%phtupx
!       if(irem.gt.0) tuprem=meqrec%phr(irem)%curd%phtupx
       if(iadd.gt.0) tupadd=meqrec%phr(iadd)%phtupix
       if(irem.gt.0) tuprem=meqrec%phr(irem)%phtupix
       if(.not.btest(meqrec%status,MMQUIET)) then
          if(iadd.gt.0) then
             phnames='+'
             call get_phasetup_name(tupadd,phnames(2:))
             if(irem.gt.0) then
                kk=len_trim(phnames)+3
                phnames(kk-1:kk-1)='-'
                call get_phasetup_name(tuprem,phnames(kk:))
             endif
          else
             phnames='-'
             call get_phasetup_name(tuprem,phnames(2:))
          endif
          addph: if(formap) then
!             if(btest(meqrec%status,MMSTEPINV)) then
! This did not work to handle invariants during STEP
! we are exiting an invariant node for a STEP calculation, allow phase change
! meq_sameset wants to ADD a phase, instead remove the last stable phase
!                write(*,*)'MM meq_phaseset invariant node',meqrec%noofits,iadd
!                do jj=1,meqrec%nstph
!                   irem=meqrec%stphl(jj)
!                   if(iadd.eq.0 .and. &
!                        meqrec%phr(irem)%curd%amfu.eq.zero) then
!                      meqrec%phr(irem)%curd%amfu=1.0D-1
!                   endif
!                   write(*,*)'MM stable: ',jj,irem,meqrec%phr(irem)%curd%amfu
!                enddo
!                if(iadd.gt.0 .and. meqrec%nstph.gt.1) then
!                   meqrec%nstph=meqrec%nstph-1
!                   meqrec%phr(irem)%curd%amfu=zero
!                   write(*,*)'MM ignore adding ',iadd,' but remove ',irem
!                   iadd=0
!                   goto 200
!                endif
!                exit addph
!             endif
! This can be too strong, we can have a tie-line betwen two stoichiometric
! phases, i.e. a new phase appears at first attempt to step in two-phase region.
! UNFINISHED handling of many exceptions during mapping 
             write(*,'(a,a)')'MM Phase change not allowed: ',trim(phnames)
             gx%bmperr=4210; goto 1000
#ifdef silent
#else
          elseif(ceq%eqno.ne.1) then
!             write(*,219)meqrec%noofits,iadd,irem,' at equil: ',ceq%eqno
!219          format('Phase change: its/add/remove: ',3i5,a,i5)
             if(.not.btest(meqrec%status,MMQUIET)) &
                  write(*,219)ceq%eqno,meqrec%noofits,trim(phnames)
219          format('Phase change (equil: ',i3,') iteration: ',i5,', phase: ',a)
#endif
          else
             if(iadd.gt.0) then
                phnames='+'
                call get_phasetup_name(tupadd,phnames(2:))
                if(irem.gt.0) then
                   kk=len_trim(phnames)+3
                   phnames(kk-1:kk-1)='-'
                   call get_phasetup_name(tuprem,phnames(kk:))
                endif
             else
                phnames='-'
                call get_phasetup_name(tuprem,phnames(2:))
             endif
#ifdef silent
#else             
             if(.not.btest(meqrec%status,MMQUIET)) &
                  write(*,281)meqrec%noofits,trim(phnames)
281          format('Phase change iteration: ',i5,2x,a)
#endif
          endif addph
       endif
    endif
222 continue
    remove: if(irem.gt.0) then
! remove a phase ---------------------------
       if(ocv()) write(*,223)'Phase to be removed: ',meqrec%phr(irem)%iph,&
            meqrec%phr(irem)%ics,meqrec%phr(irem)%curd%amfu,meqrec%noofits
       if(meqrec%nstph.eq.1) then
          if(.not.REPLACE) then
! we must be able to REPLACE the only stable phase for a unary system
             write(*,*)'Attempt to remove the only stable phase!!!'
             goto 200
          endif
!          write(*,*)'MM replacing one stable phase with another',irem,iadd
       else
! make sure replace is false unless explitly set below
          replace=.FALSE.
       endif
       if(meqrec%noofits-meqrec%phr(irem)%itadd.lt.minrem) then
! if phase was just added do not remove before minrem iterations
          if(ocv()) write(*,*)'Too soon to remove phase',&
               meqrec%phr(irem)%curd%phtupx,meqrec%noofits,&
               meqrec%phr(irem)%itadd
          if(phloopaddrem1.gt.0) then
             if(phloopaddrem2.eq.meqrec%phr(irem)%curd%phtupx) then
                phloopaddrem1=phloopaddrem1+1
             else
                phloopaddrem2=0
                phloopaddrem1=0
             endif
          else
             phloopaddrem2=meqrec%phr(irem)%curd%phtupx
             phloopaddrem1=1
          endif
          goto 200
       endif
! shift phases after irem down in meqrec%stphl
! irem is index to meqrec%phr(), meqrec%stphl(jph) is index to meqrec%phr
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
       meqrec%phr(irem)%itrem=meqrec%noofits
       meqrec%phr(irem)%prevam=zero
       meqrec%phr(irem)%stable=0
       meqrec%phr(irem)%curd%amfu=zero
! save irem as it is used to restore a phase if massbalance problem
       iremsave=irem
       irem=0
       lastchange=meqrec%noofits
! one can remove and add a phase at the same time !!!
       if(iadd.eq.0) then
          toomanystable=0
          jrem=0
          goto 200
       endif
    endif remove
!------------------------------------------- 
    add: if(iadd.gt.0) then
! add a phase.  This can be tricky
! NOTE it must be added so meqrec%stphl in ascending order
       if(ocv()) write(*,223)'Phase to be added:   ',meqrec%phr(iadd)%iph,&
            meqrec%phr(iadd)%ics,meqrec%phr(iadd)%curd%dgm,meqrec%noofits
223    format(a,2x,2i4,1pe15.4,i7)
       if(meqrec%noofits-meqrec%phr(iadd)%itrem.lt.minadd .and. .not.force) then
! if phase was just removed, do not add it before minadd iterations
!          if(.not.btest(meqrec%status,MMQUIET))write(*,224)
          if(ocv()) write(*,224)meqrec%phr(iadd)%curd%phtupx,&
               meqrec%noofits,meqrec%phr(iadd)%itrem,phloopaddrem1,&
               phloopaddrem2,minadd
224       format('Too soon to add phase: ',i3,2x,i4,2x,5i5)
          if(phloopaddrem1.gt.0) then
             if(phloopaddrem2.eq.meqrec%phr(iadd)%curd%phtupx) then
                phloopaddrem1=phloopaddrem1+1
             else
                phloopaddrem2=0
                phloopaddrem1=0
             endif
          else
             phloopaddrem2=meqrec%phr(iadd)%curd%phtupx
             phloopaddrem1=1
          endif
          goto 200
       endif
!       if(iadd.eq.abs(iremsave)) then
!          write(*,*)'Phase just removed, do not add: ',iadd
!          iremsave=0
!          goto 200
!       endif
! make sure iremsave is zero
       iremsave=0
       if(meqrec%nstph.eq.meqrec%maxsph) then
! No more phases allowed, we must see if  some other phase may be removed
          if(toomanystable.ge.3) then
!             write(*,*)'Attempt to set too many phases stable',meqrec%maxsph
!             gx%bmperr=4201; goto 1000
! During mapping do not replace phases ...
             if(formap) then
                gx%bmperr=4201; goto 1000
             endif
! UNFINISHED code below
             if(jrem.eq.0) then
! try to remove a stable phase ... which? Replace the one that does not
! disturb the order of phases in meqrec%stphl by adding iadd
                do iph=1,meqrec%nstph
                   if(iadd.gt.meqrec%stphl(iph)) cycle
                   jrem=meqrec%stphl(iph); exit
                enddo
! if jrem zero here replace the last
                if(jrem.eq.0) jrem=meqrec%stphl(meqrec%nstph)
                krem=jrem
                irem=jrem
                if(.not.btest(meqrec%status,MMQUIET)) &
                     write(*,241)meqrec%noofits,irem,iadd,ceq%tpval(1)
241             format('MM Too many stable phases at iter ',i3,', phase ',i3,&
                     ' replaced by ',i3,', T= ',F8.2)
!                write(*,240)meqrec%noofits,irem,iadd,ceq%tpval(1),&
!                     (meqrec%stphl(iph),iph=1,meqrec%nstph)
!240             format('Too many stable phases at iter ',i3,', phase ',i3,&
!                     ' replaced by ',i3,', T= ',F8.2/3x,15(i3))
                replace=.TRUE.
                goto 222             
             else
                write(*,*)'MM setting too many phases stable',meqrec%maxsph
                gx%bmperr=4201; goto 1000
             endif
          else
! try ignore adding 3 times
!             write(*,*)'Ignoring attempt to set too many phases stable',&
!                  meqrec%maxsph,toomanystable
             toomanystable=toomanystable+1
             goto 200
          endif
       endif
! the phase must be added in sequential order of phase and composition set no
       findplace: do jph=1,meqrec%nstph
          jj=meqrec%stphl(jph)
          if(meqrec%phr(iadd)%iph.gt.meqrec%phr(jj)%iph) then
             cycle
          endif
          if(meqrec%phr(iadd)%iph.lt.meqrec%phr(jj)%iph) then
             exit
          endif
! if same phase number compare composition set numbers
          if(meqrec%phr(iadd)%iph.eq.meqrec%phr(jj)%iph) then
             if(meqrec%phr(iadd)%ics.gt.meqrec%phr(jj)%ics) then
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
!       write(*,*)'Phase added: ',jph,meqrec%nstph,meqrec%maxsph
! phase added at jph, (note jph may be equal to nstph+1)
       meqrec%stphl(jph)=iadd
       meqrec%nstph=meqrec%nstph+1
       meqrec%phr(iadd)%itadd=meqrec%noofits
       meqrec%phr(iadd)%curd%dgm=zero
       lastchange=meqrec%noofits
! maybe some more variables should be set?
       meqrec%phr(iadd)%curd%amfu=addedphase_amount
       meqrec%phr(iadd)%stable=1
       iadd=0
       toomanystable=0
       jrem=0
       goto 200
    endif add
!---------------------------------------------------
! found stable phase set or error
1000 continue
    if(gx%bmperr.eq.0) then
! equilibrium calculation converged, one should add check on stability
!
! >> add calculate eigenvalues of phase matrix to check stability, 
! >> a negative eigenvalue means inside spinodal
! >> Note charge problems for metastable phases, phase must be neutral ...
!
!------------------------------------------------------------
! clear bits: no equilibrium calculated/ inconsistent conditions and result/
! equilibrium calculation failed/ only gridcal
       ceq%status=ibclr(ceq%status,EQNOEQCAL)
       ceq%status=ibclr(ceq%status,EQINCON)
       ceq%status=ibclr(ceq%status,EQFAIL)
       ceq%status=ibclr(ceq%status,EQGRIDCAL)
! set stable bit in stable phases and clear it in all others
       kk=1
       do jj=1,mph
          if(jj.eq.meqrec%stphl(kk)) then
             meqrec%phr(jj)%curd%status2=&
                  ibset(meqrec%phr(jj)%curd%status2,CSABLE)
! the stable phase list should be ordered in increasing phase number
             kk=min(kk+1,meqrec%nstph)
!             write(*,*)'mm max kk: ',kk,meqrec%nstph
          else
             meqrec%phr(jj)%curd%status2=&
                  ibclr(meqrec%phr(jj)%curd%status2,CSABLE)
          endif
       enddo
!-----------------------
! loop through all phases and if their status is entered set it as PHENTUNST
! unless stablestable phases and set the PHENTST for phases in stable set
! That is important for extracting values later ...
       do jph=1,meqrec%nphase
          if(meqrec%phr(jph)%curd%phstate.ge.PHENTUNST .and. &
               meqrec%phr(jph)%curd%phstate.le.PHENTSTAB) then
             meqrec%phr(jph)%curd%phstate=PHENTUNST
          endif
       enddo
       do jph=1,meqrec%nstph
          jj=meqrec%stphl(jph)
          if(meqrec%phr(jj)%curd%phstate.lt.PHFIXED) then
             meqrec%phr(jj)%curd%phstate=PHENTSTAB
          endif
       enddo
!-----------------------------------------
    else
!       write(*,*)'MM cleaning up due to error'
! set some failure bits
       ceq%status=ibset(ceq%status,EQINCON)
       ceq%status=ibset(ceq%status,EQFAIL)
       ceq%status=ibclr(ceq%status,EQGRIDCAL)
! even when not converged save the current chemical potentials
       do jj=1,meqrec%nrel
          ceq%complist(jj)%chempot(1)=ceq%cmuval(jj)*ceq%rtn
       enddo
    endif
! restore phases set dormant
    jph=0
    if(gx%bmperr.ne.0) then
! save any error already set and clear error code
       saverr=gx%bmperr; gx%bmperr=0
    else
       saverr=0
    endif
    jj=meqrec%dormlink
1200 continue
    if(jj.ne.0) then
!       if(.not.btest(meqrec%status,MMQUIET)) &
!            write(*,*)'Restore from dormant: ',jj,meqrec%phr(jj)%iph,&
!            meqrec%phr(jj)%ics
       kk=meqrec%phr(jj)%phtupix
       phnames=' '
       call get_phasetup_name(kk,phnames)
       if(gx%bmperr.ne.0) then
          write(*,*)'MM cannot find phasetup name: ',jj,kk,gx%bmperr
          gx%bmperr=0
       endif
       if(.not.btest(meqrec%status,MMQUIET)) then
          if(meqrec%phr(jj)%curd%dgm.gt.zero) then
             write(*,1220)jj,kk,trim(phnames),meqrec%phr(jj)%curd%dgm
1220         format('MM Restoring phase:  ',2i5,2x,a,5x,1pe12.4)
          else
             write(*,1220)jj,kk,trim(phnames)
          endif
       endif
       if(meqrec%phr(jj)%curd%dgm.gt.1.0D-2) jph=jj
! do I have two places for suspendeded ?? YES!!
       meqrec%phr(jj)%phasestatus=PHENTUNST
! below is in the phase_varres record, previous is temporary equilibrium data
       meqrec%phr(jj)%curd%phstate=PHENTUNST
       jj=meqrec%phr(jj)%dormlink
       goto 1200
    endif
    if(jph.gt.0) then
       if(.not.btest(meqrec%status,MMQUIET)) &
            write(*,*)'MM warning, a restored phase wants to be stable:',jph
       gx%bmperr=4363
    endif
! we may already have had an error ...
    if(saverr.ne.0) gx%bmperr=saverr
! try to find problem with listed chemical potential    
! chempot(2) should be value with user defined reference state,
    if(gx%bmperr.eq.0) then
       do jj=1,meqrec%nrel
          xxx=zero
          lokph=ceq%complist(jj)%phlink
          if(lokph.gt.0) then
! we must also handle reference state at fix T !!
! lokph is index of phase in phlista, calcg_endmember want index in phases ....
!             write(*,*)'Component has defined reference state: ',jj,lokph
             tpvalsave=ceq%tpval
! modified calcg_endmember to convert negative phase index to phase number ...
!             write(*,*)'MM calling calcg_endmember 1: ',-lokph
! MUS same as TC MU
             call calcg_endmember(-lokph,ceq%complist(jj)%endmember,xxx,ceq)
             if(gx%bmperr.ne.0) then
                write(*,68)'MM error calculating reference state',gx%bmperr,&
                     -lokph,jj,xxx,tpvalsave(1),ceq%complist(jj)%endmember
68              format(a,3i5,2(1pe12.4),2x,10i3)
                ceq%tpval=tpvalsave
!                stop
                goto 998
             endif
          endif
! MU same as TC MUR
          ceq%complist(jj)%chempot(2)=ceq%complist(jj)%chempot(1)+xxx*ceq%rtn
       enddo
!    else
!       write(*,69)'Unable to calculate reference states due to errors'
!69     format(a)
    endif
!    write(*,37)'mu1: ',(ceq%complist(jj)%chempot(1),jj=1,meqrec%nrel)
!    write(*,37)'mu2: ',(ceq%complist(jj)%chempot(2),jj=1,meqrec%nrel)
!37  format(a,6(1pe12.4))
!-------------
998 continue 
    if(.not.formap) then
! if called during mapping keep phr
       deallocate(meqrec%phr)
    endif
! >>>> here one can allow new composition set in parallelization
    return
  end subroutine meq_phaseset

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine meq_sameset
!\begin{verbatim}
  recursive subroutine meq_sameset(irem,iadd,mapx,meqrec,phr,inmap,ceq)
! iterate until phase set change, converged or error (incl too many its)
! iadd = -1 indicates called from calculating a sequence of equilibria
! mapx is used when calling meq_sameset from step/map
    implicit none
    integer irem,iadd,inmap,mapx
    TYPE(meq_setup) :: meqrec
    TYPE(meq_phase), dimension(*), target :: phr
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer increase,ioff,ik,jj,jph,ie,ierr,jmaxy
    integer kk,kkz,level3,mph,negam,negamph,nj,nk,nl
    integer nz1,nz2
    TYPE(meq_phase), pointer :: pmi
! Using SAVE not possible for parallel calculations here once is just warning
    logical, save :: once=.true.
    double precision, dimension(5) :: qq
    double precision, dimension(maxconst) :: ycormax
    double precision, dimension(:,:), allocatable :: smat
    double precision, dimension(:), allocatable :: svar
! these arrays should maybe be allocated ....
    double precision, dimension(maxconst) :: ycorr,yarr
    integer converged,jz
    double precision chargefact,chargerr
    double precision dgm,summ,dgmmax,gsurf,phf,phs
    double precision prevmaxycorr,pv,signerr
    double precision xxx,ycormax2,yprev,ys,ysmm,ysmt,yss,yst
    double precision, parameter :: ylow=1.0D-3,ymin=1.0D-12,ymingas=1.0D-30
    double precision yvar1,yvar2
    double precision maxphch
    double precision sum
    double precision, dimension(:), allocatable :: cit
    double precision deltat,deltap,deltaam,yfact
! this is an emergecy fix to improve convergence for ionic liquid
    double precision, parameter :: ionliqyfact=3.0D-1
! to check if we are calculating a single almost stoichiometric phase ...
    integer iz,tcol,pcol,nophasechange,notagain
    double precision maxphasechange,molesofatoms,factconv
    double precision lastdeltat,deltatycond,phfmin,value
    integer notf,dncol,iy,jy,iremsave,phasechangeok,nextch,iremax,srem,errall
    character phnames*50
    double precision, dimension(:), allocatable :: lastdeltaam
    logical vbug,stoikph
! NOTE using save cannot be reconciled with parallel calculations
    save notagain
!
! do not allow return unless meqrec%noofits greater or equal to nextch
    mapx=0
    nextch=meqrec%noofits+4
    stoikph=.true.
    nophasechange=0
    maxphasechange=zero
! this is set each time the set of phases changes, controls change in T
! when there is a condition on y
    deltaTycond=2.5D1
    if(iadd.eq.-1 .or. ocv()) then
       write(*,*)'Debug output in meq_sameset'
       vbug=.TRUE.; iadd=0
    else
       vbug=.FALSE.
    endif
!    vbug=.TRUE.
    if(vbug)write(*,*)'entering meq_sameset',meqrec%nphase,irem
    iremsave=irem
! this is max correction of constituent fraction for each phases
    ycormax=zero
! magic trying to force decreasing step in fractions
!    ymagic=one
!    nmagic=0
! this is an attempt to decrease variation in phase amount corrections
    allocate(lastdeltaam(meqrec%nstph),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 9: ',errall
       gx%bmperr=4370; goto 1000
    endif
    lastdeltaam=zero
! dimension matrix for conditions, components+stable phases
    nz1=meqrec%nrel-meqrec%nfixmu+meqrec%nstph-meqrec%nfixph
    if(meqrec%tpindep(1)) nz1=nz1+1
    if(meqrec%tpindep(2)) nz1=nz1+1
    if(ocv()) write(*,11)meqrec%nrel,meqrec%nfixmu,meqrec%nstph,&
         meqrec%nfixph,meqrec%tpindep,nz1,ceq%tpval(1)
11  format('In meq_sameset, sysmat: ',4i7,2l2,i5,1pe12.4)
    nz2=nz1+1
    if(vbug) write(*,*)'Allocating smat: ',nz1
    allocate(smat(nz1,nz2),stat=errall)
    allocate(svar(nz1),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 10: ',errall
       gx%bmperr=4370; goto 1000
    endif
! check if constituent fraction correction in stable phases increases
! for each iteration.  Needed for the Re-V case ....
    prevmaxycorr=zero
    increase=0
    level3=0
! this is set TRUE after 3 iterations
    phasechangeok=meqrec%noofits
    if(phasechangeok.eq.1) then
       notagain=0
    endif
! debugging problem with changing axis in mapping
    if(ocv() .and. meqrec%tpindep(1)) write(*,*)'variable T: ',ceq%tpval(1)
!-------------------------------------------------------------
! return here until converged or phase set change
100 continue
    meqrec%noofits=meqrec%noofits+1
    cerr%flag=0
! nonzero flag means error output below
!    cerr%flag=1
    if(nophasechange.gt.100) then
       if(maxphasechange.lt.1.0E-10) then
! if we have not changed the set of stable phases for many iterations
! and the changes in phase amounts is small maybe we are calculationg an
! almost stoichiometric phase?  Changes in MU can be large!
          if(stoikph .and. meqrec%nphase.gt.1) then
! write this message if VERBOSE is set
             if(btest(globaldata%status,GSVERBOSE)) write(*,30)nophasechange,&
                  converged,cerr%nvs,ceq%tpval(1)
30           format('Slow converge at ',3i3,F10.2)
             if(cerr%flag.ne.0) then
                write(*,31)(cerr%typ(iz),cerr%val(iz),cerr%dif(iz),&
                     iz=1,cerr%nvs)
31              format('MM 31: ',3(i3,1pe12.4,e10.2))
             endif
! write message only (once for each minimization)
             stoikph=.false.
! if this happends during step/map give error message to force smaller steps
             if(inmap.eq.1 .and. meqrec%noofits.eq.ceq%maxiter) then
                gx%bmperr=4359; goto 1000
             endif
          endif
!+          converged=0
!+          goto 1000
!       else
! maybe use this to improve concergence??
!          if(.not.allocated(loopfact)) then
!             allocate(loopfact(meqrec%nrel))
!          endif
       endif
    endif
    nophasechange=nophasechange+1
    cerr%nvs=0
    cerr%mconverged=0
! this is magic ....
!    nmagic=nmagic+1
!    if(mod(nmagic,5).eq.0) ymagic=0.5*ymagic
!    if(mod(nmagic,25).eq.0) ymagic=one
! end of magic
!101 format(a)
!    write(*,*)'Iteration: ',meqrec%noofits,' ----------------------------- '
    if(ocv()) write(*,199)meqrec%noofits,ceq%tpval(1),meqrec%nstph,&
!    write(*,199)meqrec%noofits,ceq%tpval(1),meqrec%nstph,&
         (meqrec%stphl(jz),jz=1,meqrec%nstph)
199 format(/'Equil iter: ',i3,f8.2,', stable phases: ',i3,2x,100i3)
    if(meqrec%noofits.gt.ceq%maxiter) goto 1200
    converged=0
    if(vbug) write(*,*)'Iteration: ',meqrec%noofits,converged
! loop for all phases and composition sets, loop over phr
!    if(meqrec%tpindep(1)) write(*,*)'variable T: ',meqrec%noofits,ceq%tpval(1)
!
! >>>>>>>>>>>> here we can parallelize 
!
!-$omp parallel do private(pmi) shared(meqrec)
! nullify liquid pointer
    nullify(meqrec%pmiliq)
    parallel: do mph=1,meqrec%nphase
       pmi=>phr(mph)
! this routine calculates G and derivatives, the phase matrix and inverts it.
! it also calculates the amounts of moles of components in the phase
!-$     write(*,*)'Phase and tread: ',mph,omp_get_thread_num()
! to set correct pmiliq we must calculate all liquids first!!
!       write(*,*)'MM call onephase: ',pmi%iph,pmi%ics
       call meq_onephase(meqrec,pmi,ceq)
       if(gx%bmperr.ne.0) then
! using LAPCK gives severe problems if we do not stop
          goto 1000
          if(pmi%stable.eq.0) then
! if this happends for an unstable phase just continue but ensure it will
! not be stable (in a very crude way)
!             write(*,*)'Matrix inversion error for unstable phase',pmi%iph
             pmi%curd%gval(1,1)=one
             gx%bmperr=0
          else
! Inversion error for stable phase is fatal, error code already set
             if(once) then
                write(*,*)'Warning, matrix inversion problem: ',pmi%iph
                once=.false.
             else
                goto 1000
             endif
             gx%bmperr=0
          endif
       endif
!107       format(a,6(1pe12.3))
! end of pmi% scope
    enddo parallel
!-$omp end parallel do
!
!=======================================================================
! step 2: calculation of equil matrix
! Solve for chemical potentials and conditions using all stable phases
! The EQUIL MATRIX (smat) has one row for each stable phase and
! one row for each component representing a condition
! (If a fix phase condition or chem.pot. condition slightly different??)
!----------------------------------------
300 continue
!    if(vbug) write(*,301)'Calculating general equil matrix',meqrec%nfixmu,&
!    write(*,301)'Calculating general equil matrix',meqrec%nfixmu,&
!         meqrec%nfixph,meqrec%tpindep,meqrec%noofits
301 format(a,2i2,2l2,i5)
! some arguments here are redundant but kept for some
    call setup_equilmatrix(meqrec,phr,nz1,smat,tcol,pcol,&
         dncol,converged,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'Back from setup_equilmatrix',tcol
!=====================================================================
! debug output of equil matrix, last column is right hand side
!380 continue
!    open(33,file='eqmat.dat ',access='sequential',status='unknown')
!    write(33,*)'Equilibrium matrix',nz1
!    do iz=1,nz1
!       write(33,112)iz,(smat(iz,jz),jz=1,nz2)
!112 format('>',i4,1x,4(1pe15.6))
!    enddo
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> debug
! debug output to follow the minimization: all mu_i, and 
! for all stable phases np^alpha, G^alpha, and x^alpha_i
!    call calc_molmass(xdum,wdum,tmdum,wmdum,ceq)
!    write(*,116)'MM mu:',meqrec%nstph,(ceq%cmuval(iz),iz=1,meqrec%nrel),&
!         (xdum(iz),iz=1,meqrec%nrel)
!116 format(a,i3,6(1pe12.4))
!    do iz=1,meqrec%nstph
!       jj=meqrec%stphl(iz)
!       call calc_phase_molmass(phr(jj)%iph,phr(jj)%ics,&
!            xdum,wdum,tmdum,wmdum,dumdum,ceq)
!       if(gx%bmperr.ne.0) stop 'debug'
! amount of phase, G of phase, x_i of phase
!       write(*,116)'MM ph:',jj,phr(jj)%curd%amfu,smat(iz,nz2),&
!            (xdum(ioff),ioff=1,meqrec%nrel)
!    enddo
! end debug output
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if(vbug) then
! when convergence problem list smat here and (and svar below) and study!!!
       call list_conditions(kou,ceq)
       do iz=1,nz1
          write(*,228)'smat1:',(smat(iz,jz),jz=1,nz2)
       enddo
    endif
228 format(a,6(1pe12.4),(8x,6e12.4))
! This is an emergecy check that the smat matrix does not contain
! values >1E+50.  We should test for Infinity and NaN but how??
    do iz=1,nz1
       do jz=1,nz2
          if(abs(smat(iz,jz)).gt.1.0D+50) then
             write(*,118)iz,jz
118          format('meq_sameset has illegal values in equilibrium matrix',2i4)
             gx%bmperr=4354; goto 990
          endif
       enddo
    enddo
! HERE new values of chemical potentials and and amount of phases
    call lingld(nz1,nz2,smat,svar,nz1,ierr)
    if(ierr.ne.0) then
       if(vbug) write(*,*)'Error solving equil matrix 1',meqrec%noofits,ierr,&
            iremsave
       if(iremsave.gt.0) then
! parallel2 goes into a loop here when phase iremsave has been suspended
! after at has been set suspended .... fixed by not returning nonzero irem 
! equil matrix wrong at first iteration after removing a phase
! This can be caused by having no phase with solubility of an element
! (happened in Fe-O-U-Zr calculation with just C1_MO2 stable and C1 does not
! dissolve Fe).  Try to set back the last phase removed!!
          if(.not.btest(meqrec%status,MMQUIET)) then
             kk=meqrec%phr(iremsave)%phtupix
             phnames=' '
             call get_phasetup_name(kk,phnames)
             write(*,*)'Error, restoring previously removed phase: ',&
                  trim(phnames)
          endif
! NOTE: it should also be removed from the dormant list!!
          iadd=iremsave
          notagain=iremsave
          goto 1100
       endif
       if(vbug) then
          do iz=1,nz1
             write(*,228)'smat2:',(smat(iz,jz),jz=1,nz2)
          enddo
       endif
! debug output ...
!       write(*,229)'ce:',meqrec%noofits
!       call list_conditions(kou,ceq)
!       do iz=1,nz1
!          write(*,228)'smat2:',(smat(iz,jz),jz=1,nz2)
!       enddo
!       gx%bmperr=4203; goto 1000
    endif
! when problems output svar here !! (and smat1: above)
!    write(33,*)'Solution'
!    write(*,228)'PHMAT: ',(svar(jz),jz=1,nz1)
!    close(33)
!    write(*,228)'svar1:',(svar(jz),jz=1,nz1)
    if(vbug) write(*,228)'svar1:',(svar(jz),jz=1,nz1)
!
! if no error at first calculation after phase set change iremsave=0
    iremsave=0
    if(vbug) write(*,229)'pm: ',meqrec%noofits,(svar(iz),iz=1,nz1)
!    write(*,229)'pm: ',meqrec%noofits,(svar(iz),iz=1,nz1)
229 format(a,i3,6(1pe12.4))
!---------
! copy the chemical potentials, take care of fixed values ....
! new potentials are in svar(1..meqrec%nrel-meqrec%nfixmu)
    iz=1
    notf=1
    setmu: do ik=1,meqrec%nrel
       if(notf.le.meqrec%nfixmu) then
          if(ik.eq.meqrec%mufixel(notf)) then
! this potential is fixed, no incrementing "iz", ceq%cmuval(ik) is a condition
             ceq%complist(ik)%chempot(1)=meqrec%mufixval(1)*ceq%rtn
             notf=notf+1
             cycle setmu
          endif
       endif
!       if(abs(svar(iz)-ceq%cmuval(ik)).gt.ceq%xconv) then
       if(abs(svar(iz)-ceq%cmuval(ik)).gt.abs(ceq%xconv*ceq%cmuval(ik))) then
!          if(vbug) write(*,387)'Unconverged pot: ',iz,ik,&
          if(nophasechange.gt.100) then
! Attempt to improve convergence for a 15 component system ... failed
!             xxx=0.25D0*(3.0D0*svar(iz)+1.0D0*ceq%cmuval(ik))
!             write(*,387)'Uncnv pot: ',iz,ik,&
!                  svar(iz),ceq%cmuval(ik),xxx,abs(svar(iz)-ceq%cmuval(ik)),&
!                  abs(ceq%xconv*ceq%cmuval(ik))
!387          format(a,2i3,3(1pe14.5),2(1pe10.2))
! take mean value ... DO NOT TRY THIS IF IT IS NOT ALMOST CONVERGED!!!
!             svar(iz)=xxx
          endif
          converged=7
          cerr%mconverged=converged
       endif
       ceq%cmuval(ik)=svar(iz)
! svar(iz) is mu/RT, chemput is mu
       ceq%complist(ik)%chempot(1)=svar(iz)*ceq%rtn
       iz=iz+1
    enddo setmu
    ioff=meqrec%nrel-meqrec%nfixmu+1
!------------
! update T and P if variable
    if(meqrec%tpindep(1)) then
       xxx=ceq%tpval(1)
! check convergence
!       write(*,*)'Delta T: ',svar(ioff),1.0D2*ceq%xconv
!       if(abs(svar(ioff)).gt.1.0D2*ceq%xconv) then
! this convergece criteria needed for the CHO-gas calculation!!!
! but causes problem calculating phase diagrams ... inmap=1 for step/map
! OBS svar(ioff) is Delta T, not absolute value
       if(inmap.eq.0 .and. abs(svar(ioff)).gt.1.0D1*ceq%xconv) then
          converged=8
          cerr%mconverged=converged
       endif
! limit changes in T to +/- 20% of current value
       if(abs(svar(ioff)/ceq%tpval(1)).gt.0.2D0) then
          svar(ioff)=sign(0.2D0*ceq%tpval(1),svar(ioff))
       endif
! limit change in T when there is condition on y
       if(ycondTlimit) then
          deltat=svar(ioff)
! Suck it happend that svar(ioff) changed sign each iteration ....
          if(lastdeltat*deltat.lt.zero) then
             deltatycond=max(deltatycond-one,one)
! never increase during one minimization ...
!          else
!             deltatycond=2.5D1
          endif
          if(abs(svar(ioff)).gt.deltatycond) then
             if(svar(ioff).gt.zero) then
                svar(ioff)=deltatycond
             else
                svar(ioff)=-deltatycond
             endif
             write(*,*)'MM ycondTlimit: ',deltat,svar(ioff)
             lastdeltat=svar(ioff)
          endif
       endif
       deltat=svar(ioff)
! limit maximum change in deltat
       if(abs(deltat).gt.meqrec%tpmaxdelta(1)) then
          deltat=sign(meqrec%tpmaxdelta(1),deltat)
          if(ocv()) write(*,386)'limit the change in T: ',&
               ceq%tpval(1),deltat,svar(ioff)
386       format(a,3(1pe12.4))
       endif
       ceq%tpval(1)=ceq%tpval(1)+deltat
! problems here when -finit-local-zero is removed
       if(vbug) write(*,*)'T and deltaT:',ceq%tpval(1),deltat
       if(ceq%tpval(1).le.1.0D-2) then
          write(*,*)'Attempt to set a temperature less than 0.01 K !!!'
          gx%bmperr=4187; goto 1000
       endif
       ioff=ioff+1
    endif
    if(meqrec%tpindep(2)) then
! if pressure variable
       xxx=ceq%tpval(2)
! check convergence
! ??? svar(ioff) much too small!! why? add a factor ...
!       svar(ioff)=1.0D2*svar(ioff)
       if(abs(svar(ioff)).gt.1.0D4*ceq%xconv) then
          converged=8
          cerr%mconverged=converged
       endif
!       write(*,389)'HMS pv: ',ioff,converged,svar(ioff),ceq%tpval(2)
!389    format(a,2i3,4(1pe12.4))
       if(abs(svar(ioff)/ceq%tpval(2)).gt.0.2D0) then
          svar(ioff)=sign(0.2D0*ceq%tpval(2),svar(ioff))
       endif
       deltap=svar(ioff)
! limit the changes in P
       if(abs(deltap).gt.meqrec%tpmaxdelta(2)) then
          deltap=sign(meqrec%tpmaxdelta(2),deltap)
          if(ocv()) write(*,386)'limit the change in P: ',&
               ceq%tpval(2),deltap,svar(ioff)
       endif
       ceq%tpval(2)=ceq%tpval(2)+svar(ioff)
       if(ceq%tpval(2).le.1.0D-2) then
          write(*,*)'Attempt to set pressure lower than 0.01 Pa!!!'
          gx%bmperr=4187; goto 1000
       endif
       ioff=ioff+1
    endif
!------------
! update phase amounts, take care of fixed phases ....
! the change in amounts are in svar(ioff+...)
    negamph=0
    negam=0
    irem=0
    iremax=0
    phfmin=zero
! dncol+1 should be the first Delta_phase-amount
    ioff=dncol+1
! scale all changes in phase amount with total number of atoms. At present
! assume this is unity.  Without scaling phase changes can be +/-1E+11 or more
! which creates instabilities
    maxphch=zero
!    normphchange: do jph=1,meqrec%nstph
    normphchange: do jph=1,meqrec%nstph-meqrec%nfixph
       if(abs(svar(ioff+jph-1)).gt.maxphch) maxphch=abs(svar(ioff+jph-1))
    enddo normphchange
    if(maxphch.gt.one) then
       ioff=dncol+1
       do jph=1,meqrec%nstph-meqrec%nfixph
          svar(ioff+jph-1)=svar(ioff+jph-1)/maxphch
       enddo
    endif
!
    ioff=dncol+1
! do not change phase amounts the first iteration
!    write(*,554)svar
!554 format('MM svar: ',6(1pe12.4))
!    if(meqrec%noofits.eq.1) then
!       goto 555
!    endif
    phamount2: do jph=1,meqrec%nstph
! loop for all stable phases
       jj=meqrec%stphl(jph)
!       phr(jj)%curd%damount=zero
!       kkz=test_phase_status(phr(jj)%iph,phr(jj)%ics,xxx,ceq)
       kkz=phr(jj)%phasestatus
! new -4=hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fixed
       if(kkz.ge.PHENTUNST .and. kkz.le.PHENTSTAB) then
! phase is entered so its amount can change, -svar(ioff) is the change
          phs=phr(jj)%curd%amfu
          if(ioff.gt.size(svar)) then
! error here calculating Fe-Si-C with 2 phases set fix zero
! setting w(si)=w(c)=none and fix T; should have w(si) fix and T=none
             write(*,42)'MM Too many phases with variable amount',ioff,&
                  size(svar),meqrec%nstph,phr(jj)%iph
42           format(a,10i4)
            gx%bmperr=4193; goto 1000
          endif
          deltaam=svar(ioff)
! Sigli convergence problem, bad guess of start amount of phases??
! NOTE sign! -deltaam is the change in amount of phase, 
!          write(*,43)'Deltaam: ',meqrec%noofits,jj,deltaam,lastdeltaam(jph),&
!               phr(jj)%curd%amfu,phr(jj)%curd%amfu-deltaam
!43        format(a,2i3,6(1pe12.4))
! tried to avoid too large changes in phase amount, just made things worse
!          if(meqrec%noofits.lt.3 .and. &
!               abs(deltaam).gt.0.5D0*phr(jj)%curd%amfu) then
!             deltaam=sign(0.1D0*phr(jj)%curd%amfu,deltaam)
!             write(*,43)'Modified: ',meqrec%noofits,jj,deltaam
!          endif
! limit change in amount of phase
          if(abs(deltaam).gt.ceq%xconv) then
! For the equil O-U with conditions on N(O) and N(U) there is no problem
! with the amount of C1 but with N= and x(O)= the phase amount change varies
! with sign and converges very slowly.  Probably an interference with the
! charge balance criteria.
             if(lastdeltaam(jph)*deltaam.lt.zero) then
! wow, this seems to work ... other attmepts interfere directly with the
! charge balance so one should carefully check how they are connected...
!                deltaam=5.0D-1*deltaam
! The half worked to C1+tetragonal, it did not work for ionic liquid misc. gap
! and in that case there is no charge balance criteria ... suck
!                deltaam=5.0D-1*deltaam
! Dubbelt wow ... 0.2 works for both cases ... why?? More iterations though .. 
                deltaam=2.0D-1*deltaam
                if(ocv()) write(*,3)'Phase amount sign change: ',&
                     meqrec%noofits,jph,jj,phs,lastdeltaam(jph),deltaam
!                write(*,3)'Phase amount sign change: ',&
!                     meqrec%noofits,jph,jj,phs,lastdeltaam(jph),deltaam
3               format(a,3i3,6(1pe12.4))
             endif
             if(converged.lt.6) then
                converged=6
                cerr%mconverged=converged
             endif
             if(vbug) write(*,381)'Phase amount change: ',meqrec%noofits,jj,&
                  phs,deltaam
381          format(a,2i3,4(1pe12.4))
          endif
          lastdeltaam(jph)=deltaam
          if(phr(jj)%curd%amfu-deltaam.le.zero) then
             if(meqrec%nstph.eq.1) then
! this is the only stable phase!  cannot have negative or zero amount!
                deltaam=phr(jj)%curd%amfu-1.0D-2
             endif
          endif
!          if(-deltaam.gt.one) then
          if(abs(deltaam).gt.one) then
! try to prevent too large increase/decrease in phase amounts. 
! Should be related to total amount of components.
             if(.not.btest(meqrec%status,MMQUIET)) &
                  write(*,*)'Large change in phase amount: ',deltaam
!             deltaam=-one
             deltaam=sign(0.5D0,deltaam)
          endif
          if(abs(deltaam).gt.maxphasechange) then
! to allow checks when phase set does not change and amount changes are small
! like when calculating an almost stoichiometric composition like UO2 with
! n(o)=2 and n(u)=1 at low T
             maxphasechange=abs(deltaam)
          endif
! special test for Al-Ni fcc/fcc#2 two-phase
! Calculations with Al-Ni T=1000, x(al)=.2 gives just a single FCC phase
! possible problems that we change the amounts of the wrong composition set?
! HOWEVER, I found the error is the second derivatives are wrong!!
!          if(meqrec%noofits.lt.10) deltaam=0.1*deltaam
!          write(*,383)'MM phase change: ',meqrec%noofits,jj,&
!               phr(jj)%iph,phr(jj)%ics,phr(jj)%curd%amfu,deltaam,svar(ioff)
!383       format(a,2i3,2x,2i3,3(1pe12.4))
          phf=phr(jj)%curd%amfu-deltaam
          if(phs.gt.0.2D0 .and. phf.le.zero) then
! violent change of phase fractions in Siglis case, liquid change from 1 to 0
! Prevent changes larger than 0.1 if value larger than 0.5
! old value of amfu in phs
             phf=0.1D0
          endif
!          write(*,363)' >>>> Stable phase: ',jj,phr(jj)%iph,&
!               phr(jj)%ics,phf,phs,deltaam,sum
363          format(a,3i3,6(1pe12.4))
!          phr(jj)%curd%damount=deltaam
          ioff=ioff+1
       elseif(kkz.eq.PHFIXED) then
! phase is fix, there is no change in its amounts
          phf=phr(jj)%curd%amfu
!          write(*,*)'Fixed phase: ',jj,phf
       else
! phase is dormant or suspended, must not be stable!!!!
          call get_phase_name(phr(jj)%iph,phr(jj)%ics,phnames)
          if(gx%bmperr.ne.0) goto 1000
!          write(*,373)phr(jj)%iph,phr(jj)%ics,kkz
          write(*,373)trim(phnames),kkz
373       format('MM The phase ',a,' cannot vary its amount:',3i7)
          gx%bmperr=4194; goto 1000
       endif
! problem with Fe-O-U-Zr convergence, all phases disappear ??
!       write(*,364)'Stable phase: ',meqrec%noofits,jj,phr(jj)%iph,&
!       phr(jj)%ics,phf,phs,phr(jj)%prevam
!364    format(a,4i3,6(1pe12.4))
! make sure the driving force of stable phases to zero
       phr(jj)%curd%dgm=zero
       if(phf.lt.zero) then
! phase has negative amount, NOT ALLOWED if it is the only stable phase 
          if(meqrec%nstph-meqrec%nfixph.eq.1) then
!             write(*,367)'Trying to remove the only stable phase ',jj,&
!                  phr(jj)%curd%amfu
367          format(a,i3,1pe14.6)
             phf=0.5D0*phr(jj)%curd%amfu
             gx%bmperr=4195; goto 1000
          else
! select phase with most negative amount
             if(phf.lt.phfmin) then
                phfmin=phf
                iremax=jj
             endif
! trying to improve convergence by allowing phases to be removed quicker
!             write(*,363)'Phase with negative amount: ',jj,meqrec%noofits,0,&
!                  phf,phs,phr(jj)%prevam
!             if(phf.lt.-1.0D-2) phf=zero
             if(jj.ne.notagain .and. phr(jj)%prevam.lt.zero) then
! remove this phase if negative amount previous iteration also
                irem=jj
!                write(*,376)'meq_sameset remove: ',meqrec%noofits,nextch,&
!                     jj,notagain
376             format(a,4i4)
! jumping to 1000 here means constitutions not changed in this iteration
                goto 1000
             else
! mark this phase had negative amount this iteration
! PROBLEM removing one of two composition sets of the same phase,
! (miscibility gap), they may change which have negative amount each iteration
                phr(jj)%prevam=-one
                phf=zero
             endif
          endif
       else ! phase has positive amount, mark in prevam
          phr(jj)%prevam=one
       endif
! store the new phase fraction (moles formula units)
       phr(jj)%curd%amfu=phf
    enddo phamount2 ! end of loop for jph=1,meqrec%nstph
!555 continue
!
!    if(iremax.gt.0) then
!       write(*,*)'meq_sameset remove?',meqrec%noofits,iremax,phfmin
!    endif
    if(vbug) write(*,*)'finished updating phase amounts: ',&
         meqrec%noofits,phasechangeok,irem
!    if(meqrec%nfixmu.gt.0) then
!       write(*,33)'mu1: ',(ceq%cmuval(nj),nj=1,meqrec%nrel)
!       write(*,33)'mu2: ',(ceq%complist(nj)%chempot(1),nj=1,meqrec%nrel)
!       write(*,33)'mu3: ',(ceq%complist(nj)%chempot(2),nj=1,meqrec%nrel)
!       write(*,33)'mu4: ',(svar(nj),nj=1,meqrec%nrel)
!33     format(a,6(1pe12.4))
!    endif
!-------------------------------------------------------
! After solving the equil matrix and updating the chemical potentials,
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
!    chargefact=1.0D-1
! kk is used to check if a charged phase is stable,
! it is incremented for each stable phase
    kk=1
! iadd is set to the unstable phase with largest positive driving force
! dgmmax is the largest psoitive driving force
    iadd=0
    dgmmax=zero
    ysmm=zero
!-----------------------------------------------------
! Update the constitutions.  If irem>0 remove this phase unless
! we have made at least 3 iterations with the current phase set
    if(irem.gt.0 .and. meqrec%noofits-phasechangeok.gt.3) goto 1000
!--------------------------
! These are needed to avoid several phases have exactly the same fracions
! if the start guess is very bad and limitations are used
       yvar1=1.0D-4
       yvar2=1.0D-13
!-----------------------------------------
    lap: do jj=1,meqrec%nphase
! The current chemical potentials are in ceq%cmuval(i)
!       if(vbug) write(*,*)'Phase: ',phr(jj)%iph,phr(jj)%ics,&
!              phr(jj)%curd%amfu
       if(jj.eq.meqrec%stphl(kk)) then
! jj is stable, increment kk but do not make it larger than meqrec%nstph
! save index in meqrec%stphl in jph !!!!!!!!!!! kk never used !!!!!!!!!
          jph=kk
          kk=min(kk+1,meqrec%nstph)
!          if(meqrec%noofits.le.2) write(*,83)'dy1: ',jj,jph,kk
!83        format(a,3i3,6(1pe12.4))
       else ! phase is not stable
! calculate driving force for unstable phases. First calculate the sum
! of the current phase composition and the calculated chemical potentials
          jph=0
          gsurf=zero; summ=zero
          do ie=1,meqrec%nrel
! fatal parallel execution error once here
! index '1' of dimension 1 of array 'phr' above upper bound of 0
             gsurf=gsurf+phr(jj)%xmol(ie)*ceq%cmuval(ie)
             summ=summ+phr(jj)%xmol(ie)
          enddo
          gsurf=gsurf/summ
! calculate G_m plus any deltat and deltap terms
          dgm=phr(jj)%curd%gval(1,1)
          if(meqrec%tpindep(1)) then
             dgm=dgm+phr(jj)%curd%gval(2,1)*deltat
          endif
          if(meqrec%tpindep(2)) then
             dgm=dgm+phr(jj)%curd%gval(3,1)*deltap
          endif
! scale dgm per mole atoms
          molesofatoms=phr(jj)%curd%abnorm(1)
          if(molesofatoms.lt.0.5D0) then
! problem when BCC becomes just vacancies
!             write(*,*)'MM Phase: ',jj,' moles of atoms: ',molesofatoms
             molesofatoms=0.5D0
          endif
!          dgm=gsurf-dgm/phr(jj)%curd%abnorm(1)
          dgm=gsurf-dgm/molesofatoms
          if(phr(jj)%phasestatus.gt.0) then
! we should be here only for UNSTABLE phases, phr(jj)%phasestatus<=0
! For some reason a phase has entered/fixed status (>0) THAT IS AN ERROR
! It happened in SMP2A when mapping Al-Ni and correcting too long step in T
             write(*,'(a,i4,i3)')'MM phase status reset:',jj,phr(jj)%phasestatus
             phr(jj)%phasestatus=0
          endif
          if(dgm.gt.dgmmax) then
             if(phr(jj)%phasestatus.ge.PHENTUNST .and. &
                phr(jj)%phasestatus.le.PHENTERED) then
! phase is entered, can have status changed
! if this is another constitution set of an already stable phase then check
! below if the constitution of this phase is very similar to the stable one
                iadd=jj
                dgmmax=dgm
!                write(*,379)'meq_sameset add: ',meqrec%noofits,nextch,&
!                     iadd,dgmmax
379             format(a,3i4,4(1pe12.4))
             endif
          endif
! The difference between previous and current DGM is used to check for
! convergence below.  Very important to check if continue iterating!!
          phr(jj)%prevdg=phr(jj)%curd%dgm
          phr(jj)%curd%dgm=dgm
       endif
! Update constituent fractions for ALL phases, stable or not
! if phr(jj)%xdone=1 then phase has no composition variation
       if(phr(jj)%xdone.eq.1) cycle
!----------------------------------------------------
       allocate(cit(phr(jj)%idim),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 11: ',errall
          gx%bmperr=4370; goto 1000
       endif
       cit=zero
       if(meqrec%tpindep(1)) then
! variable T, code copied from calc_dgdyterms, cit(nj) used below
!          write(*,44)'index 1: ',jj,phr(jj)%ncc,phr(jj)%idim,&
!               size(phr(jj)%invmat)
          do jy=1,phr(jj)%ncc
             sum=zero
             do iy=1,phr(jj)%ncc
                sum=sum+phr(jj)%invmat(iy,jy)*&
                     phr(jj)%curd%dgval(2,iy,1)
             enddo
             cit(iy)=sum*deltat
!             write(*,44)'index 2: ',jj,jy,iy,0,sum
!44           format(a,4i3,6(1pe12.4))
          enddo
!! end copy
!          write(*,*)'Adding contribution from variable T to delta-y',&
!               phr(jj)%ncc
! missing code for correction due to variable P?????
       endif
! These are used to introduce some variation in fractions when the values
! exceed limits.  Otherwise one can as Sigli found have two stable phases
! with exactly the same fractions and have a crash
!
       moody: do nj=1,phr(jj)%ncc
          ys=zero
          do nk=1,phr(jj)%ncc
             pv=zero
             do nl=1,meqrec%nrel
! ceq%cmuval(nl) is the chemical potential of element nl (divided by RT)
! phr(jj)%dxmol(nl,nk) is the derivative of component nl
! wrt constituent nk
!                write(*,*)'ycorr: ',nl,ceq%complist(nl)%chempot(1)/ceq%rtn
!                write(*,612)'MM y1: ',nk,nl,&
!                     ceq%complist(nl)%chempot(1)/ceq%rtn,ceq%cmuval(nl)
!612             format(a,2i4,6(1pe12.4))
                pv=pv+ceq%complist(nl)%chempot(1)/ceq%rtn*phr(jj)%dxmol(nl,nk)
!                write(*,111)'pvx: ',nj,pv,ceq%complist(nl)%chempot(1),&
!                     ceq%rtn,phr(jj)%dxmol(nl,nk)
!                pv=pv+ceq%cmuval(nl)*phr(jj)%dxmol(nl,nk)
!                pv=pv+svar(nl)*phr(jj)%dxmol(nl,nk)
             enddo
             pv=pv-phr(jj)%curd%dgval(1,nk,1)
             ys=ys+phr(jj)%invmat(nj,nk)*pv
!             write(*,111)'pvx: ',nj,ys,pv,phr(1)%curd%dgval(1,nk,1),&
!                  phr(1)%invmat(nj,nk)
!111          format(a,i2,6(1pe12.4))
          enddo
          if(phr(jj)%chargebal.eq.1) then
! For charged phases add a term 
! phr(jj)%invmat(phr(jj)%idim,phr(jj)%idim)*Q
             ys=ys-chargefact*phr(jj)%invmat(nj,phr(jj)%idim)*&
                  phr(jj)%curd%netcharge
!             ys=ys-chargefact*phr(jj)%invmat(nj,phr(jj)%idim)*&
!                  phr(jj)%charge
! jph is nonzero only for stable phases
             if(jph.gt.0 .and. &
!             if(jj.eq.meqrec%stphl(kk) .and. &
! Hm, is this check correct?  kk is updated above to be the next stable phase..
!                  abs(phr(jj)%charge).gt.chargerr) then
!                chargerr=abs(phr(jj)%charge)
!                signerr=phr(jj)%charge
                  abs(phr(jj)%curd%netcharge).gt.chargerr) then
                chargerr=abs(phr(jj)%curd%netcharge)
                signerr=phr(jj)%curd%netcharge
             endif
!             write(*,*)'Charge: ',jj,phr(jj)%netcharge
          endif
! when T is variable
          ycorr(nj)=ys+cit(nj)
          if(abs(ycorr(nj)).gt.ycormax2) then
             ycormax2=ycorr(nj)
          endif
! Sigli converge problem, fixed by changing stable phases in different order
!          write(*,111)converged,jj,nj,ys
!111       format('Y corr: cc/ph/cons/y: ',i2,2i4,1pe12.4)
! should possibly be ycorr(nj) instead of ys (ycorrmax)
          if(abs(ys).gt.ceq%xconv) then
! if the change in any constituent fraction larger than xconv continue iterate
!             write(*,*)'Convergence criteria, phase/const: ',jj,nk
             if(phr(jj)%stable.eq.0) then
! Phase is not stable
! Handle convergence criteria different if inmap=1 or not
                mapping7: if(inmap.eq.0) then
! we are NOT in STEP/MAP, increase convergence criteria to handle
! the Mo-Ni-Re 3 phase equilibria
                   if(abs(ys).gt.1.0D1*phr(jj)%curd%yfr(nj)) then
! for unstable phases the corrections must be smaller than ...????
                      if(converged.lt.3) then
                         converged=3
                         cerr%mconverged=converged
                         yss=ys
                         yst=phr(jj)%curd%yfr(nj)
                      endif
                   elseif(abs(ys).gt.1.0D2*ceq%xconv) then
!212                   format(a,3i3,i4,4(1pe12.4))
                      if(converged.lt.4) then
                         if(phr(jj)%ncc.gt.10) then
! Calculation with the COST507 database and 20 elements too many iterations
! ... allow larger gdconv(1) 
                            factconv=1.0D1
                         else
                            factconv=one
                         endif
                         if(phr(jj)%curd%dgm-phr(jj)%prevdg.gt.&
                              factconv*ceq%gdconv(1)) then
! Must be less than this  if(phr(jj)%curd%dgm-phr(jj)%prevdg.gt.5.0E-3) then
                            converged=4
                            cerr%mconverged=converged
                            yss=ys
                            yst=phr(jj)%curd%yfr(nj)
                         endif
                      endif
                   else
                      if(converged.eq.0) then
                         converged=1
                         cerr%mconverged=converged
                         yss=ys
                         yst=phr(jj)%curd%yfr(nj)
                      endif
                   endif
                else
! we are doing step/map NO CHANGE, use old convergence criteria
! otherwise step1 and mmap4 are uncomplete with those above ...
                   if(abs(ys).gt.1.0D1*phr(jj)%curd%yfr(nj)) then
! for unstable phases the corrections must be smaller than ...????
                      if(converged.lt.3) then
                         converged=3
                         cerr%mconverged=converged
                         yss=ys
                         yst=phr(jj)%curd%yfr(nj)
                      endif
                   elseif(abs(ys).gt.1.0D2*ceq%xconv) then
! maybe accept 100 times larger correction than for stable phases
!                   write(*,107)'metast ph ycorr: ',ys,&
!                        phr(jj)%curd%yfr(nj)
                      if(converged.lt.2) then
                         converged=2
                         cerr%mconverged=converged
                         yss=ys
                         yst=phr(jj)%curd%yfr(nj)
                      endif
                   else
                      if(converged.eq.0) then
                         converged=1
                         cerr%mconverged=converged
                         yss=ys
                         yst=phr(jj)%curd%yfr(nj)
                      endif
                   endif
                endif mapping7
             elseif(converged.lt.4) then
! large correction in fraction of constituent fraction of stable phase
!                write(*,*)'mm converged 4A: ',jj,nj,ys
                converged=4
                cerr%mconverged=converged
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
! check if the change in any fraction is larger than the fraction ...
             if(ycorr(nj).gt.phr(jj)%curd%yfr(nj)) then
!                write(*,612)'MM y2: ',jj,nj,ycorr(nj),phr(jj)%curd%yfr(nj)
                if(converged.lt.4) then
                   converged=4
                   cerr%mconverged=converged
                endif
             endif
          endif
       enddo moody
! end of correction of y fractions
!---------------------------------
! Limit change in fractions .... all ycorr(nj) multiplied with same factor
! keeping the sum of corrections in all sublattices as zero
       if(vbug) write(*,74)'maximum corr: ',&
            meqrec%noofits,jj,ycormax2,ycormax(jj)
74     format(a,2i3,2(1pe12.4))
       if(ycormax(jj)*ycormax2.le.zero) then
! the condition is zero at first step, limit that
          yfact=one/(2.0D0+abs(ycormax2))
          ycormax2=yfact*ycormax2
       elseif(phr(jj)%ionliq.gt.0 .and. ycormax2.lt.1.0D-4) then
! step seems to be very small ... try to decrease number of iteration
          yfact=2.0d0
       else
          yfact=one
       endif
       moody2: do nj=1,phr(jj)%ncc
! all corrections of constituent fractions in ycorr(1..phr(jj)%ncc)
! ymagic is halfed every 5th iteration when same phase set, after 5 times reset
          yprev=phr(jj)%curd%yfr(nj)
!          yarr(nj)=yprev+ycorr(nj)
          if(phr(jj)%ionliq.gt.0) then
! For ionic liquids, an even smaller step is allowed ...
! The O-Pu-U test case converged up to 2800 without any particular factor
! with a factor 0.4 it converged up to 3000K (~150 its), yfact does not
! has any significant influence. 
!             yarr(nj)=yprev+4.0D-1*ycorr(nj)*yfact
! tafidbug, 0.2 created problems
!             yarr(nj)=yprev+2.0D-1*ycorr(nj)*yfact
!             yarr(nj)=yprev+3.0D-1*ycorr(nj)*yfact
             yarr(nj)=yprev+ionliqyfact*ycorr(nj)*yfact
!             yarr(nj)=yprev+ycorr(nj)*yfact
!             write(*,281)'ycorr: ',nj,yfact,yprev,yarr(nj)
!281           format(a,i3,6(1pe12.4))
          else
             yarr(nj)=yprev+ycorr(nj)*yfact
          endif
!          if(vbug) then
! output to check reasons for bad convergence
!             write(*,57)'MM y&dy ',phr(jj)%iph,phr(jj)%ics,&
!                  phr(jj)%stable,nj,&
!                  ys,cit(nj),phr(jj)%curd%yfr(nj),yarr(nj),ycorr(nj)
!57           format(a,3i2,i3,5(1pe12.4))
!          endif
          if(yarr(nj).lt.ymin) then
! this added to avoid too drastic jumps in small fractions
! The test case ccrfe1.OCM needs this
             if(yprev.gt.ylow) then
!                write(*,*)'Applying fraction change limitation 4 ',jj
                yarr(nj)=0.9*ylow
             elseif(test_phase_status_bit(phr(jj)%iph,PHGAS)) then
! for gas phase one must allow smaller constituent fractions
                if(yarr(nj).lt.ymingas) then
                   yarr(nj)=ymingas
                endif
             else
!                write(*,*)'Applying fraction change limitation 5 ',jj
                yarr(nj)=ymin+yvar2
                yvar2=2.0D0*yvar2
                if(yvar2.gt.1.0D-11) yvar2=1.0D-13
             endif
          endif
          if(yarr(nj).gt.one) then
!             write(*,*)'Applying fraction change limitation 6 ',jj
             yarr(nj)=one-yvar1
             yvar1=2.0D0*yvar1
             if(yvar1.gt.1.0D-3) yvar1=1.0D-4
          endif
       enddo moody2 ! end loop for all constituents nj in phase jj
!
       ycormax(jj)=ycormax2
! >>>>>>>>>>>>>>>>>> HERE the new constitution is set <<<<<<<<<<<<<<<<<<<<<
!       if(meqrec%noofits.le.2) write(*,83)'dy2: ',jj,phr(jj)%iph,kk,&
!            (yarr(nj),nj=1,phr(jj)%ncc)
!       write(*,114)'YARR: ',jj,phr(jj)%ics,(yarr(nj),nj=1,phr(jj)%ncc)
!114       format(a,2i3,8(F7.4))
!       write(*,*)'MM calling set_constitution 1:',phr(jj)%iph,phr(jj)%ics
       call set_constitution(phr(jj)%iph,phr(jj)%ics,yarr,qq,ceq)
       if(gx%bmperr.ne.0) goto 1000
!  >>>>>>>>>>>>>>>>>> for all phases <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       deallocate(cit)
    enddo lap
! finished correction of all constituent fractions in all phases
!-------------------------------------------------------
!    do jph=1,meqrec%nstph
!       jj=meqrec%stphl(jph)
!       write(*,393)'Stable phase: ',phr(jj)%iph,phr(jj)%ics,&
!            phr(jj)%curd%amfu
!    enddo
!393 format(a,2i4,6(1pe12.4))
! check if fraction corrections in stable phases increases
! it solved a problem in ReV when fractions initially changed very little
! but the change increased each iteration
    if(meqrec%noofits.gt.8) then
! this means minimum 8 iterations!!
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
! The request for 100 times better than ceq%xconv is OK with conditions 
! N(U)= N(O)= but not with N= x(O)=
!    if(chargerr.gt.1.0D-2*ceq%xconv) then
! strengthen charge balance convergence criteria
    if(chargerr.gt.ceq%xconv) then
       if(ocv()) write(*,654)'Charge error: ',signerr,chargerr,ceq%xconv
654    format(a,6(1pe12.4))
       if(converged.lt.6) then
          converged=6
          cerr%mconverged=converged
       endif
    endif
!-------------------------------------------------------
    if(converged.eq.3) then
! force one extra iterations with large fraction variations in unstable phases
!       write(*,267)'End of iteration: ',meqrec%noofits,converged,&
!            increase,yss,yst
       level3=level3+1
    elseif(converged.eq.4) then
! this means large fraction variations in stable phases
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
!    call get_state_var_value('X(O) ',value,phnames,ceq)
! trying to understand how STEP/MAP sets fix phases ....
!    write(*,*)'MM Fraction of O: ',value
    if(iadd.gt.0) then
! check if phase to be added is already stable as another composition set
! This check should maybe be above as maybe another phase want to be stable??
       if(same_composition(iadd,phr,meqrec,ceq,dgm)) iadd=0
    endif
! check if phase iadd is stoichiometric and if so check of any stable phase
! phase that is stoichiometric has the same composition!!  IF SO
! remove that phase at the same time ...
    srem=0
    if(meqrec%nrel.gt.1 .and. iadd.gt.0) then
! skip this for unary system!!!
       jy=meqrec%phr(iadd)%phtupix
       samestoi: do nj=1,meqrec%nstph
! loop through all stable phases for other phase with same stoichiometry
          jj=meqrec%stphl(nj)
          if(jj.ne.iadd) then
             iy=meqrec%phr(jj)%phtupix
! check if same composition ... how? same_stoik in gtp3Y.F90
             if(same_stoik(jy,iy)) then
                srem=jj
                exit samestoi
             endif
          endif
       enddo samestoi
    endif
    if(srem.gt.0) then
       jy=meqrec%phr(iadd)%phtupix
       call get_phasetup_name(jy,phnames)
       iz=len_trim(phnames)+2
       call get_phasetup_name(iy,phnames(iz:))
!       write(*,*)'MM Same stoichiometry: ',trim(phnames),inmap,value
! try to handle this by calculating the T when the two stochiometric phases
! has the same Gibbs energy.  Use this only if maping and T is not a condition
       if(inmap.ne.0) then
! inmap=0 if we are not in a step/map calculation
! I do not understand why iy and jy here ?? I think iadd and srem ...
          call two_stoich_same_comp(iy,jy,mapx,meqrec,inmap,ceq)
       endif
       iadd=iy; irem=jy
!       write(*,*)'Phases: ',iadd,irem
! after this routine set the error code to return to mapping
!       stop 'same stoichimetries'

! to be handelled either by map/step routines or meq_phaseset
       gx%bmperr=4364; goto 1000
    endif
!    if(meqrec%noofits.gt.2 .and. (irem.gt.0 .or. iadd.gt.0)) then
    if(irem.gt.0 .or. iadd.gt.0) then
       goto 1100
    endif
!--------------------------------------------------------------------
!    write(*,*)'Iterations and convergence: ',meqrec%noofits,converged
!--------------------------------------------------------------------
! check convergence
!    if(meqrec%noofits.gt.400) then
!       write(*,778)'Test converged: ',meqrec%noofits,converged
!778    format(a,2i4)
!    endif
!------------------------------------------------------------
! This output gives a good indication for convergence problem
    if(vbug) write(*,*)'Convergence criteria: ',converged,level3
! converged=1 or 2 means constituent fraction in metastable phase not converged
    if(converged.gt.3) goto 100
! converged 3 means large change conts. fraction of unstable phase change a lot
! level3 is nuber of previous iteration with converged=3
! with allcost I had the correct equilibrium but occational converged=4
! probably because a metastable liquid with almost identical composition
! as the stable interfeared. Accept converged=3 twice in a row as correct!!
!    if(converged.eq.3 .and. level3.lt.4) goto 100
    if(converged.eq.3 .and. level3.lt.2) goto 100
! converged 4 means a constituent fraction of a stable phase change a lot
! converged=5 means a condition not fullfilled
! converged=6 means charge balance not converged or large phase fraction change
! converged=7 means large change in chemical potentials
! converged=8 means large change T or P
! always force 4 iterations, there is a minimum above forcing 9 iterations.
    if(meqrec%noofits.lt.4) goto 100
    if(increase.ne.0) then
! continue if corrections in constituent fractions in stable phases increases
! This is needed to change fractions in a gas from 1E-20 to some significant
! value
       goto 100
    endif
!------------------------
! equilibrium calculation converged, do some common thing
!    write(*,*)'Converged: ',converged
    goto 800
!
!==============================================================
! equilibrium calculation converged, save chemical potentials (svar*RT)
800 continue
!------------------------------------------------------
! do not save system matrix but save -dimension for use with derivatives
    ceq%sysmatdim=-nz1
! but save components with fix mu and fix phases
    ceq%nfixmu=meqrec%nfixmu
    if(allocated(ceq%fixmu)) deallocate(ceq%fixmu)
    if(ceq%nfixmu.gt.0) then
       allocate(ceq%fixmu(ceq%nfixmu),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 12: ',errall
          gx%bmperr=4370; goto 1000
       endif
       do ie=1,ceq%nfixmu
          ceq%fixmu(ie)=meqrec%mufixel(ie)
       enddo
    endif
    ceq%nfixph=meqrec%nfixph
    if(allocated(ceq%fixph)) deallocate(ceq%fixph)
    if(ceq%nfixph.gt.0) then
       allocate(ceq%fixph(2,ceq%nfixph),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 13: ',errall
          gx%bmperr=4370; goto 1000
       endif
       do ie=1,ceq%nfixph
! phase and composition set numbers
          ceq%fixph(1,ie)=meqrec%fixph(1,ie)
          ceq%fixph(2,ie)=meqrec%fixph(2,ie)
       enddo
    endif
!-------------------------------------
    if(vbug) write(*,*)'At 800 in meq_sameset: ',meqrec%nrel
    ceq%rtn=globaldata%rgas*ceq%tpval(1)
    do ie=1,meqrec%nrel
       ceq%complist(ie)%chempot(1)=ceq%cmuval(ie)*ceq%rtn
!       write(*,*)'Chempot/RT: ',cea%cmuval(ie),svar(ie)
    enddo
! list stable phases on exit
!    do jph=1,meqrec%nstph
!       jj=meqrec%stphl(jph)
!       write(*,393)'Stable phase Z: ',phr(jj)%iph,phr(jj)%ics,&
!            phr(jj)%curd%amfu
!    enddo
! set status of the stable phases on exit
    do jph=1,meqrec%nstph
       jj=meqrec%stphl(jph)
       call mark_stable_phase(phr(jj)%iph,phr(jj)%ics,ceq)
!       write(*,393)'Stable phase Z: ',phr(jj)%iph,phr(jj)%ics,&
!            phr(jj)%curd%amfu
    enddo
!----------------------
! save inverted phase matrix and more for future use when calculating H.T etc
! If already allocated then dealloc/alloc as number of constituents can change
!    if(vbug) write(*,*)'allocate/deallocate in meq_sameset: ',meqrec%nphase
    do jj=1,meqrec%nphase
       if(allocated(phr(jj)%curd%cinvy)) then
          deallocate(phr(jj)%curd%cinvy)
          deallocate(phr(jj)%curd%cxmol)
          deallocate(phr(jj)%curd%cdxmol)
       endif
! why is the dimension if invmat so different???
       ie=phr(jj)%idim
       if(vbug) write(*,*)'Save inverted phase matrix in meq_sameset: ',jj,ie
!       ie=int(sqrt(real(size(phr(jj)%invmat)))+0.1)
!       write(*,*)'Size: ',ie,phr(jj)%ncc
       allocate(phr(jj)%curd%cinvy(ie,ie),stat=errall)
       allocate(phr(jj)%curd%cxmol(meqrec%nrel),stat=errall)
       allocate(phr(jj)%curd%cdxmol(meqrec%nrel,phr(jj)%ncc),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 14: ',errall
          gx%bmperr=4370; goto 1000
       endif
       phr(jj)%curd%cinvy=phr(jj)%invmat
       phr(jj)%curd%cxmol=phr(jj)%xmol
       phr(jj)%curd%cdxmol=phr(jj)%dxmol
!----------------------
    enddo
    goto 1000
! output of equilibrium matrix when error return
990 continue
    do iz=1,nz1
       write(*,228)'smat1:',(smat(iz,jz),jz=1,nz2)
    enddo
!
1000 continue
    if(gx%bmperr.ne.0) then
       ceq%status=ibset(ceq%status,EQFAIL)
!      write(*,*)'minimization error: ',gx%bmperr
!   elseif(irem.eq.0 .and. iadd.eq.0) then
    endif
! jump here if phase change
1100 continue
    if(vbug) write(*,*)'Deallocating smat and svar'
    deallocate(smat)
    deallocate(svar)
    if(vbug) write(*,*)'Final return from meq_sameset'
!    if(gx%bmperr.ne.0) write(*,*)'Error return from meq_sameset',gx%bmperr
!    if(irem*iadd.gt.0) write(*,*)'Leaving meq_sameset: ',irem,iadd
!    write(*,*)'Exit meq_sameset'
    return
! too many iterations
1200 continue
!    write(*,*)'Too many iterations: ',meqrec%noofits,ceq%maxiter
    gx%bmperr=4204
    goto 1000
  end subroutine meq_sameset
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable setup_comp2cons
!\begin{verbatim}
  subroutine setup_comp2cons(meqrec,phr,nz1,smat,tval,xknown,converged,ceq)
! calculate internal equilibrium in a phase for given overall composition
! meqrec and phr contains data for phases, nz1 is dimension of equlibrium
! matrix, smat is the equilibrium matrix, tval is fixed T and P
! xknown is the overall composition
    TYPE(meq_setup) :: meqrec
    TYPE(meq_phase), dimension(*), target :: phr
    double precision smat(nz1,*),tval(*),xknown(*)
    integer nz1,converged
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
!    TYPE(gtp_condition), pointer :: condition,lastcond
    TYPE(meq_phase), pointer :: pmi
! cmix dimensioned for 2 terms ...
    integer tcol,pcol,dncol
    integer sel,jph,jj,ie,je,ncol
    integer nz2,nrow,errall
    double precision cvalue,totam,pham,mag,mat,map,xxx
! the next line of values are a desperate search for a solution
!    double precision amount
!    double precision hmval
    double precision, dimension(:), allocatable :: xcol,mamu
!    double precision, allocatable :: xxmm(:),wwnn(:),hval(:)
!    logical :: calcmolmass
!    character encoded*32
!-------------------------------------------------------------------
! Formulating the equil equation in general:
! Variables (one column per variable):
! - The chemical potentials of the components:     MEQREC%NREL
!   minus the number of fixed chemical potentials: -MEQREC%NFIXMU
! - The variation in T if not fixed                +1
! - The variation in P if not fixed                +1
! - The variation of the amounts of stable phases: MEQREC%NSTPH
!   minus those that have fixed amount:            -MEQREC%NFIXPH
!
! The variables will be ordered: MU, DeltaT, DeltaP, Delta Phase amounts
! this is important for the order of columns in the equil matrix
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
!    write(*,*)'MM: in comp2cons'
    allocate(mamu(meqrec%nrel),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 15: ',errall
       gx%bmperr=4370; goto 1000
    endif
!    goto 1000
! zero all values in equil matrix, dimension (nz1)x(nz1)
    nz2=nz1+1
    tcol=0
    pcol=0
    dncol=0
!-----------------------------------------------------------
! step 2.1 the Gibbs energies for the phases, we have just one !!
!    allstableph: do jph=1,meqrec%nstph
    jph=1
    jj=meqrec%stphl(jph)
! one column with amount of each component to be multiplied with the
! chemical potential
    ncol=1
    xxx=zero
    gloop: do je=1,meqrec%nrel
! I cannot understand how smat changes columns and rows !!!!
       smat(1,ncol)=phr(jj)%xmol(je)
!       smat(ncol,1)=phr(jj)%xmol(je)
       ncol=ncol+1
    enddo gloop
! column nz2 is the right hand side of the equation, the molar G
!?    smat(jph,nz2)=phr(jj)%curd%gval(1,1)
    smat(1,nz2)=phr(jj)%curd%gval(1,1)
!?    write(*,11)'MM smat1: ',1,(smat(1,ncol),ncol=1,nz2)
!    do nrow=1,nz1
!       write(*,11)'MM smat1: ',nrow,(smat(nrow,ncol),ncol=1,nz2)
!    enddo
!11  format(a,i2,6(1pe12.4))
!------------------------------------------------------------
! insert code to calculate N(A)=fix for all elements in this phase
!
!    case(11) ! N or X with or without indices and normalization
!1100   continue
    nrow=1
! conditions are N(A)=fix for all elements
    elloop1: do sel=1,meqrec%nrel
! Formulate equation for total amount N:
! rhs:  N-N+\sum_alpha N^a + \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j
! \sum_B \sum_alpha N^a \sum_i \sum_j dM^_A/dy_i dM^a_B/dy_j*z^a_ij  *mu(B)
!        \sum_alpha N^a \sum_i d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j      *deltaT
!        \sum_alpha N^a \sum_i d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j      *deltaP
!        \sum_A M^a_A                                    *deltaN^a
       allocate(xcol(nz2),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 16: ',errall
          gx%bmperr=4370; goto 1000
       endif
       xcol=zero
!       totam=zero
!       nallph: do jj=1,meqrec%nphase
! we have just one phase
       jj=1
       pmi=>phr(jj)
! moles formula units of phase ??
       pham=one
! multiply terms with the inverse phase matrix
       ie=sel
! MAYBE use calc_dgdyterms1X ??
       call calc_dgdyterms1(meqrec%nrel,ie,meqrec%tpindep,&
            mamu,mag,mat,map,pmi,ceq%cmuval,meqrec%noofits)
       if(gx%bmperr.ne.0) goto 1000
! the call above calculates (A is "ie", z_ij is the inverted phase matrix): 
! mamu_A(B=1..nrel) = \sum_i \sum_j dM^a_A/dy_i dM^a_B/dy_j z^a_ij
! mag_A             = \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j
! mat_A             = \sum_i \sum_j d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j
! map_A             = \sum_i \sum_j d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j
! calculate a term for each column to be multiplied with chemical potential
! if the potential is fixed add the term to the rhs
       ncol=1
       elloop2: do je=1,meqrec%nrel
! mamu(B) = \sum_i \sum_j \sum_A dM^a_B/dy_i dM^a_A z^a_ij
          xcol(ncol)=xcol(ncol)-pham*mamu(je)
          ncol=ncol+1
       enddo elloop2
! last columns on lhs are amounts of element ie for all stable non-fix phases
! dncol should indicate last column with potential, can be different for
! derivative, notf is set above
! Amount of component in phase
!       totam=totam+pham*pmi%xmol(sel)
! pmi%xmol(sel) is the M per formula unit, not mole fraction !!!!
       jj=size(pmi%xmol)
!       write(*,411)'xmol: ',jj,pmi%sumxmol,(pmi%xmol(ncol2),ncol2=1,jj)
!411    format(a,i2,6(1pe12.4))
       totam=pham*pmi%xmol(sel)/pmi%sumxmol
       xcol(ncol)=pham*pmi%xmol(sel)
! right hand side (rhs) contribution is
! - NP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij
       xxx=xcol(nz2)
!       write(*,11)'MM xxx: ',nrow+1,(xcol(je),je=1,nz2)
       xcol(nz2)=xcol(nz2)-pham*mag
!
! in xcol are values summed over all phases and components
! then copy summed columns to row nrow in matrix smat
       nrow=nrow+1
       if(nrow.gt.nz1) then
          write(*,*)'MM too many equations 11A',nrow
          gx%bmperr=4212; goto 1000
       endif
       do ncol=1,nz2
          smat(nrow,ncol)=xcol(ncol)
       enddo
       deallocate(xcol)
! add N^prescribed - N^current to rhs (right hand side)
! cvalue is the prescibed composition assuming one F.U. of phase ...??
       cvalue=xknown(sel)
       smat(nrow,nz2)=smat(nrow,nz2)-cvalue+totam
!       write(*,11)'MM row: ',nrow,cvalue,totam,(smat(nrow,ncol),ncol=1,nz2)
! relative check for convergence if cvalue>1.0
       conv: if(abs(totam-cvalue).gt.ceq%xconv*max(1.0d0,abs(cvalue)))then
          if(converged.lt.5) then
             converged=5
!             write(*,*)'1: converged=5',cerr%nvs
             cerr%mconverged=converged
             if(cerr%nvs.lt.10) then
                cerr%nvs=cerr%nvs+1
                cerr%typ(cerr%nvs)=5
                cerr%val(cerr%nvs)=cvalue
                cerr%dif(cerr%nvs)=totam-cvalue
!             write(*,266)'Unconverged condition N or N(A): ',sel,cvalue,totam
!266          format(a,i3,4(1pe12.4))
             endif
          endif
       endif conv
    enddo elloop1
!----------------------------------------------------------
! all conditions set
!380 continue
! there was a strange error that the matrix had been changed on return ...
!    do nrow=1,nz1
!       write(*,11)'MM smat2: ',nrow,(smat(nrow,ncol),ncol=1,nz2)
!    enddo
    goto 1000
1000 continue
    return
  end subroutine setup_comp2cons

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine setup_equilmatrix
!\begin{verbatim}
  subroutine setup_equilmatrix(meqrec,phr,nz1,smat,tcol,pcol,&
       dncol,converged,ceq)
! handels external conditions on extensive variables in the equil matrix
! meqrec and phr contains data for phases, nz1 is dimension of equlibrium
! matrix, smat is the equilibrium matrix, tcol and pcol are columns for
! variable T or P, dncol is the column with phase amount variables.
! converged is used to indicate calling routine and set if not converged
! external variable.
    TYPE(meq_setup) :: meqrec
    TYPE(meq_phase), dimension(*), target :: phr
    double precision smat(nz1,nz1+1)
    integer nz1,tcol,pcol,converged,dncol
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    TYPE(gtp_condition), pointer :: condition,lastcond
    TYPE(meq_phase), pointer :: pmi
! cmix dimensioned for 2 terms ...
    integer cmix(22),cmode,stvix,stvnorm,sel,sph,scs,jph,jj,ie,je,ke,ncol
    integer notf,nz2,nrow,nterms,mterms,moffs,ncol2,iph
    integer xterm,yindex,jy,errall
    double precision cvalue,totam,pham,mag,mat,map,xxx,zval,xval,ccf(5),evalue
! the next line of values are a desperate search for a solution
    double precision totalmol,totalmass,check1,check2,amount,mag1,mat1,map1
    double precision hmval,gref,tpvalsave(2),cib
    double precision, dimension(:), allocatable :: xcol,mamu,mamu1,zcol,qmat
    double precision, allocatable :: xxmm(:),wwnn(:),hval(:)
    logical :: vbug=.FALSE.,calcmolmass,notdone,nosave
    double precision bbug,dvalue
    character encoded*32,name*32
! For saving calculated terms in calc_dgdyterms
!    type(saveddgdy), target :: savedrec
!    type(saveddgdy), pointer :: saved
!-------------------------------------------------------------------
! Formulating the equil equation in general:
! Variables (one column per variable):
! - The chemical potentials of the components:     MEQREC%NREL
!   minus the number of fixed chemical potentials: -MEQREC%NFIXMU
! - The variation in T if not fixed                +1
! - The variation in P if not fixed                +1
! - The variation of the amounts of stable phases: MEQREC%NSTPH
!   minus those that have fixed amount:            -MEQREC%NFIXPH
!
! The variables will be ordered: MU, DeltaT, DeltaP, Delta Phase amounts
! this is important for the order of columns in the equil matrix
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
! A serious bug concerning mole fraction condition was fixed 2014.09.30
!
!-------------------------------------------------------------------
! zero all values in equil matrix, dimension (nz1)x(nz1)
    nz2=nz1+1
    smat=zero
    ycondTlimit=.false.
! CCI Bugfixes by Clemnet Introini indicated by CCI    2018.02.20
    evalue=zero
!    dncol=0
!    write(*,*)'in setup_equil: ',converged,nz1,meqrec%tpindep
    if(converged.ge.0) then
! converged < 0 means called from dot derivative, then tcol or pcol set
! otherwise set them to zero
       tcol=0
       pcol=0
       dncol=0
!    else
!       write(*,11)meqrec%nstph,dncol
!11     format('setup: ',10i5)
    endif
!-----------------------------------------------------------
! step 2.1 the Gibbs energies for the stable phases (incl fixed)
    allstableph: do jph=1,meqrec%nstph
       jj=meqrec%stphl(jph)
!       if(meqrec%noofits.le.2) &
!            write(*,12)'pha: ',jph,meqrec%nstph,jj,&
!            phr(jj)%iph,phr(jj)%ics,&
!            phr(jj)%curd%amfu,phr(jj)%curd%gval(1,1)
!12     format(a,5i3,6(1pe12.4))
! column nz2 is the right hand side of the equation, to molar G
       smat(jph,nz2)=phr(jj)%curd%gval(1,1)
!       write(*,313)'Gm: ',0,0,jph,nz2,smat(jph,nz2),ceq%tpval(1)
! one column with amount of component A for each variable chemical potential
! components with fixed chemical potential are automatically skipped
       ncol=1
       xxx=zero
       gloop: do je=1,meqrec%nrel
          do ke=1,meqrec%nfixmu
             if(meqrec%mufixel(ke).eq.je) then
! meqrec%mufixel(ke) is the component number with fix chemical potential
! DONE: reference state must be handelled (may depend on T) ??
!
!---------------------------------------------------------
! handling of user defined reference states for components
                iph=ceq%complist(je)%phlink
                if(iph.gt.0) then
! lokph is index of phase record, to get phase index use phlink ....
!                   iph=ceq%phase_varres(lokph)%phlink
!                   write(*,34)'MM refst: ',je,ke,iph,ceq%complist(je)%endmember
34                 format(a,3i4,4x,10i3)
! we must also handle reference state at fix T !!
                   tpvalsave=ceq%tpval
!                   write(*,*)'MM calling calcg_endmember 2: ',-iph
                   call calcg_endmember(-iph,ceq%complist(je)%endmember,&
                        gref,ceq)
                   if(gx%bmperr.ne.0) then
                      write(*,*)'MM error calculating reference state'
                      ceq%tpval=tpvalsave
                      goto 1000
                   endif
! this is only place where we need to use %mufixvalref
! mufixval should be referred to SER, mufixvalref prescribed value for user ref
                   meqrec%mufixval(ke)=meqrec%mufixvalref(ke)+gref
!                   write(*,35)'MM gref: ',ke,meqrec%mufixvalref(ke),gref,&
!                        meqrec%mufixval(ke)
35                 format(a,i3,6(1pe12.4))
! also copy to cmuval !!?? YES !!!
                   ceq%cmuval(je)=meqrec%mufixval(ke)
!                else
!                   write(*,*)'No userdefined reference state'
                endif
!---------------------------------------------------------
!
                xxx=smat(jph,nz2)
                smat(jph,nz2)=smat(jph,nz2)-&
                     phr(jj)%xmol(je)*meqrec%mufixval(ke)
!                write(*,312)'fix mu G: ',jj,je,ke,xxx,smat(jph,nz2),&
!                     phr(jj)%xmol(je),meqrec%mufixval(ke)
312             format(a,3i3,6(1pe12.4))
                cycle gloop
             endif
          enddo
          smat(jph,ncol)=phr(jj)%xmol(je)
          ncol=ncol+1
       enddo gloop
!       write(*,*)'MM dncol: ',ncol,dncol,meqrec%tpindep
! variable T and P?       
       if(meqrec%tpindep(1)) then
! column for variable T, value is -dG/dT ??
          if(tcol.eq.0) then
             tcol=ncol
             dncol=ncol
             ncol=ncol+1
          endif
          smat(jph,tcol)=-phr(jj)%curd%gval(2,1)
       endif
       if(meqrec%tpindep(2)) then
! column for variable P, value is +dG/dP ??
          if(pcol.eq.0) then
             pcol=ncol
             dncol=ncol
             ncol=ncol+1
          endif
! PVARIABLE in G
          smat(jph,pcol)=-phr(jj)%curd%gval(3,1)
       endif
!       if(meqrec%noofits.le.2) &
!            write(*,13)'Row: ',jph,jj,(smat(jph,je),je=1,nz2)
!13     format(a,2i2,7(1pe10.2))
    enddo allstableph
! we have generated meqrec%nstph rows with ncol columns and rhs in column nz2
! The columns for delta_phase-amounts should be zero
! dncol is number of variable potentials (including T or P if variable)
    if(dncol.eq.0) dncol=ncol-1
!    do iz=1,dncol
!       write(*,228)'smat 1: ',(smat(iz,jz),jz=1,nz2)
!    enddo
!228    format(a,6(1pe12.4))
!    nrow=meqrec%nstph
!-------------------------------------------------------------------
! step 2.2 equations due to user conditions on extensive/normalizzed properties
! nz2 is number of columns, last column is right hand side (rhs)
! nrow is number of nows already filled (G for stable ph)
!    nz2=nz1+1
!
! >>>>>>>>>>> THIS IS UNFINISHED, ONLY A FEW STATE VARIABLES ALLOWED
! expressions only for N and x and H ... added V mm ... y 190720
!
    nrow=meqrec%nstph
    lastcond=>ceq%lastcondition
    condition=>lastcond
    allocate(mamu(meqrec%nrel),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 17: ',errall
       gx%bmperr=4370; goto 1000
    endif
! for saving partial dgdyterms, set nosave=.TRUE. to use old calc_dgdyterms1
!    nosave=.TRUE.
! nosave always FALSE as there are places to save results in phase_varres
    nosave=.FALSE.
!    savedrec%sameit=0
!    saved=>savedrec
350 continue
! cmode=0 means calculate and return current value
    cmode=0
    cmix=0
    condition=>condition%next
! This is the condition, cvalue is the prescibed value
! cmode and cmix contain information how to calculate its current value
!    write(*,*)'MM calling apply',condition%noofterms
! apply_condition in gtp3X.F90 ??
    call apply_condition_value(condition,cmode,cvalue,cmix,ccf,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,71)'MM apply 2: ',cmode,cvalue,cmix,ccf(1)
71  format(a,i3,1pe12.4,22i4,1pe12.4)
!    if(condition%noofterms.gt.1) write(*,351)nrow,cmode,cmix,nterms,cvalue,&
!         (ccf(jj),jj=1,condition%noofterms)
! Only cmix(1)=5 is interesting here. potentials already cared for
    if(cmix(1).ne.5) then
! loop if not the last condition
!       write(*,*)'Taking next condition: ',cmix(1)
       if(.not.associated(condition,lastcond)) goto 350
       goto 380
    endif
! check if several terms
    mterms=1
    nterms=condition%noofterms
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
! Q       10x   "                            14     Internal stability
! N       11x  (component/phase#set,component) 16  moles of components
! X       111   "                            17     mole fraction of components
! B       12x   "                            18     mass of components
! W       122   "                            19     mass fraction of components
! Y       13    phase#set,constituent#subl   20     constituent fraction
!----- model variables <<<< these now treated differently
    stvix=cmix(2)/10
! stvnorm is normalization, 0, 1, 2, 3 or 4
! 0=none; 1=per mole; 2=per mass; 3=per volume; 4=per formula unit
    stvnorm=mod(cmix(2),10)
    select case(stvix)
    case default
       write(*,*)'not a condition:',stvix,stvnorm,cmix(1),cmix(2),cmix(3)
       gx%bmperr=4208; goto 1000
    case(1,5) 
! stvix=1..6: U, S, V, H, A, G, some conditions not implemented
!             1  2  3  4  5  6
       write(*,*)'Not implemented yet: ',stvix,stvnorm
       gx%bmperr=4207; goto 1000
!------------------------------------------------------------------
! Entropy for system or phase(s)
    case(2) ! S entropy condition
       write(*,*)'MM entropy condition testing: ',nterms,nrow,nz1
       if(stvnorm.eq.0) then
! not normallized          
          if(cmix(3).eq.0) then
! condition is S=value
             sph=0
          else
! condition is S(phase#set)=value
             sph=cmix(3); scs=cmix(4)
          endif
          write(*,*)'MM not normallized entropy conditions not implemented'
          gx%bmperr=4207; goto 1000
       else
! entropy difference: to use the condition SM(solid)-SM(liquid)=0
! for equientropy lines ...
! s1-s2=0: delta-s = ds/dT dT + ds/dy dy + ... = 0
          xterm=3
          xxx=zero
! calculate and store the drivatives in xcol
! How to know wich is the independent for each index?
! SEE HOW A V condition is calculated below!!
          allocate(xcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 18: ',errall
             gx%bmperr=4370; goto 1000
          endif
220       continue
          if(mterms.le.nterms) then
! loop over ALL phases
             sph=cmix(xterm); scs=cmix(xterm+1)
             do jph=1,meqrec%nphase
                if(phr(jph)%iph.eq.sph .and. phr(jph)%ics.eq.scs) then
! extract the value of SM for phase in mterms
!                   write(*,*)'MM: sph,scs: ',mterms,xterm,sph,scs
                   xxx=xxx+ccf(mterms)*phr(sph)%curd%gval(2,1)/&
                        phr(sph)%curd%abnorm(1)
                   write(*,230)'MM S: ',ccf(mterms),phr(sph)%curd%gval(2,1),&
                        phr(sph)%curd%abnorm(1),phr(sph)%curd%amfu,xxx
230                format(a,6(1pe12.4))
                   xterm=xterm+4
                   mterms=mterms+1
                   goto 220
                endif
             enddo
             write(*,*)'MM cannot find phase for EEC',mterms
             gx%bmperr=4399; goto 1000
          endif
       endif
       nrow=nrow+1
       write(*,230)'MM equientropy: ',ceq%tpval(1),cvalue,xxx
       smat(nrow,1)=xxx
!       gx%bmperr=4207; goto 1000
!------------------------------------------------------------------
    case(3) ! V volume condition, almost the same a H condition
! Volume for system or phase, NOT normallized
       if(stvnorm.eq.0) then
! not normallized
          if(cmix(3).eq.0) then
! condition is V=value
             sph=0
          else
! condition is V(phase#set)=value
             sph=cmix(3); scs=cmix(4)
          endif
! FU(alpha) is formula units of alpha phase, V=\sum_alpha VM(alpha) VM(alpha)
! dVM(alpha) = d2GM/dPdy_i*c_iA*\mu_A+
!     \sum_i dGM/dP*dP + ??
!     \sum_alpha ???
! UNFINISHED ??
          allocate(xcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 19: ',errall
             gx%bmperr=4370; goto 1000
          endif
          xcol=zero
          totam=zero
          notf=0
          check1=zero
          check2=zero
          notdone=.TRUE.
          vallph: do jph=1,meqrec%nstph
! sum over all stable phases
             jj=meqrec%stphl(jph)
             pmi=>phr(jj)
! if phase is not fixed there is a column in xcol for variable amount
! This has to be done before loop of elements
             if(pmi%phasestatus.ne.PHFIXED) notf=notf+1
             if(sph.gt.0) then
! if a phase is specified, skip all other phases
                if(.not.(sph.eq.phr(jj)%iph .and. scs.eq.phr(jj)%ics)) &
                     cycle vallph
             endif
! moles formula unit of phase
             pham=pmi%curd%amfu
             allocate(hval(pmi%ncc),stat=errall)
             if(errall.ne.0) then
                write(*,*)'MM Allocation error 20: ',errall
                gx%bmperr=4370; goto 1000
             endif
             notdone=.FALSE.
             if(.not.allocated(mamu1)) then
! it will be deallocated when leaving this subroutine ??
                allocate(mamu1((meqrec%nrel)),stat=errall)
                if(errall.ne.0) then
                   write(*,*)'MM Allocation error 21: ',errall
                   gx%bmperr=4370; goto 1000
                endif
             endif
             ncol=1
             if(stvix.eq.3) then
! V condition, calculate the terms d2G/dPdy_i for all constituents
                do ie=1,pmi%ncc
                   hval(ie)=pmi%curd%dgval(3,ie,1)
                enddo
!                write(*,*)'Volume condition: ',pcol,pmi%ncc,hval(1)
             endif
!             write(*,75)'hval: ',hval
!             write(*,75)'cmuvamanyl: ',(ceq%cmuval(ie),ie=1,meqrec%nrel)
! calculate the terms to be multiplied with the unknown mu(ie)
             vallel: do ie=1,meqrec%nrel
! multiply terms with the inverse phase matrix and hval()
! but also return values without this in mamu1,mag1,mat1 and map1 needed
! for normalization and if there is a condition on chemical potentials
                call calc_dgdytermshm(meqrec%nrel,ie,meqrec%tpindep,hval,&
                     mamu,mag,mat,map,mamu1,mag1,mat1,map1,&
                     pmi,ceq%cmuval,meqrec%noofits)
                if(gx%bmperr.ne.0) goto 1000
! calculate a term for each column to be multiplied with chemical potential
! if the potential is fixed add the term to the rhs
                do ke=1,meqrec%nfixmu
                   if(meqrec%mufixel(ke).eq.ie) then
! components with fix chemical potential added to rhs, do not increment ncol!!!
                      xcol(nz2)=xcol(nz2) + pham*meqrec%mufixval(ke)*mamu(ie)
!                      write(*,102)'fix mu V:',nz2,ie,pham,&
!                           meqrec%mufixval(ke),mamu1(ie),&
!                           pham*meqrec%mufixval(ke)*mamu1(ie),xcol(nz2)
                      cycle vallel
                   endif
                enddo
                xcol(ncol)=xcol(ncol) - pham*mamu(ie)
                ncol=ncol+1
             enddo vallel
! vallel loop should end here as mat and map are element independent
! If T or P are variable, mat and map include \sum_j hval(j)
             if(tcol.gt.0) then
                xxx=xcol(tcol)
! gval(2,1) is dG/dT, gval(4,1) is d2G/dT2, gval(5,1) is d2G/dTdP=dV/dT
                xcol(tcol)=xcol(tcol)+&
                     2.0D-3*pham*(ceq%tpval(1)*pmi%curd%gval(5,1)-mat)
! Why is d2G/dTdP multiplied by T??
! >500 its            1.0D-3*pham*(ceq%tpval(1)*pmi%curd%gval(5,1)-mat)
! 27 its              2.0D-3*pham*(ceq%tpval(1)*pmi%curd%gval(5,1)-mat)
! 80 its              5.0D-3*pham*(ceq%tpval(1)*pmi%curd%gval(5,1)-mat)
! 158 its             1.0D-2*pham*(ceq%tpval(1)*pmi%curd%gval(5,1)-mat)
! slow                 pham*(ceq%tpval(1)*pmi%curd%gval(5,1)-mat)
! bad                  pham*ceq%tpval(1)*(ceq%tpval(1)*pmi%curd%gval(5,1)-mat)
! slow                 pham*(ceq%tpval(1)*pmi%curd%gval(5,1)-mat)
! wrong                pham*(mat-ceq%tpval(1)*pmi%curd%gval(5,1))
!                write(*,*)'VCONDT: ',tcol,xcol(tcol)
             endif
! PVARIABLE for condition on V
             if(pcol.gt.0) then
                xxx=xcol(pcol)
! gval(3,1) is dG/dP, gval(6,1) is d2G/dP2, sign???
                xcol(pcol)=xcol(pcol)+pham*(pmi%curd%gval(6,1)-map)
!                xcol(pcol)=xcol(pcol)+pham*(map-pmi%curd%gval(6,1))
!                     pmi%curd%gval(3,1)-ceq%tpval(1)*pmi%curd%gval(5,1))
!                write(*,*)'VCONDP: ',pcol,xcol(pcol)
             endif
! uncertain if enddo hallel here or after label 7000 above ...
             deallocate(hval)
             if(stvix.eq.3) then
! sum the total volume (or for a single phase its volume)
! slow                totam=totam+pham*pmi%curd%gval(3,1)
!                totam=totam+pham*pmi%curd%gval(3,1)
                totam=totam+pham*pmi%curd%gval(3,1)
! wrong                totam=totam+pham*pmi%curd%gval(3,1)*ceq%rtn
             endif
! Now the term multipled with change of the amount of the phase
             if(pmi%phasestatus.ne.PHFIXED) then
                xcol(dncol+notf)=pmi%curd%gval(3,1)
             endif
! term to the RHS, sign???
             xcol(nz2)=xcol(nz2)+pham*mag
! slow            xcol(nz2)=xcol(nz2)+pham*mag
! as slow         xcol(nz2)=xcol(nz2)-pham*mag
          enddo vallph
          if(sph.gt.0 .and. notdone) then
! if sph.ne.0 it is possible that the specified phase is not stable, check that
! the vallph loop has beed done at least once
             write(*,*)'Unnormalized volume condition of unstable phase'
! These values are most probably all zero making system matrix singular
             write(*,177)'xcol: ',nz2,(xcol(jj),jj=1,nz2)
             gx%bmperr=4196; goto 1000
          endif
! Add difference to the RHS.  Totam is summed above, cvalue is prescribed value
!          write(*,74)'Volume: ',nrow+1,ceq%tpval(1),ceq%rtn,&
!               xcol(nz2),totam,cvalue/ceq%rtn
! sign?   xcol(nz2)=xcol(nz2)+totam-cvalue/ceq%rtn
          xcol(nz2)=xcol(nz2)-totam+cvalue/ceq%rtn
!          xcol(nz2)=xcol(nz2)-totam+cvalue
!          write(*,75)'RHS: ',xcol(nz2),totam,cvalue/ceq%rtn,ceq%rtn,&
!               totam*ceq%rtn,ceq%tpval(1)
! test if condition converged, use relative error 
! slow          if(abs(totam-cvalue/ceq%rtn).gt.ceq%xconv*abs(cvalue)) then
          if(abs(totam-cvalue/ceq%rtn).gt.ceq%xconv*abs(cvalue)) then
!          if(abs(totam-cvalue).gt.ceq%xconv*abs(cvalue)) then
!                  write(*,75)'Unconverged volume: ',ceq%tpval(1),&
!             if(vbug) write(*,75)'Unconverged volume: ',ceq%tpval(1),&
!             write(*,75)'Unconverged volume: ',ceq%tpval(1),&
!                  totam,cvalue,totam-cvalue,ceq%xconv*abs(cvalue)
!                  totam,cvalue/ceq%rtn,totam-cvalue/ceq%rtn
             if(converged.lt.5) then
                converged=5
!                write(*,*)'2: converged=5',cerr%nvs
                cerr%mconverged=converged
                if(cerr%nvs.lt.10) then
                   cerr%nvs=cerr%nvs+1
                   cerr%typ(cerr%nvs)=5
                   cerr%val(cerr%nvs)=cvalue
                   cerr%dif(cerr%nvs)=totam-cvalue
                endif
             endif
          endif
! we have one more equation to add to the equilibrium matrix
          nrow=nrow+1
          if(nrow.gt.nz1) stop 'MM too many equations 5A'
          do ncol=1,nz2
             smat(nrow,ncol)=xcol(ncol)
          enddo
          deallocate(xcol)
       else
! volume is normalized
          write(*,*)'Normalized volume condition not implemented yet'
          gx%bmperr=4207; goto 1000
       endif
!------------------------------------------------------------------
    case(4) ! Enthaly condition (Heat balance). 
! Enthalpy for system or phase, normallized or not
!       gx%bmperr=4207; goto 1000
       if(stvnorm.eq.0) then
! not normallized
          if(cmix(3).eq.0) then
! condition is H=value or V=value
             sph=0
          else
! condition is H(phase#set)=value or V(phase#set)=value
             sph=cmix(3); scs=cmix(4)
          endif
! FU(alpha) is formula units of alpha phase
! dH=\sum_alpha FU(alpha)(dG/y_i-Td2G/dTdy_i)*c_iA*\mu_A + 
!   (-Td2G/dT2 + \sum_i (dG/dy_i - Td2G/dTdY_i)*c_iT)*dT + ...
!   +\sum_alpha (G-TdG/dT)*\delta FU(alpha) =
!    \sum_alpha FU(alpha)\sum_i(dG/dy_i-Td2G/dTdy_i)*c_iG + H\tilde - H
!          write(*,*)'Condition on H: ',pmi%ncc,dncol
! dV = \sum_alpha FU(alpha)(d2G/dPdy_i)*c_iA*\mu_A+
!     \sum_i dG/dP*dP + ??
!     \sum_alpha ???
! Condition H=value and H(phase)=value are OK, HM=value is NOT OK  Why??
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          allocate(xcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 22: ',errall
             gx%bmperr=4370; goto 1000
          endif
          xcol=zero
          totam=zero
          notf=0
          check1=zero
          check2=zero
          notdone=.TRUE.
          hallph: do jph=1,meqrec%nstph
! sum over all stable phases
             jj=meqrec%stphl(jph)
             pmi=>phr(jj)
! if phase is not fixed there is a column in xcol for variable amount
! This has to be done before loop of elements
             if(pmi%phasestatus.ne.PHFIXED) notf=notf+1
             if(sph.gt.0) then
! if a phase is specified, skip all other phases
                if(.not.(sph.eq.phr(jj)%iph .and. scs.eq.phr(jj)%ics)) &
                     cycle hallph
             endif
! moles formula unit of phase
             pham=pmi%curd%amfu
             allocate(hval(pmi%ncc),stat=errall)
             notdone=.FALSE.
             if(.not.allocated(mamu1)) then
! it will be deallocated when leaving this subroutine ??
                allocate(mamu1((meqrec%nrel)),stat=errall)
             endif
             if(errall.ne.0) then
                write(*,*)'MM Allocation error 23: ',errall
                gx%bmperr=4370; goto 1000
             endif
             ncol=1
             if(stvix.eq.3) then
! V condition, calculate the terms d2G/dPdy_i for all constituents
! THIS IS REDUNDANT, V HAS ITIS OWN CASE NOW
                do ie=1,pmi%ncc
                   hval(ie)=pmi%curd%dgval(3,ie,1)
                enddo
!                write(*,*)'Volume condition: ',pcol,pmi%ncc,hval(1)
             else
! H condition, calculate the terms dG/dy_i - T*d2G/dTdy_i for all constituents
                do ie=1,pmi%ncc
                   hval(ie)=pmi%curd%dgval(1,ie,1)-&
                        ceq%tpval(1)*pmi%curd%dgval(2,ie,1)
                enddo
!                write(*,*)'Enthalpy condition: ',tcol,hval(1)
             endif
!             write(*,75)'hval: ',hval
!             write(*,75)'cmuvamanyl: ',(ceq%cmuval(ie),ie=1,meqrec%nrel)
! calculate the terms to be multiplied with the unknown mu(ie)
             hallel: do ie=1,meqrec%nrel
! multiply terms with the inverse phase matrix and hval()
! but also return values without this in mamu1,mag1,mat1 and map1 needed
! for normalization and if there is a condition on chemical potentials
                call calc_dgdytermshm(meqrec%nrel,ie,meqrec%tpindep,hval,&
                     mamu,mag,mat,map,mamu1,mag1,mat1,map1,&
                     pmi,ceq%cmuval,meqrec%noofits)
                if(gx%bmperr.ne.0) goto 1000
!                write(*,99)'hfix 1: ',ceq%tpval(1),mag,mat,map,mamu
99              format(a,6(1pe12.4))
! calculate a term for each column to be multiplied with chemical potential
! if the potential is fixed add the term to the rhs
                do ke=1,meqrec%nfixmu
                   if(meqrec%mufixel(ke).eq.ie) then
! components with fix chemical potential added to rhs, do not increment ncol!!!
                      xcol(nz2)=xcol(nz2) + pham*meqrec%mufixval(ke)*mamu(ie)
!                      write(*,102)'fix mu H6:',nz2,ie,pham,&
!                           meqrec%mufixval(ke),mamu1(ie),&
!                           pham*meqrec%mufixval(ke)*mamu1(ie),xcol(nz2)
102                   format(a,2i3,6(1pe12.4))
                      cycle hallel
                   endif
                enddo
                xcol(ncol)=xcol(ncol) - pham*mamu(ie)
                ncol=ncol+1
             enddo hallel
! I think hallel loop should end here as mat and map are element independent
! If T or P are variable, mat and map include \sum_j hval(j)
             if(tcol.gt.0) then
                xxx=xcol(tcol)
! gval(2,1) is dG/dT, gval(4,1) is d2G/dT2, sign????
                xcol(tcol)=xcol(tcol)+&
                     pham*(ceq%tpval(1)*pmi%curd%gval(4,1)-mat)
             endif
! PVARIABLE condition on H
             if(pcol.gt.0) then
                xxx=xcol(pcol)
! gval(3,1) is dG/dP, gval(5,1) is d2G/dTdP, sign???
                xcol(pcol)=xcol(pcol)+pham*(pmi%curd%gval(3,1)-map)
!                xcol(pcol)=xcol(pcol)+pham*(map-&
!                     pmi%curd%gval(3,1)-ceq%tpval(1)*pmi%curd%gval(5,1))
!>>                xcol(pcol)=xcol(pcol)-pham*(map-&
!                     pmi%curd%gval(3,1)+ceq%tpval(1)*pmi%curd%gval(5,1))
!                write(*,363)'d2G/dPdy: H',nrow+1,ie,pcol,&
!                     xxx,xcol(pcol),pham,mat
             endif
! uncertain if enddo hallel here or after label 7000 above ...
!             enddo hallel
! hval no longer needed
             deallocate(hval)
             if(stvix.eq.3) then
! sum the total volune (or for a single phase its volume)
                totam=totam+pham*pmi%curd%gval(3,1)
!                write(*,211)'HMS total volume:',totam,ceq%rtn*totam,cvalue
211             format(a,5(1pe12.4))
             else
! Sum the total enthalpy (for a single phase just one value)
                totam=totam+pham*(pmi%curd%gval(1,1)-&
                     ceq%tpval(1)*pmi%curd%gval(2,1))
!             write(*,73)'pham:  ',sph,jj,pham,totam,ceq%cmuval(1),ceq%cmuval(2)
             endif
! Now the term multipled with change of the amount of the phase
             if(pmi%phasestatus.ne.PHFIXED) then
                xcol(dncol+notf)=pmi%curd%gval(1,1)-&
                     ceq%tpval(1)*pmi%curd%gval(2,1)
             endif
! term to the RHS, sign???
!             xcol(nz2)=xcol(nz2)-pham*mag
             xcol(nz2)=xcol(nz2)+pham*mag
!             write(*,76)'Check2: ',jj,pham,mag,pham*mag
          enddo hallph
          if(sph.gt.0 .and. notdone) then
! if sph.ne.0 it is possible that the specified phase is not stable, check that
! the hallph loop has beed done at least once
             write(*,*)'Unnormalized enthalpy condition of unstable phase'
! These values are most probably all zero making system matrix singular
             write(*,177)'xcol: ',nz2,(xcol(jj),jj=1,nz2)
177          format(a,i2,6(1pe10.2))
             gx%bmperr=4196; goto 1000
          endif
!          write(*,177)'xcol: ',nz2,(xcol(jj),jj=1,nz2)
! Add difference to the RHS.  Totam is summed above, cvalue is prescribed value
!          write(*,74)'Enthalpy: ',nrow+1,ceq%tpval(1),ceq%rtn,&
!               xcol(nz2),totam,cvalue/ceq%rtn
          xcol(nz2)=xcol(nz2)+totam-cvalue/ceq%rtn
!          write(*,75)'RHS: ',xcol(nz2),totam,cvalue,ceq%rtn,cvalue/ceq%rtn
! test if condition converged, use relative error 
          if(abs(totam-cvalue/ceq%rtn).gt.ceq%xconv*abs(cvalue)) then
!             write(*,75)'Unconverged enthalpy: ',ceq%tpval(1),&
!                  totam,cvalue/ceq%rtn,totam-cvalue/ceq%rtn
             if(converged.lt.5) then
                converged=5
!                write(*,*)'3: converged=5',cerr%nvs
                cerr%mconverged=converged
                if(cerr%nvs.lt.10) then
                   cerr%nvs=cerr%nvs+1
                   cerr%typ(cerr%nvs)=5
                   cerr%val(cerr%nvs)=cvalue/ceq%rtn
                   cerr%dif(cerr%nvs)=totam-cvalue/ceq%rtn
                endif
             endif
          endif
! we have one more equation to add to the equilibrium matrix
          nrow=nrow+1
          if(nrow.gt.nz1) stop 'MM too many equations 5A'
          do ncol=1,nz2
             smat(nrow,ncol)=xcol(ncol)
          enddo
!          write(*,*)'H conv: ',ceq%tpval(1)
!          write(*,74)'hline: ',nrow,xcol
75        format(a,6(1pe12.4))
74        format(a,i2,6(1pe11.3))
73        format(a,2i3,6(1pe11.3))
! check1 and check2 should be equal if we set H as current value and release T
!          write(*,75)'Check: ',check1,check2
          deallocate(xcol)
! ..........................................................
       else
! normallized HM (per mole, 1), HW (per mass, 2) or HV (per volume, 3)
!          write(*,*)'*** Normallized enthalpy not yet implemented as condition'
!          gx%bmperr=4207; goto 1000
! UNFINISHED 
          if(stvnorm.ne.1) then
             write(*,*)'Only normallizing per mole implemented'
             gx%bmperr=4207; goto 1000
          endif
! ie=0 means no element specification
          ie=0
          if(cmix(3).eq.0) then
! condition is HM=value
             sph=0
          else
! condition is HM(phase#set)=value
! UNFINISHED, does not converge 
!             gx%bmperr=4207; goto 1000
             sph=cmix(3); scs=cmix(4)
          endif
! dH=\sum_alpha FU(alpha)(dG/y_i-Td2G/dTdy_i)c_iA\mu_A + 
!   (-Td2G/dT2 + \sum_i (dG/dy_i - Td2G/dTdY_i)c_iT)dT + ...
!   +\sum_alpha (G-TdG/dT)\delta FU(alpha) =
!    \sum_alpha FU(alpha)\sum_i(dG/dy_i-Td2G/dTdy_i)c_iG + H-\tilde H
          allocate(xcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 24: ',errall
             gx%bmperr=4370; goto 1000
          endif
          xcol=zero
          totam=zero
          notf=0
          check1=zero
          check2=zero
          notdone=.TRUE.
          if(.not.allocated(mamu1)) then
! it will be deallocated when leaving this subroutine ??
             allocate(mamu1((meqrec%nrel)),stat=errall)
             if(errall.ne.0) then
                write(*,*)'MM Allocation error 25: ',errall
                gx%bmperr=4370; goto 1000
             endif
          endif
! current value of molar enthalpy
          if(sph.eq.0) then
             call get_state_var_value('HM ',hmval,encoded,ceq)
             totalmol=one
          else
! current value of molare enthalpy for a phase
             call get_phase_name(sph,scs,encoded)
             name='HM('//encoded
             jj=len_trim(name)
             name(jj+1:)=')'
             call get_state_var_value(name,hmval,encoded,ceq)
          endif
          call get_state_var_value('N ',totalmol,encoded,ceq)
          hmval=hmval/ceq%rtn
! this is not yet implemented
          write(*,*)'hmval, totalmol: ',hmval, totalmol
          if(gx%bmperr.ne.0) goto 1000
          hmallph: do jph=1,meqrec%nstph
! sum over all stable phases
             jj=meqrec%stphl(jph)
             pmi=>phr(jj)
! if phase is not fixed there is a column in xcol for variable amount
! This has to be done before loop of elements
             if(pmi%phasestatus.ne.PHFIXED) notf=notf+1
             if(sph.gt.0) then
! if a phase is specified, skip all other phases
                if(.not.(sph.eq.phr(jj)%iph .and. scs.eq.phr(jj)%ics)) &
                     cycle hmallph
                pham=one
             else
                pham=pmi%curd%amfu
             endif
! moles formula unit of phase
             allocate(hval(pmi%ncc),stat=errall)
             if(errall.ne.0) then
                write(*,*)'MM Allocation error 26: ',errall
                gx%bmperr=4370; goto 1000
             endif
             notdone=.FALSE.
! calculate the terms dG/dy_i - T*d2G/dTdy_i for all constituents
             do ie=1,pmi%ncc
                hval(ie)=pmi%curd%dgval(1,ie,1)-&
                     ceq%tpval(1)*pmi%curd%dgval(2,ie,1)
             enddo
             write(*,73)'hmval: ',sph,ie,hmval
!             write(*,75)'cmuvamanyl: ',(ceq%cmuval(ie),ie=1,meqrec%nrel)
! ncol is increemented for each variable chemical potential
             ncol=1
! calculate the terms to be multiplied with the unknown mu(ie)
             hmallel: do ie=1,meqrec%nrel
! multiply terms with the inverse phase matrix and hval
! but also return values without this in mamu1,mag1,mat1 and map1 needed
! for normalization ...
                call calc_dgdytermshm(meqrec%nrel,ie,meqrec%tpindep,hval,&
                     mamu,mag,mat,map,mamu1,mag1,mat1,map1,&
                     pmi,ceq%cmuval,meqrec%noofits)
                if(gx%bmperr.ne.0) goto 1000
! In this loop we subtract H/N*\sum_B \Delta M_B for all terms
                ncol2=1
                hmloop1: do je=1,meqrec%nrel
                   do ke=1,meqrec%nfixmu
                      if(meqrec%mufixel(ke).eq.je) then
! components with fix chemical potential added to rhs, do not increment ncol2!!!
                         xcol(nz2)=xcol(nz2)+&
                              pham*hmval*mamu1(je)*meqrec%mufixval(ke)
!                         write(*,102)'fix mu 1: ',nz2,je,pham,mamu1(je),&
!                              meqrec%mufixval(ke)
                         cycle hmloop1
                      endif
                   enddo
! mamu(B) = \sum_i \sum_j \sum_A dM^a_B/dy_i dM^a_A z^a_ij
                   xcol(ncol2)=xcol(ncol2)-pham*hmval*mamu1(je)
!                   write(*,102)'HM jel:',je,ncol2,pham,&
!                        mamu(je),hmval,mamu1(je),xcol(ncol2)
                   ncol2=ncol2+1
                enddo hmloop1
! calculate a term for each column to be multiplied with chemical potential
! if the potential is fixed add the term to the rhs
                do ke=1,meqrec%nfixmu
                   if(meqrec%mufixel(ke).eq.ie) then
! components with fix chemical potential added to rhs, do not increment ncol!!!
                      xcol(nz2)=xcol(nz2) + pham*meqrec%mufixval(ke)*mamu(ie)
!                      write(*,102)'fix mu HM 3:',nz2,ke,pham,&
!                           meqrec%mufixval(ke),mamu1(ie),xcol(nz2)
                      cycle hmallel
                   endif
                enddo
! mamu(ie) = \sum_i hval(i) \sum_j \sum_B dM^a_B/dy_j z^a_ij
                xcol(ncol)=xcol(ncol) - pham*mamu(ie)
!                write(*,102)'HM col:',ie,ncol,pham,mamu(ie),xcol(ncol)
                ncol=ncol+1
!                check1=check1-pham*mamu(ie)*ceq%cmuval(ie)
!                write(*,76)'check1: ',ie,check1,pham*mamu(ie)*ceq%cmuval(ie)
76              format(a,i2,6(1pe12.4))
             enddo hmallel
! UNFINSHED: problems converging with normallized enthalpy condition 
! If T or P are variable, mat and map include \sum_j hval(j)
             if(tcol.gt.0) then
                xxx=xcol(tcol)
! gval(2,1) is dG/dT, gval(4,1) is d2G/dT2, sign????
! the equation above should be better but ....
                xcol(tcol)=xcol(tcol)+&
                     pham*(ceq%tpval(1)*pmi%curd%gval(4,1)-mat)
!                     pham*(ceq%tpval(1)*pmi%curd%gval(4,1)-mat+hmval*mat1)
!                write(*,102)'HM dt: ',0,tcol,pham,&
!                     ceq%tpval(1)*pmi%curd%gval(4,1),mat,hmval,mat1,xcol(tcol)
             endif
             if(pcol.gt.0) then
! condition on H and variable P
                xxx=xcol(pcol)
! gval(3,1) is dG/dP, gval(5,1) is d2G/dTdP, sign??? UNFINISHED TEST
                xcol(pcol)=xcol(pcol)+pham*(map-hmval*map1-&
                     pmi%curd%gval(3,1)-ceq%tpval(1)*pmi%curd%gval(5,1))
             endif
! Now the term multipled with change of the amount of the phase, not pham
             if(pmi%phasestatus.ne.PHFIXED) then
                xcol(dncol+notf)=xcol(dncol+notf)+pmi%curd%gval(1,1)-&
                     ceq%tpval(1)*pmi%curd%gval(2,1)
!                     ceq%tpval(1)*pmi%curd%gval(2,1)-hmval
!                write(*,102)'HM dn: ',ie,dncol+notf,0.0,&
!                     pmi%curd%gval(1,1)-ceq%tpval(1)*pmi%curd%gval(2,1),&
!                     hmval,xcol(dncol+notf)
             endif
! term to the RHS
!             xcol(nz2)=xcol(nz2)+pham*(mag-hmval*mag1)
             xcol(nz2)=xcol(nz2)+pham*mag
!             write(*,102)'HM rhs:',ie,nz2,pham,mag,hmval,mag1,xcol(nz2)
! hval can be differnt for next phase
             deallocate(hval)
          enddo hmallph
          if(sph.gt.0 .and. notdone) then
! if sph.ne.0 it is possible that the specified phase is not stable, check that
! the hallph loop has beed done at least once
             write(*,*)'Normalized enthalpy condition of unstable phase'
! These values are most probably all zero making system matrix singular
             write(*,177)'xcol: ',nz2,(xcol(jj),jj=1,nz2)
             gx%bmperr=4196; goto 1000
          endif
!          write(*,177)'xcol: ',nz2,(xcol(jj),jj=1,nz2)
! Add difference to the RHS.  Totam is summed above, cvalue is prescribed value
!          write(*,74)'Enthalpy: ',nrow+1,ceq%tpval(1),ceq%rtn,&
!               xcol(nz2),totam,cvalue/ceq%rtn
!          xcol(nz2)=xcol(nz2)+totam-cvalue/ceq%rtn
          xcol(nz2)=xcol(nz2)/totalmol-hmval+cvalue/ceq%rtn
!          write(*,75)'RHS: ',xcol(nz2),hmval,cvalue/ceq%rtn,totalmol,&
!               ceq%tpval(1)
! test if condition converged, use relative error 
          if(abs(hmval-cvalue/ceq%rtn).gt.ceq%xconv*abs(cvalue)) then
             write(*,75)'Unconverged enthalpy: ',&
                  hmval*ceq%rtn,cvalue,hmval-cvalue/ceq%rtn
             if(converged.lt.5) then
                converged=5 
!                write(*,*)'4: converged=5',cerr%nvs
                cerr%mconverged=converged
                if(cerr%nvs.lt.10) then
                   cerr%nvs=cerr%nvs+1
                   cerr%typ(cerr%nvs)=5
                   cerr%val(cerr%nvs)=hmval
                   cerr%dif(cerr%nvs)=hmval-cvalue/ceq%rtn
                endif
             endif
          endif
! we have one more equation to add to the equilibrium matrix
          nrow=nrow+1
          if(nrow.gt.nz1) stop 'MM too many equations 5B'
! we must divide all terms in the LHS with totalmol
          do ncol=1,nz1
             smat(nrow,ncol)=xcol(ncol)/totalmol
          enddo
          smat(nrow,nz2)=xcol(nz2)
!          write(*,*)'H conv: ',ceq%tpval(1)
!          write(*,74)'hline: ',nrow,xcol
! check1 and check2 should be equal if we set H as current value and release T
!          write(*,75)'Check: ',check1,check2
          deallocate(xcol)
       endif
! already calculated above
!------------------------------------------------------------------
    case(6) ! G
! Gibbs energy, for system or a phase
       gx%bmperr=4207; goto 1000
       if(stvnorm.eq.0) then
! not normallized
          if(cmix(3).eq.0) then
! condition is G=value
             sph=0
          else
! condition is G(phase#set)=value
             gx%bmperr=4207; goto 1000
             sph=cmix(3); scs=cmix(4)
          endif
! current value of dG=\sum_A dM_A \mu_A + G -\tilde G=0
          allocate(xcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 27: ',errall
             gx%bmperr=4370; goto 1000
          endif
          xcol=zero
!...UNFINISHED
          gx%bmperr=4207; goto 1000
          nrow=nrow+1
          if(nrow.gt.nz1) stop 'MM too many equations 6A'
          do ncol=1,nz2
             smat(nrow,ncol)=xcol(ncol)
          enddo
! set rhs to G^prescribed - G^current 
!          smat(nrow,nz2)=cvalue
          deallocate(xcol)
       else
! normallizing can be M (per mole, 1), W (per mass, 2) or V (per volume, 3?)
          gx%bmperr=4207; goto 1000
       endif
!------------------------------------------------------------------
    case(7) ! NP
! Amount of phase in moles, use fix phase instead
       write(*,352)stvix,stvnorm
352    format('Not implemented yet, use set status phase=fix: ',2i5)
       gx%bmperr=4207; goto 1000
       nrow=nrow+1
       if(nrow.gt.nz1) stop 'MM too many equations 7A'
!------------------------------------------------------------------
    case(8) ! BP
! Amount of phase in mass, use fix phase instead
       write(*,352)stvix,stvnorm
       gx%bmperr=4207; goto 1000
       nrow=nrow+1
       if(nrow.gt.nz1) stop 'MM too many equations 8A'
!------------------------------------------------------------------
! 9 and 10 (DG and Q) not allowed as conditions
!------------------------------------------------------------------
    case(11) ! N or X with or without indices and normalization
! 160818: adding possibility to have several terms a*N(A)-b*N(B)=cvalue
1100   continue
       if(stvnorm.eq.0) then
          moffs=0
!          write(*,*)'MM condition for N: ',nterms,sph,sel
! return here for second term
1107      continue
          if(cmix(3).eq.0) then
! condition is N=fix
             sel=0; sph=0
          elseif(cmix(4+moffs).eq.0) then
! condition is N(A)=fix
             sel=cmix(3+moffs); sph=0
          else
! condition is N(phase#set,A)=fix;  how to handle if phase#set not stable?
!             write(*,*)'Condition N(phase#set,A)=fix not allowed'
!             gx%bmperr=4208; goto 1000
             sel=cmix(5+moffs); sph=cmix(3+moffs); scs=cmix(4+moffs)
          endif
!          write(*,*)'Condition on N, N(A) or N(phase,A)',sph,sel
! Formulate equation for total amount N:
! rhs:  N-N+\sum_alpha N^a + \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j
! \sum_B \sum_alpha N^a \sum_i \sum_j dM^_A/dy_i dM^a_B/dy_j*z^a_ij  *mu(B)
!        \sum_alpha N^a \sum_i d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j      *deltaT
!        \sum_alpha N^a \sum_i d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j      *deltaP
!        \sum_A M^a_A                                    *deltaN^a
          allocate(xcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 28: ',errall
             gx%bmperr=4370; goto 1000
          endif
          xcol=zero
          totam=zero
! notf keeps track on entered non-fixed phases with variable amount
          notf=0
! THE CALCULATION FOR N= and N(A)= seems OK
! sum over all phases to handle conditions like N(phase#set,A)=fix
! as the phase#set may not be stable
!          write(*,*)'Loop for all all phases for condition N='
          nallph: do jj=1,meqrec%nphase
             pmi=>phr(jj)
             if(sph.eq.0) then
! skip if not stable
                if(phr(jj)%stable.eq.0) cycle nallph
             else
! condition is for a specific phase#compset, N(phase#compset,comp)=A
                if(phr(jj)%iph.ne.sph .or. phr(jj)%ics.ne.scs) cycle nallph
                write(*,*)'N(phase#set,component) not implemented'
                gx%bmperr=4207; goto 1000
             endif
! moles formulat unit of phase set above
             pham=pmi%curd%amfu
!             write(*,*)'MM pham: ',phr(jj)%iph,pham
! if phase is not fixed there is a column in xcol for variable amount
! This has to be done before loop of elements
             if(pmi%phasestatus.ne.PHFIXED) notf=notf+1
             ncol=1
!             write(*,*)'Loop for elements: ', jj,phr(jj)%iph,phr(jj)%ics,ncol
             nallel: do ie=1,meqrec%nrel
! if sel=/=0 then skip all components except sel
                if(sel.gt.0 .and. ie.ne.sel) cycle nallel
! multiply terms with the inverse phase matrix
! This is called for each condition, maybe try to save values ...
                if(nosave) then
                   call calc_dgdyterms1(meqrec%nrel,ie,meqrec%tpindep,&
                        mamu,mag,mat,map,pmi,ceq%cmuval,meqrec%noofits)
                else
! this routine which should work in parallel ...
                   call calc_dgdyterms1X(meqrec%nrel,ie,meqrec%tpindep,&
                        mamu,mag,mat,map,pmi,meqrec%noofits)
                endif
                if(gx%bmperr.ne.0) goto 1000
! the call above calculates (A is "ie", z_ij is the inverted phase matrix): 
! mamu_A(B=1..nrel) = \sum_i \sum_j dM^a_A/dy_i dM^a_B/dy_j z^a_ij
! mag_A             = \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j
! mat_A             = \sum_i \sum_j d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j
! map_A             = \sum_i \sum_j d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j
! calculate a term for each column to be multiplied with chemical potential
! if the potential is fixed add the term to the rhs
!                 goto 8000  .... skipping nloop1 with je fails ....
!??????????????????? is this loop needed ?????????????????????? YES !!!
                ncol=1
                nloop1: do je=1,meqrec%nrel
                   do ke=1,meqrec%nfixmu
                      if(meqrec%mufixel(ke).eq.je) then
! components with fix chemical potential added to rhs, do not increment ncol!!!
                         xcol(nz2)=xcol(nz2)+pham*mamu(je)*meqrec%mufixval(ke)
!                         write(*,102)'fix mu N: ',sel,je,pham,&
!                              meqrec%mufixval(ke),mamu(je),&
!                              pham*mamu(je)*meqrec%mufixval(ke),xcol(nz2)
                         cycle nloop1
                      endif
                   enddo
! mamu(B) = \sum_i \sum_j \sum_A dM^a_B/dy_i dM^a_A z^a_ij
                   xcol(ncol)=xcol(ncol)-pham*mamu(je)
                   ncol=ncol+1
                enddo nloop1
!                goto 9000
!9000            continue
! If T or P are variable
                if(tcol.gt.0) then
                   xxx=xcol(tcol)
                   xcol(tcol)=xcol(tcol)+pham*mat
!                   write(*,363)'d2G/dTdy 2: ',nrow+1,ie,tcol,&
!                        xxx,xcol(tcol),pham,mat
                endif
! condition on N and variable P
                if(pcol.gt.0) then
                   xxx=xcol(pcol)
                   xcol(pcol)=xcol(pcol)+pham*map
!                   write(*,363)'MM d2G/dPdyi: ',nrow+1,ie,pcol,&
!                        xxx,xcol(pcol),pham,map
                endif
! last columns on lhs are amounts of element ie for all stable non-fix phases
! dncol should indicate last column with potential, can be different for
! derivative, notf is set above
                if(pmi%phasestatus.ne.PHFIXED) then
! notf indicates the column for amount of a component in stable nonfixed phase
! sum of moles in phase will be multiplied with delta-phase_amount
                   if(sel.gt.0 .and. sel.eq.ie) then
                      xcol(dncol+notf)=pmi%xmol(ie)
                   else
                      xcol(dncol+notf)=xcol(dncol+notf)+pmi%xmol(ie)
                   endif
                endif
! Maybe this should be included also for fixed phases ....?? YES
! right hand side (rhs) contribution is
! - NP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij
                xxx=xcol(nz2)
                xcol(nz2)=xcol(nz2)-pham*mag
             enddo nallel
! this is to used on the RHS for compare with prescribed value
             if(sel.gt.0) then
                totam=totam+pham*pmi%xmol(sel)
             else
                totam=totam+pham*pmi%sumxmol
             endif
! tafidbug
!             write(*,665)xxx,pham,mag,cvalue,totam,&
!                  xxx-pham*mag+cvalue-totam
665          format('RHS: ',6(1pe12.4))
          enddo nallph
!
! 160818: adding code to have several terms ... same as for x below
          nmany: if(mterms.lt.nterms) then
! this branch if 2 or more terms
             if(mterms.eq.1) then
! allocate arry to save intermediate results
! -Wuninitialized gave a warning: qmat.dim[0].ubound may be uninitilzed
! when used a few lines below but adding this removed this ... 
                if(.not.allocated(qmat)) then
                   allocate(qmat(nz2),stat=errall)
                   if(errall.ne.0) then
                      write(*,*)'MM Allocation error 29: ',errall
                      gx%bmperr=4370; goto 1000
                   endif
                endif
                qmat=zero
                evalue=zero
             endif
! save xcol and then go back and calculate next term
! maybe ccf should be included ??? YES!!! must correct also xterms!!!
             do ncol=1,nz2
                qmat(ncol)=qmat(ncol)+ccf(mterms)*xcol(ncol)
             enddo
             evalue=evalue+ccf(mterms)*totam
!             write(*,664)'MM nsel1:',moffs,sel,sph,totam,ccf(mterms),&
!                  cvalue,evalue
!             write(*,666)'MM evalue1: ',mterms,evalue,ccf(mterms),totam
!             write(*,666)'MM q:',mterms,evalue,(qmat(ncol),ncol=1,nz2)
666          format(a,i2,6(1pe12.4))
! prepare for next term by incrementing mterms and moffs
             mterms=mterms+1
             moffs=moffs+4
             deallocate(xcol)
!             deallocate(zcol)
             goto 1107
          elseif(nterms.gt.1) then
! for last term when more than 1
             nrow=nrow+1
             if(nrow.gt.nz1) then
                write(*,*)'MM too many equations 11A0',nrow
                gx%bmperr=4209; goto 1000
             endif
             do ncol=1,nz2
                smat(nrow,ncol)=qmat(ncol)+ccf(mterms)*xcol(ncol)
             enddo
             evalue=evalue+ccf(mterms)*totam
             smat(nrow,nz2)=smat(nrow,nz2)-cvalue+evalue
!             write(*,664)'MM nsel2:',moffs,sel,sph,totam,ccf(mterms),&
!                  cvalue,evalue
664          format(a,3i3,6(1pe12.4))
!             write(*,666)'MM evalue: ',mterms,evalue,ccf(mterms),totam
!             write(*,666)'MM s:',mterms,evalue,(smat(nrow,ncol),ncol=1,nz2)
! 160818: end code added for N(A)-N(B)
          else
! only one terms (original code unchanged)
! in xcol are values summed over all phases and components
! then copy summed columns to row nrow in matrix smat
             nrow=nrow+1
             if(nrow.gt.nz1) then
                write(*,*)'MM too many equations 11A',nrow
                gx%bmperr=4212; goto 1000
             endif
             do ncol=1,nz2
                smat(nrow,ncol)=xcol(ncol)
             enddo
! add N^prescribed - N^current to rhs (right hand side)
             xxx=smat(nrow,nz2)
! convergence problems using condition fix phase with amount >0, change sign ...
             smat(nrow,nz2)=smat(nrow,nz2)-cvalue+totam
             evalue=totam
          endif nmany
! tafidbug
!          smat(nrow,nz2)=smat(nrow,nz2)+cvalue-totam
!          write(*,355)'MM N: ',cvalue,totam,(smat(nrow,jj),jj=1,nz2)
355          format(a,6(1pe12.4))
!          write(*,363)'RHSN: ',nrow,nz2,0,smat(nrow,nz2),xxx,cvalue,totam,&
!               cvalue-totam
          deallocate(xcol)
! relative check for convergence if cvalue>1.0
!          if(abs(totam-cvalue).gt.ceq%xconv*max(1.0d0,abs(cvalue))) then
          if(abs(evalue-cvalue).gt.ceq%xconv*max(1.0d0,abs(cvalue))) then
             if(converged.lt.5) then
                converged=5
                cerr%mconverged=converged
                if(cerr%nvs.lt.10) then
                   cerr%nvs=cerr%nvs+1
                   cerr%typ(cerr%nvs)=5
                   cerr%val(cerr%nvs)=cvalue
                   cerr%dif(cerr%nvs)=evalue-cvalue
                endif
!                write(*,*)'5: converged=5',cerr%nvs
             endif
!          endif
             if(vbug) then
                if(sel.eq.0) then
                   write(*,266)'Unconverged condition N or N(A): ',sel,&
                        cvalue,evalue,evalue-cvalue
                else
                   write(*,266)'Unconverged condition N or N(A): ',sel,&
                        cvalue,evalue,evalue-cvalue
                endif
             endif
          endif
!----------------------------------------------------------
       elseif(stvnorm.gt.1) then
! only normallizing of N with respect to amount of moles (M) is allowed
          write(*,*)'N can only be normalled with M',stvix,stvnorm,cmix(2)
          gx%bmperr=4208; goto 1000
       else
!------------------------------------------------------------
! condition is x(A)=fix or x(phase,A)=fix or several terms for x(...)
! return here if several terms, value of xxmm ???
          moffs=0
1120      continue
! x(A)=fix and x(phase#set,A)=fix conditions. x(A)=N(A)/N; x(ph,A)=N(ph,A)/N(ph)
! above N=fix and N(A)=fix are treated as they have a "simple" summation, 
! We must sum over all phases and constituents for the normallizing factor
! definition: X(A)=N(A)/N; 
! derivative: dX(A)=dN(A)/N - N(A)/N**2 *dN
! sum dN(A) and dN at the same time and multiply the sums with 1/N 
! and -N(A)/N**2 in the end.
          if(cmix(3+moffs).eq.0) then
             write(*,*)'Condition NM=fix is illegal'
             gx%bmperr=4208; goto 1000
          elseif(cmix(4+moffs).eq.0) then
! condition is x(A)=fix
             sel=cmix(3+moffs); sph=0
          else
! condition is x(phase#set,A)=fix
!             write(*,33)cmix
33           format('Condition x(phase#set,A)=fix?',10i4)
             sel=cmix(5+moffs); sph=cmix(3+moffs); scs=cmix(4+moffs)
          endif
          if(.not.allocated(xxmm)) then
! this call returns the current fractions and total amounts.  We need
! to do it only once inside this subroutine. xxmm are deallocated at exit
             allocate(xxmm(meqrec%nrel),stat=errall)
             allocate(wwnn(meqrec%nrel),stat=errall)
             calcmolmass=.FALSE.
             if(errall.ne.0) then
                write(*,*)'MM Allocation error 30: ',errall
                gx%bmperr=4370; goto 1000
             endif
          endif
          if(.not.calcmolmass) then
             call calc_molmass(xxmm,wwnn,totalmol,totalmass,ceq)
             if(gx%bmperr.ne.0) goto 1000
             calcmolmass=.TRUE.
          endif
! two summations, zcol sums the term dN(A); xcol sums dN (as above)
          allocate(xcol(nz2),stat=errall)
          allocate(zcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 31: ',errall
             gx%bmperr=4370; goto 1000
          endif
          xcol=zero
          zcol=zero
          totam=zero
          zval=zero
          xval=zero
! LOOP FOR ALL PHASES (why not all stable??)
! dncol+notf indicate column for the amount of phases with variable amount
          notf=0
! sum over all phases to handle conditions like x(phase#set,A)=fix
! as the phase#set may not be stable
          bbug=zero
          xallph: do jj=1,meqrec%nphase
             pmi=>phr(jj)
             if(sph.eq.0) then
! skip this phase if not stable and condition not on a specific phase (sph)
! WOW COMPLICATION, I have another test for stability ... suck
                if(phr(jj)%stable.eq.0) cycle xallph
                pham=pmi%curd%amfu
             else
! condition on specific phase, skip this phase if not the right one
                if(phr(jj)%iph.ne.sph .or. &
                     phr(jj)%ics.ne.scs) cycle xallph
! note this destroys calculated values from calc_molmass above ...
                call calc_phase_molmass(sph,scs,xxmm,wwnn,&
                     totalmol,totalmass,amount,ceq)
                calcmolmass=.FALSE.
                pham=one
                totalmol=one
!                write(*,355)'MM cpm: ',totalmol,amount,pham,xxmm
! totalmol depend on amout of phase stable, irrelevant here
                if(gx%bmperr.ne.0) goto 1000
             endif
! notf indicates the column for the variable amount of the phase
             if(pmi%phasestatus.ne.PHFIXED) notf=notf+1
             xallel: do ie=1,meqrec%nrel
! we cannot skip summation over all element as that is needed for normallizing
! calculate a term for each column to be multiplied with chemical potential
! we must sum xcol for all elemenets and add to zcol for element sel
! if sel=/=0 then we sum also zcol(sel) for all phases
                if(nosave) then
                   call calc_dgdyterms1(meqrec%nrel,ie,meqrec%tpindep,&
                        mamu,mag,mat,map,pmi,ceq%cmuval,meqrec%noofits)
                else
                   call calc_dgdyterms1X(meqrec%nrel,ie,meqrec%tpindep,&
                        mamu,mag,mat,map,pmi,meqrec%noofits)
                endif
                if(gx%bmperr.ne.0) goto 1000
!                write(*,355)'MM dgdy: ',mamu
                ncol=1
                xloop2: do je=1,meqrec%nrel
!---------------------------------------------------------------------
! BIG TROUBLE HERE FOR FIXED CHEMICAL POTENTIAL !!!!! FIXED NOW ... NO!! ??
! but still problems combining with other conditions on H etc ...
! it works when we have N(A)=fix (code above) but not with x(A)=fix
! Calculate one column for each component to be multiplied with chem.pot.
! components with fix chemical potential added to rhs, do not increment ncol!!!
                   do ke=1,meqrec%nfixmu
! check for elements with fixed chemical potentials, they go to RHS
                      if(meqrec%mufixel(ke).eq.je) then
! the sign here should be opposite from xcol(ncol)= below
!                         write(*,*)'In xloop2: ',ie,ke,je,sel,nrow
!                         xcol(nz2)=xcol(nz2)-&
                         xcol(nz2)=xcol(nz2)+&
                              pham*mamu(je)*meqrec%mufixval(ke)
!                         bbug=bbug-pham*mamu(je)
!                         write(*,102)'fix mu xall: ',sel,je,pham,&
!                              meqrec%mufixval(ke),mamu(je),&
!                              pham*mamu(je)*meqrec%mufixval(ke),xcol(nz2)
! zcol needed because we have a normallized property (mole fraction)
! NOTE it should be ie here and NOT je ??? and opposite sign from xcol(nz2)
                         if(ie.eq.sel) then
                            zcol(nz2)=zcol(nz2)+&
                                 pham*mamu(je)*meqrec%mufixval(ke)
!                            write(*,102)'fix mu xsel: ',ie,je,pham,&
!                                 meqrec%mufixval(ke),mamu(je),&
!                                 pham*mamu(je)*meqrec%mufixval(ke),zcol(nz2)
!                            abug=-pham*mamu(je)
                         endif
                         cycle xloop2
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
                enddo xloop2
!-----------------------------------------------------------------------
! If T or P are variable, mat is \sum_i d2G/dy_idT, map is \sum_i d2G/dy_idP
                if(tcol.gt.0) then
                   xcol(tcol)=xcol(tcol)+pham*mat
                   if(sel.eq.ie) then
                      zcol(tcol)=zcol(tcol)+pham*mat
                   endif
!                   write(*,363)'d2G/dTdy 3: ',nrow+1,ie,tcol,&
!                        xxx,xcol(tcol),pham,mat
363                format(a,3i3,6(1pe12.4))
                endif
                if(pcol.gt.0) then
                   xcol(pcol)=xcol(pcol)+pham*map
                   if(sel.eq.ie) then
                      zcol(pcol)=zcol(pcol)+pham*map
                   endif
                endif
! columns for phase amounts
                if(pmi%phasestatus.ne.PHFIXED) then
!                   write(*,*)'MM 363A: ',dncol,notf,ie,sel
                   if(sph.eq.0) then
                      xcol(dncol+notf)=xcol(dncol+notf)+pmi%xmol(ie)
!                   write(*,*)'MM 363B: ',dncol,notf,ie,xcol(dncol+notf)
                      if(ie.eq.sel) then
                         zcol(dncol+notf)=zcol(dncol+notf)+pmi%xmol(ie)
                      endif
                   endif
                endif
! right hand side (rhs) contribution is (normallized below)
! - NP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij 
                xcol(nz2)=xcol(nz2)-pham*mag
                if(sel.eq.ie) then
                   zcol(nz2)=zcol(nz2)-pham*mag
                endif
             enddo xallel
! totam and zval not used !!??
             totam=totam+pham*pmi%sumxmol
! UNFINISHED: if sph nonzero next line must be changed to be for sph
             zval=zval+pham*pmi%xmol(sel)
!             sel=cmix(5); sph=cmix(3); scs=cmix(4)
!             write(*,*)'MM x(p,c): ',sph,scs,sel,zval
          enddo xallph
!-------------- new code begin
! can handle the case of several terms like x(liquid,S)-x(pyrrh,S)=0
!                                       x(Mg)-2*x(Si)=0
          xterms: if(mterms.lt.nterms) then
! this branch if 2 or more terms
             if(mterms.eq.1) then
! allocate array for saving intermediate results
                if(.not.allocated(qmat)) then
                   allocate(qmat(nz2),stat=errall)
                   if(errall.ne.0) then
                      write(*,*)'MM Allocation error 32: ',errall
                      gx%bmperr=4370; goto 1000
                   endif
                endif
                qmat=zero
                evalue=zero
             endif
! save zcol and xcol then go back and calculate next term
! corrected by adding ccf factor!! (not needed for x(liq,a)-x(sol,a)=0 ....
             do ncol=1,nz2
                qmat(ncol)=qmat(ncol)+ccf(mterms)*&
                     (zcol(ncol)-xcol(ncol)*xxmm(sel))/totalmol
             enddo
             evalue=evalue+ccf(mterms)*xxmm(sel)
! prepare for next term by incrementing mterms and moffs
             mterms=mterms+1
             moffs=moffs+4
!             write(*,1117)'MM 2nd indices: ',moffs,(cmix(jj+moffs),jj=3,6)
!1117         format(a,i3,2x,4i3)
!             write(*,1118)'MM xxmm:',mterms,sel,xxmm(sel)
             deallocate(xcol)
             deallocate(zcol)
             goto 1120
          elseif(nterms.gt.1) then
! for last term of expression
             nrow=nrow+1
             if(nrow.gt.nz1) then
                write(*,*)'MM too many equations 11B: ',nrow,nz1,meqrec%nfixph
                gx%bmperr=4209; goto 1000
             endif
! insert results in smat
!             write(*,1118)'MM endofexp:',mterms,sel,evalue,xxmm(sel)
1118         format(a,2i3,6(1pe12.4))
             do ncol=1,nz2
                smat(nrow,ncol)=qmat(ncol)+&
                     ccf(mterms)*(zcol(ncol)-xcol(ncol)*xxmm(sel))/totalmol
             enddo
             evalue=evalue+ccf(mterms)*xxmm(sel)
! add x^prescribed - x^current to rhs (right hand side)
             smat(nrow,nz2)=smat(nrow,nz2)-cvalue+evalue
!------------------new code end
          else
! use this else branch when nterms=1, just a single x(a)=value
             nrow=nrow+1
!             if(bbug.ne.zero) then
! looking for bug with activity conditions
!                write(*,16)'abug: ',sel,abug,bbug,xxmm(sel),&
!                     abug-bbug*xxmm(sel),meqrec%mufixval(1),&
!                     (abug-bbug*xxmm(sel))*meqrec%mufixval(1)
!16              format(a,i3,6(1pe12.4))
!             else
!                write(*,16)'nomy : ',sel,zcol(1),xcol(1),&
!                     xxmm(sel),zcol(1)-xcol(1)*xxmm(sel)
!             endif
             if(nrow.gt.nz1) then
                write(*,*)'MM too many equations 11B: ',nrow,nz1,meqrec%nfixph
                gx%bmperr=4209; goto 1000
             endif
! in xcol is dN and in zcol dN(A) summed over all phases and components
! calculate the normallized values now
! xmat=1/N*(dN(A) - (N(A)/N)*dN)
! sum zcol and xcol to nrow in smat multiplying xcol with current amount
! and normallizing with total amount, including the RHS (column nz2)
             do ncol=1,nz2
                smat(nrow,ncol)=(zcol(ncol)-xcol(ncol)*xxmm(sel))/totalmol
             enddo
! subract x^prescribed - x^current to rhs (right hand side)
             smat(nrow,nz2)=smat(nrow,nz2)-cvalue+xxmm(sel)
             evalue=xxmm(sel)
          endif xterms
          deallocate(xcol)
          deallocate(zcol)
! phase composition problem
!          write(*,355)'MM X: ',cvalue,xxmm(sel),totalmol,pham,&
!               (smat(nrow,jj),jj=1,nz2)
! check on convergence
!          if(abs(xxmm(sel)-cvalue).gt.ceq%xconv) then
          if(abs(evalue-cvalue).gt.ceq%xconv) then
             if(converged.lt.5) then
                converged=5
!                write(*,*)'6: converged=5',cerr%nvs
                cerr%mconverged=converged
                if(cerr%nvs.lt.10) then
                   cerr%nvs=cerr%nvs+1
                   cerr%typ(cerr%nvs)=5
                   cerr%val(cerr%nvs)=xxmm(sel)
                   cerr%dif(cerr%nvs)=xxmm(sel)-cvalue
                endif
             endif
!             write(*,266)'Unconverged condition x(A): ',sel,cvalue,evalue
!             if(vbug) write(*,266)'Unconverged condition x(A): ',sel,&
!                  cvalue,evalue
          endif
       endif
! finished conditions on N and X with indices
       if(allocated(xxmm)) then
          deallocate(xxmm)
          deallocate(wwnn)
       endif
!
!------------------------------------------------------------------
  case(12) ! B or W
! Amount of component in mass, can have indices and normallization
! code copied from the case(11) for N and X and modified for mass
1200   continue
       if(stvnorm.eq.0) then
          if(cmix(3).eq.0) then
! condition is B=fix
             if(bwarning) then
                write(*,491)
491             format(' *** WARNING, using B=value as condition can disable',&
                     ' the gridminimizer'/&
                     ' and cause convergence problem. Use N=value instead.')
! Issue this message only once for each calculation
                bwarning=.FALSE.
             endif
!             write(*,*)'MM condition B=fix: ',stvnorm,cmix(3)
             sel=0; sph=0
          elseif(cmix(4).eq.0) then
! condition is B(A)=fix
             sel=cmix(3); sph=0
          else
! condition is B(phase#set,A)=fix;  how to handle if phase#set not stable?
             write(*,*)'Condition B(phase#set,A)=fix not implemented'
             gx%bmperr=4208; goto 1000
             sel=cmix(5); sph=cmix(3); scs=cmix(4)
          endif
! Formulate equation for total amount B: each M_A multiplied with mass_A
! rhs:  B-B+\sum_alpha N^a + \sum_i \sum_j dM^a_A/dy_i z^a_ij dG/dy_j 
! \sum_B \sum_alpha N^a \sum_i \sum_j dM^_A/dy_i dM^a_B/dy_j*z^a_ij  *mu(B)
!        \sum_alpha N^a \sum_i d2M^a_A/dTdy_i z^a_ij d2G/dTdy_j      *deltaT
!        \sum_alpha N^a \sum_i d2M^a_A/dPdy_i z^a_ij d2G/dPdy_j      *deltaP
!        \sum_A M^a_A                                    *deltaN^a
          allocate(xcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 33: ',errall
             gx%bmperr=4370; goto 1000
          endif
          xcol=zero
          totam=zero
!          write(*,222)'MM xcol 1',totam,xcol
222       format(a,10(1pe11.3))
! notf keeps track on entered non-fixed phases with variable amount
          notf=0
! not used          zval=zero
          ballph: do jph=1,meqrec%nstph
! sum over all stable phases
             jj=meqrec%stphl(jph)
             pmi=>phr(jj)
! if phase is not fixed there is a column in xcol for variable amount
             if(pmi%phasestatus.ne.PHFIXED) notf=notf+1
! amount of phase, amfu is moles formula units, abnorm(2) is mass per form.unit
             pham=pmi%curd%amfu
             ballel: do ie=1,meqrec%nrel
! if sel=/=0 then skip all components except sel
                if(sel.gt.0 .and. ie.ne.sel) cycle
! multiply terms with the inverse phase matrix
                if(nosave) then
                   call calc_dgdyterms1(meqrec%nrel,ie,meqrec%tpindep,&
                        mamu,mag,mat,map,pmi,ceq%cmuval,meqrec%noofits)
                else
                   call calc_dgdyterms1X(meqrec%nrel,ie,meqrec%tpindep,&
                        mamu,mag,mat,map,pmi,meqrec%noofits)
                endif
                if(gx%bmperr.ne.0) goto 1000
!                write(*,*)'Calculated dgdyterms 3: ',mat
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
!                         write(*,98)'fix mu b:',sel,je,&
!                              pham*mamu(je),meqrec%mufixval(ke),mass_of(ie,ceq)
                         xcol(nz2)=xcol(nz2)+&
                              pham*mamu(je)*meqrec%mufixval(ke)*mass_of(ie,ceq)
                         cycle bloop1
                      endif
                   enddo
! mamu(B) = \sum_i \sum_j \sum_A dM^a_B/dy_i dM^a_A z^a_ij mass_A mass_B ???
                   xcol(ncol)=xcol(ncol)-pham*mamu(je)*mass_of(ie,ceq)
                   ncol=ncol+1
                enddo bloop1
! If T or P are variable
                if(tcol.gt.0) then
!                   xxx=xcol(tcol)
                   xcol(tcol)=xcol(tcol)+pham*mat*mass_of(ie,ceq)
!                   write(*,363)'d2G/dTdy 4: ',nrow-1,ie,tcol,&
!                        xxx,xcol(tcol),pham,mat
                endif
                if(pcol.gt.0) then
!                   xxx=xcol(pcol)
                   xcol(pcol)=xcol(pcol)+pham*map*mass_of(ie,ceq)
!                   write(*,363)'d2G/dPdy: ',nrow-1,ie,pcol,&
!                        xxx,xcol(pcol),pham,mat
                endif
! last columns are amounts of element ie for all stable non-fix phases
! for all stable (non fixed) phases we have the mass multiplied with deltaaleph
                if(pmi%phasestatus.ne.PHFIXED) then
! ??                    zval=zval+pmi%xmol(ie)*mass_of(ie,ceq)
                   if(sel.gt.0 .and. sel.eq.ie) then
                      xcol(dncol+notf)=&
                           pmi%xmol(ie)*mass_of(ie,ceq)
!                     write(*,363)'xcola: ',ncol,ie,0,xcol(ncol),mass_of(ie,ceq)
                   else
                      xcol(dncol+notf)=xcol(dncol+notf)+&
                           pham*pmi%xmol(ie)*mass_of(ie,ceq)
                   endif
                endif
! right hand side (rhs) contribution is
! - BP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij
                xcol(nz2)=xcol(nz2)-pham*mag*mass_of(ie,ceq)
!                write(*,222)'MM xcol 2',totam,xcol
             enddo ballel
! sum of mass in phase will be multiplied with delta-phase_amount
!             write(*,202)'sumxmol mm:  ',sel,pham,pmi%sumxmol,pmi%sumwmol
             if(sel.gt.0) then
                totam=totam+pham*pmi%xmol(sel)*mass_of(sel,ceq)
             else
                totam=totam+pham*pmi%sumwmol
             endif
          enddo ballph
!          write(*,222)'MM xcol 3',totam,xcol
!......debug
          if(.not.allocated(xxmm)) then
! this call returns the current fractions and total amounts.  We need
! to do it only once inside this subroutine. xxmm are deallocated at exit
             allocate(xxmm(meqrec%nrel),stat=errall)
             allocate(wwnn(meqrec%nrel),stat=errall)
             if(errall.ne.0) then
                write(*,*)'MM Allocation error 34: ',errall
                gx%bmperr=4370; goto 1000
             endif
             call calc_molmass(xxmm,wwnn,totalmol,totalmass,ceq)
             if(gx%bmperr.ne.0) goto 1000
          endif
!          if(sel.eq.0) write(*,*)'totalmass: ',totalmass,totam
!
! in xcol are values summed over all phases and components
! copy summed columns to smat nrow
          nrow=nrow+1
          if(nrow.gt.nz1) then
             write(*,*)'MM too many equations 12A',nrow
             gx%bmperr=4209; goto 1000
          endif
          do ncol=1,nz2
             smat(nrow,ncol)=xcol(ncol)
          enddo
!          write(*,97)'Totalmass B: ',sel,totam,cvalue,totalmass,wwnn(sel)
97        format(a,i4,6(1pe12.4))
! add B^prescribed - B^current to rhs (right hand side)
          xxx=smat(nrow,nz2)
          smat(nrow,nz2)=smat(nrow,nz2)-cvalue+totam
!          write(*,363)'RHSB: ',nrow,nz2,0,smat(nrow,nz2),xxx,cvalue,totam,&
!               cvalue-totam
          deallocate(xcol)
! check convergence
          if(abs(totam-cvalue).gt.ceq%xconv) then
!            write(*,266)'Unconverged condition B(A): ',sel,cvalue,zval
             if(converged.lt.5) then
                converged=5
!                write(*,*)'7: converged=5',cerr%nvs
                cerr%mconverged=converged
                if(cerr%nvs.lt.10) then
                   cerr%nvs=cerr%nvs+1
                   cerr%typ(cerr%nvs)=5
                   cerr%val(cerr%nvs)=cvalue
                   cerr%dif(cerr%nvs)=totam-cvalue
                endif
             endif
          endif
!          write(*,222)'MM xcol 3',totam,xcol
          if(vbug) then
             if(sel.eq.0) then
                write(*,363)'Condition B=fix',0,0,0,cvalue,totam
             else
                write(*,363)'Condition B(a)=fix',sel,0,0,cvalue,totam
             endif
          endif
!          write(*,223)'MM smat 1',nrow,(smat(nrow,ncol),ncol=1,nz2)
223       format(a,i2,10(1pe11.3))
       elseif(stvnorm.ne.2) then
! only normallizing of B with respect to mass (W) is allowed
          write(*,*)'Allowed normallizing with W only',stvix,stvnorm,cmix(2)
          gx%bmperr=4208; goto 1000
       else
!-------------------------------
! Conditions like w(A)=fix, w(phase#set,A)=fix
! B=fix and B(A)=fix treated above as they have a "simple" summation, 
! We must sum over all phases and constituents for the normallizing factor
! definition: W(A)=B(A)/B; 
! derivative: dW(A)=dB(A)/B - B(A)/N**2 *dB
! sum dB(A) and dB at the same time and multiply the sums with 1/B
! and -B(A)/B**2 in the end.
          if(cmix(3).eq.0) then
             write(*,*)'Condition BW=fix is illegal'
             gx%bmperr=4208; goto 1000
          elseif(cmix(4).eq.0) then
! condition is x(A)=fix
             sel=cmix(3); sph=0
          else
             sel=cmix(5); sph=cmix(3); scs=cmix(4)
          endif
          if(.not.allocated(xxmm)) then
! this call returns the current fractions and total amounts.  We need
! to do it only once inside this subroutine. xxmm are deallocated at exit
             allocate(xxmm(meqrec%nrel),stat=errall)
             allocate(wwnn(meqrec%nrel),stat=errall)
             if(errall.ne.0) then
                write(*,*)'MM Allocation error 35: ',errall
                gx%bmperr=4370; goto 1000
             endif
             calcmolmass=.FALSE.
          endif
          if(.not.calcmolmass) then
             call calc_molmass(xxmm,wwnn,totalmol,totalmass,ceq)
             if(gx%bmperr.ne.0) goto 1000
             calcmolmass=.TRUE.
          endif
!          write(*,267)'wwnn: ',(wwnn(ncol),ncol=1,noel())
!          write(*,267)'xxmm: ',(xxmm(ncol),ncol=1,noel())
! two summations, zcol sums the term dN(A); xcol sums dN (as above)
          allocate(xcol(nz2),stat=errall)
          allocate(zcol(nz2),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 36: ',errall
             gx%bmperr=4370; goto 1000
          endif
          xcol=zero
          zcol=zero
          totam=zero
          zval=zero
          xval=zero
          notf=0
!          wallph: do jph=1,meqrec%nstph
!             jj=meqrec%stphl(jph)
! sum over all phases to handle conditions like x(phase#set,A)=fix
! as the phase#set may not be stable
          wallph: do jj=1,meqrec%nphase
             pmi=>phr(jj)
             if(sph.eq.0) then
! skip this phase if not stable and condition not on a specific phase
                if(phr(jj)%stable.eq.0) cycle wallph
                pham=pmi%curd%amfu
             elseif(sph.gt.0) then
! condition on a composition of a phase
                if(phr(jj)%iph.ne.sph .or. &
                     phr(jj)%ics.ne.scs) cycle wallph
! We need the phase comoposition
                call calc_phase_molmass(sph,scs,xxmm,wwnn,&
                     totalmol,totalmass,amount,ceq)
                pham=one
!                totalmol=one
                totalmass=one
             endif
!             pmi=>phr(jj)
! amount formula units of phase, set above
!             pham=pmi%curd%amfu
             if(pmi%phasestatus.ne.PHFIXED) notf=notf+1
             wallel: do ie=1,meqrec%nrel
! calculate a term for each column to be multiplied with chemical potential
! we must sum xcol for all elemenets and add to zcol for element sel
! if sel=/=0 then we sum also zcol(sel) for all phases
                if(nosave) then
                   call calc_dgdyterms1(meqrec%nrel,ie,meqrec%tpindep,&
                        mamu,mag,mat,map,pmi,ceq%cmuval,meqrec%noofits)
                else
                   call calc_dgdyterms1X(meqrec%nrel,ie,meqrec%tpindep,&
                        mamu,mag,mat,map,pmi,meqrec%noofits)
                endif
                if(gx%bmperr.ne.0) goto 1000
!                write(*,*)'Calculated dgdyterms 4: ',mat
                ncol=1
! BUG TROUBLE WITH MIXED FIX CHEMICAL POT AND MASS FRACTION CONDITION !!!
                wloop2: do je=1,meqrec%nrel
! Calculate one column for each component to be multiplied with chem.pot.
! components with fix chemical potential added to rhs, do not increment ncol!!!
! modified in accordance with condition on x
                   do ke=1,meqrec%nfixmu
                      if(meqrec%mufixel(ke).eq.je) then
!                         write(*,98)'fix mu w:',sel,je,&
!                              pham*mamu(je),meqrec%mufixval(ke),mass_of(ie,ceq)
98                       format(a,2i3,6f12.4)
                         xcol(nz2)=xcol(nz2)+&
                              pham*mamu(je)*meqrec%mufixval(ke)*mass_of(ie,ceq)
                         if(ie.eq.sel) then
!                            write(*,98)'fix mu u:',sel,ie,&
!                                 pham*mamu(je),meqrec%mufixval(ke),&
!                                 mass_of(ie,ceq)
! VERY STRANGE ... zcol and xcol have both the term added here but
! when calculating with mole frac and fix chem.pot they have different signs!!!
                            zcol(nz2)=zcol(nz2)+&
                                 pham*mamu(je)*meqrec%mufixval(ke)*&
                                 mass_of(ie,ceq)
                         endif
                         cycle wloop2
                      endif
                   enddo
! mamu(B) = \sum_i \sum_j dM^a_B/dy_i dM^a_A z^a_ij
                   xcol(ncol)=xcol(ncol)-pham*mamu(je)*mass_of(ie,ceq)
                   if(sel.eq.ie) then
                      zcol(ncol)=zcol(ncol)-pham*mamu(je)*mass_of(ie,ceq)
                   endif
! problem that this reurn whatever for 2nd and higher equilibria
!                   write(*,*)'mass of: ',ie,mass_of(ie,ceq)
                   ncol=ncol+1
                enddo wloop2
! If T or P are variable
                if(tcol.gt.0) then
                   xcol(tcol)=xcol(tcol)+pham*mat*mass_of(ie,ceq)
                   if(sel.eq.ie) then
                      zcol(tcol)=zcol(tcol)+pham*mat*mass_of(ie,ceq)
                   endif
!                   write(*,363)'d2G/dTdy 5: ',nrow-1,ie,tcol,&
!                        xxx,xcol(tcol),pham,mat
                endif
                if(pcol.gt.0) then
                   xcol(pcol)=xcol(pcol)+pham*map*mass_of(ie,ceq)
                   if(sel.eq.ie) then
                      zcol(pcol)=zcol(pcol)+pham*map*mass_of(ie,ceq)
                   endif
                endif
! last columns are amounts of element ie for all stable non-fix phase,
                if(pmi%phasestatus.ne.PHFIXED) then
                   if(sph.eq.0) then
! all phases with variable amount, sum over all components
                      xcol(dncol+notf)=xcol(dncol+notf)+&
                           pmi%xmol(ie)*mass_of(ie,ceq)
                      if(ie.eq.sel) then
                         zcol(dncol+notf)=zcol(dncol+notf)+&
                              pmi%xmol(ie)*mass_of(ie,ceq)
                      endif
!                   else
! no coefficint for phase amount if phase specific composition!!
                   endif
                endif
! right hand side (rhs) contribution is
! - NP(phase)*\sum_i \sum_j dM(ie)/dy_i * dG/dy_j * z_ij * mass_ie
                xcol(nz2)=xcol(nz2)-pham*mag*mass_of(ie,ceq)
                if(sel.eq.ie) then
                   zcol(nz2)=zcol(nz2)-pham*mag*mass_of(ie,ceq)
                endif
             enddo wallel
! totam never used ???
             if(sel.gt.0) then
                totam=totam+pham*pmi%xmol(sel)*mass_of(sel,ceq)
             else
                totam=totam+pham*pmi%sumwmol
             endif
! UNFINISHED: if sph=/=0 next line must be changed
!             zval=zval+pham*pmi%xmol(sel)*mass_of(sel,ceq)
          enddo wallph
! in xcol is dB and in zcol dB(A) summed over all phases and components
! calculate the normallized values now
! xmat=dB(A)/B - B(A)*dB/B**2
          nrow=nrow+1
          if(nrow.gt.nz1) then
             write(*,*)'MM too many equations 12B',nrow,nz1
             gx%bmperr=4209; goto 1000
          endif
!          write(*,97)'Totalmass W: ',sel,wwnn(sel),cvalue,totalmass,totam
! copy to smat row nrow.  totalmass=1 if phase specific composition
          do ncol=1,nz2
             smat(nrow,ncol)=(zcol(ncol)-xcol(ncol)*wwnn(sel))/totalmass
          enddo
! add W^prescribed - W^current to rhs (right hand side)
          smat(nrow,nz2)=smat(nrow,nz2)-cvalue+wwnn(sel)
          deallocate(xcol)
          deallocate(zcol)
! check on convergence
!          write(*,266)'massbalance condition w(A): ',sel,cvalue,wwnn(sel)
          if(abs(wwnn(sel)-cvalue).gt.ceq%xconv) then
             if(converged.lt.5) then
                converged=5
!                write(*,*)'8: converged=5',cerr%nvs
                cerr%mconverged=converged
                if(cerr%nvs.lt.10) then
                   cerr%nvs=cerr%nvs+1
                   cerr%typ(cerr%nvs)=5
                   cerr%val(cerr%nvs)=wwnn(sel)
                   cerr%dif(cerr%nvs)=wwnn(sel)-cvalue
                endif
!                write(*,*)'8B: converged=5',cerr%nvs
             endif
!             write(*,266)'Unconverged condition w(A): ',sel,cvalue,wwnn(sel)
266          format(a,i3,3(1pe14.6))
!             write(*,267)'wwnn: ',(wwnn(ncol),ncol=1,noel())
!             write(*,267)'xxmm: ',(xxmm(ncol),ncol=1,noel())
!267          format(a,8F9.5)
          endif
!          if(sph.eq.0) then
!             write(*,363)'Condition w(A)=fix',sel,0,0,cvalue,wwnn(sel)
!          else
! this is not implemented yet
!             write(*,363)'Condition w(phase#set,A)=fix',sph,sel,0,cvalue,zval
!          endif
       endif
! finished conditions on B and W with indices
       if(allocated(xxmm)) then
          deallocate(xxmm)
          deallocate(wwnn)
       endif
!
!------------------------------------------------------------------
    case(13) ! Y ycond
! Constituent fraction: phase#set, (subl.,) constituent index (over all subl)
! NOTE differences also interesting y(B2,A)-y(B2,A#2) is 2nd order transf
! nterms is number of terms, mterms=1 here
       moffs=3
!       write(*,*)'MM stvix, mterms & nterms: ',stvix,mterms,nterms,nz2
! xcol not needed as we have no sums over several phases
!       allocate(xcol(nz2))
!       xcol=zero
! we do not use calc_dgdyterms as we have a single constituent yindex
!             call calc_dgdyterms1X(meqrec%nrel,ie,meqrec%tpindep,&
!                  mamu,mag,mat,map,pmi,meqrec%noofits)
! mamu is an array, normally set to zero in calc_dgdyterms
! also mag, mat and map
       mamu=zero
       mag=zero
       mat=zero
       map=zero
       dvalue=zero
! this is executed for each iteration, this value must be set earlier
!       deltaTycond=2.5d1
       yterms: do mterms=1,nterms
! loop for all terms in constion, we may have y_i-y_j =fix
! cmix(3,4,5,6) are for first term, cmix(7,8,9,10) for second etc
! for each term ccf(1..5) gives the factor in front of y
! constituent is cmix(3) (sequental for all sublattices?)
!          write(*,*)'MM phase and compset   :',cmix(moffs),cmix(moffs+1)
!          write(*,'(a,i3,2(1pe12.4))')'MM constituent & value :',&
!               cmix(moffs+2),cvalue,ccf(mterms)
! cmix(moffs+4) NOT USED          
          sph=cmix(moffs); scs=cmix(moffs+1)
          yindex=cmix(moffs+2)
          findphase: do jj=1,meqrec%nphase
! phase with y condition may not be stable ... loop for all phases
!             write(*,*)'MM phase: ',meqrec%nphase,jj,phr(jj)%iph
             if(phr(jj)%iph.ne.sph .or. phr(jj)%ics.ne.scs) cycle findphase
             pmi=>phr(jj)
!             write(*,*)'MM found phase: ',jj,yindex,phr(jj)%curd%yfr(yindex)
! The equation is \Delta y_i = yknown (or \Delta y_i - \Delta y_j = dyknown)
! we should set up a row where index "i"  is known constituent
!  \sum_A \sum_k dM_A/dy_i e_ik \mu_A + \sum_k d2G/dy_i dT e_ik \Delta T = 
!                                       \sum_k dG/dy_i e_ik +ycurr - yknown
! where e_ij is the inverted phase matrix
! The values of the constituent fractions must be set before calculating e_ij
! this requires some new indicator in meq_onephase
! IF the condition is a difference y_i-y_j=a it will be assumed y_i is correct
! at the start of the calculation and we set y_j=y_i-a before each iteration
             yallel: do ie=1,meqrec%nrel
                cib=zero
!                write(*,333)'MM dy: ',(pmi%dxmol(ie,jy),jy=1,pmi%ncc)
333             format(a,10(1pe12.4))
                do jy=1,pmi%ncc
! \sum_A \sum_k  e_ik dM_A/dy_k
! suck the formula below does not work unless y_i correct, suck
                   cib=cib+pmi%invmat(jy,yindex)*pmi%dxmol(ie,jy)
!                   write(*,'(a,i3,3(1pe12.4))')'MM cib 1: ',jy,cib,&
!                        pmi%invmat(jy,yindex),pmi%dxmol(ie,jy)
                enddo
                mamu(ie)=mamu(ie)+ccf(mterms)*cib
!                write(*,*)'MM mamu: ',ie,mamu(ie),cib
             enddo yallel
             cib=zero
             do jy=1,pmi%ncc
! \sum_k e_ik dG/dy_k
                cib=cib+pmi%invmat(jy,yindex)*pmi%curd%dgval(1,jy,1)
!                write(*,'(a,i3,3(1pe12.4))')'MM cib 2: ',jy,cib,&
!                     pmi%invmat(jy,yindex),pmi%curd%dgval(1,jy,1)
             enddo
! WoW it works with correct signs!  Note: y_presc - y_calc!!!
             dvalue=-ccf(mterms)*pmi%curd%yfr(yindex)
             mag=mag+ccf(mterms)*(cib-pmi%curd%yfr(yindex))
!             write(*,373)'MM mag: ',mag,cib,&
!                  ccf(mterms),-pmi%curd%yfr(yindex),cvalue
             if(meqrec%tpindep(1)) then
! failed attempt to improve convergence
!                ycondTlimit=.true.
! add coefficient for Delta T
                cib=zero
                do jy=1,pmi%ncc
! + \sum_k e_ik d2G/dTdy_ik \Delta T 
                   cib=cib+pmi%invmat(jy,yindex)*pmi%curd%dgval(2,jy,1)
! OR: + \sum_k e_ik d2G/dTdy_i  \Delta T 
! I have not tested eithor of these
!                   cib=cib+pmi%invmat(jy,yindex)*pmi%curd%dgval(2,yindex,1)
                enddo
! When T is variable with y condition one must restrict change in T !!!
                mat=mat+ccf(mterms)*cib
             endif
!             write(*,'(a,i2,6(1pe12.4))')'MM mat: ',mterms,&
!                  ceq%tpval(1),dvalue+cvalue,ccf(mterms),mat,cib
             if(meqrec%tpindep(2)) then
! add coefficient for Delta P
                cib=zero
                do jy=1,pmi%ncc
! + \sum_k e_ik d2G/dPdy_i  \Delta P
! I have not tested this
                   cib=cib+pmi%invmat(jy,yindex)*pmi%curd%dgval(3,jy,1)
                enddo
                map=map+ccf(mterms)*cib
             endif
             exit findphase
          enddo findphase
! finished this term, any more?
          moffs=moffs+4
       enddo yterms
! add the prescribed value
       mag=mag+cvalue
! dvalue is the current value which should become cvalue at equilibrium
       dvalue=dvalue+cvalue
!       write(*,373)'MM mamu: ',mat,map,mag,mamu
!       write(*,*)'MM nrow mm: ',nrow,nz1,nz2
373    format(a,10(1pe12.4))
!-------------------
       nrow=nrow+1
       if(nrow.gt.nz1) then
          write(*,*)'MM Too many equations 13'
          gx%bmperr=4209; goto 1000
       endif
! now mamu(1..nrel) are the coefficients for \mu; mat&map is coeff for Delta T&P
! assuming no activity conditionw ...
       do jj=1,meqrec%nrel
          smat(nrow,jj)=mamu(jj)
       enddo
! after exiting loop jj=meqrec%nrel+1
       if(meqrec%tpindep(1)) then
! Failed attempt to improve convergence
!          smat(nrow,jj)=-5.0D0*mat
          smat(nrow,jj)=-mat
          jj=jj+1
       endif
       if(meqrec%tpindep(2)) then
          smat(nrow,jj)=-map
          jj=jj+1
       endif
! mag is right hand side including y-y
       smat(nrow,nz2)=mag
!    write(*,'(a,i2,6(1pe12.4))')'MM *** ycond:',nrow,(smat(nrow,jj),jj=1,nz2)
!       gx%bmperr=4207; goto 1000
! 
    end select
!
! loop if not the last condition
!    write(*,*)'Taking next condition',cmix(1)
    if(.not.associated(condition,lastcond)) goto 350
!=====================================================================
380 continue
! write whole smat
! used to find ycond ....
!    do jj=1,nz1
!       write(*,390)jj,(smat(jj,jy),jy=1,nz2)
!    enddo
390 format('#:',i2,6(1pe12.4),6(4x,1pe12.4))
1000 continue
! we must ?? deallocate all data in the savedrec
!    if(allocated(savedrec%save1)) then
!       jj=size(saved%save1)
!       deallocate(savedrec%save1)
!       write(*,*)'MM deallocated saved%save1',jj
!    endif
!    if(allocated(savedrec%save2)) deallocate(savedrec%save2)
!    if(allocated(savedrec%save3)) deallocate(savedrec%save3)
!    if(allocated(savedrec%save4)) deallocate(savedrec%save4)
!    if(allocated(savedrec%save5)) deallocate(savedrec%save5)
    return
  end subroutine setup_equilmatrix

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine meq_onephase
!\begin{verbatim}
  subroutine meq_onephase(meqrec,pmi,ceq)
! this subroutine calculates new constituent fractions for a phase iph+ics
! with given T, P and chemical potentials for the components 
! For ionic liquids the sites on the sublattices varies with composition
! THIS IS A FIRST VERSION WITHOUT ANY TRICKS FOR SPEED
! this will check if EEC set and modify G for solid phases with higher entropy
! pmi is pointer to a record in meq_phase, local to this thread
! than the liquid
    implicit none
    TYPE(meq_phase), pointer :: pmi
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(meq_setup) :: meqrec
!\end{verbatim}
    integer nrel,i2sly(2),info
    integer ik,iph,ics,jz,iz,jk,ierr,kk,kkk,ll,lokcs,ncc,loksp,ncl
!    integer nd1,nd2,neq,nochange,nsl,nspel,nv,ncon,icon,jxsym,kxsym
    integer nd1,nd2,neq,nochange,nsl,nspel,nv,ncon,icon,jxsym,errall
! needed for call to get_phase_data
    integer, dimension(maxsubl) ::  nkl
    integer, dimension(maxconst) :: knr
    double precision, dimension(5) :: qq
    double precision, dimension(maxsubl) :: sites
! needed for call to get_species_data
    integer, dimension(maxspel) :: ielno
    double precision, dimension(maxspel) :: stoi
! testing lapacl+blas inverting symmetric matrix
    double precision, allocatable, dimension(:) :: lapack
    double precision xxxx,yyyy
! minimal y, charge
    double precision, parameter :: ymin=1.0D-12,ymingas=1.0D-30,qeps=1.0D-30
! derivative of moles of component wrt y_ks
    double precision, dimension(maxel) :: addmol
! for mass balance and charge
    double precision, dimension(maxconst) :: yarr,dqsum
! phase matrix, its inverse is returned as part of pmi
    double precision, dimension(:,:), allocatable :: pmat
    double precision qsp,sumsit,ykvot,ysum,qsum,spmass,yva,fion
    double precision, dimension(:,:), allocatable :: sumion
    character name*24
!    logical nolapack
!    write(*,'(a,5i5)')'in meq_onephase: ',ceq%eqno,&
!         pmi%iph,pmi%ics,meqrec%noofits
! set eecextrapol to TRUE when entering, 
! set to FALSE inside check_eec if phase has higher entropy than liquid
! I am no longer sure this is needed??
!    eecextrapol=.TRUE.
! Maybe nolapack be removed??
!    nolapack=.TRUE.
!    nolapack=.FALSE.
    iph=pmi%iph
    ics=pmi%ics
    nrel=meqrec%nrel
!    if(mmdebug.ne.0) write(*,*)'MM meq_onephase 10: ',iph,ics
! for each phase "pmi" set eeccheck=0 at first interation
! THIS IS CURRNTLY NOT USED, will be added later
    if(meqrec%noofits.eq.1) then
       pmi%eeccheck=0
    elseif(meqrec%tpindep(1).or.meqrec%tpindep(2)) then
! if T or P not conditions set eeccheck=0 at each iteration
       pmi%eeccheck=0
    endif
! extract phase structure
!    write(*,*)'MM calling get_phase_data: ',iph
!    if(mmdebug.ne.0) then
!       write(*,*)'MM meq_onephase 12: ',iph,ics
!       gtpdebug=1
!    endif
! get_phase_data modified to ignore nonexisting composition sets
    call get_phase_data(iph,ics,nsl,nkl,knr,yarr,sites,qq,ceq)
!    if(mmdebug.ne.0) then
!       write(*,*)'MM meq_onephase 13: ',iph,ics,nsl,gx%bmperr
!       gtpdebug=0
!    endif
    if(gx%bmperr.ne.0) then
! handling of parallel by openMP
!$       if(.TRUE.) then
! this is written if parallel
!$          write(*,7)'get_phase_data error in meq_onephase: ',iph,ics,&
!$               omp_get_thread_num(),gx%bmperr
!$       else
! this is written if not parallel
          write(*,7)'get_phase_data error in meq_onephase: ',iph,ics,gx%bmperr
!$       endif
7      format(a,2i3,2x,2i5)
       goto 1000
    endif
! make sure all fractions >ymin and sums in all sublattices are equal to unity
    nochange=0
    ncc=0
!    if(mmdebug.ne.0) write(*,*)'MM meq_onephase 20: ',nsl,ncc
    do ll=1,nsl
       ysum=zero
       ncl=ncc
       do ik=1,nkl(ll)
          ncc=ncc+1
          if(yarr(ncc).lt.ymin) then
             if(test_phase_status_bit(iph,PHGAS)) then
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
!       write(*,*)'MM calling set_constitution 2:',ceq%eqno,iph,ics
       call set_constitution(iph,ics,yarr,qq,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'MM never error 17',iph,ics
! output if compiled with OpenMP
!$        write(*,*)'Thread :',ceq%eqname,omp_get_thread_num(),gx%bmperr
          goto 1000
       endif
    endif
!    if(mmdebug.ne.0) write(*,*)'MM meq_onephase 30: '
    if(test_phase_status_bit(iph,PHEXCB)) then
! If external charge balance phase matrix has one more line+column
       pmi%chargebal=1
       nd1=ncc+1
!       pmi%charge=qq(2)
       pmi%curd%netcharge=qq(2)
!       if(qq(2).gt.1.0D-8) write(*,*)'Charge: ',iph,ics,qq(2)
    else
       pmi%chargebal=0
       nd1=ncc
!       pmi%charge=zero
       pmi%curd%netcharge=zero
    endif
!--------------------------
! sublattice rows, nd2=nd1+1 because I use Lukas matrix inverter
    nd1=nd1+nsl
    nd2=nd1+1
! Allocate phase matrix, one extra dimension if external charge balance
! last column of pmat is left hand side ?? (reminicent from Lukas program)
!    allocate(pmat(nd1,nd2))
! pmat should be a square matrix
    allocate(pmat(nd1,nd1),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 37: ',errall
       gx%bmperr=4370; goto 1000
    endif
! return dimension of pmi%invmat
    if(pmi%idim.eq.0) then
       pmi%idim=nd1
       pmi%ncc=ncc
       allocate(pmi%invmat(nd1,nd1),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 38: ',errall
          gx%bmperr=4370; goto 1000
       endif
       pmi%invmat=zero
!       write(*,*)'Allocated invmat: ',nd1,ncc
! meqrec is not available in this routine ?? but meqrec%nrel passed in call
       allocate(pmi%xmol(nrel),stat=errall)
       allocate(pmi%dxmol(nrel,ncc),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 39: ',errall
          gx%bmperr=4370; goto 1000
       endif
!       write(*,*)'Allocated phase matrix: ',nd2,noel(),ncc
    endif
! value of RT should be moved before phase loop
    ceq%rtn=globaldata%rgas*ceq%tpval(1)
!--------------------------------------------------
! now treat different phase types
    call get_phase_variance(iph,nv)
    if(nv.eq.0) then
!------------------------------------- stoichiometric phase, fixed composition
! For stoichiometric phases calculate just G with T and P derivatives
! and driving force.  All pmi%dxmol=zero but one must also calculate 
! pmi%xmol and save it for all future iterations
! It must also be saved in curd%abnorm(1) ?? done in set_constitution ??
!       write(*,*)'MM xdone: ',pmi%xdone,iph,nv
!       if(mmdebug.ne.0) write(*,*)'MM meq_onephase 40: '
       if(pmi%xdone.eq.1) goto 90
! we must call set_constitution once to have correct abnorm etc
!       write(*,*)'MM calling set_constitution 3: ',iph,ics
       call set_constitution(iph,ics,yarr,qq,ceq)
       qsum=zero
       dqsum=zero
       pmi%xmol=zero
       pmi%dxmol=zero
       pmi%sumxmol=zero
       pmi%sumwmol=zero
       sumsit=zero
       do ll=1,nsl
          sumsit=sumsit+sites(ll)
       enddo
       kkk=0
       sublatt: do ll=1,nsl
          allconst: do ik=1,nkl(ll)
             kkk=kkk+1
             loksp=knr(kkk)
!             call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp,spextra)
             call get_species_component_data(loksp,nspel,ielno,&
                  stoi,spmass,qsp,ceq)
             if(gx%bmperr.ne.0) goto 1000
             addmol=zero
             do jz=1,nspel
                addmol(jz)=stoi(jz)
             enddo
             dqsum(kkk)=qsp
! 160820: forgotten to multiply with site ratio??!!
             qsum=qsum+sites(ll)*qsp
             do jz=1,nspel
                if(ielno(jz).gt.0) then
! ignore vacancies, taken care of by using sumsit=qq(1) above
                   pmi%dxmol(ielno(jz),kkk)=zero
                   pmi%xmol(ielno(jz))=pmi%xmol(ielno(jz))+&
                        sites(ll)*addmol(jz)
                endif
             enddo
          enddo allconst
       enddo sublatt
!       if(qsum.ne.zero) then
       if(abs(qsum).gt.1.0D-14) then
! if qsum not zero this phase should be suspended as it cannot be stable
          write(*,88)'Stoichiometric phase with net charge: ',iph,ics,qsum
88        format(a,2i4,2(1pe12.4))
       endif
! meqrec is not available in this routine ??
       do iz=1,nrel
          pmi%sumxmol=pmi%sumxmol+pmi%xmol(iz)
          pmi%sumwmol=pmi%sumwmol+pmi%xmol(iz)*mass_of(iz,ceq)
       enddo
! phase_varres(lokcs)%abnorm already set by set_constitution
       pmi%xdone=1
!
90     continue
! lokcs is set inside this subroutine
       call calcg(iph,ics,2,lokcs,ceq)
       if (gx%bmperr.ne.0) then
!          write(*,91)'calcg error in meq_onephase ',iph,gx%bmperr,ceq%eqno
91        format(a,3i5)
          goto 1000
       endif
!       if(mmdebug.ne.0) write(*,*)'MM meq_onephase 45: '
       eec1: if(globaldata%sysreal(1).gt.one) then
! EEC check for stoichiometric phases
! gval(1:6,1) are G, G.T, G.P, G.T.T, G.T.P, G.P.P
          yyyy=zero
          if(associated(meqrec%pmiliq)) then
! NOTE gval(2,1) is dG/dT i.e. the negative of entropy!!
!             write(*,*)'MM eec1A: DS',meqrec%seecliq,&
!                  pmi%curd%gval(2,1)/pmi%curd%abnorm(1)
             if(pmi%curd%gval(2,1)/pmi%curd%abnorm(1).lt.meqrec%seecliq) then
! too high entropy, set G=1.0 (avoid 0.0 ...)
                yyyy=pmi%curd%gval(1,1)
                pmi%curd%gval(1,1)=one
!                write(*,*)'MM eec1B: new G:',pmi%curd%gval(1,1),yyyy
             endif
          else
!             write(*,*)'MM eec1 No liquid entropy for stoichiometric phase!'
          endif
       endif eec1
! set the inverted phase matrix to zero !!!
       pmi%invmat=zero
!       do ik=1,ncc
!          pmi%invmat(ik,ik)=one
!       enddo
! maybe some common ending
       goto 900
    endif
!--------------------------------------------- zero some arrays, ideal phase
    pmi%xmol=zero
    pmi%dxmol=zero
    pmi%sumxmol=zero
    pmi%sumwmol=zero
    pmi%xdone=-1
!    if(phase_model(iph,ics,PHID,ceq)) then
!    write(*,*)'MM test ideal: ',test_phase_status_bit(iph,PHID)
    if(test_phase_status_bit(iph,PHID)) then
!--------------------------------------------- ideal phase (subst, no excess)
!       write(*,*)'Phase is ideal'
       if(test_phase_status_bit(iph,PHLIQ)) then
!          write(*,*)'MM liquid ideal: ',pmi%iph,pmi%ics
          meqrec%pmiliq=>pmi
       endif
!       if(mmdebug.ne.0) write(*,*)'MM meq_onephase 50: ideal'
! special treatment of ideal phase (gas), sites assumed to be unity
! 1. Calculate M_i and dM_i/dy^s_k and the net charge charge Q and dQ/dy^s_k
       pmi%xmol=zero
       pmi%dxmol=zero
       qsum=zero
       dqsum=zero
       ncon=0
       do ik=1,nkl(1)
          loksp=knr(ik)
!          call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp,spextra)
          call get_species_component_data(loksp,nspel,ielno,stoi,spmass,&
               qsp,ceq)
          if(gx%bmperr.ne.0) goto 1000
          addmol=zero
          do jk=1,nspel
             addmol(jk)=stoi(jk)
          enddo
          dqsum(ik)=qsp
          qsum=qsum+qsp*yarr(ik)
! It seems dxmol(element,constituent) is equal to the stoichiometry
! i.e. for a molecule H2O dM_H/dy_H2O=2; dM_O/dy_H2O=1, not 2/3 and 1/3
          do jk=1,nspel
             pmi%dxmol(ielno(jk),ik)=addmol(jk)
             pmi%xmol(ielno(jk))=pmi%xmol(ielno(jk))+addmol(jk)*yarr(ik)
          enddo
          ncon=ncon+1
       enddo
! meqrec is not available in this routine ??
       do ik=1,nrel
          pmi%sumxmol=pmi%sumxmol+pmi%xmol(ik)
!          write(*,*)'sumwmol 2: ',pmi%xmol(ik),mass_of(ik,ceq)
          pmi%sumwmol=pmi%sumwmol+pmi%xmol(ik)*mass_of(ik,ceq)
       enddo
! now calculate G and all 1st and 2nd derivatives
! This can be speeded up as all 2nd derivatives of constituents are RT/y
! The calculated values are used also in other parts of the code 
       call calcg(iph,ics,2,lokcs,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error calculating ideal gas',gx%bmperr,iph,ics
          goto 1000
       endif
       eec2: if(globaldata%sysreal(1).gt.one) then
! EEC check for ideal phases except gas
! gval(1:6,1) are G, G.T, G.P, G.T.T, G.T.P, G.P.P
          xxxx=pmi%curd%gval(2,1)/pmi%curd%abnorm(1)
          if(test_phase_status_bit(iph,PHLIQ)) then
             if(associated(meqrec%pmiliq)) then
! this is a second liquid
                if(xxxx.lt.meqrec%seecliq) then
                   meqrec%pmiliq=>pmi
                   meqrec%seecliq=xxxx
                endif
             else
! this is the first (or maybe only) liquid composition set
                meqrec%pmiliq=>pmi
                meqrec%seecliq=xxxx
             endif
          elseif(.not.test_phase_status_bit(iph,PHGAS)) then
             if(associated(meqrec%pmiliq)) then
! NOTE gval(2,1) is dG/dT i.e. the negative of entropy!!
                if(xxxx.lt.meqrec%seecliq) then
! G is set to -RT*ideal entropy/RT
!                   write(*,*)'MM eec2A: ',pmi%curd%gval(1,1)
                   pmi%curd%gval(1,1)=-pmi%curd%gval(2,1)
! no need to set other derivatives
                endif
             else
                write(*,*)'MM eec2 no liquid entropy to test!'
             endif
          endif
!          write(*,*)'MM eec2B: ',pmi%curd%gval(1,1)
       endif eec2
! calculate phase matrix elements
! temporarely ignore that the phase matrix is symmetric
! ceq%phase_varres(lokcs)%...
! gval(1:6,1) are G, G.T, G.P, G.T.T, G.T.P, G.P.P
! dgval(1,1:N,1) are first derivatives of G wrt constituent 1:N
! dgval(2,1:N,1) are second derivatives of G wrt constituent 1:N and T
! dgval(3,1:N,1) are second derivatives of G wrt constituent 1:N and P
! d2gval(ixsym(N*(N+1)/2),1) are 2nd derivatives of G wrt constituents N and M
! Last index is other properties than G like TC, BMAGN etc.
!       if(.not.nolapack) then
!          if(pmi%chargebal.eq.1) then
!             neq=ncon+ll+1
!             allocate(lapack(neq*(neq+1)/2))
!          else
!             neq=ncon+ll
!             allocate(lapack(neq*(neq+1)/2))
!          endif
!          lapack=zero
!       endif
       pmat=zero
! this is for an ideal phase with no excess
       do ik=1,nkl(1)
          do jk=ik,nkl(1)
!             ll=ixsym(ik,jk)
             ll=kxsym(ik,jk)
             pmat(ik,jk)=ceq%phase_varres(lokcs)%d2gval(ll,1)
             if(jk.gt.ik) pmat(jk,ik)=pmat(ik,jk)
!             if(.not.nolapack) lapack(ll)=ceq%phase_varres(lokcs)%d2gval(ll,1)
          enddo
       enddo
       neq=nkl(1)
!       write(*,770)(yarr(ik),ik=1,nkl(1))
!770    format('yfrac: ',4(1pe16.8))
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
! invert the phase matrix (faster routine should be used) IDEAL PHASE
! removed second argument
!       call mdinv(nd1,nd2,pmat,pmi%invmat,neq,ierr)
       call mdinv(nd1,pmat,pmi%invmat,neq,ierr)
       if(ierr.eq.0) then
          write(*,*)'MM Numeric problem 1, phase/set: ',iph,ics
          write(*,*)'Phase matrix singular 1:',pmi%iph,pmi%ics,pmi%ncc,ierr
          do jk=1,neq
             write(*,73)(pmat(ik,jk),ik=1,neq)
          enddo
73        format(1x,6(1pe12.4))
          gx%bmperr=4205; goto 1000
       endif
       goto 900
    endif
!---------------------------------------------- no analytical 2nd derivatives
! phases with models with no analytical second derivatives ....
!    if(phase_model(iph,ics,PHNODGDY2,ceq)) then
!    if(test_phase_status_bit(iph,PHNODGDY2,ceq)) then
    if(test_phase_status_bit(iph,PHNODGDY2)) then
!       write(*,*)'Models without 2nd derivatives not implemented'
       gx%bmperr=4206; goto 1000
    endif
!----------------------------------------------- ionic liquid phase
!    write(*,*)'MM test I2SL: ',test_phase_status_bit(iph,PHIONLIQ)
    ionliq: if(test_phase_status_bit(iph,PHIONLIQ)) then
!       write(*,*)'Warning; ionic liquid model not fully implemented'
! Calculate M_A and dM_A/dy_i taking into account that P and Q varies 
!   call get_phase_data(iph,ics,nsl,nkl,knr,yarr,sites,qq,ceq)
       if(test_phase_status_bit(iph,PHLIQ)) then
          meqrec%pmiliq=>pmi
!          write(*,*)'MM liquid ionic: ',pmi%iph,pmi%ics
       endif
!       if(mmdebug.ne.0) write(*,*)'MM meq_onephase 55: '
       pmi%ionliq=nkl(1)
       pmi%xmol=zero
       pmi%dxmol=zero
       qsum=zero
       dqsum=zero
       pmi%sumxmol=zero
       pmi%sumwmol=zero
       allocate(sumion(nrel,2),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 40: ',errall
          gx%bmperr=4370; goto 1000
       endif
!       pmi%sumiliq=zero
! end extra
       ncon=0
       sumion=zero
       yva=zero
!       write(*,217)'y:  ',ncc,(yarr(ik),ik=1,ncc)
       i2sly=nkl(1)+nkl(2)+1
       do ll=1,nsl
          do ik=1,nkl(ll)
             ncon=ncon+1
             loksp=knr(ncon)
!             pmi%ikon(ncon)=loksp
! if only neutrals we can have a single wildcard in first sublattice ...
             if(loksp.lt.0) then
                if(ll.eq.1 .and. nkl(1).eq.1) cycle
                write(*,*)'Illegal wildcard constituent in ionic liquid'
                gx%bmperr=4197; goto 1000
             endif
             if(btest(pmi%curd%constat(ncon),CONVA)) then
! This is the nypothetical vacancy .... its charge is sites(2) = Q
                yva=yarr(ncon)
! save its index in isly(1), otherwise that is number of constit+1
                i2sly(1)=ncon
!                pmi%valency(ncon)=sites(2)
!                write(*,*)'Va: ',ncon,yva
             else
!               call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp,spextra)
                call get_species_component_data(loksp,nspel,ielno,stoi,&
                     spmass,qsp,ceq)
                if(gx%bmperr.ne.0) goto 1000
! i2sly is index of first neutral (if any) otherwise number of constit+1
                if(qsp.eq.zero .and. i2sly(2).gt.ncon) i2sly(2)=ncon
!                write(*,*)'Species: ',ncon,i2sly,qsp
!                if(qsp.eq.zero .and. i2sly(2).eq.0) i2sly(2)=ncon
!                pmi%valency(ncon)=abs(qsp)
!                write(*,*)'charge: ',ncon,qsp
                do jk=1,nspel
                   notva: if(ielno(jk).gt.0) then
! ignore vacancies in species
                      qsp=sites(ll)*stoi(jk)
                      pmi%dxmol(ielno(jk),ncon)=qsp
                      pmi%xmol(ielno(jk))=pmi%xmol(ielno(jk))+qsp*yarr(ncon)
                      sumion(ielno(jk),ll)=sumion(ielno(jk),ll)+&
                           stoi(jk)*yarr(ncon)
! take into account that the site ratios depend on constitition in corrion_..
!                      write(*,21)'ddMA:',jk,ielno(jk),ncon,ll,&
!                           pmi%dxmol(ielno(jk),ncon),qsp,sites(ll),stoi(jk)
!21                    format(a,4i3,4(1pe12.4))
! sums used in calc_dgdyterms1 to handle that sites(ll) depend on constitition
!                      pmi%sumiliq(ielno(jk),ll)=pmi%sumiliq(ielno(jk),ll)+&
!                              stoi(jk)*yarr(ncon)
! Hm, the statement above not necessary as below it is already included ....
                   endif notva
                enddo
             endif
          enddo
       enddo
! save these as needed in calc_dgdyterms
! i2sly(1) is index of vacancy, if no vacancy equal to #of constituents+1
! i2sly(2) is index if first neutral, if no neutal equal to #of constituents+1
       pmi%i2sly=i2sly
       pmi%yva=yva
! zero matrix
       pmat=zero
!...........................................
!       goto 261
! now handle that site ratios depend on constituent fractions
! (maybe also that the formula unit depend on composition)
! phlista(lokph)%i2slx; lokph=pmi%curd%phlink
! BUT: phlista is private ....
! M_A = P*M'_A + Q*M"_A           M'_A and M"_A are in sumion(A,1:2))
! P=\sum_j (-v_j)y_j + Qy_Va      j is anion
! Q=\sum_i v_iy_i                 i is cation
!       if(mmdebug.ne.0) write(*,*)'MM meq_onephase 60: '
       icon=0
       do ik=1,nkl(1)
          icon=icon+1
          do jk=1,nrel
! for cations: extra dM_A/dyi = v_i*y_Va*M'_A + v_i*M"_A where i is cation
             qsp=pmi%curd%dpqdy(icon)*(yva*sumion(jk,1)+sumion(jk,2))
! note dxmol(jk,icon) has been multiplied with sites(1) above ....
             pmi%dxmol(jk,icon)=pmi%dxmol(jk,icon)+qsp
          enddo
       enddo
! i2sly(1) is index of vacancy, i2sly(2) index of first neutral
! If no vacancy or no neutral the corresponding i2sly is ncc+1
       do ik=1,nkl(2)
          icon=icon+1
          if(icon.lt.min(i2sly(1),i2sly(2))) then
             do jk=1,nrel
! for anions: extra dM_A/dyj = (-v_j)*M'_A where j is anion
                qsp=pmi%curd%dpqdy(icon)*sumion(jk,1)
                pmi%dxmol(jk,icon)=pmi%dxmol(jk,icon)+qsp
!              write(*,654)'Extra term anjon:  ',jk,icon,pmi%dxmol(jk,icon),qsp
             enddo
          else
! note icon not updated correctly if neutrals, use ncon below
             exit
          endif
       enddo
! take care of a vacancy
       if(icon.eq.i2sly(1)) then
          do jk=1,nrel
! for Va: extra dM_A/dyi = Q*M'_A where i is vacancy
             pmi%dxmol(jk,icon)=sites(2)*sumion(jk,1)
!             write(*,654)'Extra term for Va: ',jk,icon,&
!                  pmi%dxmol(jk,icon),sites(2)*sumion(jk,1)
          enddo
       endif
! Derivatives with respect to neutrals have no extra term
!       do jk=1,nrel
!          write(*,217)'dMA2:',jk,(pmi%dxmol(jk,ik),ik=1,ncc)
!       enddo
! one may exit loop above with different values of ncon and icon, 
! ncon is the total number of constituents
       icon=ncon
!......................................... end handling P and Q variation
261    continue
! meqrec is not available in this routine ??
       do ik=1,nrel
          pmi%sumxmol=pmi%sumxmol+pmi%xmol(ik)
          pmi%sumwmol=pmi%sumwmol+pmi%xmol(ik)*mass_of(ik,ceq)
       enddo
! now calculate G and all 1st and 2nd derivatives
! The calculated values are used also in other parts of the code 
       call calcg(iph,ics,2,lokcs,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'MM Error calculating G 1: ',iph,ics,lokcs
          goto 1000
       endif
! correction of I2SL second derivatives due to variation of P and Q
       if(meqrec%noofits.gt.1) then
! NOTE pmat is dimensioned pmat(nd1,nd2)
          call corriliq_d2gdyidyj(nkl,knr,ceq%cmuval,pmi,ncon,nd1,pmat,ceq)
          if(gx%bmperr.ne.0) goto 1000
       endif
       eec3: if(globaldata%sysreal(1).gt.one) then
! EEC check for ionic liquid phase (no need to test for PHLIQ)
! gval(1:6,1) are G, G.T, G.P, G.T.T, G.T.P, G.P.P
          xxxx=pmi%curd%gval(2,1)/pmi%curd%abnorm(1)
          if(associated(meqrec%pmiliq)) then
! we already have a liquid 
             if(xxxx.gt.meqrec%seecliq) then
! this liquid has higher entropy                
                meqrec%pmiliq=>pmi
                meqrec%seecliq=xxxx
             endif
          else
! save link to liquid with higest entropy
             meqrec%pmiliq=>pmi
             meqrec%seecliq=xxxx
          endif
!          write(*,*)'MM eec3: ',meqrec%seecliq,associated(meqrec%pmiliq)
       endif eec3
!       write(*,17)'pots: ',(ceq%cmuval(ik),ik=1,3)
!       do ll=1,nd1
!          write(*,17)'cion: ',(pmat(ll,ik),ik=1,nd1)
!       enddo
! calculate phase matrix elements, the second derivatives
! note pmat has some contributions above ??
       neq=icon
       fion=one
       do ik=1,icon
          do jk=ik,icon
             pmat(ik,jk)=fion*pmat(ik,jk)+&
                  ceq%phase_varres(lokcs)%d2gval(kxsym(ik,jk),1)
!                  ceq%phase_varres(lokcs)%d2gval(ixsym(ik,jk),1)
! remove next line when using a routine inverting a symmetric matrix
             if(jk.gt.ik) pmat(jk,ik)=pmat(ik,jk)
          enddo
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
       enddo
       neq=neq+nsl
!       write(*,65)'pdim: ',nd1,nd2,neq,ncon,icon,nsl,(nkl(ll),ll=1,nsl)
!65     format(a,6i4,10i3)
!       do ll=1,nd1
!          write(*,17)'pmat: ',(pmat(ll,ik),ik=1,nd1)
17        format(a,6(1pe12.4))
!       enddo
! invert the phase matrix (faster routine should be used) IONIC LIQUID MODEL
!       call mdinv(nd1,nd2,pmat,pmi%invmat,nd1,ierr)
!       write(*,*)'Value 2 of nolapsck: ',nolapack,.not.nolapack
! removed 2nd argument
       call mdinv(nd1,pmat,pmi%invmat,nd1,ierr)
       if(ierr.eq.0) then
!          write(*,*)'MM Numeric problem 2, phase/set: ',iph,ics
          write(*,*)'Phase matrix singular 2:',pmi%iph,pmi%ics,pmi%ncc,ierr
          gx%bmperr=4205; goto 1000
       endif
!       do ll=1,nd1
!          write(*,17)'pinv: ',(pmi%mat(ll,ik),ik=1,nd1)
!       enddo
! maybe some common ending
       goto 900
    endif ionliq
!------------------------------------------------- all other phase models (CEF)
! For all other phases calculate G and all first and second derivatives
! for current composition
300 continue
!    write(*,*)'MM CEF phase',ceq%eqno
! Calculate M_i and dM_i/dy^s_k and the net charge charge Q and dQ/dy^s_k
!   call get_phase_data(iph,ics,nsl,nkl,knr,yarr,sites,qq,ceq)
! how to normalize xmol?  use qq(1)!!, it handels vacancies .... ????
!    write(*,*)'MM Phase 1: ',pmi%iph,pmi%ics
!    if(test_phase_status_bit(iph,PHLIQ)) then
!       write(*,*)'MM liquid other: ',pmi%iph,pmi%ics
!       meqrec%pmiliq=>pmi
!    endif
!    if(mmdebug.ne.0) write(*,*)'MM meq_onephase 70: '
    sumsit=one
    pmi%xmol=zero
    pmi%dxmol=zero
    qsum=zero
    dqsum=zero
    ncon=0
    pmi%sumxmol=zero
    pmi%sumwmol=zero
    subll: do ll=1,nsl
       constll: do ik=1,nkl(ll)
          ncon=ncon+1
          loksp=knr(ncon)
!          call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp,spextra)
          call get_species_component_data(loksp,nspel,ielno,stoi,spmass,qsp,ceq)
          if(gx%bmperr.ne.0) goto 1000
          addmol=zero
          do jk=1,nspel
             addmol(jk)=stoi(jk)
          enddo
          dqsum(ncon)=sites(ll)*qsp
          qsum=qsum+sites(ll)*qsp*yarr(ncon)
          do jk=1,nspel
!             write(*,963)'xmol: ',ncon,ik,jk,ielno(jk),sites(ll)
!963          format(a,4i3,6(1pe12.4))
             if(ielno(jk).gt.0) then
! ignore vacancies
! addmol(jk) can be replaced by stoi(jk) when I know it works ....
                pmi%dxmol(ielno(jk),ncon)=sites(ll)*addmol(jk)
                pmi%xmol(ielno(jk))=pmi%xmol(ielno(jk))+&
                     sites(ll)*addmol(jk)*yarr(ncon)
             endif
          enddo
       enddo constll
    enddo subll
!    write(*,*)'MM segmentation fault test 1',nrel
! meqrec is not available in this routine ??
    do ik=1,nrel
       pmi%sumxmol=pmi%sumxmol+pmi%xmol(ik)
!       write(*,*)'sumwmol 3:',pmi%xmol(ik),mass_of(ik,ceq)
       pmi%sumwmol=pmi%sumwmol+pmi%xmol(ik)*mass_of(ik,ceq)
    enddo
!    write(*,*)'MM segmentation fault test 2'
!    write(*,92)'onephase 3: ',pmi%iph,nsl,pmi%xdone,pmi%sumxmol,qq(1)
!92  format(a,3i3,6(1pe12.4))
!    write(*,17)'Vacanies: ',qq
!       do i=1,noel()
!          write(*,17)'xm: ',pmi%xmol(i)
!          write(*,17)'dxm: ',(pmi%dxmol(i,j),j=1,ncon)
!       enddo
! now calculate G and all 1st and 2nd derivatives
! The calculated values are stored and used also in other parts of the code 
!    write(*,*)'MM segmentation fault test 3',iph,ics
    call calcg(iph,ics,2,lokcs,ceq)
    if(gx%bmperr.ne.0) then
       write(*,11)'MM Error calculating G 2: ',iph,ics,lokcs,gx%bmperr
11     format(a,5i5)
       goto 1000
    endif
!    if(mmdebug.ne.0) write(*,*)'MM meq_onephase 80: '
!    write(*,*)'MM segmentation fault 10',globaldata%sysreal(1)
    eec4: if(globaldata%sysreal(1).gt.one) then
! check of EEC for a CEF phase
!       if(pmi%eeccheck.eq.0) then
! This is first iteration or we have variable T or P
       xxxx=pmi%curd%gval(2,1)/pmi%curd%abnorm(1)
!       write(*,*)'MM eec4A: ',meqrec%noofits,associated(meqrec%pmiliq)
       if(test_phase_status_bit(iph,PHLIQ)) then
! This is a liquid phase
          if(associated(meqrec%pmiliq)) then
! We have several liquids, take the highest entropy (note xxx is -entropy!)
             if(xxxx.lt.meqrec%seecliq) then
                meqrec%pmiliq=>pmi
                meqrec%seecliq=xxxx
             endif
!             write(*,*)'MM eec4B: second liquid'
          else
! this is the first (or maybe only) liquid composition set
             meqrec%pmiliq=>pmi
             meqrec%seecliq=xxxx
          endif
!          write(*,'(a,l2,5(1pe12.4))')'MM eec4C: liq:',&
!               associated(meqrec%pmiliq),meqrec%seecliq,&
!               pmi%curd%gval(2,1),pmi%curd%abnorm(1)
       elseif(.not.test_phase_status_bit(iph,PHGAS)) then
! this is a condensed phase which should have its entropy checked
! NOTE gval(2,1) is dG/dT i.e. the negative of entropy!!
          if(xxxx.lt.meqrec%seecliq) then
!             write(*,*)'MM eec4D S(solid)>S(liquid)',-xxxx,-meqrec%seecliq
! replace G and all derivates with a phase with just configurational entropy
! in the pmi%curd%gval, pmi%curd%dgval and  pmi%curd%d2gval
             yyyy=pmi%curd%gval(1,1)
             call calc_eec_gibbsenergy(pmi%curd,ceq)
             if(gx%bmperr.ne.0) goto 1000
          endif
       endif
!       write(*,'(a,5(1pe12.4))')'MM eec4F: ',pmi%curd%gval(1,1),yyyy,&
!            -meqrec%seecliq,-xxxx
!    else
!       write(*,*)'MM no EEC'
    endif eec4
! calculate phase matrix elements, first and second derivatives
!    write(*,*)'MM segmentation fault 19'
    pmat=zero
    neq=ncon
!    write(*,*)'MM segmentation fault 20'
! here we are calculating CEF models
    do ik=1,ncon
! OK       jxsym=ixsym(ik,ik); kxsym=0
       jxsym=ixsym(ik,ik)
       do jk=ik,ncon
! fatal parallel execution frequently here ... why?? Error message:
! index '0' of dimension 1 of array 'ceq' below lower bound of 1
!          pmat(ik,jk)=ceq%phase_varres(lokcs)%d2gval(ixsym(ik,jk),1)
! modified code:
          ll=kxsym(ik,jk)
!          ll=ixsym(ik,jk)
! OK          jxsym=jxsym+kxsym; kxsym=jk
! increment jxsym at the end of the loop ...
! testing replacing ixsym .... too complicated ...
          if(ll.ne.jxsym) then
!             write(*,*)'Problems: ',ik,jk,ll,jxsym
             stop "Problemns with ixsym"
!          else
!             write(*,*)'No problems: ',ik,jk,ll,jxsym,kxsym
          endif
! attempt to avoid a crash
!$          if(lokcs.le.0 .or. ll.le.0) then
!$             write(*,491)'meq_onephase error: ',lokcs,ll,omp_get_thread_num()
491          format(' *** ',a,4i5)
!$             goto 1000
!$          endif
          pmat(ik,jk)=ceq%phase_varres(lokcs)%d2gval(ll,1)
!          if(.not.nolapack) lapack(ll)=ceq%phase_varres(lokcs)%d2gval(ll,1)
! remove next line when using an inversion for symmetric matrix
          if(jk.gt.ik) pmat(jk,ik)=pmat(ik,jk)
! this is an attempt to avoid calling ixsym ... it works
          jxsym=jxsym+jk
       enddo
!       write(*,17)'row2A: ',(pmat(ik,jj),jj=1,nd1)
    enddo
! Then set the sublattice elements
!    write(*,*)'MM segmentation fault 20'
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
! if external charge balance add one column and one row
! It causes problem to invert the phase matrix below for a phase like
! M2O3 with cations CE+3 and LA+3 as the phase is always neutral 
! and the charge balance not needed.
       neq=neq+1
       do jk=1,ncon
! this is the row
          pmat(jk,neq)=dqsum(jk)
! this is the column
          pmat(neq,jk)=dqsum(jk)
       enddo
    endif
! write the phase matrix on a file
!    open(33,file='phasemat.dat ',access='sequential',status='unknown')
!    write(33,*)'Phase matrix',nd1
!    do jk=1,nd1
!       write(33,111)jk,(pmat(jk,ll),ll=1,nd1)
111    format('>',i4,1x,4(1pe15.6))
!    enddo
! debug output
!    write(*,*)'Phase matrix',nd1,neq,pmi%chargebal
!    do j=1,neq
!       write(*,17)'pmat: ',(pmat(i,j),i=1,neq)
!    enddo
! invert the phase matrix (using LAPACK+BLAS ... 50% faster than with Leo)
! removed 2nd argument
!    call mdinv(nd1,nd2,pmat,pmi%invmat,neq,ierr)
!    write(*,*)'MM segmentation fault 30'
    call mdinv(nd1,pmat,pmi%invmat,neq,ierr)
    if(ierr.eq.0) then
       write(*,*)'MM Numeric problem 3, phase/set:',iph,ics
!       if(ocv()) write(*,556)'Phase matrix singular 3:',meqrec%noofits,&
       if(pmi%chargebal.eq.1) then
! can be problem with external chargebalance not needed ...
          call get_phase_name(pmi%iph,1,name)
          write(*,553)'Try to suspend phase: ',trim(name)
553       format(a,a)
       endif
556    format(a,6i5)
! emergency fix does not work ...
       pmi%invmat=zero
       do jk=1,neq
          pmi%invmat(jk,jk)=one/neq
       enddo
!       do jk=1,neq
!          write(*,18)'3Y mat:',jk,(pmat(ik,jk),ik=1,neq)
!       enddo
!       do jk=1,neq
!          write(*,18)'3Y inv:',jk,(pmi%invmat(ik,jk),ik=1,neq)
!       enddo
18     format(a,i3,7(1pe10.2))
!       do jk=1,neq
!          write(*,73)(pmat(ik,jk),ik=1,neq)
!       enddo
       gx%bmperr=4205; goto 1000
    endif
    goto 900
!-------------------------------------------
900 continue
!
!    if(mmdebug.ne.0) write(*,*)'MM meq_onephase exit: '
    goto 1000
!--------------------------------------------
1000 continue
!    write(*,*)'MM exit meq_onephase'
    return
  end subroutine meq_onephase !ixsym
 
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine corriliq_d2gdyidyj
!\begin{verbatim}
  subroutine corriliq_d2gdyidyj(nkl,knr,curmu,pmi,ncc,nd1,pmat,ceq)
! correction of d2G/dy1dy2 for ionic liquid because the formula unit is
! not fixed.  This contributes ONLY to the second derivaties of G and
! is not really part of the model itself, only needed when minimizing G
    implicit none
    type(gtp_equilibrium_data), pointer :: ceq
    TYPE(meq_phase), pointer :: pmi
    integer ncc,nd1,nkl(*),knr(*)
    double precision curmu(*),pmat(nd1,*)
!\end{verbatim}
! corr = \sum_A \mu_A*d2(N_A)/dy_i/dy_k ; i cation, k cation, anion, Va
! N_A  = P*\sum_i b_Ai y_i + Q(\sum_j b_Aj y_j + ... ) b_Ai stoich.fact. of A
! P    = \sum_j v_j y_j + y_Va Q
! Q    = \sum_i v_i y_i 
!
! Derivativs of P and Q
! dP/dy_i = y_Va v_i;   dP/dy_j = v_j;  dP/dy_Va = Q
! dQ/dy_i = v_i         dQ/dy_j = zero  dQ/dy_Va = zero
! d2P/dy_idy_Va = v_i
! 
! d(N_A\mu_A)/dy_i = dP/dy_i\sum_jb_Aj + v_i
!
    integer icon,jcon,loksp,nspel,ielno(10),el,allions,nobug
    double precision stoi(10),spmass,qsp1,qsp2,add1,add2,yva,sumcat,bug
    double precision bugfix
!tafidbug
!    write(*,*)'Skipping liquid correction'
!    goto 1000
! this correction term affects only second derivatives and thus convergence 
! speed and stability.  But it seems just to mess up everything.
!
! dpqdy(1..ncc) is the absolute value of the charge of the species
! It is not used as we must get species data, better not to use ...
! i2sly(1) is index of vacancy, i2sly(2) is index of first neutral
! If either is missing it is equal to number of constituents+1
    allions=min(pmi%i2sly(1),pmi%i2sly(2))
!    write(*,12)'mu: ',(curmu(i1),i1=1,noel())
12  format(a,6(1pe12.4))
    if(nkl(1).eq.0) then
! no cations (bor anions), only neutrals, no need to calculate anything      
!       write(*,*)'Liquids without cations have fixed stoichiometry 1.0
       goto 1000
    endif
! If there are vacancies we save its fraction here, if not set to zero
!    if(pmi%i2sly(1).lt.ncc) then
    if(pmi%i2sly(1).le.ncc) then
       yva=pmi%curd%yfr(pmi%i2sly(1))
    else
       yva=zero
    endif
!    write(*,11)'corrion 1: ',yva,pmi%i2sly,nkl(1)+nkl(2),allions,ncc
11  format(a,1pe12.4,10i5)
! to simplify testing, 0 means include contribution from pairs of cations
    nobug=0
    bugfix=one
    sumcat=zero
! just loop for all cations here. Inside this loop we step jcon
!  for all constituents up to vacancies or last anion.
    do icon=1,nkl(1)
!    icon=0
!    do i1=1,nkl(1)
!    do i1=1,allions-1
! loop for all cations and anions
!       icon=icon+1
       loksp=knr(icon)
!       call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp1,spextra)
       call get_species_component_data(loksp,nspel,ielno,stoi,spmass,qsp1,ceq)
       if(gx%bmperr.ne.0) goto 1000
       add2=zero
       do el=1,nspel
! skip any vacancy in a species, they have zero chemical potential anyway
          if(ielno(el).gt.0) add2=add2+stoi(el)*curmu(ielno(el))
       enddo
       add1=add2
!       write(*,13)'first cat: ',icon,0,qsp1,add1
13     format(a,2i3,6(1pe12.4))
!-------------------------2nd derivatives wrt two cations
       jcon=icon
       do while(jcon.le.nkl(1))
! loop for all pairs of cations incl twins, nkl(1) is number of cations
! A smart but messy solution is to skip this loop for jcon=icon ...
          loksp=knr(jcon)
!          call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp2,spextra)
          call get_species_component_data(loksp,nspel,ielno,stoi,&
               spmass,qsp2,ceq)
          if(gx%bmperr.ne.0) goto 1000
          add2=zero
          do el=1,nspel
             if(ielno(el).ne.0) add2=add2+stoi(el)*curmu(ielno(el))
          enddo
          bug=add2
! sumcat is used below for derivative wrt cation and vacancy
          if(icon.eq.1) then
             sumcat=sumcat+pmi%curd%yfr(jcon)*add2
!                write(*,13)'sumcat:    ',0,jcon,yva,pmi%curd%yfr(jcon),&
!                     add2,sumcat
          endif
! if there are no vacancies the derivative of P is zero wrt two cations
! this is \sum_A dP/dy_icon*b_Ajcon*mu_A+\sum_A dP/dy_jcon*b_Aicon*mu_A
          if(nobug.eq.0 .and. yva.gt.zero) then
             add2=bugfix*yva*(qsp1*add2+qsp2*add1)
!             if(abs(yva*(add2)).gt.1.0D2) then
! This is a sensitive point for convergence, values of 1.0D+33 found !!!
! But bad converge also when small values, less than 100
!                add2=-1.0D2
!             endif
!             write(*,13)'pmat caca: ',icon,jcon,qsp1,yva,bug,add2
! store value in pmat as correction to d2G/dyidyj
             pmat(icon,jcon)=-add2
! tafidbug 2
!             pmat(icon,jcon)=add2
          endif
          jcon=jcon+1
       enddo
! ------------------------ 2nd derivative wrt to cation and anion
       do while(jcon.lt.allions)
! loop for all anions, allions-1 is last anion
          loksp=knr(jcon)
!          call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp2,spextra)
          call get_species_component_data(loksp,nspel,ielno,stoi,&
               spmass,qsp2,ceq)
          if(gx%bmperr.ne.0) goto 1000
          add2=zero
          do el=1,nspel
             if(ielno(el).ne.0) add2=add2+stoi(el)*curmu(ielno(el))
          enddo
          bug=add2
! This is \sum_A dP/dy_jcon*b_Aicon*mu_A+\sum_A dQ/dy_icon*b_Ajcon*mu_A
! Note dP/dy = -qsp2 as qsp2 is negative
          add2=qsp1*add2-qsp2*add1
!          write(*,13)'pmat caan: ',icon,jcon,qsp2,bug,add2
! store value in pmat as correction to d2G/dyidyj
          pmat(icon,jcon)=-add2
! tafidbug 2
!          pmat(icon,jcon)=add2
          jcon=jcon+1
       enddo
!------------- second derivative wrt cation and vacancy
!       if(icon.le.nkl(1) .and. jcon.eq.pmi%i2sly(1)) then
       if(jcon.le.ncc .and. jcon.eq.pmi%i2sly(1)) then
! if no vacancy then i2sly(1)=ncc+1
! This is \sum_A d2P/dy_icon dy_Va*\sum_k y_k*b_Ak*\mu_A + Q * b_Aicon*\mu_A
          add2=qsp1*sumcat+pmi%curd%sites(2)*add1
! It think the line above is correct but the one below works better ...
!          add2=qsp1*sumcat
!          write(*,13)'pmat cava: ',icon,jcon,qsp1,&
!               sumcat,pmi%curd%sites(2),add1,add2
! store value in pmat as correction to d2G/dyidyj
          pmat(icon,jcon)=-add2
! tafidbug 2
!          pmat(icon,jcon)=add2
          jcon=jcon+1
       endif
!------------- second derivative wrt cation and neutral
! is this really correct??
       do while(jcon.le.ncc)
          loksp=knr(jcon)
!          call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp2,spextra)
          call get_species_component_data(loksp,nspel,ielno,stoi,spmass,&
               qsp2,ceq)
          if(gx%bmperr.ne.0) goto 1000
          add2=zero
          do el=1,nspel
             if(ielno(el).ne.0) add2=add2+stoi(el)*curmu(ielno(el))
          enddo
          bug=add2
! This is \sum_A dQ/dy_icon * b_Ajcon * mu_A, icon is cation and jcon neutal
          add2=qsp1*add2
!          write(*,13)'pmat cane: ',icon,jcon,qsp1,bug,add2
          pmat(icon,jcon)=-add2
! tafidbug 2
!          pmat(icon,jcon)=add2
          jcon=jcon+1
       enddo
!------------- no other terms
    enddo
!    write(*,*)'Correction to phase matrix from corriliq: ',&
!         pmi%curd%phtupx,nobug
!    do icon=1,ncc
!       write(*,1100)(pmat(icon,jcon),jcon=1,ncc)
!    enddo
1100 format(6(1pe12.4))
1000 continue
    return
  end subroutine corriliq_d2gdyidyj

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function same_composition
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
    TYPE(meq_phase), dimension(*) :: phr
    TYPE(meq_setup) :: meqrec
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer jp,jy
! If the difference is larger than xdiff then the compositions are not the same
!    double precision, parameter :: xdiff=0.01D0
! FINETUNING: a large value of xdiff may mean you miss a miscibility gap
! a small value may create bad convergence
! 0.05 fails to find L1_2/A1/L1_0 in Au-Cu ...
!    double precision, parameter :: xdiff=0.05D0
! 0.01 works better for Au-Cu ... maybe other problems ...
    double precision, parameter :: xdiff=0.01D0
    double precision, dimension(maxel) :: xmol1,xmol2,wmass
    double precision amount,totmol,totmass,xdiffm,xdiffc
! CCI
    same_composition=.FALSE.
! check if any other compset of the phase stable with same composition
    call calc_phase_molmass(phr(jj)%iph,phr(jj)%ics,xmol1,wmass,&
         totmol,totmass,amount,ceq)
    if(gx%bmperr.ne.0) goto 1000
    xdiffm=one
!    write(*,*)'MM testing same composition',jj,phr(jj)%iph,phr(jj)%ics
! ?? strange loop limits ??
!    do jp=jj-1,1,-1
    do jp=1,meqrec%nphase
       if(phr(jp)%iph.eq.phr(jj)%iph) then
          if(phr(jp)%stable.eq.1) then
             call calc_phase_molmass(phr(jp)%iph,phr(jp)%ics,xmol2,wmass,&
                  totmol,totmass,amount,ceq)
             if(gx%bmperr.ne.0) goto 1000
             do jy=1,meqrec%nrel
                xdiffc=abs(xmol1(jy)-xmol2(jy))
                if(xdiffc.lt.xdiffm) then
                   xdiffm=xdiffc
                endif
                if(xdiffc.gt.xdiff) goto 110
!                if(abs(xmol1(jy)-xmol2(jy)).gt.xdiff) goto 110
             enddo
! we have found another stable composition set with same composition
             goto 300
          endif
       elseif(phr(jp)%iph.lt.phr(jj)%iph) then
          cycle
       else
          exit
       endif
110    continue
    enddo
    same_composition=.FALSE.
    goto 1000
!--------------------------------------------------------
! we found a stable composition set with the same composition
300 continue
    same_composition=.TRUE.
    if(ocv()) write(*,117)'Not added comp.set phase: ',phr(jj)%iph,&
         phr(jj)%ics,phr(jp)%ics,xdiffm
117 format(a,i3,2i4,2x,1pe12.4)
! One cannot have two composition sets with same composition.
! try to reset this composition set to default constition
    call set_default_constitution(phr(jj)%iph,phr(jj)%ics,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    goto 1000
!
1000 continue
    return
  end function same_composition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_dgdyterms1
!\begin{verbatim}
  subroutine calc_dgdyterms1(nrel,ia,tpindep,mamu,mag,mat,map,pmi,&
       curmux,noofits)
! THIS SUBROUTINE IS NO LONGER USED!!  Cannot be used in parallel
! any change must also be made in subroutine calc_dyterms2 and calc_dgdytermsh
! calculate the terms in the deltay expression for amounts of component ia
!
! DM_A = \sum_B mu_B*MAMU(B) - MAG - MAT*dt - MAP*dp
!
! where MAMU=\sum_i dM_A/dy_i*\sum_j invmat(i,j)*dM_B/dy_j
!       c_iB=\sum_j invmat(i,j)*dM_B/dy_j etc etc
!
! it may not be very efficient but first get it right ....
! tpindep(1) is TRUE if T variable, tpindep(2) is TRUE if P are variable
!
! >>> ATTENTION, there is a FASTER VERSION calc_dgdyterms1X
! >>> ATTENTION not safe for parallelization ....
!
    implicit none
    integer ia,nrel,noofits
    logical tpindep(2)
    double precision, dimension(*) :: mamu
    double precision mag,mat,map
    double precision curmux(*)
! pmi is the phase data record for this phase
    type(meq_phase), pointer :: pmi
!\end{verbatim} %+
! these are to be multiplied with mu(ib), nothing, deltaT, deltaP
    integer iy,jy,ib,nocon,errall
! initial values for saved results OLD VERSION
    integer :: sameit=0,big1p=0,big2p=0,big1n=0,big2n=0
    double precision cig,cit,cip,haha
    double precision morr
    double precision, allocatable, dimension(:) :: zib
! ATTENTION, see calc_dgdyterms1X !!!! for better routine
    double precision, allocatable, dimension(:,:) :: maybesave
    double precision, allocatable, dimension(:,:) ::  save1
    double precision, allocatable, dimension(:,:) ::  save2
! NOTE THIS SUBROUTINE IS NO LONGER USED!!
    save sameit,big1p,big1n,big2p,big2n
    save save1,save2
    logical big
!
!-----------
! \sum_i \sum_j e_ij*dM_A/dy_i dG/dy_j
! skip code for saving as that is implemented in calc_dgdytermes1X
!    write(*,*)'Using calc_dgdyterms1 without saving'
    goto 100
! code below to be ignored
!
    if(noofits.ne.sameit) then
! new iteration, discard saved values
       big1p=0; big1n=0
       big2p=0; big2n=0
       sameit=noofits
       goto 100
    endif
! use save values for the phases with many constituents
!                if(test_phase_status_bit(phasetuple(phr(jj)%iph)%ixphase,&
    if(10*pmi%iph+pmi%ics.eq.big1p) then
!       write(*,13)'MM using saved values 1:',noofits,sameit,big1p,big1n,ia
13     format(a,2i5,5x,2i5,5x,3i5)
       mag=zero
       mat=zero
       map=zero
       do ib=1,nrel
          mamu(ib)=zero
       enddo
       do iy=1,big1n
          morr=pmi%dxmol(ia,iy)
          do ib=1,nrel
             mamu(ib)=mamu(ib)+save1(ib,iy)*morr
          enddo
          mag=mag+save1(nrel+1,iy)*morr
          if(tpindep(1)) mat=mat+save1(nrel+2,iy)*morr
          if(tpindep(2)) map=map+save1(nrel+3,iy)*morr
       enddo
       goto 1000
    elseif(10*pmi%iph+pmi%ics.eq.big2p) then
!       write(*,13)'MM using saved values 2:',noofits,sameit,big2p,big2n,ia
       mag=zero
       mat=zero
       map=zero
       do ib=1,nrel
          mamu(ib)=zero
       enddo
       do iy=1,big2n
          morr=pmi%dxmol(ia,iy)
          do ib=1,nrel
             mamu(ib)=mamu(ib)+save2(ib,iy)*morr
          enddo
          mag=mag+save2(nrel+1,iy)*morr
          if(tpindep(1)) mat=mat+save2(nrel+2,iy)*morr
          if(tpindep(2)) map=map+save2(nrel+3,iy)*morr
       enddo
       goto 1000
    endif
!------------------------------------ calculate as usual
100 continue
!----------------------------------
    mag=zero
    mat=zero
    map=zero
!    if(tpindep(2)) then
!       write(*,99)'MM d2G/dPdy: ',(pmi%curd%dgval(3,jy,1),jy=1,pmi%ncc)
!99     format(a,6(1pe11.3))
!    endif
! noofits=1 means phase is ideal, use only diagonal
    nocon=pmi%ncc
!    if(allocated(zib)) deallocate(zib)
    allocate(zib(nrel),stat=errall)
    if(nocon.gt.nrel) then
       big=.TRUE.
       if(allocated(maybesave)) deallocate(maybesave)
       allocate(maybesave(nrel+3,nocon),stat=errall)
    else
       big=.FALSE.
    endif
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 41: ',errall
       gx%bmperr=4370; goto 1000
    endif
    do ib=1,nrel
       mamu(ib)=zero
    enddo
    do iy=1,nocon
       zib=zero
       cig=zero; cit=zero; cip=zero
       do jy=1,nocon
          haha=pmi%invmat(jy,iy)
          do ib=1,nrel
             zib(ib)=zib(ib)+haha*pmi%dxmol(ib,jy)
          enddo
          cig=cig+haha*pmi%curd%dgval(1,jy,1)
! always calculate cit because cp debug ?? dgval(2,jy,1) is d2G/dTdy_j
          if(tpindep(1)) cit=cit+haha*pmi%curd%dgval(2,jy,1)
          if(tpindep(2)) cip=cip+haha*pmi%curd%dgval(3,jy,1)
       enddo
       morr=pmi%dxmol(ia,iy)
       do ib=1,nrel
          mamu(ib)=mamu(ib)+zib(ib)*morr
          if(big) maybesave(ib,iy)=zib(ib)
       enddo
       mag=mag+morr*cig
       if(tpindep(1)) mat=mat+morr*cit
       if(tpindep(2)) map=map+morr*cip
       if(big) then
          maybesave(nrel+1,iy)=cig
          maybesave(nrel+2,iy)=cit
          maybesave(nrel+3,iy)=cip
       endif
    enddo
    goto 1000
!
! Ignore the code for saveing below, use calc_dgdyterms1X
! To speed up calculations we save same values
! what must be saved is what should be multiplied with pmi%dxmol(ia,iy)
!    write(*,*)'Checking for saving ',noofits,10*pmi%iph+pmi%ics,nocon
    if(nocon.le.nrel) goto 1000
! ATTENTION this is not really used any longer, see calc_dgdyterms1X !!!
    if(nocon.gt.big1n) then
! save all data for this phase with a large number of constituents
       big1p=10*pmi%iph+pmi%ics
       big1n=nocon
       if(allocated(save1)) deallocate(save1)
       allocate(save1(nrel+3,nocon),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 42: ',errall
          gx%bmperr=4370; goto 1000
       endif
       do iy=1,nocon
          do ib=1,nrel+3
             save1(ib,iy)=maybesave(ib,iy)
          enddo
       enddo
!       write(*,*)'Saved 1 values for ',noofits,big1p,big1n
    elseif(nocon.gt.big2n) then
! save all data for this phases with a large number of constituents
       big2p=10*pmi%iph+pmi%ics
       big2n=nocon
       if(allocated(save2)) deallocate(save2)
       allocate(save2(nrel+3,nocon),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 43: ',errall
          gx%bmperr=4370; goto 1000
       endif
       do iy=1,nocon
          do ib=1,nrel+3
             save2(ib,iy)=maybesave(ib,iy)
          enddo
       enddo
!       write(*,*)'Saved 2 values for ',noofits,big2p,big2n
!    else
!       write(*,*)'dgdy not saved: ',noofits,10*pmi%iph+pmi%ics,nocon
    endif
1000 continue
    return
  end subroutine calc_dgdyterms1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_dgdyterms1X
!\begin{verbatim}
  subroutine calc_dgdyterms1X(nrel,ia,tpindep,mamu,mag,mat,map,pmi,noofits)
! THIS SUBROUTINE using allocatable arrays in phase_varres!!
! any change must also be made in subroutine calc_dyterms2 and calc_dgdytermsh
! calculate the terms in the deltay expression for amounts of component ia
!
! DM_A = \sum_B mu_B*MAMU(B) - MAG - MAT*dt - MAP*dp
!
! where MAMU=\sum_i dM_A/dy_i*\sum_j invmat(i,j)*dM_B/dy_j
!       c_iB=\sum_j invmat(i,j)*dM_B/dy_j etc etc
!
! it may not be very efficient but first get it right ....
! tpindep(1) is TRUE if T variable, tpindep(2) is TRUE if P are variable
!
! >>> THIS IS THE PRINCIPAL VERSION of calc_dgdyterms WITH SAVE
!
    implicit none
    integer ia,nrel,noofits
    logical tpindep(2)
    double precision, dimension(*) :: mamu
    double precision mag,mat,map
! no longer used ...
!    type(saveddgdy), pointer :: saved
! pmi is the phase data record for this phase
    type(meq_phase), pointer :: pmi
!\end{verbatim} %+
! THIS IS THE ONE CURRENTLY USED IN THE MINIMIZATIONS
! these are to be multiplied with mu(ib), nothing, deltaT, deltaP
    integer iy,jy,ib,nocon,errall
! initial values for saved results
!    integer :: sameit=0,big1p=0,big2p=0,big1n=0,big2n=0
    double precision cig,cit,cip,haha
    double precision morr
    double precision, allocatable, dimension(:) :: zib
!
!-----------
! \sum_i \sum_j e_ij*dM_A/dy_i dG/dy_j
!
    nocon=pmi%ncc
    mag=zero
    mat=zero
    map=zero
    do ib=1,nrel
       mamu(ib)=zero
    enddo
! the logic here is a bit complicated ...
! At the first iteration the pmi%curd%invsaved is deallocated
!    and the pmi%curd%invsavediter set to 0
!    but at first iteration no values are saved
!    so all terms calculated at each call
! At the second iteration a new pmi%curd%invsaved is allocated
!    and values are calculated and saved and pmi%curd%invsavediter set to 2
!    and these saved values are used in second and later calls
! At later iterations new values are calculated and saved in pmi%curd%invsaved
!    at first call if pmi%curd%invsavediter is less than current iteration
!    otherwise the saved values are used. 
! The first iteration could be improved slightly but I am not sure
!    pmi%curd%invsavediter can be trusted at the first iteration.
!---------------------------------
    if(noofits.le.1) then
! At the first iteration deallocate as we may have new conditions
       if(allocated(pmi%curd%invsaved)) deallocate(pmi%curd%invsaved)
       pmi%curd%invsavediter=0
!       write(*,17)'MM dgdycalc1X: ',noofits,pmi%iph,pmi%ics,nocon,ia,&
!            pmi%curd%invsavediter,allocated(pmi%curd%invsaved)
17     format(a,6i7,l2,4i4)
       goto 100
! UNFINISHED: VALGRIND indicates unititial variable ...
    elseif(pmi%curd%invsavediter.ne.noofits) then
! no values saved for this phase and iteration, recalcute
!                  123456789.12345
!       write(*,17)'MM new iter:   ',noofits,pmi%iph,pmi%ics,nocon,ia,&
!            pmi%curd%invsavediter,allocated(pmi%curd%invsaved)
       goto 100
    elseif(.not.allocated(pmi%curd%invsaved)) then
!       write(*,17)'MM Not allocated?',noofits,pmi%iph,pmi%ics,nocon,ia,&
!            pmi%curd%invsavediter,allocated(pmi%curd%invsaved)
       goto 100
    endif
! use save values for the phase
!               123456789.12345
!    write(*,17)'MM using save: ',noofits,pmi%iph,pmi%ics,nocon,ia,&
!         pmi%curd%invsavediter,allocated(pmi%curd%invsaved),&
!         size(pmi%curd%invsaved)
    if(allocated(pmi%curd%invsaved)) then
       do iy=1,nocon
          morr=pmi%dxmol(ia,iy)
          do ib=1,nrel
             mamu(ib)=mamu(ib)+pmi%curd%invsaved(ib,iy)*morr
          enddo
          mag=mag+pmi%curd%invsaved(nrel+1,iy)*morr
          if(tpindep(1)) mat=mat+pmi%curd%invsaved(nrel+2,iy)*morr
          if(tpindep(2)) map=map+pmi%curd%invsaved(nrel+3,iy)*morr
       enddo
       goto 1000
    else
       write(*,*)'MM ERROR: not allocated!',noofits,pmi%iph,pmi%ics,nocon,ia
       gx%bmperr=4399; goto 1000
    endif
!------------------------------------ calculate as usual and save at the end
100 continue
!----------------------------------
! next time for same iteration use saved values for this phase
!    sameit=noofits
! allocate the pmi%curd%invsaved at first iteration
    if(noofits.gt.1 .and. .not.allocated(pmi%curd%invsaved)) then
       allocate(pmi%curd%invsaved(nrel+3,nocon),stat=errall)
!       write(*,17)'MM allocate    ',noofits,pmi%iph,pmi%ics,nocon,ia,&
!            pmi%curd%invsavediter,allocated(pmi%curd%invsaved),&
!            nrel,(nrel+3)*nocon,size(pmi%curd%invsaved)
    endif
    allocate(zib(nrel),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 44: ',errall
       gx%bmperr=4370; goto 1000
    endif
!
!    write(*,17)'MM calculate:  ',noofits,pmi%iph,pmi%ics,nocon,ia,&
!         pmi%curd%invsavediter,allocated(pmi%curd%invsaved),&
!         nrel,(nrel+3)*nocon
    do iy=1,nocon
       zib=zero
       cig=zero; cit=zero; cip=zero
       do jy=1,nocon
          haha=pmi%invmat(jy,iy)
          do ib=1,nrel
             zib(ib)=zib(ib)+haha*pmi%dxmol(ib,jy)
          enddo
          cig=cig+haha*pmi%curd%dgval(1,jy,1)
! always calculate cit because cp debug ?? dgval(2,jy,1) is d2G/dTdy_j
          if(tpindep(1)) cit=cit+haha*pmi%curd%dgval(2,jy,1)
          if(tpindep(2)) cip=cip+haha*pmi%curd%dgval(3,jy,1)
       enddo
       morr=pmi%dxmol(ia,iy)
       do ib=1,nrel
          mamu(ib)=mamu(ib)+zib(ib)*morr
          if(noofits.gt.1) pmi%curd%invsaved(ib,iy)=zib(ib)
       enddo
       mag=mag+morr*cig
       if(tpindep(1)) mat=mat+morr*cit
       if(tpindep(2)) map=map+morr*cip
       if(noofits.gt.1) then
          pmi%curd%invsaved(nrel+1,iy)=cig
          pmi%curd%invsaved(nrel+2,iy)=cit
          pmi%curd%invsaved(nrel+3,iy)=cip
       endif
    enddo
    pmi%curd%invsavediter=noofits
!    write(*,17)'MM saveing:    ',noofits,pmi%iph,pmi%ics,nocon,ia,&
!         pmi%curd%invsavediter,allocated(pmi%curd%invsaved),&
!         size(pmi%curd%invsaved)
!
1000 continue
    return
  end subroutine calc_dgdyterms1X

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_dgdyterms2
!\begin{verbatim} %-
  subroutine calc_dgdyterms2(iy,nrel,mamu,mag,mat,map,pmi)
! Called only by meq_calc_phase_derivative
! for the contribution to G for a single phase
! it should be similar to calc_dgdyterms1
    implicit none
    integer iy,nrel
    double precision mag,mat,map,mamu(*)
    type(meq_phase), pointer :: pmi
!\end{verbatim} %+
! these are to be multiplied with mu(ib), nothing, deltaT, deltaP
! I am not sure if this is used ...
    integer jy,ib
    double precision sum,cig,cit,cip
!
!    write(*,*)'entering calc_dgdyterms2: ',iy,nrel,allocated(pmi%invmat)
    mag=zero
    do ib=1,nrel
       sum=zero
       do jy=1,pmi%ncc
          sum=sum+pmi%invmat(iy,jy)*pmi%dxmol(ib,jy)
       enddo
       mamu(ib)=sum
    enddo
!-----------
! \sum_i \sum_j e_ij*dM_A/dy_i dG/dy_j
    cig=zero
    cit=zero
    cip=zero
    do jy=1,pmi%ncc
       cig=cig+pmi%invmat(jy,iy)*pmi%curd%dgval(1,jy,1)
       cit=cit+pmi%invmat(jy,iy)*pmi%curd%dgval(2,jy,1)
       cip=cip+pmi%invmat(jy,iy)*pmi%curd%dgval(3,jy,1)
    enddo
    mag=cig
    mat=cit
    map=cip
1000 continue
    return
  end subroutine calc_dgdyterms2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_dgdytermsh
!\begin{verbatim} %-
  subroutine calc_dgdytermsh(nrel,ia,tpindep,hval,mamu,mag,mat,map,pmi,&
       curmux,noofits)
! This is a variant of dgdyterms1 including a term multiplied with each
! term (hval) in the summation over the comstituents as needed when calculating
! an equation for fix V or H.  If hval(i)=1.0 it should give the same
! results as dgdyterms1
!
! calculate the terms in the deltay expression for amounts of component ia
!
! DM_A = \sum_B mu_B*MAMU(B) - MAG - MAT*dt - MAP*dp
!
! where MAMU=\sum_i dM_A/dy_i*\sum_j invmat(i,j)*dM_B/dy_j
!       c_iB=\sum_j invmat(i,j)*dM_B/dy_j etc etc
!
! it may not be very efficient but first get it right ....
! tpindep(1) is TRUE if T variable, tpindep(2) is TRUE if P are variable
    implicit none
    integer ia,nrel,noofits
    logical tpindep(2)
    double precision, dimension(*) :: hval,mamu
    double precision mag,mat,map
    double precision curmux(*)
! pmi is the phase data record for this phase
    type(meq_phase), pointer :: pmi
!\end{verbatim} %+
! THIS IS MODIFIED FOR CONDITIONS ON H and related properties
! these are to be multiplied with mu(ib), nothing, deltaT, deltaP
! CHARRGE BALANCE TERM ADDED 150610!!!
    integer iy,jy,ib,neq
    double precision sum,cig,cit,cip,cib
    double precision morr,curmu(maxel),maq
!
!    write(*,9)'in calc_dgdytermsh: ',ia,0,0,pmi%chargebal
9   format(a,4i3,6(1pe12.4))
    mag=zero
    do ib=1,nrel
       sum=zero
       do iy=1,pmi%ncc
          cib=zero
          do jy=1,pmi%ncc
             cib=cib+pmi%invmat(iy,jy)*pmi%dxmol(ib,jy)
          enddo
          sum=sum+cib*hval(iy)
!          write(*,11)'termsh mu: ',ib,iy,0,hval(iy),sum
11        format(a,3i2,6(1pe12.4))
       enddo
       mamu(ib)=sum
    enddo
!-----------
!    if(noofits.eq.1) then
!       curmu=zero
!    else
    do iy=1,nrel
       curmu(iy)=curmux(iy)
    enddo
!    endif
!-----------
! \sum_i \sum_j e_ij*dM_A/dy_i dG/dy_j and other terms
! for phases with extrenal chargebalance we have one more row with index
! number of constituents+sublattices+1
    if(pmi%chargebal.eq.1) neq=pmi%ncc+size(pmi%curd%sites)+1
    maq=zero
    mag=zero
    mat=zero
    map=zero
    do iy=1,pmi%ncc
       cig=zero
       cit=zero
       cip=zero
       do jy=1,pmi%ncc
! I inversed order of iy, jy, does it still converge??
          cig=cig-pmi%invmat(jy,iy)*pmi%curd%dgval(1,jy,1)
!          write(*,11)'termsh g: ',ia,iy,jy,pmi%invmat(jy,iy),&
!               pmi%curd%dgval(1,jy,1),cig
! always calculate cit because cp debug!!
! hval(j)=dG/dy_j-Td2G/dTdy_j or something similar
          if(tpindep(1)) then
             cit=cit-pmi%invmat(jy,iy)*pmi%curd%dgval(2,jy,1)
!             write(*,11)'termsh t: ',ia,iy,jy,pmi%curd%dgval(2,jy,1),cit
          endif
          if(tpindep(2)) cip=cip-pmi%invmat(jy,iy)*pmi%curd%dgval(3,jy,1)
       enddo
!       morr=pmi%dxmol(ia,iy)
       morr=hval(iy)
       mag=mag+morr*cig
       mat=mat+morr*cit
       map=map+morr*cip
!       if(pmi%chargebal.eq.1) maq=maq+morr*pmi%invmat(neq,iy)
       if(pmi%chargebal.eq.1) maq=maq+morr*pmi%invmat(iy,neq)
    enddo
!    if(pmi%chargebal.eq.1) then
! Looking for the reason of bad convergence with enthalpy condition this
! was investigated but the correction is so small it is ignored.
! For phases with external charge balance there is one more term, e_ig*Q
! number of equations are constituents+sublattices+1
!       neq=pmi%ncc+size(pmi%curd%sites)+1
!       qscale=one
!       qscale=1.0D12
!       maq=maq*pmi%curd%netcharge*qscale
!       write(*,911)'eiq> ',pmi%curd%phtupx,pmi%chargebal,neq,pmi%ncc,&
!            pmi%curd%netcharge,mag,maq,(pmi%invmat(jy,neq),jy=1,neq)
!            pmi%curd%netcharge,mag,maq,(pmi%invmat(neq,jy),jy=1,neq)
911    format(a,4i4,3(1pe12.4),/6(1pe12.4))
! The contribution \sum_i e_iq*Q should be added (or subtracted) from mag
!       mag=mag+maq
!    endif
!    write(*,11)'termsh: ',ia,0,0,mag,mat,map,(mamu(jy),jy=1,nrel)
1000 continue
    return
  end subroutine calc_dgdytermsh

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_dgdytermshm
!\begin{verbatim} %-
  subroutine calc_dgdytermshm(nrel,ia,tpindep,hval,mamu,mag,mat,map,&
       mamu1,mag1,mat1,map1,pmi,curmux,noofits)
! This is a variant of dgdyterms1 including a term multiplied with each
! term (hval) in the summation over the comstituents as needed when calculating
! an equation for fix V or H.  If hval(i)=1.0 it should give the same
! results as dgdyterms1
!
! Probably only one of calc_dgdytermshm or calc_dgdytermsh is needed
! calculate the terms in the deltay expression for amounts of component ia
!
! DM_A = \sum_B mu_B*MAMU(B) - MAG - MAT*dt - MAP*dp
!
! where MAMU=\sum_i dM_A/dy_i*\sum_j invmat(i,j)*dM_B/dy_j
!       c_iB=\sum_j invmat(i,j)*dM_B/dy_j etc etc
!
! it may not be very efficient but first get it right ....
! tpindep(1) is TRUE if T variable, tpindep(2) is TRUE if P are variable
    implicit none
    integer ia,nrel,noofits
    logical tpindep(2)
    double precision, dimension(*) :: hval,mamu,mamu1
    double precision mag,mat,map,mag1,mat1,map1
    double precision curmux(*)
! pmi is the phase data record for this phase
    type(meq_phase), pointer :: pmi
!\end{verbatim}
! THIS IS MODIFIED FOR CONDITIONS ON H and related properties
! these are to be multiplied with mu(ib), nothing, deltaT, deltaP
! CHARGE BALANCE TERM ADDED 150610!!!
    integer iy,jy,ib,neq
    double precision sum,sum1,cig,cit,cip,cib
! these variables are probably redundant
    double precision morr,curmu(maxel),maq,maq1
!
!    write(*,9)'in calc_dgdytermsh: ',ia,nrel,pmi%ncc,pmi%chargebal
9   format(a,4i3,6(1pe12.4))
    mag=zero
    do ib=1,nrel
       sum=zero
       sum1=zero
       do iy=1,pmi%ncc
          cib=zero
          do jy=1,pmi%ncc
             cib=cib+pmi%invmat(iy,jy)*pmi%dxmol(ib,jy)
          enddo
          sum=sum+cib*hval(iy)
          sum1=sum1+cib*pmi%dxmol(ia,iy)
!          write(*,11)'termsh mu: ',ib,iy,0,hval(iy),pmi%dxmol(ia,iy),sum,sum1
11        format(a,3i2,6(1pe12.4))
       enddo
       mamu(ib)=sum
       mamu1(ib)=sum1
!       write(*,11)'dgdyhm: ',ia,ib,0,mamu(ib),mamu1(ib)
    enddo
!-----------
!    if(noofits.eq.1) then
!       curmu=zero
!    else
    do iy=1,nrel
       curmu(iy)=curmux(iy)
    enddo
!    endif
!-----------
! \sum_i \sum_j e_ij*dM_A/dy_i dG/dy_j and other terms
! for phases with extrenal chargebalance we have one more row with index
! number of constituents+sublattices+1
    if(pmi%chargebal.eq.1) neq=pmi%ncc+size(pmi%curd%sites)+1
    maq1=zero
    mag1=zero
    mat1=zero
    map1=zero
    maq=zero
    mag=zero
    mat=zero
    map=zero
    do iy=1,pmi%ncc
       cig=zero
       cit=zero
       cip=zero
       do jy=1,pmi%ncc
! I inversed order of iy, jy, does it still converge??
          cig=cig-pmi%invmat(jy,iy)*pmi%curd%dgval(1,jy,1)
!          write(*,11)'termsh g: ',ia,iy,jy,pmi%invmat(jy,iy),&
!               pmi%curd%dgval(1,jy,1),cig
! always calculate cit because cp debug!!
! hval(j)=dG/dy_j-Td2G/dTdy_j or something similar
          if(tpindep(1)) then
             cit=cit-pmi%invmat(jy,iy)*pmi%curd%dgval(2,jy,1)
!             write(*,11)'termsh t: ',ia,iy,jy,pmi%curd%dgval(2,jy,1),cit
          endif
          if(tpindep(2)) cip=cip-pmi%invmat(jy,iy)*pmi%curd%dgval(3,jy,1)
       enddo
       morr=pmi%dxmol(ia,iy)
       mag1=mag1+morr*cig
       mat1=mat1+morr*cit
       map1=map1+morr*cip
       if(pmi%chargebal.eq.1) maq1=maq1+morr*pmi%invmat(iy,neq)
!
       morr=hval(iy)
       mag=mag+morr*cig
       mat=mat+morr*cit
       map=map+morr*cip
!       if(pmi%chargebal.eq.1) maq=maq+morr*pmi%invmat(neq,iy)
       if(pmi%chargebal.eq.1) maq=maq+morr*pmi%invmat(iy,neq)
    enddo
1000 continue
    return
  end subroutine calc_dgdytermshm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine meq_evaluate_all_svfun
!\begin{verbatim}
 subroutine meq_evaluate_all_svfun(kou,ceq)
! evaluate (and list if kou>0) the values of all state variable functions
   implicit none
   integer kou
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! THIS SUBROUTINE MOVED FROM gtp3D
! if kou<0 no output
   character actual_arg(10)*24,star*2
   integer kf,nsvfun
   double precision val
   nsvfun=nosvf()
   if(kou.gt.0) write(kou,75)
75 format('No  Name ',12x,'Value')
   all: do kf=1,nsvfun
! functions with bit SVFVAL set will be ignored by meq_evaluate_svfun      
!      write(*,*)'MM meq_svfun: ',kf,svflista(kf)%name,&
!           btest(svflista(kf)%status,SVFVAL),ceq%svfunres(kf)
      star='  '
      if(btest(svflista(kf)%status,SVFVAL)) star='**'
      if(btest(svflista(kf)%status,SVFEXT)) star='<>'
!      if(btest(svflista(kf)%status,SVFVAL)) then
!         write(*,*)'MM only explit evaluation of: ',trim(svflista(kf)%name)
!         if(kou.gt.0) write(kou,77)kf,svflista(kf)%name,svflista(kf)%value,'*'
!         if(kou.gt.0) write(kou,77)kf,svflista(kf)%name,ceq%svfunres(kf),'*'
!         if(kou.gt.0) write(kou,78)kf,svflista(kf)%name,ceq%svfunres(kf),'**'
78       format(i3,1x,a,1x,1PE15.7,1x,a)
!78       format(i3,1x,a,1x,1PE15.8,a,' SVFVAL set')
!         cycle all
!      endif
! actual arguments needed if svflista(kf)%nactarg>0
!      write(*,*)'MM meq_svfun evaluate ',kf,svflista(kf)%name
      if(btest(svflista(kf)%status,SVFVAL)) then
! I am not really sure where the last calculated value is ???
! better to return zero than some arbitrary value
         val=zero
      else
         val=meq_evaluate_svfun(kf,actual_arg,0,ceq)
      endif
!      write(*,*)'MM meq_svfun evaluated: ',val
      if(gx%bmperr.ne.0) then
         if(kou.gt.0) then
            write(kou,76)kf,svflista(kf)%name,gx%bmperr
76          format(i3,1x,a,'  cannot be calculated due to error ',i5)
            if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
               write(kou,992)trim(bmperrmess(gx%bmperr))
992            format('Meaning: ',a/)
            endif
         endif
         gx%bmperr=0
      elseif(kou.gt.0) then
         write(kou,77)kf,svflista(kf)%name,val,star
77       format(i3,1x,a,1x,1PE15.7,' ',a)
      endif
! save the value in current equilibrium ... probably already done ...
      ceq%svfunres(kf)=val
   enddo all
1000 continue
   return
 end subroutine meq_evaluate_all_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine meq_get_state_varorfun_value
!\begin{verbatim}
 subroutine meq_get_state_varorfun_value(statevar,value,dummy,ceq)
! used in OCPLOT to extact value of state variable of symbol
! NOTE if a specific function is given only this function evaluated
   implicit none
   character statevar*(*),dummy*(*)
   double precision value
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character encoded*64,actual_arg(2)*16
   integer lrot,mode,olderr
!
!   write(*,*)'In meq_get_state_varofun: ',trim(statevar)
! if not derivative this will work
   call get_state_var_value(statevar,value,encoded,ceq)
!   write(*,*)'MM meq_get_state_varofun: ',gx%bmperr,value
   if(gx%bmperr.ne.0) then
! if error try using meq_evaluate_svfun
      olderr=gx%bmperr
      gx%bmperr=0
      encoded=statevar
      call capson(encoded)
!      call find_svfun(encoded,lrot,ceq)
      call find_svfun(encoded,lrot)
!      write(*,*)'In meq_get_state_varofun 2: ',&
!           trim(statevar),lrot,gx%bmperr,olderr
      if(gx%bmperr.ne.0) then
! if error here return previous error code
!         write(*,*)'In meq_get_state_varofun 3: ',gx%bmperr
         value=zero
         gx%bmperr=olderr; goto 1000
      else
         mode=1
         actual_arg=' '
! segmentation fault in this routine call from smp2B for a Cp value
! after shifting to a new maptop record (several STEP/MAP)
         value=meq_evaluate_svfun(lrot,actual_arg,mode,ceq)
      endif
   endif
! return calculated state variable symbol and always set special_circumstances=0
   dummy=encoded
! always reset to zero
   special_circumstances=0
1000 continue
   return
 end subroutine meq_get_state_varorfun_value

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable double precision function meq_evaluate_svfun
!\begin{verbatim}
 double precision function meq_evaluate_svfun(lrot,actual_arg,mode,ceq)
! evaluates all funtions as they may depend on each other
! actual_arg are names of phases, components or species as @Pi, @Ci and @Si
! needed in some deferred formal parameters  (NOT IMPLEMENTED YET)
! if mode=1 always evaluate, if mode=0 several options
   implicit none
   integer lrot,mode
   character actual_arg(*)*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! THIS SUBROUTINE MOVED FROM gtp3D
!   character encoded*60
   double precision argval(20)
   type(gtp_state_variable), target :: tsvr,tsvr2
   type(gtp_state_variable), pointer :: svr,svr2
   integer jv,jt,istv,ieq,nsvfun,ii
   double precision value
!
! modified here to handle symbols that can be used as conditions
!   write(*,*)'MM: --------- start meq_evaluate_svfun',trim(svflista(lrot)%name)
   value=zero
   argval=zero
   nsvfun=nosvf()
   ieq=0
   istv=0
! FIRST ALL SYMBOLS ARE EVALUATED HERE
!   write(*,*)'MM meq_evaluate_svfun 1 ',lrot,mode,svflista(lrot)%narg
! locate function
   if(lrot.le.0 .or. lrot.gt.nsvfun) then
      gx%bmperr=4140; goto 1000
   endif
! this seems OK
!   write(*,17)'meq_evaluate_svfun 2',lrot,trim(svflista(lrot)%name),&
!        svflista(lrot)%narg,&
!        btest(svflista(lrot)%status,SVFVAL),&
!        btest(svflista(lrot)%status,SVFEXT),&
!        btest(svflista(lrot)%status,SVCONST),&
!        btest(svflista(lrot)%status,SVFTPF),&
!        btest(svflista(lrot)%status,SVFDOT),&
!        btest(svflista(lrot)%status,SVNOAM)
17 format(a,i3,2x,a,i3,6l2)
   if(svflista(lrot)%narg.eq.0) goto 300
! get values of arguments
   jv=0
   jt=0
100 continue
      jt=jt+1
      istv=svflista(lrot)%formal_arguments(1,jt)
!      write(*,*)'MM meq_evaluate_svfun 3A',jt,istv
      if(istv.gt.-1000 .and. istv.lt.0) then
! istv values between -1000 and -1 are (negative) indices to functions
! istv values less than -1000 are parameter identication symbols
! if eqnoval nonzero it indicates from which equilibrium to get its value
         ieq=svflista(lrot)%eqnoval
!********************************************************************
! Note!! it should be evaluated!! Not implemented ... ???
!********************************************************************
         if(ieq.eq.0) then
            value=ceq%svfunres(-istv)
         else
            value=eqlista(ieq)%svfunres(-istv)
         endif
!         write(*,*)'in meq_evaluate_svfun 3X',ieq,istv,value
      else
! the need for 1:10 was a new bug discovered in GNU fortran 4.7 and later
         svr=>tsvr
! inside make_stvrec istv values less than -1000 are converted
         call make_stvrec(svr,svflista(lrot)%formal_arguments(1:10,jt))
         if(gx%bmperr.ne.0) goto 1000
         if(svflista(lrot)%formal_arguments(10,jt).eq.0) then
! get state variable or symbol value ... NOTE writing a TYPE wariable !!!
!            write(*,*)'MM meq_evaluate_svfun 3D: ',svr
            call state_variable_val(svr,value,ceq)
! error check at the end of if...
         else
! if special_circumstances=1 return with error code (supress a value to plot)
            if(special_circumstances.eq.1) then
!               write(*,*)'MM special_circumstances: ',special_circumstances
! error code 4373 means value supressed due to special_circumstances
               gx%bmperr=4373; goto 1000
            endif
! state variable derivative, the denominator is the next variable
            jt=jt+1
            svr2=>tsvr2
!            write(*,*)'MM meq_evaluate_svfun 3W: ',jt,svr
            call make_stvrec(svr2,svflista(lrot)%formal_arguments(1:10,jt))
!            write(*,77)'MM meq_eval: ',jt,&
!                 (svflista(lrot)%formal_arguments(ii,jt),ii=1,10)
77          format(a,i2,':',20i5)
! This routine need access to subroutines in the minimizer !!!
            call meq_state_var_dot_derivative(svr,svr2,value,ceq)
            if(gx%bmperr.ne.0) goto 1000
!            write(*,*)'MM back from meq_state_var_dot_derivative',value
         endif
      endif
      if(gx%bmperr.ne.0) goto 1000
      jv=jv+1
      argval(jv)=value
!      write(*,*)'in meq_evaluate_svfun 3B: ',jv,jt,argval(jv)
      if(jt.lt.svflista(lrot)%narg) goto 100
! all arguments evaluated (or no arguments needed)
300 continue
!      write(*,'(a,5i5,2l2)')'MM in meq_evaluate_svfun 300: ',lrot,mode,ieq,&
!           svflista(lrot)%eqnoval,istv,&
!           btest(svflista(lrot)%status,SVFVAL),&
!           btest(svflista(lrot)%status,SVFEXT),&
!           btest(svflista(lrot)%status,SVCONST)
   modeval: if(mode.eq.0 .and. btest(svflista(lrot)%status,SVFEXT)) then
! if mode=0 and SVFEXT=TRUE use value from equilibrium svflista(lrot)%eqnoval
!      write(*,*)'MM symbol mode=0 SVFEXT=TRUE: ',lrot,ieq,istv,argval(1)
      ieq=svflista(lrot)%eqnoval
      if(ceq%eqno.eq.ieq) then
         value=evalf(svflista(lrot)%linkpnode,argval)
!         write(*,*)'MM symbol calculated: ',lrot,ieq,istv,argval(1)
         if(pfnerr.ne.0) then
            write(*,*)'MM evaluate_svfun putfunerror ',pfnerr
            gx%bmperr=4141; goto 1000
         endif
! why store value in svfunres(-istv) ??? THIS MUST BE WRONG AND UNECESSARY
! we should store the value in the function restult for this equilibrium
         ceq%svfunres(lrot)=value
!         write(*,350)'MM evaluated here: ',ieq,lrot,value
      else
         value=eqlista(ieq)%svfunres(lrot)
         ceq%svfunres(lrot)=value
!         write(*,350)'MM value from equilbrium: ',ieq,lrot,value
      endif
   elseif(mode.eq.0 .and. btest(svflista(lrot)%status,SVFVAL)) then
! If mode=0 and SVFVAL set then return the stored value
! do not evaluate, just return the stored value in svfv(lrot) !!!
! copy to current ceq!!
         value=svflista(lrot)%svfv
         ceq%svfunres(lrot)=value
!         write(*,*)'MM in meq_evaluate_svfun 19:',lrot,ieq,value
!      write(*,350)'HMS evaluate svfun 2: ',0,lrot,value,svflista(lrot)%svfv
350   format(a,2i4,4(1pe13.5))
!      write(*,*)'MM in meq_evaluate_svfun  20: ',lrot,ieq,ceq%eqno,value
   elseif(btest(svflista(lrot)%status,SVCONST)) then
! symbol is a constant, just return value
      value=svflista(lrot)%linkpnode%value
      ceq%svfunres(lrot)=value
!      write(*,*)'MM symbol is a constant',lrot,value
   else
! if mode=1 always evaluate except if wrong eqilibrium!!
!      write(*,*)'in meq_evaluate_svfun 5',argval(1)
      if(svflista(lrot)%eqnoval.eq.0) then
         value=evalf(svflista(lrot)%linkpnode,argval)
         if(pfnerr.ne.0) then
            write(*,*)'evaluate_svfun putfunerror ',pfnerr
            gx%bmperr=4141; goto 1000
         endif
         ceq%svfunres(lrot)=value
      elseif(svflista(lrot)%eqnoval.eq.ceq%eqno) then
         value=evalf(svflista(lrot)%linkpnode,argval)
!         write(*,350)'HMS evaluate svfun 8: ',ieq,lrot,value,ceq%tpval(1)
         if(pfnerr.ne.0) then
            write(*,*)'evaluate_svfun putfunerror ',pfnerr
            gx%bmperr=4141; goto 1000
         endif
         ceq%svfunres(lrot)=value
      else
         ieq=svflista(lrot)%eqnoval
         value=eqlista(ieq)%svfunres(lrot)
         write(*,360)trim(svflista(lrot)%name),ieq,ceq%eqno
360      format('Attempt to evaluate symbol ',a,&
              ' for the wrong equilibrium:',2i5)
         ceq%svfunres(lrot)=value
      endif
   endif modeval
1000 continue
   meq_evaluate_svfun=value
   return
 end function meq_evaluate_svfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine initiate_meqrec
!\begin{verbatim}
  subroutine initiate_meqrec(svr,svar,meqrec,ceq)
! this is to setup data for a state var derivative calculation
! taken from the normal initialization of an equilibrium calculation
! it also solves a modified equil matrix once to get delta-amounts and mu
    TYPE(meq_setup), pointer :: meqrec
    TYPE(gtp_state_variable), pointer :: svr
    double precision, allocatable :: svar(:)
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    TYPE(meq_phase), pointer :: pmi
    integer iph,ics,kst,ie,mph,lokph,lokcs,nz1,tcol,pcol,dncol,converged
    integer ierr,nz2,jel,ztableph1,ztableph2,ztableph3,errall
    double precision, allocatable :: smat(:,:)
    double precision xxx
!
!    if(mmdebug.ne.0) write(*,*)'MM Entering initiate_meqrec 1'
! NOTE svar must be allocated!!
!    svar=zero
!    write(*,*)'MM Entering initiate_meqrec 2'
    if(btest(ceq%status,EQNOEQCAL)) then
! error if no sucessful equilibrium calculation or a failed one
!       write(*,*)'No equilibrium calculated, no derivatives'
!       allocate(svar(1)); svar(1)=zero
       gx%bmperr=4198; goto 1000
    elseif(btest(ceq%status,EQFAIL)) then
!       write(*,*)'Last equilibrium calculation failed, no derivatives'
!       allocate(svar(1)); svar(1)=zero
       gx%bmperr=4198; goto 1000
    elseif(btest(ceq%status,EQINCON)) then
! give warning if conditions have changed
       write(*,15)
15     format('Conditions changed since last equilibrium calc,',&
            ' values may be wrong.')
! EQNOACS is not used at present but means probably "no automatic comp.set"
!       allocate(svar(1)); svar(1)=zero
!       gx%bmperr=4198; goto 1000
    endif
! meqrec is a pointer to an allocated record!
!    allocate(meqrec)
! we must enter data into meqrec here, some set outside ...
!    meqrec%typesofcond=2
    meqrec%nrel=noel()
    meqrec%maxsph=noel()+2
    meqrec%nfixph=ceq%nfixph
    meqrec%nfixmu=ceq%nfixmu
! this returns total number of phases including composition sets
!    call sumofphcs(meqrec%nphase,ceq)
!    meqrec%nphase=totalphcs(ceq)
! if we are calculating a dot_derivative the number of phases in the dynamic
! memory may be different from that in the static memory!!
!    if(mmdotder.ne.0) write(*,*)'MM inititate_meqrec, mmdotder nonzero!'
    meqrec%nphase=nonsusphcs(ceq)
    if(gx%bmperr.ne.0) goto 1000
    allocate(meqrec%phr(meqrec%nphase),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 45: ',errall
       gx%bmperr=4370; goto 1000
    endif
! this means T and P are fixed (not independent)
    meqrec%tpindep=.FALSE.
    mph=0
    ztableph1=0
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 20:',noph()
! loop for all phases, we must set values of phase number etc
! meqrec%phr is later called "pmi"
    meqrec%nstph=0
    do iph=1,noph()
       do ics=1,noofcs(iph)
          call get_phase_compset(iph,ics,lokph,lokcs)
          if(mmdotder.ne.0) then
! test if composition set exists!! lokcs is always nonzero
! but if the sites are not allocated the composition set does not exist
             if(allocated(ceq%phase_varres(lokcs)%sites)) then
                kst=test_phase_status(iph,ics,xxx,ceq)
             else
!                write(*,*)'MM composition set does not exist!'
                kst=PHSUS
             endif
          else
             kst=test_phase_status(iph,ics,xxx,ceq)
          endif
!          meqrec%nv=meqrec%nv+1
          if(kst.ge.PHDORM) then
             mph=mph+1
             meqrec%phr(mph)%iph=iph
!             write(*,*)'phases: ',mph,iph
             meqrec%phr(mph)%ics=ics
! set number of constituents, DO NOT USE size(...curd%size(yfr)!!!
             meqrec%phr(mph)%ncc=noconst(iph,ics,ceq)
             meqrec%phr(mph)%phasestatus=kst
             meqrec%phr(mph)%ionliq=-1
             meqrec%phr(mph)%i2sly=0
             if(test_phase_status_bit(iph,PHIONLIQ)) meqrec%phr(mph)%ionliq=1
! set link to calculated values of G etc.
!             call get_phase_compset(iph,ics,lokph,lokcs)
             meqrec%phr(mph)%curd=>ceq%phase_varres(lokcs)
             if(kst.ge.PHENTSTAB) then
! this phase has the stable bit set
                ztableph1=ztableph1+1
                ztableph2=lokcs
                ztableph3=iph
                meqrec%phr(mph)%stable=1
                meqrec%nstph=meqrec%nstph+1
! store the index of the phase in phr, not the phase number 
                meqrec%stphl(meqrec%nstph)=mph
             else
! unstable phase
                meqrec%phr(mph)%stable=0
             endif
             meqrec%phr(mph)%idim=0
! valgrind found one case xdone was not initiated ....
             meqrec%phr(mph)%xdone=0
!          else
! nothing to do for suspended or hidden phase
          endif
       enddo
    enddo
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 30:'
    if(ztableph1.eq.1) then
! if there is a single stable phase, does it have fixed composition?
!       write(*,*)'MM a single stable phase',ztableph2
       if(size(ceq%phase_varres(ztableph2)%sites)-&
            size(ceq%phase_varres(ztableph2)%yfr).eq.0) then
!          write(*,*)'MM fixed composition: ',ztableph2
          xxx=-ceq%tpval(1)*ceq%phase_varres(ztableph2)%gval(4,1)
! The problem here was created somewhere else when the function for a phase
! to be optimized were changed, probably when trying to create
! already existing MAPNODE records.  That error not found !!
! Calculate G for this phase !!!
!          call calcg(ztableph3,1,2,ztableph2,ceq)
          allocate(svar(1),stat=errall)
          if(errall.ne.0) then
             write(*,*)'MM Allocation error 46: ',errall
             gx%bmperr=4370; goto 1000
          endif
! ATTENTION: THIS IS A VERY TEMPORARY FIX!!!!
! gval(4,1) is the CP of a stoichiometric compound
          svar(1)=-ceq%tpval(1)*ceq%phase_varres(ztableph2)%gval(4,1)
!          write(*,321)'MM fixed composition: ',ztableph2,&
!               lokcs,xxx,svar(1),ceq%tpval(1)
!321       format(a,2i5,4(1pe12.4))
          goto 1000
       endif
    endif
    meqrec%nphase=mph
! keep memory of adding/removing phases
!    write(*,*)'MM total number of phases: ',mph
! copy current values of ceq%complist%chempot(1) to ceq%cmuval, why??
    do ie=1,meqrec%nrel
       ceq%cmuval(ie)=ceq%complist(ie)%chempot(1)/ceq%rtn
    enddo
    meqrec%dormlink=0
! This can be done in PARALLEL for all phases
! nullify liquid pointer
    nullify(meqrec%pmiliq)
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 40:',meqrec%nphase
    do mph=1,meqrec%nphase
! loop to calculte and invert the phase matrices
       pmi=>meqrec%phr(mph)
!       write(*,*)'Inverting phase matrix ',mph
! This will calculate all G, dG/dZ1 and d2G/dZ1dZ2 and the inverted phase matrix
!       if(mmdebug.ne.0) write(*,*)'MM calling meq_onephase: ',mph
       call meq_onephase(meqrec,pmi,ceq)
       if(gx%bmperr.ne.0) then
          write(*,*)'Error calculating phase matrix'
          gx%bmperr=4199; goto 1000
       endif
    enddo
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 41:'
! now we will solve a modified phase matrix and calculate svar
! copy part of it from ceq%savesmat, copy also any fix mu and phase
! no problem to allocate as meqrec just allocated
    if(ceq%nfixmu.gt.0) then
       meqrec%nfixmu=ceq%nfixmu
       allocate(meqrec%mufixel(meqrec%nfixmu),stat=errall)
       do mph=1,ceq%nfixmu
          meqrec%mufixel(mph)=ceq%fixmu(mph)
       enddo
    endif
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 42:'
    if(ceq%nfixph.gt.0) then
       meqrec%nfixph=ceq%nfixph
       allocate(meqrec%fixph(2,meqrec%nfixph),stat=errall)
       do mph=1,ceq%nfixph
          meqrec%fixph(1,mph)=ceq%fixph(1,mph)
          meqrec%fixph(2,mph)=ceq%fixph(2,mph)
       enddo
    endif
! negative value of ceq%sysmatdim means no matrix saved
    nz1=abs(ceq%sysmatdim)+1
    allocate(smat(nz1,nz1+1),stat=errall)
    smat=zero
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 43:'
    allocate(svar(nz1),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 47: ',errall
       gx%bmperr=4370; goto 1000
    endif
! savesysmat not used, all equations calculated again
!    do mph=1,nz1-1
!       do ie=1,nz1-1
!          smat(mph,ie)=ceq%savesysmat(mph,ie)
!       enddo
!    enddo  
!    write(*,*)'Saved equil matrix',nz1
!    do jel=1,nz1
!       write(*,86)(smat(jel,nz2),nz2=1,nz1+1)
!    enddo
!86  format(6(1pe12.4))
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 44:'
    tcol=0
    pcol=0
! dncol is number of variable potentials (including T or P if variable)
    dncol=meqrec%nrel-meqrec%nfixmu
    converged=-1
    if(svr%statevarid.eq.1) then
       tcol=nz1
       meqrec%tpindep(1)=.TRUE.
    elseif(svr%statevarid.eq.2) then
       pcol=nz1
       meqrec%tpindep(2)=.TRUE.
    else
       write(*,*)'Derivatives with respect to T and P allowed only'
       gx%bmperr=4213; goto 1000
    endif
!-------------------------------------------------------------------
!    if(mmdebug.ne.0) write(*,854)'dncol mm: ',tcol,pcol,dncol,converged,nz1
854 format(a,10i5)
    call setup_equilmatrix(meqrec,meqrec%phr,nz1,smat,tcol,pcol,&
         dncol,converged,ceq)
! set all terms in the RHS to zero
    nz2=nz1+1
    do mph=1,nz1
       smat(mph,nz2)=zero
    enddo
!
! Add extra variable Delta-T for all stable phases: this is dG/dT
! This is redundant now??
    do mph=1,meqrec%nstph
       jel=meqrec%stphl(mph)
       smat(mph,nz1)=-meqrec%phr(jel)%curd%gval(2,1)
    enddo
! this is the line for Delta T or Delta P, all terms zero except last
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 48:'
    smat(nz1,nz1)=one
    smat(nz1,nz2)=one
! check matrix and rhs
!    write(*,*)'Equil matrix and solution in MM initiate_meqrec'
!    do jel=1,nz1
!       write(*,89)jel,(smat(jel,nz2),nz2=1,nz1+1)
!    enddo
89  format('MM qq: ',i2,6(1pe12.4))
! solve equil matrix 
!    if(mmdebug.ne.0) write(*,*)'MM initiate_meqrec 50:',nz1
    call lingld(nz1,nz1+1,smat,svar,nz1,ierr)
    if(ierr.ne.0) then
       write(*,*)'MM initiate_meqrec: error in lingld',ierr,nz1,ceq%eqno
!       do jel=1,nz1
!          write(*,89)jel,(smat(jel,nz2),nz2=1,nz1+1)
!       enddo
!       write(*,89)0,(svar(jel),jel=1,nz1)
       gx%bmperr=4214; goto 1000
!    else
    endif
!    write(*,89)0,(svar(jel),jel=1,nz1)
1000 continue
    return
  end subroutine initiate_meqrec ! allocated svar ??

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine meq_state_var_dot_derivative
!\begin{verbatim}
  subroutine meq_state_var_dot_derivative(svr1,svr2,value,ceq)
! calculates a state variable value, dot derivative, (in some cases)
! svr1 and svr2 identifies the state variables in (dstv1/dstv2)
! check that svr2 2 is a condition
! value is calculated value
! ceq is current equilibrium
! NOTE that when plotting after a STEP/MAP the number of phases in the
! dynamic memory (ceq) may be different that the static memory
! this is indicated by setting mmdotder nonzero
!
    implicit none
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(gtp_state_variable), pointer :: svr1,svr2
    double precision value
!\end{verbatim}
! variables needed to calculate phase inverse
    TYPE(meq_setup), allocatable, target :: meqrec1
    TYPE(meq_setup), pointer :: meqrec
 !   TYPE(meq_phase), pointer :: pmi
    TYPE(gtp_condition), pointer :: pcond
    integer iel,mph,jj,nterm,errall
    double precision xxx,sumam,summass,sumvol,s298
    double precision, allocatable :: svar(:)
    character dum*128,elsym*2
!
! this indicate that the number of phases in ceq may be different from
! the static number of phases
    mmdotder=1
    value=zero
!    write(*,*)'MM meq_state_var_dot_derivative 1'
!    if(svr2%statevarid.ne.1) then
! This if statement added trying to avoid spurious error (caused by -O2??)
!       write(dum,*)'In meq_state_var_value_derivative:',&
!            svr2%statevarid,ceq%tpval(1)
!    endif
! we must check if there is a condition on svr2
    pcond=>ceq%lastcondition
    if(.not.associated(pcond)) then
!       write(*,*)'There are no conditions at all!'
       gx%bmperr=4143; goto 1000
    endif
! all conditions have just one term at present
    nterm=1
    call get_condition(nterm,svr2,pcond)
    if(gx%bmperr.ne.0) then
       write(*,71)
71     format('To calculate a derivative the state variable after the dot',&
            ' must be a condition')
       goto 1000
    elseif(pcond%active.eq.1) then
! active=1 means not active
       write(*,71)
       goto 1000
    endif
! Currently only implemented H.T and HM.T
    if(.NOT.(svr2%statevarid.eq.1 .or. svr2%statevarid.eq.2)) then
       write(*,*)'Derivatives with respect to T and P only'
       gx%bmperr=4213; goto 1000
    endif
!------------
!    write(*,17)'minimzer: meq_state_var_value_derivative: ',&
!         svr1%statevarid,svr1%oldstv,svr1%argtyp,&
!         svr2%statevarid,svr2%oldstv,svr2%argtyp
!17 format(a,10i4)
! meqrec creates the data structure for the equilibrium data
! this routine also calculated Delta-amount of phases and delta-mu
    allocate(meqrec1,stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 48: ',errall
       gx%bmperr=4370; goto 1000
    endif
! problems with map7 ???
    meqrec1%status=0
    meqrec=>meqrec1
!    write(*,88)'MM calling initiate_meqrec',svr2%statevarid,ceq%eqno
88  format(a,2i4)
! indicate this is not an iteration by setting iteration number to -1
    meqrec%noofits=-1
! looking for segmentation fault from CP calculation when plotting
!    mmdebug=1
! initiate_meqrec will ignore nonexisting compostion sets
    call initiate_meqrec(svr2,svar,meqrec,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'MM back from initiate_meqrec in meq_state_var_dot_derivative'
!    mmdebug=0
!   addremloop
    iel=size(svar)
!    write(*,18)(svar(jj),jj=1,iel)
18  format('svar: ',6(1pe12.4)(6x,6e12.4))
    if(iel.eq.1) then
! iel=1 means a single stoichiometrc phase stable, svar(1) is CP/RT/T ??
! There can be a phase specification ...       
       if(svr1%statevarid.ge.6 .and. svr1%statevarid.lt.15 .and. &
            svr1%argtyp.eq.2) then
!          write(*,*)'MM Single stoichiometric phase stable',iel,svr1%argtyp
! nothing done?
          continue
       else
          if(svr1%norm.eq.2) then
! it is HW.T ....! where is mass of the phase?  Which phase?  Which element?
!             write(*,*)'MM single phase stable which phase?'
             if(noel().eq.1) then
                call get_element_data(1,elsym,dum(1:24),dum(25:48),summass,&
                     xxx,s298)
                if(gx%bmperr.ne.0) goto 1000
             endif
             value=svar(1)*ceq%rtn/summass
          else
             value=svar(1)*ceq%rtn
          endif
          goto 1000
       endif
    endif
!---------------
!100 continue
! if no phase specified loop over all stable phases
!    write(*,*)'We have initiad meqrec: ',svr1%statevarid
    if(svr1%statevarid.eq.3 .and. svr2%statevarid.eq.1) then
! This is MU(X).T
! it should simply be svar(svr1.%component) !!
!       write(*,*)'MM: MU(A).T: ',svr1%argtyp,svr1%component
       iel=svr1%component
       call meq_calc_phase_derivative(svr1,svr2,meqrec,mph,iel,&
            svar,jj,xxx,ceq)
       value=xxx*ceq%rtn
! there can be a suffix S ??
!       gx%bmperr=4215; goto 1000
! CCI already corrected
    elseif(svr1%statevarid.ge.6 .and. svr1%statevarid.lt.15) then
! This is derivatives of U, S, etc, H has svr1%statevarid=9, oldstv=40
! Partly DONE: implement H(phase).T and normalizing 
       iel=0
       jj=1
       sumam=zero
       summass=zero
       sumvol=zero
!       write(*,*)'MM Calculating H.T',svr1%argtyp,meqrec%nphase
! This "if" statement should be included in the loop below
       if(svr1%argtyp.eq.2) then
! if argtyp=2 then it is a value for a single phase
!          write(*,*)'MM svr1%argtyp 1: ',svr1%argtyp,svr1%phase,svr1%compset
          fph: do mph=1,meqrec%nphase
!             write(*,*)'fphloop: ',mph,meqrec%phr(mph)%iph,meqrec%phr(mph)%ics
! what is meqrec%iphl(mph) ???
!             if(meqrec%iphl(mph).eq.svr1%phase .and.&
!                  meqrec%icsl(mph).eq.svr1%compset) exit fph
             if(meqrec%phr(mph)%iph.eq.svr1%phase .and.&
                  meqrec%phr(mph)%ics.eq.svr1%compset) exit fph
          enddo fph
66        if(mph.gt.meqrec%nphase) then
             gx%bmperr=4050; goto 1000
          endif
! dummy statement to avoid some strange unknown error calculating Cp
          write(dum,*)'MM svr1%argtyp 2: ',svr1%argtyp,mph,iel
          call meq_calc_phase_derivative(svr1,svr2,meqrec,mph,iel,&
               svar,jj,xxx,ceq)
          if(gx%bmperr.ne.0) goto 1000
!          write(*,*)'MM HM(phase).T: ',xxx,meqrec%phr(mph)%curd%abnorm(1),&
!               meqrec%phr(mph)%curd%amfu
! normalized? svr1%norm>0
          if(svr1%norm.ne.0) then
             if(svr1%norm.eq.1) then
! xxx is HM for one formula unit.  if %norm=1 return HM.T 
                value=xxx/meqrec%phr(mph)%curd%abnorm(1)
             elseif(svr1%norm.eq.2) then
! xxx is HW per mass HW.T
                value=xxx/meqrec%phr(mph)%curd%abnorm(2)
             elseif(svr1%norm.eq.3) then
! norm=3 per volume HV.T
                write(*,*)'Normalizing per volume not implemented: ',svr1%norm
                gx%bmperr=4399; goto 1000
             elseif(svr1%norm.eq.4) then
! norm=4 per formula unit, HF.T
                value=xxx
             else
! no other normallizing
                write(*,*)'Unown normalizing: ',svr1%norm
                gx%bmperr=4399; goto 1000
             endif
          else
! not normalized value for a single phase, if amount zero return zero
             if(meqrec%phr(mph)%curd%amfu.eq.zero) then
! elseif %amfu=0 return H.T=0
                value=zero
             else
! else returm HM.T*NP(alpha) ???
             value=xxx*meqrec%phr(mph)%curd%amfu/meqrec%phr(mph)%curd%abnorm(1)
             endif
          endif
          goto 77
       endif
! sum over all stable phases
       do mph=1,meqrec%nphase
! ignore phases with zero amount
          if(meqrec%phr(mph)%curd%amfu.gt.zero) then
! the hope is that the phase amounts in svar are in the same order as
! in svar as ordered in meqrec%phr ...
! SEGMENTATION FAULT on LINUX with -O2 unless write statement at 69 is there
! It happends in macro step1 if you run all macros.  No error if just step1
! STRANGE !!!
             call meq_calc_phase_derivative(svr1,svr2,meqrec,mph,iel,&
                  svar,jj,xxx,ceq)
             if(gx%bmperr.ne.0) goto 1000
             sumam=sumam+&
                  meqrec%phr(mph)%curd%amfu*meqrec%phr(mph)%curd%abnorm(1)
             summass=summass+&
                  meqrec%phr(mph)%curd%amfu*meqrec%phr(mph)%curd%abnorm(2)
             jj=jj+1
! this dummy write statement is to avoid SEGMENTATION FAULT when -O2
! The segmentation fault persists also in oc5P if this write removed!!
!             write(dum,69)'MM der: ',mph,ceq%tpval(1),value,xxx,&
!                  meqrec%phr(mph)%curd%amfu
69           format(a,i3,6(1pe14.6))
          else
             xxx=zero
          endif
          value=value+xxx
       enddo
       if(svr1%norm.eq.1) then
! normallize with respect to number of moles of atoms
!          write(*,*)'MM sumam: ',value,sumam
          value=value/sumam
       elseif(svr1%norm.eq.2) then
! xxx is HW per mass HW.T
          value=value/summass
       elseif(svr1%norm.ne.0) then
! no other normallizing implemented
          write(*,*)'Illegal normalizing: ',svr1%norm
          gx%bmperr=4399; goto 1000
       endif
77     continue
    elseif(svr1%statevarid.eq.17) then
! This should be x(phase,element).T
!       write(*,*)'MM: X(PHASE,A).T not implemented',svr1%argtyp,svr1%phase,&
!            svr1%compset,svr1%component
       do mph=1,meqrec%nphase
          if(svr1%phase.eq.meqrec%phr(mph)%iph .and. &
               svr1%compset.eq.meqrec%phr(mph)%ics) then
             call meq_slope(mph,svr1,meqrec,value,ceq)
             write(*,*)'meq_slope: Not implemented x(phase,element).T'
             gx%bmperr=4215; goto 1000
          endif
       enddo
!       write(*,*)'No such phase'
!       gx%bmperr=4050
    else
       write(*,900)svr1%statevarid,svr1%argtyp,svr1%phase,svr1%compset,&
            svr1%component
900    format('MM: this dot derivative not implemented',6i5)
       gx%bmperr=4215; goto 1000
    endif
1000 continue
! meqrec1 deallocated automatically?
    if(allocated(meqrec1)) deallocate(meqrec1)
!    if(svr2%statevarid.ne.1) then
! This if statement added trying to avoid spurious error (caused by -O2??)
!       write(dum,*)'MM exit meq_state_var_value_derivative',value
!    endif
! reset mmdotder to zero
!    write(*,*)'MM exit meq_state_var_dotderivative'
    mmdotder=0
    return
  end subroutine meq_state_var_dot_derivative

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
    
!\addtotable subroutine meq_calc_phase_derivative
!\begin{verbatim}
  subroutine meq_calc_phase_derivative(svr1,svr2,meqrec,iph,iel,&
       svar,jj,value,ceq)
! Calculate contribution for one phase, one or all elements
! svr1 and svr2 identifies the state variables in (dstv1/dstv2)
! value is calculated value returned
! iph and iel indicate possible phase or element
! svar is solution to equil matrix, potentials and phase amounts
! jj is an attempt to index phases in svar, starting with 1
! ceq is current equilibrium
!
! THIS IS UNFINISHED can only handle H.T
!
    implicit none
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(gtp_state_variable), target :: svr1,svr2
    TYPE(meq_setup), pointer :: meqrec
    integer iph,iel,jj
    double precision value,svar(*)
!\end{verbatim}
! variables needed to calculate phase inverse
    TYPE(meq_phase), pointer :: pmi
    integer jy,jel,jz,phncc,errall
    double precision x1,x2,x3
    double precision mag,mat,map,dpham,musum,dy,hconfig
    double precision, allocatable :: mamu(:)
!
! THE MASTER VERSION OF THIS TABLE in GTP3C.F90
! symb cmix(2) indices                   statevarid Property
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
! DG      9x    "                            15 ?   Driving force
! Q       19x   "                            14 ?   Internal stability
! N       11x  (component/phase#set,component) 16  moles of components
! X       111   "                            17     mole fraction of components
! B       12x   "                            18     mass of components
! W       122   "                            19     mass fraction of components
! Y       13    phase#set,constituent#subl   20     constituent fraction
! statevarid=1 is T, 2 is P, 3 is MU, 4 is AC, 5 is LNAC
!------------------------------------------------------------
!    write(*,*)'MM meq_calc_phase_derivative',iph,iel
    allocate(mamu(meqrec%nrel),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 49: ',errall
       gx%bmperr=4370; goto 1000
    endif
    pmi=>meqrec%phr(iph)
    value=zero
! CCI
    hconfig=zero
    if(iel.lt.0) then
! sum for all elements
       write(*,*)'sum over elements not implemented'
       gx%bmperr=4216
    elseif(iel.eq.0) then
! independent of element, return for phase
       musum=zero
! pmi%ncc here is not set correctly ... WHEN?
!       phncc=size(pmi%curd%yfr)
       phncc=pmi%ncc
! PROBLEM 181208: step5 macro: size(pmi%curd%yfr) is 1000 (default)
! but pmi%ncc is 8 (total number of constituents)
! I do not remember why this was changed
!       write(*,*)'MM derivative: ',iph,pmi%ncc,phncc
!       do jy=1,pmi%ncc
       do jy=1,phncc
! The loop to handle the contribution from fractions in each phase dZ/dyi
          dy=zero
! special if just a single element ...
          if(allocated(pmi%invmat)) then
             call calc_dgdyterms2(jy,meqrec%nrel,mamu,mag,mat,map,pmi)
             if(gx%bmperr.ne.0) goto 1000
          else
! we have a stoichiometric phase with a single component ??
!             write(*,*)'MM No inverted phase matrix allocated',jy
             mamu=zero; mag=zero
!             gx%bmperr=4399; goto 1000
          endif
          jz=1
          if(meqrec%nfixmu.gt.0) then
! if there are fixed potentials such elements should be ignored here
! as there is no value in svar (value is zero as fixed)
           write(*,*)'Dot derivatives with potential condition not implemented'
             goto 1000
          endif
! sum the contribution for the potentials
          do jel=1,meqrec%nrel
             jz=jz+1
             dy=dy+mamu(jel)*svar(jel)
!             write(*,666)'dy: ',mamu(jel),svar(1),dy
          enddo
          dy=dy-mat
!          write(*,666)'dy: ',mat,dy
! here we check which state variable we take derivative of, H is 9
!          write(*,*)'MM svr1: ',svr1%statevarid,svr1%norm
          select case(svr1%statevarid)
          case default
! state variables 1..5 are potentials, 14-15 not possible to derivate
             write(*,*)'Illegal state variable id:',svr1%statevarid
             gx%bmperr=4188; goto 1000
          case(6) !U = G + TS - PV = G - T G.T - P G.P
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(7) !S = -G.T
             hconfig=-pmi%curd%dgval(2,jy,1)
          case(8) !V = G.P
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(9) !H = G + TS = G - T G.T
! this gives contribution also when plotting H(liq).T and HM(liq).T in step1
! but it is identical to Thermo-Calc .... thus correct
             hconfig=pmi%curd%dgval(1,jy,1)-ceq%tpval(1)*pmi%curd%dgval(2,jy,1)
          case(10) !A = G - PV = G - P G.P
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(11) !G itself, dG/dy
             hconfig=pmi%curd%dgval(1,jy,1)
          case(12) !NP phase amount
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(13) !BP phase mass
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(16) !N
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(17) !X
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(18) !B
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(19) !W
             write(*,*)'Not implemented yet: ',svr1%statevarid
          case(20) !Y
             write(*,*)'Not implemented yet: ',svr1%statevarid
          end select
!          if(svr1%statevarid.eq.9) then
!             hconfig=pmi%curd%dgval(1,jy,1)-ceq%tpval(1)*pmi%curd%dgval(2,jy,1)
!          endif
          musum=musum+hconfig*dy
!          write(*,*)'musum: ',musum,dy
       enddo
       x1=zero; x2=zero
!       write(*,765)'x3= ',ceq%rtn,pmi%curd%amfu,musum
       if(svr1%norm.eq.1 .and. svr1%argtyp.eq.2) then
! for HM(phase).T the change in of phase amount should be ignored
          dpham=zero
          x3=musum*ceq%rtn
       else
! extract the change in phase amount (for stable phases!!!)
! we have to take care of fixed chemical potentials, the number of
! elements+1-(#fixed mu) should be the index of dpham,
! the change in phase amount
! The way of indexing with jj is dangerous ...
          dpham=svar(meqrec%nrel+jj)
! for current amount
          x3=musum*ceq%rtn*pmi%curd%amfu
       endif
!       write(*,665)'dpham: ',meqrec%nrel,jj,svar(meqrec%nrel+jj-1),&
!            svar(meqrec%nrel+jj)
665    format(a,2i3,6(1pe14.6))
! here we again select the state variable we take derivative of, H is 9
!       write(*,*)'MM svr2: ',svr1%statevarid,svr1%norm
       select case(svr1%statevarid)
       case default
! state variables 1..5 are potentials, 14-15 not possible to derivate
          write(*,*)'Illegal state variable id:',svr1%statevarid
          gx%bmperr=4188; goto 1000
       case(6) !U = G + TS - PV
          write(*,*)'Not implemented yet: ',svr1%statevarid
          gx%bmperr=4215
       case(7) !S = -dG/dT
          x1=-ceq%rtn*dpham*pmi%curd%gval(2,1)
          x2=-ceq%rtn*pmi%curd%amfu*ceq%tpval(1)*pmi%curd%gval(4,1)
          write(*,*)'Not implemented yet: ',svr1%statevarid
          gx%bmperr=4215
       case(8) !V = dG/dP
          write(*,*)'Not implemented yet: ',svr1%statevarid
          gx%bmperr=4215
       case(9) !H = G + TS = G - T G.T
! x1 is change in phase amount times H.  Skip this if svr1%norm.eq.1 
          x1=-ceq%rtn*dpham*&
                  (pmi%curd%gval(1,1)-ceq%tpval(1)*pmi%curd%gval(2,1))
!          write(*,666)'x1: ',ceq%rtn,dpham,pmi%curd%gval(1,1),&
!               ceq%tpval(1)*pmi%curd%gval(2,1),x1
! x2 is phase_amount * dH/dT = .. -T*d2G/dT2 = -T
! CCI changed order of tests, does not work for step1
          if(dpham.ne.zero) then
! there is a change in phase amounts
             x2=-ceq%rtn*pmi%curd%amfu*ceq%tpval(1)*pmi%curd%gval(4,1)
          elseif(svr1%norm.eq.1) then
!xCCI          if(svr1%norm.eq.1) then
! compared with Thermo-Calc this seems correct, it is just HM(phase).T
             x2=-ceq%rtn*ceq%tpval(1)*pmi%curd%gval(4,1)
!             write(*,444)'Phase: ',iph,x1,x2,x3
!444          format(a,i3,3(1pe14.6))
!xCCI          else
!xCCI             x2=-ceq%rtn*pmi%curd%amfu*ceq%tpval(1)*pmi%curd%gval(4,1)
          else
! This is H.T or H(phase).T, should be (amount of phase)*HM.T
! when there is no change of amount of phase
             x2=-ceq%rtn*pmi%curd%amfu*ceq%tpval(1)*pmi%curd%gval(4,1)
          endif
! CCI end of correction
       case(10) !A = G - PV
          write(*,*)'Not implemeneted yet: ',svr1%statevarid
       case(11) !G itself
          x1=-ceq%rtn*dpham*pmi%curd%gval(1,1)
          x2=ceq%rtn*pmi%curd%amfu*pmi%curd%gval(2,1)
!          write(*,*)'G.T: ',x1,x2
       case(12) !NP phase amount
          write(*,*)'Not implemeneted yet: ',svr1%statevarid
       case(13) !BP phase mass
          write(*,*)'Not implemeneted yet: ',svr1%statevarid
       case(16) !N moles
          write(*,*)'Not implemeneted yet: ',svr1%statevarid
       case(17) !X mole fraction
          write(*,*)'Not implemeneted yet: ',svr1%statevarid
       case(18) !B mass
          write(*,*)'Not implemeneted yet: ',svr1%statevarid
       case(19) !W mass fraction
          write(*,*)'Not implemeneted yet: ',svr1%statevarid
       case(20) !Y constituent fraction
          write(*,*)'Not implemeneted yet: ',svr1%statevarid
       end select
!       if(svr1%statevarid.eq.9) then
! x1 is change in phase amount times H
!         x1=-ceq%rtn*dpham*(pmi%curd%gval(1,1)-ceq%tpval(1)*pmi%curd%gval(2,1))
! x2 is phase_amount * dH/dT = .. -T*d2G/dT2 = -T
!         x2=-ceq%rtn*pmi%curd%amfu*ceq%tpval(1)*pmi%curd%gval(4,1)
!       endif
! x3 is phase amount times change in configuration
!       x3=ceq%rtn*pmi%curd%amfu*musum
! only derivativs wrt T are allowed!!
!       write(*,666)'CP= ',svr1%norm,x1,x2,x3,x1+x2+x3,dpham,pmi%curd%amfu
666    format(a,i3,6(1pe12.4))
! just to show the error
!       value=x2
       value=x1+x2+x3
    else
! the derivative of the chemical potential of iel wrt T
       value=svar(iel)
!       write(*,*)'Chemical potential: ',iel,value
!       gx%bmperr=4215
    endif
!
1000 continue
    return
  end subroutine meq_calc_phase_derivative

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine meq_slope
!\begin{verbatim}
  subroutine meq_slope(mph,svr,meqrec,value,ceq)
! Test subroutine for x(phase,A).T   UNFINISHED
    TYPE(meq_setup) :: meqrec
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(gtp_state_variable) :: svr
    double precision value
    integer mph
!\end{verbatim}
!    TYPE(meq_phase), pointer :: pmi
    integer nsl,nkl(10),knr(maxconst)
    double precision yarr(maxconst),sites(10),qq(5)
!    
    call get_phase_data(svr%phase,svr%compset,nsl,nkl,knr,yarr,sites,qq,ceq)
    if(gx%bmperr.ne.0) goto 1000
! UNFINISHED
!
1000 continue
    return
  end subroutine meq_slope

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

  subroutine calfun(m,n,x,f,info,niter)
! This is called by the LMDIF1 routines and calls an OC subroutine
! as I had problems using EXTERNAL
! M is number of errors
! N is number of variables
! NOTE order of X and F switched in CALFUN and ASSESSMENT_CALFUN !!!
    integer m,n,info,i,niter
    double precision f(m),x(n),sum
!    write(*,*)'MM enter calfun',info,niter,m,n
    if(info.eq.0) then
       sum=zero
       do i=1,m
          sum=sum+f(i)**2
       enddo
       if(niter.lt.0) then
! this marks end of optimization output of the individual errors
          write(*,15)-niter,sum
15        format(/'Final results after ',i3,&
               ' iteration, the sum of squares',1pe14.6)
          write(*,16)x
16        format('Scaled param: ',4(1pe14.6)/5(1pe14.6))
          write(*,17)f
17        format('Errors: '/6(1pe13.5))
          write(*,*)
       elseif(niter.ge.0) then
          write(*,18)niter,sum
18        format(/'After ',i4,' iterations the sum of squares',1pe14.6)
          write(*,19)x
19        format('Scaled param:   ',4(1pe16.8)/5(1pe16.8))
       endif
    else
! This routine is in the matsmin.F90 file
! it returns the calculated value of the property to fit
! This call removed and the whole subroutine is probably redundant
       call assessment_calfun(m,n,f,x)
    endif
    return
  end subroutine calfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine assessment_calfun
!\begin{verbatim}
  subroutine assessment_calfun(nexp,nvcoeff,errs,xyz)
! nexp is number of experiments, nvcoeff number of coefficients
! errs is the differences between experiments and value calculated by model
! returned by this subroutine
! xyz are the scaled current model parameter values
    implicit none
    integer nexp,nvcoeff
    double precision errs(*),XYZ(*)
!    type(gtp_assessmenthead), pointer :: ash
!\end{verbatim}
! firstash is the data structure for assessment head (globally declared) 
    integer i1,i2,iexp,symsym,mode,jj,savix,next
    double precision xxx,yyy,zzz
    type(gtp_equilibrium_data), pointer :: equil
    type(gtp_condition), pointer :: experiment
    type(gtp_state_variable), pointer :: svrrec
    character text*24
    double precision xa(100)
!
!    write(*,*)'MM in assessment_calfun',nexp,nvcoeff
!    if(allocated(calcexp)) write(*,*)'Calculating Jacobian'
! 1. copy values of X to the TP coefficinets, loop through all
    i2=1
    do i1=0,size(firstash%coeffstate)-1
!       write(*,*)'MM2 Testing value of firstash%coeffstate',i1
       if(firstash%coeffstate(i1).ge.10) then
!          write(*,*)'MM3 coefficient ',i1,i2,xyz(i2)
! Attempt to handle that I divide coef with scaling factor ...
          zzz=xyz(i2)*firstash%coeffscale(i1)
          xxx=xyz(i2)*firstash%coeffscale(i1)
          call get_value_of_constant_index(firstash%coeffindex(i1),zzz)
!          write(*,16)i2,i1,xyz(i2),firstash%coeffscale(i1),xxx,zzz
16        format('MM4 Opt coeff ',2i4,4(1pe12.4))
          savix=i1
          call change_optcoeff(firstash%coeffindex(i1),xxx)
          if(gx%bmperr.ne.0) goto 1000
          xa(i2)=xxx
          i2=i2+1
!       else
!          write(*,*)'MM5 coefficient not variable',i1
       endif
    enddo
! 2. calculate all differences, skipping equilibria with weight zero
! the array firstash%eqlista contain pointers to equilibria with experiments
700 continue
    if(.not.allocated(firstash%eqlista)) then
       write(kou,*)' *** Warning: no experimental data!'
       do i1=1,nexp
          errs(i1)=zero
       enddo
       goto 1000
!    else
!       write(*,*)'MM6 First equilibrium number: ',firstash%firstexpeq
!       write(*,17)size(firstash%eqlista),firstash%firstexpeq
17     format('MM Number of equilibra with experiments: ',i5,', first is ',i3)
!       do i1=1,size(firstash%eqlista)
!          write(*,21)i1,firstash%eqlista(i1)%p1%eqname
!21        format('MM Equilibrium number ',i3,' and name: ',a)
!       enddo
    endif
! Seach for any symbol that should be calculated at a particulat equilibrium
! For example a reference enthalpy.  This equilibrium must be calculated
! before any parallel calculation of the others
    next=-1
    do while(next.ne.0)
       call find_symbol_with_equilno(next,i1)
!       write(*,*)' ******* checking for equilibrium to be calculated first'
       if(i1.gt.0) then
          if(firstash%eqlista(i1)%p1%weight.gt.zero) then
             equil=>firstash%eqlista(i1)%p1
!             write(*,*)' ******* equilibrium to be calculated first: ',i1
! Force recalculation of all TP functions and parameters by changing saved T
! This does not change the value of T used for the equilibrium
             equil%eq_tpres%tpused(1)=equil%tpval(1)+one
! calculate the equilibria without grid minimizer
             mode=-1
             call calceq3(mode,.FALSE.,equil)
             if(gx%bmperr.ne.0) then
                write(kou,33)gx%bmperr,equil%eqno,trim(equil%eqname)
                gx%bmperr=0
             endif
             text=' '
! evaluate symbol "next" (which is current!!) with force 
             xxx=evaluate_svfun_old(next,text,1,equil)
             if(gx%bmperr.ne.0) then
                gx%bmperr=0
                xxx=meq_evaluate_svfun(next,text,1,equil)
             endif
! we do not need the value here, it is stored at the symbol
!             write(*,*)'MM Symbol at equil: ',next,i1,gx%bmperr,xxx
          endif
       endif
    enddo
! loop through all equilibria with experiments
! each can be calculated in parallel
    iexp=0
    if(gx%bmperr.ne.0) then
       write(*,*)'In assessment_calfun: resting error code: ',gx%bmperr
       gx%bmperr=0
    endif
    eqloop: do i1=1,size(firstash%eqlista)
       if(firstash%eqlista(i1)%p1%weight.eq.zero) then
!          write(*,29)i1,firstash%eqlista(i1)%p1%eqname
29        format('MM Skipping equilibrium number ',i3,' and name: ',a)
          cycle eqloop
       endif
!       write(*,30)i1,trim(firstash%eqlista(i1)%p1%eqname)
30     format('MM Assessment_calfun equilibrium number ',i3,' and name: ',a)
       equil=>firstash%eqlista(i1)%p1
! Force recalculation of all TP functions and parameters by changing saved T
       equil%eq_tpres%tpused(1)=equil%tpval(1)+one
! calculate the equilibria without grid minimizer
!       write(*,*)'MM calculating equil: ',equil%eqno
! mode=-1 do not use gridmin and check after ...
       mode=-1
       call calceq3(mode,.FALSE.,equil)
       if(gx%bmperr.ne.0) then
          write(kou,33)gx%bmperr,equil%eqno,trim(equil%eqname)
33        format(' *** Error ',i5,' calculating equilibrium no: ',i5,&
               ' with name ',a)
          gx%bmperr=0
          cycle
!       else
!          write(*,*)'Equilibrium calculated for ',equil%eqname
       endif
! loop through all experiments, pointer set to first
       if(.not.associated(equil%lastexperiment)) then
!          write(*,*)'No experiments for equilibrium ',equil%eqno
          cycle eqloop
       endif
       experiment=>equil%lastexperiment%next
! current value of the experiment
500    continue
          iexp=iexp+1
!          write(*,*)'MM Setting pointer to experiment ',&
!               allocated(experiment%statvar),iexp
          nostv: if(.not.allocated(experiment%statvar)) then
             symsym=experiment%statev
             text=' '
! WE MUST EVALUATE ALL SYMBOLS!!!
             call meq_evaluate_all_svfun(-1,equil)
!             write(*,*)'MM symsym: ',symsym
             xxx=evaluate_svfun_old(symsym,text,1,equil)
             if(gx%bmperr.ne.0) then
                gx%bmperr=0
!                write(*,*)'MM using meq_evaluate_svfun',gx%bmperr
                xxx=meq_evaluate_svfun(symsym,text,1,equil)
             endif
!             write(*,*)'MM value: ',iexp,xxx
          else
             svrrec=>experiment%statvar(1)
!             write(*,*)'MM exp: ',svrrec%statevarid,svrrec%argtyp
! svrrec%statevarid = 0 means symbol ...
! this can handle state variable symbols also !!??
             call state_variable_val(svrrec,xxx,equil)
          endif nostv
          if(gx%bmperr.ne.0) then
             write(kou,*)' *** Error calculating experiment ',&
                  equil%eqno,': ',trim(equil%eqname),symsym,gx%bmperr
             gx%bmperr=0
             errs(iexp)=zero
             goto 590
          endif
          if(experiment%symlink2.gt.0) then
! added check if uncertainity is a symbol
!            xxx=evaluate_svfun_old(istv,'  ',mode,ceq)
!             xxx=evaluate_svfun_old(symsym,text,1,equil)
             experiment%uncertainty=&
                  evaluate_svfun_old(experiment%symlink2,' ',1,equil)
          endif
!          write(*,510)'MM errs',iexp,experiment%prescribed,xxx,&
!               experiment%uncertainty,equil%weight
510       format(a,i4,6(1pe12.4))
          if(allocated(calcexp)) then
! this is to enable calculating RSD at the end of an assessment
! normally calcexp is not allocated!!
             calcexp(iexp)=xxx
!             write(*,555)'Jacobian: ',iexp,(xa(jj),jj=1,i2-1),xxx
555          format(a,i3,6(1pe12.4))
          endif
          if(experiment%experimenttype.eq.0) then
! take the difference between prescribed value
             errs(iexp)=(experiment%prescribed-xxx)*equil%weight/&
                  experiment%uncertainty
!             write(*,*)'MM least.sq: ',iexp,f(iexp)
          elseif(experiment%experimenttype.eq.100) then
! relative error
             yyy=1.0D-2*experiment%uncertainty*experiment%prescribed
             errs(iexp)=(experiment%prescribed-xxx)*equil%weight/yyy
          elseif(experiment%experimenttype.eq.-1) then
! less than, uncertainty is penalty function factor
             if(xxx.gt.experiment%prescribed) then
                errs(iexp)=(xxx-experiment%prescribed)*equil%weight/&
                     experiment%uncertainty
             else
                errs(iexp)=zero
             endif
          elseif(experiment%experimenttype.eq.1) then
! larger than, uncertainty is penalty function factor
             if(xxx.lt.experiment%prescribed) then
                errs(iexp)=(xxx-experiment%prescribed)*equil%weight/&
                     experiment%uncertainty
             else
                errs(iexp)=zero
             endif
          endif
590       if(.not.associated(experiment,equil%lastexperiment)) then
! if more experiments jump back to 500
             experiment=>experiment%next
             goto 500
          endif
! done all experiments for this equilibrium
    enddo eqloop
!    write(*,*)'MM experiments: ',iexp,nexp
! We have to restore the last value of the last coefficient
    if(allocated(calcexp)) then
!       write(*,*)'MM restore savix: ',savix,zzz
       call change_optcoeff(firstash%coeffindex(savix),zzz)
    endif
1000 continue
!    write(*,*)'Exit assessment_calfun'
    return
  end subroutine assessment_calfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine meq_list_experiments
!\begin{verbatim}
 subroutine meq_list_experiments(lut,ceq)
! list all experiments into text, special to handle derivatives ...
   implicit none
   integer lut
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer seqz,ip
   character text*72
   seqz=0
100 continue
      seqz=seqz+1
      ip=1
      text=' '
      call meq_get_one_experiment(ip,text,seqz,ceq)
!      write(*,*)'MM Back from get_one'
      if(gx%bmperr.ne.0) then
! error code for no more experiments or inactive experiment
!         write(*,*)'MM error line 3117: ',gx%bmperr,seqz,text(1:ip)
! speciel error code meaning experiment is not active
         if(gx%bmperr.eq.7654) then
            gx%bmperr=0; goto 100
         endif
         gx%bmperr=0; goto 1000
      else
         write(lut,120)seqz,text(1:ip)
120      format('Experiment ',i2,2x,a)
      endif
      goto 100
!------------
1000 continue
   gx%bmperr=0
   return
 end subroutine meq_list_experiments

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine meq_get_one_experiment
!\begin{verbatim} %-
 subroutine meq_get_one_experiment(ip,text,seqz,ceq)
! list the experiment with the index seqz into text
! It lists also experiments that are not active ??
! UNFINISHED current value should be appended
   implicit none
   integer ip,seqz
   character text*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer iterm,symsym,mode
   TYPE(gtp_condition), pointer :: last,current
   type(gtp_state_variable), pointer :: svrrec
   double precision xxx
   character actual_arg*16
!
   if(ip.le.0) ip=1
   text(ip:)=' '
   if(.not.associated(ceq%lastexperiment)) then
      write(*,*)'MM No experiments'
      gx%bmperr=4249; goto 1000
   endif
   last=>ceq%lastexperiment
   current=>last
!   write(*,*)'MM index of last experiment: ',current%seqz
70 continue
!   write(*,*)'MM experiment number: ',seqz,current%seqz
   if(current%seqz.eq.seqz) goto 100
   current=>current%next
   if(.not.associated(current,last)) goto 70
! no experiment with this index found or it is inactivated
   gx%bmperr=4131; goto 1000
!
100 continue
   if(current%active.eq.1) then
!      write(*,*)'MM Experiment not active '
      gx%bmperr=4218; goto 1000
   endif
   iterm=1
150 continue
!   write(*,*)'MM Testing is symbol or state variable record',&
!        allocated(current%statvar)
   nostv: if(.not.allocated(current%statvar)) then
! an experiment is a symbol!!! Then statvar is not allocated
      symsym=current%statev
!      write(*,*)'MM A symbol, not a state variable for this experiment',symsym
! we must evaluate all state variable functions!!
      call meq_evaluate_all_svfun(-1,ceq)
! get the symbol name
      text=svflista(symsym)%name
      ip=len_trim(text)+1
!      text(ip-1:ip-1)='='
!      write(*,*)'MM experiment: ',text(1:ip),ip
   else
!      write(*,*)'MM This experiment has a state variable record',&
!           allocated(current%statvar),allocated(current%indices),iterm
      symsym=0
      svrrec=>current%statvar(1)
      call encode_state_variable(text,ip,svrrec,ceq)
      if(iterm.lt.current%noofterms) then
         iterm=iterm+1; goto 150
      endif
   endif nostv
!   write(*,*)'MM ok here',symsym
   if(current%experimenttype.eq.0 .or. current%experimenttype.eq.100) then
! write = followed by the value 
!      if(text(ip:ip).ne.' ') ip=ip+1
      text(ip:)='='
      ip=ip+1
   elseif(current%experimenttype.eq.-1) then
!      if(text(ip:ip).ne.' ') ip=ip+1
      text(ip:)='<'
      ip=ip+1
   elseif(current%experimenttype.eq.1) then
!      if(text(ip:ip).ne.' ') ip=ip+1
      text(ip:)='>'
      ip=ip+1
   endif
!   write(*,*)'MM experiment line 2: ',text(1:ip),ip
   if(current%symlink1.gt.0) then
! the value is a symbol
      text(ip:)=svflista(current%symlink1)%name
      ip=len_trim(text)+1
   else
!      call wrinum(text,ip,10,0,current%prescribed)
      call wrinum(text,ip,8,0,current%prescribed)
   endif
! uncertainty can also be a symbol
   text(ip:ip)=':'
   ip=ip+1
!   write(*,*)'MM experiment line 3: ',text(1:ip),ip,current%symlink2
   if(current%symlink2.gt.0) then
! the value is a symbol
      text(ip:)=svflista(current%symlink2)%name
      ip=len_trim(text)+1
   else
!      call wrinum(text,ip,10,0,current%uncertainty)
      call wrinum(text,ip,8,0,current%uncertainty)
   endif
!   write(*,*)'MM ok here 2',symsym,text(1:ip)
!   write(*,*)'MM experiment line 2: ',text(1:ip),ip
   if(current%experimenttype.eq.100) then
      text(ip:ip)='%'
      ip=ip+1
   endif
!   write(*,*)'MM ok here 3',symsym
! add the current value of the experiment after a $ sign
! TROUBLE GETTING WRONG VALUE HERE WHEN USER DEFINED REFERENCE STATES
   if(symsym.eq.0) then
      call state_variable_val(svrrec,xxx,ceq)
   else
!      write(*,*)'MM ok here 4',symsym
      actual_arg=' '
      xxx=evaluate_svfun_old(symsym,actual_arg,1,ceq)
   endif
   if(gx%bmperr.ne.0) then
! it is maybe a derivative ... 
!      write(*,*)'MM we cannot evaluate a derivative here ...',gx%bmperr
! but meq_evaluate_svfun not available here ... it is part of the minimizer
      gx%bmperr=0
      actual_arg=' '
      mode=1
      xxx=meq_evaluate_svfun(symsym,actual_arg,mode,ceq)
!      write(*,*)'MM meq_evaluate_svfun, mode=1: ',xxx
   endif
   if(gx%bmperr.ne.0) then
      write(*,*)'MM Error evaluating symbol: ',gx%bmperr
      text(ip:)=' $ ?? '
      ip=ip+5
      gx%bmperr=0
   else
!      write(*,*)'MM experimental state variable value: ',ip,xxx
      text(ip:)=' $'
      ip=ip+3
!      call wrinum(text,ip,12,0,xxx)
      call wrinum(text,ip,8,0,xxx)
!      write(*,*)'MM experiment line 3: ',text(1:ip),ip
   endif
!   write(*,*)'MM ok here 5'
1000 continue
!   write(*,*)'MM experiment line 4: ',text(1:ip),ip,gx%bmperr
   return
 end subroutine meq_get_one_experiment

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_equilibrium_extra
!\begin{verbatim}
 subroutine list_equilibrium_extra(lut,ceq,pun)
! list the extra character variables for calculate symboles and
! list characters (if any),  It is used in pmon and is part of matsmin
! because it calls subroutines which need access to calculated results
! If the first non-blank character of ceq%eqextra(3) is 0 (zero) then pun
! will be used as a file number to generate a plotfile with calculated values
   implicit none
   integer lut,pun
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ip,slen,jj,last,kk
   character tval*24,symbol*24,encoded*24,date*12
   double precision xxx,xarr(6)
!   write(*,*)'MM calc/list extra: ',ceq%eqname
!
   tval=' '
   symbol=' '
   extra: if(allocated(ceq%eqextra)) then
      ip=1
      if(eolch(ceq%eqextra(1),ip)) goto 190
!      write(*,*)'calc "',ceq%eqextra(1)(ip:len_trim(ceq%eqextra(1))),'"',lut
      calcs: if(ceq%eqextra(1)(ip:ip).ne.' ') then
! this line contains symbols to be calculated
!         write(*,*)'MM calc extra: ',ceq%eqextra(1)(1:len_trim(ceq%eqextra(1)))
         ip=ip-1
100      continue
! Third argument 2 means terminate at a space, not at a comma ","
! because some symbols may contain a comma
         call getext(ceq%eqextra(1),ip,2,tval,' ',slen)
         if(tval(1:1).ne.' ') then
! This is for a symbol that is not a dot derivative ...
!            call find_svfun(tval,istv,ceq)
            call meq_get_state_varorfun_value(tval,xxx,symbol,ceq)
!            mode=1
!            call meq_evaluate_svfun(tval,'  ',mode,ceq)
            if(gx%bmperr.ne.0) then
!               write(*,*)'MM Cannot find symbol: ',tval,' Error reset'
               gx%bmperr=0
            else
!            mode=1
! meq_evaluate_svfun is declared in matsmin
!            xxx=meq_evaluate_svfun(istv,'  ',mode,ceq)
!            xxx=evaluate_svfun_old(istv,'  ',mode,ceq)
!            if(gx%bmperr.ne.0) then
!               write(*,*)'MM Cannot calculate symbol: ',tval,' Error reset'
!               gx%bmperr=0; goto 100
!            endif
! symbol empty??
!               write(lut,110)symbol(1:len_trim(symbol)),xxx
               write(lut,110)tval(1:len_trim(tval)),xxx
110            format(3x,a,'=',1pe16.8)
            endif
            goto 100
!         else
!            write(*,*)'Found a space at position',ip
         endif
      endif calcs
190   continue
      ip=1
      if(eolch(ceq%eqextra(2),ip)) goto 290
      lists: if(ceq%eqextra(2)(ip:ip).ne.' ') then
! this line contains state variables and related things to be listed
!         write(*,*)'MM list extra: ',ceq%eqextra(2)(1:30)
         ip=ip-1
200      continue
! Third argument 2 means terminate at a space, not at a comma ","
! because some symbols contains a comma.
         call getext(ceq%eqextra(2),ip,2,tval,' ',slen)
!         write(*,*)'MM variable: ',tval,slen
         if(tval(1:1).ne.' ') then
            if(index(tval,'*').gt.0) then
               write(*,*)'MM Not implemented wildcards'
!            call get_many_svar(tval,...
            else
               symbol=' '
               call get_state_var_value(tval,xxx,symbol,ceq)
! This checks that the phase is stable ...
!               call get_stable_state_var_value(tval,xxx,symbol,ceq)
               if(gx%bmperr.ne.0) then
!                  write(*,*)'MM Cannot list variable: ',tval,' Error reset'
                  gx%bmperr=0
               else
                  write(lut,110)trim(symbol),xxx
               endif
            endif
            goto 200
         endif
      endif lists
290   continue
      ip=1
      if(eolch(ceq%eqextra(3),ip)) goto 390
      plots: if(ceq%eqextra(3)(ip:ip).eq.'0') then
! this creates a plot file for calculated values
! This is for plot_data set 0, calculated values'
! next value must be number of columns with data to be plotted!!
         last=ip+1
         call getint(ceq%eqextra(3),last,ip)
         if(buperr.ne.0) then
            write(*,*)'MM Cannot extract number of columns'
            gx%bmperr=4399; goto 1000
         endif
         if(ip.gt.6) then
            write(*,*)'MM Too many columns in plot_data 0. Max: 6',ip
            gx%bmperr=4399; goto 1000
         endif
         if(pun.eq.0) then
            pun=30
!            plotdatafile='oc_many0'
!            write(*,*)'Opening oc_many0.plt '
            open(pun,file='oc_many0.plt',access='sequential',status='unknown')
!
! extract state variable symbols, first is x axis variable
            kk=last
            call getext(ceq%eqextra(3),kk,2,tval,' ',slen)
            call date_and_time(date)
            write(pun,305)date(1:4),date(5:6),date(7:8),&
                 trim(tval),trim(ceq%eqextra(3))
305         format('# GUNPLOT file generated by enter many_equilibria '/&
                 'set title "Open Calphad 4.0 prerelease: ',a,'-',a,'-',a,&
                 ' with GNUPLOT"'/&
                 '# set terminal pdf color'/&
                 '# set output "whatever"'/&
                 'set xlabel "',a,'"'/&
                 'set ylabel "whatever"'/&
                 'set key bottom right'/&
                 '# ',a/&
                 '# THE DATA LINES MUST BE REPEATED AS MANY TIMES AS',&
                 ' THERE ARE PLOT COMMANDS!')
            call getext(ceq%eqextra(3),kk,2,tval,' ',slen)
            if(ip.eq.2) then
! with just two columns
               write(pun,310)trim(tval)
310            format('plot "-" using 1:2 with points pt 5 ',&
                    'ps 1.5 title "',a,'"')
            else
! this first line if 3 or more columns
               write(pun,311)trim(tval)
311            format('plot "-" using 1:2 with points pt 5 ',&
                    'ps 1.5 title "',a,'",\')
            endif
! if ip>4 this for second and further lines until jj is ip-1
            do jj=3,ip-1
               call getext(ceq%eqextra(3),kk,2,tval,' ',slen)
               write(pun,312)jj,jj+3,trim(tval)
312            format('"" using 1:',i2,' with points pt ',i2,&
                    ' ps 1.5 title "',a,'",\')
            enddo
! if ip=3 this is second line, otherwise the last line
            if(ip.gt.3) then
               call getext(ceq%eqextra(3),kk,2,tval,' ',slen)
               write(pun,313)ip,ip+3,trim(tval)
313            format('"" using 1:',i2,' with points pt ',i2,&
                    ' ps 1.5 title "',a,'"')
            endif
! the line consists of several state variables to be calculated and listed
            jj=0
320         continue
!            write(*,321)trim(ceq%eqextra(3)),last
321         format('3B extract: ',a,i5,' "',a,'"')
! 3rd argument 2 means skipping , only space separators
            call getext(ceq%eqextra(3),last,2,tval,' ',slen)
!            write(*,321)trim(ceq%eqextra(3)),last,trim(tval)
            if(tval(1:1).eq.' ') then
               goto 350
            elseif(buperr.ne.0) then
               write(kou,*)'Error reading symbol: ',trim(ceq%eqextra(3))
               goto 350
            else
               jj=jj+1
               call get_state_var_value(tval,xarr(jj),encoded,ceq)
               if(gx%bmperr.ne.0) then
                  write(*,*)'Error getting: ',tval
                  goto 350
               endif
            endif
            goto 320
! no more values
350         continue
         else
! This is another line with values for plot_data set 0, file is open
            jj=0
360         continue
            call getext(ceq%eqextra(3),last,2,tval,' ',slen)
            if(tval(1:1).eq.' ') then
               goto 370
            elseif(buperr.ne.0) then
               write(kou,*)'Error reading symbol: ',trim(ceq%eqextra(3))
               goto 370
            else
               jj=jj+1
               call get_state_var_value(tval,xarr(jj),encoded,ceq)
               if(gx%bmperr.ne.0) then
                  write(*,*)'Error getting: ',tval
                  goto 370
               endif
            endif
            goto 360
! no more values
370         continue
         endif
! write the line on the plot_data file
         if(jj.ne.ip) then
            write(*,*)'Wrong number of columns',jj,ip
         endif
         write(pun,380)(xarr(jj),jj=1,ip)
380      format(6(1pe12.4))
      endif plots
!   else
!      write(*,*)'No extra lines found'
   endif extra
390 continue
1000 continue
 end subroutine list_equilibrium_extra

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine equilph1a
!\begin{verbatim}
  subroutine equilph1a(phtup,tpval,ceq)
! equilibrates the constituent fractions of a phase using its current comp.
! phtup is phase tuple
! tpval is T and P
! ceq is a datastructure with all relevant thermodynamic data
    implicit none
    double precision tpval(*)
    TYPE(gtp_phasetuple), pointer :: phtup
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(meq_setup) :: meqrec
!\end{verbatim} %+
    integer nel,ii,errall
    double precision, allocatable :: xknown(:),wmass(:),cpot(:)
    double precision totmol,totmass,amount
    nel=noel()
    allocate(xknown(nel),stat=errall)
    allocate(wmass(nel),stat=errall)
    allocate(cpot(nel),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 50: ',errall
       gx%bmperr=4370; goto 1000
    endif
! find the current molefractions
!    call calc_phase_molmass(phtup%phaseix,phtup%compset,xknown,wmass,&
    call calc_phase_molmass(phtup%ixphase,phtup%compset,xknown,wmass,&
         totmol,totmass,amount,ceq)
    if(gx%bmperr.ne.0) goto 1000
! extract the current chemical potentials
    do ii=1,nel
       cpot=ceq%cmuval(ii)
    enddo
    if(gx%bmperr.ne.0) goto 1000
! create the meqrec structure
    call equilph1_meqrec(phtup,meqrec,.FALSE.,ceq)
    if(gx%bmperr.ne.0) goto 1000
    ceq%rtn=globaldata%rgas*tpval(1)
! iterate until equilibrium found for this phase
    call equilph1c(meqrec,meqrec%phr,tpval,xknown,cpot,ceq)
    deallocate(xknown)
    deallocate(wmass)
1000 continue
    return
  end subroutine equilph1a

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine equilph1b
!\begin{verbatim} %-
  subroutine equilph1b(phtup,tpval,xknown,gval,cpot,tyst,ceq)
! equilibrates the constituent fractions of a phase for mole fractions xknown
! phtup is phase tuple
! tpval is T and P
! ceq is a datastructure with all relevant thermodynamic data
! gval is the Gibbs energy calculated as xknown(i)*cpot(i)
! cpot are the (calculated) chemical potentials
! tyst is TRUE means no outut
    implicit none
!    integer mode
    TYPE(meq_setup) :: meqrec
    double precision tpval(*),xknown(*),cpot(*),gval
    TYPE(gtp_equilibrium_data), pointer :: ceq
    logical tyst
!\end{verbatim} %+
    TYPE(gtp_phasetuple), pointer :: phtup
    integer ii
! extract the current chemical potentials as start values
    do ii=1,noel()
       cpot(ii)=ceq%cmuval(ii)
    enddo
    if(gx%bmperr.ne.0) goto 1000
! create the meqrec structure
    call equilph1_meqrec(phtup,meqrec,.FALSE.,ceq)
    if(gx%bmperr.ne.0) goto 1000
! mabe we need RT ?
    ceq%rtn=globaldata%rgas*tpval(1)
! iterate until equilibrium found for this phase
    call equilph1c(meqrec,meqrec%phr,tpval,xknown,cpot,ceq)
!    write(*,*)'We are in equilph1b',gx%bmperr
    gval=zero
    if(gx%bmperr.eq.0) then
       do ii=1,noel()
          gval=gval+xknown(ii)*cpot(ii)
!          write(*,*)'We are in equilph1b',gval
       enddo
    endif
1000 continue
    return
  end subroutine equilph1b

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine equilph1c
!\begin{verbatim}
  subroutine equilph1c(meqrec,phr,tpval,xknown,ovar,ceq)
! iterate constituent fractions of a phase for mole fractions xknown
! tpval is T and P
! xknown are mole fractions
! ceq is a datastructure with all relevant thermodynamic data
! ovar are the chemical potentials
    implicit none
!    integer phase
    double precision tpval(*),xknown(*),ovar(*)
    TYPE(meq_setup) :: meqrec
    TYPE(meq_phase), dimension(*), target :: phr
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer nz1,nz2,converged,ierr,jj,nj,nk,nl,errall
    TYPE(meq_phase), pointer :: pmi
    double precision, allocatable :: smat(:,:),svar(:),yarr(:),ycorr(:)
    double precision chargefact,chargerr,pv,qq(5),ys,ycormax2
! number of variables is number of components + one stable phase
    nz1=meqrec%nrel+1
    nz2=nz1+1
    allocate(smat(nz1,nz2),stat=errall)
    allocate(svar(nz1),stat=errall)
!    allocate(ovar(nz1))
! current values of chemical potentials
!    do jj=1,meqrec%nrel
!       ovar(jj)=ceq%cmuval(jj)
!    enddo
    allocate(ycorr(phr(1)%ncc),stat=errall)
    allocate(yarr(phr(1)%ncc),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 51: ',errall
       gx%bmperr=4370; goto 1000
    endif
    chargefact=one
    chargerr=one
!    write(*,*)'We are in equilph1c: ',phr(1)%iph,phr(1)%ics,gx%bmperr
! we have just one phase in phr, phr must be TARGET
    pmi=>phr(1)
100 continue
    converged=0
    smat=zero
! invert the phase matrix for pmi
    call meq_onephase(meqrec,pmi,ceq)
    if(gx%bmperr.ne.0) goto 1000
! all ok to here ???
! setup mass balance equations, note some components may be missing
! This is a simplified setup_equilmatrix using xknown as composition
!    call setup_equilmatrix(meqrec,phr,nz1,smat,tcol,pcol,dncol,converged,ceq)
    call setup_comp2cons(meqrec,phr,nz1,smat,tpval,xknown,converged,ceq)
    if(gx%bmperr.ne.0) goto 1000
! debug output as the matrix had changed efter return from subroutine ...
!    do nk=1,nz1
!       write(*,111)'smat4: ',nk,(smat(nk,jj),jj=1,nz2)
!    enddo
!    goto 1000
! solve the equilibrium matrix, some chemical potentials may be missing
    call lingld(nz1,nz2,smat,svar,nz1,ierr)
    if(ierr.ne.0) then
       write(*,*)'Error solving equilibrium matrix 2',ierr
       gx%bmperr=4203; goto 1000
    endif
! check that svar(1..meqrec%nrel) has converged
    do jj=1,meqrec%nrel
       if(abs(svar(jj)-ovar(jj)).gt.1.0D1*ceq%xconv) then
!          write(*,103)'chempot7: ',svar(jj),ovar(jj),svar(jj)-ovar(jj)
103       format(a,3(1pe12.4))
          converged=7
       endif
! use ovar below to correct constitutions.  Note ovar is chem.pot/RT
       ovar(jj)=svar(jj)
    enddo
!    write(*,111)'svar4: ',0,(svar(jj),jj=1,nz1)
111 format(a,i2,6(1pe12.4))
! check dxmol ... seems OK
!    do nk=1,phr(1)%ncc
!       write(*,111)'dxmol: ',nk,(phr(1)%dxmol(nl,nk),nl=1,meqrec%nrel)
!    enddo
! update constituent fractions in just one phase
!    lap: do jj=1
    jj=1
! The current chemical potentials are in ceq%cmuval(i) svar(1..n)
! jj is stable, increment kk but do not make it larger than meqrec%nstph
! save index in meqrec%stphl in jph !!!!!!!!!!! kk never used !!!!!!!!!
!    jph=kk
!    kk=min(kk+1,meqrec%nstph)
! if phr(jj)%xdone=1 then phase has no composition variation
    if(phr(jj)%xdone.eq.1) goto 1000
!----------------------------------------------------
    ycormax2=zero
!    write(*,*)'cc: ',jj,phr(jj)%ncc
! loop for all constituents
    moody: do nj=1,phr(jj)%ncc
       ys=zero
       do nk=1,phr(jj)%ncc
          pv=zero
          do nl=1,meqrec%nrel
! ceq%cmuval(nl) is the chemical potential of element nl (divided by RT)
! USE values in svar(nl)
! phr(jj)%dxmol(nl,nk) is the derivative of component nl
! wrt constituent nk
!             pv=pv+ceq%complist(nl)%chempot(1)/ceq%rtn*phr(jj)%dxmol(nl,nk)
!             write(*,111)'pv1: ',nj,pv,ceq%complist(nl)%chempot(1),&
! ovar(nl) is used instead of complist(nl)%chempot(1) as we do not want to
! change the global values of the chemical potential
             pv=pv+ovar(nl)*phr(jj)%dxmol(nl,nk)
!             write(*,111)'pv1: ',nj,pv,ovar(nl),&
!                  ceq%rtn,phr(jj)%dxmol(nl,nk)
          enddo
!          write(*,119)'cph1: ',jj,nj,nk,ys,pv,phr(jj)%curd%dgval(1,nk,1),&
!               phr(jj)%invmat(nj,nk)
119       format(a,3i3,6(1pe12.4))
          pv=pv-phr(jj)%curd%dgval(1,nk,1)
          ys=ys+phr(jj)%invmat(nj,nk)*pv
!          write(*,111)'pv2: ',nj,ys,pv,phr(1)%curd%dgval(1,nk,1),&
!               phr(1)%invmat(nj,nk)
       enddo
       if(phr(jj)%chargebal.eq.1) then
! For charged phases add a term 
! phr(jj)%invmat(phr(jj)%idim,phr(jj)%idim)*Q
          ys=ys-chargefact*phr(jj)%invmat(nj,phr(jj)%idim)*&
               phr(jj)%curd%netcharge
       endif
       ycorr(nj)=ys
       if(abs(ycorr(nj)).gt.ycormax2) then
          ycormax2=ycorr(nj)
       endif
       if(abs(ys).gt.ceq%xconv) then
! if the change in any constituent fraction larger than xconv continue iterate
          if(converged.lt.4) then
! large correction in fraction of constituent fraction of stable phase
!             write(*,*)'mm converged 4B: ',jj,nj,ys
             converged=4
             cerr%mconverged=converged
             if(cerr%nvs.lt.10) then
                cerr%nvs=cerr%nvs+1
                cerr%typ(cerr%nvs)=4
                cerr%val(cerr%nvs)=zero
                cerr%dif(cerr%nvs)=abs(ys)
             endif
!             yss=ys
!             yst=phr(jj)%curd%yfr(nj)
          endif
!       elseif(phr(jj)%stable.eq.1) then
! check to find good convergence criteria in Re-V test case
!          if(abs(ycorr(nj)).gt.ysmm) then
!             jmaxy=jj
!             ysmm=abs(ycorr(nj))
!             ysmt=phr(jj)%curd%yfr(nj)
!           endif
       endif
       yarr(nj)=phr(jj)%curd%yfr(nj)+ycorr(nj)
!       write(*,119)'ycorr4: ',jj,nj,phr(jj)%chargebal,&
!            yarr(nj),phr(jj)%curd%yfr(nj),ycorr(nj),ys
    enddo moody
! >>>>>>>>>>>>>>>>>> HERE the new constitution is set <<<<<<<<<<<<<<<<<<<<<
!    write(*,112)'YC: ',jj,(ycorr(nj),nj=1,phr(jj)%ncc)
!    write(*,112)'YZ: ',meqrec%noofits,(yarr(nj),nj=1,phr(jj)%ncc)
112 format(a,i3,8F8.5)
!    write(*,*)'MM calling set_constitution 4: ',phr(jj)%iph,phr(jj)%ics
    call set_constitution(phr(jj)%iph,phr(jj)%ics,yarr,qq,ceq)
    if(gx%bmperr.ne.0) goto 1000
!  >>>>>>>>>>>>>>>>>> for all phases <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    meqrec%noofits=meqrec%noofits+1
    if(converged.gt.3) then
       if(meqrec%noofits.le.ceq%maxiter) goto 100
       write(*,*)'MM Too many iterations',ceq%maxiter
    elseif(meqrec%noofits.lt.6) then
       goto 100
    else
       if(.not.btest(meqrec%status,MMQUIET)) write(*,202)meqrec%noofits
202 format('Calculation required ',i4,' its')
    endif
1000 continue
    return
  end subroutine equilph1c

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine equilph1d
!\begin{verbatim}
  subroutine equilph1d(phtup,tpval,xknown,cpot,tyst,nend,mugrad,mobval,ceq)
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
    integer nend
    logical tyst
    double precision tpval(*),xknown(*),cpot(*),mugrad(*),mobval(*)
    TYPE(gtp_phasetuple), pointer :: phtup
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
    TYPE(meq_setup) :: meqrec
    integer ii
! extract the current chemical potentials as start values
    do ii=1,noel()
       cpot(ii)=ceq%cmuval(ii)
    enddo
    if(gx%bmperr.ne.0) goto 1000
! create the meqrec structure
!    write(*,17)'MM equilph1d calling equilph1e',(xknown(ii),ii=1,noel())
17  format(a,10(F6.3))
!    call equilph1_meqrec(phtup,meqrec,.FALSE.,ceq)
    call equilph1_meqrec(phtup,meqrec,tyst,ceq)
    if(gx%bmperr.ne.0) goto 1000
! mabe we need RT ?
    ceq%rtn=globaldata%rgas*tpval(1)
! iterate until equilibrium found for this phase
    call equilph1e(meqrec,meqrec%phr,tpval,xknown,cpot,tyst,&
         nend,mugrad,mobval,ceq)
1000 continue
    return
  end subroutine equilph1d

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine equilph1e
!\begin{verbatim} %-
  subroutine equilph1e(meqrec,phr,tpval,xknown,ovar,tyst,&
       noofend,mugrad,mobval,ceq)
! iterate constituent fractions of a phase for mole fractions xknown
! and calculate derivatives of MU and diffusion coefficients
! tpval is T and P
! xknown are mole fractions
! nrel is the number of components (elements)
! ovar are the chemical potentials
! tyst is TRUE if no output
! mugrad is the derivatives of the chemical potentials wrt mole fractions??
! mobval are the mobilities
! ceq is a datastructure with all relevant thermodynamic data
    implicit none
    integer noofend
    double precision tpval(*),xknown(*),ovar(*),mugrad(*),mobval(*)
    logical tyst
    TYPE(meq_setup) :: meqrec
    TYPE(meq_phase), dimension(*), target :: phr
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer nz1,nz2,converged,ierr,jj,nj,nk,nl,is,jt
    integer lokph,nkl(maxsubl),first(maxsubl+1),current(maxsubl),nsl,nend
    integer deriv(maxsubl),ql,mend
    TYPE(meq_phase), pointer :: pmi
    double precision, allocatable :: smat(:,:),svar(:),yarr(:),delta(:)
! dmuenddy is derivatives of mu for endmembers wrt all constituents
    double precision, allocatable :: dmuenddy(:,:),muend(:)
    double precision, allocatable :: py(:)
    double precision chargefact,chargerr,pv,qq(5),ys,ycormax2,muall
    double precision sumsum
! ************** change in MODEL_PARAMETER_IDENTIFIER: MQ is now 1300!!
! 800 + cs where cs is the constituent index counted over all sublattices ??
! can be REDEFINED when new model parameter identifiers was added!!! 
! we get the current value (set in gtp3A.F90) by calling getmqindex below
    integer mqindex,errall
! mqindex is a constant set in gtpini in models/gtp3A.F90
! number of variables is number of components + one stable phase
    nz1=meqrec%nrel+1
    nz2=nz1+1
    allocate(smat(nz1,nz2),stat=errall)
    allocate(svar(nz1),stat=errall)
!    allocate(ovar(nz1))
! current values of chemical potentials
!    do jj=1,meqrec%nrel
!       ovar(jj)=ceq%cmuval(jj)
!    enddo
    allocate(delta(phr(1)%ncc),stat=errall)
    allocate(yarr(phr(1)%ncc),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 52: ',errall
       gx%bmperr=4370; goto 1000
    endif
    chargefact=one
    chargerr=one
! we have just one phase in phr, phr must be TARGET
    pmi=>phr(1)
100 continue
    converged=0
    smat=zero
! invert the phase matrix for pmi
    call meq_onephase(meqrec,pmi,ceq)
    if(gx%bmperr.ne.0) goto 1000
! all ok to here ???
! setup mass balance equations, note some components may be missing
! This is a simplified setup_equilmatrix using xknown as composition
!    call setup_equilmatrix(meqrec,phr,nz1,smat,tcol,pcol,dncol,converged,ceq)
    call setup_comp2cons(meqrec,phr,nz1,smat,tpval,xknown,converged,ceq)
    if(gx%bmperr.ne.0) goto 1000
!    write(*,*)'after setup_comp2cons: ',converged
! debug output as the matrix had changed efter return from subroutine ...
!    do nk=1,nz1
!       write(*,111)'smat3: ',nk,(smat(nk,jj),jj=1,nz2)
!    enddo
!    goto 1000
! solve the equilibrium matrix, some chemical potentials may be missing
    call lingld(nz1,nz2,smat,svar,nz1,ierr)
    if(ierr.ne.0) then
       write(*,*)'Error solving equilibrium matrix 3',ierr
       gx%bmperr=4203; goto 1000
    endif
! check that svar(1..meqrec%nrel) has converged
    do jj=1,meqrec%nrel
       if(abs(svar(jj)-ovar(jj)).gt.1.0D1*ceq%xconv) then
!          write(*,103)'chempot: ',svar(jj),ovar(jj),svar(jj)-ovar(jj)
103       format(a,3(1pe12.4))
          converged=7
          cerr%mconverged=converged
          if(cerr%nvs.lt.10) then
             cerr%nvs=cerr%nvs+1
             cerr%typ(cerr%nvs)=7
             cerr%val(cerr%nvs)=svar(jj)
             cerr%dif(cerr%nvs)=ovar(jj)
          endif
       endif
       ovar(jj)=svar(jj)
    enddo
!    write(*,111)'svar3: ',0,(svar(jj),jj=1,nz1)
111 format(a,i2,6(1pe12.4))
! check dxmol ... seems OK
!    do nk=1,phr(1)%ncc
!       write(*,111)'dxmol: ',nk,(phr(1)%dxmol(nl,nk),nl=1,meqrec%nrel)
!    enddo
! update constituent fractions in just one phase
!    lap: do jj=1
    jj=1
! The current chemical potentials are in ceq%cmuval(i) svar(1..n)
! jj is stable, increment kk but do not make it larger than meqrec%nstph
! save index in meqrec%stphl in jph !!!!!!!!!!! kk never used !!!!!!!!!
!    jph=kk
!    kk=min(kk+1,meqrec%nstph)
! if phr(jj)%xdone=1 then phase has no composition variation
    if(phr(jj)%xdone.eq.1) goto 1000
!----------------------------------------------------
    ycormax2=zero
!    write(*,*)'cc: ',jj
! loop for all constituents
!    write(*,112)'Y0: ',meqrec%noofits,converged,(yarr(nj),nj=1,phr(jj)%ncc)
    moody: do nj=1,phr(jj)%ncc
       ys=zero
       do nk=1,phr(jj)%ncc
          pv=zero
          do nl=1,meqrec%nrel
! ceq%cmuval(nl) is the chemical potential of element nl (divided by RT)
! When a chemical potential is fixed use meqrec%mufixval
! phr(jj)%dxmol(nl,nk) is the derivative of component nl
! wrt constituent nk
!?             pv=pv+ceq%complist(nl)%chempot(1)/ceq%rtn*phr(jj)%dxmol(nl,nk)
!             pv=pv+ceq%cmuval(nl)*phr(jj)%dxmol(nl,nk)
!             pv=pv+svar(nl)*phr(jj)%dxmol(nl,nk)
             pv=pv+ovar(nl)*phr(jj)%dxmol(nl,nk)
          enddo
          pv=pv-phr(jj)%curd%dgval(1,nk,1)
          ys=ys+phr(jj)%invmat(nj,nk)*pv
!          write(*,111)'pv: ',nj,ys,pv,phr(1)%curd%dgval(1,nk,1),&
!               phr(1)%invmat(nj,nk)
       enddo
       if(phr(jj)%chargebal.eq.1) then
! For charged phases add a term 
! phr(jj)%invmat(phr(jj)%idim,phr(jj)%idim)*Q
          ys=ys-chargefact*phr(jj)%invmat(nj,phr(jj)%idim)*&
               phr(jj)%curd%netcharge
       endif
       delta(nj)=ys
       if(abs(delta(nj)).gt.ycormax2) then
          ycormax2=delta(nj)
       endif
       if(abs(ys).gt.ceq%xconv) then
! if the change in any constituent fraction larger than xconv continue iterate
          if(converged.lt.4) then
! large correction in fraction of constituent fraction of stable phase
!             write(*,*)'mm converged 4C: ',jj,nj,ys
             converged=4
!             yss=ys
!             yst=phr(jj)%curd%yfr(nj)
          endif
!       elseif(phr(jj)%stable.eq.1) then
! check to find good convergence criteria in Re-V test case
!          if(abs(delta(nj)).gt.ysmm) then
!             jmaxy=jj
!             ysmm=abs(delta(nj))
!             ysmt=phr(jj)%curd%yfr(nj)
!           endif
       endif
       yarr(nj)=phr(jj)%curd%yfr(nj)+delta(nj)
    enddo moody
! >>>>>>>>>>>>>>>>>> HERE the new constitution is set <<<<<<<<<<<<<<<<<<<<<
!    write(*,112)'YC: ',jj,(delta(nj),nj=1,phr(jj)%ncc)
!    write(*,112)'YY: ',meqrec%noofits,converged,(yarr(nj),nj=1,phr(jj)%ncc)
112 format(a,2i3,8F8.5)
!    write(*,*)'MM calling set_constitution 5:',phr(jj)%iph,phr(jj)%ics
    call set_constitution(phr(jj)%iph,phr(jj)%ics,yarr,qq,ceq)
    if(gx%bmperr.ne.0) goto 1000
!-------------------------- end of iteration
! check convergence
    meqrec%noofits=meqrec%noofits+1
    if(converged.gt.3) then
       if(meqrec%noofits.le.ceq%maxiter) goto 100
       gx%bmperr=4204
!       write(*,*)'MM Too many iterations',ceq%maxiter
       goto 1000
    elseif(meqrec%noofits.lt.6) then
       goto 100
    else
       if(.not.btest(meqrec%status,MMQUIET)) write(*,202)meqrec%noofits
202 format('Calculation required ',i4,' its')
    endif
    do is=1,meqrec%nrel
       ovar(is)=svar(is)
    enddo
!    goto 1000
!----------------------------------------------------------
! When the calculation converged we calculate mugrad and interdiffusivites
! A nontrival expression:
!
! dmu_i/dx_j = 1/N (d2G/dx_i/dx_j - \sum_k x_k (d2G/dx_k/dx_i + d2G/dx_k/dx_j)+
!                              \sum_k\sum_m x_k x_m d2G/dx_k/dx_m )
!
! NOTE THIS IS SYMMETRICAL, dmu_i/dx_j = dmu_j/dx_i.
! If the phase is ideal then d2G/dx_i/dx_j = RT/x_i if i=j, otherwise zero
! This gives for 
! dmu_i/dx_i = RT/N * (1-x_i)/x_i
! dmu_i/dx_j = - RT/N                  (i not equal to j)
!
! We calc             sum_k (x_k*d2G/dx_k/dx_i)   in delta(i)
!         sum_m x_m ( sum_k (x_k*d2G/dx_k/dx_m))  in sumsum
!
! new use of delta !!!
    delta=zero
    muall=pmi%curd%gval(1,1)
    sumsum=zero
! Here we calculate delta(is) =           \sum_jt y(jt)*d2G/dy_jt/dy_is and
!                   sumsum = \sum_m y(is) \sum_jt y(jt)*d2G/dy_jt/dy_is
! The loop of is is for all constituents
    do is=1,phr(1)%ncc
! The loop for jt are for all constituents in all sublattices
       do jt=1,phr(1)%ncc
! STRANGE that d2G/dy_Va/dy_Va is zero ... should be 1 (*RT) ...does not matter
!          if(is.gt.jt) stop "wrong order 1"
! keep ixsym here as I do not know if jt>is or not          
          delta(is)=delta(is)+pmi%curd%yfr(jt)*pmi%curd%d2gval(ixsym(is,jt),1)
!          write(*,*)'d2G/dy/dy: ',is,jt,pmi%curd%d2gval(ixsym(is,jt),1)
       enddo
       sumsum=sumsum+pmi%curd%yfr(is)*delta(is)
       muall=muall-pmi%curd%dgval(1,is,1)*pmi%curd%yfr(is)
    enddo
! muall    = G_m - \sum_i y_i dG/dy_i
! delta(i) = \sum_j y_j d2G/dy_i/dy_j             sum for all y_j for one y_i
! sumsum   = \sum_i \sum_j y_i y_j d2G/dy_i/dy_j  sum for all y_i and y_j
!-------------------- summations over all constituents in all sublattices
! now we must generate the endmembers, loop over all sublattices
! but sublattics and number of constituents in each are in the phase record
! and protected ... use a subroutine ...
    lokph=pmi%curd%phlink
    call get_phase_structure(lokph,nsl,nkl)
    if(gx%bmperr.ne.0) goto 1000
! ---------------------------------------------------------------------
    substitutional: if(nsl.eq.1) then
! specially simple if nsl=1 (substitutional)
       noofend=nkl(1)
       allocate(muend(noofend),stat=errall)
! calculate just mu(endmember)
!       loop1: do nend=1,noofend
!          muend(nend)=muall+pmi%curd%dgval(1,nend,1)
!          loop2: do jt=1,noofend
! the chemical potential has the derivative of the constituent
!             muend(nend)=muend(nend)+pmi%curd%dgval(1,jt,1)
!          enddo loop2
!       enddo loop1
! now we calculate dmu(end)/dy_is (just for substitutional)
       allocate(dmuenddy(noofend,pmi%ncc),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 53: ',errall
          gx%bmperr=4370; goto 1000
       endif
       dmuenddy=zero
! For a substitutional solution:
! dmu_i/dx_j = 1/N ( d2G/dx_i/dx_j -
!                    \sum_k x_k d2G/dx_k/dx_i - \sum_k x_k d2G/dx_k/dx_j+
!                    \sum_k\sum_m x_k x_m d2G/dx_k/dx_m )
! NOTE THIS SHOULD BE SYMMETRICAL, dmu_i/dx_j = dmu_j/dx_i.
! use delta(i) and sumsum calculated above
!       write(*,*)'Derivatives of chemical potentials',noofend
       nl=0
       loop3: do is=1,noofend
          muend(is)=muall+pmi%curd%dgval(1,is,1)
          loop4: do jt=1,noofend
!             if(is.gt.jt) stop "wrong order 2"
! keep using ixsym here as I do not know if jt>is
             dmuenddy(is,jt)=pmi%curd%d2gval(ixsym(is,jt),1)-&
                  delta(is)-delta(jt)+sumsum
!             write(*,775)'dd1:',1,is,jt,&
!                  dmuenddy(is,jt),pmi%curd%d2gval(ixsym(is,jt),1),&
!                  delta(is),delta(jt),sumsum
             nl=nl+1
             mugrad(nl)=dmuenddy(is,jt)*ceq%rtn
          enddo loop4
!          write(*,777)'dd: ',(ceq%rtn*dmuenddy(is,jt),jt=1,noofend)
!777       format(a,6(1pe12.4))
       enddo loop3
! UNFINISHED ?? I do not divide by N
!       write(*,777)'mu: ',(muend(is),is=1,noofend)
!-------------------
    else ! not substitutional below (2 or more sublattices)
! now we have to handle sublattices and endmembers
! nsl is number of sublattices and nkl(1..nsl) the number of const in each
       noofend=1
       is=1
       first=0
       do nl=1,nsl
! nend is number of endmembers
! here first and current are set to first constituent index in each sublattice
          noofend=noofend*nkl(nl)
          first(nl)=is
          current(nl)=is
          deriv(nl)=is
          is=is+nkl(nl)
       enddo
! we need this to indicate when we reached the end
       first(nsl+1)=is
       allocate(muend(noofend),stat=errall)
       allocate(py(noofend),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 54: ',errall
          gx%bmperr=4370; goto 1000
       endif
       py=one
!       write(*,611)'first: ',nsl,(first(nj),nj=1,nsl)
!611    format(a,i2,2x,10i3)
! all partials have this term
       muend=muall
!       write(*,*)'MM muall: ',muall,pmi%curd%gval(1,1)
! The partial Gibbs energy, for each sublattice add one dG/dy_is
       nend=0
       nj=0
       allpg: do while(nj.le.nsl)
          nend=nend+1
! the partials constituents, G_I, are in current(1..nsl)
          nlloop: do nl=1,nsl
             is=current(nl)
! endmembers like 1:1:1, 1:1:2, 1:2:1, 1:2:2, 2:1:1, 2:1:2, 2:2:1, 2:2:2 =8
! constituents are in current(1..nsl)
             muend(nend)=muend(nend)+pmi%curd%dgval(1,is,1)
          enddo nlloop
! generate a new set of constituents in current
          nj=1
888       continue
          current(nj)=current(nj)+1
          if(current(nj).eq.first(nj+1)) then
! note first(nsl+1) is the end of all constituents
             current(nj)=first(nj)
             nj=nj+1
             if(nj.le.nsl) goto 888
          endif
       enddo allpg
       if(.not.tyst) then
          write(*,881)(muend(jt),jt=1,noofend)
881       format('Calculated potentials for all endmembers/RT: '/6(1x,1pe12.4))
       endif
!-----------------------------------------------------------------------
! the part below is messy and unfinished
!---------------- now the derivative of the partial Gibbs energy
! The partial Gibbs energy, for each sublattice add one dG/dy_is
! the derivative of the partial Gibbs energy wrt all other endmembers ....
! dG_i/dn_J = 1/N_J( \sum_s (d2G/dy_is/dy_js - delta(is) - delta(js)) + sumsum )
! delta(is) = \sum_s \sum_k y_k d2G/dy_is/dy_k
! sumsum    = \sum_k \sum_m y_k y_m d2G/dy_k/dy_m (already added above)
!---------------------------------------------------
! all derivatives of the partial has the sumsum term
       allocate(dmuenddy(noofend,noofend),stat=errall)
       if(errall.ne.0) then
          write(*,*)'MM Allocation error 55: ',errall
          gx%bmperr=4370; goto 1000
       endif
       dmuenddy=sumsum
       nj=0
       nend=0
       allpartg: do while(nj.le.nsl)
! loop for all partial Gibbs energies G_I
          nend=nend+1
          mend=0
!          write(*,773)'Partial:   ',nend,(current(nl),nl=1,nsl)
          allql: do while(nj.le.nsl)
! loop for all constituent endmembers n_J
             mend=mend+1
!             write(*,773)'Endmember: ',mend,(deriv(nl),nl=1,nsl)
!773          format(a,i3,2x,10i3)
             lattloop: do nl=1,nsl
! loop for all sublattices, skip sublattices with a single constituent??
!                if(nkl(nl).eq.1) cycle lattloop
                is=current(nl)
! the 2nd derivative of G for constituents in same sublattice
                jt=deriv(nl)
                dmuenddy(nend,mend)=dmuenddy(nend,mend)-delta(is)-delta(jt)
! add second derivatives wrt is and all constituents in deriv
                suckloop: do ql=1,nsl
! keep using ixsym here as I do not know if is<deriv(ql)
!                   if(is.gt.deriv(ql)) stop "wrong order 3"
                   dmuenddy(nend,mend)=dmuenddy(nend,mend)+&
                        pmi%curd%d2gval(ixsym(is,deriv(ql)),1)
                enddo suckloop
! the amount of this endmember, should be calculated only once ...
                if(mend.eq.1) py(nend)=py(nend)*pmi%curd%yfr(is)
             enddo lattloop
!             dmuenddy(nend,mend)=dmuenddy(nend,mend)+sumsum
! update the derivative endmember
             nj=1
887          continue
             deriv(nj)=deriv(nj)+1
             if(deriv(nj).eq.first(nj+1)) then
! note first(nsl) is the end of all constituents
                deriv(nj)=first(nj)
                nj=nj+1
                if(nj.le.nsl) goto 887
             endif
          enddo allql
! update the partitial Gibbs energy endmember
          nj=1
886       continue
          current(nj)=current(nj)+1
          if(current(nj).eq.first(nj+1)) then
! note first(nsl) is the end of all constituents
             current(nj)=first(nj)
             nj=nj+1
             if(nj.le.nsl) goto 886
          endif
       enddo allpartg
       nl=0
       loop7: do is=1,noofend
          loop8: do jt=1,noofend
             nl=nl+1
             mugrad(nl)=dmuenddy(is,jt)*ceq%rtn
          enddo loop8
!          write(*,705)'dmu: ',(dmuenddy(is,jt),jt=1,noofend)
!705       format(a,6(1pe12.4))
       enddo loop7
    endif substitutional
!-------------------
! D_kj = \sum_i\sum_s (delta_ki - y_ks) y_is M_i dmu_i/dy_j
! this should be calculated for the components ... but I have just endmembers 
!-------------------
! UNFINISHED calculation of diffusivities
! I can calculate D(end,jt) = py(end)*exp(mob(end))*dmudy(end,jt)
! in the database is stored mq&constituent#sublattice
! I will calculate the mob(end) as \sum_s \sum_c mq&c#s taking
! those values missing as zero ... ???
! the values of mq&c#s are in pmi%curd%gval(1,itp) where itp is
! ************** change in MODEL_PARAMETER_IDENTIFIER: MQ is now 1300!!
! 800 + cs where cs is the constituent index counted over all sublattices ??
! 1300 + cs where cs is the constituent index counted over all sublattices ??
! list additional properties:
!    write(*,400)'props: ',(pmi%curd%listprop(jt),jt=2,pmi%curd%listprop(1)-1)
!400 format(/a,12i6)
! instead of 800 use the function mqindex('MQ  ')
    ql=0
!    write(*,*)'In equi1ph1d: ',mqindex
!    jt=getmqindex()
    mqindex=get_mpi_index('MQ  ')
    if(gx%bmperr.ne.0) then
       write(*,*)'MM mqindex error: ',gx%bmperr,mqindex
       goto 1000
    endif
! note that MQ has a composition index so it must be multiplied by 100
    mqindex=100*mqindex
!
    do jt=2,pmi%curd%listprop(1)
       is=pmi%curd%listprop(jt)
       if(is.gt.mqindex .and. is.lt.mqindex+100) then
! there is a mobility for constituent (is-mqindex) stored in pmi%curd%gval(1,jt)
          jj=is-mqindex
          ql=ql+1
!          mobval(jj)=exp(pmi%curd%gval(1,jt)/ceq%rtn)/ceq%rtn
          mobval(jj)=pmi%curd%gval(1,jt)
!          write(*,410)is,jj,jt,pmi%curd%gval(1,jt)
!410       format('MM Mobility for ',2i4,' in pos ',i2,', value: ',3(1pe14.6))
       endif
    enddo
    if(ql.ne.meqrec%nrel) then
       write(*,*)'MM: WARNING found ',ql,' mobilities values out of',&
            meqrec%nrel
    endif
! we do not have mobility values for all endmembers, only for the number
! of components
!    if(ql.lt.noofend) then
!       write(*,411)noofend-ql,noofend
!411    format(' *** Warning EQUILPH1E: Missing mobility data for ',i2,&
!            ' endmembers: ',i3)
!       goto 1000
!    endif
    goto 1000
! NO CALCULATION OF DIFFUSIVITIES HERE, JUST RETURN MOBILITY VALUES
! list T and x for current values
!    write(*,412)tpval(1),(pmi%curd%yfr(jt),jt=1,3)
!412 format(/'Unreduced diffusion matrix for T= ',f8.2,' and x= ',3F8.4)
!
! TC gives for MU(i).x(j) ....:
!
! The loop below is adapted to the FCC phase in the AlCuSi system
! 2 sublattices but only substitutional diffusion
! Diffs are D_kj
!    allocate(diffs(3,3))
!
! I calculate D(is,jt) = x_is * exp(M_is/RT) * (dmu_is/dx_jt) 
!                      - x_is * \sum_nl x_nl * exp(M_nl/RT)* (dmu_nl/dy_jt)
!
!    diffs=zero
!    nend=0
!    do is=1,noofend
!       do jt=1,noofend
!          sumsum=zero
!          do nl=1,noofend
! note yarr(nl) is exp(M_nl/RT)/RT
! dumenddy is also divided by RT ...
!             sumsum=sumsum-pmi%curd%yfr(nl)*yarr(nl)*dmuenddy(nl,jt)*ceq%rtn
!          enddo
!          diffs(is,jt)=pmi%curd%yfr(is)*(yarr(is)*dmuenddy(is,jt)*ceq%rtn+&
!               sumsum)
!          nend=nend+1
!          intdiv(nend)=diffs(is,jt)
!       enddo
!       write(*,414)is,(diffs(is,jt),jt=1,3)
!    enddo
!414 format('D_kj, k=',i2,' j=1..3 ',4(1pe14.6))
!
!    write(*,415)
!415 format(/'Taking Al as reference we ger D^Al_kj = D_ij - D_1j ')
!    do is=2,3
!       write(*,416)is,(diffs(is,jt)-diffs(1,jt),jt=2,3)
!    enddo
!416 format('D_kj, k=',i2,' j=2..3 ',2(1pe16.6))
!
1000 continue
    return
  end subroutine equilph1e

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine equilph1_meqrec
!\begin{verbatim}
  subroutine equilph1_meqrec(phtup,meqrec,tyst,ceq)
!  subroutine equilph1b(phtup,tpval,xknown,cpot,tyst,ceq)
! equilibrates the constituent fractions of a phase for mole fractions xknown
! phtup is phase tuple
! tpval is T and P
! ceq is a datastructure with all relevant thermodynamic data
! cpot are the (calculated) chemical potentials
! tyst is TRUE means keep quiet
    implicit none
    integer mode,errall
    TYPE(meq_setup) :: meqrec
!    double precision tpval(*),xknown(*),cpot(*)
    TYPE(gtp_equilibrium_data), pointer :: ceq
    logical tyst
!\end{verbatim}
    TYPE(gtp_phasetuple), pointer :: phtup
! setup equilibrium calculation for a single phase, set all others as suspended
! store values in meqrec
    meqrec%nrel=noel()
    meqrec%nfixph=0
    meqrec%nfixmu=0
    meqrec%tpindep=.FALSE.
    meqrec%nphase=1
    allocate(meqrec%phr(1),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 56: ',errall
       gx%bmperr=4370; goto 1000
    endif
    meqrec%nstph=1
! wrong?? phaseix is index in phases, ixphase is index in phlista
!    meqrec%phr(1)%iph=phtup%phaseix
    meqrec%phr(1)%iph=phtup%ixphase
    meqrec%phr(1)%ics=phtup%compset
    meqrec%phr(1)%itadd=0
    meqrec%phr(1)%itrem=0
    meqrec%phr(1)%xdone=0
    meqrec%phr(1)%phasestatus=1
    meqrec%phr(1)%ionliq=-1
    meqrec%phr(1)%i2sly=0
    meqrec%stphl(1)=1
!    if(test_phase_status_bit(phtup%phaseix,PHIONLIQ)) meqrec%phr(1)%ionliq=1
    if(test_phase_status_bit(phtup%ixphase,PHIONLIQ)) meqrec%phr(1)%ionliq=1
! set link to calculated values of G etc.
!    call get_phase_compset(iph,ics,lokph,lokcs)
! link to results
    meqrec%phr(1)%curd=>ceq%phase_varres(phtup%lokvares)
! set phase stable
    meqrec%phr(1)%stable=1
    meqrec%phr(1)%prevam=one
    meqrec%phr(1)%prevdg=zero
    meqrec%phr(1)%idim=0
! number of constituents !!!
    meqrec%phr(1)%ncc=size(ceq%phase_varres(phtup%lokvares)%yfr)
    meqrec%dormlink=0
    meqrec%status=0
    if(tyst) then
       meqrec%status=ibset(meqrec%status,MMQUIET)
    else
       meqrec%status=ibclr(meqrec%status,MMQUIET)
    endif
!
    meqrec%noofits=0
! this replaces call to meq_sameset as we will never change stable phase
!    call equilph1c(meqrec,meqrec%phr,tpval,xknown,cpot,ceq)
1000 continue
    return
  end subroutine equilph1_meqrec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!-\addtotable subroutine set_eec_check
!-\begin{verbatim}
!  subroutine set_eec_check(tval)
! This set values for EEC check, called from user i/f or application software
! ceq is a datastructure with all relevant thermodynamic data
!    implicit none
!    double precision tval
!-\end{verbatim} %+
!    if(tval.gt.1.0D1) then
!       globaldata%sysreal(1)=tval
!    else
!       globaldata%sysreal(1)=zero
!    endif
!    return
!  end subroutine set_eec_check

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine check_eec
!\begin{verbatim}
  subroutine check_eec_old(pmisol,pmiliq,meqrec,ceq)
! This checks EEC after calculating all phases if the solid phase has S > S^liq
! it is called if T>globaldata%sysreal(1) (set in user i/f)
! pmisol is pointer to solid data
! pmiliq is pointer to liquid data
! ceq is a datastructure with all relevant thermodynamic data
    implicit none
    type(meq_phase), pointer :: pmiliq,pmisol
    TYPE(gtp_equilibrium_data), pointer :: ceq
    TYPE(meq_setup) :: meqrec
!\end{verbatim}
    integer sel
    double precision newg,ssol,sliq,fact,kvot
    logical :: once=.TRUE.
    save once
! check if T<globaldata%sysreal(1) already made
! Calculate:  -S^sol_m - (-S^liq_m):
    write(*,*)'MM we should never call check_eec! '
    if(.not.associated(pmiliq)) then
       if(.not.associated(pmisol)) then
          write(*,*)'MM check_eec called without any phases!'
       elseif(once) then
! This message written only once
          write(*,*)' *** WARNING EEC method fails as no liquid'
          once=.FALSE.
       endif
       goto 1000
    endif
! abnorm(1) is the number of atoms per formula units
    ssol=-pmisol%curd%gval(2,1)/pmisol%curd%abnorm(1)
    sliq=-pmiliq%curd%gval(2,1)/pmiliq%curd%abnorm(1)
    fact=sliq/ssol
    if(fact.lt.one) then
! fact<0 means solid has higher entropy than liquid
! set G.T and G.T.T for solid to those for liquid, do not care about G.Y.T
! note that values must be adjusted to number of real atoms in the phases
       kvot=pmisol%curd%abnorm(1)/pmiliq%curd%abnorm(1)
       pmisol%curd%gval(2,1)=pmiliq%curd%gval(2,1)*kvot
       pmisol%curd%gval(4,1)=pmiliq%curd%gval(4,1)*kvot
! The Gibbs energy for solid is set to G = \sum_i x_i \mu_i plus
! 10000*(number of moles of atoms of solid) to ensure it is unstable
! The mole fractions of the solid are in pmisol => meqrec%phr(ij)
       newg=zero
       do sel=1,meqrec%nrel
          newg=newg + pmisol%xmol(sel)*ceq%complist(sel)%chempot(1)
       enddo
! For Al-50%Cr a pure Al-bcc becomes almost stable at 3500 K using 10000
!       pmisol%curd%gval(1,1)=(newg+1.0E4)*pmisol%curd%abnorm(1)/ceq%rtn
       pmisol%curd%gval(1,1)=(newg+1.0E3)*pmisol%curd%abnorm(1)/ceq%rtn
! but tested for 5000 it still does not become stable.
!       pmisol%curd%gval(1,1)=(newg+5.0E3)*pmisol%curd%abnorm(1)/ceq%rtn
!    else
! Entropy of solid less than liquid, all is OK
    endif
1000 continue
!    write(*,*)'MM leaving check_eec'
    return
  end subroutine check_eec_old
  
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine tzero
!\begin{verbatim}
  subroutine tzero(iph1,iph2,icond,value,ceq)
! calculates the value of condition "icond" for two phases to have same G
    implicit none
    integer iph1,iph2,icond
    double precision value
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer, parameter :: lwa=100
    integer info,nv,ja,jb
    type(gtp_condition), pointer :: first
    type(gtp_phase_varres), pointer :: cps1,cps2
    double precision xv(5),fvec(5),tol,wa(lwa)
!    external tzcalc NOT NEEDED
!
!    write(*,'(a,3i4," and ",3i4)')'In tzero!',iph1,phasetuple(iph1)%ixphase,&
!         phasetuple(iph1)%lokph,&
!         iph2,phasetuple(iph2)%ixphase,phasetuple(iph2)%lokph
! in some way ceq, iph1, iph2 and icond must be transferred to tzcalc
    tzceq=>ceq
! always use first composition set, iph is also index in phasetuple
!    tzph1=phasetuple(iph1)%lokph; tzph2=phasetuple(iph2)%lokph
    tzph1=iph1; tzph2=iph2
! find the condition
    first=>ceq%lastcondition
    tzcond=>first%next
    ja=0
    do while(.not.associated(first,tzcond))
! we should only count ACTIVE conditions
       if(tzcond%active.eq.0) then
          ja=ja+1
          if(ja.eq.icond) goto 100
          tzcond=>tzcond%next
       endif
    enddo
! the loop above does not find the last condition !!! SUCK
    if(icond.ne.ja+1) then
       write(*,*)'No such condition'
       gx%bmperr=4399; goto 1000
!    else
! the last condition was the selected one
    endif
!
100 continue
! Set status of all phases except iph1 and iph2 as suspended
    call change_many_phase_status('* ',-3,zero,ceq)
    call change_phtup_status(iph1,1,one,ceq)
    call change_phtup_status(iph2,1,one,ceq)
    if(gx%bmperr.ne.0) goto 1000
! start value of condition to vary
    xv(1)=tzcond%prescribed
!    write(*,*)'Found condition, current value ',xv(1)
! do we need to think about parallelization?
! calculate the zero
!    write(*,*)'Calling hybrd1',xv(1)
    nv=1
! testing tzero calculation with larger composition difference in the phases?
!    tol=1.0D-2  this is max difference in G, maybe relative??
    tol=1.0D-6
    call hybrd1(tzcalc,nv,xv,fvec,tol,info,wa,lwa)
    if(info.ne.1) then
! info=0 Improper input parameters
!     =2 Too many iterations 
!     =3 tol variable too small
!     =4 Too slow progress
!       write(*,*)'HYBRD solver return error: ',info
       if(gx%bmperr.eq.0) gx%bmperr=4371
    else
    endif
    if(gx%bmperr.ne.0) goto 1000
    tzcond%prescribed=xv(1)
    value=xv(1)
1000 continue
! restore suspeded phases and set no equilibrium
    ceq%status=ibset(ceq%status,EQINCON)
    call change_many_phase_status('* ',0,zero,ceq)
    if(gx%bmperr.eq.0) then
! set amount of the two phases
       cps1=>tzceq%phase_varres(phasetuple(tzph1)%lokvares)
       cps2=>tzceq%phase_varres(phasetuple(tzph2)%lokvares)
       cps1%amfu=one
       cps2%amfu=one
    endif
    return
  end subroutine tzero

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine tzcalc
!\begin{verbatim}
  subroutine tzcalc(nv,xv,fvec,iflag)
! calculates the value of a condition for two phases to have same G
    implicit none
    integer nv,iflag
    double precision xv(*),fvec(*)
!\end{verbatim}
    type(gtp_phase_varres), pointer :: cps1,cps2
    integer mode,lokph1,lokph2,lokvares1,lokvares2
    double precision gm1,gm2
! we transfer the data needed by tzph1,tzph2,tzceq and tzcond !!
!    write(*,*)'In tzcalc: ',tzph1,tzph2
!    write(*,*)'Current value of condition: ',tzcond%prescribed,xv(1)
!
    lokph1=phasetuple(tzph1)%lokph
    lokph2=phasetuple(tzph2)%lokph
    lokvares1=phasetuple(tzph1)%lokvares
    lokvares2=phasetuple(tzph2)%lokvares
    cps1=>tzceq%phase_varres(lokvares1)
    cps2=>tzceq%phase_varres(lokvares2)
!
    mode=0
! we have to calculate each phase separately and compare G values (per atom)
! Set current value of condition
    tzcond%prescribed=xv(1)
!    write(*,*)'Prescribed condition: ',tzcond%prescribed
! on entry both phases 1 and 2 are entered, suspend phase 2
    tzcond%prescribed=xv(1)
!    write(*,*)'Prescribed condition: ',tzcond%prescribed
! on entry both phases 1 and 2 are entered, suspend phase 2
    call change_phtup_status(tzph2,-3,one,tzceq)
    call calceq3(mode,.FALSE.,tzceq)
    if(gx%bmperr.ne.0) goto 1100
    gm1=cps1%gval(1,1)/cps1%abnorm(1)
!    write(*,*)'Phase 1: ',gm1
! suspend phase 1 and restore 2
    call change_phtup_status(tzph1,-3,one,tzceq)
    call change_phtup_status(tzph2,1,one,tzceq)
    call calceq3(mode,.FALSE.,tzceq)
    if(gx%bmperr.ne.0) goto 1100
    gm2=cps2%gval(1,1)/cps2%abnorm(1)
!    write(*,*)'Phase 2: ',gm2
! restore phase 1
    call change_phtup_status(tzph1,1,one,tzceq)
    fvec(1)=gm1-gm2
! maybe relative? No as gm1 and gm2 are divided by RT and around 1.0
!    write(*,'(a,4(1pe12.4))')'tzcalc: ',xv(1),gm1,gm2,fvec(1)
!    fvec(1)=cps1%gval(1,1)/cps1%abnorm(1)-cps2%gval(1,1)/cps2%abnorm(1)
1000 continue
    return
1100 continue
! error quit also calling routine by setting value to zero
    write(*,*)'Quit tzcalc due to error: ',gx%bmperr
    fvec(1)=zero; goto 1000
  end subroutine tzcalc

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine tzcalc_stoich
!\begin{verbatim}
  subroutine tzcalc_stoich(nv,xv,fvec,iflag)
! calculates the value of a condition for two phases to have same G
! called from smp2A during mapping
! Both phases have the same composition (stoichiometric constraits)
    implicit none
    integer nv,iflag
    double precision xv(*),fvec(*)
!\end{verbatim}
    type(gtp_phase_varres), pointer :: cps1,cps2
    integer mode,lokph1,lokph2,lokvares1,lokvares2,moded
    double precision gm1,gm2
! we transfer the data needed by tzph1,tzph2,tzceq and tzcond !!
!    write(*,*)'In tzcalc: ',tzph1,tzph2
!    write(*,*)'Current value of condition: ',tzcond%prescribed,xv(1)
!
! Tested for U-O with 3 phase ORTHO_A20, TETRA_U C1_MO2:  T=942.02 K
    lokph1=phasetuple(tzph1)%lokph
    lokph2=phasetuple(tzph2)%lokph
    lokvares1=phasetuple(tzph1)%lokvares
    lokvares2=phasetuple(tzph2)%lokvares
    cps1=>tzceq%phase_varres(lokvares1)
    cps2=>tzceq%phase_varres(lokvares2)
!
    mode=0
! only G values, no derivatives
    moded=0
! we have two phases with fixed composition and search T for the same value
! of the Gibbs energy. 
! We can directly calculate the Gibbs energy of each phase
! Set current value of condition
!    tzcond%prescribed=xv(1)
!    write(*,*)'Prescribed condition: ',tzcond%prescribed
! on entry both phases 1 and 2 are entered, suspend phase 2
!    call change_phtup_status(tzph2,-3,one,tzceq)
!    call calceq3(mode,.FALSE.,tzceq)
!    if(gx%bmperr.ne.0) goto 1100
    if(xv(1).lt.1.0D-1) then
!       write(*,*)'Attempt to calculate for T less than 1'
       iflag=-1; gx%bmperr=4187; goto 1000
    endif
    tzceq%tpval(1)=xv(1)
    call calcg_internal(lokph1,moded,cps1,tzceq)
    if(gx%bmperr.ne.0) goto 1200
    call calcg_internal(lokph2,moded,cps2,tzceq)
    if(gx%bmperr.ne.0) goto 1200
    gm1=cps1%gval(1,1)/cps1%abnorm(1)
!    write(*,*)'Phase 1: ',gm1
! suspend phase 1 and restore 2
!    call change_phtup_status(tzph1,-3,one,tzceq)
!    call change_phtup_status(tzph2,1,one,tzceq)
!    call calceq3(mode,.FALSE.,tzceq)
!    if(gx%bmperr.ne.0) goto 1100
    gm2=cps2%gval(1,1)/cps2%abnorm(1)
    fvec(1)=gm1-gm1
!    write(*,*)'Phase 2: ',gm2
! restore phase 1
    call change_phtup_status(tzph1,1,one,tzceq)
    fvec(1)=gm1-gm2
!    write(*,'(a,4(1pe12.4))')'tzcalc_stoich: ',xv(1),gm1,gm2,fvec(1)
1000 continue
    return
1100 continue
! error quit also calling routine by setting value to zero
    write(*,*)'Quit tzcalc_stoich due to error: ',gx%bmperr
    iflag=-1; goto 1000
1200 continue
    write(*,1210)gx%bmperr
1210 format('Error calculating Gibbs energy ',i5)
    iflag=-1; goto 1000
  end subroutine tzcalc_stoich

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_paraeq
!\begin{verbatim}
 subroutine calc_paraeq(tupix,icond,xcond,ceq)
! calculates a paraequilibrium between two phases tupix(1&2)
! icond is the index of the fast diffusing element
! xcond are the fractions of the element in the two phases at paraequilibrium
   implicit none
   integer tupix(2),icond
   double precision xcond(2)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! at paraequilibrium the two phases has the same composition (set as conditions)
! except for one element fastel which is a fast diffusion element (such as C)
! It requires solving a nonlinear equation to find a "tie-line" between
! the two phases which have the same composition except for fastel
! We have two variables, the composion of "fastel" in each phase
! We calculate each phase separately with different fractions of fastel, x(C)
! and extract the chemical potential of fastel, mu(C)
! and a "combined" chemical potential for all other elements.
! This is calculated as (G-x(C)*mu(C))/(1-x(C)) where G is the Gibbs energy
! The two function values are the difference of these two potentials
! calculated for each phase
!
   integer nv,info,ja
   integer, parameter :: lwa=20,minus1=-1
   double precision fracs(2),fvec(2),wa(lwa),muval,xsave,ntot,nalpha,nbeta,xtest
   double precision, parameter :: tol=1.0D-10
   type(gtp_phasetuple), pointer :: ph1,ph2
   type(gtp_condition), pointer :: first,pcond
   type(gtp_state_variable), target :: p1svr,p2svr
   type(gtp_state_variable), pointer :: svr
   character encoded*24,fractions*64,elname*24
   logical verbose
! We must passing links and info to paraeqfun, the subroutine called by hybrd1
! THIS DOES NOT WORK IF CALCULATIONS ARE MADE IN PARALLEL
! tzceq is pointer to equilibrium; tzcond pointer to first condition
   verbose=.FALSE.
   xcond=zero
   tzceq=>ceq
   tzph1=tupix(1); tzph2=tupix(2)
!
!   write(*,*)'MM calc_paraeq',tzph1,tzph2,icond
! Calculate an equilibrium with the two phases
   call calceq3(minus1,verbose,tzceq)
   if(gx%bmperr.ne.0) goto 1000
! exctract various values
   call get_state_var_value('N ',ntot,encoded,tzceq)
   if(gx%bmperr.ne.0) goto 1000
! this is quite clumsy ... can it be fixed?
   call get_component_name(icond,elname,tzceq)
! do not use this routine, requires a tuple record
!   call get_phasetuple_name(tzph1,encoded)
   call get_phasetup_name(tzph1,encoded)
   fractions='X('//trim(encoded)//','//trim(elname)//') '
   call get_state_var_value(fractions,fracs(1),encoded,tzceq)
!   call get_state_var_value('X(FCC,C) ',fracs(1),encoded,tzceq)
   if(gx%bmperr.ne.0) goto 1000
   call get_phasetup_name(tzph2,encoded)
   fractions='X('//trim(encoded)//','//trim(elname)//') '
   call get_state_var_value(fractions,fracs(2),encoded,tzceq)
!   write(*,*)'MM fraction composition: ',trim(fractions)
!   call get_state_var_value('X(BCC,C) ',fracs(2),encoded,tzceq)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'MM Initial fractions: ',fracs(1),fracs(2)
! tzcond should be the condition for the element tzel
   first=>ceq%lastcondition%next
   pcond=>first
   ja=0
   findxcond: do while(.true.)
      if(pcond%active.eq.1) cycle findxcond
      ja=ja+1
      if(pcond%statvar(1)%argtyp.eq.1) then
         if(pcond%statvar(1)%component.eq.icond) then
            tzcond=>pcond
            p1svr=pcond%statvar(1)
            p2svr=pcond%statvar(1)
! this is the condition on the total amount of fast diffusing element ??
            xsave=pcond%prescribed
!            write(*,'(a,F10.6,3i4)')'MM fraction condition',pcond%prescribed,&
!                 pcond%statvar(1)%statevarid,pcond%statvar(1)%oldstv
         endif
      endif
!      write(*,*)'MM other conditions: ',ja,pcond%statvar%statevarid,&
!           pcond%statvar%component
      pcond=>pcond%next
      if(associated(pcond,first)) exit findxcond
      if(ja.gt.100) then
         write(*,*)'Eternal loop exit 1',ja
         gx%bmperr=4399; goto 1000
      endif
   enddo findxcond
! musvr and xsvr are module global variables used by hybrid subroutine
! musvr is typically
! 3 0 0 0 0 1 0 0 1 0 1.0 3
   musvr%statevarid=3; musvr%norm=0; musvr%unit=0; musvr%phref=0
   musvr%argtyp=1; musvr%phase=0; musvr%compset=0; musvr%component=icond
   musvr%constituent=0; musvr%coeff=one; musvr%oldstv=3
   svr=>musvr
   call state_variable_val(svr,muval,ceq)
   if(gx%bmperr.ne.0) goto 1000
! we should use mole fractions to calculate alloy potential
! The %oldstv is important!  But completely undocumented
!=========================================================
! This assumes fraction is mole fraction
   xsvr=musvr
   xsvr%statevarid=17; xsvr%oldstv=111
! DO NOT CHANGE xsvr, IT IS USED IN THE CALCULATING ROUTINE
!=========================================================
   nv=2
!   call list_conditions(kou,tzceq)
!
!==================================
!
! solve the non-linear equation (this is the simplified call ....)
   call hybrd1(paraeqfun,nv,fracs,fvec,tol,info,wa,lwa)
! nv number of variables and functions; fracs(nv) values of the fractions
! fvec(nv) returned values of the functions; tol required tolerance
! info returned information of result
! ws is workspace with dimension lwa; lwa integer > nv*(3*n+13)/2 (=2*19/2)
   if(info.ne.1) then
      if(info.eq.0) write(*,*)'MM hybrd1 called with illegal arguments'
      if(info.eq.2) write(*,*)'MM hybrd1 fails too many iterations'
      if(info.eq.3) write(*,*)'MM hybrd1 fails too high tolerance required'
      if(info.eq.4) write(*,*)'MM hybrd1 fails too slow progress'
      gx%bmperr=4399; goto 1000
   endif
! the phase amounts should be adjusted to a composition in the middle
   xsave=0.5*(fracs(1)+fracs(2))
! return solution:
!   write(*,*)'MM conditions at the solution:'
!   call list_conditions(kou,tzceq)
   xcond(1)=fracs(1)
   xcond(2)=fracs(2)
! We should set the phase amounts to reproduce the overall condition
   if(xcond(1).gt.xcond(2)) then
      nalpha=(xsave-xcond(2))/(xcond(1)-xcond(2))
   else
      nalpha=(xsave-xcond(1))/(xcond(2)-xcond(1))
   endif
   nbeta=ntot-nalpha
   if(nalpha.lt.zero .or. nbeta.lt.zero) then
      write(*,'(a,5(1x,F10.6))')'Paraequil error:',xsave,xcond,nalpha,nbeta
      gx%bmperr=4399
   else
! set amounts of phases correspondng to condition
!      write(*,'(a,2F10.6)')'calc_paraeq: NP(*): ',nalpha,nbeta
      call change_phase_status(phasetuple(tzph1)%ixphase,&
           phasetuple(tzph1)%compset,PHENTSTAB,nalpha,tzceq)
      call change_phase_status(phasetuple(tzph2)%ixphase,&
           phasetuple(tzph2)%compset,PHENTSTAB,nbeta,tzceq)
   endif
!
1000 continue
! restore original condition
   tzcond%prescribed=xsave
   return
 end subroutine calc_paraeq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine paraeqfun
!\begin{verbatim}
 subroutine paraeqfun(nv,fracs,fvec,iflag)
! called by hydrid1 to solve a nonlinear system of equations setup
! by calc_paraeq to calculate the difference in chemical potential 
! for a tow-phase paraequilibrium.  Arguments are:   
! nv number of variables, fracs the variable values, fvec the functions
! calculated by this routine
   implicit none
   integer nv,iflag
   double precision fracs(*),fvec(*)
!\end{verbatim}
   integer, parameter :: minus1=-1
   double precision gm,mucmat,muamat,mucgro,muagro,xcmat,xcgro,mutest,xtest,val
   type(gtp_state_variable), pointer :: svr
   integer ip
   character encoded*24
   logical verbose
! 
! The 2 variables are the fractions of the fast diffusing element in 2 phases
! The functions are the chemical potential of the fast diffusing element
! and the "extrapolated" chemical potential of an alloy with zero fraction
! of the fast diffusing element, calculated for each phases as
!    (G-x(C)*mu(C))/(1-x(C))
! where G is the Gibbs energy of the phase and x(C) the fraction of C
! The difference of these potentials calculated for each element in each phase
! should be zero at paraequilibrium.
! THIS ROUTINE DOES NOT WORK IF CALCULATIONS IN PARALLEL
!
! At paraequilibrium the two phases has the same composition (set as conditions)
! except for one element fastel which is a fast diffusion element (such as C)
! It requires solving a nonlinear equation to find a "tie-line" between
! the two phases which have the same composition except for fastel
! We have two variables, the composion of "fastel" in each phase
! We calculate each phase separately with different fractions of fastel, x(C)
! and extract the chemical potential of fastel, mu(C)
! and a "combined" chemical potential for all other elements.
! This is calculated as (G-x(C)*mu(C))/(1-x(C)) where G is the Gibbs energy
! The two function values are the difference of these two potentials
! calculated for each phase
!
! iflag should not be changed except to force termination by setting iflag=-1
! NOTE tzceq, tzcond, tzph1 and tzph2 global variables in this module!
! fractions must be betwee 1E-12 and 1
   if(fracs(1).lt.1.0D-12) fracs(1)=1.0D-12 
   if(fracs(1).gt.1.0D0) fracs(1)=1.0D0
   if(fracs(2).lt.1.0D-12) fracs(2)=1.0D-12 
   if(fracs(2).gt.1.0D0) fracs(2)=1.0D0
!   write(*,'(a,2(1pe12.4))')'>>>> Paraeqfun 1: ',fracs(1),fracs(2)
! calculate the Gibbs energy and the partial Gibbs energies of each phase 
! for the current set of conditions
! It is possible to calculate each phase separately ignoring conditions
! but then it is not trivial to obtain the chemical potentials
! or call calceq2/calceq7 with all phses except one suspended
!
   verbose=.FALSE.
!   write(*,*)'Matrix and growing phases: ',tzph1,tzph2
! suspend growing phase and calculate normal equilibrium for matrix
   call change_phase_status(phasetuple(tzph2)%ixphase,&
        phasetuple(tzph2)%compset,PHSUS,zero,tzceq)
   if(gx%bmperr.ne.0) goto 1000
! set condition on composition equal to fracs(1)
   tzcond%prescribed=fracs(1)
!   do ip=1,nooftup()
!      if(test_phase_status(phasetuple(ip)%ixphase,phasetuple(ip)%compset,&
!           val,tzceq).ge.0) then
!         write(*,*)'Stable phase (matrix)',phasetuple(ip)%ixphase,&
!              phasetuple(ip)%compset,val
!      endif
!   enddo
!   write(*,*)'Calculating with matrix phase',phasetuple(tzph1)%ixphase,&
!        phasetuple(tzph1)%compset
! calceq3 will modify the fraction of components not set as conditions (Fe)
   call calceq3(minus1,verbose,tzceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'Failed calculation for matrix phase',fracs(1),gx%bmperr
      goto 1000
   endif
! extract value of G and MU(C) and calculate M(X)=(G-x(c)*mu(C))/(1-x(C))   
!   write(*,*)'Extracting values for matrix phase'
   call get_state_var_value('GM ',gm,encoded,tzceq)
   if(gx%bmperr.ne.0) goto 1000
   svr=>musvr
   call state_variable_val(svr,mucmat,tzceq)
   if(gx%bmperr.ne.0) goto 1000
   svr=>xsvr
   call state_variable_val(svr,xcmat,tzceq)
   if(gx%bmperr.ne.0) goto 1000
! check
!   call get_state_var_value('MU(C) ',mutest,encoded,tzceq)
!   if(gx%bmperr.ne.0) goto 1000
!   call get_state_var_value('X(C) ',xtest,encoded,tzceq)
!   if(gx%bmperr.ne.0) goto 1000
!   write(*,'(a,4(1pe12.4))')'Matrix test: ',mucmat,mutest,xcmat,xtest
!
   muamat=(gm-xcmat*mucmat)/(one-xcmat)
!   write(*,'(a,4(1pe12.4))')'Matrix G, x and mu:  ',gm,xcmat,mucmat,muamat
!
! suspend matrix phase and calculate normal equilibrium for growing!
   call change_phase_status(phasetuple(tzph2)%ixphase,&
        phasetuple(tzph2)%compset,PHENTSTAB,one,tzceq)
   if(gx%bmperr.ne.0) goto 1000
   call change_phase_status(phasetuple(tzph1)%ixphase,&
        phasetuple(tzph1)%compset,PHSUS,zero,tzceq)
   if(gx%bmperr.ne.0) goto 1000
!
!   set condition on composition equal to fracs(2)
   tzcond%prescribed=fracs(2)
!   do ip=1,nooftup()
!      if(test_phase_status(phasetuple(ip)%ixphase,phasetuple(ip)%compset,&
!           val,tzceq).ge.0) then
!         write(*,*)'Stable phase (growing)',phasetuple(ip)%ixphase,&
!              phasetuple(ip)%compset,val
!      endif
!   enddo
!   write(*,*)'Calculating with growing phase: ',phasetuple(tzph2)%ixphase,&
!        phasetuple(tzph2)%compset
! calceq3 will modify the fraction of components not set as conditions (Fe)
   call calceq3(minus1,verbose,tzceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'Failed calculation for growing phase',fracs(2),gx%bmperr
      goto 1000
   endif
! extract value of G and MU(C) and calculate M(X)=(G-x(c)*mu(C))/(1-x(C))   
!   write(*,*)'Extracting values for growing phase'
!   call get_stable_state_var_value('GM ',gm,encoded,tzceq)
   call get_state_var_value('GM ',gm,encoded,tzceq)
   if(gx%bmperr.ne.0) goto 1000
   svr=>musvr
   call state_variable_val(svr,mucgro,tzceq)
   if(gx%bmperr.ne.0) goto 1000
   svr=>xsvr
   call state_variable_val(svr,xcgro,tzceq)
   if(gx%bmperr.ne.0) goto 1000
! test
!   call get_state_var_value('MU(C) ',mutest,encoded,tzceq)
!   if(gx%bmperr.ne.0) goto 1000
!   call get_state_var_value('X(C) ',xtest,encoded,tzceq)
!   if(gx%bmperr.ne.0) goto 1000
!   write(*,'(a,4(1pe12.4))')'Matrix test: ',mucgro,mutest,xcgro,xtest
! we have to use mole fraction, to calculate muamat
   muagro=(gm-xcgro*mucgro)/(one-xcgro)
!   write(*,'(a,4(1pe12.4))')'Growing G, x and mu: ',gm,xcgro,mucgro,muagro

   fvec(1)=muamat-muagro
   fvec(2)=mucmat-mucgro
!   write(*,'(a,4(1pe12.4))')'>>>> Paraeqfun 9: ',fracs(1),fracs(2),&
!        fvec(1),fvec(2)
! restore matrix as entered
   call change_phase_status(phasetuple(tzph1)%ixphase,&
        phasetuple(tzph1)%compset,PHENTSTAB,one,tzceq)
   if(gx%bmperr.ne.0) goto 1000
!
1000 continue
   if(gx%bmperr.ne.0) then
      write(*,*)'MM Error inside paraeqfun',gx%bmperr
      iflag=-1
   endif
!   iflag=-1
   return
 end subroutine paraeqfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calculate_carefully
!\begin{verbatim}
  subroutine calculate_carefully(mode,ceq)
    implicit none
! calculate an equilirium carefully (bosses_method)
! step 1: Calculate with gridminimizer and merge (already done)
!         Alternatively enter with a set of stable phases
!         which has converged at another calculation.
!      2: suspend unstable phases
!      3: calculate with iterative method, 
!      4: set all suspended as dormant, 
!      5: calculate iterative again to see if any dormant has dgm>0
!      6: if so set it entered and goto 5
!      7: set all phases entered
! mode 0 means all step, nonzero may change some of these steps
! mode 1 means set phases entered one by one, largest driving force first
    integer mode
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer, allocatable, dimension(:) :: phcsstat
    integer ntups,mtups,itup,lokcs,ns,naa,phmax,saverr,errall
    double precision dgmax
    character phname*24
    logical again
! 2. loop for all phases to suspend all not stable
! 3. calculate equilibrium with remaining phases without gridminimizer
! 4. loop for all phases to set suspended to dormant
! 5. calculate equilibrium with current set of  phases without gridminimizer
! 6. If a dormant phase has dgm>0 set it entered and go back to 5
!    Alternatively do this one by one
! 7. set all dormant phases entered
!    ntups=noofphasetuples()
! if error 4363 then reset error code and continue
    if(gx%bmperr.ne.0) then
! this error means a phase has been restored with positive dgm
!       write(*,*)'MM error code set 1:',gx%bmperr
       if(gx%bmperr.ne.4363) goto 1000
       gx%bmperr=0
    endif
    ntups=nooftup()
    allocate(phcsstat(ntups),stat=errall)
    if(errall.ne.0) then
       write(*,*)'MM Allocation error 57: ',errall
       gx%bmperr=4370; goto 1000
    endif
    phcsstat=0
    ns=0
    do itup=1,ntups
       lokcs=phasetuple(itup)%lokvares
       if(ceq%phase_varres(lokcs)%dgm.lt.zero) then
! suspend phases with negative dgm
! if already suspended ignore
          if(ceq%phase_varres(lokcs)%phstate.gt.PHSUS) then
             ceq%phase_varres(lokcs)%phstate=PHSUS
             phcsstat(itup)=PHSUS
             ns=ns+1
          endif
       endif
    enddo
    write(*,12)ns
12  format('Phases set suspended except those found stable by gridmin',i5)
! calculate equilibrium with just these phases without grid minimizer
    call calceq2(0,ceq)
    if(gx%bmperr.ne.0) then
! this error means a phase has been restored with positive dgm
!       write(*,*)'MM error code set 2:',gx%bmperr
       if(gx%bmperr.ne.4363) goto 900
       gx%bmperr=0
    endif
    write(*,17)ntups-ns,ns,' suspended '
17  format('MM Equilibrium calculated with ',i4,' entered and ',i4,a,'phases')
! set all suspended phases as dormant, maybe some has disapperard
! some composition sets may have disappeared
    mtups=nooftup()
    if(ntups-mtups.gt.0) write(*,18)ntups-mtups
18  format('MM deleted ',i3,' composition sets')
    do itup=1,mtups
       if(phcsstat(itup).eq.PHSUS) then
! set a suspended phase as dormant
          lokcs=phasetuple(itup)%lokvares
          ceq%phase_varres(lokcs)%phstate=PHDORM
          phcsstat(itup)=PHDORM
       endif
    enddo
    write(*,19)
19  format('Calculating again with all suspended phases set as dormant')
! calculate equilibrium with entered and dormant phases
100 continue
    again=.false.
    call calceq2(0,ceq)
    if(gx%bmperr.ne.0) then
! this error means a phase has been restored with positive dgm
       write(*,*)'MM error code set 3:',gx%bmperr
       if(gx%bmperr.ne.4363) goto 900
       gx%bmperr=0
    endif
    write(*,17)mtups-ns,ns,' dormant '
! if mode=0 set all dormant phases with dgm>0 entered
! if mode=1 set the dormant phase with largest dgm>0 as entered
    ntups=nooftup()
    naa=0
    dgmax=zero
    phmax=0
    do itup=1,ntups
       if(phcsstat(itup).eq.PHDORM) then
          lokcs=phasetuple(itup)%lokvares
! maybe enter phases one by one ...
          if(ceq%phase_varres(lokcs)%dgm.gt.zero) then
             if(mode.eq.0) then
                ceq%phase_varres(lokcs)%phstate=PHENTERED             
                phcsstat(itup)=PHENTERED
                ns=ns-1
                naa=naa+1
                again=.true.
             elseif(ceq%phase_varres(lokcs)%dgm.gt.dgmax) then
                dgmax=ceq%phase_varres(lokcs)%dgm
                phmax=itup
             endif
          endif
       endif
    enddo
    if(mode.eq.1 .and. phmax.gt.0) then
!       write(*,*)'MM entering phase with largest driving force ',phmax,dgmax
       lokcs=phasetuple(phmax)%lokvares
       ceq%phase_varres(lokcs)%phstate=PHENTERED             
       phcsstat(phmax)=PHENTERED
       call get_phasetup_name(phmax,phname)
       ns=ns-1
       naa=naa+1
       again=.true.
       write(*,200)trim(phname),dgmax
200    format('MM Setting ',a,' with dgm= ',1pe12.4,' as entered')
    endif
    if(again) then
       if(mode.eq.0) write(*,*)'MM set ',naa,' dormant phases as entered'
       goto 100
    endif
! we have found a solution, set all phases as entered
! or we have an error so restore suspended phases
900 continue
    if(gx%bmperr.ne.0) then
       write(*,*)'Calculation not converged, some phases remain as dormant'
       goto 1000
    endif
    saverr=gx%bmperr
    gx%bmperr=0
    ns=0
    ntups=nooftup()
!    ntups=noofphasetuples()
    do itup=1,ntups
       if(phcsstat(itup).le.PHDORM) then
          lokcs=phasetuple(itup)%lokvares
          ceq%phase_varres(lokcs)%phstate=PHENTERED             
          ns=ns+1
       endif
    enddo
    gx%bmperr=saverr
    if(ns.gt.0) write(*,'(a,i4,a)')'MM Remaining ',ns,' phases set as entered'
1000 continue
    return
  end subroutine calculate_carefully

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine listoptshort
!\begin{verbatim}
  subroutine listoptshort(lut,mexp,nvcoeff,errs)
! short listing of optimizing variables and result
    integer lut,mexp,nvcoeff
    double precision errs(*)
!    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: neweq
    integer i1,i2,j1,j2,j3,neq
    character name1*24,line*80
    double precision xxx,sum
    type(gtp_condition), pointer :: experiment
!
! list all experiments, only possible if there are experiments
    if(mexp.eq.0) then
       write(lut,666)
666    format(/'No experiments so no results'/)
       goto 1000
    endif
! list experiments, mexp is number of EXPERIMENTS, not equilibria!!
    write(lut,620)size(firstash%eqlista),mexp
620 format(/'List of ',i5,' equilibria with ',i5,&
         ' experimental data values'/&
         '  No Equil name      Weight Experiment $ calculated',18x,&
         'Error')
    j3=0
    allequil: do i1=1,size(firstash%eqlista)
! skip equilibria with zero weight
       neweq=>firstash%eqlista(i1)%p1
       if(neweq%weight.eq.zero) cycle allequil
       name1=neweq%eqname(1:12)
! LOOP for all experiments for this equilibrium (maybe none??)
       if(.not.associated(neweq%lastexperiment)) cycle allequil
       experiment=>neweq%lastexperiment%next
       if(.not.associated(experiment)) cycle allequil
700    continue
          i2=neweq%lastexperiment%seqz
!          write(*,*)'number of experiments: ',i2
          neq=neweq%eqno
          do j2=1,i2
! j1 is position in line to write experiment
             j1=1
             line=' '
! this subroutine returns experiment and calculated value: "H=1000:200 $ 5000"
             call meq_get_one_experiment(j1,line,j2,neweq)
             j3=j3+1
             if(neq.gt.0) then
!                write(lut,622)neq,name1(1:12),neweq%weight,line(1:40),errs(j3)
                write(lut,622)neq,name1(1:15),neweq%weight,line(1:39),errs(j3)
622             format(i4,1x,a,2x,F5.2,1x,a,1x,1pe12.4)
                neq=0
             else
                write(lut,623)line(1:40),errs(j3)
623             format(28x,a,1x,1pe12.4)
             endif
! list the equilibrium name just for the first (or only) experiment
!             name1=' '
          enddo
          experiment=>experiment%next
!590       if(.not.associated(experiment,neweq%lastexperiment)) then
!             experiment=>experiment%next
!             goto 700
          if(j2.lt.i2 .and. .not.associated(experiment)) then
             write(*,*)'Missing experiment in equilibrium ',neweq%eqno
             cycle allequil
          endif
    enddo allequil
! list sum of squares
    sum=zero
    do j1=1,mexp
       sum=sum+errs(j1)**2
    enddo
! same as PARROT
    j1=mexp-nvcoeff
    if(j1.gt.0) then
       write(lut,621)sum,mexp,nvcoeff,j1,sum/j1
    else
       write(lut,621)sum,mexp,nvcoeff,0,zero
    endif
621 format(/'Final sum of squared errors: ',1pe16.5,&
         ' using ',i4,' experiments and'/&
         i3,' coefficient(s).  Degrees of freedom: ',i4,&
         ', normalized error: ',1pe13.4/)
1000 continue
    return
  end subroutine listoptshort

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine calctrans
!\begin{verbatim}
  subroutine calctrans(cline,last,ceq)
! calculate a phase transition
    character cline*(*)
    integer last
    type(gtp_equilibrium_data), pointer :: ceq
!\end
    character name1*30
    integer j1,iph,ics
    double precision xxx
    type(gtp_condition), pointer :: pcond
    type(gtp_state_variable), pointer :: stvr
!
    write(kou,2090)
2090 format('To calculate when a phase will appear/disappear',&
          ' by releasing a condition.')
    if(btest(ceq%status,EQNOEQCAL)) then
       write(kou,2095)
2095   format('You must make an equilibrium calculation before using',&
            ' this command.')
       goto 1000
    endif
    call gparcx('Phase name: ',cline,last,1,name1,' ','?Calculate transform')
    call find_phase_by_name(name1,iph,ics)
    if(gx%bmperr.ne.0) goto 1000
    j1=test_phase_status(iph,ics,xxx,ceq)
    if(j1.eq.PHFIXED) then
       write(kou,*)'Phase status already fixed'
       goto 1000
    endif
    call list_conditions(kou,ceq)
    write(kou,2097)
2097 format('You must release one condition, give its number')
    call gparidx('Condition number',cline,last,j1,1,'?CALCULATE transform')
    if(j1.le.0 .or. j1.gt.noel()+2) then
       write(kou,*)'No such condition'
       goto 1000
    endif
! this finds condition with given number
    call locate_condition(j1,pcond,ceq)
    if(gx%bmperr.ne.0) goto 1000
    if(pcond%active.eq.0) then
! the condition is active, deactivate it!
       pcond%active=1
    else
       write(kou,*)'This condition is not active!'
       goto 1000
    endif
! Condition released, now set the phase as fix with zero moles
    call change_phase_status(iph,ics,PHFIXED,xxx,ceq)
    if(gx%bmperr.ne.0) goto 1000
! Calculate equilibrium
    call calceq2(1,ceq)
    if(gx%bmperr.ne.0) goto 1000
! get the value of the released condition and set it to the new value
    stvr=>pcond%statvar(1)
    call state_variable_val(stvr,xxx,ceq)
    if(gx%bmperr.ne.0) goto 1000
    write(kou,2099)xxx
2099 format('The transition occurs at ',1pe16.8,', set as condition')
    pcond%prescribed=xxx
    pcond%active=0
! set phase back as entered and stable
!    write(*,*)'Set phase back as entered'
    call change_phase_status(iph,ics,PHENTSTAB,zero,ceq)
1000 continue
    return
  end subroutine calctrans

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\addtotable subroutine list_stable_phases
!\begin{verbatim}
  subroutine list_stable_phases(text,its,iadd,irem,meqrec,ceq)
! debug listing of stable phases
! meqrec contains all necessary data ...
    character*(*) text
    integer its,iadd,irem
    type(meq_setup) :: meqrec
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer ii,ij,ik
    double precision gsum,xmol(50),wmass(50),totmol,totmass
! NOTE: phases in stphl should always be in increasing order!!
    call calc_molmass(xmol,wmass,totmol,totmass,ceq)
    if(gx%bmperr.ne.0) stop 'Error when calling list_stable_phases'
    gsum=zero
    do ii=1,meqrec%nrel
       gsum=gsum+xmol(ii)*ceq%complist(ii)%chempot(1)
    enddo
    write(*,100)text,its,iadd,irem,meqrec%nstph,gsum,totmol,&
         (meqrec%stphl(ii),ii=1,meqrec%nstph)
    do ii=2,meqrec%nstph
       if(meqrec%stphl(ii-1).gt.meqrec%stphl(ii)) then
          stop 'phases in wrong order!!'
       endif
    enddo
100 format(a,i3,2i4,i3,2(1pe12.4),20i4)
!    do ii=1,meqrec%nstph,5
!       iph1=meqrec%stphl(ii=1,meqrec%nstph)
!       write(*,200)meqrec%phr(
!    enddo
    return
  end subroutine list_stable_phases

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine two_stoich_same_comp
!\begin{verbatim}
  recursive subroutine two_stoich_same_comp(irem,iadd,mapx,meqrec,inmap,ceq)
! we have found two  phases stable with same composition
! ONLY USED WHEN MAPPING with tie-lines in plane
! ceq is equilibrium record
    implicit none
    integer irem,iadd,inmap,mapx
    type(meq_setup) :: meqrec
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
!    type(map_node), pointer :: mapnode,newnode,oldnext
!    type(map_line), pointer :: nodexit
    type(gtp_equilibrium_data), pointer :: newceq
    integer nrel,nel,iph,ics,jj,seqx,phfix,lokph,lokcs,phstable,phfixtupix
    type(gtp_condition), pointer :: pcond,lastcond
! needed for solving a nonlinear equation
    integer, parameter :: lwa=100
    type(gtp_state_variable), target :: axstv1
    type(gtp_state_variable), pointer :: axstv
    integer nv,info
    double precision newphfra,fvec(5),tol,wa(lwa),value,xv(5),tinit
    character phases*48
!    logical isotherm
    integer idum,jdum,savefix(2),saveent
!    
!    write(*,*)'In two_stoich_same_comp'
    if(meqrec%nrel.ne.2) then
! How to check if I should use this routine? Only 2 components?
! If we have an activity condition one could have 3 components ....
       write(*,*)'MM This routine should  be used only when tie-lines in plane'
       gx%bmperr=4399; goto 1000
    endif
!    call get_state_var_value('X(O) ',value,phases,ceq)
!    write(*,806)meqrec%fixph(1,1),meqrec%fixph(2,1),mapx,iadd,value
806 format('MM why fix phase/set: ',i3,i2,' entered: ',i3,', new fix: ',i3,&
         1pe12.4)
!    write(*,*)'MM two stoichiometric phases: ',iadd,irem
! new T calculated in this routine should be close to current value
    tinit=ceq%tpval(1)
!    write(*,22)'MM in two_stoich_same_comp: ',irem,iadd,ceq%tpval(1)
!22  format(/20('-')/a,2i5,F8.2)
!    call list_conditions(kou,ceq)
! We cannot calculate an equilibrium with two phases with exactly the same
! composition.  But we can calculate the T where the two stoichiometric
! phases have the same Gibbs energy using the calc_tzero routine!
! Assuming the conditions are not too involved ... but we are dealing with a
! system with tie-lines in the plane, binary or ternary.
! use the vaiables tzph1 and tzph2 (in matsmin) to specify the phases involved
! DOES NOT WORK IN PARALLEL!!
    tzph1=irem; tzph2=iadd
    phases=' '
    call get_phasetup_name(tzph1,phases)
    nv=len_trim(phases)
    call get_phasetup_name(tzph2,phases(nv+2:))
!
!    write(*,27)'MM two compounds: ',tzph1,tzph2,trim(phases)
27  format(a,2i4,2x,a)
    nv=1
    tol=1.0D-6
! hybrid1 can solve a system of nonlinear equations by calling
! subroutine tzcalc_stoich(nv,xv,fvec,iflag) is in matsmin.F90
! the tzceq is a pointer declared in matmin and used in tzcalc_stoich
    tzceq=>ceq
    xv(1)=tzceq%tpval(1)
    call hybrd1(tzcalc_stoich,nv,xv,fvec,tol,info,wa,lwa)
    if(info.ne.1) then
! info=0 means improper input parameters
!     =2 Too many calls to tzcalc_stoich
!     =3 tol is too small
!     =4 Convergence too slow
!       write(*,*)'HYBRD solver return error: ',info
       if(gx%bmperr.eq.0) gx%bmperr=4371
    endif
    if(gx%bmperr.ne.0) goto 1000
    if(abs(ceq%tpval(1)-tinit).gt.2.0D1) then
       write(*,654)ceq%tpval(1),tinit
654    format('MM Error, too large change in T: ',2F10.2)
       gx%bmperr=4399; goto 1000
    endif
! To have correct chemical potentials we must call meq_sameset again
! But with T fix and phase iadd set dormant
! Now set current value of T as condition
!    call list_conditions(kou,ceq)
! loop all conditions until we find T and set it active.
! Maybe remove some other condition??
    lastcond=>ceq%lastcondition
    pcond=>lastcond
    nv=0
    jdum=0
    condloop1: do while(.TRUE.)
! loop for all conditions
       nv=nv+1
!       write(*,*)'State variable: ',nv,pcond%statev,pcond%prescribed
       if(pcond%statev.eq.1) then
! This is T, the axis condition, set as active with calculated value of T
          pcond%prescribed=xv(1)
          if(pcond%active.ne.0) then
             write(*,*)'Error, the condition on T not inactivated!'
             gx%bmperr=4399; goto 1000
          endif
          pcond%active=0
          jdum=nv
       else
          if(gx%bmperr.ne.0) then
             write(*,*)'Error extraction axis state variable value',gx%bmperr
             goto 1000
          endif
       endif
       pcond=>pcond%next
       if(associated(pcond,lastcond)) exit condloop1
    enddo condloop1
    if(jdum.eq.0) then
       write(*,*)'Error, no condition on T!'
       gx%bmperr=4399; goto 1000
    endif
! extract which phase is fixed (only one)
    savefix(1)=meqrec%fixph(1,1)
    savefix(2)=meqrec%fixph(1,2)
! and which is entered
    jdum=0
    do jj=1,meqrec%nphase
       if(meqrec%phr(jj)%stable.eq.1) then
          if(meqrec%phr(jj)%iph.eq.savefix(1) .and. &
               meqrec%phr(jj)%ics.eq.savefix(2)) then
!             write(*,*)'MM Fix phase: ',meqrec%fixph(1,1),meqrec%fixph(1,2)
             cycle
          endif
          if(jdum.eq.0) then
!             write(*,*)'MM Entered phase',jdum,jj
             jdum=jj
!          elseif(jj.ne.irem) then
!             write(*,*)'MM More than one entered phase',jdum,jj
          endif
       endif
    enddo
! we must keep saveent to return the entered phase when generating exits!
    saveent=jdum
!    write(*,*)'MM old fix phase/set and entered: ',meqrec%fixph(1,1),&
!         meqrec%fixph(2,1),saveent
!meq_sameset and ignore any change of the set of stable phases
! We must call meq_sameset again to have correct chemical potential at this T
!    write(*,*)'MU(*) before meq_sameset: ',ceq%cmuval(1),ceq%cmuval(2)
! Now we have calculated T when both stoichiometric phases are stable
! and set this T as condition. 
! set the phase iadd as suspend to avoid it will try to be stable
    meqrec%phr(iadd)%phasestatus=PHSUS
    meqrec%noofits=0
!    call list_conditions(kou,ceq)
! Strange here we have one degree of freedom! how can we calculate?  No check!!
! But we must have a condition on the amount
! mapx set to zero inside this routine.  Make sure no error code set!!
    if(gx%bmperr.ne.0) gx%bmperr=0
!    write(*,*)'MM calling meq_sameset from two_stoich_same_comp'
!   write(*,*)'This is a recursive call as we call two_stoich from meq_sameset!'
    call meq_sameset(idum,jdum,mapx,meqrec,meqrec%phr,inmap,ceq)
!    write(*,*)'MU(*) after  meq_sameset: ',ceq%cmuval(1),ceq%cmuval(2)
    if(gx%bmperr.ne.0) then
!       write(*,*)'MM Error calling meq_sameset from two_stoich',gx%bmperr
       goto 1000
    endif
! return the entered phase in mapx (maybe not needed?)
!    call get_state_var_value('X(O) ',value,phases,ceq)
    mapx=saveent
!    write(*,807)meqrec%fixph(1,1),meqrec%fixph(2,1),mapx,iadd,value
807 format('MM old fix phase/set: ',i3,i2,' entered: ',i3,', new fix: ',i3,&
         1pe12.4)
! restore status of new phase found at nodepoint as entered
    meqrec%phr(iadd)%phasestatus=PHENTERED
!    write(*,*)'Conditions for the invariant:'
!    call list_conditions(kou,ceq)
!    write(*,*)'Exiting two_stoich_same_comp'
! we must set this error code to return to mapping routines
! This means two stoichiometric phases stable an node point
    gx%bmperr=4364
1000 continue
! Make sure status of new phase found at nodepoint as set as entered
    meqrec%phr(iadd)%phasestatus=PHENTERED
    return
  end subroutine two_stoich_same_comp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

end MODULE liboceqplus

