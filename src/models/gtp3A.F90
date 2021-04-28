!
! gtp3A.F90 included in gtp3.F90
!
!****************************************************
! general subroutines for creating and handling elements, species, phases etc
! accessable externally

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     1. Initialization
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine init_gtp
!\begin{verbatim}
 subroutine init_gtp(intvar,dblvar)
! initiate the data structure
! create element and species record for electrons and vacancies
! the allocation of many arrays should be provided calling this routne
! intvar and dblvar will eventually be used for allocations and defaults
   implicit none
   integer intvar(*)
   double precision dblvar(*)
!\end{verbatim}
   character tpname*16,tpfun*80
   integer jl,ieq,ip,lrot,npid
!
   noofel=0; noofsp=0; noofph=0; nooftuples=0
!   write(*,3)'In init_gtp',maxel,maxsp,maxph
3  format(a,10i5)
! allocate records for elements
   allocate(ellista(-1:maxel))
   allocate(elements(-1:maxel))
! allocate records for species
   allocate(splista(maxsp))
   allocate(species(maxsp))
! allocate records for phases
   allocate(phlista(0:maxph))
   allocate(phases(0:maxph))
   phases=0
   allocate(phasetuple(0:2*maxph))
   do jl=1,2*maxph
      phasetuple%nextcs=0
   enddo
! phases(0) is refrence phase, evidently this index is never set
   phases(0)=0
!---------------------------
! create  special element /-
   ellista(-1)%symbol='/-'
   ellista(-1)%name='Electron'
   ellista(-1)%ref_state='Electron_gas'
   ellista(-1)%mass=zero
   ellista(-1)%h298_h0=zero
   ellista(-1)%s298=zero
   ellista(-1)%status=0
   ellista(-1)%alphaindex=-1
! The electron does not have any corresponing species
   ellista(-1)%splink=-1
   elements(-1)=-1
! create  special elements VA
   ellista(0)%symbol='VA'
   ellista(0)%name='Vacancy'
   ellista(0)%ref_state='Vacuum'
   ellista(0)%mass=zero
   ellista(0)%h298_h0=zero
   ellista(0)%s298=0.0D0
   ellista(0)%status=0
   ellista(0)%alphaindex=0
! splink set below
!    ellista(0)%splink=0
! allocate element link array
   allocate(splista(1)%ellinks(1))
   allocate(splista(1)%stoichiometry(1))
   splista(1)%symbol='VA'
   splista(1)%mass=zero
   splista(1)%charge=zero
   splista(1)%status=0
! set status bits that is is also an element and it is the vacancy
   splista(1)%status=ibset(splista(1)%status,SPEL)
   splista(1)%status=ibset(splista(1)%status,SPVA)
   splista(1)%alphaindex=1
   splista(1)%noofel=1
   splista(1)%ellinks(1)=0
   splista(1)%stoichiometry(1)=one
   elements(0)=0
   noofsp=1
   species(1)=1
! link from element Va to species Va
   ellista(0)%splink=1
!   write(*,3)'more allocate: ',maxrefs,maxprop,maxeq,maxtpf,maxsvfun
   allocate(bibrefs(maxrefs))
   allocate(propid(maxprop))
! first free data reference record (static)
   reffree=1
   addrecs=0
!---------------------------------------
   noofem=0
   noofint=0
   noofprop=0
!----------------------------------------
! initiate equilibrium record list
! dimension arrays for in first equilibrium record including phase_varres
   allocate(eqlista(maxeq))
   do jl=1,maxeq-1
! new 2019.12.17 zero status word!!
      eqlista(jl)%status=0
      eqlista(jl)%nexteq=jl+1
   enddo
   eqlista(maxeq)%nexteq=-1
   eqfree=1
! create first equilibrium record incl complist
   call enter_equilibrium('DEFAULT_EQUILIBRIUM ',ieq)
   if(gx%bmperr.ne.0) then
      write(*,*)' error in first enter_equilibrium',gx%bmperr
      goto 1000
   endif
   firsteq=>eqlista(1)
! nullify some pointers because of error entering first
   nullify(firsteq%lastcondition,firsteq%lastexperiment)
! set phase_varres free list in firsteq. These are always allocated together
   do jl=1,2*maxph-1
      firsteq%phase_varres(jl)%nextfree=jl+1
   enddo
! NOTE last phase_varres record used for copy in shiftcompsets
   firsteq%phase_varres(2*maxph)%nextfree=-1
! csfree and highcs are declared in gtp3.F90
   csfree=1; highcs=0
! convergence criteria for constituent fractions, 1e-6 works most often
! But one should take care to equilibrate fractions smaller than xconv!!!
   firsteq%xconv=1.0D-6
   firsteq%maxiter=500
   firsteq%gdconv(1)=4.0D-3
   firsteq%gdconv(2)=zero
! initiate tp functions
!   write(*,*)'init_gtp: initiate TP fuctions'
   jl=maxtpf
   call tpfun_init(jl,firsteq%eq_tpres)
!------------------------------------
! Property records define what can be used as "id" for parameters, the first
! must be G for the "chemical" part.  The others are connected to various
! additions or are simply properties that may depend on composition and is
! needed in other contexts, like mobilities, viscosities etc.
! create property id records for G
   npid=1
   propid(npid)%symbol='G '
   propid(npid)%note='Energy '
   propid(npid)%status=0
! This indicates if there are unkown or undefined MPI in a TDB file
   nundefmpi=0
!============================================================
! VERY IMPORTANT: The properties defined below must not be equal to state
! variables, or abbreviation of state variables.
! If so they cannot be listed and other errors may occur
! IMPORTANT any changes must be propagated to gtp3F: state_variable_val3 !!!
! after label 200!!
!
! ANY CHANGES HERE MUST BE MADE ALSO IN SUBROUTINE state_variable_val3
! IN THE RESULTS THE TYPE OF VARIABLE WILL BE STORED USING THE npid INDEX HERE
! OLD SAVE FILES MAY HAVE OTHER MEANING OF npid !!
!
!============================================================
! Mixed Curie/Neel Temperature, set bits that TC and BM cannot depend on T 2
   npid=npid+1
   propid(npid)%symbol='TC '
   propid(npid)%note='Combined Curie/Neel T' 
   propid(npid)%status=0
! TC cannot depend on T but on P
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Average Bohr magneton number 3
   npid=npid+1
   propid(npid)%symbol='BMAG '
   propid(npid)%note='Average Bohr magneton numb'
   propid(npid)%status=0
! BM cannot depend on either T or P ??
   propid(npid)%status=ibset(propid(npid)%status,IDNOTP)
!.......................................
! Specific Curie temperature 4
   npid=npid+1
   propid(npid)%symbol='CTA '
   propid(npid)%note='Curie temperature'
   propid(npid)%status=0
! CTA cannot depend on either T or P ??
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Specific Neel temperature 5
   npid=npid+1
   propid(npid)%symbol='NTA '
   propid(npid)%note='Neel temperature'
   propid(npid)%status=0
! NTA cannot depend on T but on P
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Individual Bohr magneton number 6 SPECIAL THIS HAS CONSTITUENT INDEX
   npid=npid+1
   propid(npid)%symbol='IBM '
   propid(npid)%note='Individual Bohr magneton numb'
!                     123456789.123456789.12345678-
   propid(npid)%status=0
! IBM cannot depend on either T or P and it is individual
   propid(npid)%status=ibset(propid(npid)%status,IDCONSUFFIX)
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Debye or Einstein temperature 7
   npid=npid+1
   propid(npid)%symbol='LNTH '
   propid(npid)%note='LN(Debye or Einstein temp)'
   propid(npid)%status=0
! LNTH cannot depend on T but on P
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!....................................... 8
! Molar volume at T=298.15, 1 bar
   npid=npid+1
   propid(npid)%symbol='V0 '
   propid(npid)%note='Volume at T0, P0 '
   propid(npid)%status=0
! Constant independent on temperature or pressure
   propid(npid)%status=ibset(propid(npid)%status,IDNOTP)
!....................................... 9
! Thermal expansion at 1 bar
   npid=npid+1
   propid(npid)%symbol='VA '
   propid(npid)%note='Thermal expansion '
   propid(npid)%status=0
! Not P dependent, only T dependent
   propid(npid)%status=ibset(propid(npid)%status,IDONLYT)
!....................................... 10
! Bulk modulus as function of T and P
   npid=npid+1
   propid(npid)%symbol='VB '
   propid(npid)%note='Bulk modulus '
   propid(npid)%status=0
!....................................... 11
! Extra volume parameter
   npid=npid+1
   propid(npid)%symbol='VC '
   propid(npid)%note='Alternative volume parameter'
   propid(npid)%status=0
!....................................... 12
! Diffusion volume parameter, suffix S on V create confusion with S as SER?
   npid=npid+1
   propid(npid)%symbol='VD '
   propid(npid)%note='Diffusion volume parameter '
   propid(npid)%status=0
!.......................................
! Activation energy of mobility 13
   npid=npid+1
   propid(npid)%symbol='MQ '
   propid(npid)%note='Mobility activation energy'
   propid(npid)%status=0
! MQ is specific for a constituent
   propid(npid)%status=ibset(propid(npid)%status,IDCONSUFFIX)
! indicate this parameter must not have wildcard constituents
   nowildcard(1)=npid
! in subroutine equilph1e we use the index of MQ to find mobility values
   mqindex=npid*100
!.......................................
! RT*ln(Frequency factor of mobility)  14
   npid=npid+1
   propid(npid)%symbol='MF '
   propid(npid)%note='RT*LN(mobility freq.fact.)'
   propid(npid)%status=0
! MF is specific for a constituent
   propid(npid)%status=ibset(propid(npid)%status,IDCONSUFFIX)
! indicate this parameter must not have wildcard constituents
   nowildcard(2)=npid
!.......................................
! Magnetic mobility factor 15
   npid=npid+1
   propid(npid)%symbol='MG '
   propid(npid)%note='Magnetic mobility factor'
   propid(npid)%status=0
! MG is specific for a constituent
   propid(npid)%status=ibset(propid(npid)%status,IDCONSUFFIX)
! indicate this parameter must not have wildcard constituents
   nowildcard(3)=npid
!....................................... 13 fd 11
! Liquid two-state model     16
   npid=npid+1
   propid(npid)%symbol='G2   '
   propid(npid)%note='Liquid two state parameter'
   propid(npid)%status=0
!.......................................
! Smooth unit step function (or second Einstein function) 17
   npid=npid+1
   propid(npid)%symbol='THT2 '
   propid(npid)%note='LN(Smooth step function Tcrit)'
   propid(npid)%status=0
! THT2 cannot depend on T but on P
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Second Einstein delta CP 18
   npid=npid+1
   propid(npid)%symbol='DCP2 '
   propid(npid)%note='Smooth step function increm.'
   propid(npid)%status=0
! DXP2 cannot depend on T but on P
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Lattice parameter in direction X 19
   npid=npid+1
   propid(npid)%symbol='LPX '
   propid(npid)%note='Lattice param X axis'
   propid(npid)%status=0
! lattice parameters may depend on T and P
!.......................................
! Lattice parameter in direction Y 20
   npid=npid+1
   propid(npid)%symbol='LPY '
   propid(npid)%note='Lattice param Y axis'
   propid(npid)%status=0
! lattice parameters may depend on T and P
!.......................................
! Lattice parameter in direction Z 21
   npid=npid+1
   propid(npid)%symbol='LPZ '
   propid(npid)%note='Lattice param Z axis'
   propid(npid)%status=0
! lattice parameters may depend on T and P
!.......................................
! This is an angle for non-cubic lattices 22
   npid=npid+1
   propid(npid)%symbol='LPTH '
   propid(npid)%note='Lattice angle TH'
   propid(npid)%status=0
! Angle may depend on T and P 
!.......................................
! This is an elastic "constant" 23
   npid=npid+1
   propid(npid)%symbol='EC11 '
   propid(npid)%note='Elastic const C11'
   propid(npid)%status=0
! The elastic constant may depend on T and P
!.......................................
! This is another elastic "constant" 24
   npid=npid+1
   propid(npid)%symbol='EC12 '
   propid(npid)%note='Elastic const C12'
   propid(npid)%status=0
! The elastic constant may depend on T and P
!.......................................
! This is yet another elastic "constant" 25
   npid=npid+1
   propid(npid)%symbol='EC44 '
   propid(npid)%note='Elastic const C44'
   propid(npid)%status=0
! The elastic constant may depend on T and P
!.......................................
! VERY SPECIAL this model parameter identifier has no addition
! thus no check for addition in enter_parameter subroutine (gtp3B.F90)
! UNIQUAC interaction parameter 26
! IF THIS IS CHANGED TO ANOTHER NUMBER CHANGES NEEDED IN GTP3B: mpiwarning
   npid=npid+1
   propid(npid)%symbol='UQT '
   propid(npid)%note='UNIQUAC residual parameter '
   propid(npid)%status=0
! UQT is specific for a constituent, 2600+constituent index
   propid(npid)%status=ibset(propid(npid)%status,IDCONSUFFIX)
!.......................................
! Electrical resistivity 27
   npid=npid+1
   propid(npid)%symbol='RHO '
   propid(npid)%note='Electric resistivity'
   propid(npid)%status=0
!....................................... f.d. 18 now 28
! Viscosity 28
   npid=npid+1
   propid(npid)%symbol='VISC '
   propid(npid)%note='Viscosity'
   propid(npid)%status=0
!....................................... 
! Thermal conductivity as function of T and P: 29
   npid=npid+1
   propid(npid)%symbol='LAMB '
   propid(npid)%note='Thermal conductivity '
   propid(npid)%status=0
!.......................................
! From MatCalc databases 30
   npid=npid+1
   propid(npid)%symbol='HMVA '
   propid(npid)%note='Enthalpy of vacancy form. '
   propid(npid)%status=0
! this parameter does not depend on T ??
!   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Schottky anomaly T 31
   npid=npid+1
   propid(npid)%symbol='TSCH '
   propid(npid)%note='Schottky anomaly T '
   propid(npid)%status=0
! this parameter does not depend on T ??
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Schottky anomaly CP/R 32
   npid=npid+1
   propid(npid)%symbol='CSCH '
   propid(npid)%note='Schottky anomaly Cp/R. '
   propid(npid)%status=0
! this parameter does not depend on T ??
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Modified Quasichemical model coordination factor 33
   npid=npid+1
   propid(npid)%symbol='QCZ'
   propid(npid)%note='MQMQA cluster coord factor'
   propid(npid)%status=0
! This parameter does not depend on T and P
   propid(npid)%status=ibset(propid(npid)%status,IDNOTP)
!.......................................
!.......................................
! This IF statement should be at the last parameter identifier, maxprop=50 ?
   if(npid.gt.maxprop) then
      write(*,*)'Too many parameter identifiers, increase maxprop'
      gx%bmperr=4250; goto 1000
   endif
!   write(*,*)'3A number of model parameter identifiers: ',npid
! IMPORTANT any changes must be propagated to gtp3F: state_variable_val3 !!!
!.......................................
! IMPORTRANT: When adding more parameter identifiers one should NEVER
! NEVER USE A NAME ENDING IN D as that will be taken as a "disordered" part
! The number of defined properties, should be less than maxprop (=50?)
! IMPORTANT: In the addition records one must use the parameter identifier
! to extract the calculated composition dependent values
! IMPORTANT: in gtp3F new variables must be added to be able to list/plot them
   ndefprop=npid
!-------------------------------------------------
!CCI : GSVIRTUAL enables to do calculation with virtual elements
   globaldata%status=ibclr(globaldata%status,GSVIRTUAL)
!CCI
! globaldata record; set gas constant mm
   globaldata%status=0
! set beginner, no data, no phase, no equilibrium calculated
   globaldata%status=ibset(globaldata%status,GSBEG)
!   globaldata%status=ibset(globaldata%status,GSADV)
   globaldata%status=ibset(globaldata%status,GSNODATA)
   globaldata%status=ibset(globaldata%status,GSNOPHASE)
   firsteq%status=ibset(firsteq%status,EQNOEQCAL)
! set that dense grid is used by default
!   globaldata%status=ibset(globaldata%status,GSXGRID)
! set gas constant and some default values
   globaldata%name='current'
   globaldata%rgas=8.31451D0
! more recent value not used as all TDB file used the old
!   globaldata%rgas=8.3144621D0
! old value of gas constant
   globaldata%rgasuser=8.31451D0
   globaldata%pnorm=one
! zero sysparam and sysreal
   globaldata%sysparam=0
   globaldata%sysreal=zero
!   write(*,*)'init_gtp: enter R and RTLNP'
! enter R as TP function
   tpname='R'
!   write(tpfun,777)' 10 8.31451; 20000 N '
!777 format(a)
!   call enter_tpfun(tpname,tpfun,lrot,.FALSE.)
   call store_tpconstant(tpname,globaldata%rgas)
   if(gx%bmperr.ne.0) goto 1000
   tpname='RTLNP'
   tpfun=' 1 R*T*LN(1.0D-5*P); 20000 N '
!   call store_tpfun(tpname,tpfun,lrot,.FALSE.)
   call store_tpfun(tpname,tpfun,lrot,-1)
   if(gx%bmperr.ne.0) goto 1000
! default minimum fraction
   bmpymin=ymind
! putfun error code .... should use buperr at least
   pfnerr=0
!------------------------------------
! allocate array for state variable function
!   write(*,*)'init_gtp: allocate array for state variable functions'
   allocate(svflista(maxsvfun))
! number of state variable functions
   nsvfun=0
! zero the array with equilibrium index for functions, not used aywhere??
!   pflocal=0
! enter some useful state variable function
   tpfun=' R=8.31451;'
   ip=1
!   write(*,*)'init_gtp: entering function R'
   call enter_svfun(tpfun,ip,firsteq)
! mark it cannot be amended
   svflista(1)%status=ibset(svflista(1)%status,SVNOAM)
! mark it is a constant
   svflista(1)%status=ibset(svflista(1)%status,SVCONST)
!   if(gx%bmperr.ne.0) then
!      write(*,*)'Error entering R',gx%bmperr
!      goto 1000
!   endif
!   write(*,*)'Entered symbol R'
   tpfun=' RT=R*T;'
   ip=1
!   write(*,*)'init_gtp: entering function RT'
   call enter_svfun(tpfun,ip,firsteq)
! mark it cannot be amended
   svflista(2)%status=ibset(svflista(2)%status,SVNOAM)
!   if(gx%bmperr.ne.0) then
!      write(*,*)'Error entering symbol RT'
!      goto 1000
!   endif
!   write(*,*)'Entered symbol RT'
   tpfun=' T_C=T-273.15;'
   ip=1
!   write(*,*)'init_gtp: entering function T_C'
   call enter_svfun(tpfun,ip,firsteq)
! mark it cannot be amended
   svflista(3)%status=ibset(svflista(3)%status,SVNOAM)
!   if(gx%bmperr.ne.0) then
!      write(*,*)'Error entering symbol T_C'
!      goto 1000
!   endif
! we evaluate all symbols to avoid some problems ... no output
!  call meq_evaluate_all_svfun(-1,ceq) cannot be used as it is in minimizer ...
   call evaluate_all_svfun_old(-1,firsteq)
! set working directory (decleared in metlib, used now and again ...)
   call getcwd(workingdir)
! assessment initiallizing
!   write(*,*)'3A Initiallizing firstash', firstash is a pointer ...
   call assessmenthead(firstash)
!   firstash%status=0
!   write(*,*)'firstash allocated: ',firstash%status
!   nullify(firstash%prevash)
!   nullify(firstash%nextash)
! create the beginnings of a circular list
   firstash%nextash=>firstash
   firstash%prevash=>firstash
! set that dense grid used by default
!   globaldata%status=ibset(globaldata%status,GSXGRID)
! removed line above as that caused crash in parallel2 WHY????
! finished initiating
1000 continue
!   write(*,*)'exit from init_gtp'
   return
 END subroutine init_gtp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine assessmenthead
!\begin{verbatim}
 subroutine assessmenthead(ash)
! create an assessment head record and do more (later)
   type(gtp_assessmenthead), pointer :: ash
!   type(gtp_assessmenthead), allocatable :: ash
!\end{verbatim}
! it is not good to allocate a pointer, memory loss!!
   allocate(ash)
   ash%status=0
   return
 end subroutine assessmenthead

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     2. Number of things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function noel
!\begin{verbatim}
 integer function noel()
! number of elements because noofel is private
! should take care if elements are suspended
!\end{verbatim} %+
   noel=noofel
 end function noel

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function nosp
!\begin{verbatim} %-
 integer function nosp()
! number of species because noofsp is private
! should take care if species are suspended
!\end{verbatim} %+
   nosp=noofsp
 end function nosp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function noph
!\begin{verbatim} %-
 integer function noph()
! number of phases because noofph is private
! should take care if phases are hidden
!\end{verbatim} %+
   noph=noofph
 end function noph

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function noofcs
!\begin{verbatim} %-
 integer function noofcs(iph)
! returns the number of compositions sets for phase iph
   implicit none
   integer iph
!\end{verbatim} %+
   if(iph.le.0 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   endif
   noofcs=phlista(phases(iph))%noofcs
1000 continue
   return
 end function noofcs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function noconst
!\begin{verbatim} %-
 integer function noconst(iph,ics,ceq)
! number of constituents for iph (include single constituents on a sublattice)
! It tests if a constituent is suspended which can be different in each ics.
   implicit none
   integer iph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,lokcs,noc,jl
   if(iph.gt.0 .and. iph.le.noofph) then
      lokph=phases(iph)
      if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
!         write(*,*)'noconst 1 error 4072'
         gx%bmperr=4072; goto 1000
      elseif(ics.eq.0) then
         ics=1
      endif
      lokcs=phlista(lokph)%linktocs(ics)
      if(btest(ceq%phase_varres(lokcs)%status2,CSCONSUS)) then
! some constituents suspended ?? NOT POSSIBLE as not implemented 190923 !!!
!         write(*,*)'3A suspended constituents!!',iph,ics,lokph,lokcs
         write(*,*)'3A suspended constituents: ',phlista(lokph)%tnooffr,&
              allocated(ceq%phase_varres(lokcs)%constat)
         if(.not.allocated(ceq%phase_varres(lokcs)%constat)) then
            noconst=phlista(lokph)%tnooffr
         else   
            noc=phlista(lokph)%tnooffr
            do jl=1,phlista(lokph)%tnooffr
               if(btest(ceq%phase_varres(lokcs)%constat(jl),CONSUS)) then
                  noc=noc-1
               endif
            enddo
            noconst=noc
         endif
      else
         noconst=phlista(lokph)%tnooffr
      endif
   else
      gx%bmperr=4050
   endif
1000 continue
   return
 end function noconst

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function nooftup
!\begin{verbatim} %-
 integer function nooftup()
! number of phase tuples
!\end{verbatim} %+
   implicit none
   nooftup=nooftuples
   return
 end function nooftup

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function noofphasetuples
!\begin{verbatim} %-
 integer function noofphasetuples_old()
! number of phase tuples REDUNDANT !!
!\end{verbatim}
   noofphasetuples_old=nooftuples
   return
 end function noofphasetuples_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function nosvf
!\begin{verbatim}
 integer function nosvf()
! number of state variable functions
!\end{verbatim}
   implicit none
   nosvf=nsvfun
   return
 end function nosvf

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable integer function noeq
!\begin{verbatim}
 integer function noeq()
! returns the number of equilibria entered
!\end{verbatim}
   implicit none
   noeq=eqfree-1
1000 continue
   return
 end function noeq

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable integer function nosusphcs
!\begin{verbatim}
 integer function nonsusphcs(ceq)
! returns the total number of unhidden phases+composition sets
! in the system.  Used for dimensioning work arrays and in loops
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer tphic,iph,ics,lokph
   double precision xxx
   tphic=0
   do iph=1,noofph
      lokph=phases(iph)
      ics=1
      if(test_phase_status(iph,ics,xxx,ceq).ne.PHHIDDEN) then
! phase is not hidden
         do ics=1,phlista(lokph)%noofcs
!         if(test_phase_status(iph,ics,xxx,ceq).eq.4) goto 400
            if(test_phase_status(iph,ics,xxx,ceq).ne.PHSUS) then
               tphic=tphic+1
            endif
! composition set not suspended
!         tphic=tphic+phlista(lokph)%noofcs
         enddo
      endif
   enddo
1000 continue
!   write(*,*)'25 A nonsusphcs: ',tphic
   nonsusphcs=tphic
   return
 end function nonsusphcs
 
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     3. Find things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_element_by_name
!\begin{verbatim}
 subroutine find_element_by_name(name,iel)
! find an element index by its name, exact fit required
   implicit none
   character name*(*)
   integer iel
!\end{verbatim} %+
   integer lokel
   character symbol*2
   symbol=name
   call capson(symbol)
   do lokel=-1,noofel
!       write(*,*)'find_element 1: ',lokel,symbol,' ',ellista(lokel)%symbol
      if(symbol.eq.ellista(lokel)%symbol) then
         iel=ellista(lokel)%alphaindex
         goto 1000
      endif
   enddo
   iel=-100
   gx%bmperr=4042
1000 continue
   return
 end subroutine find_element_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_component_by_name
!\begin{verbatim} %-
 subroutine find_component_by_name(name,icomp,ceq)
! BEWARE: one may in the future have different components in different
! equilibria. components are a subset of the species
   implicit none
   character*(*) name
   integer icomp
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer loksp
   call find_species_record_noabbr(name,loksp)
   if(gx%bmperr.ne.0) then
      gx%bmperr=4052; goto 1000
   endif
! check that species actually is component
   do icomp=1,noofel
      if(ceq%complist(icomp)%splink.eq.loksp) goto 1000
   enddo
   gx%bmperr=4052
1000 continue
   return
 end subroutine find_component_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_species_by_name
!\begin{verbatim} %-
 subroutine find_species_by_name(name,isp)
! locates a species index from its name, unique abbreviation 
! or exact match needed
   implicit none
   character name*(*)
   integer isp
!\end{verbatim} %+
   character symbol*24
   integer loksp,lensym
   logical exact
   exact=.FALSE.
   symbol=name
   call capson(symbol)
   isp=0
   do loksp=1,noofsp
!      write(*,*)'3A find species 2: ',symbol,splista(loksp)%symbol,loksp
      if(compare_abbrev(symbol,splista(loksp)%symbol)) then
         if(isp.eq.0) then
            isp=splista(loksp)%alphaindex
            lensym=len_trim(splista(loksp)%symbol)
!            write(*,*)'3A abbr match: ',lensym,' <',symbol(1:lensym),'><',&
!                 splista(loksp)%symbol(1:lensym+1),'>'
            if(symbol(1:lensym+1).eq.splista(loksp)%symbol(1:lensym+1)) then
!               write(*,*)'3A exact match with species name'
               exact=.TRUE.
               goto 1000
            endif
         else
! abbreviation is not unique
            isp=0
            exit
         endif
      endif
   enddo
   if(isp.eq.0) then
!      write(*,*)'in find_species_by_name'
      gx%bmperr=4051
      loksp=0
   endif
1000 continue
   return
 end subroutine find_species_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_species_record
!\begin{verbatim} %-
 subroutine find_species_record(name,loksp)
! locates a species record allowing abbreviations
   implicit none
   character name*(*)
   integer loksp
!\end{verbatim} %+
   character symbol*24
   integer isp,lensp
   logical exact
   exact=.FALSE.
   symbol=name
   isp=0
   call capson(symbol)
   do loksp=1,noofsp
!      write(*,17)'3A find species: ',loksp,splista(loksp)%symbol,name
17    format(a,i3,' "',a,'" "',a,'"')
      if(compare_abbrev(symbol,splista(loksp)%symbol)) then
         if(isp.eq.0) then
            isp=loksp
! it would be enough to compare lengths of species ...
            lensp=len_trim(splista(loksp)%symbol)
            if(symbol(1:lensp+1).eq.splista(loksp)%symbol(1:lensp+1)) then
!               write(*,*)'3A exact match'
               exact=.TRUE.
               goto 1000
            endif
         else
! ambiguous species name but we may find an exact later ...
            isp=-1
         endif
      endif
   enddo
   if(isp.le.0) then
!      write(*,*)'Error in find_species_record "',name,'"'
      gx%bmperr=4051
      loksp=0
   else
      loksp=isp
   endif
1000 continue
   return
 end subroutine find_species_record

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_species_record_noabbr
!\begin{verbatim} %-
 subroutine find_species_record_noabbr(name,loksp)
! locates a species record no abbreviations allowed
   implicit none
   character name*(*)
   integer loksp
!\end{verbatim} %+
   character symbol*24
   symbol=name
   call capson(symbol)
   do loksp=1,noofsp
!       write(*,17)'find species 17B: ',loksp,splista(loksp)%symbol,name
!17     format(a,i3,' "',a,'" "',a,'"')
      if(symbol.eq.splista(loksp)%symbol) goto 1000
   enddo
!   write(*,*)'Error in find_species_record_noabbr "',name,'"'
   gx%bmperr=4051
   loksp=0
1000 continue
   return
 end subroutine find_species_record_noabbr

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_species_record_exact
!\begin{verbatim} %-
 subroutine find_species_record_exact(name,loksp)
! locates a species record, exact match needed
! for parameters, V must not be accepted as abbreviation of VA or C for CR
   implicit none
   integer loksp
   character name*(*)
!\end{verbatim}
   integer quad
   character symbol*24
   symbol=name
   call capson(symbol)
! special for quadrupoles ... they can have a trailing -Qij which may be
! different each time ...
   quad=index(symbol,'-Q')
   if(quad.gt.0) then
      quad=quad-1
   else
      quad=0
   endif
   do loksp=1,noofsp
!       write(*,17)'find species 17: ',loksp,splista(loksp)%symbol,name
!17     format(a,i3,' "',a,'" "',a,'"')
      if(quad.gt.0) then
         if(symbol(1:quad).eq.splista(loksp)%symbol(1:quad)) goto 1000
      else
! problem that V was not read from database ...
         if(symbol.eq.splista(loksp)%symbol) goto 1000
      endif
   enddo
! This message cannot be written as it is used when reading a TDB file ...
!   write(kou,*)'Exact match to species name requited'
   gx%bmperr=4051
   loksp=0
1000 continue
   return
 end subroutine find_species_record_exact

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_phasetuple_by_name
!\begin{verbatim}
 subroutine find_phasetuple_by_name(name,phcsx)
! finds a phase with name "name", returns phase tuple index
! handles composition sets either with prefix/suffix or #digit
! When no pre/suffix nor # always return first composition set
   implicit none
   character name*(*)
   integer phcsx
!\end{verbatim} %+
   integer iph,ics
   iph=0
   ics=0
   phcsx=0
   call find_phasex_by_name(name,phcsx,iph,ics)
1000 continue
   return
 end subroutine find_phasetuple_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_phase_by_name
!\begin{verbatim} %-
 subroutine find_phase_by_name(name,iph,ics)
! finds a phase with name "name", returns address of phase, first fit accepted
! handles composition sets either with prefix/suffix or #digit
! When no pre/suffix nor # always return first composition set
   implicit none
   character name*(*)
   integer iph,ics
!\end{verbatim} %+
   integer phcsx
   phcsx=0
   call find_phasex_by_name(name,phcsx,iph,ics)
1000 continue
   return
 end subroutine find_phase_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function find_phasetuple_by_indices
!\begin{verbatim} %-
 integer function find_phasetuple_by_indices(iph,ics)
! subroutine find_phasetuple_by_indices(iph,ics)
! find phase tuple index given phase index and composition set number
   integer iph,ics
!\end{verbatim} %+
   integer ij
   ij=iph
   if(ij.gt.0 .and. ij.le.nooftuples) then
      do while(ij.gt.0)
         if(ics.eq.phasetuple(ij)%compset) then
            find_phasetuple_by_indices=ij
            goto 1000
         else
            ij=phasetuple(ij)%nextcs
         endif
      enddo
   endif
   write(*,*)'Wrong arguments to find_phasetuple_by_indices: ',iph,ics,ij
   gx%bmperr=4073
1000 continue
 return
end function find_phasetuple_by_indices

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_phasex_by_name
!\begin{verbatim} %-
 subroutine find_phasex_by_name(name,phcsx,iph,zcs)
! finds a phase with name "name", returns index and tuplet of phase.
! All phases checked and error return if name is ambiguous
! handles composition sets either with prefix/suffix or #digit or both
! if no # check all composition sets for prefix/suffix
! special if phcsx = -1 and there are several composition sets then
! zcs is set to -(number of composition sets).  Used when changing status
! phcsx, iph and zcs are values to return!
   implicit none
   character name*(*)
   integer phcsx,iph,zcs
!\end{verbatim} %+
   character name1*36,csname*36,name2*24,name3*24
   TYPE(gtp_phase_varres), pointer :: csrec
   integer kp,kcs,lokph,jcs,lokcs,first1,fcs,lcs,ics,lenam,allsets
! set ics to an illegal value
   ics=-1
   allsets=phcsx
! convert to upper case locally
   name1=name
   call capson(name1)
! composition set as #digit
   kp=index(name1,'#')
   if(kp.gt.0) then
      ics=ichar(name1(kp+1:kp+1))-ichar('0')
! negative ics should give error, 0 should be the same as 1
      if(ics.eq.0) ics=1
      if(ics.lt.1 .or. ics.gt.9) then
         gx%bmperr=4093; goto 1000
      endif
      allsets=ics
      name1(kp:)=' '
      kcs=ics
   else
      ics=1
      kcs=0
   endif
!   write(*,17)trim(name),ics,kcs,kp,noofph
17 format('3A find_phase 3: ',a,2x,10i4)
   first1=0
   loop1: do lokph=1,noofph
      if(kcs.eq.0) then
! no composition set specified explicitly, all sets must be checked
         fcs=2; lcs=phlista(lokph)%noofcs
!      elseif(kcs.eq.1) then
!         fcs=1; lcs=1
      elseif(kcs.le.phlista(lokph)%noofcs) then
! we shoud check pre and suffix ...
         fcs=max(2,kcs); lcs=kcs
      else
! this phase does not have a composition set kcs
         cycle loop1
      endif
      name2=phlista(lokph)%name
      if(kcs.le.1) then
         if(compare_abbrev(name1,name2)) then
            if(first1.eq.0) then
               first1=lokph
               if(len_trim(name1).eq.len_trim(name2)) then
! exact match, we already know there is a composition set
!                  write(*,*)'3A exact match',name1(1:len_trim(name1)),lokph
                  goto 300
               endif
            else
! another phase with same abbreviation, phase name is ambiguous
               gx%bmperr=4121
               goto 1000
            endif
         endif
      endif
! if composition set specified check only that set, otherwise all from 2
!      write(*,*)'3A first1: ',first1,fcs,lcs
      loop2: do jcs=fcs,lcs
         lokcs=phlista(lokph)%linktocs(jcs)
         csrec=>firsteq%phase_varres(lokcs)
         kp=len_trim(csrec%prefix)
         if(kp.gt.0) then
            csname=csrec%prefix(1:kp)//'_'//name2
         else
            csname=name2
         endif
         kp=len_trim(csrec%suffix)
         if(kp.gt.0) csname=csname(1:len_trim(csname))//'_'//&
              csrec%suffix(1:kp)
!         write(*,244)ics,kcs,jcs,kp,fcs,lcs,first1,name1(1:len_trim(name1)),&
!              csname(1:len_trim(csname))
244      format('3A: find_phase: ',7i3,'<',a,'>=?=<',a,'>')
         if(compare_abbrev(name1,csname)) then
            if(first1.eq.lokph) then
! match already with first composition set, that is OK
               cycle loop2
            elseif(first1.eq.0) then
               first1=lokph
               ics=jcs
               allsets=ics
            else
! ambiguous phase name
               gx%bmperr=4121; goto 1000
            endif
         elseif(kcs.gt.1) then
! No mach with phase name including pre/suffix but if user has specified #
! accept also match with original name without pre/suffix
            if(compare_abbrev(name1,name2)) then
               if(first1.eq.0) then
                  first1=lokph
                  ics=jcs
               else
! another phase with same abbreviation, phase name is ambiguous
                  gx%bmperr=4121
                  goto 1000
               endif
            endif
         endif
      enddo loop2
   enddo loop1
!   write(*,*)'3A first1: ',first1
   if(first1.eq.0) then
! no phase found
      gx%bmperr=4050
      goto 1000
   endif
300 continue
! first1 is lokph for phase
   iph=phlista(first1)%alphaindex
   if(allsets.eq.-1) then
! special to set status: return -(number of composition sets) in zcs if >1
! DO NOT CHANGE PHCSX
      lcs=phlista(first1)%noofcs
      if(lcs.gt.1) then
         ics=-lcs; zcs=-lcs
      else
         ics=1; zcs=1
      endif
   else
! ics set above, return it in zcs
      zcs=ics
      phcsx=firsteq%phase_varres(phlista(first1)%linktocs(ics))%phtupx
   endif
   gx%bmperr=0
1000 continue
   return
1100 continue
   gx%bmperr=4073
   goto 1000
 END subroutine find_phasex_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_phase_by_name_exact
!\begin{verbatim} %-
 subroutine find_phase_by_name_exact(name,iph,ics)
! finds a phase with name "name", returns address of phase. exact match req.
! handles composition sets either with prefix/suffix or #digit
! no pre/suffix nor # gives first composition set
   implicit none
   character name*(*)
   integer iph,ics
!\end{verbatim}
   character name1*36,csname*36,name2*24
   TYPE(gtp_phase_varres), pointer :: csrec
   integer kp,kcs,iphfound,lokph,jcs,lokcs
! convert to upper case locally
   name1=name
   call capson(name1)
! composition set as #digit
   kp=index(name1,'#')
   if(kp.gt.0) then
      ics=ichar(name1(kp+1:kp+1))-ichar('0')
! negative ics should give error, 0 should be the same as 1
      if(ics.eq.0) ics=1
      if(ics.lt.1 .or. ics.gt.9) then
         gx%bmperr=4093; goto 1000
      endif
      name1(kp:)=' '
      kcs=ics
   else
      ics=1
      kcs=0
   endif
!    write(*,17)ics,kcs
17  format('find_phase 3: ',2i4)
!    write(*,11)'fpbne 1: ',name,noofph
11  format(a,a,'; ',2i3)
   iphfound=0
   loop1: do lokph=1,noofph
      name2=phlista(lokph)%name
!       write(*,*)'find_phase 2: ',name1,name2
      if(compare_abbrev(name1,name2)) then
         if(ics.le.phlista(lokph)%noofcs) then
! possible phase, if exact match no more checks
            if(trim(name1).eq.trim(name2)) then
               iphfound=lokph
               goto 300
            endif
            if(iphfound.ne.0) then
               if(trim(name1).eq.trim(name2)) then
                  iphfound=lokph
                  goto 300
               else
                  iphfound=-lokph
               endif
            else
               iphfound=lokph
            endif
         else
!            write(*,18)ics,phlista(lokph)%noofcs
18  format('find_phase 4: ',2i4)
            gx%bmperr=4072; goto 1000
         endif
      endif
   enddo loop1
!    write(*,*)'find_phase ',iphfound
   if(iphfound.lt.0) then
! several phases found
      gx%bmperr=4121; goto 1000
   elseif(iphfound.le.0) then
! no phase found
      gx%bmperr=4050; goto 1000
   else
      lokph=iphfound
      goto 300
   endif
! if there are composition sets check name including prefix/suffix
   write(*,*)'find_phase 5: ',lokph,phlista(lokph)%noofcs
   do jcs=2,phlista(lokph)%noofcs
      lokcs=phlista(lokph)%linktocs(jcs)
      csrec=>firsteq%phase_varres(lokcs)
      kp=len_trim(csrec%prefix)
      if(kp.gt.0) then
         csname=csrec%prefix(1:kp)//'_'//name2
      else
         csname=name2
      endif
      kp=len_trim(csrec%suffix)
      if(kp.gt.0) csname=csname(1:len_trim(csname))//'_'//&
           csrec%suffix(1:kp)
      if(compare_abbrev(name1,csname)) then
! if user has provided both #<digit> and pre/suffix these must be consistent
         if(kcs.gt.0 .and. kcs.ne.jcs) goto 1100
         ics=jcs
         goto 300
      endif
   enddo
250 continue
! no phase with this name
   gx%bmperr=4050
   goto 1000
300 continue
   iph=phlista(lokph)%alphaindex
   gx%bmperr=0
1000 continue
   return
1100 continue
! composition set index and pre/suffix does not match
   gx%bmperr=4073
   goto 1000
 END subroutine find_phase_by_name_exact

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine findeq
!\begin{verbatim}
 subroutine findeq(name,ieq)
! finds the equilibrium with name "name" and returns its index
! ieq should be the current equilibrium
   implicit none
   character name*(*)
   integer ieq
!\end{verbatim} %+
   character name2*64
   integer jeq
   name2=name
   call capson(name2)
! Accept abbreviations of PREVIOUS and FIRST (DEFAULT is the same as the first)
   jeq=0
   if(compare_abbrev(name2,'PREVIOUS ')) then
      jeq=max(1,ieq-1); goto 200
   elseif(compare_abbrev(name2,'FIRST ')) then
      jeq=1; goto 200
   elseif(compare_abbrev(name2,'DEFAULT ')) then
      jeq=1; goto 200
!   elseif(compare_abbrev(name2,'LAST ')) then
!      jeq=1; goto 200
   endif
100 jeq=jeq+1
!    write(*,*)'findeq 2: ',jeq,name2
   if(jeq.ge.eqfree) then
      gx%bmperr=4124
      goto 1000
   endif
!    write(*,*)'findeq 3: ',jeq,eqlista(jeq)%eqname
   if(.not.compare_abbrev(name2,eqlista(jeq)%eqname)) goto 100
!    if(eqlista(jeq)%eqname.ne.name2) goto 100
200 continue
   ieq=jeq
1000 continue
 end subroutine findeq

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine selecteq
!\begin{verbatim} %-
 subroutine selecteq(ieq,ceq)
! checks if equilibrium ieq exists and if so set it as current
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer ieq
!\end{verbatim}
   if(ieq.lt.0 .or. ieq.ge.eqfree) then
      gx%bmperr=4124
      goto 1000
   endif
   ceq=>eqlista(ieq)
1000 continue
 end subroutine selecteq

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     4. Get things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_phase_record
!\begin{verbatim}
 subroutine get_phase_record(iph,lokph)
! given phase index iph this returns the phase location lokph
   implicit none
   integer iph,lokph
!\end{verbatim} %+
   if(iph.lt.1 .or. iph.gt.noofph) then
!      write(*,*)'gpr: ',iph,noofph
      gx%bmperr=4050
   else
      lokph=phases(iph)
   endif
   return
 end subroutine get_phase_record

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_phase_variance
!\begin{verbatim} %-
 subroutine get_phase_variance(iph,nv)
! returns the number of independent variable fractions in phase iph
   implicit none
   integer iph,nv
!\end{verbatim} %+
   integer lokph
   call get_phase_record(iph,lokph)
   nv=phlista(lokph)%tnooffr-phlista(lokph)%noofsubl
   return
 end subroutine get_phase_variance

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_constituent_location
!\begin{verbatim} %-
 subroutine get_constituent_location(lokph,cno,loksp)
! returns the location of the species record of a constituent
! requred for ionic liquids as phlista is private
   implicit none
   integer lokph,loksp,cno
!\end{verbatim} %+
   loksp=phlista(lokph)%constitlist(cno)
   return
 end subroutine get_constituent_location

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_phase_compset
!\begin{verbatim} %-
 subroutine get_phase_compset(iph,ics,lokph,lokcs)
! Given iph and ics the phase and composition set locations are returned
! Checks that ics and ics are not outside bounds.
   implicit none
   integer iph,ics,lokph,lokcs
!\end{verbatim} %+
   if(iph.le.0 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   endif
   lokph=phases(iph)
! find composition set
   if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
      gx%bmperr=4072; goto 1000
   elseif(ics.eq.0) then
      ics=1
   endif
   lokcs=phlista(lokph)%linktocs(ics)
1000 continue
   return
 end subroutine get_phase_compset

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine find_constituent
!\begin{verbatim}
 subroutine find_constituent(iph,spname,mass,icon)
! find the constituent "spname" of a phase. spname can have a sublattice #digit
! Return the index of the constituent in icon.  Additionally the mass
! of the species is returned.
   implicit none
   character*(*) spname
   double precision mass
   integer iph,icon
!\end{verbatim}
! BUG found, asking for a constituent N it returned the constituent NB !!!
! Must search for exact match!!!
   character spname1*24
   integer lokph,kp,ll,kk,loksp,ls,first,jabbr
   lokph=phases(iph)
   kp=index(spname,'#')
   if(kp.gt.0) then
      ls=ichar(spname(kp+1:kp+1))-ichar('0')
      spname1=spname(1:kp-1)
   else
      ls=0
      spname1=spname
   endif
   call capson(spname1)
   icon=0
   jabbr=0
   first=0
   lloop: do ll=1,phlista(lokph)%noofsubl
      sploop: do kk=1,phlista(lokph)%nooffr(ll)
         icon=icon+1
         if(ls.eq.0 .or. ls.eq.ll) then
            loksp=phlista(lokph)%constitlist(icon)
! constituent icon is the requested one ??
!            write(*,55)ll,kk,icon,trim(spname1),trim(splista(loksp)%symbol)
55          format('find_const 7: ',3i3,1x,a,2x,a)
            if(compare_abbrev(spname1,splista(loksp)%symbol)) then
!               write(*,*)'3A abbreviation OK: ',trim(spname1),'?',&
!                    trim(splista(loksp)%symbol),icon
               if(trim(spname1).eq.trim(splista(loksp)%symbol)) then
! if exact match accept
                  first=loksp; goto 90
               elseif(first.eq.0) then
! constituent name is an abbreviation, if only one accept
                  first=loksp
                  jabbr=icon
               else
                  gx%bmperr=4121
                  goto 1000
               endif
            endif
         endif
!         write(*,*)'3A current: ',icon,first,loksp
      enddo sploop
   enddo lloop
90 continue
!   write(*,*)'3A current: ',icon,first,loksp,' "',trim(spname1),'"'
   if(first.eq.0) then
! no such constituent
      gx%bmperr=4096
   else
      if(jabbr.gt.0) then
! accept unique abbreviation
!         write(*,*)'3A abbreviation: ',icon,jabbr,loksp
         icon=jabbr
      endif
      mass=splista(first)%mass
   endif
1000 continue
   return
 end subroutine find_constituent

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_constituent_name
!\begin{verbatim} %-
 subroutine get_constituent_name(iph,iseq,spname,mass)
! find the constituent with sequential index iseq in phase iph
! return name in "spname" and mass in mass
   implicit none
   character*(*) spname
   integer iph,iseq
   double precision mass
!\end{verbatim}
   integer lokph,loksp
   if(iph.gt.0 .and. iph.le.noofph) then
      lokph=phases(iph)
   else
      gx%bmperr=4050
      goto 1000
   endif
   if(iseq.gt.0 .and. iseq.le.phlista(lokph)%tnooffr) then
      loksp=phlista(lokph)%constitlist(iseq)
      spname=splista(loksp)%symbol
      mass=splista(loksp)%mass
   else
!      write(*,*)'No such constituent'
      gx%bmperr=4096
   endif
!   write(*,*)'3A get_constituent_name: ',iph,iseq,' "',trim(spname),'"'
1000 continue
   return
 end subroutine get_constituent_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_element_data
!\begin{verbatim}
 subroutine get_element_data(iel,elsym,elname,refstat,mass,h298,s298)
! return element data as that is stored as private in GTP
   implicit none
   character elsym*2, elname*(*),refstat*(*)
   double precision mass,h298,s298
   integer iel
!\end{verbatim}
   integer lokel
   if(iel.le.noofel) then
      lokel=elements(iel)
      elsym=ellista(lokel)%symbol
      elname=ellista(lokel)%name
      refstat=ellista(lokel)%ref_state
      mass=ellista(lokel)%mass
      h298=ellista(lokel)%h298_h0
      s298=ellista(lokel)%s298
   else
      gx%bmperr=4042
   endif
 end subroutine get_element_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine new_element_data
!\begin{verbatim}
 subroutine new_element_data(iel,elsym,elname,refstat,mass,h298,s298)
! set new values in an element record, only mass allowed to change ...
   implicit none
   character elsym*2, elname*(*),refstat*(*)
   double precision mass,h298,s298
   integer iel
!\end{verbatim}
   integer lokel
   if(iel.gt.0 .and. iel.le.noofel) then
      lokel=elements(iel)
!      ellista(lokel)%symbol=elsym
!      ellista(lokel)%name)=elname
!      ellista(lokel)%ref_state=refstate
      ellista(lokel)%mass=mass
!      ellista(lokel)%h298_h0=h298
!      ellista(lokel)%s298=s298
   else
      gx%bmperr=4042
   endif
 end subroutine new_element_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_component_name
!\begin{verbatim}
 subroutine get_component_name(icomp,name,ceq)
! return the name of component icomp
   implicit none
   character*(*) name
   integer icomp
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   if(icomp.gt.noofel) then
      gx%bmperr=4052
   else
! strange error buperr set here when plotting q(phase) in step2.OCM ??
      if(buperr.ne.0) then
         write(*,*)'3A buperr set entering get_component_name',buperr
         buperr=0
      endif
      name=splista(ceq%complist(icomp)%splink)%symbol
! no reason buperr should be set here
!      if(buperr.ne.0) then
!         write(*,*)'3A gcn buperr: ',trim(name),buperr
!         gx%bmperr=buperr
!      endif
   endif
1000 continue
   return
 end subroutine get_component_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_species_name
!\begin{verbatim} %-
 subroutine get_species_name(isp,spsym)
! return species name, isp is species number
   implicit none
   character spsym*(*)
   integer isp
!\end{verbatim} %+
   if(isp.le.0 .or. isp.gt.noofsp) then
!      write(*,*)'in get_species_name'
      gx%bmperr=4051; goto 1000
   endif
!   loksp=species(isp)
!   spsym=splista(loksp)%symbol
   spsym=splista(species(isp))%symbol
1000 return
 end subroutine get_species_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_species_location
!\begin{verbatim} %-
 subroutine get_species_location(isp,loksp,spsym)
! return species location and name, isp is species number
   implicit none
   character spsym*(*)
   integer isp,loksp
!\end{verbatim}
   if(isp.le.0 .or. isp.gt.noofsp) then
!      write(*,*)'in get_species_name'
      gx%bmperr=4051; goto 1000
   endif
   loksp=species(isp)
   spsym=splista(loksp)%symbol
!   spsym=splista(species(isp))%symbol
1000 return
 end subroutine get_species_location

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_species_data
!\begin{verbatim}
 subroutine get_species_data(loksp,nspel,ielno,stoi,smass,qsp,nextra,extra)
! return species data, loksp is from a call to find_species_record
! nspel: integer, number of elements in species
! ielno: integer array, element indices
! stoi: double array, stoichiometric factors
! smass: double, mass of species
! qsp: double, charge of the species
! nextra, integer, number of additional values
! extra: double, some additional values like UNIQUAC volume and area
   implicit none
   integer, dimension(*) :: ielno
   double precision, dimension(*) :: stoi,extra
   integer loksp,nspel,nextra
   double precision smass,qsp
!\end{verbatim} %+
   integer jl,iel
   if(loksp.le.0 .or. loksp.gt.noofsp) then
!      write(*,*)'in get_species_data'
      gx%bmperr=4051; goto 1000
   endif
   nspel=splista(loksp)%noofel
   elements: do jl=1,nspel
      iel=splista(loksp)%ellinks(jl)
      ielno(jl)=ellista(iel)%alphaindex
      stoi(jl)=splista(loksp)%stoichiometry(jl)
   enddo elements
   smass=splista(loksp)%mass
   qsp=splista(loksp)%charge
! extraproperties for UNIQUAC model (and maybe others)
   nextra=0
   if(allocated(splista(loksp)%spextra)) then
      nextra=size(splista(loksp)%spextra)
      do jl=1,nextra
         extra(jl)=splista(loksp)%spextra(jl)
      enddo
   endif
1000 return
 end subroutine get_species_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_species_component_data
!\begin{verbatim} %-
 subroutine get_species_component_data(loksp,nspel,compnos,stoi,smass,qsp,ceq)
! return species data, loksp is from a call to find_species_record
! Here we return stoichiometry using components 
! nspel: integer, number of components in species
! compno: integer array, component (species) indices
! stoi: double array, stoichiometric factors
! smass: double, mass of species
! qsp: double, charge of the species
   implicit none
   integer, dimension(*) :: compnos
   double precision, dimension(*) :: stoi(*)
   integer loksp,nspel
   double precision smass,qsp
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer jl,iel,jk,ncomp,locomp,nspx
   integer, allocatable :: components(:)
   double precision, allocatable :: compstoi(:)
! this can be UNIQUAC parameters: area, volume
   double precision qextra(10)
!
! if the components are the elements then use get_species_data
   if(.not.btest(globaldata%status,GSNOTELCOMP)) then
      call get_species_data(loksp,nspel,compnos,stoi,smass,qsp,nspx,qextra)
      goto 1000
!   else
!      write(*,11)globaldata%status,GSNOTELCOMP
!11    format('3A using other components than elements',Z8,i4)
   endif
   allocate(components(noofel))
   allocate(compstoi(noofel))
   components=0
   compstoi=zero
   if(loksp.le.0 .or. loksp.gt.noofsp) then
!      write(*,*)'in get_species_data'
      gx%bmperr=4051; goto 1000
   endif
   nspel=splista(loksp)%noofel
   elements: do jl=1,nspel
! splista(loksp)%ellinks is the location of the element record in ellista
! To find the element index in alphabetical order use the %alphaindex
      iel=ellista(splista(loksp)%ellinks(jl))%alphaindex
! ignore vacancies
      if(iel.le.0) cycle elements
      allcomp: do jk=1,noofel
! this is a loop for all components
! locomp is the species record of the component
         if(abs(ceq%invcompstoi(jk,iel)).gt.1.0D-12) then
! the stoichiometry of this component is nonzero for this element
! add to compstoi(jk)
! convert the element to components using the inverted stoichiometry matrix
! for example elements Ca O Si
! components CaO SiO2 O
! matrix  components/elemenets    Ca   O    Si
!         CaO                     1    1    0
!         SiO2                    0    2    1
!         O                       0    1    0
! inverted matrix                 CaO SiO2  O
!                          Ca     1    0    -1  invmat(1,1) (2,1) (3,1)
!                          O      0    1    0
!                          Si     0    1    -2
! for Ca return 2 components,  1 * CaO -1 * O
! for SiO return 2 components  1*SiO   -1 * O
            compstoi(jk)=compstoi(jk)+&
                 splista(loksp)%stoichiometry(jl)*ceq%invcompstoi(jk,iel)
            qsp=splista(loksp)%charge
         endif
      enddo allcomp
   enddo elements
! return components with nonzero stoichiometry.  
! Note stoichiometry can be negative
! There are always as many components as elements
   smass=zero
   nspel=0
   reduce: do jk=1,noofel
      if(abs(compstoi(jk)).gt.1.0D-12) then
         nspel=nspel+1
         compnos(nspel)=jk
         stoi(nspel)=compstoi(jk)
         smass=smass+stoi(nspel)*ceq%complist(jk)%mass
! maybe save species charge in the component record??
! the lines below needed only if a component is charged !! hopefully never ...
         locomp=ceq%complist(jk)%splink
         qsp=qsp+stoi(nspel)*splista(locomp)%charge
         if(splista(locomp)%charge.ne.zero) then
            write(*,*)'3A charge: ',loksp,qsp,stoi(nspel),splista(locomp)%charge
         endif
      endif
   enddo reduce
1000 return
 end subroutine get_species_component_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_new_stoichiometry
!\begin{verbatim}
 subroutine set_new_stoichiometry(loksp, new_stoi, ispel)
! provided by Clement Instoini
! Change the stoichiometric coefficient of the ispel-th element of loksp-th
! species (the last one when ispel is not given)
! loksp: index of the species (input integer)
! new_stoi: new value of the stoichiometric coefficient (input double precision)
! ispel: index of the element (optional, input integer)
   implicit none
   integer, intent(in):: loksp
   integer, intent(in), optional :: ispel
   double precision, intent(in):: new_stoi
!\end{verbatim}
   character el_name*12,spe_name*24
   integer iel,jl,nspel
   double precision :: old_stoi
   ! number of elements in species
   nspel=splista(loksp)%noofel
   spe_name = trim(splista(loksp)%symbol)
!
   if( .not. present(ispel) ) then
    iel = nspel
    !change the stoichiometric coefficient of the last element
    old_stoi=splista(loksp)%stoichiometry(iel)
    splista(loksp)%stoichiometry(iel)=new_stoi
   else
     iel = ispel
     if (iel.gt.0) then
       ! Change the stoichiometric coefficient of the ispel-th element
       old_stoi=splista(loksp)%stoichiometry(iel)
       splista(loksp)%stoichiometry(iel)=new_stoi
!     else
!       nothing to be done
     end if
   end if
   el_name=ellista(splista(loksp)%ellinks(iel))%name
!   if(ocv()) then
!      write(*,*)"set_new_stoichiometry: (species,element,old_stoi,new_stoi)',&
!           ' = (",spe_name,",",el_name,",",old_stoi,",",new_stoi,")"
!   endif
 end subroutine set_new_stoichiometry

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable double precision function mass_of
!\begin{verbatim}
 double precision function mass_of(component,ceq)
! return mass of component
! smass: double, mass of species
   implicit none
   integer :: component
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   if(component.le.0 .or. component.gt.noofel) then
      write(*,*)'Calling mass_of with illegal component number: ',component
      gx%bmperr=4251; goto 1000
   endif
! return in kg
   mass_of=ceq%complist(component)%mass
1000 return
 end function mass_of

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_phase_name
!\begin{verbatim} %
 subroutine get_phase_name(iph,ics,name)
! Given the phase index and composition set number this subroutine returns
! the name with pre- and suffix for composition sets added and also 
! a \# followed by a digit 1-9 if there are more than one composition sets
   implicit none
   character name*(*)
   integer iph,ics
!\end{verbatim} %+
   character phname*36
   integer lokph,lokcs,kp
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   if(ics.eq.1) then
      name=phlista(lokph)%name
      if(phlista(lokph)%noofcs.ge.2) then
! this was added 2020.04.02 because a call to change_many_phase_status
! using a phase name returned from this routine will suspend all compsets
         kp=len_trim(name)+1
         name(kp:)='#1'
      endif
   else
      kp=len_trim(firsteq%phase_varres(lokcs)%prefix)
      if(kp.gt.0) then
         phname=firsteq%phase_varres(lokcs)%prefix(1:kp)//'_'//&
              phlista(lokph)%name
      else
         phname=phlista(lokph)%name
      endif
      kp=len_trim(firsteq%phase_varres(lokcs)%suffix)
      if(kp.gt.0) then
         phname(len_trim(phname)+1:)='_'//firsteq%phase_varres(lokcs)%suffix
      endif
      phname(len_trim(phname)+1:)='#'//char(ics+ichar('0'))
      name=phname
   endif
1000 continue
   return
 end subroutine get_phase_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_phasetup_name
!\begin{verbatim} %-
 subroutine get_phasetup_name(phtupx,name)
! phasetuple(phtupx)%phase is index to phlista
! the name has pre- and suffix for composition sets added and also 
! a \# followed by a digit 2-9 for composition sets higher than 1.
   implicit none
   character name*(*)
   integer phtupx
!\end{verbatim} %+
   integer phx,phy
!   phx=phlista(phasetuple(phtupx)%phaseix)%alphaindex
   phx=phlista(phasetuple(phtupx)%lokph)%alphaindex
   call get_phase_name(phx,phasetuple(phtupx)%compset,name)
1000 continue
   return
 end subroutine get_phasetup_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_phasetuple_name
!\begin{verbatim} %-
 subroutine get_phasetuple_name(phtuple,name)
! phtuple is a phase tuple
! the name has pre- and suffix for composition sets added and also 
! a \# followed by a digit 2-9 for composition sets higher than 1.
   implicit none
   character name*(*)
   type(gtp_phasetuple) :: phtuple
!\end{verbatim} %+
!   integer phx,phy
   call get_phase_name(phtuple%ixphase,phtuple%compset,name)
1000 continue
   return
 end subroutine get_phasetuple_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!
!-\addtotable subroutine get_phasetup_name_old
!-\begin{verbatim}
! subroutine get_phasetup_name_old(phtuple,name)
! Given the phase tuple this subroutine returns the name with pre- and suffix
! for composition sets added and also a \# followed by a digit 2-9 for
! composition sets higher than 1.
!   implicit none
!   character name*(*)
!   type(gtp_phasetuple) :: phtuple
!-\end{verbatim} %+
!
! PROBABLY REDUNDANT and wrong ...
!
!   call get_phase_name(phtuple%phaseix,phtuple%compset,name)
!   call get_phase_name(phtuple%ixphase,phtuple%compset,name)
!1000 continue
!   return
! end subroutine get_phasetup_name_old
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_phasetup_record
!\begin{verbatim} %-
 subroutine get_phasetup_record(phtx,lokcs,ceq)
! return lokcs when phase tuple known
   implicit none
   integer phtx,lokcs
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   if(phtx.lt.1 .or. phtx.gt.nooftuples) then
!      write(*,*)'Wrong tuple index',phtx
      gx%bmperr=4252; goto 1000
   endif
   write(*,*)'Calling get_phasetup_record is redundant'
   stop
!   lokcs=phlista(phasetuple(phtx)%phaseix)%linktocs(phasetuple(phtx)%compset)
1000 continue
   return
 end subroutine get_phasetup_record

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 !\begin{verbatim}
 subroutine get_sublattice_number(iph,nsl,ceq)
! return the number of sublattices for phase iph
! nsl: integer, number of sublattices
! ceq: pointer, to current gtp_equilibrium_data record
   implicit none
   integer iph,nsl
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph
   nsl = 1
   if(iph.lt.1 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   else
      lokph=phases(iph)
   endif
   nsl=phlista(lokph)%noofsubl
1000 continue
   return

 end subroutine get_sublattice_number

 !/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 !\begin{verbatim}
 subroutine get_sublattice_structure(iph,ics,nsl,nkl,nsites,ceq)
! return the structure of the sublattices for phase iph (ics composition set)
! nsl: integer, number of sublattices
! nkl: integer array, number of constituents in each sublattice
! nsites: double array, number of sites in each sublattice
! ceq: pointer, to current gtp_equilibrium_data record
   implicit none
   integer, intent (in) :: iph,ics,nsl
   integer, dimension(nsl), intent (out) :: nkl, nsites
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer :: i, lokph,lokcs, ncs
!
   if(iph.lt.1 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   else
      lokph=phases(iph)
   endif
   if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
      gx%bmperr=4072; goto 1000
   else
      ncs=max(ics,1)
   endif
! extra check if using saved equilibria which may have less composition sets
   lokcs=phlista(lokph)%linktocs(ncs)
   if(lokcs.le.0) then
      write(*,*)'Index of composition set missing, maybe using a saved equil.'
      gx%bmperr=4072
      goto 1000
   endif
   do i=1,nsl
      nkl(i)=phlista(lokph)%nooffr(i)
      nsites(i)=ceq%phase_varres(lokcs)%sites(i)
   enddo
1000 continue
   return

 end subroutine get_sublattice_structure

 !/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 !\begin{verbatim}
 subroutine get_constituent_data(iph,ics,icons,yarr,charge,csname,ceq)
 ! return the constitution for phase iph (ics composition set)
 ! yarr: double, fraction of constituent
 ! charge: integer, charge of constituent
 ! consname: name of the constituent
 ! ceq: pointer, to current gtp_equilibrium_data record
   implicit none
   integer, intent (in) :: iph,ics,icons
   double precision, intent (inout) :: yarr
   integer, intent (inout) :: charge
   character*(*) , intent (inout) :: csname
   
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer :: i, lokph,lokcs,ncs,loksp

   if(iph.lt.1 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   else
      lokph=phases(iph)
   endif
   if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
      gx%bmperr=4072; goto 1000
   else
      ncs=max(ics,1)
   endif
! extra check if using saved equilibria which may have less composition sets
   lokcs=phlista(lokph)%linktocs(ncs)
   if(lokcs.le.0) then
      write(*,*)'Index of composition set missing, maybe using a saved equil.'
      gx%bmperr=4072
      goto 1000
   endif
   
   yarr=ceq%phase_varres(lokcs)%yfr(icons)
   loksp=phlista(lokph)%constitlist(icons)
   csname=splista(loksp)%symbol
   if(loksp.gt.0) then
      charge = splista(loksp)%charge
   else
      charge=0.D0
   endif
1000 continue
   return
 end subroutine get_constituent_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 !\addtotable subroutine get_phase_data
!\begin{verbatim}
 subroutine get_phase_data(iph,ics,nsl,nkl,knr,yarr,sites,qq,ceq)
! return the structure of phase iph and constituntion of comp.set ics
! nsl: integer, number of sublattices
! nkl: integer array, number of constituents in each sublattice
! knr: integer array, species location (not index) of constituents (all subl)
! yarr: double array, fraction of constituents (in all sublattices)
! sites: double array, number of sites in each sublattice
! qq: double array, (must be dimensioned at least 5) although only 2 used:
! qq(1) is number of real atoms per formula unit for current constitution
! qq(2) is net charge of phase for current constitution
! ceq: pointer, to current gtp_equilibrium_data record
   implicit none
   integer, dimension(*) :: nkl,knr
   double precision, dimension(*) :: yarr,sites,qq
   integer iph,ics,nsl
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,lokcs,kkk,ll,jj,loksp
   double precision vsum,qsum,ql,vl,yz
!
   if(iph.lt.1 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   else
      lokph=phases(iph)
   endif
!   if(gtpdebug.ne.0) write(*,*)'3A get_phase_data 1: ',iph,ics,lokph
   nsl=phlista(lokph)%noofsubl
   if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
      gx%bmperr=4072; goto 1000
   elseif(ics.eq.0) then
      ics=1
   endif
! extra check if using saved equilibria which may have less composition sets
   lokcs=phlista(lokph)%linktocs(ics)
   if(lokcs.le.0) then
      write(*,*)'Index of composition set missing, maybe using a saved equil.'
      gx%bmperr=4072
      goto 1000
   endif
!   if(gtpdebug.ne.0) write(*,*)'3A get_phase_data 10: ',lokcs
!   lokcs=phlista(lokph)%cslink
!   jcs=ics-1
!   do while(jcs.gt.0)
!      lokcs=ceq%phase_varres(lokcs)%next
!      if(lokcs.le.0) then
!         write(*,*)'get_phase_data error 4072'
!         gx%bmperr=4072; goto 1000
!      endif
!      jcs=jcs-1
!   enddo
! >>>>> get_phase_data missing: for ionic liquid sites vary with composition 
   vsum=zero
   qsum=zero
   kkk=0
   if(.not.btest(ceq%phase_varres(lokcs)%status2,CSCONSUS)) then
! CSCONSUS set if a constituent is suspended ... not implemented yet
!      if(gtpdebug.ne.0) write(*,*)'3A get_phase_data 20: ',lokph,nsl
      sublat: do ll=1,nsl
         nkl(ll)=phlista(lokph)%nooffr(ll)
!         if(gtpdebug.ne.0) then
!            write(*,*)'3A get_phase_data 21: ',lokcs,ll,nkl(ll),&
!                 allocated(ceq%phase_varres(lokcs)%sites),&
!                 size(ceq%phase_varres(lokcs)%sites)
!         endif
         if(.not.allocated(ceq%phase_varres(lokcs)%sites)) then
! This can happen for plotting when different dynamic ceq
! have different number of composition sets
            write(*,777)trim(phlista(lokph)%name),ics
777         format('3A site array for phase: ',a,' set ',i2,' not allocated')
            gx%bmperr=4399; goto 1000
         endif
! we get strange error "index 1 or array ceq above bound of 0"
         if(size(ceq%phase_varres(lokcs)%sites).lt.1) then
!            write(*,*)'Strange error when step: ',iph,ics,lokcs,ll
            gx%bmperr=4253; goto 1000
         endif
!         write(*,17)'3 A Strange error: ',iph,ics,lokcs,ll,&
!              size(ceq%phase_varres(lokcs)%sites)
! another strange error "below lower bound of 4 ..."
! I do not now how to check for a lower boundary ... 
17       format(a,10i6)
         sites(ll)=ceq%phase_varres(lokcs)%sites(ll)
!         if(gtpdebug.ne.0) write(*,*)'3A get_phase_data 25: ',sites(ll)
         ql=zero
         vl=zero
!         if(gtpdebug.ne.0) write(*,*)'3A get_phase_data 30: ',&
!              ll,nkl(ll),sites(ll)
         const: do jj=1,nkl(ll)
            kkk=kkk+1
            loksp=phlista(lokph)%constitlist(kkk)
            knr(kkk)=loksp
            yz=ceq%phase_varres(lokcs)%yfr(kkk)
            yarr(kkk)=yz
            if(loksp.gt.0) then
! loksp is -99 for wildcards.  ionic liquid can have that in first sublattice
               ql=ql+yz*splista(loksp)%charge
               if(btest(splista(loksp)%status,SPVA)) then
                  vl=yz
               endif
            endif
         enddo const
         vsum=vsum+sites(ll)*(one-vl)
         qsum=qsum+sites(ll)*ql
      enddo sublat
!      if(gtpdebug.ne.0) write(*,*)'3A get_phase_data 40: ',vsum
      qq(1)=vsum
      qq(2)=qsum
!      write(*,*)'get_phase_data: ',qq(1),qq(2)
   else
! >>>> unfinished handle the case with suspended constituents
!      write(*,*)'get_phase_data with suspended constituents not implemented'
      gx%bmperr=4080; goto 1000
   endif
!
1000 continue
!   if(gtpdebug.ne.0) write(*,*)'3A get_phase_data exit: ',iph,ics
   return
 end subroutine get_phase_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_phase_structure
!\begin{verbatim} %-
 subroutine get_phase_structure(lokph,nsl,nkl)
! return the number of sblattices and constituents in each.
! nsl: integer, number of sublattices
! nkl: integer array, number of constituents in each sublattice
! USED when calculating derivatives of chemical potentials and diffusion coef
   implicit none
   integer, dimension(*) :: nkl
   integer lokph,nsl
!\end{verbatim}
   integer ii
   if(lokph.le.0 .or. lokph.gt.noofph) then
!      write(*,*)'You are way off your head'
      gx%bmperr=4050; goto 1000
   endif
   nsl=phlista(lokph)%noofsubl
   do ii=1,nsl
      nkl(ii)=phlista(lokph)%nooffr(ii)
   enddo
1000 continue
   return
 end subroutine get_phase_structure

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function get_phtuplearray
!\begin{verbatim}
 integer function get_phtuplearray(phcs)
! copies the internal phase tuple array to external software
! function value set to number of tuples
   type(gtp_phasetuple), dimension(*) :: phcs
!\end{verbatim} %+
   integer iz
   do iz=1,nooftuples
! phasetuple(iz)%phase is lokph!!!  .... probably never used ...
      phcs(iz)=phasetuple(iz)
!      phcs(iz)%phase=phasetuple(iz)%phase
!      phcs(iz)%compset=phasetuple(iz)%compset
   enddo
1000 continue
   get_phtuplearray=nooftuples
   return
 end function get_phtuplearray

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     5. Set things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_constitution
!\begin{verbatim}
 subroutine set_constitution(iph,ics,yfra,qq,ceq)
! set the constituent fractions of a phase and composition set and the
! number of real moles and mass per formula unit of phase
! returns number of real atoms in qq(1), charge in qq(2) and mass in qq(3)
! for ionic liquids sets the number of sites in the sublattices
   implicit none
   double precision, dimension(*) :: yfra,qq
   integer iph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokph,lokcs,ll,ml,ic,loksp,jl,locva,zl,zel
   double precision charge,spat,asite,bsite,badd,yz,yva,sumat,asum,bsum,csum
!   double precision charge1,bion1,ionsites(2)
   double precision charge1,bion1,compsum,comp1
! The mass is not calculated correctly in version 2, attempt to fix
   double precision bliq1
! This is needed if we have other components than the elements
   double precision, allocatable :: compam(:),elam(:),iliqcats(:)
!   TYPE(gtp_fraction_set), pointer :: disrec
   logical ionicliq
!   write(*,*)'3A In set_constitution ...',ceq%eqno,iph,ics
   if(iph.le.0 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   endif
   lokph=phases(iph)
   if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
      gx%bmperr=4072; goto 1000
   elseif(ics.eq.0) then
      ics=1
   endif
   lokcs=phlista(lokph)%linktocs(ics)
!   write(*,*)'3A segmentation fault 1',iph,ics,lokcs
   ionicliq=btest(phlista(lokph)%status1,PHIONLIQ)
   if(ionicliq) then
! default values of i2slx
      phlista(lokph)%i2slx(1)=phlista(lokph)%tnooffr+1
      phlista(lokph)%i2slx(2)=phlista(lokph)%tnooffr+1
      yva=zero
      locva=0
   endif
!----
   if(btest(globaldata%status,GSNOTELCOMP)) then
      allocate(elam(noofel))
      allocate(compam(noofel))
      elam=zero
      compam=zero
      if(ionicliq) then
! we must save the amounts on sublattice 1 as we do not know the sites
         allocate(iliqcats(noofel))
         iliqcats=zero
      endif
   endif
!   write(*,*)'3A segmentation fault 10',lokcs
   if(ocv()) write(*,8)'3Ay:',iph,ics,&
        (yfra(ic),ic=1,phlista(lokph)%tnooffr)
8  format(a,2i2,6(1pe11.3))
   nosuscon: if(btest(ceq%phase_varres(lokcs)%status2,CSCONSUS)) then
! >>>> unfinished: handle the case when some constituents are suspended
!      write(*,*)'set_constitution with suspended constituents not implemented'
      write(*,*)'suspended const in: ',lokph,lokcs
      gx%bmperr=4080; goto 1000
   else
! no suspended constituents
! As the application program may have errors first make sure than
! the constituents fractions are correct:
! - no negative fractions
! - sum of fractions in each sublattice unity
!      if(ocv()) write(*,*)'3A 2: ',ionicliq
      ic=0
!      write(*,*)'3A segmentation fault 30',phlista(lokph)%noofsubl
      do ll=1,phlista(lokph)%noofsubl
!         write(*,*)'3A sumy 2: ',ll,ic,phlista(lokph)%noofsubl
         asite=zero
         do ml=1,phlista(lokph)%nooffr(ll)
            yz=yfra(ic+ml)
            if(yz.lt.bmpymin) yz=bmpymin
            ceq%phase_varres(lokcs)%yfr(ic+ml)=yz
            asite=asite+yz
         enddo
! make sure sum of fractions is unity in each sublattice
         do ml=1,phlista(lokph)%nooffr(ll)
            ceq%phase_varres(lokcs)%yfr(ic+ml)=&
                 ceq%phase_varres(lokcs)%yfr(ic+ml)/asite
         enddo
!         write(*,13)'3A y: ',ll,ic,asite,bmpymin,&
!              (ceq%phase_varres(lokcs)%yfr(ic+ml),&
!              ml=1,phlista(lokph)%nooffr(ll))
13       format(a,2i2,2(1pe12.4),1x,4(1pe12.4))
         ic=ic+phlista(lokph)%nooffr(ll)
      enddo
!--------
      ll=1; ml=0; asum=zero; bsum=zero; csum=zero; charge=zero
!      write(*,*)'3A segmentation fault 40'
      if(ionicliq) then
! For ionic liquid we do not know the number of sites
         asite=one
         bion1=zero
      else
         asite=ceq%phase_varres(lokcs)%sites(ll)
      endif
! what is bsite used for???
      bsite=asite; badd=zero
      spat=zero
      allcon: do ic=1,phlista(lokph)%tnooffr
         yz=ceq%phase_varres(lokcs)%yfr(ic)
!         if(ocv()) write(*,*)'3A 3: ',ic,yz
         notva: if(btest(ceq%phase_varres(lokcs)%constat(ic),CONVA)) then
! the constituent is the vacancy
! i2slx(1) should be set to the index of vacancies (if any)
            if(ionicliq) phlista(lokph)%i2slx(1)=ic
            locva=ic
            yva=yz
         else
! sum charge and for constituents with several atoms spat sum number of atoms
            loksp=phlista(lokph)%constitlist(ic)
            charge=charge+bsite*yz*splista(loksp)%charge
! derivates of sites for ionic liquid model
!            if(ocv()) write(*,*)'3A 4: ',loksp,charge
            if(ionicliq) then
               ceq%phase_varres(lokcs)%dpqdy(ic)=abs(splista(loksp)%charge)
!               if(ocv()) write(*,*)'3A dpqdy:     ',&
!                    ic,abs(splista(loksp)%charge)
! i2slx(2) should be set to the index of the first neutral (if any)
               if(splista(loksp)%charge.eq.zero .and.&
                    phlista(lokph)%i2slx(2).gt.ic) &
                    phlista(lokph)%i2slx(2)=ic
            endif
! add the mass of the constituents
            badd=badd+bsite*yz*splista(loksp)%mass
!            write(*,56)'3A badd: ',iph,loksp,splista(loksp)%mass,yz,bsite,badd
56          format(a,2i3,6(1pe12.4))
            sumat=zero
! This is summing atoms per formula unit of the phase
            do jl=1,splista(loksp)%noofel
               sumat=sumat+splista(loksp)%stoichiometry(jl)
            enddo
!--------------------------------------------------------------
            if(btest(globaldata%status,GSNOTELCOMP)) then
! When there are other components than the elements we must sum the number
! of each atom, not just the total. elam was alloctated and zeroed above
               do jl=1,splista(loksp)%noofel
! NOTE that the ellinks specify the location, not alphabetically!!
! we must use %alphaindex to have the alphabetical index of the element ?? YES
                  zel=ellista(splista(loksp)%ellinks(jl))%alphaindex
! FOR IONIC LIQUID MODEL asite is unity and must be updatated below!!
                  elam(zel)=elam(zel)+yz*splista(loksp)%stoichiometry(jl)*asite
!                  write(*,14)'3A elam: ',zel,yz,&
!                     splista(loksp)%stoichiometry(jl),(elam(zl),zl=1,noofel),&
!                       trim(splista(loksp)%symbol)
!14                format(a,i2,5(1pe11.3),2x,a)
               enddo
!               write(*,*)'3A NOTELCOMP: ',compsum,trim(splista(loksp)%symbol)
!               csum=csum+yz*compsum
            endif
            spat=spat+yz*sumat
! check sum number of atoms for ionic liquid
!            if(sumat.gt.1) then
!               write(*,7)'spat: ',lokph,splista(loksp)%noofel,sumat,yz,spat
!7              format(a,2i3,3F10.4)
!            endif
!             write(*,11)loksp,yz,splista(loksp)%mass,badd,bsum
11           format('set_const 3: ',i3,4(1PE15.7))
         endif notva
! ml is constituent number in this sublattice, ic for all sublattices
         ml=ml+1
!         if(ocv()) write(*,*)'3A 5: ',ml
         newsubl: if(ml.ge.phlista(lokph)%nooffr(ll)) then
! next sublattice
            ionliq: if(ionicliq) then
! for ioniq liquids the number of sites is the charge on opposite sublattice
               if(ll.eq.1) then
! Q=\sum_i v_i y_i = charge
!                  write(*,88)'ionliq: ',ll,badd,bion1
88                format(a,i3,6(1pe12.4))
                  ceq%phase_varres(lokcs)%sites(2)=charge
!                  write(*,*)'Ionic 2: ',ceq%phase_varres(lokcs)%sites(2)
!                  bsite=one
                  charge1=charge
                  charge=zero
! same the mass of the constituents on first sublattice
                  bliq1=badd
                  badd=zero
! initiate vacancy and neutral indices beyond last index (already done??)
                  phlista(lokph)%i2slx=phlista(lokph)%tnooffr+1
                  if(btest(globaldata%status,GSNOTELCOMP)) then
                     iliqcats=elam
                     elam=zero
                  endif
               elseif(ll.eq.2) then
! P=\sum_j (-v_j)y_j + Qy_Va. Note charge is total charge and valences 
! on 2nd sublattice is negative
! Now we know number of sites on sublattice 1, update asum and bsum
! Cryptic programming ... sumat is here set to sites on first sublattice
                  sumat=-charge+charge1*yva
                  ceq%phase_varres(lokcs)%sites(1)=sumat
!                  write(*,*)'Ionic 1: ',ceq%phase_varres(lokcs)%sites(1)
                  asum=asum*sumat
                  bsum=bion1*sumat
                  charge=zero
                  if(btest(globaldata%status,GSNOTELCOMP)) then
                     elam=elam+sumat*iliqcats
                  endif
!                  write(*,88)'3A iliq: ',ll,badd,bion1,bsum,sumat,yva
! new way to calculate mass of ionic liquid
                  bsum=sumat*bliq1+ceq%phase_varres(lokcs)%sites(2)*badd
!                  write(*,66)'3A ilmass: ',ll,ceq%phase_varres(lokcs)%sites,&
!                       bliq1,badd,bsum
66                format(a,i3,6(1pe12.4))
                  badd=zero
               else
!                  write(*,*)'Ionic liquid must have two sublattices',ll
                  gx%bmperr=4255; goto 1000
               endif
            endif ionliq
! note: for ionic liquid previous values of asum and bsum are updated 
! when fractions in sublattice 2 have been set
            asum=asum+asite*spat
            bsum=bsum+badd
!            write(*,33)'3A g:',lokcs,ll,asum,asite,spat
33          format(a,2i2,6(1pe12.4))
!            write(*,39)'set_con: ',ll,ml,asum,asite,spat
!39          format(a,2i5,3(1pe12.4))
!            write(*,12)'set_const 12: ',ll,asum,asite,bsum,badd
!12          format(a,i3,4(1pe12.4))
            if(ll.lt.phlista(lokph)%noofsubl) then
               ll=ll+1; ml=0
!               asite=phlista(lokph)%sites(ll); spat=zero
               asite=ceq%phase_varres(lokcs)%sites(ll)
               spat=zero; bion1=badd; badd=zero
! if ionic liquid bsite must be 1.0 when summing second sublattice. Why???
               if(.not.ionicliq) bsite=asite
            endif
         endif newsubl
      enddo allcon
!      write(*,33)'3A h:',lokcs,ll,asum,asite,spat
   endif nosuscon
!   write(*,*)'3A NO segmentation fault 100'
! save charge, number of moles and mass of real atoms per formula unit
!   write(*,33)'3A isum:',lokcs,0,charge,asum,bsum,asite,spat
   ceq%phase_varres(lokcs)%netcharge=charge
   ceq%phase_varres(lokcs)%abnorm(2)=bsum
   if(btest(globaldata%status,GSNOTELCOMP)) then
! Now we can convert the amount of atoms to amount of components
! use ceq%invcompstoi to convert to components
!      write(*,279)'3A elsm: ',iph,asum,(elam(zl),zl=1,noofel)
279   format(a,i3,6(1pe12.4))
      csum=zero
      do zl=1,noofel
         comp1=zero
!         write(*,278)'3A inv: ',(ceq%invcompstoi(zl,zel),zel=1,noofel)
278      format(a,6(1pe12.4))
         do zel=1,noofel
            comp1=comp1+ceq%invcompstoi(zl,zel)*elam(zel)
         enddo
         compam(zl)=comp1
         csum=csum+compam(zl)
      enddo
!      write(*,*)'3A segmentation fault 200'
!      write(*,277)'3A cpam: ',iph,csum,(compam(zl),zl=1,noofel)
277   format(a,i3,6(1pe12.4))
! abnorm(3) is the number of moles of user defined components
!      write(*,299)'3A comp/FU: ',iph,ics,asum,csum
299   format(a,2i3,4(1pe12.4))
      ceq%phase_varres(lokcs)%abnorm(1)=csum
      ceq%phase_varres(lokcs)%abnorm(3)=asum
   else
! if elements are constituents then set abnorm(3)=abnorm(1)
      ceq%phase_varres(lokcs)%abnorm(1)=asum
      ceq%phase_varres(lokcs)%abnorm(3)=asum
   endif
!   write(*,*)'3A sety: ',lokcs,ceq%phase_varres(lokcs)%abnorm(1)
   if(ionicliq .and. locva.gt.0) then
! the ionic liquid vacancy charge is the number of sites on second subl.
      ceq%phase_varres(lokcs)%dpqdy(locva)=ceq%phase_varres(lokcs)%sites(2)
!      if(ocv()) write(*,*)'3A dpqdy(va): ',&
!           locva,ceq%phase_varres(lokcs)%sites(2)
   endif
!   if(ionicliq) then
!      write(*,301)'3A xsc:',lokcs,asum,bsum,ceq%phase_varres(lokcs)%sites,&
!           charge1
!301 format(a,i3,6(1pe12.4))
!      write(*,301)'3A y:  ',ic,ceq%phase_varres(lokcs)%yfr
!   endif
!   write(*,*)'3A NO segmentation fault 300'
   qq(1)=asum
   qq(2)=charge
   qq(3)=bsum
!   write(*,*)'3A segmentation fault 301'
! set disordered fractions if any
   if(btest(phlista(lokph)%status1,phmfs)) then
!now set disordered fractions if any
!      write(*,*)'3A call calc_disfrac for: ',lokph,lokcs
      call calc_disfrac(lokph,lokcs,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,*)'3A segmentation fault 311',lokph,lokcs
   endif
314 format(a,8F8.3)
1000 continue
!   write(*,*)'3A no segmentation fault at exit'
!   if(ionicliq) write(*,*)'3A s_c: ',phlista(lokph)%i2slx
   return
 end subroutine set_constitution

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_reference_state
!\begin{verbatim}
 subroutine set_reference_state(icomp,iph,tpval,ceq)
! set the reference state of a component to be "iph" at tpval
   implicit none
   integer icomp,iph
   double precision, dimension(2) :: tpval
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! NOTE if elements have mixed reference state EQMIXED is set and SER used
! That applies to integral properties like G, S but not MU or AC
   integer nsl,nkl(maxsubl),knr(maxconst),splink,j1,ie,elink
   integer ll,jj,nrel,lokph,noendm,jerr,lokres,ny,endmemx,endmemxy,ics
   double precision sites(maxsubl),qq(5),yarrsave(maxconst),xsum,gmin,gval
   double precision, dimension(:), allocatable :: yarr,xcomp,xmol
   integer, dimension(:), allocatable :: maxjj,jend,jendsave
   double precision tpsave(2),molat,saveg(6)
! iph negative means remove current reference state
   if(iph.lt.0) then
      if(allocated(ceq%complist(icomp)%endmember)) then
! I do not understand the code here any longer but this gave error
! as unallocated when I tried to ser reference state back to SER
         deallocate(ceq%complist(icomp)%endmember)
!      else
!         write(*,4)icomp,ceq%complist(icomp)%phlink
!4        format('3A This component has no previous reference state: ',2i4)
      endif
      ceq%complist(icomp)%phlink=0
      ceq%complist(icomp)%tpref=zero
      ceq%complist(icomp)%refstate='SER (default)'
      goto 1000
   endif
! calculate the composition of the component in mole fractions
   nrel=noel()
   allocate(xcomp(nrel))
   splink=ceq%complist(icomp)%splink
   xcomp=zero
   xsum=zero
   do j1=1,splista(splink)%noofel
      elink=splista(splink)%ellinks(j1)
      ie=ellista(elink)%alphaindex
      xcomp(ie)=splista(splink)%stoichiometry(j1)
      xsum=xsum+xcomp(ie)
   enddo
!   write(*,17)'3A srs x1: ',iph,xsum,(xcomp(ie),ie=1,nrel)
!   do ie=1,splista(splink)%noofel changed 190710/BoS
   do ie=1,nrel
      xcomp(ie)=xcomp(ie)/xsum
   enddo
!   write(*,17)'3A srs x2: ',iph,xsum,(xcomp(ie),ie=1,nrel)
17 format(a,i3,15(f5.2))
! find suitable endmember with correct composition and lowest G
! Note that lowest G is calculated at current T, may be different at another T
! WE CAN HAVE SEVERAL SUBLATTICES ...
   call get_phase_data(iph,1,nsl,nkl,knr,yarrsave,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
   allocate(maxjj(0:nsl))
   allocate(jend(1:nsl))
   allocate(jendsave(1:nsl))
! generate all endmembers, maybe there is a better way ...
! and set unity fraction in yarr and check composition
   ny=0
   maxjj(0)=1
   do ll=1,nsl
      ny=ny+nkl(ll)
      maxjj(ll)=ny
   enddo
   allocate(yarr(ny))
   yarr=zero
   jj=1
   do ll=1,nsl
      yarr(jj)=one
      jend(ll)=jj
      jj=jj+nkl(ll)
   enddo
   allocate(xmol(nrel))
!   lokph=phases(iph)
! we must save the gval for lokres (composition set 1)
   ics=1
   call get_phase_compset(iph,ics,lokph,lokres)
   if(gx%bmperr.ne.0) goto 1000
   gmin=1.0D5
   noendm=0
   tpsave=ceq%tpval
   if(tpval(1).gt.zero) then
! negative tpval means current temperature, else use tpval(1)
      ceq%tpval(1)=tpval(1)
   endif
!   write(*,*)'3A tp: ',tpval(1),ceq%tpval(1)
   ceq%tpval(2)=tpval(2)
   do ie=1,6
      saveg(ie)=ceq%phase_varres(lokres)%gval(ie,1)
   enddo
!   write(*,912)'3G Saved G: ',lokres,ceq%phase_varres(lokres)%gval(1,1),&
!        saveg(1)
!----------------------------------------------
! return here for each endmember
   endmemx=0
200 continue
!   write(*,*)'3G endm: ',(jend(jj),jj=1,nsl)
!   write(*,17)'3G srs y: ',iph,(yarr(jj),jj=1,ny)
   call set_constitution(iph,1,yarr,qq,ceq)
   if(gx%bmperr.ne.0) goto 900
! this subroutine converts site fractions in phase iph, compset 1
! to mole fractions of components (or elements ??? )
   endmemx=endmemx+1
   call calc_phase_mol(iph,xmol,ceq)
   if(gx%bmperr.ne.0) goto 900
!   write(*,202)'3A srs xem: ',iph,endmemx,(xmol(ie),ie=1,nrel)
202 format(a,2i4,15(F5.2))
   do jj=1,nrel
      if(abs(xmol(jj)-xcomp(jj)).gt.1.0D-12) goto 250
   enddo
!--------------------------------------------------
! we have an endmember with the correct composition
   call calcg(iph,1,0,lokres,ceq)
   if(gx%bmperr.ne.0) goto 900
   gval=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
!   write(*,222)'3A srs gval: ',iph,qq(1),gval,gmin,ceq%tpval(1)
222 format(a,i3,F10.3,3(1pe12.4))
   if(gval.lt.gmin) then
! we should check if electrically neutral ??
      noendm=noendm+1
      gmin=gval
      jendsave=jend
      molat=qq(1)
      endmemxy=endmemx
!      write(*,229)'3G min: ',gmin,jendsave
229   format(a,1pe12.4,10i4)
   endif
250 continue
! change constitution .... quit when all endmembers done
   ll=nsl
! should this always be 0?
   maxjj(0)=0
260 continue
! jend is the current endmember
   jj=jend(ll)
   yarr(jj)=zero
   jj=jj+1
   if(jj.gt.maxjj(ll)) then
      jend(ll)=maxjj(ll-1)+1
      yarr(jend(ll))=one
      ll=ll-1
! if ll becomes zero here all endmemebrs have been generated (?)
      if(ll.ge.1) goto 260
   else
      jend(ll)=jj
      yarr(jj)=one
      goto 200
   endif
!----------------------------------------------
   if(noendm.eq.0) then
! if no endmember found this phase cannot be reference phase
!      write(*,*)'This phase cannot be reference state for for this component'
      gx%bmperr=4256; goto 900
   endif
!-----------------------------------------------
! Now we store the reference state and set some bits
! mark that conditions and equilibrium may not be consistent
   ceq%status=ibset(ceq%status,EQINCON)
! endmemx and endmemxy redundant
!   write(*,808)'3G reference state endmember',lokph,endmemxy,jendsave
808 format(a,i3,2x,10i3)
! If all OK then save phase location, endmember array, T and P
   ceq%complist(icomp)%phlink=lokph
   if(.not.allocated(ceq%complist(icomp)%endmember)) then
! if the user changes reference state do not allocate again
!      write(*,*)'3A Allocating endmember for this reference state'
      allocate(ceq%complist(icomp)%endmember(nsl))
   endif
!   write(*,*)'3A refendm: ',icomp,size(ceq%complist),noofel
   ceq%complist(icomp)%endmember=jendsave
!   allocate(ceq%complist(icomp)%endmember(1))
!   ceq%complist(icomp)%endmember=endmemxy
! molat is probably redundant as calcg_endmember returns for one mole component
   ceq%complist(icomp)%molat=molat
! Note tpval(1) can be negative indicating current T
   ceq%complist(icomp)%tpref=tpval
   ceq%complist(icomp)%refstate=phlista(lokph)%name
! NEW 2019.12.02 unless all elements have the same phase and T as reference
! we must set the EQMIXED bit in the CEQ record to enforce use of SER 
! for integral properties like G, H etc.  Element specific MU etc not affected
   allel: do ie=1,noofel
      if(ceq%complist(ie)%refstate.ne.ceq%complist(icomp)%refstate) exit allel
      if(ceq%complist(ie)%tpref(1).ne.ceq%complist(icomp)%tpref(1)) exit allel
      if(ceq%complist(ie)%tpref(2).ne.ceq%complist(icomp)%tpref(2)) exit allel
!      write(*,*)'3A mixed: ',ie,ceq%complist(ie)%tpref,&
!           ceq%complist(ie)%refstate
   enddo allel
! if loop finishes without exit then ie=noofel+1 (Fortran standard)
! and all elements have the same reference state
   if(ie.le.noofel) then
! different phase or T in the elements
!      write(*,*)'3A setting mixed bit'
      ceq%status=ibset(ceq%status,EQMIXED)
   else
! all elements have the same reference phase and T
!      write(*,*)'3A clearing mixed bit'
      ceq%status=ibclr(ceq%status,EQMIXED)
   endif
!-------------------------------------------------------
! restore original constitution of compset 1
!   write(*,*)'3A gval: ',gval
900 continue
   ceq%tpval=tpsave
   jerr=gx%bmperr; gx%bmperr=0
   call set_constitution(iph,1,yarrsave,qq,ceq)
   if(jerr.ne.0) then
      gx%bmperr=jerr
   endif
! restore original values of G and derivatives
   do ie=1,6
      ceq%phase_varres(lokres)%gval(ie,1)=saveg(ie)
   enddo
!   write(*,912)'3G Restored G: ',lokres,ceq%phase_varres(lokres)%gval(1,1),&
!        saveg(1)
912 format(a,i5,6(1pe12.4))
1000 continue
   return
 end subroutine set_reference_state

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine amend_components
!\begin{verbatim}
 subroutine amend_components(line,ceq)
! amend the set of components
   implicit none
   character line*(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer c1,c2,c3,i1,i2,nspel,ierr,lokph,lokcs,nspx
   integer, allocatable :: ielno(:),loksp(:)
   double precision, allocatable :: stoi(:),smass(:),yarr(:)
   double precision qsp,spextra(10),qq(5)
   double precision, allocatable :: matrix(:,:),imat(:,:)
   character name*24
   type(gtp_condition), pointer :: pcond,qcond,last
   type(gtp_equilibrium_data), pointer :: curceq
!
   allocate(loksp(noofel))
   allocate(ielno(noofel))
   allocate(stoi(noofel))
   allocate(smass(noofel))
   allocate(matrix(noofel,noofel))
   matrix=zero
   c2=1
   do c1=1,noel()
      c3=c2+index(line(c2:),' ')
      name=line(c2:c3-1)
!      write(*,*)'3A name: "',trim(name),'"',c3,c1,' "',trim(line(c3:)),'"'
      c2=c3
      call find_species_record_exact(name,loksp(c1))
      if(gx%bmperr.ne.0) goto 1000
      call get_species_data(loksp(c1),nspel,ielno,stoi,&
           smass(c1),qsp,nspx,spextra)
      if(qsp.gt.zero) then
         write(*,*)'Charged species must not be components'
         gx%bmperr=4399; goto 1000
      endif
      do i1=1,nspel
         matrix(ielno(i1),c1)=stoi(i1)
      enddo
!      do i1=1,nspel
!         matrix(c1,ielno(i1))=stoi(i1)
!      enddo
   enddo
!   do c1=1,noofel
!      write(*,70)'3A mat: ',c1,(matrix(c2,c1),c2=1,noofel)
!   enddo
70 format(a,i1,6(1pe12.4))
! check that the matrix has an inverse
   allocate(imat(noofel,noofel))
! removed second index as not used!
!   call mdinvold(noofel,noofel+1,matrix,imat,noofel,ierr)
   call mdinvold(noofel,matrix,imat,noofel,ierr)
   if(ierr.eq.0) then
!      write(*,*)'Error inverting component matrix, dependent components'
      gx%bmperr=4362; goto 1000
   endif
!   do c1=1,noofel
!      write(*,70)'3A imt: ',c1,(imat(c2,c1),c2=1,noofel)
!   enddo
!   gx%bmperr=4399
!   write(*,*)'3A *** All seems OK so far ... but only testing yet'
!   goto 1000
!----------------------------------------------------------
! We have a new set of components!!
! At present (and maybe forever) use the same components in all equilibria ...
   do c1=1,noofel
      do c2=1,noofel
         ceq%compstoi(c2,c1)=matrix(c2,c1)
         ceq%invcompstoi(c2,c1)=imat(c2,c1)
! set bit GSNOTELCOMP if there are non-zero off-diagonal terms in invcompstoi 
         if(c1.ne.c2 .and. imat(c2,c1).ne.zero) then
            globaldata%status=ibset(globaldata%status,GSNOTELCOMP)
         endif
      enddo
!   enddo
! enter the components, no alphabetical order ... ??
!   do c1=1,noofel
      ceq%complist(c1)%splink=loksp(c1)
      ceq%complist(c1)%phlink=0
      ceq%complist(c1)%tpref(1)=2.9815D2
      ceq%complist(c1)%tpref(2)=1.0D5
      ceq%complist(c1)%mass=smass(c1)
   enddo
! delete all conditions and experiments in all equilibria
! the argument 0 means only conditions and experiments deleted, 
! not the ceq itself
!   write(*,*)'3A deleting all conditions in all equilibria',eqfree-1
   do c1=1,eqfree-1
      curceq=>eqlista(c1)
      call delete_all_conditions(0,curceq)
! delete if there are some extra things
      if(allocated(curceq%eqextra)) deallocate(curceq%eqextra)
   enddo
! we must go through all (stoichiometric?) phases and set a new
! value for abnorm(1) and (3)
!   write(*,*)'3A update asum and csum for all phases'
   do c1=1,noofph
! the value stored in phases(i) is the location of phase record!!
      lokph=phases(c1)
      do c2=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(c2)
         c3=size(ceq%phase_varres(lokcs)%yfr)
         if(.not.allocated(yarr)) then
            allocate(yarr(c3))
         endif
         yarr=ceq%phase_varres(lokcs)%yfr
! this will update abnorm(1) and (3) for THIS equilibrium ... loop for all??
         call set_constitution(c1,c2,yarr,qq,ceq)
      enddo
      deallocate(yarr)
   enddo
1000 continue
! deallocate temporary things (maybe default?)
   deallocate(loksp)
   deallocate(ielno)
   deallocate(stoi)
   deallocate(smass)
   deallocate(matrix)
   return
 end subroutine amend_components

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function gettupix
!\begin{verbatim}
 integer function gettupix(iph,ics)
! convert phase and compset index to tuple index
   implicit none
   integer iph,ics
!\end{verbatim}
   integer ii,tupix
   ii=ics
   tupix=iph
   loop: do while(ii.gt.1)
      tupix=phasetuple(tupix)%nextcs
      if(tupix.le.0) then
         gx%bmperr=4072; exit loop
      endif
      ii=ii-1
   enddo loop
   return
 end function gettupix

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine suspend_somephases
!\begin{verbatim}
  subroutine suspend_somephases(mode,invph,dim1,dim2,ceq)
! This was added to handle calculating restricted equilibria during mapping
! to suspend (mode=1) or restore (mode=0) phases not involved
! in an invariant equilibrium.
! invph is array with phases that are involved, it has dimension (dim1,*)
! the current status is saved and restored 
    implicit none
    integer mode,dim1,dim2,invph(dim1,*)
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer, save, allocatable, dimension(:) :: phtupixstatus
    integer, save :: ntup
    integer ii,jj,kk,lokcs,lokph
    character phname*24
    ii=nooftup()
    kk=0
    if(mode.eq.1) then
! after saving current status suspend all phases not included in invph
!       write(*,*)'3A suspending some phases',ii
       ntup=ii
       if(allocated(phtupixstatus)) then
          write(*,*)'3A calls to suspend_somephases cannot be nested'
          gx%bmperr=4399; goto 1000
       else
          allocate(phtupixstatus(ntup))
       endif
       loop1: do ii=1,ntup
          lokcs=phasetuple(ii)%lokvares
          phtupixstatus(ii)=ceq%phase_varres(lokcs)%phstate
          do jj=1,dim2
!             write(*,*)'3A suspend? ',jj,lokcs,&
!                  phlista(invph(1,jj))%linktocs(invph(2,jj)),phtupixstatus(ii)
! invph(1,jj) is index in phases (phase and alphabetcal order)
! lokph is the order the phase were entered into phlista (arbitrary)
             lokph=phases(invph(1,jj))
             if(lokcs.eq.phlista(lokph)%linktocs(invph(2,jj))) then
!                write(*,'(a,6i5)')'3A not suspending',jj,invph(1,jj),&
!                     invph(2,jj),phlista(lokph)%linktocs(invph(2,jj))
                cycle loop1
             endif
          enddo
! this phase should be suspended
          kk=kk+1
          ceq%phase_varres(lokcs)%phstate=PHSUS
       enddo loop1
!       write(*,'(a,i3,a,i3)')'3A suspededed ',kk,' phases out of ',ntup
    elseif(mode.eq.0) then
! restore status of all phases except those in invph
!       write(*,*)'3A restoring some phases',ii
       if(ii.ne.ntup) then
          write(*,*)'3A number of phases and compsets changed!',ntup,ii
          stop
       endif
       do ii=1,ntup
          ceq%phase_varres(phasetuple(ii)%lokvares)%phstate=phtupixstatus(ii)
       enddo
!       write(*,'(a,i3,a)')'3A restored phase status for ',ntup,' phases'
       deallocate(phtupixstatus)
    else
       write(*,*)'3A mode must be 0 or 1'
       gx%bmperr=4399
    endif
1000 continue
    return
  end subroutine suspend_somephases

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine delete_unstable_compsets
!\begin{verbatim}
  subroutine delete_unstable_compsets(lokph,ceq)
! This was added to explictly delete unstable composition sets with AUTO set
! Compsets will be shifted down if a stable compset is after an unstable
! See subroutine TOTO_AFTER in gtp3Y.F90
!
    implicit none
    integer lokph
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer ii,iph,lokcs
    write(*,*)'3A delete unstable compsets for phase: ',&
         trim(phlista(lokph)%name),phlista(lokph)%noofcs
! the first composition sets cannot be deleted even if unstable
    do ii=phlista(lokph)%noofcs,2,-1
       lokcs=phlista(lokph)%linktocs(ii)
       write(*,100)ii,btest(ceq%phase_varres(lokcs)%status2,CSAUTO),&
            btest(ceq%phase_varres(lokcs)%status2,CSTEMPAR)
100    format('3A compset: ',i2,' bits: ',2l2)
    enddo
!    call remove_composition_set(iph,.FALSE.)
    write(*,*)'Not implemented yet'
1000 continue
    return
  end subroutine delete_unstable_compsets

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

