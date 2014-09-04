!
! included in pmod25.F90
!
!****************************************************
! general subroutines for creating and handling elements, species, phases etc
! accessable externally

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     1. Initialization
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine init_gtp(intvar,dblvar)
! initiate the data structure
! create element and species record for electrons and vacancies
! the allocation of many arrays should be provided calling this routne
! these will eventually be used for allocations and defaults
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
   allocate(phasetuple(0:2*maxph))
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
   ellista(0)%ref_state='Vaccum'
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
      eqlista(jl)%next=jl+1
   enddo
   eqlista(maxeq)%next=-1
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
   firsteq%phase_varres(2*maxph)%nextfree=-1
! csfree is not declared ... how can that be?? where is it declared ??
   csfree=1
! convergence criteria for constituent fractions, 1e-6 works most often
! But one should take care to equilibrate fractions smaller than xconv!!!
   firsteq%xconv=1.0D-6
   firsteq%maxiter=500
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
!============================================================
! VERY IMPORTANT: The properties defined below must not be equal to state
! variables, if so they cannot be listed and other errors may occur
!
! ANY CHANGES HERE MUST BE MADE ALSO IN SUBROUTINE state_variable_val, pmod25c
!
!============================================================
! Mixed Curie/Neel Temperature, set bits that TC and BM cannot depend on T
   npid=npid+1
   propid(npid)%symbol='TC '
   propid(npid)%note='Mix Curie/Neel T'
   propid(npid)%status=0
! TC cannot depend on T but on P
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Average Bohr magneton number
   npid=npid+1
   propid(npid)%symbol='BMAG '
   propid(npid)%note='Aver Bohr magn no'
   propid(npid)%status=0
! BM cannot depend on either T or P ??
   propid(npid)%status=ibset(propid(npid)%status,IDNOTP)
!.......................................
! Specific Curie temperature
   npid=npid+1
   propid(npid)%symbol='CTA '
   propid(npid)%note='Curie temperature'
   propid(npid)%status=0
! CTA cannot depend on either T or P ??
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Specific Neel temperature
   npid=npid+1
   propid(npid)%symbol='NTA '
   propid(npid)%note='Neel temperature'
   propid(npid)%status=0
! NTA cannot depend on T but on P
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Individual Bohr magneton number
   npid=npid+1
   propid(npid)%symbol='IBM '
   propid(npid)%note='Ind. Bohr magn no'
   propid(npid)%status=0
! IBM cannot depend on either T or P and it is individual
   propid(npid)%status=ibset(propid(npid)%status,IDCONSUFFIX)
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Debye or Einstein temperature
   npid=npid+1
   propid(npid)%symbol='THETA '
   propid(npid)%note='Debye or Einst temp'
   propid(npid)%status=0
! THETA cannot depend on T but on P
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! logarithm of individual mobility
   npid=npid+1
   propid(npid)%symbol='MQ '
   propid(npid)%note='LN mob. of const.'
   propid(npid)%status=0
! MQ is specific för a constituent
   propid(npid)%status=ibset(propid(npid)%status,IDCONSUFFIX)
!.......................................
! Electrical resistivity 
   npid=npid+1
   propid(npid)%symbol='RHO '
   propid(npid)%note='Elect resistivity'
   propid(npid)%status=0
!.......................................
! Magnetic suseptibility 
   npid=npid+1
   propid(npid)%symbol='MAGS '
   propid(npid)%note='Magn suseptibility'
   propid(npid)%status=0
!.......................................
! Glas trasition temperature 
   npid=npid+1
   propid(npid)%symbol='GTT '
   propid(npid)%note='Glas trans temperature'
   propid(npid)%status=0
! Cannot depend on temperature 
   propid(npid)%status=ibset(propid(npid)%status,IDONLYP)
!.......................................
! Viscosity 
   npid=npid+1
   propid(npid)%symbol='VISCA '
   propid(npid)%note='Viscosity'
   propid(npid)%status=0
!.......................................
! Lattice parameter in direction X
   npid=npid+1
   propid(npid)%symbol='LPX '
   propid(npid)%note='Lat par X axis'
   propid(npid)%status=0
! lattice parameters may depend on T and P
!.......................................
! Lattice parameter in direction Y
   npid=npid+1
   propid(npid)%symbol='LPY '
   propid(npid)%note='Lat par Y axis'
   propid(npid)%status=0
! lattice parameters may depend on T and P
!.......................................
! Lattice parameter in direction Z
   npid=npid+1
   propid(npid)%symbol='LPZ '
   propid(npid)%note='Lat par Z axis'
   propid(npid)%status=0
! lattice parameters may depend on T and P
!.......................................
! This is an angle for non-cubic lattices
   npid=npid+1
   propid(npid)%symbol='LPTH '
   propid(npid)%note='Lat angle TH'
   propid(npid)%status=0
! Angle may depend on T and P
!.......................................
! This is an elastic "constant"
   npid=npid+1
   propid(npid)%symbol='EC11 '
   propid(npid)%note='Elast const C11'
   propid(npid)%status=0
! The elastic constant may depend on T and P
!.......................................
! This is another elastic "constant"
   npid=npid+1
   propid(npid)%symbol='EC12 '
   propid(npid)%note='Elast const C12'
   propid(npid)%status=0
! The elastic constant may depend on T and P
!.......................................
! This is yet another elastic "constant"
   npid=npid+1
   if(npid.gt.maxprop) then
      write(*,*)'Too many parameter identifiers, increase maxprop'
      gx%bmperr=7777; goto 1000
   endif
   propid(npid)%symbol='EC44 '
   propid(npid)%note='Elast const C44'
   propid(npid)%status=0
! The elastic constant may depend on T and P
!.......................................
! IMPORTRANT: When adding more parameter identifiers one should never
! use a name ending in D as that would be taken as a "disordered"
! The number of defined properties, should be less than maxprop
! IMPORTANT: In the addition records one must use the parameter identifier
! to extract the calculated composition dependent values
   ndefprop=npid
!-------------------------------------------------
   highcs=0
! globaldata record; set gas constant mm
   globaldata%status=0
! set advanced user, no data, no phase, no equilibrium calculated
   globaldata%status=ibset(globaldata%status,GSADV)
   globaldata%status=ibset(globaldata%status,GSNODATA)
   globaldata%status=ibset(globaldata%status,GSNOPHASE)
   firsteq%status=ibset(firsteq%status,EQNOEQCAL)
! set gas constant and some default values
   globaldata%name='current'
   globaldata%rgas=8.31451D0
! more recent value
!   globaldata%rgas=8.3144621D0
! old value of gas constant
   globaldata%rgasuser=8.31451D0
   globaldata%pnorm=one
!   write(*,*)'init_gtp: enter R and RTLNP'
! enter R as TP function
   tpname='R'
!   write(tpfun,777)' 10 8.31451; 20000 N '
!777 format(a)
!   call enter_tpfun(tpname,tpfun,lrot,.FALSE.)
   call enter_tpconstant(tpname,globaldata%rgas)
   if(gx%bmperr.ne.0) goto 1000
   tpname='RTLNP'
   tpfun=' 10 R*T*LN(1.0D-5*P); 20000 N '
   call enter_tpfun(tpname,tpfun,lrot,.FALSE.)
   if(gx%bmperr.ne.0) goto 1000
! default minimum fraction
   bmpymin=ymind
! putfun error code .... should use buperr at least
   pfnerr=0
!------------------------------------
! allocate array for state variable function
!   write(*,*)'init_gtp: allocate array for state variable functions'
   allocate(svflista(maxsvfun))
! number of state variable function
   nsvfun=0
! zero the array with equilibrium index for functions, not used aywhere??
!   pflocal=0
! enter some useful state variable function
   tpfun=' R=8.31451;'
   ip=1
!   write(*,*)'init_gtp: entering function R'
   call enter_svfun(tpfun,ip,firsteq)
!   if(gx%bmperr.ne.0) then
!      write(*,*)'Error entering R',gx%bmperr
!      goto 1000
!   endif
!   write(*,*)'Entered symbol R'
   tpfun=' RT=R*T;'
   ip=1
!   write(*,*)'init_gtp: entering function RT'
   call enter_svfun(tpfun,ip,firsteq)
!   if(gx%bmperr.ne.0) then
!      write(*,*)'Error entering symbol RT'
!      goto 1000
!   endif
!   write(*,*)'Entered symbol RT'
   tpfun=' T_C=T-273.15;'
   ip=1
!   write(*,*)'init_gtp: entering function T_C'
   call enter_svfun(tpfun,ip,firsteq)
!   if(gx%bmperr.ne.0) then
!      write(*,*)'Error entering symbol T_C'
!      goto 1000
!   endif
! finished initiating
1000 continue
!   write(*,*)'exit from init_gtp'
   return
 END subroutine init_gtp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     2. Number of things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 integer function noel()
! number of elements because noofel is private
! should take care if elements are suspended
!\end{verbatim} %+
   noel=noofel
 end function noel

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 integer function nosp()
! number of species because noofsp is private
! should take care if species are suspended
!\end{verbatim} %+
   nosp=noofsp
 end function nosp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 integer function noph()
! number of phases because noofph is private
! should take care if phases are hidden
!\end{verbatim} %+
   noph=noofph
 end function noph

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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

!\begin{verbatim} %-
 integer function noconst(iph,ics,ceq)
! number of constituents for iph (include single constituents on a sublattice)
! It tests if a constituent is suspended which can be different in each ics.
   implicit none
   integer iph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
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
! some constituents suspended
         noc=phlista(lokph)%tnooffr
         do jl=1,phlista(lokph)%tnooffr
            if(btest(ceq%phase_varres(lokcs)%constat(jl),CONSUS)) then
               noc=noc-1
            endif
         enddo
         noconst=noc
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

!\begin{verbatim}
 integer function nosvf()
! number of state variable functions
!\end{verbatim}
   implicit none
   nosvf=nsvfun
   return
 end function nosvf

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

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
!>     3. Find things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine find_element_by_name(name,iel)
! find an element index by its name
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

!\begin{verbatim} %-
 subroutine find_species_by_name(name,isp)
! locates a species index from its name
   implicit none
   character name*(*)
   integer isp
!\end{verbatim} %+
   character symbol*24
   integer loksp
   symbol=name
   call capson(symbol)
   do loksp=1,noofsp
!     write(*,*)'find species 2: ',loksp,splista(loksp)%symbol
      if(compare_abbrev(symbol,splista(loksp)%symbol)) then
         isp=splista(loksp)%alphaindex
         goto 1000
      endif
   enddo
!    write(*,*)'in find_species_by_name'
   gx%bmperr=4051
   loksp=0
1000 continue
   return
 end subroutine find_species_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine find_species_record(name,loksp)
! locates a species record allowing abbreviations
   implicit none
   character name*(*)
   integer loksp
!\end{verbatim} %+
   character symbol*24
   symbol=name
   call capson(symbol)
   do loksp=1,noofsp
!       write(*,17)'find species 17: ',loksp,splista(loksp)%symbol,name
!17     format(a,i3,' "',a,'" "',a,'"')
       if(compare_abbrev(symbol,splista(loksp)%symbol)) goto 1000
   enddo
!   write(*,*)'Error in find_species_record "',name,'"'
   gx%bmperr=4051
   loksp=0
1000 continue
   return
 end subroutine find_species_record

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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

!\begin{verbatim} %-
 subroutine find_species_record_exact(name,loksp)
! locates a species record, exact match needed
! for parameters, V must not be accepted as abbreviation of VA or C for CR
   implicit none
   integer loksp
   character name*(*)
!\end{verbatim} %+
   character symbol*24
   symbol=name
   call capson(symbol)
   do loksp=1,noofsp
!       write(*,17)'find species 17: ',loksp,splista(loksp)%symbol,name
!17     format(a,i3,' "',a,'" "',a,'"')
      if(symbol.eq.splista(loksp)%symbol) goto 1000
   enddo
!    write(*,*)'in find_species_record'
   gx%bmperr=4051
   loksp=0
1000 continue
   return
 end subroutine find_species_record_exact

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine find_phase_by_name(name,iph,ics)
! finds a phase with name "name", returns address of phase, first fit accepted
! handles composition sets either with prefix/suffix or #digit
! no pre/suffix nor # gives first composition set
   implicit none
   character name*(*)
   integer iph,ics
!\end{verbatim} %+
   character name1*36,csname*36,name2*24
   TYPE(gtp_phase_varres), pointer :: csrec
   integer kp,kcs,lokph,jcs,lokcs
! convert to upper case locally
   name1=name
   call capson(name1)
!   write(*,*)'find phase: ',name1
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
!17  format('find_phase 3: ',2i4)
   loop1: do lokph=1,noofph
      name2=phlista(lokph)%name
!      write(*,*)'find_phase 2: ',name1,name2
      if(compare_abbrev(name1,name2)) then
         if(ics.le.phlista(lokph)%noofcs) then
            goto 300
         else
!            write(*,18)ics,phlista(lokph)%noofcs
!18  format('find_phase 4: ',2i4)
            gx%bmperr=4072; goto 1000
         endif
      endif
! if there are composition sets check name including prefix/suffix
!      write(*,*)'find_phase 5: ',lokph,phlista(lokph)%noofcs
      csno: do jcs=2,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(jcs)
!         write(*,*)'25A: find phase: ',jcs,lokcs,phlista(lokph)%noofcs
         csrec=>firsteq%phase_varres(lokcs)
         kp=len_trim(csrec%prefix)
         if(kp.gt.0) then
            csname=csrec%prefix(1:kp)//'_'//name2
         else
            csname=name2
         endif
!         write(*,*)'find phase: ',kp
         kp=len_trim(csrec%suffix)
         if(kp.gt.0) csname=csname(1:len_trim(csname))//'_'//&
              csrec%suffix(1:kp)
!         write(*,244)ics,kcs,jcs,kp,name1(1:len_trim(name1)),&
!              csname(1:len_trim(csname))
244      format('25A: find_phase: ',4i3,'<',a,'>=?=<',a,'>')
         if(compare_abbrev(name1,csname)) then
! if user has provided both #<digit> and pre/suffix these must be consistent
            if(kcs.gt.0 .and. kcs.ne.jcs) then
! automatically created composition sets all have the suffix _AUTO but
! can have several numbers
!               write(*,*)'25A: mix? ',jcs,phlista(lokph)%noofcs
               if(jcs.eq.phlista(lokph)%noofcs) goto 1100
               cycle csno
            endif
            ics=jcs
            goto 300
         endif
      enddo csno
   enddo loop1
! no phase found
   gx%bmperr=4050
   goto 1000
300 continue
   iph=phlista(lokph)%alphaindex
   gx%bmperr=0
1000 continue
   return
1100 continue
   gx%bmperr=4073
   goto 1000
 END subroutine find_phase_by_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
! possible phase, if iphfound>0 exact match is required
            if(iphfound.ne.0) then
               if(name1.eq.name2) then
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
   gx%bmperr=4073
   goto 1000
 END subroutine find_phase_by_name_exact

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
   character spname1*24
   integer lokph,kp,ll,kk,loksp,ls
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
   lloop: do ll=1,phlista(lokph)%noofsubl
      sploop: do kk=1,phlista(lokph)%nooffr(ll)
         icon=icon+1
         if(ls.eq.0 .or. ls.eq.ll) then
            loksp=phlista(lokph)%constitlist(icon)
! constituent icon is the requested one
!             write(*,55)ll,kk,icon,spname1(1:3),splista(loksp)%symbol(1:3)
!55           format('find_const 7: ',3i3,1x,a,2x,a)
            if(compare_abbrev(spname1,splista(loksp)%symbol)) then
               mass=splista(loksp)%mass
               goto 1000
            endif
         endif
      enddo sploop
   enddo lloop
   gx%bmperr=4096
1000 continue
   return
 end subroutine find_constituent

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
!   write(*,*)'25A equil: ',name2(1:20),ieq
   if(compare_abbrev(name2,'PREVIOUS ')) then
      jeq=max(1,ieq-1); goto 200
   elseif(compare_abbrev(name2,'FIRST ')) then
      jeq=1; goto 200
   endif
   jeq=0
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
      write(*,*)'No such constituent'
      gx%bmperr=7777
   endif
1000 continue
   return
 end subroutine get_constituent_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine get_phase_constituent_name(iph,icon,name)
! return the name of constituent icon of phase iph
! redundant?
   implicit none
   character*(*) name
   integer iph,icon
!\end{verbatim}
   integer lokph
   if(iph.le.0 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   endif
   lokph=phases(iph)
   if(icon.le.0 .or. icon.gt.phlista(lokph)%tnooffr) then
      gx%bmperr=4096; goto 1000
   endif
!   loksp=phlista(lokph)%constitlist(icon)
!   name=splista(loksp)%symbol
   name=splista(phlista(lokph)%constitlist(icon))%symbol
1000 continue
   return
 end subroutine get_phase_constituent_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
      name=splista(ceq%complist(icomp)%splink)%symbol
      if(buperr.ne.0) then
         write(*,*)'gcn buperr: ',buperr
         gx%bmperr=buperr
      endif
   endif
1000 continue
   return
 end subroutine get_component_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine get_species_name(isp,spsym)
! return species name, isp is species number
   implicit none
   character spsym*(*)
   integer isp
!\end{verbatim}
   if(isp.le.0 .or. isp.gt.noofsp) then
!       write(*,*)'in get_species_name'
      gx%bmperr=4051; goto 1000
   endif
!   loksp=species(isp)
!   spsym=splista(loksp)%symbol
   spsym=splista(species(isp))%symbol
1000 return
 end subroutine get_species_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine get_species_data(loksp,nspel,ielno,stoi,smass,qsp)
! return species data, loksp is from a call to find_species_record
! nspel: integer, number of elements in species
! ielno: integer array, element indices
! stoi: double array, stocichiometric factors
! smass: double, mass of species
! qsp: double, charge of the species
   implicit none
   integer, dimension(*) :: ielno
   double precision, dimension(*) :: stoi(*)
   integer loksp,nspel
   double precision smass,qsp
!\end{verbatim}
   integer jl,iel
   if(loksp.le.0 .or. loksp.gt.noofsp) then
!       write(*,*)'in get_species_data'
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
1000 return
 end subroutine get_species_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
      gx%bmperr=7777; goto 1000
   endif
! return in kg
   mass_of=ceq%complist(component)%mass
1000 return
 end function mass_of

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine phase_name(phtuple,name)
! Given the phase tuple this subroutine returns the name with pre- and suffix
! for composition sets added and also a \# followed by a digit 2-9 for
! composition sets higher than 1.
   implicit none
   character name*(*)
   type(gtp_phasetuple) :: phtuple
!\end{verbatim} %+
   call get_phase_name(phtuple%phase,phtuple%compset,name)
1000 continue
   return
 end subroutine phase_name

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine get_phase_name(iph,ics,name)
! Given the phase index and composition set number this subroutine returns
! the name with pre- and suffix for composition sets added and also 
! a \# followed by a digit 2-9 for composition sets higher than 1.
   implicit none
   character name*(*)
   integer iph,ics
!\end{verbatim}
   character phname*36
   integer lokph,lokcs,kp
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   if(ics.eq.1) then
      name=phlista(lokph)%name
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
!\end{verbatim}
   integer lokph,lokcs,kkk,ll,jj,loksp
   double precision vsum,qsum,ql,vl,yz
!
   if(iph.lt.1 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   else
      lokph=phases(iph)
   endif
   nsl=phlista(lokph)%noofsubl
   if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
      gx%bmperr=4072; goto 1000
   elseif(ics.eq.0) then
      ics=1
   endif
   lokcs=phlista(lokph)%linktocs(ics)
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
      sublat: do ll=1,nsl
         nkl(ll)=phlista(lokph)%nooffr(ll)
! For ionic liquid one must use this but at present values are not set
         sites(ll)=ceq%phase_varres(lokcs)%sites(ll)
!         sites(ll)=phlista(lokph)%sites(ll)
!         write(*,*)'get_phase_data 7:',lokcs,ll,sites(ll)
         ql=zero
         vl=zero
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
   return
 end subroutine get_phase_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine extract_stvr_of_condition(pcond,nterm,coeffs,statevar)
! finds a condition record with the given state variable record
! returns it as a state variable record !!!
! nterm: integer, number of terms in the condition expression
! pcond: pointer, to a gtp_condition record
   implicit none
   TYPE(gtp_condition), pointer :: pcond
! ONE CANNOT HAVE ARRAYS OF POINTERS!!! STUPID
!   TYPE(gtp_state_variable), dimension(*), pointer :: statevar
   TYPE(gtp_state_variable), dimension(*) :: statevar
   integer nterm
   double precision coeffs(*)
!\end{verbatim}
   TYPE(gtp_condition), pointer :: last,current,first
   integer, dimension(4) :: indx
   integer ncc,nac,j1,j2,istv,iref,iunit
!
   write(*,*)'not implemented!!'
   gx%bmperr=7777; goto 1000
!--------------------------------------------------------
   if(.not.associated(pcond)) goto 900
   first=>pcond%next
   current=>first
!   write(*,*)'get_condition start: ',current%statev,current%active
   ncc=1
   nac=0
!   write(*,98)'new:',0,nterm,istv,(indices(i,1),i=1,4),iref,iunit
98 format(a,2x,i2,5x,2i4,5x,4i4,5x,2i3)
100 continue
!   write(*,98)'old:' ,current%nid,current%noofterms,current%statev,&
!        (current%indices(i,1),i=1,4),current%iref,current%iunit
   if(nterm.eq.0) then
!      write(*,*)'get_condition: ',istv,ncc,nac
      if(current%active.eq.0) then
! this call just looks for active condition istv
         nac=nac+1
! why should fix phase conditions have istv=nac?? Check!!
         if(nac.eq.istv) then
! a condition specified like this must not be a phase status change
            if(current%statev.lt.0) then
            write(kou,*)'You must use "set phase status" to change fix status'
            else
               goto 150
            endif
         endif
      endif
      goto 200
   endif
   if(current%noofterms.ne.nterm .or. current%statev.ne.istv .or. &
        current%iref.ne.iref .or. current%iunit.ne.iunit) goto 200
   do j1=1,nterm
!      do j2=1,4
!         if(current%indices(j2,j1).ne.indices(j2,j1)) goto 200
!      enddo
   enddo
150 continue
! found condition
   pcond=>current
!   write(*,*)'Found condition: ',pcond%nid,ncc
   goto 1000
200 continue
   current=>current%next
   ncc=ncc+1
   if(.not. associated(current,first)) goto 100
900 continue
! no such condition
   gx%bmperr=4131; goto 1000
1000 continue
   return
 end subroutine extract_stvr_of_condition

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     5. Set things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!
!-\begin{verbatim}
! subroutine set_mean_constitution(iph,ics,ceq)
! sets a start constitution 1/ns in all sublattices where ns is the number
! of constituents in the sublattice.
!   implicit none
!   integer iph,ics
!   TYPE(gtp_equilibrium_data), pointer :: ceq
!-\end{verbatim}
!   integer, dimension(maxsubl) :: knl
!   double precision, dimension(maxsubl) :: sites
!   integer, dimension(maxconst) :: knr
!   double precision, dimension(maxconst) :: yarr
!   double precision, dimension(5) :: qq
!   double precision df
!   integer nsl,kkk,ll,jl
!   call get_phase_data(iph,ics,nsl,knl,knr,yarr,sites,qq,ceq)
!   if(gx%bmperr.ne.0) goto 1000
!   kkk=0
!   do ll=1,nsl
!      if(knl(ll).gt.1) then
!         df=one/dble(knl(ll))
!         do jl=1,knl(ll)
!            kkk=kkk+1
!            yarr(kkk)=df
!         enddo
!      endif
!   enddo
!   write(*,17)iph,ics,(yarr(j),j=1,kkk)
!17 format('Default cons: ',2i3,5(1pe12.4))
!   call set_constitution(iph,ics,yarr,qq,ceq)
!1000 continue
!   return
! end subroutine set_mean_constitution
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
   integer lokph,lokcs,ll,ml,ic,loksp,jl,nsl,locva
   double precision charge,spat,asite,bsite,badd,yz,yva,sumat,asum,bsum
!   double precision charge1,bion1,ionsites(2)
   double precision charge1,bion1
   TYPE(gtp_fraction_set), pointer :: disrec
   logical ionicliq
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
   ionicliq=btest(phlista(lokph)%status1,PHIONLIQ)
   if(ionicliq) then
! default values of i2slx
      phlista(lokph)%i2slx(1)=phlista(lokph)%tnooffr+1
      phlista(lokph)%i2slx(2)=phlista(lokph)%tnooffr+1
      yva=zero
      locva=0
   endif
!----
   if(ocv()) write(*,8)'25Ay:',iph,ics,&
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
!      if(ocv()) write(*,*)'25A 2: ',ionicliq
      ic=0
      do ll=1,phlista(lokph)%noofsubl
!         write(*,*)'25A sumy 2: ',ll,ic,phlista(lokph)%noofsubl
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
!         write(*,13)'25A y: ',ll,ic,asite,bmpymin,&
!              (ceq%phase_varres(lokcs)%yfr(ic+ml),&
!              ml=1,phlista(lokph)%nooffr(ll))
13       format(a,2i2,2(1pe12.4),1x,4(1pe12.4))
         ic=ic+phlista(lokph)%nooffr(ll)
      enddo
!--------
      ll=1; ml=0; asum=zero; bsum=zero; charge=zero
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
!         if(ocv()) write(*,*)'25A 3: ',ic,yz
         notva: if(btest(ceq%phase_varres(lokcs)%constat(ic),CONVA)) then
! i2slx(1) should be set to the index of vacancies (if any)
            if(ionicliq) phlista(lokph)%i2slx(1)=ic
            locva=ic
            yva=yz
         else
! sum charge and for constituents with several atoms spat sum number of atoms
            loksp=phlista(lokph)%constitlist(ic)
            charge=charge+bsite*yz*splista(loksp)%charge
! derivates of sites for ionic liquid model
!            if(ocv()) write(*,*)'25A 4: ',loksp,charge
            if(ionicliq) then
               ceq%phase_varres(lokcs)%dpqdy(ic)=abs(splista(loksp)%charge)
!               if(ocv()) write(*,*)'25A dpqdy:     ',&
!                    ic,abs(splista(loksp)%charge)
! i2slx(2) should be set to the index of the first neutral (if any)
               if(splista(loksp)%charge.eq.zero .and.&
                    phlista(lokph)%i2slx(2).gt.ic) &
                    phlista(lokph)%i2slx(2)=ic
            endif
            badd=badd+bsite*yz*splista(loksp)%mass
!            write(*,56)'setcon: ',iph,loksp,splista(loksp)%mass,yz,badd
56          format(a,2i3,3(1pe12.4))
            sumat=zero
! This is not adopted for other components than the elements
            do jl=1,splista(loksp)%noofel
               sumat=sumat+splista(loksp)%stoichiometry(jl)
            enddo
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
!         if(ocv()) write(*,*)'25A 5: ',ml
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
! initiate vacancy and neutral indices beyond last index (already done??)
                  phlista(lokph)%i2slx=phlista(lokph)%tnooffr+1
               elseif(ll.eq.2) then
! P=\sum_j (-v_j)y_j + Qy_Va. Note charge is total charge and valences 
! on 2nd sublattice is negative
! Now we know number of sites on sublattice 1, update asum and bsum
                  sumat=-charge+charge1*yva
                  ceq%phase_varres(lokcs)%sites(1)=sumat
!                  write(*,*)'Ionic 1: ',ceq%phase_varres(lokcs)%sites(1)
                  asum=asum*sumat
                  bsum=bion1*sumat
                  charge=zero
!                  write(*,88)'ionliq: ',ll,badd,bion1,bsum,sumat,yva
               else
                  write(*,*)'Ionic liquid must have two sublattices'
                  gx%bmperr=7777; goto 1000
               endif
            endif ionliq
! note: for ionic liquid previous values of asum and bsum are updated 
! when fractions in sublattice 2 have been set
            asum=asum+asite*spat
            bsum=bsum+badd
!            write(*,33)'25A g:',lokcs,ll,asum,asite,spat
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
!      write(*,33)'25A h:',lokcs,ll,asum,asite,spat
   endif nosuscon
! save charge, number of moles and mass of real atoms per formula unit
!   write(*,33)'25A i:',lokcs,0,asum,asite,spat
   ceq%phase_varres(lokcs)%netcharge=charge
   ceq%phase_varres(lokcs)%abnorm(1)=asum
   ceq%phase_varres(lokcs)%abnorm(2)=bsum
   if(ionicliq .and. locva.gt.0) then
! the ionic liquid vacancy charge is the number of sites on second subl.
      ceq%phase_varres(lokcs)%dpqdy(locva)=ceq%phase_varres(lokcs)%sites(2)
!      if(ocv()) write(*,*)'25A dpqdy(va): ',&
!           locva,ceq%phase_varres(lokcs)%sites(2)
   endif
!   if(ionicliq) then
!      write(*,301)'25A xsc:',lokcs,asum,bsum,ceq%phase_varres(lokcs)%sites,&
!           charge1
!301 format(a,i3,6(1pe12.4))
!      write(*,301)'25A y:  ',ic,ceq%phase_varres(lokcs)%yfr
!   endif
   qq(1)=asum
   qq(2)=charge
   qq(3)=bsum
! set disordered fractions if any
   if(btest(phlista(lokph)%status1,phmfs)) then
!now set disordered fractions if any
      call calc_disfrac(lokph,lokcs,ceq)
      if(gx%bmperr.ne.0) goto 1000
   endif
314 format(a,8F8.3)
1000 continue
!   if(ionicliq) write(*,*)'25A s_c: ',phlista(lokph)%i2slx
   return
 end subroutine set_constitution

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

