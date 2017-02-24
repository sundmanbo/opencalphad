!
! gtp3B included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     6. Enter data
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_element(symb,name,refstate,mass,h298,s298)
! Creates an element record after checks.
! symb: character*2, symbol (it can be a single character like H or V)
! name: character, free text name of the element
! refstate: character, free text name of reference state.
! mass: double, mass of element in g/mol
! h298: double, enthalpy difference between 0 and 298.14 K
! s298: double, entropy at 298.15 K
   implicit none
   CHARACTER*(*) symb,name,refstate
   DOUBLE PRECISION mass,h298,s298
!\end{verbatim}
   CHARACTER symb2*2,symb24*24
   integer knr(1),jl,jjj,kkk,nsl,loksp,lokph,nycomp
   double precision stoik(1)
   character ch1*1,model*24,phname*24,const(1)*24
   if(.not.allowenter(1)) then
      gx%bmperr=4125
      goto 1000
   endif
! check input data
100 continue
   call capson(symb)
   if(ucletter(symb(1:1))) then
      if(len(symb).ge.2) then
         if(ucletter(symb(2:2)) .or. symb(2:2).eq.' ') then
            goto 200
         endif
      else
         goto 200
      endif
   endif
! element name error, must be only letters (except /- already entered)
!   write(6,*)'new element not allowed ',symb,gx%bmperr
   gx%bmperr=4033
   goto 1000
200 continue
! check element not already entered
   symb2=symb(1:2)
!    write(*,202)'3B new element 1: ',symb,symb2
202 format(a,'"',a,'"',a,'"')
   reallynew: do jl=0,noofel
      if(symb2.eq.ellista(jl)%symbol) then
         gx%bmperr=4034
         goto 1000
      endif
   enddo reallynew
! element name is not really needed but must start with letter
!    write(*,12)symb,name,refstate,mass,h298,s298
!12  format('3B new_el: "',a,'"',a,'"',a,'"',3(1PE12.4))
   call capson(name)
   if(name(1:1).ne.' ') then
! allow empty element state
      if(.not.ucletter(name(1:1))) then
         gx%bmperr=4035
         goto 1000
      endif
   endif
300 continue
! reference state must start with letter, no other check
   call capson(refstate)
   if(refstate(1:1).ne.' ') then
! allow empty reference state
      if(.not.ucletter(refstate(1:1))) then
! error here when 1/2_MOLE_O2(G) etc ....
         refstate='GAS_'//refstate
!         gx%bmperr=4036
!         goto 1000
      endif
   endif
400 continue
! mass, h298-h0 and s298  must not be negative
   if(mass.lt.zero) then
      gx%bmperr=4037
      goto 1000
   endif
   if(h298.lt.zero) then
      gx%bmperr=4038
      goto 1000
   endif
   if(s298.lt.zero) then
      gx%bmperr=4039
      goto 1000
   endif
! All OK, increment noofel and store values in record noofel
   noofel=noofel+1
   if(noofel.gt.maxel) then
      gx%bmperr=4040
      goto 1000
   endif
! ensure that symbol has no strange characters
!    write(*,202)'3B new element 1B: ',symb,symb2
   ellista(noofel)%symbol='  '
   ellista(noofel)%symbol=symb
   ellista(noofel)%name=name
   ellista(noofel)%ref_state=refstate
   ellista(noofel)%mass=mass
   ellista(noofel)%h298_h0=h298
   ellista(noofel)%s298=s298
   ellista(noofel)%status=0
   ellista(noofel)%alphaindex=noofel
! value 0 is H298, 1 H0, 2 G
   ellista(noofel)%refstatesymbol=0
! Now create corresponding species
   noofsp=noofsp+1
   if(noofel.gt.maxsp) then
      gx%bmperr=4041
      goto 1000
   endif
   ellista(noofel)%splink=noofsp
!   write(*,202)'3B new element 1C: ',symb,symb2
   symb24=' '
   symb24=symb2
!    write(*,77)symb,symb2,symb24
!77  format('3B new element 77: ',a,'"',a,'"',a,'"')
   splista(noofsp)%symbol=symb24
   splista(noofsp)%mass=mass
   splista(noofsp)%charge=zero
   splista(noofsp)%status=0
   splista(noofsp)%status=ibset(splista(noofsp)%status,SPEL)
   splista(noofsp)%alphaindex=noofsp
   splista(noofsp)%noofel=1
! allocate
   allocate(splista(noofsp)%ellinks(1))
   allocate(splista(noofsp)%stoichiometry(1))
   splista(noofsp)%ellinks(1)=noofel
   splista(noofsp)%stoichiometry(1)=one
! return with error code 0 i.e. no error
!    gx%bmperr=0
! rearrange ELEMENTS and SPECIES to maintain these in alphabetical order
   elements(noofel)=noofel
   call alphaelorder
   species(noofsp)=noofsp
   call alphasporder
! As this is an element add the species to the component list of firsteq
!------------------------------------------------
! Beware that the alphabetical order may have changed. jjj used later
   jjj=ellista(noofel)%alphaindex
   if(jjj.lt.noofel) then
!      write(*,*)'3B TDB MUST HAVE ELEMENTS IN ALPHABETICAL ORDER!',jjj,noofel
      do kkk=noofel,jjj+1,-1
         firsteq%complist(kkk)%splink=firsteq%complist(kkk-1)%splink
         firsteq%complist(kkk)%phlink=firsteq%complist(kkk-1)%phlink
         firsteq%complist(kkk)%refstate=firsteq%complist(kkk-1)%refstate
         firsteq%complist(kkk)%tpref=firsteq%complist(kkk-1)%tpref
         firsteq%complist(kkk)%mass=firsteq%complist(kkk-1)%mass
      enddo
   else
      jjj=noofel
   endif
! %splink is location of species
   firsteq%complist(jjj)%splink=noofsp
   firsteq%complist(jjj)%phlink=0
! do not copy element reference state name here
   firsteq%complist(jjj)%refstate='SER (default)'
   firsteq%complist(jjj)%tpref(1)=2.9815D2
   firsteq%complist(jjj)%tpref(2)=1.0D5
! copy mass of component from species record
   firsteq%complist(jjj)%mass=mass
! check
!   call compmassbug(firsteq)
! NOTE jjj is used below when adding this element to reference phase
! also set the stoichiometry matrix, just the diagonal.  Also the inverse
   firsteq%compstoi(noofel,noofel)=one
   firsteq%invcompstoi(noofel,noofel)=one
!   write(*,*)'3B new_el: ',noofel,name,symb24
   nycomp=noofel
   if(noofel.eq.1) then
! create reference phase with index 0
!       phname='ELEMENT_REFERENCE_PHASE '
      phname='SELECT_ELEMENT_REFERENCE'
      nsl=1
      knr(1)=1
!      const(1)=name
      const(1)=symb24
      stoik(1)=one
      model='NON_MIXING'
      ch1='Z'
      call enter_phase(phname,nsl,knr,const,stoik,model,ch1)
      if(gx%bmperr.ne.0) goto 1000
! set phase hidden as it should never be included in calculations
      lokph=0
      phlista(lokph)%status1=ibset(phlista(lokph)%status1,phhid)
! add all additions ??
   else
! Add the element to the reference phase (phase 0) by extending the
! constituent list (and many other arrays)
      loksp=firsteq%complist(jjj)%splink
      call add_to_reference_phase(loksp)
      if(gx%bmperr.ne.0) goto 1000
   endif
   if(noofel.gt.0) then
! clear the nodata bit
      globaldata%status=ibclr(globaldata%status,GSNODATA)
   endif
!    if(gx%bmperr.ne.0) goto 1000
1000 continue
!    write(*,*)'3B created new species: ',noofsp,splista(noofsp)%symbol
   return
 END subroutine enter_element

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_species(symb,noelx,ellist,stoik)
! creates a new species
! symb: character*24, name of species, often equal to stoichiometric formula
! noelx: integer, number of elements in stoichiometric formula (incl charge)
! ellist: character array, element names (electron is /-)
! stoik: double array, must be positive except for electron.
   implicit none
   character symb*(*),ellist(*)*(*)
   integer noelx
   double precision stoik(*)
!\end{verbatim}
   double precision mass,charge
   integer elindex(10)
   integer loksp,noelxx,jl,jk
   if(.not.allowenter(1)) then
      gx%bmperr=4125
      goto 1000
   endif
   call capson(symb)
!   write(*,*)'3B Entering ',symb,noelx
   if(.not.ucletter(symb(1:1))) then
      gx%bmperr=4044
      goto 1000
   endif
   if(noelx.le.0 .or. noelx.gt.10) then
      gx%bmperr=4045
      goto 1000
   endif
! check if there is a period "." in the species, that is a common error!
   if(index(symb,'.').gt.0) then
      gx%bmperr=4044; goto 1000
   endif
! check symb is unique
!   call find_species_record(symb,loksp)
   call find_species_record_noabbr(symb,loksp)
   if(gx%bmperr.eq.0) then
! strange error reading cadarache database
      do jl=1,noofsp
         write(*,*)'3B entered species ',jl,splista(jl)%symbol
      enddo
      gx%bmperr=4049; goto 1000
   endif
   mass=zero
   charge=zero
   noelxx=noelx
!    write(*,*)'3B enter_species 1A: ',noelx
   checkel: do jl=1,noelx
      loopel: do jk=-1,noofel
         if(ellist(jl).eq.ellista(jk)%symbol) goto 200
      enddo loopel
! an unknown element
      gx%bmperr=4046
      goto 1000
200    continue
      elindex(jl)=jk
      if(jk.ge.0) then
         if(stoik(jl).lt.zero) then
            gx%bmperr=4047
            goto 1000
         else
            mass=mass+stoik(jl)*ellista(jk)%mass
         endif
      else
! this is the electron, save negative of stoick as charge negative
! the electron is not counted as "element" when storing
         charge=-stoik(jl)
         noelxx=noelxx-1
         if(jl.ne.noelx) then
! this must be the last element .... otherwise problem storing stoik
            gx%bmperr=4048
            goto 1000
         endif
      endif
!     write(6,*)'enter_species 2: ',symb,jl,mass,charge
   enddo checkel
   noofsp=noofsp+1
   if(noofsp.gt.maxsp) then
      gx%bmperr=4125
      goto 1000
   endif
! store species data
   splista(noofsp)%symbol=symb
   splista(noofsp)%mass=mass
   splista(noofsp)%charge=charge
   splista(noofsp)%alphaindex=noofsp
   splista(noofsp)%noofel=noelxx
   splista(noofsp)%status=0
! I am not sure what this will be used for ...
   splista(noofsp)%extra=one
   if(charge.ne.zero) then
      splista(noofsp)%status=ibset(splista(noofsp)%status,SPION)
   endif
! allocate
   allocate(splista(noofsp)%ellinks(noelxx))
   allocate(splista(noofsp)%stoichiometry(noelxx))
   loop2: do jl=1,noelxx
      splista(noofsp)%ellinks(jl)=elindex(jl)
      splista(noofsp)%stoichiometry(jl)=stoik(jl)
   enddo loop2
! return with no error
   gx%bmperr=0
! add species last and rearrange
   species(noofsp)=noofsp
   call alphasporder
! error: continue would be a nice use of non-digit labels ....
1000 continue
   return
 END subroutine enter_species

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_phase(name,nsl,knr,const,sites,model,phtype)
! creates the data structure for a new phase
! name: character*24, name of phase
! nsl: integer, number of sublattices (range 1-9)
! knr: integer array, number of constituents in each sublattice
! const: character array, constituent (species) names in sequential order
! sites: double array, number of sites on the sublattices
! model: character, free text
! phtype: character*1, specifies G for gas, L for liquid
   implicit none
   character name*(*),model*(*),phtype*(*)
   integer nsl
   integer, dimension(*) :: knr
   double precision, dimension(*) :: sites
   character, dimension(*) :: const*(*)
!\end{verbatim}
   type(gtp_phase_add), pointer :: addrec
   character ch1*1
   double precision formalunits,endch
   integer kconlok(maxconst),kalpha(maxconst),iord(maxconst),klok(maxconst)
   integer iva(maxconst),endm(maxsubl),endm0(maxsubl+1)
   logical externalchargebalance,tupix
   integer iph,kkk,lokph,ll,nk,jl,jk,mm,lokcs,nkk,nyfas,loksp,tuple
!   write(*,*)'3B enter enter_phase: ',trim(name),' ',trim(model)
   if(.not.allowenter(2)) then
      gx%bmperr=4125
      goto 1000
   endif
! check input
   call capson(name)
!   if(.not.ucletter(name)) then
   if(.not.proper_symbol_name(name,0)) then
      write(*,*)'3B Error for phase name: ',name(1:min(24,len(name)))
      gx%bmperr=4053; goto 1000
   endif
! name unique?
   call find_phase_by_name_exact(name,iph,kkk)
!   write(6,*)'new phase 1A ',name,nsl,gx%bmperr,const(1)
   if(gx%bmperr.eq.0) then
! if phase found then error as name not unique ... but check explicitly
      lokph=phases(iph)
      if(name.eq.phlista(lokph)%name) then
         gx%bmperr=4054
         goto 1000
      endif
! name was not exactly the same, accept this phase name also
   else
      gx%bmperr=0
   endif
   if(nsl.lt.1 .or. nsl.gt.maxsubl) then
      gx%bmperr=4056
      goto 1000
   endif
   site1: do ll=1,nsl
      if(sites(ll).le.zero) then
!        write(6,*)' new phase 1B: ',name,ll,nsl,sites(ll)
         gx%bmperr=4057
         goto 1000
      endif
   enddo site1
   nk=0
   knrtest: do ll=1,nsl
      if(knr(ll).lt.1 .or. knr(ll).gt.maxconst) then
         write(*,*)'3B enter phase error:',ll,knr(ll),maxconst
         gx%bmperr=4058; goto 1000
      endif
      if(ll.ge.2 .and. knr(ll).gt.maxcons2) then
         gx%bmperr=4059; goto 1000
      endif
      nk=nk+knr(ll)
   enddo knrtest
   nkk=nk
!  write(6,*)' enter_phase 3: ',name,nsl,nkk,noofsp
! check constituents exists as species
   externalchargebalance=.false.
   constest: do jl=1,nkk
      if(jl.eq.1 .and. model(1:13).eq.'IONIC_LIQUID ') then
! in this case * is allowed on first sublattice!!
         if(const(1)(1:2).eq.'* ') then
            kalpha(jl)=-99
            kconlok(jl)=-99
            cycle constest
         endif
      endif
      call capson(const(jl))
!      write(6,297)' enter_phase constituent: ',jl,const(jl),nkk
      findspecies: do jk=1,noofsp
         if(const(jl).eq.splista(jk)%symbol) then
!            write(*,*)'3B at new 300: ',noofsp,jk,const(jl)
            goto 300
         endif
      enddo findspecies
!      write(6,297)' enter_phase constituent error: ',jl,const(jl),jk,nkk
297 format(a,i3,'>',A,'<',2i3)
      write(kou,*)'Unknown constituent, name must be exact: ',trim(const(jl))
      gx%bmperr=4051
      goto 1000
! found species,
300   continue
! check for duplicates in same sublattice
      kalpha(jl)=splista(jk)%alphaindex
      ll=1
      mm=1
      nk=knr(1)
310   continue
      if(jl.gt.nk) then
         if(ll.lt.nsl) then
            ll=ll+1
            mm=nk+1
            nk=nk+knr(ll)
            goto 310
         else
            write(*,*)'3B Impossible: constituent index outside range!'
            gx%bmperr=4257; goto 1000
         endif
      else
         do mm=mm,jl-1
!            write(*,314)mm,jl,kalpha(mm),kalpha(jl),&
!                 const(jl)(1:len_trim(const(jl))),name(1:len_trim(name))
314         format('3B Species: ',4i4,' "',a,'" in ',a)
            if(kalpha(mm).eq.kalpha(jl)) then
               write(*,315)name(1:len_trim(name)),&
                    const(jl)(1:len_trim(const(jl))),ll
315            format(' *** Error, the ',a,' phase has constituent ',a,&
                    ' twice in sublattice ',i2)
               gx%bmperr=4258; goto 1000
            endif
         enddo
      endif
      kconlok(jl)=jk
!     write(6,73)' enter_phase 4B: ',jl,const(jl),jk,kconlok(jl),kalpha(jl)
!73   format(A,i3,1x,A6,3I3)
! mark that PHEXCB bit must be set if species has a charge 
      if(splista(jk)%charge.ne.zero) then
         externalchargebalance=.true.
      endif
   enddo constest
! reserve a new phase record and store data there and in other records
! the first phase entered is the reference phase created by init_gtp
   if(noofph.eq.0 .and. phtype(1:1).eq.'Z') then
! phtyp=Z is the reference phase
      nyfas=0
   else
! sort the phase in alphabetical order but always gas (if any) first
! then liquids specified by the phtype letter (G, L, etc)
      noofph=noofph+1
!      if(nyfas.gt.size(phlista)) then
      if(noofph.gt.size(phlista)) then
!         write(*,*)'3B Too many phases: ',noofph
         gx%bmperr=4259; goto 1000
      endif
      nyfas=noofph
   endif
   phlista(nyfas)%name=name
   phlista(nyfas)%status1=0
   ionliq: if(model(1:13).eq.'IONIC_LIQUID ') then
! the external charge balance set above, not needed
!      write(*,*)'3B  *** ionic liquid entered!!!'
      externalchargebalance=.FALSE.
! ionic liquid may have phtype='Y', change that to L
      if(phtype(1:1).eq.'Y') phtype(1:1)='L'
      if(nsl.ne.2) then
! if entered with only one sublattice then no cations and only neutrals!!
         write(*,*)'3B Ionic liquid must have 2 sublattices'
         gx%bmperr=4255; goto 1000
      endif
      phlista(nyfas)%status1=ibset(phlista(nyfas)%status1,PHIONLIQ)
! constituents in ionic liquid must be sorted in a special way
      call sort_ionliqconst(lokph,0,knr,kconlok,klok)
      if(gx%bmperr.ne.0) goto 1000
   else ! else link is for all other phases except ionic liquid
! sort the constituents in each sublattice according to alphaspindex
!  write(6,70)5,(kalpha(i),i=1,nkk)
!  write(6,70)5,(kconlok(i),i=1,nkk)
!70    format('enter_phase ',I2,': ',20I3)
      nk=1
      sort1: do ll=1,nsl
         call sortin(kalpha(nk),knr(ll),iord(nk))
         if(buperr.ne.0) then
            gx%bmperr=buperr
            goto 1000
         endif
! iord(nk+1:nk+knr(ll)) has numbers 1..knr(ll), add on nk-1 to these
! to be in parity with index of kalpha(nk+1:nk+knr(ll))
         adjust: do mm=0,knr(ll)-1
            iord(nk+mm)=iord(nk+mm)+nk-1
         enddo adjust
         nk=nk+knr(ll)
      enddo sort1
!  write(6,70)6,(kalpha(i),i=1,nkk)
!  write(6,70)6,(kconlok(iord(i)),i=1,nkk)
! in constituent record store kconlok(iord(i))
! verify we can find species name ...
!  test7: do kk=1,nkk
!    write(6,71)kk,iord(kk),kconlok(iord(kk)),splista(kconlok(iord(kk)))%symbol
!71 format('enter_phase 7: ',3I3,1x,A)
!  enddo test7
      do jl=1,nkk
         klok(jl)=kconlok(iord(jl))
      enddo
   endif ionliq
!----------------------------------------
!  write(6,79)8,name,(klok(kk),kk=1,nkk)
79    format('enter_phase ',I2,': ',A6,10I3)
   ch1=phtype(1:1)
   call capson(ch1)
! sort the phase in alphabetical but order but first gas, then liquid etc
! legal values of ch1 is G, L, S and C (gas, liquid, solution, compound)
!   write(*,*)'3B phase byte: ',ch1
   if(ch1.eq.'G') then
      phlista(nyfas)%status1=ibset(phlista(nyfas)%status1,PHGAS)
      model='ideal'
   elseif(ch1.eq.'L') then
      phlista(nyfas)%status1=ibset(phlista(nyfas)%status1,PHLIQ)
   endif
! Handle option F and B for permutations
   if(ch1.eq.'F') then
!      write(*,*)'3B Setting PHFORD bit'
      phlista(nyfas)%status1=ibset(phlista(nyfas)%status1,PHFORD)
!      call set_phase_status_bit(lokph,PHFORD)
   elseif(ch1.eq.'B') then
!      write(*,*)'3B Setting PHBORD bit'
      phlista(nyfas)%status1=ibset(phlista(nyfas)%status1,PHBORD)
!      call set_phase_status_bit(lokph,PHBORD)
   endif
! :I is used by TC to indicate charge balance needed, ignore
   if(ch1.eq.' ' .or. ch1.eq.'I') ch1='S'
!   ch1='S'
   phlista(nyfas)%phletter=ch1
   phlista(nyfas)%models=model
!   if(nyfas.eq.0) then
!      continue
!   else
   if(nyfas.gt.0) then
      call alphaphorder(tuple)
      phlista(nyfas)%nooffs=1
   else
! uninitiated below for reference phase
      tuple=0
   endif
   phlista(nyfas)%noofsubl=nsl
   allocate(phlista(nyfas)%nooffr(nsl))
! sites stored in phase_varres
!   allocate(phlista(nyfas)%sites(nsl))
   formalunits=zero
   do ll=1,nsl
      phlista(nyfas)%nooffr(ll)=knr(ll)
      formalunits=formalunits+sites(ll)
   enddo
!  write(*,*)'3B enter_phase 8x: ',nyfas,nkk
   phlista(nyfas)%tnooffr=nkk
!  write(*,*)'3B enter_phase 8y: ',nyfas,phlista(nyfas)%tnooffr
! create consituent record
   call create_constitlist(phlista(nyfas)%constitlist,nkk,klok)
! in phase_varres we will indicate the VA constituent, indicate in iva
   valoop: do jl=1,nkk
      iva(jl)=0
      loksp=phlista(nyfas)%constitlist(jl)
      if(loksp.gt.0) then
! ionic liquid can have a wildcard */-99 as constituent in first sublattice
         if(btest(splista(loksp)%status,SPVA)) iva(jl)=ibset(iva(jl),CONVA)
      endif
   enddo valoop
!    write(*,32)'3B phase14A: ',nyfas,(phlista(nyfas)%constitlist(iz),iz=1,nkk)
32  format(a,i3,50(i3))
!    write(*,33)nkk,(iva(i),i=1,nkk)
!33 format('3B enter_phase 14B: ',i3,2x,10i3)
!   nprop=10
!   write(*,*)'3B enter_phase: ',lokcs,name
   call create_parrecords(nyfas,lokcs,nsl,nkk,maxcalcprop,iva,firsteq)
!   write(*,*)'3B enter_phase 15: ',nyfas,lokcs
   if(gx%bmperr.ne.0) goto 1000
! zero array of pointer to phase_varres record, then set first
   phlista(nyfas)%linktocs=0
   phlista(nyfas)%linktocs(1)=lokcs
   phlista(nyfas)%noofcs=1
   firsteq%phase_varres(lokcs)%phlink=nyfas
   firsteq%phase_varres(lokcs)%prefix=' '
   firsteq%phase_varres(lokcs)%suffix=' '
! Initiated to total number of sites, will be updated in set_condition
   firsteq%phase_varres(lokcs)%abnorm(1)=formalunits
! ncc no longer part of this record
!   firsteq%phase_varres%ncc=nkk
! zero the phstate (means entered and not known (unknown) if stable)
   firsteq%phase_varres(lokcs)%phstate=0
! sites must be stored in phase_varres
   do ll=1,nsl
      firsteq%phase_varres(lokcs)%sites(ll)=sites(ll)
   enddo
! make sure status word and some other links are set
   firsteq%phase_varres(lokcs)%status2=0
   firsteq%phase_varres(lokcs)%phtupx=tuple
! set link to lokcs in phase tuple!
!   phasetuple(tuple)%lokvares=lokcs
!   write(*,*)'3B new phase tuple: ',nyfas,lokcs,tuple
! If one has made NEW the links are not always zero
! set some phase bits (PHGAS and PHLIQ set above)
! external charge balance etc.
   if(externalchargebalance) then
! do not set this if all endmembers have zero charge  m2o3(Ce+3,La+3)2(O-2)3
      jl=1
      endch=zero
      do ll=1,nsl
         endm0(ll)=jl
         endm(ll)=jl
         jk=phlista(nyfas)%constitlist(jl)
         endch=endch+splista(jk)%charge*sites(ll)
         jl=jl+phlista(nyfas)%nooffr(ll)
      enddo
      endm0(nsl+1)=phlista(nyfas)%tnooffr+1
400   continue
!      write(*,405)'3B netcharge ',trim(name),endch,(endm(ll),ll=1,nsl)
!405   format(a,a,1pe12.4,9i4)
      if(abs(endch).gt.1.0D-6) goto 420
      ll=nsl
410   continue
      if(endm(ll).lt.endm0(ll+1)-1) then
         jk=phlista(nyfas)%constitlist(endm(ll))
         endch=endch-splista(jk)%charge*sites(ll)
         endm(ll)=endm(ll)+1
         jk=phlista(nyfas)%constitlist(endm(ll))
         endch=endch+splista(jk)%charge*sites(ll)
         goto 400
      elseif(ll.gt.1) then
         jk=phlista(nyfas)%constitlist(endm(ll))
         endch=endch-splista(jk)%charge*sites(ll)
         endm(ll)=endm0(ll)
         jk=phlista(nyfas)%constitlist(endm(ll))
         endch=endch+splista(jk)%charge*sites(ll)
         ll=ll-1
         goto 410
      endif
!      write(*,*)'3B charge balance not needed for ',trim(name)
      goto 430
! jump here if any endmember has a net charge
420   continue
!      write(*,*)'3B charge balance needed for ',trim(name)
      phlista(nyfas)%status1=ibset(phlista(nyfas)%status1,PHEXCB)
! jump here if all neutral
430   continue
   endif
! set net charge to zero
   firsteq%phase_varres(lokcs)%netcharge=zero
   if(nsl.eq.1) then
! if no sublattices set ideal bit.  Will be cleared if excess parameter entered
      phlista(nyfas)%status1=ibset(phlista(nyfas)%status1,PHID)
   endif
   if(nkk.eq.nsl) then
! as many constiuents as sublattice
      phlista(nyfas)%status1=ibset(phlista(nyfas)%status1,PHNOCV)
   endif
! nullify links
   nullify(phlista(nyfas)%additions)
   nullify(phlista(nyfas)%ordered)
   nullify(phlista(nyfas)%disordered)
! initiate phcs, the phase composition set counter for nyfas redundant ??
! (not for reference phase 0) 
!   if(nyfas.gt.0) phcs(nyfas)=1
   if(noofph.gt.0) then
! clear the nophase bit
      globaldata%status=ibclr(globaldata%status,GSNOPHASE)
!---------------------- new code to generate phase tuple array here
! NOTE nooftuples updated in alphaphorder ... for old times sake
      do ll=1,noofph
! this is index in phlista
!         phasetuple(ll)%phaseix=phases(ll)
         phasetuple(ll)%lokph=phases(ll)
         phasetuple(ll)%compset=1
! this is alphabetical index
         phasetuple(ll)%ixphase=ll
! this is link to higher tuple of same phase
         phasetuple(ll)%nextcs=0
! this is the link to phase tuple from the phase
         jl=phlista(phases(ll))%linktocs(1)
         firsteq%phase_varres(jl)%phtupx=ll
         phasetuple(ll)%lokvares=jl
      enddo
!---------------------- new code end
   endif
! always enter volume model1, nyfas is lokph, use alphabetical index
   if(nyfas.gt.0) then
!      write(*,*)'3B enter_phase adding volume model: ',trim(name),nyfas
      call add_addrecord(nyfas,' ',volmod1)
   endif
1000 continue
   return
 END subroutine enter_phase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine sort_ionliqconst(lokph,mode,knr,kconlok,klok)
! sorts constituents in ionic liquid, both when entering phase
! and decoding parameter constituents
! order: 1st sublattice only cations
! 2nd: anions, VA, neutrals
! mode=0 at enter phase, wildcard ok in 1st sublattice if neiher anions nor Va
! mode=1 at enter parameter (wildcard allowed, i.e. some kconlok(i)=-1)
! some  parameters not allowed, L(ion,A+:B,C), must be L(ion,*:B,C), check!
   implicit none
   integer lokph,knr(*),kconlok(*),klok(*),mode
!\end{verbatim}
   integer nk,jl,jk,mm,kkk,ionva,byte
   integer, dimension(:), allocatable :: kalpha,iord,iva,anion
!
   allocate(kalpha(knr(1)+knr(2)))
   allocate(iord(knr(1)+knr(2)))
   allocate(iva(knr(1)+knr(2)))
   allocate(anion(knr(1)+knr(2)))
! check1: constituents in sublattice 1 must all have positive charge
!   if(mode.eq.1) then
!      write(*,17)'3B sl2: ',knr(1),knr(2),(kconlok(mm),mm=1,knr(1)+knr(2))
!17    format(a,2i3,2x,10i3)
!   endif
   do nk=1,knr(1)
      if(kconlok(nk).lt.0) then
! wildcard give index -99. If mode=0 more checks later
         kalpha(nk)=-99
      elseif(splista(kconlok(nk))%charge.le.zero) then
         write(*,*)'3B In ionic_liquid only cations on first sublattice'
         gx%bmperr=4260; goto 1000
      else
         kalpha(nk)=splista(kconlok(nk))%alphaindex
      endif
   enddo
!   write(*,69)'3B In 1: ',knr(1),(kconlok(mm),mm=1,knr(1))
   if(knr(1).gt.1) then
      call sortin(kalpha,knr(1),iord)
      if(buperr.ne.0) then
         gx%bmperr=buperr
         goto 1000
      endif
      if(mode.eq.0 .and. kalpha(1).lt.0) then
! when entering phase a single wildcard allowed in first sublattice
         write(*,*)'3B Illegal parameter with wildcard mixed with cations'
         gx%bmperr=4261; goto 1000
      endif
      do jl=1,knr(1)
         klok(jl)=kconlok(iord(jl))
      enddo
   else
      klok(1)=kconlok(1)
   endif
!   write(*,69)'3B 1st:  ',knr(1),(kalpha(mm),mm=1,knr(1))
! check2: constituents in sublattice 1 must be ANIONS, VA and NEUTRALS
! in that order
   kkk=knr(1)
   jl=0
   jk=0
   ionva=0
   do nk=1,knr(2)
      if(mode.eq.0 .and. kconlok(nk+kkk).lt.0) then
! when entering phase no wildcards allowed in second sublattice
         write(*,*)'3B You cannot enter phase with wildcard on 2nd sublattice'
         gx%bmperr=4262; goto 1000
      elseif(kconlok(nk+kkk).lt.0) then
! wildcard, treat as anion ?? DO NOT ALLOW, what stoichiometry??
         write(*,*)'3B Ionic_liq parameter with wildcard on 2nd sublat. illegal'
         gx%bmperr=4262; goto 1000
!         jk=jk+1
!         anion(jk)=nk
      elseif(splista(kconlok(nk+kkk))%charge.gt.zero) then
         write(*,*)'3B No cations allowed on second sublattice'
         gx%bmperr=4263; goto 1000
      elseif(btest(splista(kconlok(nk+kkk))%status,SPVA)) then
! this is the hypothetical vacancy
         ionva=nk
      elseif(splista(kconlok(nk+kkk))%charge.eq.zero) then
! neutral species allowed, use iva, must be sorted after all anions and Va
         jl=jl+1
         iva(jl)=nk
      else
! anion
         jk=jk+1
         anion(jk)=nk
      endif
   enddo
!   write(*,88)'3B at 1:  ',knr(2),(kconlok(knr(1)+mm),mm=1,knr(2))
88 format(a,i4,2x,20i3)
! There are jl neutrals and jk anions, if vacancies set it as jk+1
! if wildcard on first sublattice neither ainons nor Va allowed on 2nd
   if(klok(1).lt.0 .and. (jk.gt.0 .or. ionva.ne.0)) then
      write(*,*)'3B Only neutrals on second sublattice if wildcard on first'
      gx%bmperr=4264; goto 1000
   endif
   do nk=1,jk
      if(anion(nk).gt.nk) then
! shift the anion to position nk, kconlok must be updated
         if(ionva.eq.nk) then
            byte=kconlok(kkk+nk)
            kconlok(kkk+nk)=kconlok(kkk+anion(nk))
            ionva=anion(nk)
            kconlok(kkk+ionva)=byte
!            write(*,88)'3B byt 1: ',knr(2),(kconlok(knr(1)+mm),mm=1,knr(2))
         else
            do mm=1,jl
               if(iva(mm).eq.nk) exit
            enddo
            if(mm.gt.jl) stop 'big bug'
            byte=kconlok(kkk+nk)
            kconlok(kkk+nk)=kconlok(kkk+anion(nk))
            iva(mm)=anion(nk)
            kconlok(kkk+iva(mm))=byte
!            write(*,88)'3B byt 2: ',knr(2),(kconlok(knr(1)+mm),mm=1,knr(2))
         endif
         anion(nk)=nk
      endif
   enddo
!   write(*,88)'3B at 2:  ',knr(2),(kconlok(knr(1)+mm),mm=1,knr(2))
! now all ions should be in positions 1..jk.  Fix position of vacancy
! by moving neiutrals
   if(ionva.gt.jk+1) then
      byte=kconlok(kkk+jk+1)
      kconlok(kkk+jk+1)=kconlok(kkk+ionva)
      kconlok(kkk+ionva)=byte
      iva(ionva)=ionva
      ionva=jk+1
   endif
!   write(*,88)'3B at 3:  ',knr(2),(kconlok(knr(1)+mm),mm=1,knr(2))
!   write(*,69)'3B 2nda: ',jk,&
!        (splista(kconlok(kkk+anion(mm)))%alphaindex,mm=1,jk)
!   if(ionva.gt.0) &
!        write(*,69)'3B 2ndv: ',1,splista(kconlok(kkk+ionva))%alphaindex
!   write(*,69)'3B 2ndn: ',jl,&
!        (splista(kconlok(kkk+iva(mm)))%alphaindex,mm=1,jl)
69 format(a,i3,2x,10i3,i5,10i3)
   do mm=1,knr(2)
      if(kconlok(kkk+mm).lt.0) then
         kalpha(mm+kkk)=-99
      else
         kalpha(mm+kkk)=splista(kconlok(kkk+mm))%alphaindex
      endif
   enddo
   kkk=knr(1)+1
!   write(*,69)'3B 2ndx: ',knr(2),(kalpha(mm+kkk-1),mm=1,knr(2))
   if(jk.gt.1) then
!      write(*,69)'3B kalpha: ',jk,(kalpha(kkk+mm-1),mm=1,jk)
      call sortin(kalpha(kkk),jk,iord)
      if(buperr.ne.0) then
         gx%bmperr=buperr; goto 1000
      endif
!      write(*,69)'3B sort jk: ',jk,(iord(kkk+mm-1),mm=1,jk)
      do mm=1,jk
         klok(kkk+mm-1)=kconlok(kkk+iord(mm)-1)
      enddo
   elseif(jk.gt.0) then
      klok(kkk)=kconlok(kkk)
   endif
   kkk=kkk+jk
   if(ionva.gt.0) then
      klok(kkk)=kconlok(kkk)
      kkk=kkk+1
   endif
   if(jl.gt.1) then
      call sortin(kalpha(kkk),jl,iord)
      if(buperr.ne.0) then
         gx%bmperr=buperr; goto 1000
      endif
      do mm=1,jl
         klok(kkk+mm-1)=kconlok(kkk+iord(mm)-1)
      enddo
   elseif(jl.gt.0) then
      klok(kkk)=kconlok(kkk)
   endif
   if(mode.eq.1) then
! final check for parameters:
! if only neutrals on sublatice 2 no interaction allowed on sublattice 1
      if(jk.eq.0 .and. ionva.eq.0) then
         if(knr(1).gt.1) then
            write(*,*)'3B Illegal interaction parameter'
            gx%bmperr=4265; goto 1000
         else
! replace whatever constituent specified in sublattice 1 by wildcard
            klok(1)=-99
         endif
      endif
   endif
!   write(*,69)'3B al1: ',knr(1)+knr(2),&
!        (klok(mm),mm=1,knr(1)+knr(2))
!   write(*,69)'3B al2: ',knr(1)+knr(2),&
!        (splista(klok(mm))%alphaindex,mm=1,knr(1)+knr(2))
!----------------------------------------------------------
1000 continue
   return
 end subroutine sort_ionliqconst

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_composition_set(iph,prefix,suffix,icsno)
! adds a composition set to a phase.
! iph: integer, phase index
! prefix: character*4, optional prefix to original phase name
! suffix: character*4, optional suffix to original phase name
! icsno: integer, returned composition set index (value 2-9)
! ceq: pointer, to current gtp_equilibrium_data
!
! BEWARE this must be done in all equilibria (also during parallel processes)
! There may still be problems with equilibria saved during STEP and MAP
!
   implicit none
   integer iph,icsno
   character*(*) prefix,suffix
!   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! also update phasetuple array !! csfree,highcs
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer lokph,ncs,nsl,nkk,lokcs,lokcs1,nprop,lastcs,jl,nyttcs
   integer leq,nydis,tuple,nz,jz
   character*4 pfix,sfix
   integer iva(maxconst)
   TYPE(gtp_phase_varres), pointer :: peq,neq,ndeq
   logical once
!
!   write(*,*)'3B in enter_composition set',iph,phases(iph),nooftuples
   once=.TRUE.
   if(iph.le.0 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   endif
   lokph=phases(iph)
   ncs=phlista(lokph)%noofcs
   if(ncs.gt.8) then
! max 9 composition sets
      gx%bmperr=4092; goto 1000
   endif
   ceq=>firsteq
   icsno=ncs+1
! test if mmy is correct in all existing compsets
! OK here
!   do jl=1,ncs
!      lokcs=phlista(lokph)%linktocs(jl)
!      write(*,7)lokcs,firsteq%phase_varres(lokcs)%mmyfr
!7     format('3B mmy: ',i4,10F5.1)
!   enddo
! collect some data needed
   nsl=phlista(lokph)%noofsubl
   nkk=phlista(lokph)%tnooffr
   lokcs=phlista(lokph)%linktocs(phlista(lokph)%noofcs)
   lokcs1=lokcs
   nprop=ceq%phase_varres(lokcs)%nprop
   lastcs=phlista(lokph)%linktocs(phlista(lokph)%noofcs)
! one must set the VA bit in the constituent status array
   ivaloop: do jl=1,nkk
      iva(jl)=ceq%phase_varres(lastcs)%constat(jl)
   enddo ivaloop
! check that prefix is empty or start with a letter
   if(biglet(prefix(1:1)).ne.' ' .and. &
        (biglet(prefix(1:1)).lt.'A' .or. biglet(prefix(1:1)).gt.'Z')) then
      write(kou,*)'Prefix of composition set must start with a letter'
      gx%bmperr=4167; goto 1000
   endif
   if(biglet(suffix(1:1)).ne.' ' .and. &
        (biglet(suffix(1:1)).lt.'A' .or. biglet(suffix(1:1)).gt.'Z')) then
      write(kou,*)'Suffix of composition set must start with a letter'
      gx%bmperr=4167; goto 1000
   endif
!------------------------------------------------------------------
! begin threadprotected code >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! composition sets must be created in all equilibria
! note that indices to phase_varres same in all equilibria
! >>> beware not tested created composition sets with several equilibria 
! maybe this call can be replaced by a simple assignment????
! create_parrecord in GTP3G.F90 update csfree etc 
!   call create_parrecords(lokph,nyttcs,nsl,nkk,maxcalcprop,iva,ceq)
!   call create_parrecords(lokph,nyttcs,nsl,nkk,maxcalcprop,iva,firsteq)
   call create_parrecords(lokph,nyttcs,nsl,nkk,maxcalcprop,iva,firsteq)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3B add_cs: ',nyttcs
! add new tuple at the end and save tuple index
   tuple=nooftuples+1
!   phasetuple(tuple)%phaseix=phases(iph)
   phasetuple(tuple)%lokph=phases(iph)
   phasetuple(tuple)%compset=icsno
! New variables in phase tuple!, phase index and phase_varrres
   phasetuple(tuple)%ixphase=iph
   phasetuple(tuple)%lokvares=nyttcs
! nextcs is the index of next phasetuple for same phase
   leq=iph
! why upper bound error??
   do while(leq.le.nooftuples .and. phasetuple(leq)%nextcs.gt.0)
      leq=phasetuple(leq)%nextcs
   enddo
!   write(*,56)'3B setting nextcs in tuple: ',iph,phases(iph),nyttcs,leq,tuple
!56 format(a,10i5)
   phasetuple(leq)%nextcs=tuple
   nooftuples=tuple
!   write(*,*)'3B Adding phase tuple: ',tuple,iph,phases(iph)
! save index of tuple in new phase_varres record
   firsteq%phase_varres(nyttcs)%phtupx=tuple
!   write(*,31)'3B Phase tuple: ',nyttcs,tuple,iph,icsno,phases(iph)
31 format(a,10i5)
!   firsteq%phase_varres(lastcs)%phtupx=tuple
!   peq=>eqlista(1)%phase_varres(lastcs)
   peq=>firsteq%phase_varres(lastcs)
! sum up number of constituents!!
   nz=phlista(lokph)%tnooffr
!   write(*,*)'3B check: ',phlista(lokph)%nooffr,size(peq%yfr)
!   write(*,*)'3B added compset: ',iph,icsno,noeq()
!-------------------------------------------------------------------
! loop for all equilibria
   alleq: do leq=1,noeq()
! LOOP for all equilibria records to add this composition set to phase lokph
! lastcs is the previously last composition set, nyttcs is the new,
! same in all equilibria, also for firsteq (eqlista(1))!!
      neq=>eqlista(leq)%phase_varres(nyttcs)
!      write(*,19)'3B equil loop 1: ',leq,eqlista(leq)%eqno,lokph,icsno,&
!           phlista(lokph)%linktocs(icsno),nyttcs,tuple,neq%phlink
19    format(a,10i4)
! why is phlista updated here? It is outside the equilibrium record ...
!      phlista(lokph)%linktocs(icsno)=nyttcs
      neq%phlink=lokph
!      write(*,19)'3B equil loop 2: ',phlista(lokph)%linktocs(icsno),neq%phlink
! prefix and suffix, only letters and digits allowed but not checked ...
      pfix=prefix; sfix=suffix; call capson(pfix); call capson(sfix)
      neq%prefix=pfix
      neq%suffix=sfix
! tuple index
      neq%phtupx=tuple
! initiate the phstate as entered (value 0)
      neq%phstate=PHENTERED
! increment composition set counter when leq=1, phlista same in all equilibria
      if(leq.eq.1) then
         phlista(lokph)%linktocs(icsno)=nyttcs
         phlista(lokph)%noofcs=phlista(lokph)%noofcs+1
      endif
!      write(*,19)'3B add tupple: ',leq,nooftuples,tuple,neq%phtupx,icsno,&
!           nyttcs,phlista(lokph)%linktocs(icsno),&
!           firsteq%phase_varres(nyttcs)%phtupx
!      write(*,311)'3B sites: ',leq,iph,icsno,neq%sites
! sites, abnorm and amount formula units 
      if(.not.allocated(neq%sites)) then
!         write(*,*)'3B allocation 1: ',nsl
         allocate(neq%sites(nsl))
      endif
      neq%sites=peq%sites
      neq%abnorm=peq%abnorm
      neq%amfu=zero
!      write(*,311)'3B amfu: ',leq,iph,icsno,neq%amfu,neq%abnorm,peq%abnorm
311   format(a,3i3,6(1pe12.4))
! NOTE: these allocations below because create_parrecords does not work ...
! fractions and related
! NOTE: peq%yfr in firsteq is allocated maxconst=1000 as it is done
! before any elements entered!!! nz set above!!
!      nz=size(peq%yfr)
      if(.not.allocated(neq%yfr)) then
!         write(*,*)'3B allocation 2: ',nz
         allocate(neq%yfr(nz))
         neq%yfr=peq%yfr
      endif
! mmyfr is allocated here ...
!      write(*,*)'3B enter_compset: ',allocated(peq%mmyfr)
      if(allocated(peq%mmyfr)) then
         if(.not.allocated(neq%mmyfr)) then
!            write(*,*)'3B allocation 3: ',nz
            allocate(neq%mmyfr(nz))
            neq%mmyfr=peq%mmyfr
         endif
      endif
      if(allocated(peq%dpqdy)) then
! for ionic liquid, emergency bugfix 2017/02/16 Bo+Karl
         if(.not.allocated(neq%dpqdy)) then
            jz=size(peq%dpqdy)
            allocate(neq%dpqdy(jz))
            neq%dpqdy=peq%dpqdy
            jz=size(peq%d2pqdvay)
            allocate(neq%d2pqdvay(jz))
            neq%d2pqdvay=peq%d2pqdvay
         endif
      endif
! end bugfix
      if(.not.allocated(neq%constat)) then
! important!! constat has identification of the vacancy constituent !!
!         write(*,*)'3B allocation 4: ',nz
         allocate(neq%constat(nz))
         neq%constat=peq%constat
      endif
! copy status word but clear some bits CSDEFCON means default constitution
      neq%status2=peq%status2
      neq%status2=ibclr(neq%status2,CSDEFCON)
! set duplicate bit for auto in all equilibria
      if(len(suffix).ge.4) then
         if(suffix.eq.'AUTO') then
!            write(*,*)'3B setting bit CSTEMPAR in ',leq,nyttcs
            neq%status2=ibset(neq%status2,CSTEMPAR)
         endif
      endif
!
      if(.not.allocated(neq%gval)) then
! result arrays should have been allocated in create_parrecords ...
! but I do not call create_parrecords !!
!         write(*,83)'3B gval: ',leq,lokph,nyttcs,nprop,nz
83       format(a,10i5)
         allocate(neq%gval(6,nprop))
         allocate(neq%dgval(3,nz,nprop))
         allocate(neq%d2gval(nz*(nz+1)/2,nprop))
         allocate(neq%listprop(nprop))
      endif
!------------------- add addg ...
      if(btest(neq%status2,CSADDG)) then
         if(.not.allocated(neq%addg)) then
!            write(*,*)'3B allocation 6: ',1
            allocate(neq%addg(1))
            neq%addg(1)=peq%addg(1)
         endif
      endif
!--------------------
!      write(*,88)'3B cs: ',nz,neq%status2,neq%constat
88    format(a,i2,2x,Z16,2x,10(1x,i3))
! if there is a disordered fraction set one must copy the fraction set record
! and add a new parrecords to this. lokcs1 is first composition set
      disordered: if(btest(phlista(lokph)%status1,phmfs)) then
! copy the old fraction set record to the new
!------------------------ does this work??? disfra has a lot of data
         neq%disfra=peq%disfra
!------------------------- yes it works!!
!         write(*,*)'3B disfra 1: ',peq%disfra%ndd,neq%disfra%ndd
!         write(*,*)'disfra 2: ',peq%disfra%dxidyj(2),neq%disfra%dxidyj(2)
!--------------------------------------
         nsl=peq%disfra%ndd
         nkk=peq%disfra%tnoofxfr
!         write(*,*)'3B Creating disordered fraction set 1',lokcs1,nyttcs,nkk
         do jl=1,nkk
            iva(jl)=ceq%phase_varres(lokcs1)%constat(jl)
         enddo
         if(leq.eq.1) then
! allocate a parrecord for DISORDERED FRACTION SET for first equilibrium.
! Then use the same index: nydis, for all other equilibria.
! Maybe this can be made by a simple assignement????
            call create_parrecords(lokph,nydis,nsl,nkk,maxcalcprop,iva,firsteq)
            if(gx%bmperr.ne.0) goto 1000
         elseif(once) then
            write(kou,*)'3B creates a composition set in all equilibria'
            once=.FALSE.
!            write(kou,170)trim(eqlista(leq)%eqname),leq,lokcs1,nydis
!170         format('3B New composition set in equilibrium ',a,i4,&
!                 ' with lokcs and nydis index: ',2i4)
         endif
         ndeq=>eqlista(leq)%phase_varres(nydis)
         ndeq%phlink=lokph
         ndeq%prefix=' '
         ndeq%suffix=' '
! sites must be copied to disordered phase_varres
!         write(*,*)'3B dsites: ',size(neq%disfra%dsites),size(neq%sites)
         ndeq%disfra%dsites=peq%disfra%dsites
! some status bits must be set
         ndeq%status2=ibset(ndeq%status2,CSDFS)
         neq%status2=ibset(neq%status2,CSDLNK)
! set the link from ordered disfra record to the disordered phase_varres record
         neq%disfra%varreslink=nydis
      endif disordered
   enddo alleq
! end threadprotected code <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< dpqdy
!-------------------------------------------------
!   write(*,*)'3B Link from ordred ',lastcs,&
!        ' to disordered ',ceq%phase_varres(lastcs)%disfra%varreslink
!   next=ceq%phase_varres(lastcs)%next
!   write(*,*)'3B Link from ordred ',next,&
!        ' to disordered ',ceq%phase_varres(next)%disfra%varreslink
1000 continue
! test if mmy is correct in all existing compsets
! OK here also ...
!   do jl=1,icsno
!      lokcs=phlista(lokph)%linktocs(jl)
!      write(*,7)lokcs,firsteq%phase_varres(lokcs)%mmyfr
!   enddo
   return
 end subroutine enter_composition_set

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine suspend_composition_set(iph,parallel,ceq)
! the last composition set is suspended in all equilibria
!
! If parallel is TRUE then execution is not in parallel (threaded)
!
   implicit none
   integer iph
   logical parallel
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   TYPE(gtp_phase_varres), pointer :: varres,disvarres
   integer ics,lokph,lokcs,ncs,nsl,nkk,lastcs,nprop,idisvarres,kcs,leq
   lokph=phases(iph)
   ncs=phlista(lokph)%noofcs
! cannot remove composition set 1 or a nonexisting one
   if(ncs.le.1) goto 1000
   lokcs=phlista(lokph)%linktocs(ncs)
!   write(*,*)'3B suspend compset ',parallel
   if(parallel) then
! we have to stop all threads to do anyting with other equilibria, to
! suspend composition sets in other threads, skip that just suspend the
! last composition set of iph in this equilibrium, ceq
!$      if(omp_get_num_threads().eq.1) then
!$         write(*,*)'3B suspend ',iph,ncs
!$         if(btest(ceq%phase_varres(lokcs)%status2,CSTEMPAR)) then
!$            ceq%phase_varres(lokcs)%phstate=PHSUS
!$         endif
!-$      else
!-$        write(*,*)' *** Cannot suspend_composition_set in parallel'
!$      endif
      goto 1000
   endif
! we have many equilibria but is not running parallel
! suspend last composition set of iph in all equilibria where it is not stable
   do leq=1,noeq()
!      write(*,*)'3B suspend ',iph,ncs,&
!           eqlista(leq)%phase_varres(lokcs)%phstate,&
!           btest(eqlista(leq)%phase_varres(lokcs)%status2,CSAUTO),&
!           btest(eqlista(leq)%phase_varres(lokcs)%status2,CSTEMPAR)
      if(btest(eqlista(leq)%phase_varres(lokcs)%status2,CSTEMPAR) .and. &
           eqlista(leq)%phase_varres(lokcs)%phstate.le.PHENTERED) then
         eqlista(leq)%phase_varres(lokcs)%phstate=PHSUS
      endif
   enddo
!      
1000 continue
 end subroutine suspend_composition_set

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!
!\begin{verbatim} %-
 subroutine remove_composition_set(iph,force)
! the last composition set is deleted, update csfree and highcs
! SPURIOUS ERRORS OCCUR IN THIS SUBROUTINE
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>> NOTE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
! Not safe to remove composition sets when more than one equilibrium       !
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
!
! If force is TRUE delete anyway ... very dangerous ...
!
   implicit none
!
! BEWARE must be for all equilibria but maybe not allowed when threaded
!
   integer iph,jl,tuple
   logical force
!\end{verbatim}
   TYPE(gtp_phase_varres), pointer :: varres,disvarres
   integer ics,lokph,lokcs,ncs,nsl,nkk,lastcs,nprop,idisvarres,kcs,leq
!
   write(*,*)'In remove_compsets',iph
   if(iph.le.0 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   endif
   lokph=phases(iph)
   ncs=phlista(lokph)%noofcs
   if(ncs.eq.1) then
! cannot remove composition set 1 or a nonexisting one
      gx%bmperr=4093; goto 1000
   else
      ics=ncs
   endif
   if(btest(globaldata%status,GSNOREMCS)) then
      write(*,*)'3B Not allowed to delete composition sets'
      gx%bmperr=4211; goto 1000
   endif
!   write(*,*)'3B Delete highest composition set: ',iph,lokph,ics
   if(noeq().gt.1) then
! the deletion of composition sets when many equilibia not allowed until
! further testing
      write(*,*)' Warning, attempt to delete composition set',&
           ' with many equilibria ignored'
      goto 1000
      if(force) then
         write(*,*)' *** WARNING: deleting composition sets',&
              ' in many equilibria may cause errors'
      else
         write(*,*)'Attempt to delete composition sets when many equilibria'
         gx%bmperr=4211; goto 1000
      endif
   endif
!$   if(.TRUE.) then
!      write(*,*)'Deleting composition sets impossible when running parallel'
!       write(*,*)'This subroutine must be executed in sequential'
!$      goto 1000
!$   endif
! find the tuple for this phase+compset
   loop: do jl=1,nooftuples
!      write(*,*)'3B tuple compset: ',jl,ics,phasetuple(jl)%compset
!      if(phasetuple(jl)%phaseix.eq.lokph) then
      if(phasetuple(jl)%lokph.eq.lokph) then
         if(phasetuple(jl)%compset.eq.ics) then
            tuple=jl; exit loop
         endif
      endif
   enddo loop
!   write(*,*)'3B Delete composition set',iph,ics,lokph,tuple
   if(tuple.le.0) then
!      write(*,*)'No such tuple!!'
      gx%bmperr=4252; goto 1000
   endif
! collect some data
   nsl=phlista(lokph)%noofsubl
   nkk=phlista(lokph)%tnooffr
   lokcs=phlista(lokph)%linktocs(ics)
   lastcs=lokcs
   nprop=firsteq%phase_varres(lokcs)%nprop
!   write(*,*)'3B Removing varres record: ',lastcs
!-------------------------------------
! begin threadprotected code to remove lastcs >>>>>>>>>>>>>>>>>>>
! delete compset ics, shift higher down (not necessary)
! deallocate data in lokcs and return records to free list
!-------------------------------------
! note that the index to phase_varres is the same in all equilibria!!!!
   alleq: do leq=1,noeq()
      varres=>eqlista(leq)%phase_varres(lastcs)
      deallocate(varres%constat)
      deallocate(varres%yfr)
      deallocate(varres%mmyfr)
      deallocate(varres%sites)
! these may not be allocated ...
!      write(*,*)'3B delete varres dsitesdy: ',leq,lokcs,size(varres%dsitesdy)
!      if(size(varres%dsitesdy).gt.1) deallocate(varres%dsitesdy)
!      if(size(varres%d2sitesdy2).gt.1) deallocate(varres%d2sitesdy2)
      deallocate(varres%listprop)
      deallocate(varres%gval)
      deallocate(varres%dgval)
      deallocate(varres%d2gval)
! There is a disordered fraction record .... more to deallocate
      disordered: if(allocated(varres%disfra%y2x)) then
         deallocate(varres%disfra%dsites)
         deallocate(varres%disfra%nooffr)
         deallocate(varres%disfra%splink)
         deallocate(varres%disfra%y2x)
         deallocate(varres%disfra%dxidyj)
! now deallocate and release the phase_varres record with disordered fractions
         idisvarres=varres%disfra%varreslink
         disvarres=>eqlista(leq)%phase_varres(idisvarres)
!         write(*,*)'3B Deallocationg disordered varres record ',idisvarres
         deallocate(disvarres%constat)
         deallocate(disvarres%yfr)
         if(allocated(disvarres%mmyfr)) deallocate(disvarres%mmyfr)
         deallocate(disvarres%sites)
! these may not be allocated ...
!         write(*,*)'3B delete cs dsitesdy: ',leq,size(disvarres%dsitesdy)
!         if(size(disvarres%dsitesdy).gt.1) deallocate(disvarres%dsitesdy)
!         if(size(disvarres%d2sitesdy2).gt.1) deallocate(disvarres%d2sitesdy2)
         deallocate(disvarres%listprop)
         deallocate(disvarres%gval)
         deallocate(disvarres%dgval)
         deallocate(disvarres%d2gval)
! BOS 1401227: I do not think this is an error, just ignore ...
!         if(size(disvarres%disfra%dsites).gt.0) then
!            write(*,*)'ERROR, only one level of disordering allowed',leq,&
!                 size(disvarres%disfra%dsites)
!            stop
!         endif
      else
         idisvarres=0
      endif disordered
   enddo alleq
!   write(*,*)'3B Done all equilibrium records'
! decrement the composition set counter for this phase
! the phlista record is global, not part of the equilibria
   phlista(lokph)%noofcs=phlista(lokph)%noofcs-1
! link the released phase_varres record back to free list,
! maintained in firsteq only
   if(idisvarres.ne.0) then
! there was a disordered phase_varres record, link it into free list
!      write(*,*)'3B Free list 2: ',csfree,idisvarres
      firsteq%phase_varres(idisvarres)%nextfree=csfree
      csfree=idisvarres
! make used but released
      firsteq%phase_varres(idisvarres)%status2=&
           ibset(firsteq%phase_varres(idisvarres)%status2,CSDEL)
! UNFINISHED this is not correct ....
      idisvarres=newhighcs(.false.)
      if(idisvarres.eq.highcs) highcs=idisvarres-1
   endif
! link the free phase_varres into the free list
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! UNFINISHED: the free list for phase_varres is not updated correctly
! The use of csfree is DANGEROUS, there can be unallocated varres recored
! before the record indiceted by csfree
! and allocated after!!!
!   write(*,*)'3B Free list 1: ',csfree,lastcs
   firsteq%phase_varres(lastcs)%nextfree=csfree
   csfree=lastcs
! mark this record used but deleted
   firsteq%phase_varres(lastcs)%status2=&
        ibset(firsteq%phase_varres(lastcs)%status2,CSDEL)
! UNFINISHED this is not correct
   idisvarres=newhighcs(.false.)
   if(highcs.eq.lastcs) highcs=lastcs-1
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! finally shift all composition sets in phlista(lokph)%linktocs
! if last deleted then ics>phlista(lokph)%noofcs
   do kcs=ics,phlista(lokph)%noofcs
      phlista(lokph)%linktocs(kcs)=phlista(lokph)%linktocs(kcs+1)
   enddo
! and zero the last pointer to composition set.
   phlista(lokph)%linktocs(phlista(lokph)%noofcs+1)=0
!
! cleaning up phasetuple
   jl=phasetuple(tuple)%ixphase
!   write(*,*)
!   write(*,*)'3B cleaning up phase tuple when removing tuple: ',tuple,jl
   if(phasetuple(tuple)%compset.eq.2) then
! if the removed phasetuple has compset index 2 then zero the link in
! the original phase tuple ...
!      write(*,*)'3B link to tuple for compset 2 set to zero: ',tuple
      phasetuple(jl)%nextcs=0
   else
      jl=phasetuple(jl)%nextcs
! zero the nextcs pointer in the phase tuple pointing to tuple
      eternity: do while(phasetuple(jl)%nextcs.ne.tuple)
         if(jl.eq.phasetuple(tuple)%nextcs) then
            exit eternity
         endif
         if(phasetuple(jl)%nextcs.eq.0) then
!            write(*,*)'3B No such tuple: ',phasetuple(tuple)%compset,tuple
            gx%bmperr=4252; goto 1000
         endif
         jl=phasetuple(jl)%nextcs
      enddo eternity
      phasetuple(jl)%nextcs=0
   endif
!
!>>>>>>>>>>>>>>>>> THINK <<<<<<<<<<<<<<<<<<<<<<<
!
! The assumption is that phase tuples are always ordered in increasing
! composition set number.  One will always delete the highest number.
! The main problem is to ensure that %nextcs is correct and that the
! nextcs from the first composition set is updated correctly, also when
! phase tuples from other phases are deleted.
!
!   write(*,*)'3B Free list 1: ',csfree,highcs,lokcs
! update phasetuple array, overwrite tuple.  This means tuples may change phase
! NOTE the first tuple for a phase+compset will never change position.  Only
! those created later may be shifted ... but that may be complicated enough ...
!   write(*,*)'3B Shifting phase tuples above deleted: ',tuple,nooftuples
!   write(*,770)'3B1:',(jl,phasetuple(jl),jl=tuple-1,nooftuples)
770 format(a,3(6i4,';'),(/4x,6i4,';',6i4,';',6i4,';'))
! It is always the last compset of a phase that is removed,
! all nextcs links goes to higher tuples
   do jl=tuple+1,nooftuples
!      phasetuple(jl-1)%phaseix=phasetuple(jl)%phaseix
      phasetuple(jl-1)%lokph=phasetuple(jl)%lokph
      phasetuple(jl-1)%compset=phasetuple(jl)%compset
      phasetuple(jl-1)%ixphase=phasetuple(jl)%ixphase
      phasetuple(jl-1)%lokvares=phasetuple(jl)%lokvares
! all tuples have moved down one position ... thus nextcs decremented by one
      if(phasetuple(jl)%nextcs.gt.0) then
         phasetuple(jl-1)%nextcs=phasetuple(jl)%nextcs-1
      else
! unless it is zero in which case it keeps its value
         phasetuple(jl-1)%nextcs=0
      endif
! we must change the link to this tuple starting from ixphase ??
      if(phasetuple(jl-1)%compset.eq.2) then
!         write(*,*)'3B Changing link to compset 2: ',&
!              phasetuple(jl-1)%ixphase,jl-1
         phasetuple(phasetuple(jl-1)%ixphase)%nextcs=jl-1
      endif
!
! THERE IS SOME ERROR HERE ... macro Nestor-800 with 21 elements returned
! sometimes that a tuple did not exist.
!
! we must change the link in the phase_varres records also!!
!      lokph=phasetuple(jl-1)%phaseix
      lokph=phasetuple(jl-1)%lokph
      lokcs=phlista(lokph)%linktocs(phasetuple(jl-1)%compset)
      if(lokcs.le.0) then
         write(*,*)'3B index pf phase_varres <=0',jl-1,lokph
         gx%bmperr=4399; goto 1000
      endif
!      write(*,771)'3B Shifting down ',jl,nooftuples,phasetuple(jl-1)%phaseix,&
!           phasetuple(jl-1)%compset,lokph,lokcs
!771   format(a,10i5)
! in all equilibrium records, luckily the phase_varres record the same in all!!
      do leq=1,noeq()
         eqlista(leq)%phase_varres(lokcs)%phtupx=jl-1
      enddo
   enddo
!   write(*,770)'3B2:',(jl,phasetuple(jl),jl=tuple-1,nooftuples)
   nooftuples=nooftuples-1
! the last tuple must explicitly have its link set to zero ?? done
!   phasetuple(nooftuples)%nextcs=0
!      write(*,*)'3B Warning: phase tuples may have changed phase ...'
!   write(*,770)'3B 2: ',(phasetuple(jl),jl=tuple-4,nooftuples)
! end threadprotected code <<<<<<<<<<<<<<<<<<<<<<<<
!-------------------------
1000 continue
   return
 end subroutine remove_composition_set

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,&
      lfun,refx)
! enter a parameter for a phase from database or interactivly
! typty is the type of property, 1=G, 2=TC, ... , n*100+icon MQ&const#subl
! fractyp is fraction type, 1 is site fractions, 2 disordered fractions
! nsl is number of sublattices
! endm has one constituent index for each sublattice
! constituents in endm and lint should be ordered so endm has lowest
! (done by decode_constarr)
! nint is number of interacting constituents (can be zero)
! lint is array of sublattice+constituent indices for interactions
! ideg is degree
! lfun is link to function (integer index)
! refx is reference (text)
! if this is a phase with permutations all interactions should be in
! the first or the first two identical sublattices (except interstitals)
! a value in endm can be negative to indicate wildcard
! for ionic liquid constituents must be sorted specially
   implicit none
   integer, dimension(*) :: endm
   character refx*(*)
   integer lokph,fractyp,typty,nsl,nint,ideg,lfun
   integer, dimension(2,*) :: lint
!\end{verbatim}
   character notext*20,funexp*1024
   integer iord(maxsubl),jord(2,maxsubl)
   integer again,kkk,ll,kk1,mint,kk,lokint,iz,it,kint,ib,jl,zz
   integer lj,i1,i2,newint,ifri,lokcs,noperm,firstint,listfun
   integer, dimension(24) :: intperm
   integer, dimension(:,:), allocatable :: elinks
   integer, dimension(:,:), allocatable :: intlinks
   type(gtp_endmember), pointer :: newem,endmemrec,lastem
   type(gtp_interaction), pointer :: intrec,lastint,newintrec
!   type(gtp_interaction), allocatable, target :: newintrec
   type(gtp_property), pointer :: proprec,lastprop
   TYPE(gtp_fraction_set) :: disfra
   logical ionliq
!
   if(gx%bmperr.ne.0) then
      write(*,*)'3B Error ',gx%bmperr,' already set calling enter_parameter!'
      gx%bmperr=0
   endif
! listfun used when calling this routine just to list a parameter
   listfun=0
   if(fractyp.eq.2) goto 50
! this is for site fractions
!   write(*,6)'enter_parameter 1: ',lokph,nsl,phlista(lokph)%noofsubl,nint,ideg
6  format(a,10i5)
   if(nsl.ne.phlista(lokph)%noofsubl) then
      gx%bmperr=4065; goto 1000
   endif
   kkk=0
   jord=0
   sublloop: do ll=1,nsl
      emloop: do kk=1,phlista(lokph)%nooffr(ll)
         kk1=kkk+kk
!         write(*,12)lokph,nsl,ll,endm(ll),kk1,phlista(lokph)%constitlist(kk1)
!12       format('3B enter_parameter 2A: '4I4,5x,2i5)
         if(endm(ll).eq.phlista(lokph)%constitlist(kk1)) then
            iord(ll)=kk1
            goto 17
         endif
      enddo emloop
      if(endm(ll).eq.-99) then
! wildcard, sorted at the end
         iord(ll)=-99
      else
!         write(*,*)'3B error in enter_parameter ',endm(ll)
         gx%bmperr=4096; goto 1000
      endif
17     continue
      kkk=kkk+phlista(lokph)%nooffr(ll)
   enddo sublloop
!   write(*,*)'3B enter_parameter 2B: ',(iord(ll),ll=1,nsl)
!  if(nint.eq.2) write(*,*)'enter_parameter 2C: ************************ '
! end member constituents found, check interaction
! interactions are in sublattice order in lint
!80  continue
   mint=1
23 continue
   kkk=0
   if(mint.le.nint) then
      do ll=1,nsl
         if(lint(1,mint).eq.ll) then
            intloop: do kk=1,phlista(lokph)%nooffr(ll)
               kkk=kkk+1
!              write(*,15)mint,lint(2,mint),kkk,phlista(lokph)%constitlist(kkk)
               if(lint(2,mint).eq.phlista(lokph)%constitlist(kkk)) then
! write(*,*)'enter_parameter jord: ',mint,ll,kkk
!                  write(*,*)'3B Int no, subl, const: ',mint,ll,kkk
                  jord(1,mint)=ll
                  jord(2,mint)=kkk
                  mint=mint+1
!  write(*,*)'3B enter_parameter mint1: ',mint,ll,kkk,nint
                  if(mint.gt.nint) goto 28
                  goto 23
               endif
            enddo intloop
! a constituent does not exist in sublattice ll
!    write(*,16)ll,mint,lint(1,mint),lint(2,mint)
            gx%bmperr=4066; goto 1000
         endif
         kkk=kkk+phlista(lokph)%nooffr(ll)
      enddo
   endif
28  continue
!   write(*,*)'3B enter_parameter mint2: ',mint,nint
15  format('enter_parameter x: ',4I4)
16  format('enter_parameter y: ',4I4)
   if(mint.lt.nint) then
!      write(*,*)'3B enter_param error: ',nint,mint,lint(1,mint),lint(2,mint)
      gx%bmperr=4067; goto 1000
   endif
!   write(*,33)'3B epar 1: ',nint,((lint(i,j),i=1,2),j=1,nint)
33 format(a,i3,' : ',3(2i3,2x))
   goto 90
!----------------
! code below is for disordered fraction types, use fractset record
! one could try to handle both fraction types in the same code but
! that would just make it very very messy
50  continue
   if(.not.btest(phlista(lokph)%status1,PHMFS)) then
! there are no disordered fraction sets for this phase
      gx%bmperr=4068; goto 1000
   endif
   lokcs=phlista(lokph)%linktocs(1)
   disfra=firsteq%phase_varres(lokcs)%disfra
! number of sublattices in the disordered set
   if(nsl.ne.disfra%ndd) then
      gx%bmperr=4069; goto 1000
   endif
   kkk=0
!   write(*,*)'3B: disordered parameter: ',nsl
   do ll=1,nsl
      do kk=1,disfra%nooffr(ll)
         kk1=kkk+kk
!          write(*,12)ll,endm(ll),kk1,disfra%splink(kk1)
         if(endm(ll).eq.disfra%splink(kk1)) then
            iord(ll)=kk1
            goto 67
         endif
      enddo
      if(endm(ll).eq.-99) then
! wildcard
         iord(ll)=-99
      else
!         write(*,*)'3B in enter_parameter'
         gx%bmperr=4051; goto 1000
      endif
67     continue
      kkk=kkk+disfra%nooffr(ll)
   enddo
! check interaction constituents
   mint=1
73  continue
   kkk=0
   if(mint.le.nint) then
      do ll=1,nsl
         if(lint(1,mint).eq.ll) then
            do kk=1,disfra%nooffr(ll)
               kkk=kkk+1
               if(lint(2,mint).eq.disfra%splink(kkk)) then
                  jord(1,mint)=ll
                  jord(2,mint)=kkk
!   write(*,75)mint,lint(1,mint),lint(2,mint),kkk,ll,jord(1,mint),jord(2,mint)
75 format('ep 75: ',8i4)
                  mint=mint+1
                  if(mint.gt.nint) goto 78
                  goto 73
               endif
            enddo
! a constituent does not exist in sublattice ll
            gx%bmperr=4066; goto 1000
         endif
         kkk=kkk+disfra%nooffr(ll)
      enddo
   endif
78  continue
   if(mint.lt.nint) then
      gx%bmperr=4067; goto 1000
   endif
!---------------------------------------------------
! we have found all constituents for the end member and interactions
! now look if there are parameter records, otherwise create them
! try to keep end member records in some order of constituents ...
90 continue
!   if(fractyp.eq.2) then
!   write(*,92)'3B: endmembers: ',(iord(ii),ii=1,nsl)
!   write(*,92)'3B: interactions: ',(jord(2,ii),ii=1,nint)
!   endif
   nullify(lastem)
! check that interactions are in sublattice and alphabetical order!!
   again=0
   intcheck: do lokint=2,nint
      if(jord(1,lokint).lt.jord(1,lokint-1)) then
         corrsubl: do iz=1,2
            it=jord(iz,lokint)
            jord(iz,lokint)=jord(iz,lokint-1)
            jord(iz,lokint-1)=it
         enddo corrsubl
         again=1
      elseif(jord(1,lokint).eq.jord(1,lokint-1)) then
         if(jord(2,lokint).lt.jord(2,lokint-1)) then
            it=jord(2,lokint)
            jord(2,lokint)=jord(2,lokint-1)
            jord(2,lokint-1)=it
!            write(*,*)'3B interactions: ',jord(2,lokint),jord(2,lokint-1)
            again=1
         elseif(jord(2,lokint).eq.jord(2,lokint-1)) then
!            write(*,656)'3B Illegal with same interaction constituent twice',&
!                 phlista(lokph)%name
656         format(a/' phase: ',a)
            gx%bmperr=4266; goto 1000
         endif
      endif
   enddo intcheck
!   write(*,*)'3B Again: ',again
   if(again.eq.1) goto 90
! Make sure the endmember has the alphabetically lowest constituent
! and that the interaction is not the same as the endmember
!   write(*,92)'3B endmembers: ',(iord(i),i=1,nsl)
92 format(a,10i4)
!   write(*,92)'3B interactions: ',(jord(2,i),i=1,nint)
   placeibloop: do kint=1,nint
! ll is the sublattice with interaction
      ll=jord(1,kint)
      placeib: if(jord(2,kint).eq.iord(ll)) then
!         write(*,*)'pmod3B: Illegal with interaction with same constituent'
! subroutine enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,&
!      lfun,refx)
!         write(*,97)lokph,typty,fractyp,nsl,(endm(zz),zz=1,nsl),&
!              ideg,nint,(lint(1,zz),lint(2,zz),zz=1,nint)
97       format('pmod3B: Illegal with interaction with same constituent:'/&
              3i3,i4,2x,15(i5))
        gx%bmperr=4266; goto 1000
      elseif(jord(2,kint).lt.iord(ll)) then
! constituent in iord higher than that in jord, exchange jord and iord.  
         ib=iord(ll)
         iord(ll)=jord(2,kint)
         if(kint.eq.nint) then
! there are no more interactions, just put ib in the place of jord(2,kint)
            jord(2,kint)=ib
         else
! a bit problematic, we may have to shift constituents in jord
            moreint: do mint=kint+1,nint
               if(jord(1,mint).gt.ll) then
! next interaction in another sublattice, put ib in jord(2,mint-1)
                  jord(2,mint-1)=ib
               else
                  shiftint: if(ib.lt.jord(2,mint)) then
! next interaction is higher, put ib in jord(2,mint-1)
                     jord(2,mint-1)=ib
                  else
! interacting constituent is lower, we must shift constituents down in jord
! It can be done one at a time?? Example: user enter:
! L(fcc,D,E,C,A,B): iord(1)='D', jord(2,*)='A', 'B', 'C', 'E' (ordered above)
! kint=1 replaces iord(1)='A'; look for the place for 'D'; ninit=4
! loop mint=2 but 'D' is higher than 'B' so shift jord one step making
!    jord(2,*)='B', 'C', 'C', 'E'; 
! loop mint=3 but D is higher than 'C' so shift jord(2,*)='B', 'C', 'E', 'E'; 
! Now 'D' is lesser than 'E' so place it in jord(2,3):
! jord(2,*)='B', 'C', 'D', 'E'; 
                     jord(2,mint-1)=jord(2,mint)
                     if(mint.lt.nint .and. jord(1,mint+1).eq.ll) then
                        jord(2,mint)=jord(2,mint+1)
                     else
                        jord(2,mint)=ib
                     endif
                  endif shiftint
               endif
            enddo moreint
         endif
      endif placeib
   enddo placeibloop
! there may be permutations for ordered phases  ... implemented for fcc
   intperm=0
   ftyp1: if(fractyp.eq.1) then
      if(btest(phlista(lokph)%status1,PHFORD)) then
! These permutations may require 2 interaction records created ...
         call fccpermuts(lokph,nsl,iord,noperm,elinks,nint,jord,&
              intperm,intlinks)
         if(gx%bmperr.ne.0) goto 1000
! make sure iord is alphabtically ordered to find the correct parameter
! endm() are species index, iord() are constituent index (overal all subl)
! elinks are constituent index, iord(*,1) is identical to elinks(*,1)
         do jl=1,nsl
            iord(jl)=elinks(jl,1)
         enddo
      elseif(btest(phlista(lokph)%status1,PHBORD)) then
         call bccpermuts(lokph,nsl,iord,noperm,elinks,nint,jord,&
              intperm,intlinks)
         if(gx%bmperr.ne.0) goto 1000
! make sure iord is alphabtically ordered to find the correct parameter
! endm() are species index, iord() are constituent index (overal all subl)
         do jl=1,nsl
            iord(jl)=elinks(jl,1)
         enddo
      else
         noperm=1
      endif
   else
! fraction type 2 has no permutations
      noperm=1
   endif ftyp1
! parameters for site fractions
   if(fractyp.eq.1) then
      endmemrec=>phlista(lokph)%ordered
   else
      endmemrec=>phlista(lokph)%disordered
   endif
!   write(*,91)'3B enter_param 90: ',fractyp,nsl,(iord(ii),ii=1,nsl)
91 format(a,i2,i3,10i4)
   ionliq=btest(phlista(lokph)%status1,PHIONLIQ)
   findem: do while(associated(endmemrec))
      if(.NOT.ionliq) then
         lika:do lj=1,nsl
! iord(lj) can be negative for wildcard.  Wildcard endmedmemers at the end
            i1=iord(lj)
            i2=endmemrec%fraclinks(lj,1)
            if(i1.gt.0) then
               if(i2.lt.0 .or. i1.lt.i2) then
! The new end member record should be inserted before this record
                  goto 100
               elseif(i1.gt.i2) then
! continue searching for the end member or place to create it
                  lastem=>endmemrec
                  endmemrec=>endmemrec%nextem
                  cycle findem
               endif
! here i1<0
            elseif(i2.gt.0) then
! continue searching for the end member or place to create it
               lastem=>endmemrec
               endmemrec=>endmemrec%nextem
               cycle findem
            endif
! It is the same "wildcard" value if both i1 and i2 are negative
         enddo lika
      else
! for ionic liquids insert endmembers in order of second sublattice ...
! This is important as we want to calculate all parameters with anions
! before we come to vacancy and neutrals which should be multiplied with Q
         illika:do lj=nsl,1,-1
! iord(lj) can be negative for wildcard.  Wildcard endmedmemers at the end
            i1=iord(lj)
            i2=endmemrec%fraclinks(lj,1)
            if(i1.gt.0) then
               if(i2.lt.0 .or. i1.lt.i2) then
! The new end member record should be inserted before this record
                  goto 100
               elseif(i1.gt.i2) then
! continue searching for the end member or place to create it
                  lastem=>endmemrec
                  endmemrec=>endmemrec%nextem
                  cycle findem
               endif
! here i1<0
            elseif(i2.gt.0) then
! continue searching for the end member or place to create it
               lastem=>endmemrec
               endmemrec=>endmemrec%nextem
               cycle findem
            endif
! It is the same "wildcard" value if both i1 and i2 are negative
         enddo illika
      endif
!-------------------------------------------------
! found end member record with same constituents
      goto 200
   enddo findem
!
! if lfun=-1 we want to list the function and not create anything
   if(lfun.lt.0) goto 900
!
100 continue
! we have not found any endmember record so we have to insert a record here
! lokem may be nonzero if we exited from findem loop to this label
   call create_endmember(lokph,newem,noperm,nsl,iord,elinks)
!    write(*,*)'3B enter_par: created endmember ',new
   if(gx%bmperr.ne.0) goto 1000
! insert link to new from last end member record, lastem.
   if(.not.associated(lastem)) then
      if(fractyp.eq.1) then
         phlista(lokph)%ordered=>newem
      else
         phlista(lokph)%disordered=>newem
      endif
   else
!      emlista(lastem)%next=new
      lastem%nextem=>newem
   endif
! insert link from new to next (if lokem=0 this record is the last)
   newem%nextem=>endmemrec
   endmemrec=>newem
!---------------------------------------------------
! Here we have found or created the endmember record
! look for or create interaction record, no wildcards in interactions
! Interacting elements should be in sublattice and alphabetical order!!
200 continue
!   write(*,*)'3B enter_parameter mint3: ',mint,nint
   lokint=0
   someint: if(nint.gt.0) then
! when there are interaction records the ideal bit must be cleared
      phlista(lokph)%status1=ibclr(phlista(lokph)%status1,PHID)
! to locate interaction record,
      nullify(lastint)
      mint=1
      intrec=>endmemrec%intpointer
!      write(*,202)'3B enter_parameter 12A: ',lokph,typty,nsl,ideg,typty,lokem,&
!           (lint(1,i),i=1,nint),(lint(2,i),i=1,nint)
202   format(/a,7i4,4x,10i4)
      if(.not.associated(intrec)) then
! no interaction record for this endmember, create one unless lfun=-1
         if(lfun.eq.-1) goto 900
         call create_interaction(newintrec,mint,jord,intperm,intlinks)
         if(gx%bmperr.ne.0) goto 1000
!         write(*,*)'3B created interaction 9:',mint,nint
         endmemrec%intpointer=>newintrec
         intrec=>newintrec
         lastint=>intrec
!         mint=mint+1
         newint=1
!         write(*,*)'3B created interaction:  ',newint,mint
      else
!         write(*,*)'3B existing interaction: ',intrec%status
         newint=0
         firstint=0
      endif
300   continue
!      write(*,303)'3B  at 300A: ',lokph,newint,nint,mint,intrec%status
303   format(a,10i3)
! interaction records should be ordered according to the sublattice
! with the interaction.  For interaction with permutations use the 
! sublattice of the first permutation
      findint: do while(mint.le.nint)
!         write(*,*)'3B At findint: ',mint,nint,newint
         if(intrec%sublattice(1).eq.jord(1,mint) .and. &
              intrec%fraclink(1).eq.jord(2,mint)) then
! found an interaction with same constituent (maybe just created)
            if(mint.eq.nint) then
!               write(*,*)'3B same interaction, level: ',mint
               goto 400
            endif
            lastint=>intrec
            intrec=>intrec%highlink
            mint=mint+1
            newint=1
            if(.not.associated(intrec)) exit findint
         else
            if(mint.eq.nint) then
! error when storing permutations because newint=0 below.  Moved it to the end
! but that gave error L(liq,C,Cr,V) was stored as L(Liq,C,Cr,Fe,V)
! Add a check on mint, if mint=nint one cannot store it as higher
               newint=0
            endif
! we must store interactions in sublattice order and in order of constituent
! in jord(2,mint) otherwise we will never be able to find a permutation. 
            if(intrec%sublattice(1).gt.jord(1,mint)) then
!               write(*,*)'3B insering interaction before existing'
               exit findint
            endif
            lastint=>intrec
            intrec=>intrec%nextlink
            if(.not.associated(intrec)) exit findint
            firstint=1
! more records on this interaction level ?
! this worked for permutations but gave other errors, see above
!            newint=0
         endif
      enddo findint
! we can be here either because mint>nint or no more interaction records
! we must create at least one interactionrecord, newint=0 if same level
! If intrec is associated the nextint link should be set to this
310    continue
!      write(*,*)'3B At 310',mint,nint
      if(mint.le.nint) then
! if lfun=-1 and parameter does not exist just skip away
         if(lfun.eq.-1) goto 900
!         write(*,303)'3B  Linking at 310:',mint,nint,newint,firstint
         call create_interaction(newintrec,mint,jord,intperm,intlinks)
         if(gx%bmperr.ne.0) goto 1000
         if(newint.eq.1) then
!            write(*,*)'3B Linking as higher'
            lastint%highlink=>newintrec
         elseif(associated(intrec)) then
!            write(*,*)'3B Linking as previous'
            newintrec%nextlink=>intrec
!            write(*,*)'3B Ho ho said the sixth'
            if(associated(lastint)) then
               lastint%nextlink=>newintrec
            else
! this should be linked from the endmember or lower order interaction
!               write(*,*)'3B No previous interaction on this level'
               endmemrec%intpointer=>newintrec
            endif
!            write(*,*)'3B Ho ho said the sixth'
         else
!            write(*,*)'3B Linking as next'
            lastint%nextlink=>newintrec
         endif
! redundant as newint set to 1 below ...
!         newint=0
         intrec=>newintrec
         lastint=>intrec
         mint=mint+1
! there may be more interaction records .... but they must all be created
         newint=1
         goto 310
      endif
! Now we should have found or created the interaction record,
! check property list
400    continue
      proprec=>intrec%propointer
      if(.not.associated(proprec)) then
! do not create anything if lfun=-1
         if(lfun.eq.-1) goto 900
         call create_proprec(intrec%propointer,typty,ideg,lfun,refx)
      else
         goto 800
      endif
!     write(*,*)'3B enter_parameter 17: ',lokint,lokem,link
   else
! Found endmember and there is no interaction
! search the property list, there may not be the correct property!
      proprec=>endmemrec%propointer
      if(.not.associated(proprec)) then
! if on property record and lfun=-1 just list parameter equal to zero
         if(lfun.lt.0) goto 900
         call create_proprec(endmemrec%propointer,typty,ideg,lfun,refx)
      else
         goto 800
      endif
   endif someint
! all done
   goto 1000
!--------------------------------------------------------
! we found correct parameter record with a property, now search property list
! This loop both for endmembers and interactions
800 continue
   do while(associated(proprec))
      lastprop=>proprec
      if(proprec%proptype.eq.typty) then
! found property record, one should delete old and insert new function
! one must alse change the reference !!! And add the reference if new.
! mode=0 means no change of reference text if reference already exists
         call capson(refx)
         notext='*** Not set by user'
         call tdbrefs(refx,notext,0,ifri)
         if(ideg.le.proprec%degree) then
            if(lfun.eq.-1) then
               listfun=proprec%degreelink(ideg)
            else
               proprec%degreelink(ideg)=lfun
               proprec%reference=refx
            endif
         elseif(lfun.ge.0) then
            call extend_proprec(proprec,ideg,lfun)
            proprec%reference=refx
         endif
         if(lfun.eq.-1) goto 900
         goto 1000
      endif
      proprec=>proprec%nextpr
   enddo
! if lfun=-1 we just want to list a the parameter which is zero
   if(lfun.lt.0) goto 900
! no record for this property present, add a new property record
   call create_proprec(lastprop%nextpr,typty,ideg,lfun,refx)
! all done and go home
   goto 1000
!--------------------------------------------------------
! this is for listing parameter
900 continue
   write(*,*)'3B list parameter ',lfun,listfun
   if(listfun.gt.0) then
      call list_tpfun(listfun,0,funexp)
! for the moment use the TPFUN symbol ...
      call wrice2(kou,0,12,78,1,funexp)
   else
      write(kou,*)'Parameter is zero'
   endif
!----------------------------------------------------------
1000 continue
   if(gx%bmperr.eq.0) then
! mark that the phase has at least one parameter
      phlista(lokph)%status1=ibset(phlista(lokph)%status1,PHHASP)
   endif
!  write(*,*)'3B enter_parameter 99: ',gx%bmperr
!  write(*,1010)'enter_parameter 77: ',(phlista(lokph)%constitlist(i),i=1,6)
!1010 format(A,6I3)
   return
 end subroutine enter_parameter
! lfun

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine fccpermuts(lokph,nsl,iord,noperm,elinks,nint,jord,intperm,intlinks)
! finds all fcc/hcp permutations needed for this parameter
! The order of elements in the sublattices is irrelevant when one has F or B
! ordering as all permutations are stored in one place (with some exceptions)
! Thus the endmembers are ordered alphabetically in the sublattices and also
! the interaction parameters.  Max 2 levels of interactions allowed.
   implicit none
   integer, dimension(*) :: iord,intperm
   integer, dimension(2,*) :: jord
   integer lokph,nsl,noperm,nint
!\end{verbatim} %+
   integer l2,ll,ib,again,clink,lshift,mshift,a211
   integer odd,inz,ip,iqq1,iqq2,isp,jb,jp,jsp,l3,level1,level2
   integer level2perm,lj,loksp,lsp,niqq1,nl1,nl2,nll,np,nq,nz
   integer, dimension(4) :: elal,esame
   integer, dimension(:,:), allocatable :: elinks
   integer, dimension(:,:), allocatable :: intlinks
   logical notsame
   character carr*64
!   integer, dimension(3) :: esame
!
!-------------------------------------------------------------------
!
! This is a very long and messy subroutine and it calls others that are
! equally complicated.  It is important it is understandable and correct,
! all possible cases has not been tested.  Do not try to simplify it by making
! it more messy, this subroutine is not important for calculating speed
! but the structure it creates is important for speed.
! The corresponing routine for bcc permutations is even worse ...
!
!-------------------------------------------------------------------
!
!   if(nint.eq.2) then
!      write(*,501)'3B fccpermuts1: ',jord(1,1),jord(2,1),jord(1,2),jord(2,2)
!   endif
! I assume the ordering is in the first 4 sublattices, that could be changed
   if(nsl.lt.4) then
      write(*,*)'3B There must be at least 4 sublattices for fcc/hcp option'
      gx%bmperr=4267; goto 1000
   endif
   if(nint.gt.2) then
      write(*,*)'3B Maximum 2nd level interaction with option F'
      gx%bmperr=4268; goto 1000
   endif
! rearrange constituents in alphabetcal order in the sublattices,
! change interactions also!
!   write(*,11)'3B fp1: ',(iord(i),i=1,4),nint,((jord(j,k),j=1,2),k=1,nint)
11 format(a,4i4,' interactions: ',i2,4i4)
   do l2=1,4
      if(iord(l2).gt.0) then
         loksp=phlista(lokph)%constitlist(iord(l2))
         elal(l2)=splista(loksp)%alphaindex
      else
         elal(l2)=iord(l2)
      endif
   enddo
!   write(*,11)'3B fp2: ',(elal(i),i=1,4),nint,((jord(j,k),j=1,2),k=1,nint)
   again=1
   lagain: do while(again.ne.0)
! yet another messy sorting 
      again=0
      do l2=1,3
         do ll=l2+1,4
            equal: if(elal(ll).lt.elal(ll-1)) then
               again=1
               ib=elal(ll)
               elal(ll)=elal(ll-1)
               elal(ll-1)=ib
!               write(*,*)'3B call 1',ll-1,elal(ll-1)
               call findconst(lokph,ll-1,elal(ll-1),iord(ll-1))
               if(gx%bmperr.ne.0) goto 1000
!               write(*,*)'3B call 2',ll,elal(ll)
               call findconst(lokph,ll,elal(ll),iord(ll))
               if(gx%bmperr.ne.0) goto 1000
! if there are interacting constituents in ll or ll-1 shift them also
               do lj=1,nint
                  if(jord(1,lj).eq.ll) then
! write(*,21)'3B fpi1: ',lj,jord(1,lj),jord(2,lj)
21 format(a,i2,2i4)
                     jord(1,lj)=ll-1
                     loksp=phlista(lokph)%constitlist(jord(2,lj))
                     ib=splista(loksp)%alphaindex
!                     write(*,*)'3B call 3',ll-1,ib
                     call findconst(lokph,ll-1,ib,jord(2,lj))
                     if(gx%bmperr.ne.0) goto 1000
! write(*,21)'3B fpi2: ',lj,jord(1,lj),jord(2,lj)
                  elseif(jord(1,lj).eq.ll-1) then
! write(*,21)'3B fpi3: ',lj,jord(1,lj),jord(2,lj)
                     jord(1,lj)=ll
                     loksp=phlista(lokph)%constitlist(jord(2,lj))
                     ib=splista(loksp)%alphaindex
!                     write(*,*)'33B call 4',ll,ib
                     call findconst(lokph,ll,ib,jord(2,lj))
                     if(gx%bmperr.ne.0) goto 1000
! write(*,21)'3B fpi4: ',lj,jord(1,lj),jord(2,lj)
                  else
!                     write(*,23)'3B No interactions in sublattice: ',jord(1,lj)
23 format(a,2i3)
                  endif
               enddo
            endif equal
         enddo
      enddo
   enddo lagain
! elements are now ordered in alphabetical order over the sublattices
! find how many equal
!   if(nint.eq.2) then
!      write(*,501)'3B fccpermuts2A: ',jord(1,1),jord(2,1),jord(1,2),jord(2,2)
!   endif
   esame=0
   ib=1
   esame(ib)=1
   do ll=2,4
      if(elal(ll).eq.elal(ll-1)) then
         esame(ib)=esame(ib)+1
      else
         ib=ib+1
         esame(ib)=1
      endif
   enddo
   if(jord(1,1).ne.jord(1,2)) then
! we can have a case AX:AY:A:A and it should not be changed to AXY:A:A:A
      notsame=.true.
   else
      notsame=.false.
   endif
! we must rearrange interactions so they are in the first sublattice with
! the same endmember element for each level separately
! This is probably redundant as decode_constarr also sorts
   do l2=1,nint
      ib=elal(jord(1,l2))
      do ll=1,jord(1,l2)-1
         if(elal(ll).eq.ib) then
!            write(*,*)'3B Shifting interacting constituent to sublattice: ',ll
            nll=ll
            if(l2.eq.2 .and. notsame) then
! if interactions should not be in same sublattice but with the same element
! in the endmember, increment ll to interact in next sublattice.  It should
! be the same endmember constituent there!
               if(ll.eq.jord(1,1)) nll=ll+1
!               write(*,*)'3B nll: ',ll,nll
            endif
            jord(1,l2)=nll
            loksp=phlista(lokph)%constitlist(jord(2,l2))
            ib=splista(loksp)%alphaindex
!            write(*,*)'3B call 5',nll,ib
            call findconst(lokph,nll,ib,jord(2,l2))
            if(gx%bmperr.ne.0) goto 1000
         endif
      enddo
   enddo
!   if(nint.eq.2) then
!      write(*,501)'3B fccpermuts2B: ',jord(1,1),jord(2,1),jord(1,2),jord(2,2)
!   endif
!   write(*,11)'3B fp3: ',(elal(i),i=1,4),nint,((jord(j,k),j=1,2),k=1,nint)
!   write(*,11)'3B fp4: ',(iord(i),i=1,4)
! make sure that any interaction is connected to the first possible endmember
! for example A:A,B:B:B should be changed to A,B:A:B:B
! Also A,C:A,B:A:A should be A,B:A,C:A:A to have a unique record
   do l2=1,nint
      lj=jord(1,l2)
      do ll=1,lj-1
! ll must be less than 4 in this loop
         equalem: if(elal(ll).eq.elal(lj)) then
            if(l2.eq.1 .or. .not.notsame) then
               jord(1,l2)=ll
               loksp=phlista(lokph)%constitlist(jord(2,l2))
               ib=splista(loksp)%alphaindex
!               write(*,*)'3B call 6',ll,ib
               call findconst(lokph,ll,ib,jord(2,l2))
               if(gx%bmperr.ne.0) goto 1000
            else
! l2 must be 2 here, i.e. second order interaction
               loksp=phlista(lokph)%constitlist(jord(2,1))
               ib=splista(loksp)%alphaindex
               loksp=phlista(lokph)%constitlist(jord(2,2))
               jb=splista(loksp)%alphaindex
               if(jb.lt.ib) then
! change them so the lowest constituent comes first in sublattice order
!                  write(*,*)'3B call 7',ll,jb
                  call findconst(lokph,ll,jb,jord(2,1))
                  if(gx%bmperr.ne.0) goto 1000
!                  write(*,*)'3B call 8',lj,ib
                  call findconst(lokph,lj,ib,jord(2,2))
                  if(gx%bmperr.ne.0) goto 1000
!                  write(*,*)'3B exchange: ',ib,jb,jord(2,1),jord(2,2)
               else
! The interactions should not be in same sublattice, the next sublattice
! must have the same endmember constituent as jord(1,1), put it there
                  if(ll.eq.jord(1,1)) then
                     nll=ll+1
                  else
                     nll=ll
                  endif
                  jord(1,l2)=nll
                  loksp=phlista(lokph)%constitlist(jord(2,l2))
                  ib=splista(loksp)%alphaindex
!                  write(*,*)'3B call 9',nll,ib
                  call findconst(lokph,nll,ib,jord(2,l2))
                  if(gx%bmperr.ne.0) goto 1000
               endif
            endif
         endif equalem
      enddo
   enddo
!   if(nint.eq.2) then
!      write(*,501)'3B fccpermuts2C: ',jord(1,1),jord(2,1),jord(1,2),jord(2,2)
!   endif
!--------------------------------
! now we can calculate the number of endmember permutations
! Generate also all endmember links in elinks to be stored in endmember record
   lshift=phlista(lokph)%nooffr(1)
   if(esame(1).eq.4) then
! all 4 equal
      noperm=1
      allocate(elinks(nsl,noperm))
      do ll=1,nsl
         elinks(ll,1)=iord(ll)
      enddo
   elseif(esame(1).eq.3) then
! first 3 equal, one different: A:A:A:B; A:A:B:A; A:B:A:A; B:A:A:A
      noperm=4
      allocate(elinks(nsl,noperm))
      do np=1,noperm
         do ll=1,nsl
            elinks(ll,np)=iord(ll)
         enddo
         if(np.lt.4) then
! shift the single different element forward step by step
            ib=iord(4-np)+lshift
            iord(4-np)=iord(5-np)-lshift
            iord(5-np)=ib
         endif
      enddo
   elseif(esame(1).eq.2) then
      if(esame(2).eq.2) then
! the two first equal and also last two: A:A:B:B
! A:B:A:B; A:B:B:A; B:A:B:A; B:B;A:A; B:A:A:B
! I have no idea how to make this into a loop so I handle each separately
         noperm=6
         allocate(elinks(nsl,noperm))
         np=1
         do ll=1,nsl
            elinks(ll,np)=iord(ll)
         enddo
! shift sublattice 2 and 3: A:B:A:B
         ib=iord(2)+lshift
         iord(2)=iord(3)-lshift
         iord(3)=ib
         np=np+1
         do ll=1,nsl
            elinks(ll,np)=iord(ll)
         enddo
! shift sublattice 3 and 4: A:B:B:A
         ib=iord(3)+lshift
         iord(3)=iord(4)-lshift
         iord(4)=ib
         np=np+1
         do ll=1,nsl
            elinks(ll,np)=iord(ll)
         enddo
! shift sublattice 1 and 2: B:A:B:A
         ib=iord(1)+lshift
         iord(1)=iord(2)-lshift
         iord(2)=ib
         np=np+1
         do ll=1,nsl
            elinks(ll,np)=iord(ll)
         enddo
! shift sublattice 2 and 3: B:B:A:A
         ib=iord(2)+lshift
         iord(2)=iord(3)-lshift
         iord(3)=ib
         np=np+1
         do ll=1,nsl
            elinks(ll,np)=iord(ll)
         enddo
! shift sublattice 2 and 4 (double lenght): B:A:A:B
         ib=iord(2)+2*lshift
         iord(2)=iord(4)-2*lshift
         iord(4)=ib
         np=np+1
         do ll=1,nsl
            elinks(ll,np)=iord(ll)
         enddo
      else
! the first two equal and last 2 different: A:A:B:C
         a211=1
         noperm=12
         allocate(elinks(nsl,noperm))
         call fccpe211(1,elinks,nsl,lshift,iord)
      endif
   elseif(esame(2).eq.3) then
! first different and last 3 equal: A:B:B:B; B:A:B:B; B:B:A:B; B:B:B:A
      noperm=4
      allocate(elinks(nsl,noperm))
      do np=1,noperm
         do ll=1,nsl
            elinks(ll,np)=iord(ll)
         enddo
         if(np.lt.4) then
! shift the single different element backward step by step
            ib=iord(np)+lshift
            iord(np)=iord(np+1)-lshift
            iord(np+1)=ib
         endif
      enddo
   elseif(esame(2).eq.2) then
! two equal but first and last different
      a211=2
      noperm=12
      allocate(elinks(nsl,noperm))
      call fccpe211(2,elinks,nsl,lshift,iord)
   elseif(esame(3).eq.2) then
! first two different but last two equal
      a211=3
      noperm=12
      allocate(elinks(nsl,noperm))
      call fccpe211(3,elinks,nsl,lshift,iord)
   else
! all 4 different
      noperm=24
      allocate(elinks(nsl,noperm))
      call fccpe1111(elinks,nsl,lshift,iord)
   endif
! always skip debug output of endmembers for interaction parameters
   intperm(1)=0
!   if(nint.eq.2) then
!      write(*,501)'3B fccpermuts3: ',jord(1,1),jord(2,1),jord(1,2),jord(2,2)
!   endif
   if(nint.gt.0) goto 200
! comment next line to have debug output
!   goto 200
!--------------------
! debug output of endmembers after rearranging
   carr='fp6: '
   ib=6
   l3=1
   do ll=1,4
      if(elal(ll).gt.0) then
         l2=len_trim(splista(species(elal(ll)))%symbol)
         write(carr(ib:),16)splista(species(elal(ll)))%symbol(1:l2)
16       format(a)
         ib=ib+l2
      else
         carr(ib:)='*'
         ib=ib+1
      endif
17    continue
      if(l3.le.nint) then
         if(jord(1,l3).eq.ll) then
            loksp=phlista(lokph)%constitlist(jord(2,l3))
            l2=len_trim(splista(loksp)%symbol)
            write(carr(ib:),18)splista(loksp)%symbol(1:l2)
18          format(',',a)
            ib=ib+l2+1
            l3=l3+1
            goto 17
         endif
      endif
      if(ll.lt.4) carr(ib:ib)=':'
      ib=ib+1
   enddo
!   write(*,19)carr(1:ib)
!   write(*,19)'3B fp7: ',esame,noperm
19 format(a,4i3,i5)
! More debug output: all endmember permutations
   do np=1,noperm
! listing indices in constituent list (stored in endmember record)
!      write(*,31)np,(elinks(ll,np),ll=1,nsl)
31    format('3B elinks: ',i3,3x,10i4)
   enddo
   do np=1,noperm
! Easier to check listing of permutations using constituent names
      carr=' '
      ib=1
      do ll=1,nsl
         if(elinks(ll,np).gt.0) then
            loksp=phlista(lokph)%constitlist(elinks(ll,np))
            l2=len_trim(splista(loksp)%symbol)
            write(carr(ib:),32)splista(loksp)%symbol(1:l2)
32          format(a,':')
            ib=ib+l2+1
         else
            carr(ib:)='*:'
            ib=ib+2
         endif
      enddo
!      write(*,33)np,carr
33    format('3B emperm ',i3,': ',a)
   enddo
! debug output of endmembers end
!--------------------
200 continue
! done arranging component array and permutations of endmembers
   if(nint.eq.0) then
      goto 1000
   endif
!===============================================
! Now the 1st level interactions ... store in intlinks(1..2)
   allocate(intlinks(2,100))
! intperm(1)=number of interaction permutations on level 1 for each endmember
!   on level 1 each endmember perumtation has the same
! intperm(2)=total number of permutation links for level 1
! intperm(3..) used for 2nd level
   select case(noperm)
   case default ! error
!      write(*,*)'3B Unknown case for endmemeber permutations: ',noperm
      gx%bmperr=4269
!----------
   case(1) ! A:A:A:A
!      if(nint.eq.2) then
!         write(*,501)'3B fccpermuts4: ',jord(1,1),jord(2,1),jord(1,2),jord(2,2)
!      endif
      if(jord(1,1).ne.1) then
!         write(*,*)'3B Interaction must be in sublattice 1'
         gx%bmperr=4270; goto 1000
      endif
      intperm(1)=4
      intperm(2)=4
      clink=jord(2,1)
! set links to interaction with same element in all 4 sublattices
      do l2=1,4
         intlinks(1,l2)=l2
         intlinks(2,l2)=clink
         clink=clink+lshift
      enddo
      level1=1
!----------
   case(4) ! A:A:A:B and A:B:B:B
      if(esame(1).eq.3) then
         if(jord(1,1).eq.1) then
! the interaction must be AX:A:A:B
            call fccint31(jord,lshift,intperm,intlinks)
            level1=2
         else
! the interaction must be A:A:A:BX
            intperm(1)=1
            intperm(2)=4
            intlinks(1,1)=4
            intlinks(2,1)=jord(2,1)
            do ll=2,4
               intlinks(1,ll)=5-ll
               intlinks(2,ll)=intlinks(2,ll-1)-lshift
            enddo
            level1=3
         endif
      elseif(jord(1,1).eq.2) then
! the interaction must be A:BX:B:B
         call fccint31(jord,lshift,intperm,intlinks)
         level1=4
      else
! the interaction must be AX:B:B:B
         intperm(1)=1
         intperm(2)=4
         intlinks(1,1)=1
         intlinks(2,1)=jord(2,1)
         do ll=2,4
            intlinks(1,ll)=ll
            intlinks(2,ll)=intlinks(2,ll-1)+lshift
         enddo
         level1=5
      endif
!----------
   case(6) ! A:A:B:B
      call fccint22(jord,lshift,intperm,intlinks)
      level1=6
!----------
   case(12) ! A:A:B:C; A:B:B:C; A:B:C:C
      if(a211.eq.jord(1,1)) then
         call fccint211(a211,jord,lshift,intperm,intlinks)
         level1=7
      else
! interaction with one of the single constituents
! a single permutation follows the single different element in 4 sublattices
! starting from sublattice 1.  There are 12 enemember permutations
         intperm(1)=1
!         intperm(2)=12
         intperm(2)=noperm
         l2=jord(1,1)
         ib=phlista(lokph)%constitlist(elinks(l2,1))
         intlinks(1,1)=jord(1,1)
         intlinks(2,1)=jord(2,1)
         do ll=2,noperm
            do l3=1,4
               jb=phlista(lokph)%constitlist(elinks(l3,ll))
               if(jb.eq.ib) goto 410
            enddo
            write(*,*)'3B Cannot find endmember element for premutation ',ll,ib
            gx%bmperr=4271; goto 1000
410         continue
            intlinks(1,ll)=l3
            mshift=(intlinks(1,ll)-intlinks(1,ll-1))*lshift
            intlinks(2,ll)=intlinks(2,ll-1)+mshift
!            write(*,422)ll,l3,jord(1,1),mshift,intlinks(1,ll),intlinks(2,ll)
         enddo
         level1=8
      endif
!----------
   case(24) ! A:B:C:D
      write(*,77)
77    format(' *** CONGRATULATIONS, '/&
           '     You must be the first to enter a parameter like this!!!')
      intperm(1)=1
      intperm(2)=noperm
      l2=jord(1,1)
! species number in endmember of interacting sublattice
      ib=phlista(lokph)%constitlist(elinks(l2,1))
      intlinks(1,1)=l2
      intlinks(2,1)=jord(2,1)
      do ll=2,24
         do l3=1,4
            jb=phlista(lokph)%constitlist(elinks(l3,ll))
            if(jb.eq.ib) goto 420
!            write(*,419)'3B elinks,ib: ',ll,l3,ib,jb,elinks(l3,ll)
!419         format(a,2i4,2x,3i4)
         enddo
         write(*,*)'3B Cannot find endmember element for premutation ',ll,ib
         gx%bmperr=4271; goto 1000
420      continue
         intlinks(1,ll)=l3
         mshift=(intlinks(1,ll)-intlinks(1,ll-1))*lshift
         intlinks(2,ll)=intlinks(2,ll-1)+mshift
!         write(*,422)ll,l3,jord(1,1),mshift,intlinks(1,ll),intlinks(2,ll)
422      format('3B spec: ',3i3,2x,i10,2x,2i10)
      enddo
      level1=9
   end select
500 continue
   if(nint.eq.1) goto 900
!================================================================
! 2nd level interaction permutations
!   write(*,*)'3B First level interaction type: ',level1
!   write(*,502)'3B elinks and jord: ',elal,((jord(i,j),i=1,2),j=1,2)
501 format(a,2(2i4,2x))
502 format(a,4(i4),' : ',2(2i4,2x))
!
! The simplest 2nd level interaction is in the same sublattice as first
   if(jord(1,2).eq.jord(1,1)) then
! AXY:B:C:D where X and Y are two different constituents (not A) and B, C, D
! can be any constituents.  There are no new permutations, just add Y
!      write(*,*)'3B shortcut'
      intperm(3)=1
      intperm(4)=1
      nz=intperm(2)
      loksp=phlista(lokph)%constitlist(jord(2,2))
      isp=splista(loksp)%alphaindex
      do np=1,intperm(2)
         intlinks(1,nz+np)=intlinks(1,np)
         call findconst(lokph,intlinks(1,np),isp,intlinks(2,nz+np))
         if(gx%bmperr.ne.0) goto 1000
      enddo
! for debug output
      goto 900
   endif
!-----------------------------------------------------------
   select case(level1)
   case default !error
      write(*,*)'3B Unknown case for permutations on level 1: ',level1
      gx%bmperr=4272
!-----------------------------------------------------------
   case(1) ! AXY:A:A:A or AX:AX:A:A or AX:AY:A:A
      call fccip2A(lokph,jord,intperm,intlinks)
      if(gx%bmperr.ne.0) goto 1000
!-----------------------------------------------------------
   case(2) ! AXY:A:A:B or AX:AY:A:B or AX:A:A:BY
!      write(*,*)'3B case 2: ',jord(1,2),jord(2,2)
      if(jord(1,2).eq.4) then
! AX:A:A:BY, there should be 12 permutations, no new on second level
         intperm(3)=1
         intperm(4)=1
         intperm(5)=12
         nz=intperm(2)
         loksp=phlista(lokph)%constitlist(jord(2,2))
         isp=splista(loksp)%alphaindex
         do np=1,4
! sublattice for B the same for 3 permutations
            do nq=1,3
               nz=nz+1
               intlinks(1,nz)=5-np
               call findconst(lokph,5-np,isp,intlinks(2,nz))
               if(gx%bmperr.ne.0) goto 1000
            enddo
         enddo
      else
! AX:AY:A:B
         call fccip2B(1,lokph,lshift,jord,intperm,intlinks)
         if(gx%bmperr.ne.0) goto 1000
      endif
!-----------------------------------------------------------
   case(3) ! A:A:A:BXY
! never here as taken care by shortcut above ??
      if(jord(1,2).ne.jord(1,1)) then
!         write(*,*)'3B Thinking error, restructure!'
         gx%bmperr=4273; goto 1000
      endif
!-----------------------------------------------------------
   case(4) ! A:BXY:B:B or A:BX:BY:B; no AY:BX:B:B as that would be case 5
! A:BX:BY:B
      call fccip2B(2,lokph,lshift,jord,intperm,intlinks)
      if(gx%bmperr.ne.0) goto 1000
!-----------------------------------------------------------
   case(5) ! AX:BY:B:B
! This parameter has just 4 endmember permutations.  On this level 3 more
! AX:B:B:B  AX:BY:B:B AX:B:BY:B AX:B:B:BY
! B:AX:B:B  B:AX:BY:B B:AX:B:BY BY:AX:B:B etc
      intperm(3)=1
      intperm(4)=3
      intperm(5)=12
      loksp=phlista(lokph)%constitlist(jord(2,2))
      isp=splista(loksp)%alphaindex
      nz=intperm(2)
      do np=1,4
         nll=intlinks(1,np)
         do ip=1,3
            nz=nz+1
            nll=nll+1
            if(nll.gt.4) nll=1
            intlinks(1,nz)=nll
            call findconst(lokph,nll,isp,intlinks(2,nz))
            if(gx%bmperr.ne.0) goto 1000
         enddo
      enddo
!      endif
!-----------------------------------------------------------
! This is the important one as it includes the reciprocal excess parameter
   case(6) ! AX:A:B:B or A:A:BX:B, 6 endmem and 2 level 1 permutations = 12
! AX:A:B:B: AX:AX:B:B: 1; 0 totally 6 permutations
! AX:A:B:B: AX:AY:B:B and AY:AX:B:B; 2 additional permutations, totally 24
      loksp=phlista(lokph)%constitlist(jord(2,2))
      jsp=splista(loksp)%alphaindex
      if(abs(jord(1,2)-jord(1,1)).gt.1) then
! level 2 interaction with another endmember constituent than level 1
! AX:A:BY:B; 2 additional permutations, totally 24
! The endmember permutations will put element B in sublattices:
! 3,4; 2,4; 2,3; 1,3; 1,2; 1,4; If that changes this must be changed too ...
         intperm(3)=1
         intperm(4)=2
         intperm(5)=24
         nz=intperm(2)
         nl1=3
         nl2=4
         do ip=1,6
            nz=nz+1
            intlinks(1,nz)=nl1
            call findconst(lokph,nl1,jsp,intlinks(2,nz))
            if(gx%bmperr.ne.0) goto 1000
            nz=nz+1
            intlinks(1,nz)=nl2
            call findconst(lokph,nl2,jsp,intlinks(2,nz))
            if(gx%bmperr.ne.0) goto 1000
            nz=nz+1
            intlinks(1,nz)=nl1
            call findconst(lokph,nl1,jsp,intlinks(2,nz))
            if(gx%bmperr.ne.0) goto 1000
            nz=nz+1
            intlinks(1,nz)=nl2
            call findconst(lokph,nl2,jsp,intlinks(2,nz))
            if(gx%bmperr.ne.0) goto 1000
            select case(nl1)
            case default
!               write(*,*)'3B Error in fccpermut, case(lavel1=6), case(nl1)'
               gx%bmperr=4274; goto 1000
            case(1) ! change nl2 to 2 or 4, nl1 should be 1
               if(nl2.eq.2) nl2=4 
               if(nl2.eq.3) nl2=2
            case(2) ! change nl2 to 3
               if(nl2.eq.3) then
                  nl1=1
                  nl2=3
               else
                  nl2=3
               endif
            case(3) ! change nl1 to 2
               nl1=2
            end select
         enddo
      else
! interaction with same endmember element in 2 different sublattices
!         write(*,*)'3B smart?'
         loksp=phlista(lokph)%constitlist(jord(2,1))
         isp=splista(loksp)%alphaindex
         if(isp.eq.jsp) then
! AX:AX:B:B or A:A:BX:BX, there are 12 permutations of AX:A:B:B on level 1
! but there are only 6 second level interactions
! The endmember permutations will put element A in sublattices:
! 1,2; 1,3; 1,4; 2,4; 3,4; 2,3;  and element B in sublattices:
! 3,4; 2,4; 2,3; 1,3; 1,2; 1,4; 
            intperm(3)=2
            intperm(4)=1
            intperm(5)=0
            intperm(6)=6
            nz=intperm(2)
            if(jord(1,1).eq.1) then
               nll=2
            else
               nll=4
            endif
            odd=1
            do np=1,12
               odd=1-odd
               do jp=1,intperm(4+odd)
! this loop is done 1 or 0 times twice; nll=2,3,4; 4,4,3 // 4,4,3; 3,2,4
                  nz=nz+1
                  intlinks(1,nz)=nll
                  call findconst(lokph,nll,jsp,intlinks(2,nz))
                  if(gx%bmperr.ne.0) goto 1000
! nz=  13,14,15,16,17,18,19
! nll=  2, 3, 4, 4, 4, 3, - if jord(1,1)=1
! nll=  4, 4, 3, 3, 2, 4, - if jord(1,1)=2
                  select case(nz)
                  case default
                     write(*,*)'3B Error in fccpermut, case(lavel1=6), nz=',nz
                     gx%bmperr=4274; goto 1000
                  case(13) ! change nll to 3 if 2, else same
                     if(nll.eq.2) nll=3  ! 3 or same
                  case(14)
! the if ..,
!                  if(nll.eq.4) then
!                     nll=3
!                  else
!                     nll=4
!                  endif
! is same as nll=7-nll
                     nll=7-nll
                  case(15,18) ! no change!!
                     continue
                  case(16)
                     if(nll.eq.3) nll=2
                  case(17)
                     if(nll.eq.4) nll=3
                     if(nll.eq.2) nll=4
                  end select
               enddo
            enddo
! if the case and loops above works they are smart and easy to understand ???
         else
! AX:AY:B:B or A:A:BX:BY
! In this case we have the sume number of level2 permutations as level1
! Just add an interaction on the other sublattice with same endmember
! The endmember permutations will put element A in sublattices:
! 1,2; 1,3; 1,4; 2,4; 3,4; 2,3;  and element B in sublattices:
! 3,4; 2,4; 2,3; 1,3; 1,2; 1,4; 
! The first interaction will be with the first of the sublattices, the
! second in the second, just switch
            intperm(3)=1
            intperm(4)=1
            intperm(5)=intperm(2)
            nz=intperm(2)
            do np=1,6
! Here AX:AY:B:B and AY:AX:B:B
               nz=nz+1
               nll=intlinks(1,nz-11)
               nl2=intlinks(1,nz-12)
               intlinks(1,nz)=nll
               write(*,73)'3B loop 6B: ',np,nll,nl2,nz
73             format(a,10i4)
               call findconst(lokph,nll,jsp,intlinks(2,nz))
               if(gx%bmperr.ne.0) goto 1000
! set the second interaction in sublattice with level 1 interaction
               nz=nz+1
               intlinks(1,nz)=nl2
               call findconst(lokph,nl2,jsp,intlinks(2,nz))
               if(gx%bmperr.ne.0) goto 1000
            enddo
! if the case and loops above works they are smart and easy to understand ???
         endif
      endif
!-----------------------------------------------------------
! Maybe this can wait a little ...
   case(7) ! AX:A:B:C or A:BX:B:C or A:B:CX:C
!      write(*,*)'3B Not implemented yet 7'
      gx%bmperr=4275
!-----------------------------------------------------------
! Maybe this can wait a little ...
   case(8) ! A:A:BX;C or similar
!      write(*,*)'3B Not implemented yet 8'
      gx%bmperr=4275
!-----------------------------------------------------------
! Maybe this can wait a little ...
   case(9) ! AX:B:C:D or similar
!      write(*,*)'3B Not implemented yet 9'
      gx%bmperr=4275
   end select
!-----------------------------------------------------------
! done permutations of interactions
!   write(*,510)'3B 510: ',(intperm(j),j=1,7)
510 format(a,10i4)
!------- debug output of first level interaction permutations
900 continue
! to skip remove comment on next line
! goto 1000
   if(nint.eq.2) then
!      write(*,905)'3B Permutations of endmem and intlevel 1: ',noperm,&
!           intperm(1),intperm(2)
!      write(*,905)'3B Permutations of intlevel 2: ',intperm(3),&
!           (intperm(3+i),i=1,intperm(3))
905   format(a,i5,2x,10i4)
   endif
! these are the base pointers to first and second level permutations
   iqq1=0
   iqq2=intperm(2)+1
   inz=0
   emdmem: do np=1,noperm
! for each endmember permutation there are intperm(1) level 1 permutations
      intlev1: do niqq1=1,intperm(1)
         iqq1=iqq1+1
         if(nint.eq.2) then
            level2=1
            if(intperm(3).eq.1) then
! there is a fixed number of 2nd level permutations
               level2perm=intperm(4)
            else
! the number of 2nd level interaction varies with the first level, it can be 0
               level2perm=intperm(3+niqq1)
               if(level2perm.eq.0) cycle intlev1
            endif
         else
! no 2nd level interaction
            iqq2=0
         endif
910      continue
         carr=' '
         ib=1
         subl: do ll=1,nsl
! endmember constituent, can be wildcard
            loksp=elinks(ll,np)
            if(loksp.gt.0) then
               loksp=phlista(lokph)%constitlist(loksp)
               lsp=len_trim(splista(loksp)%symbol)
               carr(ib:)=splista(loksp)%symbol(1:lsp)
               ib=ib+lsp
            else
               carr(ib:ib)='*'
               ib=ib+1
            endif
920         continue
            if(intlinks(1,iqq1).eq.ll) then
! level 1 interaction constituent
! NOTE: For error checks output of intlinks is more important than the
! constituent name in carr as the link also indicates the sublattice!!!
!               if(nint.eq.2) &
!                    write(*,922)1,iqq1,intlinks(1,iqq1),intlinks(2,iqq1)
922            format('3B intlinks: ',2i5,2x,2i5,2x,3i5)
               loksp=phlista(lokph)%constitlist(intlinks(2,iqq1))
               lsp=len_trim(splista(loksp)%symbol)
               carr(ib:)=','//splista(loksp)%symbol(1:lsp)
               ib=ib+lsp+1
            endif
            if(iqq2.gt.0) then
               if(intlinks(1,iqq2).eq.ll) then
! level 2 interaction constituent
! NOTE: For error checks output of intlinks is more important than the
! constituent name in carr as the link also indicates the sublattice!!!
!                write(*,922)2,iqq2,intlinks(1,iqq2),intlinks(2,iqq2),jord(2,2)
                  loksp=phlista(lokph)%constitlist(intlinks(2,iqq2))
                  lsp=len_trim(splista(loksp)%symbol)
                  carr(ib:)=','//splista(loksp)%symbol(1:lsp)
                  ib=ib+lsp+1
               endif
            endif
            if(ll.lt.nsl) then
               carr(ib:)=': '
               ib=ib+2
            endif
         enddo subl
         inz=inz+1
!         write(*,925)inz,carr(1:len_trim(carr))
925      format('3B inter perm ',i3,': ',a)
         if(iqq2.gt.0) then
! there are level2perm number of 2nd order permutations
            level2=level2+1
            iqq2=iqq2+1
            if(level2.le.level2perm) goto 910
         endif
      enddo intlev1
   enddo emdmem
!------- debug output end
1000 continue
   return
 end subroutine fccpermuts

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine fccip2A(lokph,jord,intperm,intlinks)
! 2nd level interaction permutations for fcc
   implicit none
   integer, dimension(*) :: intperm
   integer, dimension(2,*) :: jord,intlinks
   integer lokph
!\end{verbatim} %+
   integer loksp,isp,jsp,ij,nll,ll,iqq,nz,ik
! AX:A:A:A, 2nd level can be AXY:A:A:A, AX:AX:A:A or AX:AY:A:A
   loksp=phlista(lokph)%constitlist(jord(2,2))
   isp=splista(loksp)%alphaindex
!   write(*,2)'3B fccip2A1: ',((jord(i,j),i=1,2),j=1,2)
!2  format(a,2(2i3,2x))
! 2nd level interaction in another sublattice, AX:AX:A:A or AX:AY:A:A
   loksp=phlista(lokph)%constitlist(jord(2,1))
   jsp=splista(loksp)%alphaindex
!      write(*,*)'3B fccip2A2: ',isp,jsp
   if(isp.eq.jsp) then
! 2nd level interacting constituent same as first level constituent:
! Level 1:  Level2:
! AX:A:A:A; AX:AX:A:A; AX:A:AX:A; AX:A:A:AX      3 permutations
! A:AX:A:A; A:AX:AX:A; A:AX:A:AX                 2 permutations
! A:A:AX:A; A:A:AX:AX                            1 permutations
! A:A:A:AX; none                                 0 permutations
!         write(*,*)'3B same interaction constituent in different sublattices'
      intperm(3)=4
      intperm(4)=3
      intperm(5)=2
      intperm(6)=1
      intperm(7)=0
      intperm(8)=24
      iqq=intperm(2)
      do ij=1,3
! loop only to 3 as there is no 2nd level permutation for ij=4
         nll=intlinks(1,ij)
         do ll=1,intperm(3+ij)
            iqq=iqq+1
            nll=nll+1
            intlinks(1,iqq)=nll
            if(nll.gt.4) then
!               write(*,*)'3B Error in 2nd level interaction of AX:AX:A:A'
               gx%bmperr=4276; goto 1000
            endif
            call findconst(lokph,intlinks(1,iqq),isp,intlinks(2,iqq))
            if(gx%bmperr.ne.0) goto 1000
!               write(*,76)'3B loop:',ij,nll,iqq,intlinks(1,iqq),intlinks(2,iqq)
76             format(a,3i3,2x,2i4)
         enddo
      enddo
! debug output
!         nc=0
!         nc1=0
!         nc2=intperm(2)
!         do lj=1,4
!            do ljj=1,intperm(3+lj)
!               nc=nc+1
!               nc1=nc1+1
!               nc2=nc2+1
!               write(*,77)nc,lj,ljj,&
!                    (intlinks(i,nc1),i=1,2),(intlinks(i,nc2),i=1,2)
77             format('3B AX:AX:A:A: ',i3,2x,2i3,2x,2(2i4,2x))
!            enddo
!         enddo
   else
! If 2nd level interacting element different
! Level 1:  Level2:
! AX:A:A:A; AX:AY:A:A; AX:A:AY:A; AX:A:A:AY      3 permutations
! A:AX:A:A; AY:AX:A:A; A:AX:AY:A; A:AX:A:AY      3 permutations
! A:A:AX:A; AY:A:AX:A; A:AY:AX:A; A:A:AX:AY      3 permutations
! A:A:A:AX; AY:A:A:AX; A:AY:A:AX; A:A:AY:AX      3 permutations
!     write(*,*)'3B different interaction constituent in different sublattices'
      intperm(3)=1
      intperm(4)=3
      intperm(5)=12
      nz=intperm(2)
      do ik=1,4
! Note that these permutations include AY:AX:A:A linked from AX:A:A:A
! A first level interaction AY:A:A:A is stored in another interaction record
! with no link to this 2nd level interaction.
         nll=intlinks(1,ik)
         do ll=1,3
            nll=nll+1
            if(nll.gt.4) nll=1
            nz=nz+1
            intlinks(1,nz)=nll
            call findconst(lokph,nll,isp,intlinks(2,nz))
            if(gx%bmperr.ne.0) goto 1000
!               write(*,88)nz,ik,ll,intlinks(1,nz),intlinks(2,nz)
88          format('3B loop: ',3i3,2x,2i5)
         enddo
      enddo
   endif
1000 continue
   return
 end subroutine fccip2A

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine fccip2B(lq,lokph,lshift,jord,intperm,intlinks)
! 2nd level interaction permutations for fcc
   implicit none
   integer lq,lokph,lshift
   integer, dimension(*) :: intperm
   integer, dimension(2,*) :: jord,intlinks
!\end{verbatim} %+
   integer loksp,isp,jsp,ny,nz,mp,isub2,nll,ip,np
! lq=1 means AX:AY:A:B or AX:AX:A:B
! lq=2 means A:BX:BY:B or A:BX:BX:B
! This parameter has 4 endmember permuts each with 3 permuts on level 1
! if X is same as Y only 2; 1; 0
   loksp=phlista(lokph)%constitlist(jord(2,1))
   isp=splista(loksp)%alphaindex
   loksp=phlista(lokph)%constitlist(jord(2,2))
   jsp=splista(loksp)%alphaindex
!   write(*,*)'3B fccip2B3: ',isp,jsp
   if(isp.eq.jsp) then
! Endmember  Level 1    Level 2   2; 1; 0;
! A:A:A:B    AX:A:A:B   AX:AX:A:B  AX:A:AX:B
!            A:AX:A:B   A:AX:AX:B
!            A:A:AX:B   none
! A:A:B:A    AX:A:B:A   AX:AX:B:A  AX:A:B:AX
!            A:AX:B:A   A:AX:B:AX
!            A:A:B:AX   none
! A:B:A:A    AX:B:A:A   AX:B:AX:A  AX:B:A:AX
!            A:B:AX:A   A:B:AX:AX
!            A:B:A:AX   none
! B:A:A:A    B:AX:A:A   B:AX:AX:A  B:AX:A:AX
!            B:A:AX:A   B:A:AX:AX
!            B:A:A:AX   none
! or the same for endmember A:B:B:B
      intperm(3)=3
      intperm(4)=2
      intperm(5)=1
      intperm(6)=0
      intperm(7)=intperm(2)
      ny=0
      nz=intperm(2)
      mp=3
! these loops are frustratingly messy .... but they seem to work ...
      nploop: do np=1,intperm(2)
         mp=mp+1
         if(lq.eq.1) then
! isub2 is the endmember sublattice occupied by the "different" constituent
!            isub2=(20-np)/4
            isub2=(15-np)/3
         else
!            isub2=(3+np)/4
            isub2=(2+np)/3
         endif
! nll is the sublattice with 1st level interaction
         ny=ny+1
         nll=intlinks(1,ny)
! np           = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
! mp           = 4, 5, 6, 4, 5, 6, 4, ...
! intperm(mp)  = 2, 1, 0, 2, 1, 0, 2, 1, 0,  2,  1,  0  
         do ip=1,intperm(mp)
            nll=nll+1
            if(nll.eq.isub2) nll=nll+1
            nz=nz+1
            intlinks(1,nz)=nll
!            write(*,13)'3B AX:AX:A:B: ',np,mp,ip,isub2,nz,nll,jsp
13          format(a,4i3,2x,i3,2i5)
            call findconst(lokph,nll,jsp,intlinks(2,nz))
            if(gx%bmperr.ne.0) goto 1000
         enddo
         if(mod(np,3).eq.0) mp=3
      enddo nploop
   else
! Endmember  Level 1    Level 2   2;
! A:A:A:B    AX:A:A:B   AX:AY:A:B  AX:A:AY:B
!            A:AX:A:B   A:AX:AY:B  AY:AX:A:B
!            A:A:AX:B   AY:A:AX:B  A:AY:AX:B
! A:A:B:A    AX:A:B:A   AX:AY:B:A  AX:A:B:AY etc
! There are 2 additional permutations for each of the 12 existing, the problem
! is mainly to know in which sublattice to add the interaction
      intperm(3)=1
      intperm(4)=2
      intperm(5)=2*intperm(2)
      ny=0
      nz=intperm(2)
      do np=1,intperm(2)
         if(lq.eq.1) then
! isub2 is the endmember sublattice occupied by the "different" constituent
            isub2=(15-np)/3
         else
! isub2 should be 1 for np=1..4, 2 for np=4..7 etc
            isub2=(np+2)/3
         endif
! nll is the sublattice with 1st level interaction
         ny=ny+1
         nll=intlinks(1,ny)
         do ip=1,2
! set 2nd interaction in sublattice after first interaction.  If that
! sublattice is >4 set it in first.  If the endmember is the single other
! constituent set it in next.  If that is >4 set it in first
            nll=nll+1
            if(nll.gt.4) nll=1
            if(nll.eq.isub2) nll=nll+1
            if(nll.gt.4) nll=1
            nz=nz+1
            intlinks(1,nz)=nll
!            write(*,13)'3B AX:AY:A:B: ',np,ip,0,isub2,nz,nll,jsp
            call findconst(lokph,nll,jsp,intlinks(2,nz))
            if(gx%bmperr.ne.0) goto 1000
         enddo
      enddo
   endif
1000 continue
   return
 end subroutine fccip2B

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine fccint31(jord,lshift,intperm,intlinks)
! 1st level interaction in sublattice l1 with endmember A:A:A:B or A:B:B:B
! set the sublattice and link to constituent for each endmember permutation
! 1st permutation of endmember: AX:A:A:B; A:AX:A:B; A:A:AX;B  4      0 1 2
! 2nd permutation of endmember: AX:A:B:A; A:AX:B:A; A:A:B:AX  3      0 1 3
! 3rd permutation of endmember: AX:B:A:A; A:B:AX:A; A:B:A:AX  3      0 2 3
! 4th permutation of endmember: B:AX:A:A; B:A:AX:A; B:A:A:AX  1 or   1 2 3
! 1st permutation of endmember: A:BX:B:B; A:B:BX:B; A:B:B:BX  4      0 1 2
! 2nd permutation of endmember: BX:A:B:B; B:A:BX:B; B:A:B:BX  1 etc -1 1 2
! 3rd -1 0 2 ; -1 0 1
! suck
   implicit none
   integer lshift
   integer, dimension(2,*) :: jord,intlinks
   integer, dimension(*) :: intperm
!\end{verbatim} %+
   integer l2,shift0,shift1,shift2,clink,idis,np
!
   intperm(1)=3
   intperm(2)=12
   l2=jord(1,1)
   clink=jord(2,1)
   idis=0
   shift0=0
   shift1=1
   shift2=2
   do np=1,4
      intlinks(1,idis+1)=l2+shift0
      intlinks(2,idis+1)=clink+shift0*lshift
      intlinks(1,idis+2)=l2+shift1
      intlinks(2,idis+2)=clink+shift1*lshift
      intlinks(1,idis+3)=l2+shift2
      intlinks(2,idis+3)=clink+shift2*lshift
      idis=idis+3
      subl: if(l2.eq.1) then
         if(np.eq.1) then
            shift2=3
         elseif(np.eq.2) then
            shift1=2
         elseif(np.eq.3) then
            shift0=1
         endif
      else
         if(np.eq.1) then
            shift0=-1
         elseif(np.eq.2) then
            shift1=0
         else
            shift2=1
         endif
      endif subl
   enddo
1000 return
 end subroutine fccint31

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine fccint22(jord,lshift,intperm,intlinks)
! 1st level for endmember A:A:B:B with interaction in sublattice jord(1,1) 
! 6 permutations of endmember, 2 permutations of interactions, 12 in total
! 1st endmemperm: AX:A:B:B; A:AX:B:B      0  1
! 2nd endmemperm: AX:B:A:B; A:B:AX:B      0  2
! 3rd endmemperm: AX:B:B:A; A:B:B:AX      0  3
! 4th endmemperm: B:AX:B:A; B:A:B:AX      1  3
! 5th endmemperm: B:B:AX:A; B:B:A:AX      2  3
! 6th endmemperm: B:AX:A:B; B:A:AX:B or   1  2
! 1th endmemperm: A:A:BX:B; A:A:B:BX      0  1
! 2nd endmemperm: A:BX:A:B; A:B:A:BX     -1  1
! 3rd endmemperm: A:BX:B:A; A:B:BX:A     -1  0
! 4th endmemperm: BX:A:B:A; B:A:BX:A     -2  0
! 5th endmemperm: BX:B:A:A; B:BX:A:A     -2 -1
! 6th endmemperm: BX:A:A:B; B:A:A:BX     -2  1
   implicit none
   integer lshift
   integer, dimension(2,*) :: jord,intlinks
   integer, dimension(*) :: intperm
!\end{verbatim} %+
   integer shift0,shift1,l2,clink,idis,np
!
   intperm(1)=2
   intperm(2)=12
   l2=jord(1,1)
   clink=jord(2,1)
   idis=0
   shift0=0
   shift1=1
   do np=1,6
      intlinks(1,idis+1)=l2+shift0
      intlinks(2,idis+1)=clink+shift0*lshift
      intlinks(1,idis+2)=l2+shift1
      intlinks(2,idis+2)=clink+shift1*lshift
      idis=idis+2
      subl: if(l2.eq.1) then
         select case(np)
         case default
            write(*,*)'3B Case error in fccint22: ',np
         case(1) !A:B:A:B is next endmember
            shift1=2
         case(2) !A:B:B:A
            shift1=3
         case(3) !B:A:B:A
            shift0=1
         case(4) !B:B:A:A
            shift0=2
         case(5) !B:A:A:B
            shift0=1
            shift1=2
         case(6) ! no more
         end select
      else
         select case(np)
         case default
            write(*,*)'3B Case error in fccint22: ',np
         case(1) !A:B:A:B is next endmember
            shift0=-1
         case(2) !A:B:B:A
            shift1=0
         case(3) !B:A:B:A
            shift0=-2
         case(4) !B:B:A:A
            shift1=-1
         case(5) !B:A:A:B
            shift1=1
         case(6) ! no more
         end select
      endif subl
   enddo
1000 continue
   return
 end subroutine fccint22

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine fccint211(a211,jord,lshift,intperm,intlinks)
! 1st level interaction in sublattice l1 with endmember like A:A:B:C
! 12 endmember permutations of AABC; ABBC; or ABCC
! 2 interaction permutations for each, 24 in total
   implicit none
   integer a211,lshift
   integer, dimension(2,*) :: jord,intlinks
   integer, dimension(*) :: intperm
!\end{verbatim} %+
   integer l2,clink,idis,shift0,shift1,np
   intperm(1)=2
   intperm(2)=24
   l2=jord(1,1)
   if(l2.ne.a211) then
!      write(*,*)'3B Error calling fccint211',a211,l2
      gx%bmperr=4276; goto 1000
   endif
   clink=jord(2,1)
   idis=0
   shift0=0
   shift1=1
! endmemeber A:A:B:C; first permutation interactions: AX:A:B:C; A:AX:B:C
! endmemeber A:B:B:C; first permutation interactions: A:BX:B:C; A;B:BX:C
! endmemeber A:B:C:C; first permutation interactions: A:B:CX:C; A:B:C:CX
   do np=1,12
      intlinks(1,idis+1)=l2+shift0
      intlinks(2,idis+1)=clink+shift0*lshift
      intlinks(1,idis+2)=l2+shift1
      intlinks(2,idis+2)=clink+shift1*lshift
      idis=idis+2
      subl: if(l2.eq.1) then
! endmember A:A:B:C
         select case(np)
         case default
            write(*,*)'3B Case error in fccint211: ',np,a211
         case(1) !A:A:C:B is next endmember 
            continue
         case(2) !A:C:A:B
            shift1=2
         case(3) !A:C:B:A
            shift1=3
         case(4) !A:B:C:A
            continue
         case(5) !A:B:A:C
            shift1=2
         case(6) !B:A:A:C
            shift0=1
         case(7) !B:A:C:A
            shift1=3
         case(8) !B:C:A:A
            shift0=2
         case(9) !C:B:A:A
            continue
         case(10) !C:A:B:A
            shift0=1
         case(11) !C:A:A:B
            shift1=2
         case(12) ! no more
         end select
      elseif(l2.eq.2) then
! endmember A:B:B:C
         select case(np)
         case default
            write(*,*)'3B Case error in fccint211: ',np,a211
         case(1) !A:B:C:B is next endmember
            shift1=2
         case(2) !C:B:A;B
            continue
         case(3) !C:B:B:A
            shift1=1
         case(4) !B:B:C:A
            shift0=-1
            shift1=0
         case(5) !B:B:A:C
            continue
         case(6) !B:A:B:C
            shift1=1
         case(7) !B:A:C:B
            shift1=2
         case(8) !C:A:B:B
            shift0=1
         case(9) !A:C:B:B
            continue
         case(10) !B:C:A:B
            shift0=-1
         case(11) !B:C:B:A
            shift1=1
         case(12) ! no more
         end select
      else
! endmember A:B:C:C
         select case(np)
         case default
            write(*,*)'3B Case error in fccint211: ',np,a211
         case(1) !A:C:B:C is next endmember
            shift0=-1
         case(2) !C:A:B:C
            shift1=0
         case(3) !C:B:A:C
            shift0=-2
         case(4) !B:C:A:C
            shift1=-1
         case(5) !B:A:C:C
            shift1=1
         case(6) !B:C:C:A
            shift1=1
         case(7) !C:B:C:A
            shift1=1
         case(8) !C:C:B:A
            shift1=1
         case(9) !C:C:A:B
            shift1=1
         case(10) !C:A:C:B
            shift1=1
         case(11) !A:C:C:B
            shift1=1
         case(12) ! no more
         end select
      endif subl
   enddo
1000 continue
   return
 end subroutine fccint211

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine fccpe211(l1,elinks,nsl,lshift,iord)
! sets appropriate links to constituents for the 12 perumations of
! A:A:B:C (l1=1), A:B:B:C (l1=2) and A:B:C:C (l1=3)
   implicit none
   integer, dimension(nsl,*) :: elinks
   integer, dimension(*) :: iord
   integer l1,nsl,lshift
!\end{verbatim} %+
   integer odd,np,ll,ib
! l1=1; keep 1 and change 3o4 and 2o3 6 times; then change 1o2 and
! loop 2 times
! changing 3o4 and 2o3; then change 1o2 and loop 2 times changing 2o3
! and 3o4
! AABC; AACB; ACAB; ACBA; ABCA; ABAC; ! BAAC; BACA; BCAA; ! CBAA;
! CABA; CAAB;
! l1=2; keep 2 and change 3o4 and 1o3 6 times; then change 2o3 and
! loop 2 times
! changing 3o4 and 1o3; then change 
! ABBC; ABCB; CBAB; CBBA; BBCA; BBAC; ! BABC; BACB; CABB; ! ACBB;
! BCAB; BCBA;
! l1=3; keep 4 and change 2o3 and 1o2 6 times; then change
! ABCC; ACBC; CABC; CBAC; BCAC; BACC; !  
!   write(*,*)'3B fccpe211: ',l1
   odd=0
   loop12: do np=0,11
      do ll=1,nsl
         if(iord(ll).lt.0) iord(ll)=-99
         elinks(ll,np+1)=iord(ll)
      enddo
! note l1 and ll are different !!!
      if(l1.eq.1) then
! AABC. Keep constituent in sublattice 1 first 6 loops; then for 3 and 3
         if(np.eq.5) then
            ib=iord(1)+lshift
            iord(1)=iord(2)-lshift
            iord(2)=ib
            odd=1-odd
         elseif(np.eq.8) then
            ib=iord(1)+lshift
            iord(1)=iord(2)-lshift
            iord(2)=ib
            odd=1-odd
         elseif(odd.eq.0) then
            ib=iord(3)+lshift
            iord(3)=iord(4)-lshift
            iord(4)=ib
            odd=1-odd
         else
            ib=iord(2)+lshift
            iord(2)=iord(3)-lshift
            iord(3)=ib
            odd=1-odd
         endif
      elseif(l1.eq.2) then
! ABBC. Keep constituent in sublattice 2 for first 6; then for 3 and 3
         if(np.eq.5) then
            ib=iord(2)+lshift
            iord(2)=iord(3)-lshift
            iord(3)=ib
            odd=1-odd
         elseif(np.eq.8) then
            ib=iord(1)+lshift
            iord(1)=iord(2)-lshift
            iord(2)=ib
            odd=1-odd
         elseif(odd.eq.0) then
            ib=iord(3)+lshift
            iord(3)=iord(4)-lshift
            iord(4)=ib
            odd=1-odd
         else
            ib=iord(1)+2*lshift
            iord(1)=iord(3)-2*lshift
            iord(3)=ib
            odd=1-odd
         endif
      else
! ABCC. Keep constituent in sublattice 4 for first 6; then for 3 and 3
         if(np.eq.5) then
            ib=iord(2)+2*lshift
            iord(2)=iord(4)-2*lshift
            iord(4)=ib
         elseif(np.eq.8) then
            ib=iord(3)+lshift
            iord(3)=iord(4)-lshift
            iord(4)=ib
            odd=1-odd
         elseif(odd.eq.0) then
            ib=iord(2)+lshift
            iord(2)=iord(3)-lshift
            iord(3)=ib
            odd=1-odd
         else
            ib=iord(1)+lshift
            iord(1)=iord(2)-lshift
            iord(2)=ib
            odd=1-odd
         endif
      endif
   enddo loop12
1000 continue
   return
 end subroutine fccpe211

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine fccpe1111(elinks,nsl,lshift,iord)
! sets appropriate links to 24 permutations when all 4 constituents different
! A:B:C:D
! The do loop keeps the same constituent in first sublattice 6 times, changing
! the other 3 sublattice, then changes the constituent in the first sublattice
! and goes on changing in the other 3 until all configurations done
   implicit none
   integer, dimension(nsl,*) :: elinks
   integer, dimension(*) :: iord
   integer nsl,lshift
!\end{verbatim}
   integer np,ll,odd,ib
! odd is either 0 or 1
   odd=1
   loop24: do np=0,23
      do ll=1,nsl
         if(iord(ll).lt.0) iord(ll)=-99
         elinks(ll,np+1)=iord(ll)
      enddo
! keep the same constituent in sublattice 1 for 6 endmembers, then shift
      if(np.eq.5) then
! shift 1 and 2, change odd
         ib=iord(2)-lshift
         iord(2)=iord(1)+lshift
         iord(1)=ib
         odd=1-odd
      elseif(np.eq.11) then
! shift 1 and 4, keep odd
         ib=iord(3)-2*lshift
         iord(3)=iord(1)+2*lshift
         iord(1)=ib
      elseif(np.eq.17) then
! shift 1 and 4, change odd
         ib=iord(4)-3*lshift
         iord(4)=iord(1)+3*lshift
         iord(1)=ib
         odd=1-odd
      elseif(odd.eq.0) then
         odd=1-odd
! shift 3 and 4
         ib=iord(4)-lshift
         iord(4)=iord(3)+lshift
         iord(3)=ib
      else
         odd=1-odd
! shift 2 and 3
         ib=iord(3)-lshift
         iord(3)=iord(2)+lshift
         iord(2)=ib
      endif
   enddo loop24
1000 continue
   return
 end subroutine fccpe1111

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine bccpermuts(lokph,nsl,iord,noperm,elinks,nint,jord,intperm,intlinks)
! finds all bcc permutations needed for this parameter
   implicit none
   integer lokph,nsl,noperm,nint
   integer, dimension(*) :: iord,intperm
   integer, dimension(2,*) :: jord
   integer, dimension(:,:), allocatable :: elinks
   integer, dimension(:,:), allocatable :: intlinks
!\end{verbatim}
! I assume the ordering is in the first 4 sublattices, that could be changed
   if(nsl.lt.4) then
      write(*,*)'3B There must be at least 4 sublattices for bcc option'
      gx%bmperr=4267; goto 1000
   endif
! unifinished
   write(*,*)'3B BCC permutations not implemented yet'
   gx%bmperr=4277
1000 continue
   return
 end subroutine bccpermuts

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine findconst(lokph,ll,spix,constix)
! locates the constituent index of species with index spix in sublattice ll
! and returns it in constix.  For wildcards spix is -99; return -99
! THERE MAY ALREADY BE A SIMULAR SUBROUTINE ... CHECK
   implicit none
   integer lokph,ll,spix,constix
!\end{verbatim}
   integer nc,l2,loksp
   if(spix.eq.-99) then
      constix=-99
      goto 1000
   endif
   nc=1
   do l2=1,ll-1
! The number of constituents in each sublattice can vary, add together
      nc=nc+phlista(lokph)%nooffr(l2)
   enddo
   constix=0
   do l2=nc,nc+phlista(lokph)%nooffr(ll)-1
      loksp=phlista(lokph)%constitlist(l2)
      if(splista(loksp)%alphaindex.eq.spix) then
         constix=l2; exit
      endif
   enddo
   if(constix.eq.0) then
!      write(*,90)spix,nc
90    format('3B No such constituent with index ',i5,' in sublattice',i3)
      gx%bmperr=4066; goto 1000
   endif
1000 continue
   return
 end subroutine findconst
 
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine tdbrefs(refid,line,mode,iref)
! store a reference from a TDB file or given interactivly
! If refid already exist and mode=1 then amend the reference text
   implicit none
   character*(*) refid,line
   integer mode,iref
!\end{verbatim}
   integer ip,ml,nr,mc,nc,jl
! make sure refid is left adjusted
   ip=0
10 ip=ip+1
   if(ip.gt.len(refid)) then
      gx%bmperr=4154; goto 1000
   endif
   if(refid(ip:ip).eq.' ') goto 10
   if(ip.gt.1) refid=refid(ip:)
! make it upper case
   call capson(refid)
! look if refid already exist
   do iref=1,reffree-1
      if(refid.eq.bibrefs(iref)%reference) then
         if(mode.eq.1) then
!            write(*,70)i,refid,bibrefs(i)%refspec
!70          format('3B tdbrefs: ',i4,a,a)
!            deallocate(bibrefs(iref)%refspec)
!            deallocate(bibrefs(iref)%nyrefspec)
            deallocate(bibrefs(iref)%wprefspec)
            goto 200
         else
! reference already exist and no changes needed
            goto 1000
         endif
      endif
   enddo
! if bibliographic reference does not exist do not create
   if(mode.eq.1) goto 1000
   iref=reffree
   reffree=reffree+1
   bibrefs(iref)%reference=refid
200 continue
   ml=len_trim(line)
!   nr=(ml+63)/64
!   allocate(bibrefs(iref)%refspec(nr))
   if(ml.gt.1024) then
      write(*,*)'Bibliographic references longer than 1024 will be truncated'
      mc=nwch(1024)+1
   else
      mc=nwch(ml)+1
   endif
   allocate(bibrefs(iref)%wprefspec(mc))
! This requires Fortran 2003/2008 standard
!   allocate(character(len=mc) :: bibrefs(iref)%nyrefspec)
!   mc=1
!   nc=0
!   bibrefs(iref)%nyrefspec=line(1:mc)
   bibrefs(iref)%wprefspec(1)=ml
   call storc(2,bibrefs(iref)%wprefspec,line(1:ml))
!   write(*,202)'3B newref: ',iref,refid,nr,line(1:min(32,len_trim(line)))
!202 format(a,i4,1x,a,i3,1x,a)
!   do jl=1,nr
! 1-64       mc=1, nc=64
! 65-122
!      bibrefs(iref)%refspec(jl)=' '
!      nc=nc+min(ml,64)
!      bibrefs(iref)%refspec(jl)=line(mc:nc)
!      mc=nc+1
!      ml=ml-64
!   enddo
1000 continue
   return
 end subroutine tdbrefs

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_equilibrium(name,number)
! creates a new equilibrium.  Allocates arrayes for conditions
! components, phase data and results etc.
! returns index to new equilibrium record
! THIS CAN PROBABLY BE SIMPLIFIED, especially phase_varres array can be
! copied as a whole, not each record structure separately ... ???
   implicit none
   character name*(*)
   integer number
!\end{verbatim}
! allocate
   TYPE(gtp_phase_varres), pointer :: cpv,cp1
   character name2*64
   integer ieq,ipv,nc,jz,iz,jl,jk,novarres
   if(.not.allowenter(3)) then
      write(*,*)'3B: not allowed enter equilibrium: ',name
      gx%bmperr=4153; goto 1000
   endif
   name2=name
   call capson(name2)
   if(ocv()) write(*,*)'3B In routine to enter equilibria: ',&
        name,noofph,eqfree
   if(.not.proper_symbol_name(name2,0)) then
! the name must start with a letter A-Z and contain letters, numbers and _
      gx%bmperr=4122
      goto 1000
   endif
   call findeq(name2,ieq)
   if(gx%bmperr.eq.0) then
! error as equilibrium with this name already exists
      gx%bmperr=4123
      goto 1000
   else
! OK to reset error code as we are creating a new equilibrium
      gx%bmperr=0
   endif
   if(eqfree.le.maxeq) then
      ieq=eqfree
      eqfree=eqfree+1
   endif
   number=ieq
   if(ocv()) write(*,*)'3B create eq',eqfree,maxeq,ieq
! allocate data arrayes in equilibrium record
   eqlista(ieq)%next=0
   eqlista(ieq)%eqname=name2
   eqlista(ieq)%eqno=ieq
   eqlista(ieq)%weight=-one
   eqlista(ieq)%comment=' '
! component list and matrix, if second or higher equilibrium copy content
   if(ocv()) write(*,*)'3B: entereq 1: ',maxel,ieq,noofel
   if(ieq.eq.1) then
      allocate(eqlista(ieq)%complist(maxel))
      allocate(eqlista(ieq)%compstoi(maxel,maxel))
      allocate(eqlista(ieq)%invcompstoi(maxel,maxel))
      allocate(eqlista(ieq)%cmuval(maxel))
      eqlista(ieq)%cmuval=zero
! this is a bit meaningless but skipping it has given raise to strange errors
      eqlista(ieq)%compstoi=zero
      eqlista(ieq)%invcompstoi=zero
      do jl=1,maxel
         eqlista(ieq)%compstoi(jl,jl)=one
         eqlista(ieq)%invcompstoi(jl,jl)=one
! valgrind complained this was not set !!
         eqlista(ieq)%complist(jl)%chempot=zero
      enddo
! Maybe valgrind complained of this ... it can have to do with -finit-local-zero
      eqlista(ieq)%status=0
   else
      eqlista(ieq)%status=0
      allocate(eqlista(ieq)%complist(noofel))
! copy mass of components, maybe other components?
      do jl=1,noofel
         eqlista(ieq)%complist(jl)%mass=firsteq%complist(jl)%mass
      enddo
      allocate(eqlista(ieq)%compstoi(noofel,noofel))
      allocate(eqlista(ieq)%invcompstoi(noofel,noofel))
      allocate(eqlista(ieq)%cmuval(noofel))
! this is a bit meaningless but skipping it has given raise to strange errors
      eqlista(ieq)%compstoi=zero
      eqlista(ieq)%invcompstoi=zero
      do jl=1,noofel
         eqlista(ieq)%compstoi(jl,jl)=one
         eqlista(ieq)%invcompstoi(jl,jl)=one
      enddo
      eqlista(ieq)%cmuval=zero
      if(ocv()) write(*,*)'3B: entereq 1B: '
      do jl=1,noofel
         eqlista(ieq)%complist(jl)%splink=firsteq%complist(jl)%splink
         eqlista(ieq)%complist(jl)%phlink=firsteq%complist(jl)%phlink
         eqlista(ieq)%complist(jl)%status=firsteq%complist(jl)%status
!         if(firsteq%complist(jl)%phlink.gt.0) then
! only if there is a defined reference state
         eqlista(ieq)%complist(jl)%refstate=firsteq%complist(jl)%refstate
         eqlista(ieq)%complist(jl)%tpref=firsteq%complist(jl)%tpref
         eqlista(ieq)%complist(jl)%chempot=zero
         do jk=1,noofel
            eqlista(ieq)%compstoi(jl,jk)=firsteq%compstoi(jl,jk)
            eqlista(ieq)%invcompstoi(jl,jk)=firsteq%invcompstoi(jl,jk)
         enddo
         if(allocated(firsteq%complist(jl)%endmember)) then
            iz=size(firsteq%complist(jl)%endmember)
            if(ocv()) write(*,*)'3B: entereq 1E: ',iz
            allocate(eqlista(ieq)%complist(jl)%endmember(iz))
            eqlista(ieq)%complist(jl)%endmember=&
                 firsteq%complist(jl)%endmember
         endif
!         endif
      enddo
   endif
!   write(*,*)'3B enter_eq 2, after this segmentation fault'
! these records keep calculated values of G and derivatives for each phase
! For phase lokph the index to phase_varres is in phlista(lokph)%cslink
! For phase lokph the index to phase_varres is in phlista(lokph)%linktocs(ics)
   if(ocv()) write(*,*)'3B: entereq 2: ',maxph
   alleq: if(ieq.eq.1) then
! %multiuse is used for axis and direction of a start equilibrium
      allocate(eqlista(ieq)%phase_varres(2*maxph))
      firsteq=>eqlista(ieq)
      firsteq%multiuse=0
! endif is at label 900, no need for goto
!      goto 900
   else
      eqlista(ieq)%multiuse=0
! UNFINISHED vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! this is not good, csfree is not the last used phase_varres
! there may be allocated records after and unallocated before !!
      if(highcs.ne.csfree-1) then
         write(*,*)'3B Beware, problems with varres records!',csfree,highcs
      endif
      novarres=highcs
! the next line should be removed when highcs correctly implemented
      novarres=csfree-1
      iz=noofph
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! for ieq>1 allocate the current number of phase_varres records plus 10
! for extra composition sets added later
!      allocate(eqlista(ieq)%phase_varres(iz+10))
      allocate(eqlista(ieq)%phase_varres(2*maxph))
!      write(*,*)'3B enter_eq 2B, after this segmentation fault'
!      write(*,*)'3B varres: ',ieq,size(eqlista(ieq)%phase_varres),iz
      if(ocv()) write(*,*)'3B varres: ',ieq,size(eqlista(ieq)%phase_varres)
! now copy the current content of firsteq%phase_varres to this equilibrium
! note, the SELECT_ELEMENT_REFERENCE phase has phase number 0
! and phase_varres index 1, the number of phase_varres records is not the
! same as number of phases ....
      copypv: do ipv=1,novarres
! note eqlista(1) is identical to firsteq
         if(.not.allocated(firsteq%phase_varres(ipv)%yfr)) then
! UNFINISHED this handels unallocated records below novarres
!            write(*,*)'3B error creating varres record',ipv
! BUT what about allocated after !!! no problem so far but .............
            cycle copypv
         endif
         cp1=>eqlista(1)%phase_varres(ipv)
         cpv=>eqlista(ieq)%phase_varres(ipv)
         cpv%nextfree=cp1%nextfree
         cpv%phlink=cp1%phlink
         cpv%status2=cp1%status2
         cpv%abnorm=cp1%abnorm
         cpv%prefix=cp1%prefix
         cpv%suffix=cp1%suffix
         cpv%phtupx=cp1%phtupx
! allocate and copy arrays
         nc=size(cp1%yfr)
! note SIZE gives rubbish unless array is allocated
         if(ocv()) write(*,*)'3B copy yfr 1: ',nc
         allocate(cpv%yfr(nc))
         cpv%yfr=cp1%yfr
! problems with phase_varres in equilibrium 2 ...
!         write(*,46)'3B 1: ',cp1%yfr
!         write(*,46)'3B v: ',cpv%yfr
46       format('yfr ',a,10(F7.3))
         allocate(cpv%constat(nc))
         cpv%constat=cp1%constat
!         write(*,*)'3B enter_eq 2C, after this segmentation fault'
         if(allocated(cp1%mmyfr)) then
! problem with mmyfr???  .... no
!            if(ocv()) write(*,*)'3B mmyfr 1: ',ipv,cpv%phlink,nc
            allocate(cpv%mmyfr(nc))
            cpv%mmyfr=cp1%mmyfr
!            write(*,34)'3B mmyfr 2: ',(cpv%mmyfr(jz),jz=1,nc)
34          format(1x,a,10(F7.3))
!         else
!            write(*,*)'3B mmyfr not allocated'
         endif
         jz=size(cp1%sites)
         allocate(cpv%sites(jz))
         cpv%sites=cp1%sites
! these are currently not allocated (ionic liquid model) Maybe not needed??
         if(allocated(cp1%dpqdy)) then
            jz=size(cp1%dpqdy)
            allocate(cpv%dpqdy(jz))
            cpv%dpqdy=cp1%dpqdy
            jz=size(cp1%d2pqdvay)
            allocate(cpv%d2pqdvay(jz))
            cpv%d2pqdvay=cp1%d2pqdvay
         endif
! the values in the following arrays are irrelevant, just allocate and zero
!         write(*,*)'3B enter_eq 2D, after this segmentation fault',ipv,novarres
         cpv%nprop=cp1%nprop
         allocate(cpv%listprop(cp1%nprop))
         allocate(cpv%gval(6,cp1%nprop))
         allocate(cpv%dgval(3,nc,cp1%nprop))
         allocate(cpv%d2gval(nc*(nc+1)/2,cp1%nprop))
         cpv%listprop=0
         cpv%amfu=zero
         cpv%dgm=zero
         cpv%phstate=PHENTERED
         cpv%netcharge=zero
         cpv%gval=zero
         cpv%dgval=zero
         cpv%d2gval=zero
! copy the disordered fraction record, that should take care of all
! array allocations inside the disfra record ???
         cpv%disfra=cp1%disfra
!         disordered: if(cpv%disfra%varreslink.gt.0) then
! if there is a disordered phase_varres record that must be taken care of
!            lokdis=cpv%disfra%varreslink
!            eqlista(ieq)%phase_varres(lokdis)%abnorm=&
!                 eqlista(1)%phase_varres(lokdis)%abnorm
! !!!! WOW it really seems to copy a whole tructure just by = !!!
!            eqlista(ieq)%phase_varres(lokdis)=eqlista(1)%phase_varres(lokdis)
! BUT THEN I HAVE TO CHANGE EVERYTHING ABOVE ... NEXT RELEASE ...
!            write(*,*)'3B copied dis: ',lokdis
!            write(*,77)eqlista(ieq)%phase_varres(lokdis)%yfr(2),&
!                 eqlista(1)%phase_varres(lokdis)%yfr(2)
!77          format('enter eq: ',2(1pe15.6))
!            continue
!         endif disordered
      enddo copypv
!      write(*,*)'3B enter_eq 2E, after this segmentation fault'
   endif alleq
! From here also for first equilibria
900 continue
!   write(*,*)'3B enter_eq 3'
   if(ocv()) write(*,*)'3B: entereq 3: '
! nullify condition links, otherwise "if(associated(..)" does not work
   nullify(eqlista(ieq)%lastcondition)
   nullify(eqlista(ieq)%lastexperiment)
   if(ocv()) write(*,*)'3B set T and P',ieq
! also set default local values of T and P (not conditions)
   eqlista(ieq)%tpval(1)=1.0D3; eqlista(ieq)%tpval(2)=1.0D5
! allocate and copy tpfun result array also for first equilibria
!   jz=size(firsteq%eq_tpres)
   jz=maxtpf
!   write(*,*)'3B enter_eq 4',jz,maxsvfun
   if(ocv()) write(*,*)'3B: entereq 4: ',jz,maxsvfun
!    write(*,*)'3B create equil tpres size ',jz,notpf()
! Valgrind wants us to initiate eq_tpres%forcenewcalc !!!
! This is probably quite messy as eq_pres are pointers???
!! eq_tpres already allocated in gtp_init???
!   allocate(eqlista(ieq)%eq_tpres(jz))
   if(.not.allocated(eqlista(ieq)%eq_tpres)) then
!      if(ieq.ne.1) then
!         write(*,*)'3B Allocating eq_tpres for equil: ',ieq,jz,freetpfun
         allocate(eqlista(ieq)%eq_tpres(jz))
!      endif
   endif
! this should be done in init_tpfun (gtp3Z.F90) ??
   do iz=1,jz
      eqlista(ieq)%eq_tpres(iz)%forcenewcalc=0
   enddo
! allocate result array for state variable functions (svfunres)
   if(ocv()) write(*,*)'3B maxsvfun: ',ieq,maxsvfun,jz
!   write(*,*)'3B Allocating svfunres for equilibrium: ',name(1:len_trim(name))
   allocate(eqlista(ieq)%svfunres(maxsvfun))
! convergence criteria PHTUPX
   eqlista(ieq)%xconv=firsteq%xconv
   eqlista(ieq)%maxiter=firsteq%maxiter
1000 continue
   if(ocv()) write(*,*)'3B finished enter equilibrium',ieq
   return
 end subroutine enter_equilibrium !allocate

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine enter_many_equil(cline,last,pun)
! executes an enter many_equilibria command
! and creates many similar equilibria from a table
! pun is file units for storing experimental dataset, pun(i)>0 if i is open
   implicit none
   character*(*) cline
   integer last,pun(9)
!\end{verbatim}
! enter many_equilibria
! by default all phases suspended
! 1 entered phases <list>
! 2 fixed phases <list>
! 3 dormant phases <list>
! 4 conditions ....
! 5 experiments ....
! 6 calculate symbols <list>
! 7 list state_variables <list> 
! 8 table_start
! <equil name> values in columns ...
! 9 table_end
!10 referece state
!11 plot_data
!12 not used
!
! values required by @<column> will appear in table in column order
! EXAMPLE:
! enter many_equilibria
! fixed 1 liquid @1
! condition T=1000 p=1e5
! experiment x(liq,cr)=@2:@3 x(@1,cr)=@4:10%
! table_start
! <equil name> bcc 0.15 0.02 0.20
! ...
! table_end
! expanded experiment line:
! experiment x(liq,cr)=0.15:0.02 x(bcc,cr)=0.20:10%
!
! ncom is numbor of command, ncol is max number of columns in tables
   integer, parameter :: ncom=12,ncol=9
   character (len=12), dimension(ncom), parameter :: commands=&
!        123456789.12...123456789.12...123456789.12...123456789.12
       ['FIXED       ','ENTERED     ','DORMANT     ','CONDITIONS  ',&
        'EXPERIMENTS ','CALCULATE   ','LIST        ','TABLE_START ',&
        'COMMENT     ','REFERENCE_S ','PLOT_DATA   ','            ']
   character*128 rowtext(ncom),text*128,dummy*128,tval*24,savetitle*24
   character*128 eqlin(ncom),eqname*24,plotdatafile*8,encoded*24
   integer dcom,kom,done(ncom),ip,jp,kp,ival,jval,neq,slen,shift,ieq,nystat
   integer iel,iph,ics,maxcol,jj
   type(gtp_equilibrium_data), pointer ::ceq
   double precision xxx,xxy,pxx,pyy,tpa(2),xarr(6)
! This is to know where to store column values from a row
   TYPE gtp_row
      integer column,position
   end type gtp_row
   type(gtp_row), dimension(ncom,ncol) :: colvar,coleq 
   logical plotfile
!
   done=0
   plotfile=.FALSE.
   do ip=1,ncom-1
      colvar(ip,1)%column=0
      rowtext(ip)=' '
   enddo
   maxcol=0
   dcom=0
100 continue
   call gparcd('Table head line: ',cline,last,5,text,' ',q1help)
   kom=ncomp(text,commands,ncom,last)
   if(kom.le.0) then
      write(kou,110)text(1:len_trim(text))
110   format('Error in subcommand to enter many: ',a)
      gx%bmperr=4278; goto 1000
   endif
! the table_start command means end of head, generate one equilibria per row
   if(kom.eq.8) goto 299
! =================================================================
! the heading is stored in character array rowtext(1..12)
! Keep the whole line, the only thing we handle now are column references
   dcom=dcom+1
   rowtext(dcom)=cline
!=====================================================================
! seach for column indicators @digit (0< digit <=9)
   ip=1
200 continue
!  write(*,*)'3B at 200: ',rowtext(kom)(ip:len_trim(rowtext(kom))),ip
   jp=index(rowtext(dcom)(ip:),'@')
   if(jp.gt.0) then
! only a single digit allowed!!
      ival=ichar(rowtext(dcom)(ip+jp:ip+jp))-ichar('0')
! maxcol is the maximal column referred to in the head
      if(ival.gt.maxcol) maxcol=ival
      if(ival.le.0 .or. ival.gt.9) then
! column 0 is name of equilibrium, not a value
         write(kou,*)ival,rowtext(dcom)(1:jp+1)
210      format('Illegal column for variable: ',i3,': ',a)
      else
         do kp=1,ncol
            if(colvar(dcom,kp)%column.eq.0) then
               if(kp.lt.ncol) colvar(dcom,kp+1)%column=0
               colvar(dcom,kp)%column=ival
               colvar(dcom,kp)%position=ip+jp-1
               goto 250
            endif
         enddo
!         write(kou,240)ncol,dcom,rowtext(dcom)(1:len_trim(rowtext(dcom)))
240      format('More than ',i2,' column variables used in row ',i3/a)
         gx%bmperr=4279; goto 1000
! no problem, continue
250      continue
      endif
      ip=ip+jp
      if(ip.lt.len_trim(rowtext(dcom))) then
         goto 200
      endif
   endif
! force reading next command line from file or keyboard
   last=len(cline)
   goto 100
!
!------------------------------------------------------------
! Now start generating one equilibrium per line in table
299 continue
   neq=0
300 continue
! we must not destroy the values in colvar and rowtext!!
   coleq=colvar
   eqlin=rowtext
   call gparc('Table row: ',cline,last,5,text,' ',q1help)
! allow empty lines
   if(len_trim(text).le.1) goto 300
! remove TAB characters
   call untab(text)
! make all upper case
   call capson(text)
!   write(*,*)'3B 300: ',cline(1:len_trim(cline))
   if(text(1:5).eq.'TABLE') then
! finish if first word on line is "TABLE" meaning TABLE_END
      write(kou,310)neq
310   format('Created ',i5,' equilibria')
      goto 1000
   endif
! values are in column order,the digit after @
   ip=0
   values: do ival=0,maxcol
! value in column ival should replace all @digit in all lines, allow "," in tval
      call getext(text,ip,2,tval,' ',slen)
!      write(*,*)'3B tval: ',tval,slen,ival
      if(slen.le.0) then
         write(kou,*)'Table row missing value in column: ',ival
         gx%bmperr=4280; goto 1000
      endif
! first value, in column 0, is equilibrium name
      if(ival.eq.0) then
         eqname=tval; cycle values
      endif
! the column value can be used in several places, also in the same row
      com2: do jp=1,ncom-1
         shift=0
         com3: do kp=1,ncol
            if(coleq(jp,kp)%column.gt.0) then
!               write(*,330)'3B replace: ',jp,kp,coleq(jp,kp)%column,ival,&
!                    shift,tval
330            format(a,2i3,i14,2i4,': ',a)
               if(coleq(jp,kp)%column.eq.ival) then
! insert column value at coleq(jp,kp)%position
                  dummy=eqlin(jp)(coleq(jp,kp)%position+2:)
                  eqlin(jp)(coleq(jp,kp)%position:)=tval
                  eqlin(jp)(coleq(jp,kp)%position+slen:)=dummy
!                  write(*,*)'3B eqlin: ',eqlin(jp)(1:len_trim(eqlin(jp)))
                  shift=shift+slen-2
               else
! we must update all following positions in coleq(jp,...)
!                  write(*,332)'3B shifting: ',jp,kp,coleq(jp,kp)%position,shift
332               format(a,2i3,2x,2i4)
               coleq(jp,kp)%position=coleq(jp,kp)%position+shift
               endif
            else
               cycle com2
            endif
         enddo com3
      enddo com2
   enddo values
! check the final equilibrium description
   neq=neq+1
!   write(*,*)'3B New equilibrium: ',eqname,dcom
!   do kom=1,dcom
!      write(kou,340)neq,kom,eqlin(kom)(1:len_trim(eqlin(kom)))
!340   format('3B cc',2i3,' :',a)
!   enddo
!========================================================================
! create the equilibrium using the row values
!   write(*,*)'3B enter equilibrium: ',eqname,ieq
   call enter_equilibrium(eqname,ieq)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3B entered equilibrium: ',eqname
   call selecteq(ieq,ceq)
!   write(kou,515)eqname,ieq
!515 format('3B Entered equilibrium: ',a,' with number ',i4)
! by default set all phases suspended
   ip=-1; jp=1; nystat=PHSUS; xxx=zero
!   write(*,*)'3B suspending all phases'
   call change_phase_status(ip,jp,nystat,xxx,ceq)
!   call change_many_phase_status(tval,nystat,xxx,ceq)
   if(gx%bmperr.ne.0) goto 1000
!========================================================================
! now set values for the equilibrium description with dcom lines
   do jval=1,dcom
      kom=ncomp(eqlin(jval),commands,ncom,last)
!      write(*,12)'3B eqlin: ',jval,trim(eqlin(jval)),last,kom
12    format(a,i3,' "',a,'" ',2i3)
      SELECT CASE(kom)
!---------------------------
      CASE DEFAULT
         write(*,*)'Error generating equilibrium: ',trim(eqlin(jval))
!---------------------------
      CASE(1,2)! fixed and entered phases
! pick up the number of moles of the phases as first argument after command
         call getext(eqlin(jval),last,1,tval,'1.0',slen)
         ip=1
         call getrel(tval,ip,xxx)
         if(buperr.ne.0) then
            write(*,11)'3B Line causing error: ',trim(eqlin(jval))
11          format(/a,' "',a,'"'/)
            gx%bmperr=4281; goto 1000
         endif
         nystat=PHFIXED
         if(kom.eq.2) nystat=PHENTERED
         if(eolch(eqlin(jval),last)) then
            write(*,*)'3B no phase name after status command'
            gx%bmperr=4282; goto 1000
         endif
         call change_many_phase_status(eqlin(jval)(last:),nystat,xxx,ceq)
         if(gx%bmperr.ne.0) goto 1000
!---------------------------
      CASE(3)! domant phases
         nystat=PHSUS
         xxx=zero
         call change_many_phase_status(eqlin(jval)(last:),nystat,xxx,ceq)
         if(gx%bmperr.ne.0) goto 1000
!---------------------------
      CASE(4)! conditions
         ip=0
         call set_condition(eqlin(jval)(last:),ip,ceq)
         if(gx%bmperr.ne.0) goto 1000
!---------------------------
      CASE(5)! experiments
         ip=0
         call enter_experiment(eqlin(jval)(last:),ip,ceq)
         if(gx%bmperr.ne.0) goto 1000
!---------------------------
      CASE(6)! calculate symbol
         if(.not.allocated(ceq%eqextra)) then
            allocate(ceq%eqextra(3))
            ceq%eqextra(2)=' '
            ceq%eqextra(3)=' '
         endif
         ceq%eqextra(1)=eqlin(jval)(last:)
!---------------------------
      CASE(7)! list state variables and modelled properties
         if(.not.allocated(ceq%eqextra)) then
            allocate(ceq%eqextra(3))
            ceq%eqextra(1)=' '
            ceq%eqextra(3)=' '
         endif
         ceq%eqextra(2)=eqlin(jval)(last:)
!---------------------------
!      CASE(8)! table start should never occur
!---------------------------
      CASE(9)! comment
         ceq%comment=eqlin(jval)(last:)
!---------------------------
      CASE(10)! reference state
         call gparc('Component name: ',eqlin(jval),last,1,tval,' ',q1help)
         call find_component_by_name(tval,iel,ceq)
         if(gx%bmperr.ne.0) goto 1000
         call gparc('Reference phase: ',eqlin(jval),last,1,tval,'SER ',q1help)
         if(tval(1:4).eq.'SER ') then
!            write(kou,*)'Reference state is stable phase at 298.15 K and 1 bar'
! this means no reference phase, SER is at 298.15K and 1 bar
            iph=-1
         else
            call find_phase_by_name(tval,iph,ics)
            if(gx%bmperr.ne.0) goto 1000
! temperature * means always to use current temperature
            xxy=-one
            call gparr('Temperature: /*/: ',eqlin(jval),last,xxx,xxy,q1help)
            if(buperr.ne.0) then
!               write(*,*)'3B buperr: ',buperr
               buperr=0
               tpa(1)=-one
            elseif(xxx.le.zero) then
               tpa(1)=-one
            else
               tpa(1)=xxx
            endif
            xxy=1.0D5
            call gparrd('Pressure: ',eqlin(jval),last,xxx,xxy,q1help)
            if(xxx.le.zero) then
               tpa(2)=xxy
            else
               tpa(2)=xxx
            endif
         endif
!         write(*,*)'3B Reference T and P: ',tpa
         call set_reference_state(iel,iph,tpa,ceq)
!---------------------------
      CASE(11)! PLOT_DATA
         call getint(eqlin(jval),last,ip)
         if(buperr.ne.0) then
            write(kou,*)'Dataset number must be 1 to 9',buperr
         elseif(ip.eq.0) then
! this is a special plotdata file for calculated values, store in eqextra(3)
            if(.not.allocated(ceq%eqextra)) then
               allocate(ceq%eqextra(3))
               ceq%eqextra(1)=' '
               ceq%eqextra(2)=' '
            endif
            ceq%eqextra(3)=' 0 '//eqlin(jval)(last:)
!            write(*,*)'3B eqextra(3): ',trim(ceq%eqextra(3))
         else
            if(ip.le.0 .or. ip.gt.9) then
               write(*,*)'3B plot_data dataset must be from 1 to 9'
               goto 1000
            else
! this is for plot datafile 1 to 9
               if(pun(ip).eq.0) then
                  pun(ip)=30+ip
                  plotdatafile='oc_many0'
                  plotdatafile(8:8)=char(ichar('0')+ip)
                  write(*,*)'3B Opening ',plotdatafile//'.plt'
                  open(pun(ip),file=plotdatafile//'.plt',access='sequential',&
                       status='unknown')
                  call getrel(eqlin(jval),last,pxx)
                  call getrel(eqlin(jval),last,pyy)
                  call getint(eqlin(jval),last,iel)
                  if(buperr.ne.0) then
                     write(*,*)'3B Incorrect values in plot_data',&
                          trim(eqlin(jval))
                     buperr=0
                  endif
                  if(eolch(eqlin(jval),last)) then
                     savetitle='Unknown'
                  else
                     savetitle=trim(eqlin(jval)(last:))
                  endif
                  write(pun(ip),600)iel,trim(eqlin(jval)(last:))
600               format('# GUNPLOT file generated by enter many_equilibria '/&
     'set title "Open Calphad 4.0 : with GNUPLOT"'/&
     'set xlabel "whatever"'/&
     'set ylabel "whatever"'/&
     'set key bottom right'/&
     '#'/'# One can use expressions to convert values'/&
     '# (1-$3) above means the x-value will be, "1-value in column 3"'/'#'/&
     '# pt pointtype 1 +, 2 x, 3 star, 4 square, 5 fill square, 6 circle',/&
     '#',14x,'7 filled circle, 8 triangle up, 9 filled triangle up'/&
     '#',14x,'10 triangle down, 11 filled triangle down, 12 romb'/&
     '#',14x,'13 filled romb, 14 pentad, 15 filled pentad, 16 same as 1 etc'/&
     '# ps pointsize'/'# color ???'/&
     '# set style line 1 lt 2 lc rgb "#000000" lw 2'/&
     '# set multiplot ?? '/'#'/&
     '# To make a nice plot with different symbols for each experimentalist'/&
     '# Move all plot lines after the first one and remove the #'/&
     '# and add a ,\ at the end except for the last.'/&
     '# Remove also the # for the line with a single e',/'#'/&
     'plot "-" using 2:3 with points pt ',i3,&
     ' ps 1.5 title "',a,' and maybe others"')
               else
! This is for plot dataset when the file is open
                  call getrel(eqlin(jval),last,pxx)
                  call getrel(eqlin(jval),last,pyy)
                  call getint(eqlin(jval),last,iel)
                  if(buperr.ne.0) then
                     write(*,*)'3B Incorrect values in plot_data',trim(eqlin(jval))
                     buperr=0
                  endif
                  if(eolch(eqlin(jval),last)) then
                     savetitle='Unknown'
                  endif
                  if(trim(savetitle).ne.trim(eqlin(jval)(last:))) then
! new title, this should preferably be with the previous plot command !!
                     write(pun(ip),605)iel,trim(eqlin(jval)(last:))
605                  format('#e'/'#"" using 2:3 with points pt ',i3,&
                          ' ps 1.5 title "',a,'"')
                     savetitle=trim(eqlin(jval)(last:))
                  endif
!                  write(*,610)'3B debug: ',ip,pxx,pyy,iel,trim(savetitle)
                  write(pun(ip),610)' ',ip,pxx,pyy,iel
610               format(a,i2,2x,2(1pe14.6),i3,5x,a)
               endif
! endif open files and extracting data
            endif
! endif data_plot type
         endif
! endif buperr
!---------------------------
!      CASE(12)! unused
!         continue
      end SELECT
   enddo
!
! force reading next row with values for another equilibrium
   last=len(cline)
   goto 300
!
1000 continue
! we can have many enter many with plot data, do not close here!
! The file(s) will be closed when the command  enter range
!   if(plotfile) then
!      write(pun,1010)
!1010  format('e'/'pause mouse'/)
!      close(pun)
!   endif
   return
 end subroutine enter_many_equil

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine delete_all_conditions(mode,ceq)
! deletes the (circular) list of conditions in an equilibrium
! it also deletes any experiments
! if mode=1 the whole equilibrium is removed, do not change phase status
! because the phase_varres records have been deallocated !!!
! I am not sure it releases any memory though ...
   implicit none
   integer mode
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   type(gtp_condition), pointer :: last,current,next
   integer iph,ics,lokcs
!
!   write(*,*)'3B deleting conditions and experiments',trim(ceq%eqname)
   last=>ceq%lastcondition
   do while(associated(last))
      next=>last%next
      do while(.not.associated(next,last))
         current=>next
         next=>current%next
! if mode=0 then the equilibrium is not deleted, just the conditions
         if(mode.eq.0 .and. current%active.eq.0) then
! if condition is active and that a phase is fix change the phase status!!
! A fix phase has a negative statevariable-id
            iph=-current%statvar(1)%statevarid
!            write(*,*)'3B Active condition: ',iph
            if(iph.gt.0) then
!               write(*,*)'3B rest status for phase: ',iph
               ics=current%statvar(1)%compset
110            continue
               if(phasetuple(iph)%compset.ne.ics) then
                  iph=phasetuple(iph)%nextcs
                  if(iph.gt.0) goto 110
! this composition set does not exist
                  gx%bmperr=4399; goto 1000
               else
                  lokcs=phasetuple(iph)%lokvares
! set the phase status to entered and unknown
!                  write(*,*)'3B remove phase condition: ',iph,ics,lokcs
                  ceq%phase_varres(lokcs)%phstate=0
               endif
            endif
!         else
!            write(*,*)'3B inactive condition: ',current%statvar(1)%statevarid
         endif
         deallocate(current)
      enddo
!      write(*,*)'3B last condition'
      if(mode.eq.0 .and. last%active.eq.0) then
! if condition is active and that a phase is fix change the phase status!!
! A fix phase has a negative statevariable-id
         iph=-last%statvar(1)%statevarid
!         write(*,*)'3B Active condition: ',iph
         if(iph.gt.0) then
!            write(*,*)'3B restore status for phase: ',iph
            ics=last%statvar(1)%compset
120         continue
            if(phasetuple(iph)%compset.ne.ics) then
               iph=phasetuple(iph)%nextcs
               if(iph.gt.0) goto 120
! this composition set does not exist
               gx%bmperr=4399; goto 1000
            else
               lokcs=phasetuple(iph)%lokvares
! set the phase status to entered and stable (not fix)
!               write(*,*)'3B change phase status: ',iph,ics,lokcs
               ceq%phase_varres(lokcs)%phstate=phentstab
!               write(*,*)'3B new phase status: ',&
!                    ceq%phase_varres(lokcs)%phstate
            endif
         endif
      endif
!      write(*,*)'3B deallocate last condition'
      deallocate(last)
!      write(*,*)'3B last condition deallocated'
   enddo
   nullify(ceq%lastcondition)
!------------------------------
! same for experiments (no fix phases)
   last=>ceq%lastexperiment
   do while(associated(last))
      next=>last%next
      do while(.not.associated(next,last))
         current=>next
         next=>current%next
         deallocate(current)
      enddo
      deallocate(last)
   enddo
   nullify(ceq%lastexperiment)
! same for experiments ...
1000 continue
! mark conditions and current result may not be compatible
   ceq%status=ibset(ceq%status,EQINCON)
   return
 end subroutine delete_all_conditions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine delete_equilibria(name,ceq)
! deletes equilibria (needed when repeated step/map)
! name can be an abbreviation line "_MAP*"
! deallocates all data.  Minimal checks ... one cannot delete "ceq"
   implicit none
   character name*(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   type(gtp_equilibrium_data), pointer :: curceq
   type(gtp_condition), pointer :: lastcond,pcond,qcond
   integer cureq,ieq,ik,novarres,ipv
!
   cureq=ceq%eqno
!   write(*,*)'In delete_equilibria ',cureq,trim(name)
   ik=index(name,'*')-1
   if(ik.lt.0) ik=min(24,len(name))
   do ieq=eqfree-1,2,-1
! we cannot have "holes" in the free list??  Delete from the end...
      if(ieq.eq.cureq) exit 
      if(eqlista(ieq)%eqname(1:ik).ne.name(1:ik)) exit
!      write(*,*)'3B Deleting equil: ',trim(eqlista(ieq)%eqname),ieq
      eqlista(ieq)%eqname=' '
      deallocate(eqlista(ieq)%complist)
      deallocate(eqlista(ieq)%compstoi)
      deallocate(eqlista(ieq)%invcompstoi)
      deallocate(eqlista(ieq)%cmuval)
!
      novarres=highcs
! the next line should be removed when highcs implemented
      novarres=csfree-1
!      write(*,*)'3B deallocationg phase_varres'
      do ipv=1,novarres
! can happen it is not allocated when previous errors
         if(.not.allocated(eqlista(ieq)%phase_varres(ipv)%yfr)) cycle
         deallocate(eqlista(ieq)%phase_varres(ipv)%yfr)
         deallocate(eqlista(ieq)%phase_varres(ipv)%constat)
         if(allocated(eqlista(ieq)%phase_varres(ipv)%mmyfr)) &
              deallocate(eqlista(ieq)%phase_varres(ipv)%mmyfr)
         deallocate(eqlista(ieq)%phase_varres(ipv)%sites)
         deallocate(eqlista(ieq)%phase_varres(ipv)%listprop)
         deallocate(eqlista(ieq)%phase_varres(ipv)%gval)
         deallocate(eqlista(ieq)%phase_varres(ipv)%dgval)
         deallocate(eqlista(ieq)%phase_varres(ipv)%d2gval)
! do not deallocate explicitly disfra as it is another phase_varres record ...
      enddo
      deallocate(eqlista(ieq)%phase_varres)
      deallocate(eqlista(ieq)%eq_tpres)
!      write(*,*)'3B Deallocating svfunres for equilibrium:',trim(name)
      deallocate(eqlista(ieq)%svfunres)
! this deletes the conditions and experiments (if any)
      curceq=>eqlista(ieq)
      call delete_all_conditions(1,curceq)
      if(gx%bmperr.ne.0) then
         write(kou,800)gx%bmperr,ieq
800      format(' *** Error ',i6,' deleting equilibrium ',i5)
         gx%bmperr=0
      endif
   enddo
! we have deleted all equilibria until ieq+1
   if(ocv()) write(*,900)ieq+1,eqfree-1
!   write(*,900)ieq+1,eqfree-1
900 format('Deleted all data in equilibria from ',i3,' to ',i3)
   eqfree=ieq+1
1000 continue
   return
 end subroutine delete_equilibria

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine copy_equilibrium(neweq,name,ceq)
! creates a new equilibrium which is a copy of ceq.  
   implicit none
   character name*(*)
   type(gtp_equilibrium_data), pointer ::neweq,ceq
!\end{verbatim} %+
   integer number
   call copy_equilibrium2(neweq,number,name,ceq)
1000 continue
   return
 end subroutine copy_equilibrium

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine copy_equilibrium2(neweq,number,name,ceq)
! creates a new equilibrium which is a copy of ceq.  
! Allocates arrayes for conditions,
! components, phase data and results etc. from equilibrium ceq
! returns a pointer to the new equilibrium record
! THIS CAN PROBABLY BE SIMPLIFIED, especially phase_varres array can be
! copied as a whole, not each record structure separately ... ???
   implicit none
   character name*(*)
   integer number
   type(gtp_equilibrium_data), pointer ::neweq,ceq
!\end{verbatim}
   type(gtp_condition), pointer :: oldcond,lastcond
   type(gtp_condition), pointer :: newcond1,newcond2
   type(gtp_condition), pointer :: bugcond
   character name2*64
   integer ieq,ipv,jz,iz,jl,jk,novarres,oldeq
   logical okname
!
!   write(*,*)'In copy_equilibrium',len_trim(name)
   nullify(neweq)
   if(.not.allowenter(3)) then
!      write(*,*)'3B Not allowed enter a copy'
      gx%bmperr=4153; goto 1000
   endif
!   write(*,*)'3B allow enter OK'
! not allowed to enter equilibria if there are no phases
!   if(btest(globaldata%status,GSNOPHASE)) then
!      write(*,*)'3B Meaningless to copy equilibria with no phase data'
!      gx%bmperr=7777; goto 1000
!   endif
! equilibrium names starting with _ are automatically created by mapping
! and in some other cases.
   if(name(1:1).eq.'_') then
      name2=name(2:)
      jk=1
   elseif(name(1:1).eq.' ') then
      write(*,*)'A name must start with a letter'
      gx%bmperr=4284; goto 1000
   else
      name2=name
      jk=0
   endif
   call capson(name2)
!   write(*,*)'3B Entering copy equilibria: ',name2,jk
! program crashed with this construction
!   if(.not.proper_symbol_name(name2,0)) then
   okname=proper_symbol_name(name2,0)
   if(.not.okname) then
! the name must start with a letter A-Z and contain letters, numbers and _
      gx%bmperr=4122
      goto 1000
   endif
!   write(*,*)'3B name check ok: ',jk
! remove initial "_" used for automatically created equilibria
   if(jk.eq.1) then
! changing this cause a lot of trouble ... but I do not understand
      name2='_'//name2
!      name2=name2(2:)
   endif
! check if name already used
!   write(*,*)'3B check if name unique: ',name2
   call findeq(name2,ieq)
   if(gx%bmperr.eq.0) then
      gx%bmperr=4123
      goto 1000
   else
! reset error code
      gx%bmperr=0
   endif
!   write(*,*)'3B check if name unique: ',eqfree
   if(eqfree.le.maxeq) then
      ieq=eqfree
      eqfree=eqfree+1
   else
!      write(*,*)'Too many equilibrium required, increase dimension',eqfree
      gx%bmperr=4283; goto 1000
   endif
   number=ieq
   if(ieq.eq.1) then
!      write(*,*)'Cannot copy to default equilibria'
      gx%bmperr=4285; goto 1000
   endif
!   write(*,*)'3B copy eq',eqfree,maxeq,ieq
! allocate data arrayes in equilibrium record
   eqlista(ieq)%next=0
   eqlista(ieq)%eqname=name2
   eqlista(ieq)%eqno=ieq
! do not copy comment but set it to blanks
   eqlista(ieq)%comment=' '
! component list and matrix, if second or higher equilibrium copy content
!   write(*,*)'3B: entereq 1A: ',maxel,noofel
   allocate(eqlista(ieq)%complist(noofel))
   allocate(eqlista(ieq)%compstoi(noofel,noofel))
   allocate(eqlista(ieq)%invcompstoi(noofel,noofel))
   allocate(eqlista(ieq)%cmuval(noofel))
!   write(*,*)'3B: entereq 1B: ',noofel
! careful here because FIRSTEQ has other dimensions than the other
   do jl=1,noofel
      eqlista(ieq)%complist(jl)=ceq%complist(jl)
      eqlista(ieq)%cmuval(jl)=ceq%cmuval(jl)
      do jk=1,noofel
         eqlista(ieq)%compstoi(jk,jl)=ceq%compstoi(jk,jl)
         eqlista(ieq)%invcompstoi(jk,jl)=ceq%invcompstoi(jk,jl)
      enddo
   enddo
   oldeq=ceq%eqno
! what about the weight?
   eqlista(ieq)%weight=ceq%weight
!   write(*,*)'3B copyeq 1: ',ceq%weight,eqlista(ieq)%weight
!   write(*,*)'3B: entereq 2: ',noofel
   do jl=1,noofel
      eqlista(ieq)%complist(jl)%splink=eqlista(oldeq)%complist(jl)%splink
      eqlista(ieq)%complist(jl)%phlink=firsteq%complist(jl)%phlink
      eqlista(ieq)%complist(jl)%status=firsteq%complist(jl)%status
      if(firsteq%complist(jl)%phlink.gt.0) then
! only if there is a defined reference state
         eqlista(ieq)%complist(jl)%refstate=firsteq%complist(jl)%refstate
         eqlista(ieq)%complist(jl)%tpref=firsteq%complist(jl)%tpref
         eqlista(ieq)%complist(jl)%chempot=zero
         do jk=1,noofel
            eqlista(ieq)%compstoi(jl,jk)=firsteq%compstoi(jl,jk)
            eqlista(ieq)%invcompstoi(jl,jk)=firsteq%invcompstoi(jl,jk)
         enddo
         if(.not.allocated(eqlista(ieq)%complist(jl)%endmember)) then
            iz=size(firsteq%complist(jl)%endmember)
            allocate(eqlista(ieq)%complist(jl)%endmember(iz))
            eqlista(ieq)%complist(jl)%endmember=firsteq%complist(jl)%endmember
         endif
      else
         eqlista(ieq)%complist(jl)%refstate=firsteq%complist(jl)%refstate
      endif
   enddo
! these records keep calculated values of G and derivatives for each phase
! For phase lokph the index to phase_varres is in phlista(lokph)%cslink
! For phase lokph the index to phase_varres is in phlista(lokph)%linktocs(ics)
! for ieq>1 allocate the current number of phase_varres records plus 10
! for extra composition sets added later
   novarres=highcs
! the next line should be deleted when highcs implemented
   novarres=csfree-1
!   write(*,*)'3B: entereq 3: ',highcs,novarres
! BEWARE: allocation: calculating with one phase with 8 composition sets
! and disordered fractions sets !!!
   iz=max(noofph,novarres)
   allocate(eqlista(ieq)%phase_varres(2*iz))
!   write(*,*)'3B eqlista%phase_varres: ',size(eqlista(ieq)%phase_varres)
! now copy the current content of ceq%phase_varres to this equilibrium
! note, the SELECT_ELEMENT_REFERENCE phase has phase number 0
! and phase_varres index 1, the number of phase_varres records is not the
! same as number of phases ....
!
! strange error here running STEP on bigfcc4: crash with message:
! Index "3" of dimension 1 of array "eqlista" above upper bound of 2
!   write(*,*)'3B 3737:',novarres,ieq,oldeq,size(eqlista(oldeq)%phase_varres)
! Ahhhh, there are 2 phase_varres records for each phase because of 
! disordered fraction set, one for the ordered with 33 y-fractions, one for
! the disordered with 8 y-fractions.  
! A simple dimensioning problem: 1 phase, 8 compsets, disordered fracset
! requires 17 phase_varres.  Before the "max" above I had dimensioned for 2
   copypv: do ipv=1,novarres
      eqlista(ieq)%phase_varres(ipv)=eqlista(oldeq)%phase_varres(ipv)
! in matsmin nprop seemed suddenly to be zero in copied equilibria ....
!      write(*,*)'3B copyeq 2: ',ieq,ipv,eqlista(ieq)%phase_varres(ipv)%nprop
   enddo copypv
900 continue
!   write(*,*)'3B To copy conditions:'
! copy conditions (and experiments) !!!
   lastcond=>eqlista(oldeq)%lastcondition
   if(associated(lastcond)) then
      jz=1
      call copy_condition(eqlista(ieq)%lastcondition,lastcond)
!      write(*,770)'3B cc1: ',jz,lastcond%prescribed,&
!           eqlista(ieq)%lastcondition%prescribed
      newcond1=>eqlista(ieq)%lastcondition
      bugcond=>newcond1
      oldcond=>lastcond%next
      do while(.not.associated(oldcond,lastcond))
         jz=jz+1
         newcond2=>newcond1
         call copy_condition(newcond1%next,oldcond)
         newcond1=>newcond1%next
!         write(*,770)'3B cc2: ',jz,oldcond%prescribed,newcond1%prescribed
770      format(a,i2,6(1pe12.4))
         newcond1%previous=>newcond2
         oldcond=>oldcond%next
      enddo
      newcond1%next=>bugcond
!      write(*,*)'3B Copied all condition',jz
   else
      nullify(eqlista(ieq)%lastcondition)
   endif
! copy experiments) ... later
!
   nullify(eqlista(ieq)%lastexperiment)
!
! copy TPfuns and symbols and current values
!   write(*,*)'3B Copy tpval arrays'
   eqlista(ieq)%tpval=ceq%tpval
   allocate(eqlista(ieq)%eq_tpres(maxtpf))
!   write(*,*)'3B allocated tpres arrays'
   eqlista(ieq)%eq_tpres=ceq%eq_tpres
   allocate(eqlista(ieq)%svfunres(maxsvfun))
!   write(*,*)'3B allocated svfunres arrays'
   eqlista(ieq)%svfunres=ceq%svfunres
! copy convergence criteria
   eqlista(ieq)%xconv=ceq%xconv
   eqlista(ieq)%maxiter=ceq%maxiter
!   write(*,*)'3B finished copy equilibrium',ieq
   eqlista(ieq)%eqno=ieq
   neweq=>eqlista(ieq)
!   write(*,*)'3B Assigned pointer to new equilibrium',neweq%eqno
1000 continue
!   write(*,*)'3B exit copy_equilibrium'
   return
 end subroutine copy_equilibrium2

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine copy_condition(newrec,oldrec)
! Creates a copy of the condition record "oldrec" and returns a link
! to the copy in newrec.  The links to "next/previous" are nullified
   implicit none
   type(gtp_condition), pointer :: oldrec
   type(gtp_condition), pointer :: newrec
!\end{verbatim}
!   write(*,*)' *** In copy_condition:         ',oldrec%prescribed
   allocate(newrec)
!   write(*,*)' *** Allocated'
   newrec=oldrec
!   write(*,*)' *** Copied old condition to new',newrec%prescribed
1000 continue
   return
 end subroutine copy_condition

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 logical function check_minimal_ford(lokph)
! some tests if the fcc/bcc permutation model can be applied to this phase
! The function returns FALSE if the user may set the FORD or BORD bit of lokph
   implicit none
   integer lokph
!\end{verbatim}
   integer nsl,nc,jl,ll,j2,loksp,lokcs
   logical notallowed
   integer, dimension(:), allocatable :: const
   double precision ss
   notallowed=.true.
   nsl=phlista(lokph)%noofsubl
   if(btest(phlista(lokph)%status1,PHHASP)) then
! The PHASP bit is set if a parameter has been entered (never cleared)
      write(kou,*)'Permutation must be set before parameters are entered'
      goto 1000
   endif
   if(nsl.lt.4) then
      write(kou,*)'Phase with permutation must have 4 or more sublattices'
      goto 1000
   else
! ordering assumed in first 4 sublattices, that is not really necessary
!      ss=phlista(lokph)%sites(1)
      lokcs=phlista(lokph)%linktocs(1)
      ss=firsteq%phase_varres(lokcs)%sites(1)
      nc=phlista(lokph)%nooffr(1)
      allocate(const(nc))
      do jl=1,nc
         loksp=phlista(lokph)%constitlist(jl)
         const(jl)=splista(loksp)%alphaindex
      enddo
      jl=nc
      do ll=2,4
!         if(abs(phlista(lokph)%sites(ll)-ss).gt.1.0D-12) then
         if(abs(firsteq%phase_varres(lokcs)%sites(ll)-ss).gt.1.0D-12) then
            write(kou,12)
12          format(' Permutation requires the same number of',&
                 ' sites in first 4 sublattices')
            goto 1000
         endif
         if(phlista(lokph)%nooffr(ll).ne.nc) then
            write(kou,13)
13          format(' Permutation requires that the number of constituents',&
                 ' are equal'/' in all 4 sublattices for ordering')
            goto 1000
         endif
! one must also check the constituents are identical
         do j2=1,nc
            loksp=phlista(lokph)%constitlist(jl+j2)
            if(splista(loksp)%alphaindex.ne.const(j2)) then
               write(kou,14)
14             format(' Permutation requires that the constituents in the',&
                    ' 4 sublattices for'/' ordering are identical')
               goto 1000
            endif
         enddo
         jl=jl+nc
      enddo
   endif
   notallowed=.false.
1000 continue
   check_minimal_ford=notallowed
   return
 end function check_minimal_ford


!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 integer function newhighcs(reserved)
! updates highcd and arranges csfree to be in sequential order
! highcs is the higest used varres record before the last reservation
! or release of a record.  release is TRUE if a record has been released 
! csfree is the beginning of the free list of varres records.
   implicit none
   logical reserved
!\end{verbatim}
   integer high,lok,free,prev
! Do not be smart, go through the whole array
! in all used varres record the %nextfree is zero
   high=0
   free=0
   do lok=1,size(firsteq%phase_varres)
      if(firsteq%phase_varres(lok)%nextfree.eq.0) then
         high=lok
      elseif(free.eq.0) then
! we have the first record belonging to the free list
         free=lok
         prev=lok
      else
         firsteq%phase_varres(prev)%nextfree=lok
         prev=lok
      endif
   enddo
! verification
   prev=2*noofph+2
!   write(*,*)'3B high and free: ',high,free,reserved,highcs,csfree
!   write(*,110)(firsteq%phase_varres(lok)%nextfree,lok=free,prev)
110 format(12(i6))   
   newhighcs=high
1000 continue
 end function newhighcs

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

