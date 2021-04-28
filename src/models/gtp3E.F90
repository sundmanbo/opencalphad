!
! gtp3E included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     9. Save and read things from files
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine gtpsave
!\begin{verbatim}
 subroutine gtpsave(filename,str)
! save all data on file, unformatted, TDB or macro
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
!
   implicit none
   character*(*) filename,str
!\end{verbatim}
! separate UNFORMATTED, DIRECT, TDB, MACRO or LaTeX 
   if(str(1:1).eq.'U') then
! unformatted
      call gtpsaveu(filename,str(3:))
   elseif(str(1:1).eq.'D') then
! direct (random access)
      call gtpsavedir(filename,str(3:))
   elseif(str(1:1).eq.'T') then
! TDB format
      write(*,*)'In gtpsave '
      call gtpsavetdb(filename,str(3:))
   elseif(str(1:1).eq.'L') then
! LaTeX format
      call gtpsavelatex(filename,str(3:))
   else
! macro and PDB formats
      call gtpsavetm(filename,str)
   endif
1000 continue
   return
 end subroutine gtpsave

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine gtpsavelatex
!\begin{verbatim}
 subroutine gtpsavelatex(filename,specification)
! save all data on LaTeX format on a file (for publishing)
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
! equilibrium record(s) with conditions, componenets, phase_varres records etc
! anything else?
   implicit none
   character*(*) filename,specification
!\end{verbatim} %+
1000 continue
   return
 end subroutine gtpsavelatex

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine gtpsavedir
!\begin{verbatim} %-
 subroutine gtpsavedir(filename,specification)
! save all data on a direct file (random access)
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
! equilibrium record(s) with conditions, componenets, phase_varres records etc
! anything else?
   implicit none
   character*(*) filename,specification
!\end{verbatim} %+
1000 continue
   return
 end subroutine gtpsavedir

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine gtpsavetm
!\begin{verbatim}
 subroutine gtpsavetm(filename,str)
! save all data on file in macro format
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
!
   implicit none
   character*(*) filename,str
!\end{verbatim} %+
   logical tdbmode
   if(str(1:1).eq.'T') then
! TDB file
      tdbmode=.true.
   else
! MACRO mode
      tdbmode=.false.
   endif
   write(*,*)'TDB and MACRO save not implemented yet'
   goto 1000
! UNFINISHED ....
! open file and write (either as TDB, MACRO or LaTeX):
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
!
! For inspiration look at the LIST subroutines in pmod25E.F90
!
1000 continue
   return
 end subroutine gtpsavetm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine gtpsavetdb
!\begin{verbatim}
 subroutine gtpsavetdb(filename,specification)
! save all data in TDB format on an file UNFINISHED
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
! equilibrium record(s) with conditions, componenets, phase_varres records etc
! anything else?
   implicit none
   character*(*) filename,specification
!\end{verbatim}
1000 continue
   return
 end subroutine gtpsavetdb

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine gtpsaveu
!\begin{verbatim}
 subroutine gtpsaveu(filename,specification)
! save all data unformatted on an file
! First move it to an integer workspace, then write that on a file
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
! equilibrium record(s) with conditions, componenets, phase_varres records etc
! anything else?
   implicit none
   character*(*) filename,specification
!\end{verbatim} %+
!
   character id*40,comment*72
! size of workspace for unformatted storage
   integer miws
   integer, allocatable, dimension(:) :: iws
   integer i,isp,jph,lokph,lut,last,lok,rsize,displace,ibug,ffun,lokeq,ccc
   integer nspx,check
! these depend on hardware, bytes/word and words/double. Defined in metlib3
!   integer, parameter :: nbpw=4,nwpr=2
! integer function nwch calculates the number of words to store a character
!   write(*,*)'3E In gtpsaveu: ',trim(specification),' version: ',trim(savefile)
!
! positions reserved in the beginning of the workspace
! 3 element list
! 4 element version
! 5 species list
! 6 species version
! 7 tpfun list 
! 8 tpfun version
! 9 phlista lista
! 10 phase version
! 11 endmember version
! 12 interaction version
! 13 property version
! 14 phase tuples lista
! 15 phase tuples version
! 16 equilibrium lista
! 17 equilibrium data version
! 18 component version
! 19 phase_varres version
! 20 global data record
! 21 global data version (not saved?)
! 22 bibliography lista
! 23 bibliography version
! 24 svfun lista
! 25 svfun version
! 26 assessment record list
! 27 assessment version
! 28
! 29
! 30
! missing: parameter_id_lista ... step/map/plot data
! range record? experiments ...
!----------------------------------------------------------------------
! allocate the workspace, words 3-102 for pointers and things listed above
!   write(*,*)'3E allocating iws',miws
   if(btest(globaldata%status,GSNOPHASE)) then
      write(*,1)
1     format('There is no data to save!')
      goto 1001
   endif
! dimension iws depending on number of equuilibria stored
! 7000 is for a 6 component system with 50 phases
! steel1 with 6 elements:         30000 static and       6000 per equilibrium
! TAFID 15 elements:              90000 for static and  30000 per equilibrium
! TAFID 41 elements, 350 phases: 300000 for static and 120000 per equilibrium
! estimate static: nel*nph*x; x=30: static=100000+nel*nph*30
! equilibrium: 40*350*10 ... too small when few elements
   ccc=max(20*noofel*noofph,10000)
! eqfree may not be the higest used equilibrium record!!
   i=eqfree
   last=eqfree
   do while(eqlista(i)%nexteq.gt.0)
      if(eqlista(i)%nexteq.ne.i+1) then
! if eqlista(i)%nexteq does not increment sequentially there are some holes!
         last=eqlista(i)%nexteq
         write(*,*)'3E Beware: unused equilibria before the last used,'//&
              ' cannot be saved'
         gx%bmperr=4399; goto 1000
      endif
      i=eqlista(i)%nexteq
   enddo
!   
   miws=100000+noofel*noofph*30+ccc*last
!   miws=2000000+50000*eqfree
!   write(*,*)'3E allocating workspace: ',miws
   write(*,7)'3E allocating workspace: ',miws,30*noofel*noofph,ccc,last
7  format(a,i10,10x,'(',i7,', ',i7,', ',i4,')')
   allocate(iws(miws))
   call winit(miws,100,iws)
   if(buperr.ne.0) goto 1100
!----------------------------------------------------------------------
! note the use of gtp_xxx_version to handle versions of data structures
!----------------------------------------------------------------------
!
!>>>>> 1: elementlist
!   write(*,*)'3E nbpw and nwpr: ',nbpw,nwpr,nwch(1)
!   rsize=1+1+12/nbpw+24/nbpw+3*nwpr+4 should be enough but ....
   rsize=1+1+nwch(12)+nwch(24)+3*nwpr+5
!   write(*,*)'3E Storing elements',noofel,rsize
   last=3
   iws(last+1)=gtp_element_version
   do i=1,noofel
! next, symbol*2, name*12, ref_state*24, mass, h298, s298,
! splink, status, alphaindex, refstatesymbol
      call wtake(lok,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving element record'
         gx%bmperr=4356; goto 1100
      endif
      call storc(lok+1,iws,ellista(i)%symbol)
      call storc(lok+2,iws,ellista(i)%name)
      displace=3+nwch(12)
      call storc(lok+displace,iws,ellista(i)%ref_state)
      displace=displace+nwch(24)
      call storr(lok+displace,iws,ellista(i)%mass)
      call storr(lok+displace+nwpr,iws,ellista(i)%h298_h0)
      call storr(lok+displace+2*nwpr,iws,ellista(i)%s298)
      displace=displace+3*nwpr
      iws(lok+displace)=ellista(i)%splink
      iws(lok+displace+1)=ellista(i)%status
      iws(lok+displace+2)=ellista(i)%alphaindex
      iws(lok+displace+3)=ellista(i)%refstatesymbol
!      write(*,*)'3E element: ',i,displace+3,rsize,ellista(i)%refstatesymbol
! link sequentially in first word
      iws(last)=lok
      last=lok
      ibug=lok+displace+3
!      write(*,*)'3E refstatesymbol 0: ',ibug,iws(ibug),iws(1)
   enddo
! bug??
! added one saved integer for size of spextra (normally zero)
   ibug=lok+displace+4
!   write(*,*)'3E refstatesymbol 1: ',ibug,iws(ibug),iws(1)
!-----------
!>>>>> 2: specieslist
   rsize=1+nwch(24)+3*nwpr+3
! next, symbol*24, mass, charge, extra, noofel, status, alphaindex
! (allocated) ellinks, stoichiometry
!   write(*,*)'3E storing species',noofsp,rsize,'+ellinks'
   last=5
   iws(last+1)=gtp_species_version
   do isp=1,noofsp
      if(allocated(splista(isp)%spextra)) then
         nspx=size(splista(isp)%spextra)
      else
         nspx=0
      endif
      check=rsize+splista(isp)%noofel*(1+nwpr)+nspx*nwpr
      call wtake(lok,rsize+splista(isp)%noofel*(1+nwpr)+nspx*nwpr,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving species record'
         gx%bmperr=4356; goto 1100
      endif
!      write(*,*)'3E refstatesymbol 2: ',ibug,iws(ibug),lok
      call storc(lok+1,iws,splista(isp)%symbol)
      displace=2+nwch(24)
      call storr(lok+displace,iws,splista(isp)%mass)
      call storr(lok+displace+nwpr,iws,splista(isp)%charge)
      iws(lok+displace+2*nwpr)=splista(isp)%noofel
      iws(lok+displace+2*nwpr+1)=splista(isp)%status
      iws(lok+displace+2*nwpr+2)=splista(isp)%alphaindex
      iws(lok+displace+2*nwpr+3)=nspx
! displace one less as the index i is added
      displace=displace+2*nwpr+3
      do i=1,splista(isp)%noofel
         iws(lok+displace+i)=splista(isp)%ellinks(i)
      enddo
      displace=displace+splista(isp)%noofel+1
! storing splista(isp)%noofel doubles in iws(lok+displace)
!      write(*,*)'3E displace store: ',lok,displace
! storrn starts storing in iws(lok+displace)
      call storrn(splista(isp)%noofel,&
           iws(lok+displace),splista(isp)%stoichiometry)
!  if nspx>0 save also all double variables in spextra
      if(nspx.gt.0) then
         displace=displace+splista(isp)%noofel*nwpr
         call storrn(nspx,iws(lok+displace),splista(isp)%spextra)
!         write(*,*)'3E species with extra data: ',isp,nspx
      endif
!      write(*,'(a,2i5)')'3E species record check: ',check,&
!           displace+nspx*nwpr
!      write(*,*)'3E refstatesymbol 3: ',ibug,iws(ibug),lok+displace
!      do i=1,splista(isp)%noofel
!         call storr(lok+displace+(i-1)*nwpr,iws,
!      enddo
!      write(*,*)'3E stored species ',isp,lok,displace+splista(isp)%noofel*nwpr
! link records sequentially in first word
      iws(last)=lok
      last=lok
   enddo
!   write(*,*)'3E last species link: ',last,lok,iws(lok),iws(1)
!
!------------- tpfuns
!>>>>> 20: tpfuns
   call wtake(lok,freetpfun,iws)
   if(buperr.ne.0) then
      write(*,*)'3E Error reserving tpfun record'
      gx%bmperr=4356; goto 1100
   endif
   iws(7)=lok
   iws(8)=tpfun_expression_version
   iws(lok)=freetpfun
!   write(*,*)'3E saving TPfuns: ',iws(7),iws(iws(7))
   do i=1,freetpfun-1
! store all TPfuns here. In parameters just store an index!
! we have to pass iws also ....
      call save0tpfun(ffun,iws,i)
      if(gx%bmperr.ne.0) goto 1100
!      write(*,*)'3E TPfun: ',i,' stored in ',iws(lok+i),iws(iws(lok+i))
      iws(lok+i)=ffun
   enddo
! write the record for TP function 3 as check
!   call wrttprec(3,iws)
! All seems OK this far
!------------- phases and parameters, static data
!>>>>> 3: phaselist, start from 0 (reference phase)
! including sublattces, endmembers, interactions, properties etc
! save version of various records
!   write(*,*)'3E saving phases',noofph
   last=9
   iws(last+1)=gtp_phase_version
   iws(last+2)=gtp_endmember_version
   iws(last+3)=gtp_interaction_version
   iws(last+4)=gtp_property_version
   call savephases(last,iws)
   if(gx%bmperr.ne.0) goto 1100
! save all phase tuples in a single reord
   last=14
   iws(last+1)=gtp_phasetuple_version
!   write(*,*)'3E Saving tuples: ',nooftuples
   if(nooftuples.gt.0) then
      call wtake(lok,1+nooftuples*5,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving phase tuple record'
         gx%bmperr=4356; goto 1100
      endif
      iws(lok)=nooftuples
      do i=0,nooftuples-1
         iws(lok+5*i+1)=phasetuple(i+1)%lokph
         iws(lok+5*i+2)=phasetuple(i+1)%compset
         iws(lok+5*i+3)=phasetuple(i+1)%ixphase
         iws(lok+5*i+4)=phasetuple(i+1)%lokvares
         iws(lok+5*i+5)=phasetuple(i+1)%nextcs
      enddo
      iws(last)=lok
   else
! no phase tuples
      iws(last)=0
   endif
!   write(*,*)'3E tuples saved: '
!------------------------------------
! save link to the global data record and version in 20-21
   last=20
! extended globaldata record 190317/BoS
   rsize=1+nwch(24)+3*nwpr+11+5*nwpr
   call wtake(lok,rsize,iws)
   if(buperr.ne.0) then
      write(*,*)'3E Error reserving globaldata record'
      gx%bmperr=4356; goto 1100
   endif
   iws(last)=lok
   iws(lok+1)=globaldata%status
   call storc(lok+2,iws,globaldata%name)
! BUG name was ovewritten by rgas etc !!!
   displace=2+nwch(24)
   call storr(lok+displace,iws,globaldata%rgas)
   call storr(lok+displace+nwpr,iws,globaldata%rgasuser)
   call storr(lok+displace+2*nwpr,iws,globaldata%pnorm)
! extended globaldata record 190317/BoS
   displace=displace+3*nwpr
! these used for testing when reading
!   globaldata%sysparam(1)=987
!   globaldata%sysparam(10)=17
   do i=0,9
      iws(lok+displace+i)=globaldata%sysparam(i+1)
   enddo
   displace=displace+10
!   globaldata%sysreal(5)=12345678.9D0
   call storrn(5,iws(lok+displace),globaldata%sysreal)
!   write(*,*)'3E globalsave:: ',rsize,displace+5*nwpr
!   write(*,*)'3E name: "',globaldata%name,'"'
!   goto 900
! unfinished
!------------- state variable functions
!>>>>> 30: svfuns
!   write(*,*)'3E saving state variable functions'
   call svfunsave(lok,iws,firsteq)
   if(gx%bmperr.ne.0) goto 1100
   iws(24)=lok
   iws(25)=gtp_putfun_lista_version
!------------- references
!>>>>> 40: bibliographic references
!   write(*,*)'3E saving bibliography'
! link to bibliography is stored in 22
   call bibliosave(lok,iws)
   if(gx%bmperr.ne.0) goto 1100
   iws(22)=lok
   iws(23)=gtp_biblioref_version
! document use of workspace
   call wrkchk(rsize,miws,iws)
   write(*,*)'3E used ',rsize,' words out of ',miws,' for storing static data'
!-------------------------------------------------------
! write the equilibrium records
! conditions, components, phase_varres for all composition sets etc
!>>>>> 50: equilibria
!   write(*,*)'3E saving equilibria'
!   write(lut)gtp_equilibrium_data_version,gtp_component_version,&
!        gtp_phase_varres_version
   lokeq=0
! all equilibria are saved here
   call saveequil(lokeq,iws)
   if(gx%bmperr.ne.0) goto 1100
! finished saving equilibria
!   write(*,*)'3E first saved equilibrium at: ',lokeq
   iws(16)=lokeq
   iws(17)=gtp_equilibrium_data_version
   iws(18)=gtp_component_version
   iws(19)=gtp_phase_varres_version
! disfra record version??
!-------------------------------------------------------
! assessment head record
   write(*,*)'3E Saving assessment record'
   if(associated(firstash)) then
      iws(27)=gtp_assessment_version
      lok=26
      call saveash(lok,iws)
      if(gx%bmperr.ne.0) goto 1100
   endif
! free list for phase_varres records
!   write(*,*)'3E Phase_varres first free/highcs: ',csfree,highcs
! UNFINISHED we should write assessment records and step/map/plot records
!-------------------------------------------------------
! finally write the workspace to the file ...
900 continue
   if(index(filename,'.').eq.0) then
      filename(len_trim(filename)+1:)='.OCU'
   endif
   lut=21
!**********************************************************
! IMPORTANT savefile
! is a variable in gtp3.F90
! which MUST BE CHANGED whenever there is a change in the unformatted
! file layout
!***********************************************************
   open(lut,file=filename,access='sequential',status='unknown',&
           form='unformatted',iostat=gx%bmperr,err=1000)
   id='This is a save file for OC version:    '
   comment=specification
   call wrkchk(rsize,miws,iws)
! NOTE: savefile is a character*8 in gtp3.F90
   last=5+nwch(40)+nwch(8)+nwch(72)
   write(lut)id,savefile,comment,noofel,noofsp,noofph,nooftuples,rsize+5
   write(lut)(iws(i),i=1,rsize+5)
   close(lut)
   if(buperr.ne.0) then
      write(kou,990)buperr
990   format(/' *** WARNING *** , workspace save error: ',i7/)
   endif
   write(kou,989)rsize+5+last,miws,1.0D2*real(rsize+5+last)/real(miws)
989 format('Used ',i8,' words out of ',i8,', ',&
         F6.2,'% for unformatted memory save')
   write(kou,991)nbpw*(rsize+5+last),trim(filename)
991 format('Written workspace ',i10,' bytes unformatted on ',a)
1000 continue
   deallocate(iws)
1001 continue
   return
1100 continue
   write(*,*)'3E Error storing record, nothing written on file',buperr,gx%bmperr
   gx%bmperr=4357
   goto 1000
 end subroutine gtpsaveu

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine savephases
!\begin{verbatim}
 subroutine savephases(phroot,iws)
! save data for all phases and store pointer in iws(phroot)
! For phases with disordered set of parameters we must access the number of
! sublattices via firsteq
   implicit none
   integer phroot,iws(*)
!\end{verbatim} %+
   integer doneord,i,j,level,lokcs,nem,noi,nop,nox,nsl,nup,noendm,fipsize
   integer iph,lok,rsize,displace,lokph,iwsph,lokem,lastem,lokpty,last
   integer phreclink
   type(gtp_endmember), pointer :: emrec
   type(gtp_interaction), pointer :: intrec
   type(gtp_property), pointer :: proprec
   character*8 dummy
   logical higher
! to keep track of interaction records
   type saveint
      type(gtp_interaction), pointer :: p1
      integer lok
   end type saveint
   type(saveint), dimension(:), pointer :: stack
   type(gtp_phase_add), pointer :: addlink
   allocate(stack(5))
! do not save the phases array, regenerate it on read
!   call wtake(lok,noofph+1,iws)
!   do i=1,noofph
!      iws(lok+i)=phases(i)
!   enddo
! store this link in last and set link to next in first word
!   iws(last)=lok
!   last=lok
! loop for all phases
   iwsph=phroot
! phases start from 0, the SER phase
   do iph=0,noofph
!      lokph=phases(iph)
      lokph=iph
!>>>>> 5: phase header
! link,name*24,model*72,phletter*1,status1,alphaindex,noofcs,nooffs,additionlink
      rsize=1+nwch(24)+nwch(72)+nwch(1)+5
! the endmember_ord, endmember_dis and endmemrecarray is not used ...
! noofsubl,tnooffr,linktocs(9),nooffr(subl),constlist(tnooffr),i2slx
      rsize=rsize+2+9+phlista(lokph)%noofsubl+phlista(lokph)%tnooffr+2
! we must also have links to two endmember lists and addtions
      rsize=rsize+3
      call wtake(lok,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving phase record',trim(phlista(lokph)%name),&
              buperr
         gx%bmperr=4356; goto 1000
      endif
! link all phase records sequentially from phroot using iwsph
      iws(iwsph)=lok
      iwsph=lok
! store phase data
!      write(*,*)'3E creating record for ',trim(phlista(lokph)%name),lok
      call storc(lok+1,iws,phlista(lokph)%name)
      displace=1+nwch(24)
      call storc(lok+displace,iws,phlista(lokph)%models)
      displace=displace+nwch(72)
! we should store at least 4 characters ...
      dummy=phlista(lokph)%phletter
      call storc(lok+displace,iws,dummy(1:4))
      displace=displace+1
      iws(lok+displace)=phlista(lokph)%status1
      iws(lok+displace+1)=phlista(lokph)%alphaindex
      iws(lok+displace+2)=phlista(lokph)%noofcs
      iws(lok+displace+3)=phlista(lokph)%nooffs
! mark there are additions, it is handled below
!      if(associated(phlista(lokph)%additions)) then
!         iws(lok+displace+4)=1
!      endif
      if(allocated(phlista(lokph)%oendmemarr)) then
         write(*,*)'3E Attention!! ignoring endmemberrec array!'
      endif
!>>>>> 6: sublattice and constituent info
      nsl=phlista(lokph)%noofsubl
      iws(lok+displace+4)=nsl
      j=phlista(lokph)%tnooffr
      iws(lok+displace+5)=j
! displace one less as loops starts from 1
      displace=displace+5
      do i=1,9
         iws(lok+displace+i)=phlista(lokph)%linktocs(i)
      enddo
      displace=displace+9
      do i=1,nsl
         iws(lok+displace+i)=phlista(lokph)%nooffr(i)
      enddo
      displace=displace+nsl
      do i=1,j
         iws(lok+displace+i)=phlista(lokph)%constitlist(i)
      enddo
      displace=displace+j+1
! saving i2sl is probably not necessary as it is calculated each time ??
      iws(lok+displace)=phlista(lokph)%i2slx(1)
      iws(lok+displace+1)=phlista(lokph)%i2slx(2)
! links to endmembers and additions to be stored here and afterwards
! iws(phreclink) ordered, iws(phreclink+1) disordered, iws(phreclink+2) addition
      phreclink=lok+displace+2
!      write(*,*)'3E phreclink 1: ',phreclink,iws(1)
!--------- endmember list, interaction tree and property records
! save all parameter data starting from the endmember list
      doneord=0
      emrec=>phlista(lokph)%ordered
!      write(*,*)'3E saving endmembers',doneord,nsl
! there can be phases without any ordered parameters ...
      if(.not.associated(emrec)) goto 400
! The start of the sequentail list of endmember records (for ordered fractions)
      lokem=phreclink
! we come back here if there are disordered parameters
200   continue
! if doneord=1 then we have listed the ordered parameters
      if(doneord.eq.1) then
         emrec=>phlista(lokph)%disordered
         if(ocv()) write(*,*)'3E Saving disordered parameters'
      endif
      emlista: do while(associated(emrec))
         proprec=>emrec%propointer
         intrec=>emrec%intpointer
!         nop=0
!         noi=0
!         nem=0
!         if(associated(proprec)) nop=1
!         if(associated(intrec)) noi=1
!         if(associated(emrec%nextem)) nem=1
!>>>>> 7: endmember record (basic or disordered)
!         write(lut)emrec%noofpermut,emrec%phaselink,emrec%antalem,nop,noi,nem
!         do j=1,emrec%noofpermut
!            write(lut)(emrec%fraclinks(i,j),i=1,nsl)
!         enddo
! In the endmember recorde we store:
! link to next endmember, link to interaction, link to property record : 3
! link to phase record, number of permutations, seq.order of creation? : 3
! for each permutation nsl indices to fractions                      : perm*nsl
!
         rsize=6+emrec%noofpermut*nsl
         call wtake(lok,rsize,iws)
         if(buperr.ne.0) then
            write(*,*)'3E Error reserving endmember record'
            gx%bmperr=4356; goto 1000
         endif
!         write(*,*)'3E emrec:    ',emrec%noofpermut,lok,rsize,emrec%antalem
! maintain the sequential link between all endmember records
         iws(lokem)=lok
         lokem=lok
! iws(lok) to next, iws(nop=lok+1) to property, iws(noi=lok+2) to intercation, 
! these are nem, noi, nop
!         write(lut)emrec%noofpermut,emrec%phaselink,emrec%antalem,nop,noi,nem
         iws(lok+3)=emrec%noofpermut
         iws(lok+4)=emrec%phaselink
         iws(lok+5)=emrec%antalem
         displace=5
         do j=1,emrec%noofpermut
            do i=1,nsl
               iws(lok+displace+i)=emrec%fraclinks(i,j)
            enddo
            displace=displace+nsl
!            write(lut)(emrec%fraclinks(i,j),i=1,nsl)
         enddo
! this is the place to store the start of property records
         nop=lok+1
         level=nop
         emproplista: do while(associated(proprec))
!            if(associated(proprec%nextpr)) nox=1
!>>>>> 8: endmember property record (loop) add place for %modelparamid
!            rsize=5+nwch(16)+proprec%degree+1
            rsize=5+nwch(20)+proprec%degree+1
            call wtake(lokpty,rsize,iws)
            if(buperr.ne.0) then
               write(*,*)'3E Error reserving endmember record'
               gx%bmperr=4356; goto 1000
            endif
! link the property recordds sequentially
            iws(nop)=lokpty
!            write(*,*)'3E endmem property record',iws(nop),lokpty,&
!                 proprec%proptype,proprec%degree
            level=nop
            nop=lokpty
!            write(lut)proprec%reference,proprec%proptype,&
!                 proprec%degree,proprec%extra,proprec%antalprop,nox
            iws(lokpty+1)=proprec%proptype
            iws(lokpty+2)=proprec%degree
            iws(lokpty+3)=proprec%extra
            iws(lokpty+4)=proprec%antalprop
            call storc(lokpty+5,iws,proprec%reference)
            displace=5+nwch(16)
!            write(*,*)'place to save modelparamid 1: ',&
!                 proprec%modelparamid,lokpty+displace
            call storc(lokpty+displace,iws,proprec%modelparamid)
            displace=displace+nwch(4)
            do i=0,proprec%degree
! store a link in iws(lokpty+displace+i) to the TP fun stored as a text
! we have to pass iws also ....
!               call save1tpfun(lut,.FALSE.,proprec%degreelink(i))
! third argument 1 means do not store function name
!               call save2tpfun(lokpty+displace+i,iws,1,proprec%degreelink(i))
!               if(gx%bmperr.ne.0) goto 1000
               iws(lokpty+displace+i)=proprec%degreelink(i)
               if(gx%bmperr.ne.0) goto 1000
!               write(*,*)'3E place of function: ',iws(lokpty+displace+i),&
!                    ' stored in ',lokpty+displace+i,iws(1)
            enddo
            proprec=>proprec%nextpr
         enddo emproplista
! at the end of the propoerty list iws(lokpty)=0 (zero)
! start interaction tree
         level=0
         noi=lokem+2
! return here for new interaction record
300      continue
         intlista: do while(associated(intrec))
! noi is next, nup is higher, nop is property
310         continue
!>>>>> 9: interaction record
! next, higher,property,status,antalint,order,fipsize  :7
! very complex for permutations ...
! look in gtp3G, create_interaction, for use of intrec%noofip
! noofip,sublattice(noofip),fraclink(noofip) 
!            write(*,*)'3E save link: ',intrec%antalint,intrec%noofip(2)
            fipsize=size(intrec%noofip)
! 2020.06.08/BoS problem handling interactions permutations ???
! should there be something separate for level=2 ???
! Sometimes error when reading a parameter level=2 WITHOUT any permutations
            if(level.eq.1) then
               rsize=7+fipsize+2*intrec%noofip(fipsize)
            else
               rsize=7+fipsize+2*intrec%noofip(2)
            endif
            call wtake(lok,rsize,iws)
            if(buperr.ne.0) then
               write(*,*)'3E Error reserving interaction record',&
                    buperr,rsize,fipsize
               gx%bmperr=4356; goto 1000
            endif
! store data
            iws(lok+3)=intrec%status
            iws(lok+4)=intrec%antalint
            iws(lok+5)=intrec%order
            iws(lok+6)=fipsize
            displace=6
            do i=1,fipsize
               iws(lok+displace+i)=intrec%noofip(i)
            enddo
            displace=displace+fipsize
! intrec%noofip(2) is OK for 1st order, for 2nd order we must use
! intrec%noofip(fipsize)
!            write(*,*)'3E fipsize: ',level,fipsize,&
!                 intrec%noofip(2),intrec%noofip(fipsize)
            if(level.ne.1) then
               do i=1,intrec%noofip(2)
                  iws(lok+displace+2*i-1)=intrec%sublattice(i)
                  iws(lok+displace+2*i)=intrec%fraclink(i)
               enddo
            elseif(level.eq.1) then
               do i=1,intrec%noofip(fipsize)
                  iws(lok+displace+2*i-1)=intrec%sublattice(i)
                  iws(lok+displace+2*i)=intrec%fraclink(i)
               enddo
            endif
!            write(*,11)'3E interaction: ',intrec%antalint,higher,lok,noi,&
!                 intrec%noofip(2),intrec%sublattice(1),intrec%fraclink(1)
11          format(a,i3,l3,10i5)
! link from previous, iws(lok+1) is link to higher, iws(lok+2) is property
            iws(noi)=lok
            noi=lok
! interaction property, link from nop
            proprec=>intrec%propointer
            nop=lok+2
            intproplista: do while(associated(proprec))
!>>>>> 10: interaction property record (loop)
!               rsize=5+nwch(16)+proprec%degree+1
               rsize=5+nwch(20)+proprec%degree+1
               call wtake(lokpty,rsize,iws)
               if(buperr.ne.0) then
                  write(*,*)'3E Error reserving inteaction property record'
                  gx%bmperr=4356; goto 1000
               endif
! link the property records sequentially
               iws(nop)=lokpty
               nop=lokpty
!               write(*,*)'3E interact property record',lokpty,&
!                    proprec%proptype
!            write(lut)proprec%reference,proprec%proptype,&
!                 proprec%degree,proprec%extra,proprec%antalprop,nox
               iws(lokpty+1)=proprec%proptype
               iws(lokpty+2)=proprec%degree
               iws(lokpty+3)=proprec%extra
               iws(lokpty+4)=proprec%antalprop
               call storc(lokpty+5,iws,proprec%reference)
               displace=5+nwch(16)
!               write(*,*)'place to save modelparamid 2: ',&
!                    proprec%modelparamid,lokpty+displace
               call storc(lokpty+displace,iws,proprec%modelparamid)
               displace=displace+nwch(4)
               do i=0,proprec%degree
! store a link in iws(lokpty+displace+i) to the TP fun stored as a text
! we have to pass iws also ....
!                  call save2tpfun(lokpty+displace+i,iws,1,&
!                       proprec%degreelink(i)) 
                  iws(lokpty+displace+i)=proprec%degreelink(i)
               enddo
               proprec=>proprec%nextpr 
            enddo intproplista
! save on stack and check if higher level
            level=level+1
            if(level.gt.5) then
!            write(*,*)'3E Too many interaction levels'
               gx%bmperr=4164; goto 1000
            endif
! save this interaction record and take link to higher
            stack(level)%p1=>intrec
            stack(level)%lok=lok
            intrec=>intrec%highlink
! link to higher should be in lok+1
            noi=lok+1
            if(associated(intrec)) higher=.true.
         enddo intlista
! we come here when there is no higher level
! pop previous intrec and take link to next interaction (on same level)
         higher=.false.
         if(level.gt.0) then
            intrec=>stack(level)%p1
            noi=stack(level)%lok
            intrec=>intrec%nextlink
            level=level-1
            goto 300
         endif
!---- next endmember
         emrec=>emrec%nextem
      enddo emlista
! no more endmembers, check if the disordered (if any) has been written
400   continue
! take link to higher higher interaction
      if(doneord.eq.0) then
         if(ocv()) write(*,*)'3E any disordered endmembers?'
         if(associated(phlista(lokph)%disordered)) then
! there are also disordered parameters
! the disfra record is written in saveequil??
! we have to change nsl ...three % vojvoj
            doneord=1
            lokcs=phlista(lokph)%linktocs(1)
            nsl=firsteq%phase_varres(lokcs)%disfra%ndd
!>>>>> 11A: write disordered endmemebers
!            write(lut)2,nsl
! emrec should already be null but for security ....
            nullify(emrec)
            lokem=phreclink+1
            goto 200
         endif
      endif
!------ save additions list, use lokpty...
500 continue
! iws error check
      addlink=>phlista(lokph)%additions
      lokpty=phreclink+2
      addition: do while(associated(addlink))
! WHEN SAVING MORE ADDITION YOU MUST ALSO CHANGE READING UNFORMATTED readphases
!  integer, public, parameter :: INDENMAGNETIC=1
!  integer, public, parameter :: XIONGMAGNETIC=2
!  integer, public, parameter :: DEBYECP=3
!  integer, public, parameter :: EINSTEINCP=4
!  integer, public, parameter :: TWOSTATEMODEL1=5
!  integer, public, parameter :: ELASTICMODEL1=6
!  integer, public, parameter :: VOLMOD1=7
!  integer, public, parameter :: UNUSED_CRYSTALBREAKDOWNMOD=8
!  integer, public, parameter :: SECONDEINSTEIN=9
!  integer, public, parameter :: SCHOTTKYANOMALY=10
!  integer, public, parameter :: DIFFCOEFS=11
! with composition independent G2 parameter
!         if(addlink%type.eq.1) then
         if(addlink%type.eq.INDENMAGNETIC) then
!>>>>> 12A: additions id, regenerate all when reading this
!            rsize=3
! also saving status
            rsize=4
            call wtake(lok,rsize,iws)
            if(buperr.ne.0) then
               write(*,*)'3E Error reserving addition record'
               gx%bmperr=4356; goto 1000
            endif
            iws(lokpty)=lok
            lokpty=lok
            iws(lok+1)=addlink%type
            iws(lok+2)=addlink%aff
            iws(lok+3)=addlink%status
!            write(*,*)'3E saving additions in: ',phreclink+2,lok,iws(lok+1),&
!                 iws(lok+2)
         elseif(addlink%type.eq.XIONGMAGNETIC) then       ! 2
!>>>>> 12A: additions id, regenerate all when reading this
!            rsize=3
! also saving status, there is a real
            rsize=4
            call wtake(lok,rsize,iws)
            if(buperr.ne.0) then
               write(*,*)'3E Error reserving addition record'
               gx%bmperr=4356; goto 1000
            endif
            iws(lokpty)=lok
            lokpty=lok
            iws(lok+1)=addlink%type
! we have no aff but for xiongmagnetic we specify -1 for BCC
!            write(*,*)'3E xionmagnetic: ',addlink%status,ADDBCCMAG
            if(btest(addlink%status,ADDBCCMAG)) then
               iws(lok+2)=-1
            else
               iws(lok+2)=0
            endif
! there is no need to save this because record is will be regenerated
            iws(lok+3)=addlink%status
! addrecord typ 3 not used
! link the property recordds sequentially
         elseif(addlink%type.eq.EINSTEINCP) then                ! 4
!            write(*,*)'Not saving Einstein addition'          
            rsize=4
            call wtake(lok,rsize,iws)
            if(buperr.ne.0) then
               write(*,*)'3E error saving addition record'
               gx%bmperr=4356; goto 1000
            endif
            iws(lokpty)=lok
            lokpty=lok
            iws(lok+1)=addlink%type
            iws(lok+3)=addlink%status
         elseif(addlink%type.eq.TWOSTATEMODEL1) then          ! 5
!            write(*,*)'Not saving liquid two-state addition'
            rsize=4
            call wtake(lok,rsize,iws)
            if(buperr.ne.0) then
               write(*,*)'3E error saving addition record'
               gx%bmperr=4356; goto 1000
            endif
            iws(lokpty)=lok
            lokpty=lok
            iws(lok+1)=addlink%type
            iws(lok+3)=addlink%status
! addrecord typ 6 not used
         elseif(addlink%type.eq.VOLMOD1) then                 ! 7  
!>>>>> 12A: additions id, regenerate all when reading this
!           rsize=3
            rsize=4
            call wtake(lok,rsize,iws)
            if(buperr.ne.0) then
               write(*,*)'3E Error reserving addition record'
               gx%bmperr=4356; goto 1000
            endif
            iws(lokpty)=lok
            lokpty=lok
            iws(lok+1)=addlink%type
! save also the status word
            iws(lok+3)=addlink%status
!            iws(lok+2)=addlink%aff
!            write(*,*)'3E saving additions in: ',phreclink+2,lok,iws(lok+1),&
!                 iws(lok+2)
         elseif(addlink%type.eq.DIFFCOEFS) then               ! 11
            write(*,*)'Not saving Diffusion addition'
         else
            write(*,99)addlink%type
99          format(78('*')/'3E *** NOT SAVED addition record type ',i3/78('*')/)
         endif
         addlink=>addlink%nextadd
      enddo addition
   enddo
!   write(*,*)'3E phreclink 2: ',phreclink,iws(phreclink),iws(phreclink+1),&
!        iws(phreclink+2)
1000 continue
   return
 end subroutine savephases

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine saveequil
!\begin{verbatim}
 subroutine saveequil(lok1,iws)
! subroutine saveequil(lok1,iws,ceq)
! save data for equilibrium record ceq including phase_varres
   implicit none
   integer lok1,iws(*),jeq
!\end{verbatim} %+
   character text*1024
   type(gtp_phase_varres), pointer :: firstvarres
   TYPE(gtp_fraction_set), pointer :: fslink
!   TYPE(gtp_condition), pointer :: condrec
   integer i,isp,j,k,kl,lokcs,lokph,mc,mc2,nsl,lokeq,rsize,displace,lokvares
   integer lokdis,disz,lok,qsize,eqdis,iws1,dcheck,lokcc,seqz,offset,dmc
   integer loklast,eqnumber,lokhighcs,ceqsize,ceqsize2
   type(gtp_equilibrium_data), pointer :: ceq
! loop to save all equilibria
   eqnumber=0
   ceqsize2=ceqrecsize()
17 continue
   eqnumber=eqnumber+1
   if(eqnumber.eq.1) then
! calculate the size of the first equilibrium record saved
      ceqsize=iws(1)
   elseif(eqnumber.eq.2) then
      ceqsize=iws(1)-ceqsize
      write(*,18)ceqsize,ceqsize2
18    format(' 3E Saving an equilibrium record requires ',2i8,' words')
   endif
   ceq=>eqlista(eqnumber)
! check if enything entered ...
   if(.not.allocated(ceq%complist)) then
      write(*,*)'3E not storing unused equilibria from: ',eqnumber
      goto 1000
!   else
!      write(*,*)'3E storing equilibrium number: ',eqnumber
   endif
!>>>>> 50:
!   write(lut)ceq%eqname,ceq%eqno,ceq%status,ceq%next
! status,multi,eqno,next,name,comment,tpval(2),rtn,weight,
! (links to cond,exper), complist(nel),(link to compstoi*(nel*nel))
! old: highcs, (link to phase_varres), mu(nel), xconc,gmind,eqextra,maxiter
! highcs, (link to phase_varres),mu(nel),xconc, gdconv(2),gmind,eqextra,maxiter
!   rsize=4+nwch(24)+nwch(72)+4*nwpr+2+2*noofel+4+3*nwpr
   rsize=4+nwch(24)+nwch(72)+4*nwpr+2+2*noofel+4+5*nwpr
   call wtake(lokeq,rsize,iws)
   if(buperr.ne.0) then
      write(*,*)'3E Error reserving equilibrium record'
      gx%bmperr=4356; goto 1000
   endif
   if(lok1.eq.0) then
! return pointer to first
      lok1=lokeq
   else
! else link from previous
!      write(*,*)'Linking equilibria: ',lok1,loklast,lokeq
      iws(loklast)=lokeq
   endif
   loklast=lokeq
! iws(lokeq) is pointer to next
!   write(*,16)lokeq,ceq%status
16 format('3E equilibrium status word: ',i8,1x,z8)
   iws(lokeq+1)=ceq%status
   iws(lokeq+2)=ceq%multiuse
   iws(lokeq+3)=ceq%eqno
   iws(lokeq+4)=ceq%nexteq
   call storc(lokeq+5,iws,ceq%eqname)
   displace=5+nwch(24)
   call storc(lokeq+displace,iws,ceq%comment)
   displace=displace+nwch(72)
   call storrn(2,iws(lokeq+displace),ceq%tpval)
   call storr(lokeq+displace+2*nwpr,iws,ceq%weight)
   displace=displace+3*nwpr
! svfunres not stored
!---- conditions, write as text and recreated when reading file
   call get_all_conditions(text,0,ceq)
   if(gx%bmperr.ne.0) goto 1000
   kl=index(text,'CRLF')-1
!   write(*,*)'3E cond: ',trim(text),kl
   if(kl.gt.1) then
      call wtake(lok,1+nwch(kl),iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving condition record'
         gx%bmperr=4356; goto 1000
      endif
      call storc(lok+1,iws,text(1:kl))
      iws(lok)=kl
      iws(lokeq+displace)=lok
   else
! no conditions
      iws(lokeq+displace)=0
   endif
!---- save experiments as text
! a bit strange one has to loop incrementing seqz until there is an error ...
   iws(lokeq+displace+1)=0
   seqz=0
   lokcc=lokeq+displace+1
133 continue
   seqz=seqz+1
   j=1
   text=' '
   call get_one_experiment(j,text,seqz,.FALSE.,ceq)
   if(gx%bmperr.ne.0) then
! no or no more experiments
      gx%bmperr=0
   else
! do not save the "current value" after the $
!      write(*,*)'3E save experiment: "',trim(text),'"'
      kl=index(text,'$')-1
      if(kl.le.0) then
         kl=len_trim(text)
      endif
      if(kl.gt.0) then
!         write(*,*)'3E experiment: ',text(1:kl),seqz
         call wtake(lok,2+nwch(kl),iws)
         if(buperr.ne.0) then
            write(*,*)'3E Error reserving experiments record'
            gx%bmperr=4356; goto 1000
         endif
         call storc(lok+2,iws,text(1:kl))
         iws(lok+1)=kl
! create a linear list
         iws(lokcc)=lok
         lokcc=lok
      endif
      goto 133
   endif
!   write(*,*)'3E buperr 1: ',buperr
!---- if components different from elements
   if(btest(globaldata%status,GSNOTELCOMP)) then
      write(*,*)'3E Not implemented saving other components than elements'
      gx%bmperr=4399; goto 1000
!      do i=1,noofel
!         isp=ceq%complist(i)%splink
!         write(lut)isp
!         write(lut)ceq%complist(i)%phlink,ceq%complist(i)%status,&
!              ceq%complist(i)%refstate,ceq%complist(i)%tpref,&
!              ceq%complist(i)%mass
!      enddo
!      do i=1,noofel
!         if(ocv()) write(*,99)'comp.matrix: ',(ceq%invcompstoi(j,i),j=1,noofel)
!      enddo
!99    format(a,7e11.3)
!      do i=1,noofel
!         write(lut)(ceq%compstoi(j,i),j=1,noofel)
!      enddo
   else
! save component records in a linked list NEEDED FOR MANY THINGS
! like reference state etc
      lokcc=lokeq+displace+2
      rsize=5+nwch(16)+1+6*nwpr
      do j=1,noofel
         if(allocated(ceq%complist(j)%endmember)) then
! this component has a user defined reference state
            kl=size(ceq%complist(j)%endmember)
         else
            kl=0
         endif
         call wtake(lok,rsize+kl,iws)
         if(buperr.ne.0) then
            write(*,*)'3E Error reserving varres record 1',j,rsize+kl
            gx%bmperr=4356; goto 1000
         endif
! sequential link
         iws(lokcc)=lok
         lokcc=lok
! data
         iws(lok+1)=ceq%complist(j)%splink
         iws(lok+2)=ceq%complist(j)%phlink
         iws(lok+3)=ceq%complist(j)%status
         call storc(lok+4,iws,ceq%complist(j)%refstate)
         disz=4+nwch(16)
         iws(lok+disz)=kl
         if(kl.gt.0) then
            do mc=1,kl
               iws(lok+disz+mc)=ceq%complist(j)%endmember(mc)
            enddo
            disz=disz+kl+1
         else
            disz=disz+1
         endif
!         write(*,*)'3E refstate 1: ',ceq%complist(j)%tpref
         call storrn(2,iws(lok+disz),ceq%complist(j)%tpref)
         disz=disz+2*nwpr
         call storrn(2,iws(lok+disz),ceq%complist(j)%chempot)
         disz=disz+2*nwpr
         call storr(lok+disz,iws,ceq%complist(j)%mass)
!         write(*,*)'3E saving component mass',lok,disz,j,ceq%complist(j)%mass
         call storr(lok+disz+nwpr,iws,ceq%complist(j)%molat)
!         write(*,*)'3e comprec size: ',lok,lok+disz+nwpr,iws(1)
      enddo
   endif
117 continue
! LINKED LIST of phase_varres records stored from lokeq+lokvares
   lokhighcs=lokeq+displace+3
!   write(*,118)'3E highcs: ',eqnumber,highcs,csfree,lokhighcs
118 format(a,3i5,i10)
   iws(lokhighcs)=highcs
   lokvares=lokhighcs+1
   eqdis=displace+5
!   write(*,*)'3E buperr 2: ',buperr
!   write(*,*)'3E link to first phase_varres in ',lokvares,highcs
!--------------------------------------------------------- below is varres
!---- varres records, one for each composition set of the phases and sometimes
! one for disordered fraction sets ....
! write them in records linked from lokvares as they can be very different
   compset: do j=1,highcs
! loop for all composition sets
      firstvarres=>ceq%phase_varres(j)
      if(.not.allocated(firstvarres%yfr)) then
! if this phase_varres is no longer used this should be unallocated
         call wtake(lok,4,iws)
         if(buperr.ne.0) then
            gx%bmperr=4356; goto 1000
         endif
         write(*,*)'3E unused phase_varres:',j,highcs,lok
! this is the free list
         iws(lok+1)=firstvarres%nextfree
! this should be phlink but set to illegal value
         iws(lok+2)=-1
! this links all phase varres records together
         iws(lokvares)=lok
         lokvares=lok
         cycle compset
      endif
      lokph=firstvarres%phlink
      if(btest(firstvarres%status2,CSDFS)) then
! this phase_varres/parres record belong to disordered fraction_set
! A bit tricky to find the number of sublattices and constituents ....
         lokcs=phlista(lokph)%linktocs(1)
         nsl=ceq%phase_varres(lokcs)%disfra%ndd
         mc=ceq%phase_varres(lokcs)%disfra%tnoofxfr
      else
!         lokcs=0
         nsl=phlista(firstvarres%phlink)%noofsubl
!         mc=phlista(firstvarres%phlink)%tnooffr
! if this phase_varres has been removed this may be unallocated
         if(.not.allocated(firstvarres%yfr)) then
            write(*,*)'3E highcs not updated when removing compset!',j,highcs
! we should update??             iws(lokeq+displace+3)=highcs
            cycle compset
         endif
! wow, firstvarres%yfr is dimensioned to 1000
         mc=phlista(firstvarres%phlink)%tnooffr
!         write(*,*)'3E mc: ',trim(phlista(lokph)%name),mc,size(firstvarres%yfr)
      endif
      if(btest(firstvarres%status2,CSDLNK)) then
! the offset here shold be the place to store the disfra record ...
         offset=6+2*nwch(4)+3*nwpr+mc*(1+2*nwpr)+nsl*nwpr
!         write(*,202)'3E offset 0: ',j,highcs,lokph,nsl,mc,offset
      endif
      mc2=mc*(mc+1)/2
! nextfree,phlink,status2,phstate,phtupx,abnorm(3),prefix*4,suffix*4
! constat(mc),yfr(mc),mmyfr(mc)+2 extra for nsl and mc 
      rsize=6+2*nwch(4)+3*nwpr+mc+2*mc*nwpr
! sites(nsl),disfralink,amfu,netcharge,dgm and link to ionliq dpqdy record!!
! also added qcbonds!!
!      rsize=rsize+nsl*nwpr+1+4*nwpr+2
      rsize=rsize+nsl*nwpr+1+5*nwpr+2
! results g, dg, d2g some exra space
      rsize=rsize+6*nwpr+3*mc*nwpr+mc2*nwpr+5+2
      qsize=rsize
      call wtake(lok,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving varres record 2',j,rsize,nsl,mc
         gx%bmperr=4356; goto 1000
      endif
      iws1=iws(1)
!      lokph=firstvarres%phlink
!      write(*,107)'3E saving: ',j,lok,rsize,mc,nsl,trim(phlista(lokph)%name)
!      write(*,107)'3E saving: ',j,phasetuple(j)%ixphase,nsl,0,0
107   format(a,i3,2i10,i4,i3,2x,a)
! link from lokvares and use iws(lok) to link to next
      iws(lokvares)=lok
      lokvares=lok
! data
      iws(lok+1)=firstvarres%nextfree
      iws(lok+2)=firstvarres%phlink
      iws(lok+3)=firstvarres%status2
      iws(lok+4)=firstvarres%phstate
      iws(lok+5)=firstvarres%phtupx
      iws(lok+6)=nsl
      iws(lok+7)=mc
      call storc(lok+8,iws,firstvarres%prefix)
      displace=8+nwch(4)
      call storc(lok+displace,iws,firstvarres%suffix)
      displace=displace+nwch(4)
      call storrn(3,iws(lok+displace),firstvarres%abnorm)
      displace=displace+3*nwpr
!      write(*,*)'3E sizes:',allocated(firstvarres%constat),&
!           size(firstvarres%constat),size(firstvarres%yfr),mc
      do i=1,mc
         iws(lok+displace+i-1)=firstvarres%constat(i)
      enddo
      displace=displace+mc
      call storrn(mc,iws(lok+displace),firstvarres%yfr)
      displace=displace+mc*nwpr
! mmyfr is just reals ... do not bother (although space for double reserved)
!      write(lut)(firstvarres%mmyfr(i),i=1,mc)
      displace=displace+mc*nwpr
!      write(*,*)'3E sites:',lok,displace,lok+displace
      call storrn(nsl,iws(lok+displace),firstvarres%sites)
      displace=displace+nsl*nwpr
! do not save the cmuval array
! dsitesdy is interesting only for ionic liquids
!      if(btest(phlista(lokph)%status1,PHIONLIQ)) then
!         call wtake(mc+mc2,iws
!         call storrn(mc,iws(lok+displace),firstvarres%dpqdy)
!         displace=displace+mc
!         call storrn(mc2,iws(lok+displace),firstvarres%d2pqdvay)
!         displace=displace+mc2
!         write(*,*)'3E odd:   ',lok,displace
!      else
!         iws(
!      endif
      fsrec: if(btest(firstvarres%status2,CSDLNK)) then
! we need a record for a disordered fraction_set record
! latd,ndd,tnoofxfr,tnoofyfr,varreslink,totdis, id*1, dsites(nsl), 
! nooffr(mc), splink(mc), y2x(mc), dxidyj(mc),fsites
         fslink=>firstvarres%disfra
         nsl=fslink%ndd
! dmc because we store G and dG/dy later for original mc
         dmc=fslink%tnoofxfr
         rsize=8+nwch(1)+nsl+dmc+1+mc*(1+nwpr)+nsl*nwpr+nwpr
         call wtake(lokdis,rsize,iws)
         if(buperr.ne.0) then
            write(*,*)'3E Error reserving disordered varres record',rsize
            gx%bmperr=4356; goto 1000
         endif
!         write(*,202)'3E disfracset 1: ',j,lok,displace,lokdis,nsl,dmc
202      format(a,10i6)
! set link from varres record
         iws(lok+displace)=lokdis
! store data
         iws(lokdis)=fslink%latd
         iws(lokdis+1)=fslink%ndd
         iws(lokdis+2)=fslink%tnoofxfr
         iws(lokdis+3)=fslink%tnoofyfr
         iws(lokdis+4)=fslink%totdis
         iws(lokdis+5)=fslink%varreslink
         call storc(lokdis+6,iws,fslink%id)
!         write(*,202)'3E disfracset 2: ',j,iws(lokdis+1),iws(lokdis+2),&
!              iws(lokdis+5)
! set disz to one less as i starts from 1
         disz=6+nwch(1)
! number of constituents in each sublattice
         do i=1,nsl
            iws(lokdis+disz+i)=fslink%nooffr(i)
         enddo
         disz=disz+nsl
!         write(*,*)'3E disfra 1: ',lokdis,disz
! species index for all constituents
         do i=1,dmc
            iws(lokdis+disz+i)=fslink%splink(i)
         enddo
         disz=disz+dmc+1
         iws(lokdis+disz)=mc
         disz=disz
! NOTE y2x and dxidy1 has dimension mc!!
!         write(*,*)'3E disfra 2: ',lokdis,disz,dmc,mc,size(fslink%y2x)
! This has to do with the fractions that should be added together
         do i=1,mc
            iws(lokdis+disz+i)=fslink%y2x(i)
         enddo
         disz=disz+mc+1
!         write(*,*)'3E disfra 3: ',lokdis,disz
! number of sites in each sublattice
         call storrn(nsl,iws(lokdis+disz),fslink%dsites)
         disz=disz+nsl*nwpr
!         write(*,*)'3E disfra 4: ',lokdis,disz,dmc,mc,size(fslink%y2x)
! converting ordered fractions to disordered fractions
         call storrn(mc,iws(lokdis+disz),fslink%dxidyj)
! formula unit factor
         disz=disz+mc*nwpr
!         write(*,*)'3E disfra 5: ',lokdis,disz
         call storr(lokdis+disz,iws,fslink%fsites)
!         write(*,*)'3E disfra: ',lokdis+disz+nwpr,iws(1)
      else
! mark no link to disordered record
!         write(*,*)'3E no disorderd record',lok+displace,iws(1),iws(2)
         iws(lok+displace)=0
      endif fsrec
!      write(*,*)'3E buperr 7: ',buperr,lok,displace
!------------------------------------- end of disorderd record
! save some results stored in the phase_varres record
      displace=displace+1
      call storr(lok+displace,iws,firstvarres%amfu)
      call storr(lok+displace+nwpr,iws,firstvarres%netcharge)
      call storr(lok+displace+2*nwpr,iws,firstvarres%dgm)
! record size increased to save qcbonds ... and increment of displace below
      call storr(lok+displace+3*nwpr,iws,firstvarres%qcbonds)
! Maybe firstvarres%nprop is not always initiated??
! it seems that additional compsets have an arbitrary value ...
      if(firstvarres%nprop.ne.20) then
         iws(lok+displace+4*nwpr)=20
      else
         iws(lok+displace+4*nwpr)=firstvarres%nprop
      endif
!      write(*,303)'3E saving nprop: ',lok,displace+3*nwpr,lok+displace+3*nwpr,&
!           iws(lok+displace+3*nwpr),trim(phlista(lokph)%name)
303   format(a,4i8,2x,a)
      displace=displace+4*nwpr+1
! only save G and derivatives
      do i=1,6
         call storr(lok+displace+nwpr*(i-1),iws,firstvarres%gval(i,1))
      enddo
! problem here with SELECT_ELEMENT_REFERENCE phase ...
!      write(*,304)'3E bug: ',trim(phlista(lokph)%name),mc,&
!           size(firstvarres%dgval)
!304   format(a,a,5i5)
! in the ENTER_EQUILIBRIUM the incorrect size of dgval was allocated !!! fixed
      displace=displace+6*nwpr
      do i=1,3
         do k=1,mc
!            write(*,*)'indices: ',i,k
            call storr(lok+displace,iws,firstvarres%dgval(i,k,1))
            displace=displace+nwpr
         enddo
      enddo
      do i=1,mc2
         call storr(lok+displace+nwpr*(i-1),iws,firstvarres%d2gval(i,1))
      enddo
!      write(*,*)'3E last values used ',j,lok+displace+mc2*nwpr,lok+qsize,iws1
   enddo compset
!----------------------------------------
! we must set csfree to highcs+1
! as new composition sets will use that as the free list pointer
   csfree=highcs+1
!-----------------------------------------
! mu(nel), xconc,gmind,eqextra,maxiter
! MODIFIED: mu(nel), xconc, gdconv(2), gmind,eqextra,maxiter
   iws(lokeq+eqdis)=ceq%maxiter
   call storrn(noofel,iws(lokeq+eqdis+1),ceq%cmuval)
   eqdis=eqdis+1+noofel*nwpr
   call storr(lokeq+eqdis,iws,ceq%xconv)
   call storr(lokeq+eqdis+nwpr,iws,ceq%gdconv(1))
   call storr(lokeq+eqdis+2*nwpr,iws,ceq%gdconv(2))
   call storr(lokeq+eqdis+3*nwpr,iws,ceq%gmindif)
! last use of lokeq !!
!   write(*,*)'3E NOT saving the character eqextra!'
!   call storc(lokeq+displace+2*nwpr,iws,ceq%eqextra)
!   write(*,*)'3E check rsize: ',rsize,eqdis+2*nwpr
!>>>>>> 64: savesysmat
!   write(*,*)'3E not saving sysmat?',ceq%sysmatdim,ceq%nfixmu,ceq%nfixph
! NOTE:: ceq%sysmatdim negative, not initiallized??
! NOTE:: phasetuples not saved !!!
!   write(lut)ceq%sysmatdim,ceq%nfixmu,ceq%nfixph
!   if(ceq%nfixmu.gt.0) write(lut)(ceq%fixmu(kl),kl=1,ceq%nfixmu)
!   if(ceq%nfixph.gt.0) write(lut)&
!        (ceq%fixph(1,kl),ceq%fixph(2,kl),kl=1,ceq%nfixph)
!   if(ceq%sysmatdim.gt.0) then
!      do mc=1,ceq%sysmatdim
!         write(lut)(ceq%savesysmat(mc,kl),kl=1,ceq%sysmatdim)
!      enddo
!   endif
! loop for all entered equilibria
   goto 17
!----------
1000 continue
   return
 end subroutine saveequil

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine svfunsave
!\begin{verbatim}
 subroutine svfunsave(loksvf,iws,ceq)
! saves all state variable functions as texts in iws
   implicit none
   integer iws(*),loksvf
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   character text*512,symbols(20)*32,afterdot*32
   integer ip,ipos,istv,js,jt,kl,ks,lrot,rsize,lok
   type(gtp_state_variable), target :: svr2
   type(gtp_state_variable), pointer :: svrrec
   rsize=nsvfun+5
   call wtake(loksvf,rsize,iws)
   if(buperr.ne.0) then
      write(*,*)'3E Error reserving state variable function record',rsize,iws(1)
      gx%bmperr=4356; goto 1000
   endif
   iws(loksvf)=nsvfun
   iws(loksvf+1)=3
! do not save the first three, R, RT and T_C
   symbols=' '
   write(*,*)'3E saving ',nsvfun,' symbols as texts'
   do lrot=4,nsvfun
      ipos=1
      text=' '
      call list_svfun(text,ipos,lrot,ceq)
      rsize=1+nwch(ipos)
      call wtake(lok,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving state variable func record',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
!      write(*,*)'3E storing svfun: ',text(1:ipos)
      iws(lok)=ipos
! NOTE position 1-7 is equilibrium number and status
      call storc(lok+1,iws,text(1:ipos))
      iws(loksvf+lrot)=lok
   enddo
1000 continue
   return
 end subroutine svfunsave

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine bibliosave
!\begin{verbatim}
 subroutine bibliosave(bibhead,iws)
! saves references on a file
   implicit none
   integer bibhead,iws(*)
!\end{verbatim} %+
   character longline*2048
   integer ir,jp,ll,nl,lok,rsize
!>>>>> 40:
!   write(*,*)'3E Saving reference version and number of:',&
!        gtp_biblioref_version,reffree-1
   rsize=3+reffree-1
   call wtake(bibhead,rsize,iws)
   if(buperr.ne.0) then
      write(*,*)'3E Error reserving biblographiic record',rsize,iws(1)
      gx%bmperr=4356; goto 1000
   endif
   iws(bibhead)=reffree-1
   do ir=1,reffree-1
! a bibliographic reference contains 16 character identifier and a variable
! characters text.  Concatinate that into a single text and save one
! reference in each record linked from bibhead
      longline=bibrefs(ir)%reference
      jp=17
! This require Fortran 2003/2008 standard, not available in GNU Fortran 4.8 
!      longline(17:)=bibrefs(ir)%nyrefspec
      ll=bibrefs(ir)%wprefspec(1)
      call loadc(2,bibrefs(ir)%wprefspec,longline(17:ll+16))
!      nl=size(bibrefs(ir)%refspec)
!      do ll=1,nl
!         longline(jp:)=bibrefs(ir)%refspec(ll)
!         jp=jp+64
!      enddo
      jp=len_trim(longline)
      rsize=1+nwch(jp)
      call wtake(lok,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving biblographiic record',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
      iws(lok)=jp
      call storc(lok+1,iws,longline(1:jp))
      iws(bibhead+ir)=lok
   enddo
1000 continue
   return
 end subroutine bibliosave

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine saveash
!\begin{verbatim}
 subroutine saveash(lok,iws)
! saving assessment records
   integer lok,iws(*)
!\end{verbatim} %+
   integer lok1,lok2,last,rsize,i1,i2,disp
   type(gtp_assessmenthead), pointer :: assrec
!   type(gtp_equilibrium_data), pointer :: ceq
!
   assrec=>firstash%nextash
   if(.not.allocated(assrec%eqlista)) then
      write(kou,*)'3E No experimental equilibrium range set'
!      iws(lok)=0
!      goto 1000
   endif
20 continue
! next, status, varcoef, first, and 8 allocatable arrays
!   rsize=4+2*nwch(64)+10
! added one location for pointer to RSD values
!   rsize=4+2*nwch(64)+11
   rsize=5+2*nwch(64)+11
   write(*,*)'3E allocating assessment head record',rsize
   call wtake(lok1,rsize,iws)
   if(buperr.ne.0) then
      write(*,*)'3E Error reserving assessment record',rsize,iws(1)
      gx%bmperr=4356; goto 1000
   endif
   if(iws(lok).eq.0) then
      iws(lok)=lok1
      last=lok1
   else
      iws(last)=lok1
      last=lok1
   endif
   iws(lok1+1)=assrec%status
   iws(lok1+2)=assrec%varcoef
   iws(lok1+3)=assrec%firstexpeq
   iws(lok1+4)=assrec%firstexpeq
!   call storc(lok1+4,iws,assrec%general)
   call storc(lok1+5,iws,assrec%general)
   disp=5+nwch(64)
   call storc(lok1+disp,iws,assrec%special)
   disp=disp+nwch(64)
! eqlista CAN BE EMPTY!
   if(allocated(assrec%eqlista)) then
      i1=size(assrec%eqlista)
      rsize=1+i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
!   write(*,*)'3E in saveash 1:',lok,lok1,lok2,i1
      iws(lok2)=i1
      if(i1.gt.0) then
! Hm assrec%eqlista(i2)%p1 is a pointer to an element in the global eqlists
!   ceq=>assrec%eqlista(1)%p1
         do i2=1,i1
            iws(lok2+i2)=assrec%eqlista(i2)%p1%eqno
         enddo
      endif
   else
! mark that no experimental records
      lok2=0
   endif
   iws(lok1+disp+1)=lok2
! coeffvalues
   if(allocated(assrec%coeffvalues)) then
      i1=size(assrec%coeffvalues)
      rsize=1+nwpr*i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
!      write(*,*)'3E in saveash 2:',lok2,i1,rsize
      iws(lok2)=i1
      call storrn(i1,iws(lok2+1),assrec%coeffvalues)
      iws(lok1+disp+2)=lok2
! relative standard deviation
      i1=size(assrec%coeffvalues)
      rsize=1+nwpr*i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
!      write(*,*)'3E in saveash 2 RSD:',lok2,i1,rsize
      iws(lok2)=i1
      call storrn(i1,iws(lok2+1),assrec%coeffrsd)
      iws(lok1+disp+3)=lok2
! coeffscale
      i1=size(assrec%coeffscale)
      rsize=1+nwpr*i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
      iws(lok2)=i1
!      write(*,*)'3E in saveash 3:',lok2,i1
      call storrn(i1,iws(lok2+1),assrec%coeffscale)
!      iws(lok1+disp+3)=lok2
      iws(lok1+disp+4)=lok2
! coeffstart
      i1=size(assrec%coeffstart)
      rsize=1+nwpr*i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
      iws(lok2)=i1
!      write(*,*)'3E in saveash 4:',lok2,i1
      call storrn(i1,iws(lok2+1),assrec%coeffstart)
!      iws(lok1+disp+4)=lok2
      iws(lok1+disp+5)=lok2
! coeffmin
      i1=size(assrec%coeffmin)
      rsize=1+nwpr*i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
      iws(lok2)=i1
!      write(*,*)'3E in saveash 5:',lok2,i1
      call storrn(i1,iws(lok2+1),assrec%coeffmin)
!      iws(lok1+disp+5)=lok2
      iws(lok1+disp+6)=lok2
! coeffmax
      i1=size(assrec%coeffmax)
      rsize=1+nwpr*i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
      iws(lok2)=i1
      call storrn(i1,iws(lok2+1),assrec%coeffmax)
!      iws(lok1+disp+6)=lok2
      iws(lok1+disp+7)=lok2
! coeffindices
      i1=size(assrec%coeffindex)
      rsize=1+i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
!      write(*,*)'3E in saveash 6:',lok2,i1
      iws(lok2)=i1
      do i2=1,i1
         iws(lok2+i2)=assrec%coeffindex(i2-1)
      enddo
!      iws(lok1+disp+7)=lok2
      iws(lok1+disp+8)=lok2
! coeffstate
      i1=size(assrec%coeffstate)
      rsize=1+i1
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
      iws(lok2)=i1
      do i2=1,i1
         iws(lok2+i2)=assrec%coeffstate(i2-1)
      enddo
!      iws(lok1+disp+8)=lok2
      iws(lok1+disp+9)=lok2
   else
! pointers are zero
      write(*,*)'3E no coefficients allocated'
   endif
! maybe work array should not be saved?
   if(allocated(assrec%wopt)) then
      i1=size(assrec%wopt)
      rsize=1+nwpr*i1
      write(*,*)'3E saving assessment record: (assrec%wopt)',lok1,rsize
      call wtake(lok2,rsize,iws)
      if(buperr.ne.0) then
         write(*,*)'3E Error reserving assessment record array',rsize,iws(1)
         gx%bmperr=4356; goto 1000
      endif
      iws(lok2)=i1
      call storrn(i1,iws(lok2+1),assrec%wopt)
!      iws(lok1+disp+9)=lok2
      iws(lok1+disp+10)=lok2
   else
      write(*,*)'3E no work array (assrec%wopt) allocated'
!      iws(lok1+disp+9)=0
      iws(lok1+disp+10)=0
   endif
! check if there are several assessment records
   if(.not.associated(assrec,firstash)) then
      assrec=>assrec%nextash
      write(*,*)'3E more than one assessment records'
      goto 20
   endif
1000 continue
   return
 end subroutine saveash

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function ceqrecsize
!\begin{verbatim}
 integer function ceqrecsize()
! calculates the number of words needed to save an equilibrium record
!\end{verbatim}
   integer rsize,jj,seqz,kl,dmc,mc,mc2,nsl
   type(gtp_equilibrium_data), pointer :: ceq
   type(gtp_phase_varres), pointer :: firstvarres
   TYPE(gtp_fraction_set), pointer :: fslink
   character text*512
!
   write(*,*)'ceqrecsize not implemented',highcs
   rsize=0
   goto 1000
   ceq=>firsteq
   rsize=4+nwch(24)+nwch(72)+4*nwpr+2+2*noofel+4+5*nwpr
   text=' '
   call get_all_conditions(text,0,ceq)
   rsize=rsize+nwch(index(text,'CRLF'))
100 continue
   text=' '
   call get_one_experiment(jj,text,seqz,.FALSE.,ceq)
   if(gx%bmperr.ne.0) then
      kl=index(text,'$')-1
      if(kl.le.0) then
         kl=len_trim(text)
      endif
      rsize=rsize+2+nwch(kl)
      goto 100
   endif
   gx%bmperr=0
! ignore if a component has a defined reference state ...
   rsize=rsize+(5+nwch(16)+1+6*nwpr)*noofel
   do jj=1,highcs
! loop for phase_varres records ..
      firstvarres=>ceq%phase_varres(jj)
      if(.not.allocated(firstvarres%yfr)) then
         rsize=rsize+4
      else
         rsize=rsize+6+2*nwch(4)+3*nwpr+mc+2*mc*nwpr
         rsize=rsize+6*nwpr+3*mc*nwpr+mc2*nwpr+5+2
         if(btest(firstvarres%status2,CSDLNK)) then
! there is a disordered fraction set ...
            fslink=>firstvarres%disfra
            nsl=fslink%ndd
            rsize=8+nwch(1)+nsl+dmc+1+mc*(1+nwpr)+nsl*nwpr+nwpr
         endif
      endif
   enddo
1000 continue
   ceqrecsize=rsize
   return
 end function ceqrecsize

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine gtpread
!\begin{verbatim}
 subroutine gtpread(filename,str)
! read unformatted all data in the following order
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! references
! first equilibrium record with conditions, componenets, phase_varres etc
! state variable functions
! equilibrium record(s) with conditions, componenets, phase_varres, experim etc
! CCI changed to use iso_fortran_env to find file unit number for C++
   use :: iso_fortran_env
! CCI
   implicit none
   character*(*) filename,str
!\end{verbatim} %+
!
   character id*40,version*8,comment*72
   integer i,i1,i2,i3,isp,jph,kontroll,nel,ivers,lin,last,lok,displace,jfun
   integer nspx,saverr
   integer, allocatable :: iws(:)
! CCI
   logical is_op
! CCI
!   type(gtp_equilibrium_data), pointer :: ceq
10  format(i8)
   if(index(filename,'.').eq.0) then
      filename(len_trim(filename)+1:)='.OCU'
   endif
!CCI The previous commented lines are removed by the following lines 
!CCI that enable to find the first available logical unit. 
!CCI Such an approach can generalized in order to enable the 
!CCI opening file by several threads in the same time. 
!CCI To do this, the following lines should in a dedicated subroutine.
!CCI   lin=21
   lunit: do lin=8,99  
      inquire(lin,opened=is_op)
      if(.not.is_op) exit lunit
   enddo lunit 
   if( lin.eq.100 ) then
      write(*,*)'3E Error, no logical unit available for opening file: ',&
           trim(filename)
      goto 1000
   endif
! CCI end change   
   open(lin,file=filename,access='sequential',status='old',&
        form='unformatted',iostat=gx%bmperr,err=1100)
!   write(*,*)'3E opening file: ',trim(filename),' for unformatted read'
!
   read(lin)id,version,comment,noofel,noofsp,noofph,nooftuples,last
!**********************************************************
! IMPORTANT savefile
! is a variable in gtp3.F90
! which MUST BE CHANGED whenever there is a change in the unformatted
! file layout
!***********************************************************
   if(version.ne.savefile) then
      write(*,11)id,version,savefile
11     format('File not same version as program: ',A/a,' : ',a)
      gx%bmperr=4299; goto 900
   endif
   write(*,12)id,version,trim(comment)
12 format(/'Read unformatted file: ',a,a/'Generated: ',a/)
   str=comment
!   write(*,*)'3E numbers: ',noofel,noofsp,noofph,nooftuples,last
!-------
   allocate(iws(last))
   read(lin)(iws(i),i=1,last)
   close(lin)
!------------------------
!>>>>> 2: elementlist, follow link from iws(3)
   if(iws(4).ne.gtp_element_version) then
      write(*,*)'3E Element data structure not same:',iws(4),gtp_element_version
      gx%bmperr=4355; goto 1000
   endif
   nel=0
   last=iws(3)
   do while(last.gt.0)
      nel=nel+1
      lok=last
      call loadc(lok+1,iws,ellista(nel)%symbol)
      call loadc(lok+2,iws,ellista(nel)%name)
      displace=3+nwch(12)
      call loadc(lok+displace,iws,ellista(nel)%ref_state)
      displace=displace+nwch(24)
      call loadr(lok+displace,iws,ellista(nel)%mass)
      call loadr(lok+displace+nwpr,iws,ellista(nel)%h298_h0)
      call loadr(lok+displace+2*nwpr,iws,ellista(nel)%s298)
      displace=displace+3*nwpr
      ellista(nel)%splink=iws(lok+displace)
      ellista(nel)%status=iws(lok+displace+1)
      ellista(nel)%alphaindex=iws(lok+displace+2)
      ellista(nel)%refstatesymbol=iws(lok+displace+3)
!      write(*,17)ellista(nel)%symbol,ellista(nel)%name,ellista(nel)%ref_state,&
!           ellista(nel)%mass,ellista(nel)%h298_h0,ellista(nel)%s298,&
!           ellista(nel)%splink,ellista(nel)%status,ellista(nel)%alphaindex,&
!           ellista(nel)%refstatesymbol
17    format('3E ',a2,2x,a12,2x,a24,2x,3(1pe12.4),4i5)
! do not forget the element array!
      elements(ellista(nel)%alphaindex)=nel
!
      last=iws(last)
!      write(*,*)'3E elloop: ',nel,lok,last,iws(1)
   enddo
   if(nel.ne.noofel) then
      write(*,*)'3E Number of elements wrong: ',nel,noofel
   endif
!   write(*,*)'3E Now the species!!'
!-------
!>>>>> 3: specieslist NOTE ADDES SPEXTRA
   if(iws(6).ne.gtp_species_version) then
      write(*,*)'3E Species version wrong: ',iws(5),gtp_species_version
      gx%bmperr=4355; goto 1000
   endif
   last=iws(5)
! VA is entered automatically at first index in splista when reinitiating 
! so keep that.  We just skip the first species in iws and extract
! its alphaindex
   splista(1)%alphaindex=iws(last+2+nwch(24)+2*nwpr+2)
   species(splista(1)%alphaindex)=1
! skip the first species (this is VA)
   last=iws(last)
   isp=1
   do while(last.gt.0)
      isp=isp+1
!      write(*,*)'3E loop: ',last,isp,splista(isp-1)%symbol
      call loadc(last+1,iws,splista(isp)%symbol)
      displace=2+nwch(24)
      call loadr(last+displace,iws,splista(isp)%mass)
      call loadr(last+displace+nwpr,iws,splista(isp)%charge)
      splista(isp)%noofel=iws(last+displace+2*nwpr)
      splista(isp)%status=iws(last+displace+2*nwpr+1)
      splista(isp)%alphaindex=iws(last+displace+2*nwpr+2)
! new spextra array
      nspx=iws(last+displace+2*nwpr+3)
!      if(nspx.ne.0) write(*,*)'3E nspx value: ',nspx
      allocate(splista(isp)%ellinks(splista(isp)%noofel))
      allocate(splista(isp)%stoichiometry(splista(isp)%noofel))
      displace=displace+2*nwpr+3
      do i=1,splista(isp)%noofel
         splista(isp)%ellinks(i)=iws(last+displace+i)
      enddo
      displace=displace+splista(isp)%noofel+1
!      write(*,*)'3E displace load: ',last,displace
      call loadrn(splista(isp)%noofel,&
           iws(last+displace),splista(isp)%stoichiometry)
      species(splista(isp)%alphaindex)=isp
! handle spextra values if any
      if(nspx.gt.0) then
!         write(*,*)'We have nonzero nxsp: ',nspx
         allocate(splista(isp)%spextra(nspx))
         displace=displace+splista(isp)%noofel*nwpr
         call loadrn(nspx,iws(last+displace),splista(isp)%spextra)
      endif
! next species
      last=iws(last)
   enddo
   if(isp.ne.noofsp) then
      write(*,*)'3E wrong number of species: ',isp,noofsp
      gx%bmperr=4399; goto 1000
   endif
!---------- component record
! read inside the equilibrium record   
!---------- tpfuns
!>>>>> 20.. inside tpfunread, skip functions already read
   last=7
   if(iws(8).ne.tpfun_expression_version) then
      write(*,*)'3E tpfun_expression_version not same',iws(8),&
           tpfun_expression_version
      gx%bmperr=4355; goto 1000
   endif
   isp=iws(last)
   i3=iws(isp)
!   write(*,*)'3E tpfuns',iws(7),iws(8),i3
   if(isp.gt.0) then
! skip first 2 (R and RTLNP)
      do i=3,i3-1
         call read0tpfun(iws(isp+i),iws,i)
         if(gx%bmperr.ne.0) then
            write(*,*)'3E Error reading TP function: ',gx%bmperr
            goto 1000
         endif
      enddo
   endif
! we cannot update freetpfun before all functions are entered ....
   freetpfun=i3
! hopefully the TP functions will keep the same index ... so for parameters
! one just store the index!
!-------
!>>>>> 5: phaselist, starting from 0, the reference phase
! zero number of phases etc
   noofph=0
   nooftuples=0
   noofem=0
   noofint=0
   noofprop=0
! link to phaselist is in 9 (+10, 11, 12, 13)
   call readphases(noofph,iws)
   if(gx%bmperr.ne.0) goto 1000
!-----------
! restore phase tuples
!   write(*,*)'3E Reading phase tuples',iws(14),noofph
   lok=iws(14)
   if(lok.gt.0) then
      if(iws(15).ne.gtp_phasetuple_version) then
         write(*,*)'3E wrong phasetuple version',gtp_phasetuple_version,iws(15)
         gx%bmperr=4355; goto 1000
      endif
      nooftuples=iws(lok)
      do i=1,nooftuples
         phasetuple(i)%lokph=iws(lok+5*i-4)
         phasetuple(i)%compset=iws(lok+5*i-3)
         phasetuple(i)%ixphase=iws(lok+5*i-2)
         phasetuple(i)%lokvares=iws(lok+5*i-1)
         phasetuple(i)%nextcs=iws(lok+5*i)
      enddo
   endif
! restore the phases lista using phase tuples!
   do jph=1,noofph
      phases(jph)=phasetuple(jph)%lokph
   enddo
!-------------------------------
! the global status word in 20-21
   lok=iws(20)
   globaldata%status=iws(lok+1)
! BUGFIX and extended
   call loadc(lok+2,iws,globaldata%name)
   displace=2+nwch(24)
   call loadr(lok+displace,iws,globaldata%rgas)
   call loadr(lok+displace+nwpr,iws,globaldata%rgasuser)
   call loadr(lok+displace+2*nwpr,iws,globaldata%pnorm)
   displace=displace+3*nwpr
   do i=0,9
      globaldata%sysparam(i+1)=iws(lok+displace+i)
   enddo
! this was used to test record read correctly
!   if(globaldata%sysparam(1).ne.987 .or. &
!        globaldata%sysparam(10).ne.17) then
!      write(*,'(a,10i4)')'3E error globaldata: ',globaldata%sysparam
!   endif
   displace=displace+10
   call loadrn(5,iws(lok+displace),globaldata%sysreal)
!   if(abs(globaldata%sysreal(5)-12345678.9D0).gt.1.0D-12) then
! this was used to test the storing
!      write(*,'(a,5(1pe12.4))')'3E error 2: ',globaldata%sysreal
!   endif
!   write(*,*)'3E name: "',globaldata%name,'"'
! partly unfinished below
!---------- bibliographic references
!>>>>> 40.. inside refread
!   write(*,*)'3E reading bibliography'
   if(iws(23).ne.gtp_biblioref_version) then
      write(*,*)'3E Bibliography version wrong ',iws(23),gtp_biblioref_version
   else
      call biblioread(iws(22),iws)
      if(gx%bmperr.ne.0) goto 1000
   endif
!---------- enter the first equilibrium record without experiments!!
   if(iws(17).ne.gtp_equilibrium_data_version) then
      write(*,*)'3E Wrong equilibrium data version',&
           iws(17),gtp_equilibrium_data_version
      gx%bmperr=4355; goto 1000
   elseif(iws(18).ne.gtp_component_version) then
      write(*,*)'3E Wrong component version',iws(18),gtp_component_version
      gx%bmperr=4355; goto 1000
   elseif(iws(19).ne.gtp_phase_varres_version) then
      write(*,*)'3E Wrong phase varres version',iws(19),gtp_phase_varres_version
      gx%bmperr=4355; goto 1000
   endif
! link to first saved in equilibrium in iws(16). firsteq is eqlista(1)
   i=iws(16)
!   call readequil(i,iws,1,firsteq)
   call readequil(i,iws,1)
   if(gx%bmperr.ne.0) goto 1000
!---------- state variable functions must be present when reading experiments
! and the equilibria must
!>>>>> 30... inside svfunread
!   write(*,*)'3E reading state variable functions',iws(24)
   if(iws(25).eq.gtp_putfun_lista_version) then
      call svfunread(iws(24),iws)
      if(gx%bmperr.ne.0) goto 1000
   else
      write(*,*)'3E state variable function version error',iws(25),&
           gtp_putfun_lista_version
   endif
! we cannot list svfun as we have no ceq ...
!   call list_all_svfun(kou,ceq)
!   call list_some_svfun(kou)
!   write(*,*)'3E Now reading equilibria',iws(16)
!--------------------------------------------------------------------
! read remaining equilibria which may contain experiments
! link to first saved in equilibrium in iws(16)
   i=iws(16)
   i3=2
   call readequil(i,iws,-1)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3E read all equilibria'
!-------------------------------------------------------------------
! read assessment head recods
   if(iws(27).ne.gtp_assessment_version) then
      write(*,*)'3E wrong assemmenst record version',iws(27),&
           gtp_assessment_version
   endif
   lok=26
   call readash(lok,iws)
   if(gx%bmperr.ne.0) goto 1000
   write(*,*)'3E Read assessment record'
!------ read all ??
800 continue
! emergency exit
900 continue
! file already closed above
!   close(lin)
!
1000 continue
!CCI free the iws memory (should be done automatically?)
   if(allocated(iws)) deallocate(iws)
   if(gx%bmperr.eq.4355) then
      write(*,*)'3E *** ERROR unformatted file wrong version'
      saverr=gx%bmperr; gx%bmperr=0
! clear errr code to reinitiate ... it may not work as datastructure bad
      call new_gtp
!      if(gx%bmperr.ne.0) then
!         write(*,*)'Failed to reinitiate',gx%bmperr
!      endif
      stop 'Cannot restore data structures'
      gx%bmperr=saverr
   endif
   return
! error opening files
1100 continue
   write(*,1110)gx%bmperr,trim(filename)
1110 format('I/O error: ',i5,', opening file; ',a)
   goto 1000
 end subroutine gtpread

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine readphases
!\begin{verbatim}
 subroutine readphases(kkk,iws)
! read data for phlista and all endmembers etc
! works for test case without disordered fraction test
   implicit none
   integer kkk,iws(*)
!\end{verbatim} %+
   integer firstendmem,i,i1,i2,i3,jph,level,nem,noi,nop,nox,nup,nsl,mult,lin
   integer lok,displace,totcon,phreclink,lokem,lokint
   type(gtp_endmember), allocatable, target :: nyemrec
   type(gtp_endmember), pointer :: emrec,lem
   type(gtp_interaction), pointer :: intrec
   type(gtp_property), pointer :: proprec
   logical ifbcc
   type saveint
      type(gtp_interaction), pointer :: p1
      integer noi
   end type saveint
   type(saveint), dimension(:), pointer :: stack
   type(gtp_phase_add), pointer :: addlink,nyaddlink
!
   allocate(stack(5))
!   write(*,*)'3E in readphase:',iws(9),iws(10),iws(11),iws(12),iws(13),&
!        iws(7),iws(8)
! as the phlista record contain pointers each item must be read separately
! the phaes are stored sequentially from iws(9)
   lok=9
   if(iws(lok+1).ne.gtp_phase_version) then
      write(*,*)'3E phase version not the same ',iws(lok+1),gtp_phase_version
      gx%bmperr=4355; goto 1000
   elseif(iws(lok+2).ne.gtp_endmember_version) then
      write(*,*)'3E endmember not the same ',iws(lok+2),gtp_endmember_version
      gx%bmperr=4355; goto 1000
   elseif(iws(lok+3).ne.gtp_interaction_version) then
      write(*,*)'3E interaction not the same ',iws(lok+3),&
           gtp_interaction_version
      gx%bmperr=4355; goto 1000
   elseif(iws(lok+4).ne.gtp_property_version) then
      write(*,*)'3E property version not the same ',iws(lok+4),&
           gtp_property_version
      gx%bmperr=4355; goto 1000
   endif
! first phase (number 0) is SER phase
   jph=-1
   lok=iws(lok)
   bigloop: do while(lok.gt.0)
      jph=jph+1
      call loadc(lok+1,iws,phlista(jph)%name)
      displace=1+nwch(24)
      call loadc(lok+1,iws,phlista(jph)%models)
      displace=displace+nwch(72)
      call loadc(lok+1,iws,phlista(jph)%phletter)
      displace=displace+1
      phlista(jph)%status1=iws(lok+displace)
      phlista(jph)%alphaindex=iws(lok+displace+1)
      phlista(jph)%noofcs=iws(lok+displace+2)
      phlista(jph)%nooffs=iws(lok+displace+3)
!   read(lin)jph,phlista(jph)%name,&
!        phlista(jph)%models,phlista(jph)%phletter,phlista(jph)%status1,&
!        phlista(jph)%alphaindex,phlista(jph)%noofcs,phlista(jph)%nooffs
!>>>>> 6: sublattice info
!   read(lin)phlista(jph)%noofsubl,phlista(jph)%linktocs,phlista(jph)%tnooffr
      nsl=iws(lok+displace+4)
      phlista(jph)%noofsubl=nsl
      totcon=iws(lok+displace+5)
      phlista(jph)%tnooffr=totcon
      allocate(phlista(jph)%nooffr(nsl))
      allocate(phlista(jph)%constitlist(phlista(jph)%tnooffr))
!   read(lin)(phlista(jph)%nooffr(i),i=1,nsl),&
!        (phlista(jph)%constitlist(i),i=1,phlista(jph)%tnooffr),nem
      displace=displace+5
      do i=1,9
         phlista(jph)%linktocs(i)=iws(lok+displace+i)
      enddo
!      write(*,*)'3E Reading phase: ',trim(phlista(jph)%name),&
!           phlista(jph)%alphaindex,phlista(jph)%linktocs(1)
      displace=displace+9
      do i=1,nsl
         phlista(jph)%nooffr(i)=iws(lok+displace+i)
      enddo
      displace=displace+nsl
      do i=1,totcon
         phlista(jph)%constitlist(i)=iws(lok+displace+i)
      enddo
      displace=displace+totcon+1
      phlista(jph)%i2slx(1)=iws(lok+displace)
      phlista(jph)%i2slx(2)=iws(lok+displace+1)
! here are stored endmember records and additions
      phreclink=lok+displace+2
!------ endmember records, these must be allocated and linked now
      nullify(phlista(jph)%ordered)
      nullify(phlista(jph)%disordered)
      nullify(emrec)
!      if(associated(emrec)) then
!         write(*,*)'3E nullify does not work'
!         stop
!      endif
! if nem=0 now there are no basic (ordered) endmember (can that happen?)
! return here when endmember list empty and there is a disordered list
      firstendmem=1
      lokem=iws(phreclink)
!      write(*,*)'3E read endmember data',nsl,phreclink,iws(phreclink),lokem
!------------------
200   continue
      newendmem: do while(lokem.gt.0)
! this could probably be made nicer ...
         if(associated(emrec)) then
! emrec is allocated and the property record is also read
!            write(*,*)'3E next endmember',lokem,iws(lokem)
            call readendmem(lokem,iws,nsl,emrec%nextem)
            emrec=>emrec%nextem
         elseif(firstendmem.eq.1) then
!            write(*,*)'3E Read first endmember',jph
            call readendmem(lokem,iws,nsl,phlista(jph)%ordered)
            emrec=>phlista(jph)%ordered
         elseif(firstendmem.eq.2) then
            call readendmem(lokem,iws,nsl,phlista(jph)%disordered)
            emrec=>phlista(jph)%disordered
         endif
! in iws(lokem+2=noi) is the location of interaction records (if any)
         lokint=iws(lokem+2)
         level=0
         inttree: if(lokint.gt.0) then
!>>>>> 9A: first interaction record
            call readintrec(lokint,iws,level,emrec%intpointer)
            intrec=>emrec%intpointer
300         continue
! push before going to higher
            level=level+1
            stack(level)%p1=>intrec
            stack(level)%noi=lokint
! iws(lokint+1) is link to higher interaction
            higher: if(iws(lokint+1).gt.0) then
               call readintrec(iws(lokint+1),iws,level,intrec%highlink)
               intrec=>intrec%highlink
! problem pushing ....
               lokint=iws(lokint+1)
!               lokint=lokint+1
               goto 300
            endif higher
! There are no higher records, pop records from stack
350         continue
            pop: if(level.le.0) then
! no more interactions, take next endmember
               goto 390
            else
! loosing parameters when comming back from higher level
               intrec=>stack(level)%p1
               lokint=iws(stack(level)%noi)
               level=level-1
               if(lokint.gt.0) then
                  call readintrec(lokint,iws,level,intrec%nextlink)
                  intrec=>intrec%nextlink
               else
                  goto 350
               endif
               if(lokint.gt.0) goto 300
               goto 350
            endif pop
         endif inttree
390      continue
         lokem=iws(lokem)
      enddo newendmem
! list endmembers
!      emrec=>phlista(jph)%ordered
!      i1=1
!      do while(associated(emrec) .and. i1.lt.20)
!         write(*,*)'3E Found endmember ',i1
!         emrec=>emrec%nextem
!         i1=i1+1
!      enddo
! make sure the list of endmember as a null ending
      if(associated(emrec)) then
         nullify(emrec%nextem)
      endif
! we come here when no more endmembers in this list
      if(firstendmem.eq.1) then
!>>>>> 11: if nem read here is zero there are no disordered endmembers
         if(ocv()) write(*,*)'3E checking for disordered endmembers'
!         read(lin)nem,nsl
! we must nullify emrec to start a new list of endmembers
         nullify(emrec)
         lokem=iws(phreclink+1)
         if(lokem.ne.0) then
            firstendmem=2
            goto 200
         endif
      endif
!------ restore additions list
!500 continue
      lokem=phreclink+2
!      write(*,*)'3E Any addition for ',trim(phlista(jph)%name),lokem
      if(iws(lokem).gt.0) then
         lokem=iws(lokem)
         nullify(addlink)
510      continue
         if(iws(lokem+1).ge.1 .and. iws(lokem+1).le.11) then
! all phases has volume addition ...
            if(iws(lokem+1).ne.7) write(*,515)iws(lokem+1),&
                 additioname(iws(lokem+1)),trim(phlista(jph)%name)
515            format('3E Addition type ',i3,', ',a,' for ',a)
         elseif(iws(lokem+1).ne.0) then
            write(*,515)iws(lokem+1),'Unknown type            ',&
                    trim(phlista(jph)%name)
         endif
         if(iws(lokem+1).eq.INDENMAGNETIC) then
            call create_magrec_inden(nyaddlink,iws(lokem+2))
            if(gx%bmperr.ne.0) goto 1000
         elseif(iws(lokem+1).eq.XIONGMAGNETIC) then
            ifbcc=.FALSE.
            if(iws(lokem+2).eq.-1) ifbcc=.TRUE.
! ibm .TRUE. not implemented
!            write(*,*)'3E creating xiomagnetic record for BCC ',ifbcc
            call create_xiongmagnetic(nyaddlink,.FALSE.,ifbcc)
            if(gx%bmperr.ne.0) goto 1000
         elseif(iws(lokem+1).eq.VOLMOD1) then
            call create_volmod1(nyaddlink)
            if(gx%bmperr.ne.0) goto 1000
! just set it as a link, do not care if there are other additions ...
! Why this? it is done below ...
!            phlista(jph)%additions=>nyaddlink
!            nullify(nyaddlink%nextadd)
         elseif(iws(lokem+1).eq.EINSTEINCP) then
            call create_einsteincp(nyaddlink)
            if(gx%bmperr.ne.0) goto 1000
         elseif(iws(lokem+1).eq.TWOSTATEMODEL1) then
            call create_twostate_model1(nyaddlink)
            if(gx%bmperr.ne.0) goto 1000
         else
            write(*,*)'3E unknown addition'
            nullify(phlista(jph)%additions)
            goto 550
         endif
! copy the old status word
         nyaddlink%status=iws(lokem+3)
! link the additions sequentially
         
         if(associated(addlink)) then
            addlink%nextadd=>nyaddlink
         else
            phlista(jph)%additions=>nyaddlink
         endif
         nullify(nyaddlink%nextadd)
         addlink=>nyaddlink
550      continue
         lokem=iws(lokem)
         if(lokem.gt.0) goto 510
      else
         nullify(phlista(jph)%additions)
      endif
900   continue
! take next phase
      lok=iws(lok)
   enddo bigloop
! all data for the phase read
1000 continue
   kkk=jph
   return
 end subroutine readphases

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine readendmem
!\begin{verbatim}
 subroutine readendmem(lokem,iws,nsl,emrec)
! allocates and reads an endmember record and its property record from iws
! emrec is an un-allocated pointer in the parameter tree structure
   implicit none
   integer lokem,nsl,iws(*)
   type(gtp_endmember), pointer :: emrec
!\end{verbatim} %+
   integer i,j,displace,lokpty
   type(gtp_property), pointer :: proprec
!
   allocate(emrec)
!   write(*,*)'3E Allocating endmember for',lokem,iws(lokem),iws(lokem+1),&
!        iws(lokem+2)
! iws(lokem) is next endmember
! iws(lokem+1) is property
! iws(lokem+2) is interaction
!   read(lin)emrec%noofpermut,emrec%phaselink,emrec%antalem,nop,noi,nem
   emrec%noofpermut=iws(lokem+3)
   emrec%phaselink=iws(lokem+4)
   emrec%antalem=iws(lokem+5)
   displace=5
   allocate(emrec%fraclinks(nsl,emrec%noofpermut))
   do j=1,emrec%noofpermut
!      read(lin)(emrec%fraclinks(i,j),i=1,nsl)
      do i=1,nsl
         emrec%fraclinks(i,j)=iws(lokem+displace+i)
      enddo
      displace=displace+nsl
   enddo
   nullify(emrec%nextem)
   nullify(emrec%intpointer)
   nullify(emrec%propointer)
! called nop when storing in iws
   lokpty=lokem+1
   if(iws(lokpty).gt.0) then
! property list loop inside readproprec
      call readproprec(lokpty,iws,emrec%propointer)
!      write(*,*)'3E Back from readproprec 1'
   endif
1000 continue
   return
 end subroutine readendmem

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine readproprec
!\begin{verbatim}
 subroutine readproprec(lokpty,iws,firstproprec)
! allocates and a property record for both endmembers and interactions
   implicit none
   integer lokpty,iws(*)
   type(gtp_property), pointer :: firstproprec
!\end{verbatim} %+
   integer i,lokfun,displace
!   type(gtp_property), allocatable, target :: prec
   type(gtp_property), pointer :: proprec
! lokpty is the location where there can be a property record pointer
   nullify(proprec)
   do while(iws(lokpty).gt.0)
      lokpty=iws(lokpty)
      if(associated(proprec)) then
         allocate(proprec%nextpr)
         proprec=>proprec%nextpr
      else
         allocate(firstproprec)
         proprec=>firstproprec
      endif
!   read(lin)proprec%reference,proprec%proptype,&
!        proprec%degree,proprec%extra,proprec%antalprop,nox
!      write(*,88)lokpty,iws(lokpty),iws(lokpty+1),iws(lokpty+2)
      proprec%proptype=iws(lokpty+1)
      proprec%degree=iws(lokpty+2)
      proprec%extra=iws(lokpty+3)
      proprec%antalprop=iws(lokpty+4)
      call loadc(lokpty+5,iws,proprec%reference)
      displace=5+nwch(16)
      call loadc(lokpty+displace,iws,proprec%modelparamid)
!      write(*,*)'3E place to find modelparamid: ',&
!           proprec%modelparamid,lokpty+displace
! check that this is the same as the proptype!!
      i=proprec%proptype
      if(i.gt.100) i=i/100
      if(proprec%modelparamid.ne.propid(i)%symbol) then
         write(*,96)i,proprec%modelparamid,propid(i)%symbol
96       format('3E Model property ',i2,' has changed from ',&
              a,' to ',a/'Please contact Bo Sundman for help!')
!      else
! debug
!         write(*,96)i,proprec%modelparamid,propid(i)%symbol
      endif
!      lokfun=lokpty+5+nwch(16)
      lokfun=lokpty+displace+nwch(4)
! links to function as stored as integer indices
      allocate(proprec%degreelink(0:proprec%degree))
      do i=0,proprec%degree
! functions already read and hopefully stored with same index!
         proprec%degreelink(i)=iws(lokfun+i)
      enddo
!      write(*,*)'3E Allocated property record ',lokpty,iws(lokpty),&
!           proprec%proptype,proprec%degree
!      nullify(proprec%nextpr)
   enddo
! make sure the list is terminated by a null pointer
   nullify(proprec%nextpr)
1000 continue
   return
 end subroutine readproprec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine readintrec
!\begin{verbatim}
 subroutine readintrec(lokint,iws,level,intrec)
! allocates and reads an interaction record
   implicit none
   integer lokint,iws(*),level
   type(gtp_interaction), pointer :: intrec
!\end{verbatim} %+
   integer fipsize,noofperm,i,displace,lokpty,lokalint
   type(gtp_property), pointer :: proprec
! the storage of permutations in interaction records is complex ... one must
! take into account the number of permutations in lower order intecations ...
! for an fcc endmember A:A:A:B (4 perm) the binary interaction A:A:A,B:B has 
! 3; 3; 3 and 3 perms and the ternary A:A,B:A,B:B has 2; 2; 2; 2
! mult may not be needed ...
! one should never allocate a pointer ... but this is more or less permanent
   allocate(intrec)
!>>>>> 9D: actually read the interaction record
!   lokalint=iws(lokint)
   lokalint=lokint
   intrec%status=iws(lokalint+3)
   intrec%antalint=iws(lokalint+4)
   intrec%order=iws(lokalint+5)
! nullify Toop-Kohler link
   nullify(intrec%tooprec)
   fipsize=iws(lokalint+6)
!   write(*,'(a,5i5)')'3E readintrec 1:',intrec%antalint,lokalint,fipsize,level
   allocate(intrec%noofip(fipsize))
!   read(lin)intrec%noofip,intrec%status,noi,nup,nop
   displace=6
   do i=1,fipsize
      intrec%noofip(i)=iws(lokalint+displace+i)
   enddo
   displace=displace+fipsize
! 2020.06.08/BoS error saving a parameter with level=2 ??? but when saveing
! an interaction record there are only level=1 separate ???
   if(level.eq.1) then
      noofperm=intrec%noofip(fipsize)
   else 
      noofperm=intrec%noofip(2)
!   else
! I do not understand this error ...
!      write(*,*)'3E too many interaction levels for permutations',level
!      gx%bmperr=4399; goto 1000
   endif
!   if(level.eq.0) then
!      noofperm=intrec%noofip(2)
!   elseif(level.eq.1) then
!      noofperm=intrec%noofip(fipsize)
!   else
!      write(*,*)'3E too many interaction levels for permutations',level
!      gx%bmperr=4399; goto 1000
!   endif
! end of code changes 2020.06.08/BoS
   allocate(intrec%sublattice(noofperm))
   allocate(intrec%fraclink(noofperm))
!   write(*,*)'3E allocate link: ',intrec%antalint,intrec%noofip(2)
   do i=1,noofperm
      intrec%sublattice(i)=iws(lokalint+displace+2*i-1)
      intrec%fraclink(i)=iws(lokalint+displace+2*i)
   enddo
   nullify(intrec%nextlink)
   nullify(intrec%highlink)
   nullify(intrec%propointer)
! link to property record in lokalint+2
   lokpty=lokalint+2
   if(iws(lokpty).gt.0) then
      call readproprec(lokpty,iws,intrec%propointer)
!      write(*,*)'3E Back from readproprec 2'
! if there are no property records proprec is still nullified
   endif
1000 continue
   return
 end subroutine readintrec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine readequil
!\begin{verbatim}
 subroutine readequil(lokeq,iws,elope)
! subroutine readequil(lokeq,iws,elope,ceq)
! lokeq is index for equilibrium record in iws
   implicit none
   integer lokeq,iws(*),elope
!\end{verbatim} %+
   type(gtp_equilibrium_data), pointer :: ceq
   character text*2048,dum16*16,line*72,ctext*72
   type(gtp_phase_varres), pointer :: firstvarres
   TYPE(gtp_fraction_set), pointer :: fslink
   integer i,ierr,ip,isp,ivar,j,jp,k,lokcs,lokph,mc,mc2,nprop,nsl,kp,kl
   integer displace,llen,lok,lokvares,lokdiz,eqdis,lokcc,disz,conditionplace
   integer offset,lokd,dmc,eqnumber,fixph
   double precision, dimension(:,:), allocatable :: ca,ci
   double precision xxx
! containing conditions, components and phase varres records for wach compset
!>>>>> 50:
!   read(lin)ceq%eqname,ceq%eqno,ceq%status,ceq%next
!   write(*,*)'3E In readequil ',lokeq,elope
! constat
! elope=1 to read first equilibrium, -1 to read second or later
   eqnumber=1
17 continue
   if(elope.lt.0) then
! take next equilibrium and increment eqnumber
      lokeq=iws(lokeq)
      eqnumber=eqnumber+1
   endif
   ceq=>eqlista(eqnumber)
   if(lokeq.le.0) then
      if(elope.gt.0) then
! if this is the first equilibrium this is an error, otherwise just end of list
         write(*,*)'3E Not an equilibrium record'
         gx%bmperr=4399
      endif
      goto 1000
   endif
!   write(*,12)'3E Reading equilibrium ',lokeq,eqnumber,iws(lokeq+3),&
!        iws(lokeq+1)
12 format(a,3i5,1x,z8)
   ceq%status=iws(lokeq+1)
! set that no calculation is made in status word to prevent listing ?? why ??
!   ceq%status=ibset(ceq%status,EQNOEQCAL)
   ceq%multiuse=iws(lokeq+2)
! Hm, eqno should not be changed?  By default arbitrary value!!
   if(eqnumber.ne.iws(lokeq+3)) then
      write(*,*)'3E Should be same equilibrium number ',eqnumber,iws(lokeq+3)
   endif
   ceq%eqno=iws(lokeq+3)
   ceq%nexteq=iws(lokeq+4)
   call loadc(lokeq+5,iws,ceq%eqname)
!   write(*,*)'3E Reading equilibrium with name: ',ceq%eqname
   displace=5+nwch(24)
   call loadc(lokeq+displace,iws,ceq%comment)
!   write(*,*)'3E comment: "',trim(ceq%comment),'" ',len_trim(ceq%comment)
   displace=displace+nwch(72)
! values of T and P and weight
   call loadrn(2,iws(lokeq+displace),ceq%tpval)
   call loadr(lokeq+displace+2*nwpr,iws,ceq%weight)
   displace=displace+3*nwpr
!----- components (must be elements).  Must be entered before conditions
! complist already allocated for 20
   if(allocated(ceq%complist)) then
      deallocate(ceq%complist)
   endif
   allocate(ceq%complist(noofel))
   if(eqnumber.gt.1) then
      allocate(ceq%compstoi(noofel,noofel))
      allocate(ceq%invcompstoi(noofel,noofel))
   endif
   ceq%compstoi=zero
   ceq%invcompstoi=zero
   do kl=1,noofel
! when the elements are components ...
      ceq%compstoi(kl,kl)=one
      ceq%invcompstoi(kl,kl)=one
   enddo
   llen=0
   lokcc=iws(lokeq+displace+2)
   do while(lokcc.gt.0)
      llen=llen+1
      if(llen.gt.noofel) then
         write(*,*)'3E Too many components'
         gx%bmperr=4399; goto 1000
      endif
      ceq%complist(llen)%splink=iws(lokcc+1)
      ceq%complist(llen)%phlink=iws(lokcc+2)
      ceq%complist(llen)%status=iws(lokcc+3)
      call loadc(lokcc+4,iws,ceq%complist(llen)%refstate)
      disz=4+nwch(16)
      kl=iws(lokcc+disz)
      if(kl.gt.0) then
         allocate(ceq%complist(llen)%endmember(kl))
         do mc=1,kl
            ceq%complist(llen)%endmember(mc)=iws(lokcc+disz+mc)
         enddo
!         write(*,*)'3E endmem: ',kl,(ceq%complist(llen)%endmember(mc),mc=1,kl)
         disz=disz+kl+1
      else
         disz=disz+1
      endif
      call loadrn(2,iws(lokcc+disz),ceq%complist(llen)%tpref)
!      write(*,*)'3E refstate 2: ',ceq%complist(llen)%tpref
      disz=disz+2*nwpr
      call loadrn(2,iws(lokcc+disz),ceq%complist(llen)%chempot)
      disz=disz+2*nwpr
      call loadr(lokcc+disz,iws,ceq%complist(llen)%mass)
!      write(*,*)'3E loading component mass',lokcc,disz,llen,&
!           ceq%complist(llen)%mass
      call loadr(lokcc+disz+nwpr,iws,ceq%complist(llen)%molat)
      lokcc=iws(lokcc)
   enddo
!----- conditions (note that inactive conditions not set)
! conditions cannot be entered before the phase_varres for all phases
   conditionplace=displace
!----------- phase_varres record
!>>>>> 54:
   highcs=iws(lokeq+displace+3)
   if(ocv()) then
      write(*,*)'3E Number of phase_varres records: ',highcs
      write(*,*)'phase_varres size: ',size(ceq%phase_varres)
   endif
! link to first varres record stored here
   lokvares=iws(lokeq+displace+4)
!   write(*,*)'3E lokvares: ',lokvares,highcs,lokeq,displace+4
   eqdis=displace+5
! for equilibria 2 and higher phase_varees must be allocated!!
   if(eqnumber.gt.1) then
!      write(*,*)'3E allocating phase_varres for equilibrium: ',eqnumber
      allocate(ceq%phase_varres(highcs+5))
! we should also allocate a few other things
      allocate(ceq%eq_tpres(maxtpf))
      allocate(ceq%svfunres(maxsvfun))
      do j=1,maxtpf
         ceq%eq_tpres(j)%forcenewcalc=0
      enddo
      ceq%tpval(1)=1.0D3
      ceq%tpval(2)=1.0D5
   endif
   compset: do j=1,highcs
      if(lokvares.le.100) then
         write(*,*)'3E error linking phase_varres records ...',lokvares,j
         goto 1000
      endif
!------------------------------------------
! DEBUGPROBLEM BEWARE, using = instead of => below took 2 days to find
!------------------------------------------
! >>>      firstvarres=ceq%phase_varres(j)    <<< error
      firstvarres=>ceq%phase_varres(j)
!>>>>> 55:
      firstvarres%nextfree=iws(lokvares+1)
      lokph=iws(lokvares+2)
      if(lokph.lt.0) then
! this means this phase_varres record is not used
! we have already save the free list link, just skip the rest
         write(*,*)'3E found unused phase_varres record: ',j,lokvares
         lokvares=iws(lokvares)
         cycle compset
      endif
      firstvarres%phlink=lokph
      firstvarres%status2=iws(lokvares+3)
      firstvarres%phstate=iws(lokvares+4)
      firstvarres%phtupx=iws(lokvares+5)
      nsl=iws(lokvares+6)
      mc=iws(lokvares+7)
!      write(*,*)'3E read mc ',trim(phlista(lokph)%name),nsl,mc,j
      call loadc(lokvares+8,iws,firstvarres%prefix)
      displace=8+nwch(4)
      call loadc(lokvares+displace,iws,firstvarres%suffix)
      displace=displace+nwch(4)
      call loadrn(3,iws(lokvares+displace),firstvarres%abnorm)
      displace=displace+3*nwpr
! we need these values here! but now they are stored in iws!!
!      nsl=phlista(lokph)%noofsubl
!      mc=phlista(lokph)%tnooffr
      mc2=mc*(mc+1)/2
      if(btest(firstvarres%status2,CSDLNK)) then
! varres record with link to disordered varres record, some data to be stored
! NOTE necessary data for nsl and mc stored later ...
! we need these values here!
!         write(*,*)'3E varres with link to disordered fraction varres'
         offset=6+2*nwch(4)+3*nwpr+mc*(1+2*nwpr)+nsl*nwpr
!         write(*,202)'3E offset:',j,lokvares,displace,iws(lokvares+displace),&
!              nsl,mc,offset,iws(lokvares+offset),iws(lokvares+26)
202      format(a,10i6)
!         stop
!      elseif(btest(firstvarres%status2,CSDFS)) then
!         write(*,*)'3E varres for disordered fraction set OK',j,nsl,mc
! this phase_varres/parres record belong to disordered fraction_set
! we should use nsl and mc from disordered fraction set!
! but they are not yet created...
      endif
!      write(*,88)'3E reading phase_varres ',j,highcs,lokvares,nsl,mc,&
!           trim(phlista(lokph)%name)
!88    format(a,5i7,2x,a)
!
!      write(*,*)'3E allocating constat: ',j,mc
      allocate(firstvarres%constat(mc))
      do i=1,mc
         firstvarres%constat(i)=iws(lokvares+displace+i-1)
      enddo
      displace=displace+mc
      allocate(firstvarres%yfr(mc))
      call loadrn(mc,iws(lokvares+displace),firstvarres%yfr)
      displace=displace+mc*nwpr
!      write(*,*)'3E not allocating mmyfr'
      displace=displace+mc*nwpr
      allocate(firstvarres%sites(nsl))
!      write(*,*)'3E sites: ',lokvares,displace,lokvares+displace
      call loadrn(nsl,iws(lokvares+displace),firstvarres%sites)
      displace=displace+nsl*nwpr
!-----------------------------------
! BEWHERE allocation of the dpqdy and d2pqdvay!!! 
! They are not saved but should be allocated here! need lokph
      if(btest(phlista(lokph)%status1,PHIONLIQ)) then
!         write(*,*)'3E ionic liquid',lokph,eqnumber
         allocate(firstvarres%dpqdy(mc))
         allocate(firstvarres%d2pqdvay(mc))
         firstvarres%dpqdy=zero
         firstvarres%d2pqdvay=zero
      endif
!-------------------------------------
!      write(*,*)'3E odd:   ',lokvares,displace
      fsrec: if(btest(firstvarres%status2,CSDLNK)) then
!         write(*,*)'3E disfra record:',lokvares,displace,iws(lokvares+displace)
! disfra record
         lokd=iws(lokvares+displace)
         fslink=>firstvarres%disfra
         fslink%latd=iws(lokd)
         nsl=iws(lokd+1)
         fslink%ndd=nsl
         dmc=iws(lokd+2)
         fslink%tnoofxfr=dmc
         fslink%tnoofyfr=iws(lokd+3)
         fslink%totdis=iws(lokd+4)
         fslink%varreslink=iws(lokd+5)
         call storc(lokd+6,iws,fslink%id)
         disz=6+nwch(1)
         allocate(fslink%nooffr(nsl))
         allocate(fslink%splink(dmc))
         allocate(fslink%y2x(mc))
         allocate(fslink%dsites(nsl))
         allocate(fslink%dxidyj(mc))
         disz=6+nwch(1)
         do i=1,nsl
            fslink%nooffr(i)=iws(lokd+disz+i)
         enddo
         disz=disz+nsl
!         write(*,202)'3E disfra 1: ',lokd,disz
         do i=1,dmc
            fslink%splink(i)=iws(lokd+disz+i)
         enddo
         disz=disz+dmc+1
! we must use the ordered number of constituents here!!
         if(mc.ne.iws(lokd+disz)) then
            write(*,*)'3E constituent number error: ',mc,iws(lokd+disz)
            mc=iws(lokd+disz)
         endif
!         write(*,202)'3E disfra 2: ',lokd,disz
         do i=1,mc
            fslink%y2x(i)=iws(lokd+disz+i)
         enddo
         disz=disz+mc+1
!         write(*,202)'3E disfra 3: ',lokd,disz
         call loadrn(nsl,iws(lokd+disz),fslink%dsites)
         disz=disz+nsl*nwpr
!         write(*,202)'3E disfra 4: ',lokd,disz
         call loadrn(mc,iws(lokd+disz),fslink%dxidyj)
         disz=disz+mc*nwpr
!         write(*,202)'3E disfra 5: ',lokd,disz
         call loadr(lokd+disz,iws,fslink%fsites)
!         write(*,*)'3E disfra last: ',lokd+disz+nwpr
      else
         firstvarres%disfra%varreslink=0
      endif fsrec
      displace=displace+1
      call loadr(lokvares+displace,iws,firstvarres%amfu)
      call loadr(lokvares+displace+nwpr,iws,firstvarres%netcharge)
      call loadr(lokvares+displace+2*nwpr,iws,firstvarres%dgm)
! NEW value of qcbonds for quasichemical model, increment of displace!!
      call loadr(lokvares+displace+3*nwpr,iws,firstvarres%qcbonds)
      displace=displace+4*nwpr
      nprop=iws(lokvares+displace)
      if(nprop.lt.20) then
!         write(*,303)'3E get nprop: ',lokvares,displace,lokvares+displace,&
!              nprop,trim(phlista(lokph)%name)
         nprop=20
      endif
303   format(a,4i7,2x,a)
      firstvarres%nprop=nprop
      allocate(firstvarres%listprop(nprop))
! calculated results, only G saved
      allocate(firstvarres%gval(6,nprop))
      displace=displace+1
! we have saved only G values
      do i=1,6
         call loadr(lokvares+displace+nwpr*(i-1),iws,firstvarres%gval(i,1))
      enddo
      displace=displace+6*nwpr
      allocate(firstvarres%dgval(3,mc,nprop))
      do i=1,3
         do k=1,mc
            call loadr(lokvares+displace,iws,firstvarres%dgval(i,k,1))
            displace=displace+nwpr
         enddo
      enddo
      allocate(firstvarres%d2gval(mc2,nprop))
      do i=1,mc2
         call loadr(lokvares+displace+nwpr*(i-1),iws,firstvarres%d2gval(i,1))
      enddo
! link to next stored phase_varres record
      lokvares=iws(lokvares)
   enddo compset
!   if(elope.lt.0) then
      csfree=highcs+1
!   endif
!   write(*,*)'3E csfree: ',highcs,csfree,elope
!   write(*,*)'3E All phase_varres records created for ',ceq%eqno
!----- conditions (note that inactive conditions not set)
!   lok=iws(lokeq+displace)
   lok=iws(lokeq+conditionplace)
   nullify(ceq%lastcondition)
   nullify(ceq%lastexperiment)
   if(lok.gt.0) then
      llen=iws(lok)
      call loadc(lok+1,iws,text(1:llen))
!      write(*,*)'3E Conditions: "',text(1:llen),'"',llen
      if(llen.gt.0) then
! set the conditions, kp will be incremented by 1 in enter_condition
! the text contains " number: variable expression=value, "
! we have to set each condition separately.  There can be , but no :
! in the variable expressions.
         jp=1; ip=llen
         cloop: do while(jp.lt.ip)
            k=index(text(jp:ip),':')
            if(k.le.0) exit cloop
            line=text(jp+k:ip)
            jp=jp+k+2
! remove any commma followed by space ", " as that indicates there are more 
! conditions on the same line
            kp=index(line,', ')
            if(kp.gt.0) then
               line(kp:)=' '
            else
               kp=index(line,' ')
               line(kp:)=' '
            endif
! We must handle fix phases :: <phase>=value transforms to fix=phase == value
            if(line(1:1).eq.'<') then
               kp=index(line,'>')
               fixph=kp+1
               call getrel(line,fixph,xxx)
               if(buperr.ne.0) then
                  buperr=0; xxx=zero
               endif
               ctext=' FIX='//line(2:kp-1)//' == '//line(fixph+1:)
!               write(*,*)'3E fixph: ',trim(ctext)
               line=ctext
            endif
            kp=0
!            write(*,*)'3E set condition "',trim(line),'"',jp,ip
            call set_condition(line,kp,ceq)
!            write(*,*)'3E back from set condition "',gx%bmperr
            if(gx%bmperr.ne.0) then
               write(*,*)'3E Error setting conditions'
               write(*,*)'3E condition "',trim(line),'"',kp
               goto 1000
            endif
         enddo cloop
      endif
!   else
!      write(*,*)'3E no conditions on unformatted file'
   endif
!----- experiments
   lok=iws(lokeq+conditionplace+1)
733 continue
   kp=0
   if(lok.gt.0) then
! experiments are stored individually in a linked list
      kp=kp+1
      llen=iws(lok+1)
      text=' '
      call loadc(lok+2,iws,text(1:llen))
!      write(*,*)'3E found experiment: "',trim(text),'"'
      llen=0
      call enter_experiment(text,llen,ceq)
!      write(*,*)'3E Back from enter_experiment'
      if(gx%bmperr.ne.0) then
         write(*,*)'3E error entering experiment ',gx%bmperr,' continuing'
         gx%bmperr=0
      endif
      lok=iws(lok)
      goto 733
   endif
   if(kp.gt.0) write(*,*)'3E Found ',kp,' experiments'
!-------------------------- a few remaining things
   ceq%maxiter=iws(lokeq+eqdis)
   if(.not.allocated(ceq%cmuval)) then
      allocate(ceq%cmuval(noofel))
   endif
   call loadrn(noofel,iws(lokeq+eqdis+1),ceq%cmuval)
   eqdis=eqdis+1+noofel*nwpr
   call loadr(lokeq+eqdis,iws,ceq%xconv)
! modifed 2018.05.28 by adding gdconv(2)
   call loadr(lokeq+eqdis+nwpr,iws,ceq%gdconv(1))
   call loadr(lokeq+eqdis+2*nwpr,iws,ceq%gdconv(2))
   call loadr(lokeq+eqdis+3*nwpr,iws,ceq%gmindif)
! if elope negative continue reading next equilibrium
   if(elope.lt.0) then
!      write(*,*)'3E read the next equilibrium'
! increment the index of first free equilibrium
      eqfree=eqfree+1
      goto 17
   endif
!
1000 continue
   if(eqfree.gt.2) write(*,1010)eqfree-1
1010 format('3E Read ',i4,' equilibria')
   return
 end subroutine readequil

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine svfunread
!\begin{verbatim}
 subroutine svfunread(loksvf,iws)
! read a state variable function from save file and store it.
! by default there are some state variable functions, make sure
! they are deleted.  Done here just by setting nsvfun=0
   implicit none
   integer loksvf,iws(*)
!\end{verbatim} %+
   integer nsvfun,i,ip,lok,eqno
   character*512 text
   nsvfun=iws(loksvf)
! first 3 symbols are R, RT and T_C
   do i=iws(loksvf+1)+1,nsvfun
      lok=iws(loksvf+i)
      ip=iws(lok)
      text=' '
      call loadc(lok+1,iws,text(1:ip))
!      write(*,*)'3E Entering saved svf: "',text(1:ip),'"'
! NOTE: position 1-7 are equilibrium number and status
      ip=7
      call enter_svfun(text,ip,firsteq)
      if(gx%bmperr.ne.0) then
         write(*,*)'3E Error entering saved svf',gx%bmperr
         if(gx%bmperr.ne.4136) goto 1000
         gx%bmperr=0
      endif
! if this function should be evaluated at a particular equilibrium that is
! in position 1-5.  Extra status in position 6 and 7
!      write(*,*)'3E read symbol: ',i,': ',text(1:ip),ip
! Letters used for the status bits:
! A SVNOAM a function or constant that cannot be amended
! C SVCONST a constant that can be amended
! D SVFDOT a dot derivative (also SVFVAL set)
! N SVFVAL symbol evaluated only if explitly referenced)
! X SVFEXT only evaluated for a specific equilibrium (preceeded by eq.number)
! I SVIMPORT import value from TP function (preceeded by TP index)
! E SVEXPORT expert value to TP function constant (preceeded by TP index)
! check if symbol is a constant (can be amended)
      if(text(5:5).eq.'C') then
         svflista(i)%status=ibset(svflista(i)%status,SVNOAM)
      elseif(text(5:5).eq.'C') then
         svflista(i)%status=ibset(svflista(i)%status,SVCONST)
! check if symbol should only be evaluated when explicitly requested
      elseif(text(5:5).eq.'D') then
! D means the symbol is a dot variable, evaluates only when explitly refered
         svflista(i)%status=ibset(svflista(i)%status,SVFDOT)
         svflista(i)%status=ibset(svflista(i)%status,SVFVAL)
      elseif(text(5:5).eq.'V') then
         svflista(i)%status=ibset(svflista(i)%status,SVFVAL)
      endif
! extract any number before postion 5
      ip=0
! ip is incremented in getint
      call getint(text,ip,eqno)
      if(buperr.ne.0) then
         buperr=0
      else
         if(text(5:5).eq.'X') then
! symbol should be evaluated at a specific equilibrium (eqno)
            svflista(i)%status=ibset(svflista(i)%status,SVFEXT)
            svflista(i)%eqnoval=eqno
         elseif(text(5:5).eq.'I') then
! symbol should be imported from TP function
            svflista(i)%status=ibset(svflista(i)%status,SVIMPORT)
            svflista(i)%tplink=eqno
         elseif(text(5:5).eq.'E') then
! symbol should be exported to TP constant
            svflista(i)%status=ibset(svflista(i)%status,SVEXPORT)
            svflista(i)%tplink=eqno
         else
! a number with no meaning!
            write(*,*)trim(text)
100         format(' *** Warning, error reading symbol:'/a)
         endif
      endif
   enddo
1000 continue
   return
 end subroutine svfunread

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine biblioread
!\begin{verbatim}
 subroutine biblioread(bibhead,iws)
! read references from save file
   implicit none
   integer bibhead,iws(*)
!\end{verbatim} %+
   character text*2048
   integer i,iref,jp,nrefs,lok,kk,ir,nr
!>>>>> 40: number of references
!   write(*,*)'3E Reading reference version and nummer of'
   nrefs=iws(bibhead)
   do i=1,nrefs
      lok=iws(bibhead+i)
      jp=iws(lok)
      call loadc(lok+1,iws,text(1:jp))
      call tdbrefs(text(1:16),text(17:jp),0,iref)
   enddo
1000 continue
   return
 end subroutine biblioread

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\addtotable subroutine readash
!\begin{verbatim}
 subroutine readash(lok,iws)
! reading assessment records
   integer lok,iws(*)
!\end{verbatim}
   integer lok1,lok2,last,rsize,i1,i2,disp,kk
   double precision xxx
   type(gtp_assessmenthead), pointer :: assrec
   type(gtp_equilibrium_data), pointer :: ceq
!
   lok1=lok
   assrec=>firstash%nextash
20 continue
   if(iws(lok1).eq.0) goto 1000
   lok1=iws(lok1)
   assrec%status=iws(lok1+1)
   assrec%varcoef=iws(lok1+2)
   assrec%firstexpeq=iws(lok1+3)
   assrec%lwam=iws(lok1+4)
   call loadc(lok1+5,iws,assrec%general)
   disp=5+nwch(64)
   call loadc(lok1+disp,iws,assrec%special)
   disp=disp+nwch(64)
   lok2=iws(lok1+disp+1)
   if(lok2.gt.0) then
! eqlista
!      lok2=iws(lok2)
      i1=iws(lok2)
      if(i1.gt.0) then
!         write(*,'(a,4i10)')'3E In readash 1: ',lok,lok1,lok2,i1
         allocate(assrec%eqlista(i1))
! in iws(lok2+i2) the index to eqlista is stored, 
! assrec%eqlista(i2)%p1 is a pointer to this equilibrium
         do i2=1,i1
            ceq=>eqlista(iws(lok2+i2))
            assrec%eqlista(i2)%p1=>ceq
         enddo
      endif
   else
      write(*,*)'3E no experimental data'
   endif
   lok2=iws(lok1+disp+2)
   if(lok2.le.0) then
      write(*,*)'3E no coefficient values saved'
      goto 777
   else
! coeffvalues
!      lok2=iws(lok2)
      i1=iws(lok2)
!      write(*,*)'3E In readash 2: ',lok2,i1
      allocate(assrec%coeffvalues(0:i1-1))
      call loadrn(i1,iws(lok2+1),assrec%coeffvalues)
   endif
   lok2=iws(lok1+disp+3)
! coeffrsd
!      lok2=iws(lok2)
   if(lok2.gt.0) then
      i1=iws(lok2)
!      write(*,*)'3E In readash RSD: ',lok2,i1
      allocate(assrec%coeffrsd(0:i1-1))
      call loadrn(i1,iws(lok2+1),assrec%coeffrsd)
   endif
!   lok2=iws(lok1+disp+3)
   lok2=iws(lok1+disp+4)
   if(iws(lok2).gt.0) then
! coeffscale
!      lok2=iws(lok2)
      i1=iws(lok2)
!      write(*,*)'3E In readash 3: ',lok2,i1
      allocate(assrec%coeffscale(0:i1-1))
      call loadrn(i1,iws(lok2+1),assrec%coeffscale)
   endif
!   lok2=iws(lok1+disp+4)
   lok2=iws(lok1+disp+5)
   if(iws(lok2).gt.0) then
! coeffstart
!      lok2=iws(lok2)
      i1=iws(lok2)
!      write(*,*)'3E In readash 4: ',lok2,i1
      allocate(assrec%coeffstart(0:i1-1))
      call loadrn(i1,iws(lok2+1),assrec%coeffstart)
   endif
!   lok2=iws(lok1+disp+5)
   lok2=iws(lok1+disp+6)
   if(iws(lok2).gt.0) then
! coeffmin
!      lok2=iws(lok2)
      i1=iws(lok2)
!      write(*,*)'3E In readash 5: ',lok2,i1
      allocate(assrec%coeffmin(0:i1-1))
      call loadrn(i1,iws(lok2+1),assrec%coeffmin)
   endif
!   lok2=iws(lok1+disp+6)
   lok2=iws(lok1+disp+7)
   if(iws(lok2).gt.0) then
! coeffmax
!      lok2=iws(lok2)
      i1=iws(lok2)
!      write(*,*)'3E In readash 6: ',lok2,i1
      allocate(assrec%coeffmax(0:i1-1))
      call loadrn(i1,iws(lok2+1),assrec%coeffmax)
   endif
!   lok2=iws(lok1+disp+7)
   lok2=iws(lok1+disp+8)
   if(iws(lok2).gt.0) then
! coeffindices
!      lok2=iws(lok2)
      i1=iws(lok2)
!      write(*,*)'3E In readash 7: ',lok2,i1
      allocate(assrec%coeffindex(0:i1-1))
      do i2=1,i1
         assrec%coeffindex(i2-1)=iws(lok2+i2)
      enddo
! store these values in tpfun ...
      do kk=0,i1-1
!         write(*,333)'3E storing as TP funs ',kk,assrec%coeffindex(kk),&
!              assrec%coeffvalues(kk),assrec%coeffscale(kk)
333      format(a,2i4,6(1pe12.4))
! firstash or assrec??
         xxx=assrec%coeffvalues(kk)*assrec%coeffscale(kk)
         call change_optcoeff(assrec%coeffindex(kk),xxx)
         if(gx%bmperr.ne.0) goto 1000
      enddo
   endif
!   lok2=iws(lok1+disp+8)
   lok2=iws(lok1+disp+9)
   if(iws(lok2).gt.0) then
! coeffstate
!      lok2=iws(lok2)
      i1=iws(lok2)
!      write(*,*)'3E In readash 8: ',lok2,i1
      allocate(assrec%coeffstate(0:i1-1))
      do i2=1,i1
         assrec%coeffstate(i2-1)=iws(lok2+i2)
      enddo
   endif
777 continue
! maybe work array has been daved also?
!   lok2=iws(lok1+disp+9)
   lok2=iws(lok1+disp+10)
   if(lok2.gt.0) then
      if(iws(lok2).gt.0) then
!         lok2=iws(lok2)
         i1=iws(lok2)
!         write(*,*)'3E In readash 9: ',lok2,i1
         allocate(assrec%wopt(i1))
         call loadrn(i1,iws(lok2+1),assrec%wopt)
      endif
   endif
! check if there are several assessmentheads
   if(iws(lok1).gt.0) then
! There are more records, try to create a circular list in both directions
      write(*,*)'3E In readash 10: ',lok1,iws(lok1)
      allocate(assrec%nextash)
      assrec%nextash%prevash=>assrec
      assrec=>assrec%nextash
      firstash%prevash=>assrec
      write(*,*)'3E more assessment records'
      goto 20
   endif
1000 continue
   return
 end subroutine readash

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\addtotable subroutine new_gtp
!\begin{verbatim}
 subroutine new_gtp
!
! DELETES ALL DATA so a new TDB file can be read
!
! this is needed before reading a new unformatted file (or same file again)
! we must go through all records and delete and deallocate each
! separately.  Very similar to gtpread
   implicit none 
!\end{verbatim}
   integer isp,j,nel,intv(10),k
   double precision dblv(10)
   TYPE(gtp_equilibrium_data), pointer :: ceq
   type(gtp_phase_varres), pointer :: phdyn
!   TYPE(gtp_fraction_set) :: fslink
!   write(*,*)'3E Testing segmentation error in new_gtp'
   if(ocv()) write(*,*)'3E Removing current data'
!---------- elementlist, no need to delete, just deallocate below
!>>>>> 2:
!---------- specieslist, we have to deallocate ?? maybe not ??
!>>>>> 3:
   if(btest(globaldata%status,GSNODATA)) then
      if(ocv()) write(*,*)'3E No thermodynamic data to delete'
      goto 600
   endif
   if(gtp_species_version.ne.2) then
      write(*,17)'3E *** ERROR species',1,gtp_species_version
17    format(a,' record version error: ',2i4)
      gx%bmperr=4300; goto 1000
   endif
   ceq=>firsteq
!   write(*,*)'3E No segmentation error A'
   do isp=1,noofsp
      nel=splista(isp)%noofel
      deallocate(splista(isp)%ellinks)
      deallocate(splista(isp)%stoichiometry)
      if(allocated(splista(isp)%spextra)) deallocate(splista(isp)%spextra)
   enddo
!---------- phases, many records, here we travese all endmembers etc
!>>>>> 4
!   write(*,*)'3E No segmentation error B'
!   if(gtp_phase_version.ne.1) then
!      write(*,17)'3E **** ERROR phase',1,gtp_phase_version
!      gx%bmperr=4302; goto 1000
!   endif
!   if(gtp_endmember_version.ne.1) then
!      write(*,17)'3E **** ERROR endmember',1,gtp_endmember_version
!      gx%bmperr=4302; goto 1000
!   endif
!   if(gtp_interaction_version.ne.1) then
!      write(*,17)'3E **** ERROR interaction',1,gtp_interaction_version
!      gx%bmperr=4302; goto 1000
!   endif
!   if(gtp_property_version.ne.1) then
!      write(*,17)'3E **** ERROR property',1,gtp_property_version
!      gx%bmperr=4302; goto 1000
!   endif
   do j=0,noofph
      call delphase(j)
      if(gx%bmperr.ne.0) goto 1000
   enddo
!   write(*,*)'3E No segmentation error C1'
!----------- jump here if no thermodynamic data
600 continue
!---------- equilibrium records
!>>>>> 50: equilibrium records
!   call delete_equil(ceq)
!   do j=1,noofeq
! this loop was added in an attempt to get rid of an error occuring with
! 64 bit version, the TP functions was not cleared correctly
   do j=1,eqfree-1
      ceq=>eqlista(j)
      deallocate(ceq%svfunres)
!      write(*,*)'3E No segmentation error C2',j
      deallocate(ceq%eq_tpres)
!      write(*,*)'3E No segmentation error C3',j
      deallocate(ceq%complist)
!      write(*,*)'3E No segmentation error C4',j
      deallocate(ceq%compstoi)
!      write(*,*)'3E No segmentation error C5',j
      deallocate(ceq%invcompstoi)
!      write(*,*)'3E No segmentation error C6',j
! remove valgrind memory leak for conditions
      call delete_all_conditions(0,ceq)
! clean upp phase_varres records
      do k=1,size(ceq%phase_varres)
         phdyn=>ceq%phase_varres(k)
         if(allocated(phdyn%gval)) then
            deallocate(phdyn%gval)
            deallocate(phdyn%dgval)
            deallocate(phdyn%d2gval)
!            write(*,*)'3E No segmentation error C7',j,k
         endif
! set phstate and phlink to zero to avoid segmentation fault when plotting
! after several MAP or STEP commands with different composition sets
         phdyn%phstate=0
         phdyn%phlink=0
      enddo
!      write(*,*)'3E No segmentation error C8',j
!      deallocate(ceq%phase_varres)
   enddo
!   write(*,*)'3E No segmentation error D1'
! I am not sure if this really releases all memory, how to check .... ???
!   call deallocate_gtp(intvar,dblvar)
   deallocate(eqlista)
!   write(*,*)'3E No segmentation error D2'
!------- deallocate elements, species and phases, will be allocated in init_gtp
   deallocate(ellista)
   deallocate(elements)
!   do k=1,noofsp
!      deallocate(splista(k)%ellinks)
!   enddo
   deallocate(splista)
   deallocate(species)
   deallocate(phlista)
   deallocate(phases)
   deallocate(phasetuple)
!   write(*,*)'3E No segmentation error E'
!------ tpfunction expressions and other lists
!>>>>> 20: delete tpfuns
!   write(*,*)'3E Delete TP funs, just deallocate??',freetpfun
!   call delete_all_tpfuns
   call tpfun_deallocate
   if(gx%bmperr.ne.0) then
      write(*,*)'3E **** ERROR deleting TP functions'
   endif
!   write(*,*)'3E Back from deleting all TP funs, this is fun!!'
!------ tpfunction expressions and other lists
!>>>>> 30: delete state variable functions
   deallocate(svflista)
!   write(*,*)'3E No segmentation error F'
!   call delete_svfuns
!---------- delete bibliographic references
!>>>>> 40: references
   deallocate(bibrefs)
!   call delete_biblio
!------ parameter property records
   deallocate(propid)
!------ other things such as mqmq_data arrays
   if(allocated(mqmqa_data%contyp)) then
      deallocate(mqmqa_data%contyp)
      deallocate(mqmqa_data%constoi)
      mqmqa_data%nconst=0
   endif
!   write(*,*)'3E No segmentation error G'
!------ map results are deleted separately
!   call delete_mapresults(maptop)
!    deallocate( .... any more ???
!---------------------------
! now initiate all lists and a little more
   if(ocv()) write(*,*)'3E All data structures will be reinitiated'
! intv(1) negative means reinititate with same values as before
!   intv(1)=-1
!   write(*,*)'3E No segmentation error H', moved to pmon
!   call init_gtp(intv,dblv)
! after return firsteq must be initiated ... maybe it should be done here ??
!
! Problem when adding EEC, initialization does not work, why?
! Maybe these need initiating?
   globaldata%sysreal=zero
   globaldata%sysparam=0
! initiate Toop/Kohler record counter
   uniqid=0
!
1000 continue
   return
 end subroutine new_gtp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine deallocate_gtp
!\begin{verbatim}
 subroutine deallocate_gtp(intvar,dblvar)
! deallocate the data structure
   implicit none
   integer allocatestatus
   integer intvar(*)
   double precision dblvar(*)
!\end{verbatim}
!   integer jl
   write(*,*)'3A in deallocate_gtp'
   deallocate(ellista, STAT = allocateStatus)
   if (allocateStatus /= 0) then
     write(kou,*) 'Error during deallocation of ellista'
     goto 1000
   else
     write(kou,*) 'Deallocation of data ',  allocateStatus
   endif
!   flush(6)
   deallocate(elements)
! deallocate records for species
   deallocate(splista)
   deallocate(species)
! deallocate records for phases
   deallocate(phlista)
   deallocate(phases)
   deallocate(phasetuple)
   deallocate(bibrefs)
   deallocate(propid)
   deallocate(eqlista)
   deallocate(svflista)
   write(*,*)'3A Deallocate TP funs'
   call tpfun_deallocate
!CCI added this
   deallocate(firstash)
1000 continue
   return
 END subroutine deallocate_gtp

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\addtotable subroutine delphase
!\begin{verbatim}
 subroutine delphase(lokph)
! save data for phase at location lokph (except data in the equilibrium record)
! For phases with disordered set of parameters we must access the number of
! sublattices via firsteq
   implicit none
   integer lokph
!\end{verbatim}
   integer level,nsl,noendm
   type(gtp_endmember), pointer :: emrec,nextem
   type(gtp_interaction), pointer :: intrec,nextint
   type(gtp_property), pointer :: proprec,nextprop
! to keep track of interaction records
   type saveint
      type(gtp_interaction), pointer :: p1
   end type saveint
   type(saveint), dimension(:), pointer :: stack
   type(gtp_phase_add), pointer :: addlink,nextadd
!   write(*,*)'3E In delphase',lokph
   allocate(stack(5))
   nsl=phlista(lokph)%noofsubl
!>>>>> 6:
   deallocate(phlista(lokph)%nooffr)
   deallocate(phlista(lokph)%constitlist)
   emrec=>phlista(lokph)%ordered
   noendm=0
!>>>>> 6: sublattice info
! we come back here if there are disordered parameters
200 continue
! there can be phases without any parameters ...
   emlista: do while(associated(emrec))
      proprec=>emrec%propointer
      intrec=>emrec%intpointer
      nextem=>emrec%nextem
!>>>>> 7: after saving links deallocate endmember record with all its content
!      write(*,*)'3E deallocate endmember record'
      deallocate(emrec)
! nextem do not need to be declared as target??
      emrec=>nextem
      emproplista: do while(associated(proprec))
         nextprop=>proprec%nextpr
!>>>>> 8: endmember property records
! functions and references deallocated separately
!         write(*,*)'3E deallocate endmember property record'
         deallocate(proprec)
         proprec=>nextprop
      enddo emproplista
! interaction tree
      level=0
300   continue
      intlista: do while(associated(intrec))
!>>>>> 9: interaction record
         level=level+1
         if(level.gt.5) then
            gx%bmperr=4164; goto 1000
         endif
!         write(*,*)'3E Pushing ',level
         stack(level)%p1=>intrec%nextlink
         nextint=>intrec%highlink
         proprec=>intrec%propointer
!         write(*,*)'3E deallocate interaction record'
         deallocate(intrec)
         intproplista: do while(associated(proprec))
            nextprop=>proprec%nextpr
!>>>>> 10: interaction properties
!            write(*,*)'3E deallocate interaction property record'
            deallocate(proprec)
            proprec=>nextprop
         enddo intproplista
         intrec=>nextint
      enddo intlista
! pop the link to next interaction if any
      pop: if(level.gt.0) then
!         write(*,*)'3E popping interaction record',level
         intrec=>stack(level)%p1
         nullify(stack(level)%p1)
         level=level-1
         goto 300
      endif pop
!---- next endmember
      emrec=>nextem
   enddo emlista
! no more endmembers, check if the disordered (if any) has been written
   if(noendm.eq.0) then
! we do not have to care about that nsl is different ....
!>>>>> 11: disordered endmembers
!      write(*,*)'3E disordered endmembers'
      emrec=>phlista(lokph)%disordered
      noendm=1
      goto 200
   endif
!   write(*,*)'3E finished parameter records'
!------ additions list
500 continue
   addlink=>phlista(lokph)%additions
   addition: do while(associated(addlink))
!>>>>> 12: additions
      nextadd=>addlink%nextadd
      if(addlink%type.eq.1) then
!>>>>> 12A: delete magnetic addition ...
         deallocate(addlink%explink)
         deallocate(addlink)
      elseif(addlink%type.eq.7) then
!>>>>> 12A: delete volume addition ...
         deallocate(addlink)
      else
         write(*,*)'3E Cannot delete unknown addition type ',addlink%type
      endif
      addlink=>nextadd
   enddo addition
!   write(*,*)'3E phase location: ',lokph,size(phlista(lokph)%nooffr),&
!        size(phlista(lokph)%constitlist)
!   if(lokph.ne.0) then
! problem with phases, cannot deallocate these arrays, why??
!      deallocate(phlista(lokph)%nooffr)
!      deallocate(phlista(lokph)%constitlist)
!   endif
   phlista(lokph)%noofcs=0
   phlista(lokph)%nooffs=0
!   write(*,*)'all done'
1000 continue
! remove valgrind leak
   deallocate(stack)
   return
 end subroutine delphase

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable logical function iskeyword
!\begin{verbatim}
 logical function iskeyword(text,keyword,nextc)
! compare a text with a given keyword. Abbreviations allowed
! but the keyword and abbreviation must be surrounded by spaces
! nextc set to space character in text after the (abbreviated) keyword
   implicit none
   character text*(*),keyword*(*),key*64
   integer nextc
!\end{verbatim} %+
   character word*64
   logical ok
   integer kl,ks,kt
! extract the first word of text
   ks=1
   if(eolch(text,ks)) then
! if empty line, just exit
      ok=.false.; goto 1000
   else
! find the space after the first word
      kt=ks+index(text(ks:),' ')-1
! the abbreviation of the keyword must be at least 3 character !!!
      if(kt-ks.lt.3 .or. kt-ks.ge.64) then
         ok=.false.; goto 1000
      endif
   endif
   word=text(ks:kt)
   kt=kt-ks
   key=keyword
   kl=len_trim(key)
! check if word is an abbreviation of key
   if(word(1:kt).eq.key(1:kt)) then
! found keyword at start of line, set nextc to be positioned at the final space
      nextc=ks+kt
      ok=.true.
   else
      ok=.false.
   endif
!   write(*,100)ok,text(1:15),word(1:15),key(1:15),nextc,ks,kt,kl
!100 format('iskeyword: ',l1,' >',a,'<>',a,'<>',a,'<',5i3)
1000 continue
   iskeyword=ok
   return
 end function iskeyword

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable integer function istdbkeyword
!\begin{verbatim} %-
 integer function istdbkeyword(text,nextc)
! compare a text with a given keyword. Abbreviations allowed (not within _)
! but the keyword and abbreviation must be surrounded by spaces
! nextc set to space character in text after the (abbreviated) keyword
   implicit none
   character text*(*)
   integer nextc
!\end{verbatim} %+
! only those currently implemented ... rest ignored
   integer, parameter :: kwl=20
   integer, parameter :: nkw=14
   character (len=kwl), dimension(nkw), parameter :: keyword=&
        ['ELEMENT             ','SPECIES             ',&
         'PHASE               ','CONSTITUENT         ',&
         'FUNCTION            ','PARAMETER           ',&
         'TYPE_DEFINITION     ','LIST_OF_REFERENCES  ',&
         'ADD_REFERENCES      ','ASSESSED_SYSTEMS    ',&
         'DATABASE_INFORMATION','VERSION             ',&
         'DEFAULT_COMMAND     ','DEFINE              ']
!   
   character word*64
   integer j,ks,kt
! extract the first word of text
   ks=1
   if(eolch(text,ks)) then
! if empty line, just exit
      j=0; goto 1000
   else
! find the space after the first word
      kt=ks+index(text(ks:),' ')-1
! the abbreviation of the keyword must be at least 3 character, max kwl
      if(kt-ks.lt.3 .or. kt-ks.ge.kwl) then
!         write(*,*)'3E too long keyword: ',trim(text),kt-ks,kwl
         j=0; goto 1000
      endif
   endif
   word=text(ks:kt)
   kt=kt-ks
   call capson(word)
! replace - by _
90 continue
   j=index(word,'-')
   if(j.gt.0) then
      word(j:j)='_'
      goto 90
   endif
! check if word is an abbreviation of a keyword
!   write(*,*)'abbreviation: ',kt,'>',word(1:kt),'<'
!   do j=1,10
   do j=1,nkw
      if(word(1:kt).eq.keyword(j)(1:kt)) goto 100
   enddo
   j=0
!   write(*,99)j,nextc,text(1:nextc),trim(text)
99 format('3E Not a keyword: ',2i3,'>',a,'<'/1x,a)
   goto 1000
! found keyword at start of line, set nextc to be positioned at the final space
100 continue
   if(j.eq.11 .and. kt.lt.8) then
! we found 'DATA' at the start of several lines that is not DATABASE_INFO
!      write(*,*)'3E why? ',trim(text),kt
      j=0
      goto 1000
   endif
   nextc=ks+kt
!   write(*,101)j,nextc,text(1:nextc),trim(text)
101 format('3E Found keyword: ',2i3,'>',a,'<'/1x,a)
1000 continue
   istdbkeyword=j
   return
 end function istdbkeyword

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable integer function ispdbkeyword
!\begin{verbatim} %-
 integer function ispdbkeyword(text,nextc)
! compare a text with a given keyword. Abbreviations allowed (not within _)
! but the keyword and abbreviation must be surrounded by spaces
! nextc set to space character in text after the (abbreviated) keyword
   implicit none
   character text*(*)
   integer nextc
!\end{verbatim} %+
! only those currently implemented ... rest ignored
   integer, parameter :: kwl=20
   integer, parameter :: nkw=16
   character (len=kwl), dimension(nkw), parameter :: keyword=&
        ['ELEMENT             ','SPECIES             ',&
         'PHASE               ','CONSTITUENT         ',&
         'FUNCTION            ','PARAMETER           ',&
         'SPECIAL             ','BIBLIOGRAPHY        ',&
         'TABLE_OF_MODELS     ','TABLE_OF_IDENTIFIERS',&
         'DATABASE_INFORMATION','ASSESSED_SYSTEMS    ',&
         'DEFAULTS            ','VERSION             ',&
         'INCLUDE_FILE        ','CHECKSUM            ']
!   
   character word*64
   integer j,ks,kt
! extract the first word of text
   ks=1
   if(eolch(text,ks)) then
! if empty line, just exit
      j=0; goto 1000
   else
! find the space after the first word
      kt=ks+index(text(ks:),' ')-1
! the abbreviation of the keyword must be at least 3 character, max kwl
      if(kt-ks.lt.3 .or. kt-ks.ge.kwl) then
         j=0; goto 1000
      endif
   endif
   word=text(ks:kt)
   kt=kt-ks
   call capson(word)
90 continue
   j=index(word,'-')
   if(j.gt.0) then
      word(j:j)='_'
      goto 90
   endif
! check if word is an abbreviation of a keyword
!   write(*,*)'abbreviation: ',kt,'>',word(1:kt),'<'
!   do j=1,10
   do j=1,nkw
      if(word(1:kt).eq.keyword(j)(1:kt)) goto 100
   enddo
   j=0
   goto 1000
! found keyword at start of line, set nextc to be positioned at the final space
100 continue
   nextc=ks+kt
!   write(*,101)j,nextc,text(1:nextc)
!101 format('Found keyword: ',2i3,'>',a,'<')
1000 continue
   ispdbkeyword=j
   return
 end function ispdbkeyword

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine replacetab
!\begin{verbatim} %-
 subroutine replacetab(line,nl)
! replaces TAB by space in line
   implicit none
   character line*(*)
   integer nl
!\end{verbatim}
   integer ip
100 continue
   ip=index(line,char(9))
   if(ip.gt.0) then
      line(ip:ip)=' '
!      write(*,*)'Replaced TAB by space on line ',nl
      goto 100
   endif
1000 continue
   return
 end subroutine replacetab

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine readtdb
!\begin{verbatim}
 subroutine readtdb(filename,nel,selel)
! reading data from a TDB file with selection of elements, read_tdb
!-------------------------------------------------------
! Not all TYPE_DEFS implemented
!-------------------------------------------------------
   implicit none
   integer nel
   character filename*(*),selel(*)*2
!\end{verbatim} %+
   integer, parameter :: maxrejph=30,maxorddis=10,maxtypedefs=40
   character line*128,elsym*2,name1*24,name2*24,elsyms(10)*2
!   character longline*10000,reftext*512
! to read references in MatCalc TDB files
   character longline*40000,reftext*512
   character phtype*1,ch1*1,const(maxsp)*24,name3*24,funname*60,name4*60
   character refx*16
   character (len=1), dimension(maxtypedefs) :: typedefchar
   integer, dimension(maxtypedefs) :: typedefaction
   integer, dimension(5) :: addphasetypedef
   double precision mass,h298,s298
   integer, dimension(10) :: knr,endm
! lint(1,*) is sublattice, lint(2,*) is species
   double precision stoik(10),xsl,xxx
   integer lint(2,3),TDthisphase,nytypedef,nextc,keyw,tdbv,rewindx
   integer typty,fractyp,lp1,lp2,ix,jph,kkk,lcs,nint,noelx,idum,jdum
   logical onlyfun,nophase,ionliq,notent,mqmqa
   integer norew,newfun,nfail,nooftypedefs,nl,ipp,jp,jss,lrot,ip,jt
   integer nsl,ll,kp,nr,nrr,mode,lokph,lokcs,km,nrefs,ideg,iph,ics,ndisph
! disparttc and dispartph to handle phases with disordered parts
   integer nofunent,disparttc,dodis,jl,nd1,thisdis,cbug,nphrej,never,always
   character*24 dispartph(maxorddis),ordpartph(maxorddis),phreject(maxrejph)*24
!   character*24 disph(20)
   integer orddistyp(maxorddis),suck,notusedpar,totalpar,reason,zz,dismag
   integer enteredpar,loop,emodel
   type(gtp_phase_add), pointer :: addrec
   logical warning,dbcheck
! set to TRUE if element present in database
   logical, allocatable :: present(:)
! to prevent any output
   logical silent,thisphaserejected
!  mmyfr noofph
! if warning is true at the end pause before listing bibliography
!   write(*,*)'3E readtdb',allocated(seltdbph),nselph
   emodel=0
   nsl=0
   dbcheck=.FALSE.
   warning=.FALSE.
   silent=.FALSE.
   nphrej=0
   nytypedef=0
   totalpar=0
   notusedpar=0
   enteredpar=0
! this counts number of undefined/unused model-parameter-identifiers
   nundefmpi=0
   if(btest(globaldata%status,GSSILENT)) then
      silent=.TRUE.
!      write(*,*)'3E reading database silent'
   endif
   write(*,*)'3E reading a TDB file'
   if(ocv()) write(*,*)'3E reading a TDB file'
   if(.not.(index(filename,'.tdb').gt.0 &
       .or. index(filename,'.TDB').gt.0)) then
! no extention provided
      filename(len_trim(filename)+1:)='.TDB'
   endif
   if(nel.gt.0) then
      allocate(present(nel))
      present=.FALSE.
   endif
! disparttc counts the number of disordered phases to read, the
! disordered phase names are in dispartph(1..disparttc)
! dodis is nonzero only when reading the disordered part of phases.
   disparttc=0
   dodis=0
   open(21,file=filename,access='sequential',form='formatted',&
        err=1010,iostat=gx%bmperr,status='old')
! read whole TDB file to extract TYPE_DEFS with DIS_PART so disordered parts
! are not entered
!   call any_disordered_part(21,ndisph,disph)
   call any_disordered_part(21,ndisph,dispartph,ordpartph,orddistyp)
   if(ndisph.gt.0) then
!      write(*,*)'3E ndisph: ',ndisph
!      write(*,11)(trim(ordpartph(ip)),trim(dispartph(ip)),orddistyp(ip),&
!           ip=1,ndisph)
11    format('3E ord/dis: "',a,'"+"',a,'" ',i2)
   endif
   onlyfun=.FALSE.
   tdbv=1
   norew=0
   newfun=0
   nfail=0
   nrefs=0
! always is a dummy variable
   always=0
   nooftypedefs=0
! nophase set false after reading a PHASE keyword, 
! expecting next keyword to be CONSTITUENT
   nophase=.TRUE.
   rewindx=0
! return here after rewind
90  continue
   nl=0
! return here to look for a new keyword, end-of-file OK here
100 continue
   read(21,110,end=2000)line
110 format(a)
   nl=nl+1
   if(len_trim(line).gt.80) then
      if(.not.silent) write(*,121)nl
121   format(' *** Warning: line ',i5,' has characters beyond position 80,'/&
           'some information may be lost')
   endif
! One should remove TAB characters !! ?? YES !!
!   if(line(1:1).eq.'$') goto 100
   ipp=1
   if(eolch(line,ipp)) goto 100
   if(line(ipp:ipp).eq.'$') goto 100
! replace TAB by space
   call replacetab(line,nl)
!   goto 120
!---------------------------------------------------------
! handle all TDB keywords except function
120 continue
   keyw=istdbkeyword(line,nextc)
   if(.not.onlyfun) then
!      write(*,71)'3E back from istdbkeyword',keyw
      if(keyw.eq.0) then
         if(trim(line).eq.' DEFINE_SYSTEM_DEFAULT ELEMENT 2 !') then
            goto 100
         elseif(trim(line).eq.'DEFINE_SYSTEM_DEFAULT ELEMENT 2 !') then
            goto 100
         elseif(dodis.ne.1) then
! do not give this warning when reading disordered phases ...
            write(*,122)trim(line)
122         format(/' *** Warning, ignoring line: "',a,'"'/)
         endif
      endif
   endif
   if(keyw.eq.0) then
      ip=1
      if(.not.eolch(line,ip)) then
         if(ocv()) write(*,*)'3E Ignoring line: ',nl,ip,trim(line)
      endif
      goto 100
   elseif(onlyfun) then
      if(keyw.eq.5) goto 800
      goto 100
   endif
   if(.not.nophase .and. keyw.ne.4) then
! after a PHASE keyword one should have a CONSTITUENT
      if(.not.silent) write(kou,*)'WARNING expeciting CONSTITUENT: ',line(1:30)
      warning=.TRUE.
   endif
! check there is a ! in line, otherwise read until we find an exclamation mark
   ip=1
   longline(ip:)=line
   ip=len_trim(longline)+1
!   write(*,71)'3E line 1 ',ip,trim(longline)
!   write(*,*)'3E new keyword ',ip,'>',longline(1:40)
   do while(index(longline,'!').le.0)
      read(21,110,err=2200)line
      nl=nl+1
      if(line(1:1).ne.'$') then
         call replacetab(line,nl)
         longline(ip:)=line
         ip=len_trim(longline)+1
         if(ip.ge.len(longline)-100) then
            if(.not.silent) write(kou,69)nl,ip,longline(1:72)
69          format('Overflow in longline ',2i8,' for line starting:'/a)
            gx%bmperr=4304; goto 1000
         endif
      endif
   enddo
!   write(*,71)'3E line 2 ',ip,trim(longline)
71 format(a,i4,1x,a)
   if(dodis.eq.1) then
! if dodis=1 only read data for disordred phases
! PHASE=3, CONSTITUENT=4, PARAMETER=6 ... BIBLIOGRAPHIC REFERENCES=8,9
!      if(keyw.lt.3 .or. keyw.eq.5 .or. keyw.gt.6) goto 100
      if(.not.(keyw.eq.3 .or. keyw.eq.4 .or. keyw.eq.6 &
           .or. keyw.eq.8 .or. keyw.eq.9)) goto 100
   endif
!
! we have 13 keywords
  select case(keyw)
   case default
      if(ocv()) write(*,*)'3E default case: ',keyw,line(1:30)
!---------------------------------------------------------------------
!101 format('readtdb 1: ',i3,'>',a,'<')
!   if(line(2:9).eq.'ELEMENT ') then
   case(1) !element ------------------------------------------------
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! ELEMENT CR   BCC_A2                    5.1996E+01  4.0500E+03  2.3560E+01!
      ip=nextc
      if(eolch(longline,ip)) then
         if(.not.silent) &
              write(kou,*)'No element name after ELEMENT keyword on line ',nl
         gx%bmperr=4305; goto 1000
      endif
      elsym=longline(ip:ip+1)
      if(elsym.eq.'/-' .or. elsym.eq.'VA') goto 100
! allow lower case in TDB file ...
      call capson(elsym)
      if(nel.gt.0) then
! check if element among selected, if nel=0 accept all
         do jt=1,nel
            if(elsym.eq.selel(jt)) goto 76
         enddo
! ignore this element as not selected
         if(ocv()) write(*,*)'3E Skipping database element: ',elsym
!         write(*,*)'Skipping database element: ',elsym
!         write(*,*)'Select: ',nel,(selel(jt),jt=1,nel)
         goto 100
      endif
! mark we found a selected element
76    continue
      if(allocated(present)) then
         present(jt)=.TRUE.
      endif
! we seem to miss the first letter of the reference state below ??
      ip=ip+len_trim(elsym)-1
      if(eolch(longline,ip)) then
         name1='DUMMY'
         mass=one
         h298=zero
         s298=zero
      else
! extract the reference phase, third argument is 1 meaning until next space
! ix is the length of the reference phase (irrelevant here)
! ip is updated to character after the name extracted
         call getext(longline,ip,1,name1,' ',ix)
!         write(*,*)'3E longline: ',ip,longline(1:ip+10)
!         write(*,*)'3E element ref: ',name1
!         name1=longline(ip:)
!         ip=ip+len_trim(name1)
! after the name should be mass, H298-H0 and S298, ignore errors
         call getrel(longline,ip,mass)
         if(buperr.ne.0) then
            mass=one; buperr=0
         endif
         call getrel(longline,ip,h298)
         if(buperr.ne.0) then
            h298=zero; buperr=0
         endif
         call getrel(longline,ip,s298)
         if(buperr.ne.0) then
            s298=zero; buperr=0
         endif
         name2=elsym
      endif
      call store_element(elsym,name2,name1,mass,h298,s298)
      if(gx%bmperr.ne.0) goto 1000
   case(2) !SPECIES -------------------------------------------------
!   elseif(line(2:9).eq.'SPECIES ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! SPECIES O3PU2                       O3PU2!
      ip=nextc
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'Line after SPECIES keyword empty'
         gx%bmperr=4306; goto 1000
      endif
      name1=longline(ip:)
! find first space after non-space
      jp=index(name1,' ')
      name1(jp:)=' '
      ip=ip+jp
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'WARNING No stoichiometry for species: ',&
              trim(name1)
         warning=.TRUE.
         goto 100
      endif
      name2=longline(ip:)
      jp=index(name2,' ')
      name2(jp:)=' '
      call decode_stoik(name2,noelx,elsyms,stoik)
      if(gx%bmperr.ne.0) goto 1000
! check elements exist
      call enter_species(name1,noelx,elsyms,stoik)
!      write(*,*)'3E: entering species error: ',gx%bmperr
      if(gx%bmperr.ne.0) then
! if element not selected just skip the species
         if(gx%bmperr.eq.4046) then
            gx%bmperr=0; goto 100
         else
            if(.not.silent) write(kou,*)'Error entering species: ',name1,name2
            goto 1000
         endif
      endif
!-----------------------------------------------------------------------
   case(5) ! function
!   elseif(line(2:10).eq.'FUNCTION ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! FUNCTION GHSERCR    2.98150E+02  -8856.94+157.48*T-26.908*T*LN(T)
!      name1=line(11:18)
!      longline=' '
!      longline=line(20:)
!300    continue
!      jp=len_trim(longline)
!      if(longline(jp:jp).eq.'!') then
!          write(*,*)'3E Skipping function: ',name1
! all functions entered at the end, skip until !
!      do while(index(longline,'!').le.0)
      if(index(longline,'!').le.0) then
         if(.not.silent) &
              write(*,*)'3E Error, terminating ! not found for function!!',nl
         gx%bmperr=4307; goto 1000
      endif
!-------------------------------------------------------------------------
!   elseif(line(2:7).eq.'PHASE ') then
   case(3) ! PHASE
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! PHASE LIQUID:L %  1  1.0  !
      if(nophase) then
         nophase=.false.
! give a warning if any selected element is not present
         if(allocated(present)) then
            funname=' '
            kkk=1
            do jt=1,nel
               if(.not.present(jt)) then
                  funname(kkk:)=selel(jt)
                  kkk=len_trim(funname)+2
               endif
            enddo
            if(kkk.gt.1) then
               if(.not.silent) write(kou,68)funname(1:kkk)
68             format(/' *** Warning, elements not present in database: ',a/)
            endif
            deallocate(present)
         endif
      else
         if(.not.silent) write(kou,*) &
              'Error, a PHASE keyword must be followed by its CONSTIT'
         gx%bmperr=4308; goto 1000
      endif
! number of TYP_DEFS for this phase
      TDthisphase=0
      ip=nextc
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'line after PHASE empty'
         goto 100
      endif
      name1=longline(ip:)
! convert phase name to upp case
      call capson(name1)
      jp=index(name1,' ')
      ip=nextc+jp
      if(jp.gt.0) then
         name1(jp:)=' '
      endif
      jp=index(name1,':')
!      write(*,*)'3E readtdb 11: ',name1,ip,jp
! phytype, a letter after the phase name separated by a :, for example GAS:G
      if(jp.gt.0) then
         phtype=name1(jp+1:jp+1)
         name1(jp:)=' '
      else
         phtype=' '
      endif
! we must know if we have the mqmqa model before reading constituents!!
! tested below also.
      if(phtype.eq.'Q') then
         mqmqa=.TRUE.
      else
         mqmqa=.FALSE.
      endif
! check if phase rejected
!      write(*,*)'3E number of phases rejected: ',nphrej
      do jt=1,nphrej
         if(name1.eq.phreject(jt)) then
            thisphaserejected=.TRUE.
!            write(*,*)'3E skipping rejected phase: ',name1
! why is nophase set true? If I comment it away nothing read!!
            nophase=.true.
            goto 100
         endif
      enddo
! SELECTED_PHASES
! if seltdbph is allocated check if this phase selected when dodis=0
      if(dodis.eq.0 .and. allocated(seltdbph) .and. nselph.gt.0) then
!         write(*,*)'3E Calling isabbr: ',trim(name1),nselph
         jt=isabbr(name1,seltdbph,nselph)
!         write(*,*)'3E return from isabbr: ',jt
         if(jt.eq.0) then
!            write(*,*)'3E not selected phase: ',trim(name1)
            thisphaserejected=.TRUE.
            nophase=.true.
            goto 100
         else
            write(*,'(a,a,a,a)')'3E Phase ',trim(name1),' fits selection ',&
                 trim(seltdbph(jt))
         endif
      endif
! end elected phases
      thisphaserejected=.false.
!      write(*,*)'3E nophase set to false, phase: ',name1
      ip=ip+1
      jp=ip
      name2=longline(ip:jp)
      thisdis=0
      phdis: if(dodis.eq.1) then
! special when reading disordered parts, check phase name equail
!         write(*,*)'3E Check if disordered part: ',dodis,name1
         do jt=1,disparttc
            if(name1.eq.dispartph(jt)) goto 307
         enddo
! not a disordered part
         goto 100
307      continue
         thisdis=jt
         write(*,'(3a)')'3E ',trim(name1),' is a disordered part'
! check if disordered phase is magnetic!, we have to step though type_defs
         dismag=0
         if(.not.eolch(longline,jp)) then
            ch1=longline(jp:jp)
!            write(*,312)trim(longline),jp,trim(longline(jp+1:))
312         format('3E distypes: ',a,i3,' "',a,'"')
            dmag: do while(ch1.ne.' ')
               do jt=1,nooftypedefs
                  if(ch1.eq.typedefchar(jt) .and. &
                       ((typedefaction(jt).eq.-1 .or. &
                          typedefaction(jt).eq.-3))) then
                     dismag=typedefaction(jt)
!                     write(*,*)'3E disordered part is magnetic',dismag
                     exit dmag
                  endif
               enddo
               jp=jp+1
               ch1=longline(jp:jp)
            enddo dmag
         endif
! we skip the rest of the phase line ...
         goto 100
      elseif(dodis.eq.0 .and. ndisph.gt.0) then
! make use of initial read of TDB file to skip phases that are disordered parts
!      write(*,*)'3E comparing "',trim(name1),'" with "',trim(disph(1)),'" etc'
         do jt=1,ndisph
!            if(name1.eq.disph(jt)) then
            if(name1.eq.dispartph(jt)) then
!               write(*,*)'3E Phase ',trim(name1),' is a disordered part of ',&
!                    trim(ordpartph(jt)),jt,nphrej
! if the phase ordpartph(jt) is rejected, enter the disordered phase!!
               do zz=1,nphrej
!                  write(*,*)'3E check "',trim(ordpartph(jt)),'" and "',&
!                       trim(phreject(zz)),'"'
                  if(ordpartph(jt).eq.phreject(zz)) then
                     write(*,'(a,a,a,a,a)')'3E Keeping ',trim(name1),&
                          ' because phase ',trim(phreject(zz)),' is rejected'
                     goto 310
                  endif
               enddo
! do not enter this phase as it is a disordered part
! all these must be set ...
               thisdis=-1
               nophase=.true.
               thisphaserejected=.TRUE.
               goto 100
            endif
         enddo
      elseif(dodis.eq.0 .and. disparttc.gt.0) then
! we must not enter phases that are disordered parts
         do jt=1,disparttc
            if(name1.eq.dispartph(jt)) then
!               write(*,*)'3E Skip phase that is a disordered part: ',name1
               thisdis=-1
               goto 100
            endif
         enddo
      endif phdis
!      write(*,*)'3E Entering phase: ',name1
!      write(*,*)'3E Checking phase types for phase: ',name1,jp
! skip blanks, then read type code, finished by a blank
      if(eolch(longline,jp)) then
         if(.not.silent) &
              write(kou,*)'3E WARNING no phase typecode: ',trim(name1)
         warning=.TRUE.
      endif
      jp=jp-1
! WE MUST CHECK IF TYPE_DEFS appear after phases have been entered!!
!      write(*,311)'3E TDs: ',nooftypedefs,&
!           (typedefchar(jt),jt=1,nooftypedefs)
! return here to check for different TYPE_DEFS
310   jp=jp+1
! check which type_defs that has been entered
!      write(*,*)'3E typedefs: ',trim(name1),': ',trim(longline(jp:)),jp
! NOTE and FIX: type code expected to be after a single space: be flexible ??
      typedefcheck: if(longline(jp:jp).ne.' ') then
         ch1=longline(jp:jp)
         if(always.eq.3) then
! this code an attempt to fool -O2 compiler switch
!            write(*,*)'3E typedef for ',trim(name1),': ',ch1,TDthisphase
!            write(*,311)'3E TDs: ',nooftypedefs,&
!                 (typedefchar(jt),jt=1,nooftypedefs)
311      format(a,i3,': ',10('"',a,'", '))
            always=always+1
         endif
         do jt=1,nooftypedefs
            if(ch1.eq.typedefchar(jt)) goto 320
         enddo
! ignore typedef % meaning sequantial read ...
         if(ch1.eq.'%') goto 310
! WARNING that unknown TYPE_DEF has been used!!
         write(kou,313)trim(name1),ch1
313      format(' *** WARNING: phase ',a,' has unknown TYPE_DEF: ',a/&
              ' *** Move all TYPE_DEFS before used in any phase!')
         warning=.TRUE.
         goto 310
320      continue
         if(typedefaction(jt).eq.99) then
! ignore TYPE_DEF SEQ
            continue
         elseif(typedefaction(jt).eq.-1 .or. &
              typedefaction(jt).eq.-3) then
! Inden magnetic addition, save for after phase created
            TDthisphase=TDthisphase+1
            addphasetypedef(TDthisphase)=typedefaction(jt)
         elseif(typedefaction(jt).ge.25 .and. &
              typedefaction(jt).le.37) then
! Qing-Xiong magnetic addition
            TDthisphase=TDthisphase+1
            addphasetypedef(TDthisphase)=typedefaction(jt)
         elseif(typedefaction(jt).eq.1905) then
! Einstein
            TDthisphase=TDthisphase+1
            addphasetypedef(TDthisphase)=typedefaction(jt)
         elseif(.not.(typedefaction(jt).eq.100.or.typedefaction(jt).eq.0)) then
! give an alert if typedefaction is not 100
            write(*,*)'3E Unknown typedefaction: ',typedefaction(jt)
         endif
         goto 310
      endif typedefcheck
!      write(*,*)'3E typedefs for ',trim(name1),': ',TDthisphase,&
!           (addphasetypedef(ll),ll=1,TDthisphase)
      name2='TDB file model: '//name2
! number of sublattices
!      write(*,*)'3E buperr: ',buperr ,jp
      call getrel(longline,jp,xsl)
      if(buperr.ne.0) then
         if(.not.silent) write(kou,*)'3E tdb: "',longline(1:jp),'"',buperr
         gx%bmperr=buperr; goto 1000
      endif
! dummy statement to fool -O2 optimization
      if(nsl.lt.0) jt=1
      nsl=int(xsl)
      do ll=1,nsl
         call getrel(longline,jp,stoik(ll))
         if(buperr.ne.0) then
            gx%bmperr=buperr; goto 1000
         endif
      enddo
!      write(*,*)'3E readtdb 3A: ',nsl,(stoik(ll),ll=1,nsl)
!---------------------------------------------------------------------
! The constituent line must follow PHASE before any new phase
   case(4) !    CONSTITUENT LIQUID:L :CR,FE,MO :  !
! the phase must have been defined
      if(nophase) then
         if(thisphaserejected) goto 100
         if(.not.silent) write(kou,327)trim(longline)
327      format('A CONSTITUENT keyword not directly preceeded by PHASE!'/a)
         gx%bmperr=4308; goto 1000
      endif
      nophase=.true.
      condis1: if(dodis.eq.1) then
! searchin why sigma in TAFID does not have c disordered fraction set
!         write(*,*)'3E sigma 17:',trim(longline),thisdis
!         write(*,*)'3E sigma 17:',thisdis
         if(thisdis.eq.0) goto 100
! we skip the constituent line and go directly to create disordered fractions
         goto 395
      elseif(disparttc.gt.0 .and. thisdis.lt.0) then
! this is a disordered part, skip
         goto 100
      endif condis1
!360    continue
      jp=len_trim(longline)
!      write(*,*)'3E readtdb gas1: ',nl,jp,longline(1:jp)
! eliminate all after the exclamation mark
!      longline(jp+1:)=' '
! 
      ip=index(longline,' :')+2
      if(mqmqa) then
! this is a FactSage MQMQA model for liquids
! entering constituents as quadrupoles
!         write(*,'(a,a,a,2i5)')'3E mqmqa const: "',trim(longline(ip:jp)),&
!              '"',ip,jp
         loop=0
! MQMQA constituents created "on the fly" as quadrupols using existing species
! and additional coordination numbers n1..n4. A  / separate sublattices
! a , separate species in same sublattice. If any A B X Y species not entered
! the quadrupole is ignored (not an error)
! A/X n1 n2 A,B/X n1 n2 n3 B/X,Y n1 n2 n3 A,B/X,Y n1 n2 n3 n4 ...
         call mqmqa_constituents(longline(ip:jp),const,loop)
!         write(*,*)'3E back from entering constituents',gx%bmperr
         if(gx%bmperr.ne.0) then
            write(*,*)'3E error entering quadrupoles'
            goto 1000
         endif
         call mqmqa_rearrange(const)
!         write(*,*)'3E back from rearranging constituents',gx%bmperr
         if(gx%bmperr.ne.0) then
            write(*,*)'3E error rearranging quadrupoles'
            goto 1000
         endif
! skip the rest below except entering the phase
! stoik(1) is bonds/atom, just for output, never used explicitly
         stoik(1)=2.0D0
         knr(1)=mqmqa_data%nconst
!         write(*,*)'3E enter_p: ',trim(name1),' ',knr(1),stoik(1),' ',phtype
         name2='MQMQA '
         call enter_phase(name1,1,knr,const,stoik,name2,phtype,warning,emodel)
!         write(*,*)'3E back from entering phase',gx%bmperr
         if(gx%bmperr.ne.0) then
            write(*,*)'3E failed to enter the MQMQA phase',gx%bmperr
         endif
         goto 100
      endif
!      write(*,*)'3E readtdb gas2: ',jp,longline(1:jp)
      ll=0
      nr=0
      nrr=0
!      write(*,*)'3E readtdb 3C: ',ll,nr,nsl,longline(ip:jp)
! mode=1 indicates to getname that / + - are allowed in species names
      mode=1
370   continue
      if(ll.ge.1) then
         knr(ll)=nr
         if(nr.le.0) then
            if(ocv()) then
               write(*,*)'3E Skipping phase due to missing constituents: ',name1
!              write(*,378)name1,ll
378            format('Phase ',a,' has no constituents in sublattice ',i2)
! Not a fatal error when elements have been selected but skip this phase
            endif
            goto 100
         endif
      endif
      ll=ll+1
!      write(*,*)'3E start sublat ',ll,nsl,nr,ip
      if(ll.gt.nsl) goto 390
      nr=0
380   continue
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'Error extracting constituents 1'
         gx%bmperr=4309; goto 1000
      endif
      nr=nr+1
      nrr=nrr+1
!      write(*,379)'readtdb 3CXX: ',ip,nr,longline(ip:ip+10)
379   format(a,2i4,' >',a,'< >',a,'< >',a,'<')
      call getname(longline,ip,name3,mode,ch1)
!      write(*,379)'readtdb 3CY: ',ip,nr,longline(ip:ip+10),name3,ch1
      if(buperr.ne.0) then
!         write(*,381)'readtdb 3E: ',ll,nr,longline(1:ip+5),ip,name3
381      format(a,2i4,' "',a,'" ',i5,1x,a,'"',a)
         gx%bmperr=buperr; goto 1000
      endif
!      write(*,381)'readtdb 3E: ',ll,nr,longline(1:ip+5),ip,name3,ch1
      const(nrr)=name3
! bypass any "major" indicator %
      if(ch1.eq.'%') ip=ip+1
      if(eolch(longline,ip)) then
!         if(.not.silent) write(kou,*)'Error extracting constituents 2'
         gx%bmperr=4309; goto 1000
      endif
! check that const(nrr) among the selected elements ...
!      write(*,*)'Testing constituent: ',name3,nr
      call find_species_record_noabbr(name3,lp1)
      if(gx%bmperr.ne.0) then
! this species is not present, not a fatal error, skip it and continue
!         write(*,*)'Skipping constituent: ',name3
         gx%bmperr=0; nrr=nrr-1; nr=nr-1
      endif
      ch1=longline(ip:ip)
      if(ch1.eq.',') then
         ip=ip+1; goto 380
      elseif(ch1.eq.':') then
         ip=ip+1; goto 370
      endif
      if(ch1.ne.'!') goto 380
! when an ! found the list of constutents is finished.  But we
! should have found a : before the !
      if(.not.silent) write(kou,*)'Found "!" before terminating ":"'
      gx%bmperr=4310; goto 1000
!      write(*,*)'Species terminator error: ',ch1,nl
!      gx%bmperr=4157; goto 1000
390    continue
! name2 is model, ignored on reading TDB
      ionliq=.FALSE.
      mqmqa=.FALSE.
      if(phtype.eq.'Y') then
         name2='IONIC_LIQUID '
         ionliq=.TRUE.
      elseif(phtype.eq.'Q') then
         name2='MQMQA '
         mqmqa=.TRUE.
      else
         name2='CEF-TDB-RKM? '
      endif
      if(ocv()) write(*,*)'readtdb 9: ',name1,nsl,knr(1),knr(2),phtype
395   continue
!
! THE CODE HERE IS A MESS .... NEW PDB FORMAT NEEDED
!
!    write(*,*)'3E sigma4 label 395 add disordered fraction set: ',dodis,nphrej
      condis2: if(dodis.eq.1) then
! if we have a disordered part do not enter the phase, add disordered fracs!
! the ordered phase name is ordpart(thisdis)
!         write(*,*)'3E sigma19: ',trim(ordpartph(disparttc)),disparttc,&
!              trim(ordpartph(thisdis)),thisdis
         do jt=1,nphrej
!            if(ordpartph(disparttc).eq.phreject(jt)) then
! why disparttc?
            if(ordpartph(thisdis).eq.phreject(jt)) then
               write(*,'(a,a,a)')'3E ordered part ',trim(phreject(jt)),&
                    ' is rejected, keep disordered part '
               goto 100
            endif
         enddo
!         write(*,*)'3E sigma20: ',trim(ordpartph(thisdis))
         call find_phase_by_name(ordpartph(thisdis),iph,ics)
         if(gx%bmperr.ne.0) then
! NOTE THE ORDERED PHASE MAY NOT BE ENTERED DUE TO COMPONENTS!!
            if(.not.silent) write(kou,396)trim(ordpartph(thisdis))
396         format('3E and disordered part ',a,' has not been selected')
            warning=.TRUE.
            gx%bmperr=0
            goto 100
         else
!            if(.not.silent) write(kou,*) &
            write(kou,'(a,a,3i3)')'3E Adding disordered part to ',&
                 trim(ordpartph(thisdis)),orddistyp(thisdis),thisdis,dismag
            if(dismag.ne.0) then
! disordered phase magnetic, check if ordered is also ...
               lokph=phases(iph)
               nullify(addrec)
               addrec=>phlista(lokph)%additions
               write(*,*)'3E check if ordered phase has  magnetic model'
!   type(gtp_phase_add), pointer :: addrec
               do while(associated(addrec))
!               write(*,*)'3E addrec: ',addrec%type,INDENMAGNETIC,XIONGMAGNETIC
                  if(addrec%type.eq.INDENMAGNETIC) goto 798
                  if(addrec%type.eq.XIONGMAGNETIC) goto 798
                  addrec=>addrec%nextadd
               enddo
!               write(*,*)'3E adding magnetic model to ordered phase'
! ordered not magnetic, set the same as disordered               
               if(dismag.eq.-1) then
! Inden magnetic for BCC
                  call add_addrecord(lokph,'Y',indenmagnetic)
               elseif(dismag.eq.-3) then
! Inden magnetic for FCC/HCP
                  call add_addrecord(lokph,'N',indenmagnetic)
               endif
            endif
         endif
798      continue
! we are creating the phase, there is only one composition set, iph is ordered
!         write(*,*)'3E sigma18: get_phase_compset'
         call get_phase_compset(iph,1,lokph,lokcs)
         if(gx%bmperr.ne.0) goto 1000
! ch1 is suffix for parameters, always D
         ch1='D'
! jl=0 if NDM (sigma)
! jl=1 if phase can be totally disordered (but can have interstitials)
! nd1 is the number of sublattices to sum into disordered set
         if(orddistyp(thisdis).eq.1) then
            jl=1
            if(phlista(lokph)%noofsubl.le.5) nd1=4
            if(phlista(lokph)%noofsubl.le.3) nd1=2
!            if(.not.silent) write(kou,397) trim(ordpartph(thisdis)),nd1
            write(kou,397) trim(ordpartph(thisdis)),nd1,thisdis
397         format('3E Phase ',a,' has order/disorder partition model',&
                 ' adding first ',i2,'; thisdis: ',i2)
         else
            jl=0
            nd1=phlista(lokph)%noofsubl
         endif
!         goto 402
402      continue
         if(jl.eq.0 .and. .not.silent) write(kou,398)trim(ordpartph(thisdis))
398      format('3E The phase ',a,' cannot be completely disordered')
! add DIS_PART from TDB
!         write(*,*)'3E adding disordered fraction set',csfree,highcs
         call add_fraction_set(iph,ch1,nd1,jl)
         if(gx%bmperr.ne.0) then
            if(.not.silent) write(kou,*) &
                 '3E Error entering disordered fraction set: ',gx%bmperr
            goto 1000
         endif
!         suck= newhighcs(.true.)
!         write(*,*)'3E added disordered fraction set 1: ',csfree,highcs,suck
         if(jl.eq.0) then
! we must set the correct formula unit of the disordered phase, on the
! TDB file it is unity.  Sum up the sites for the ordered phase in lokcs
            xxx=zero
            do ll=1,nd1
               xxx=xxx+firsteq%phase_varres(lokcs)%sites(ll)
            enddo
            firsteq%phase_varres(lokcs)%disfra%fsites=xxx
         else
            xxx=one
         endif
!         if(.not.silent) write(kou,601) &
!              dispartph(thisdis)(1:len_trim(dispartph(thisdis))),ch1,nd1,jl,xxx
601      format('3E Add parameters from disordered part: ',a,5x,a,2x,2i3,F12.4)
      else
!         write(*,*)'3E enter phase: ',name1
!         call enter_phase(name1,nsl,knr,const,stoik,name2,phtype)
         call enter_phase(name1,nsl,knr,const,stoik,name2,phtype,warning,emodel)
!      write(*,*)'readtdb 9A: ',gx%bmperr
         if(gx%bmperr.ne.0) then
            if(gx%bmperr.eq.4121) then
               if(.not.silent) write(kou,*) &
                    '3E Phase ',name1(1:len_trim(name1)),&
                    ' is ambiguous or short for another phase'
            endif
            goto 1000
         endif
! any typedefs? only magnetic handelled at present
         call find_phase_by_name(name1,iph,lcs)
!         write(*,*)'readtdb 9X: ',gx%bmperr
         if(gx%bmperr.ne.0) then
            if(.not.silent) write(kou,*)'Phase ',name1,' is ambiguous'
            goto 1000
         endif
         lokph=phases(iph)
!         write(*,*)'3E typedefs for ',trim(name1),lokph,TDthisphase
         phasetypes: do jt=1,TDthisphase
!            write(*,*)'3E manage typedef ',jt,addphasetypedef(jt)
            if(addphasetypedef(jt).eq.-1) then
! Inden magnetic for BCC
               call add_addrecord(lokph,'Y',indenmagnetic)
!               call add_magrec_inden(lokph,1,-1)
            elseif(addphasetypedef(jt).eq.-3) then
! Inden magnetic for FCC and other phases
               call add_addrecord(lokph,'N',indenmagnetic)
!               call add_magrec_inden(lokph,1,-3)
            elseif(addphasetypedef(jt).eq.1905) then
! Einstein lowt model
               call add_addrecord(lokph,' ',einsteincp)
            else
! Assumed Xiong magnetic, the factor 0.37 (BCC) or 0.25 (FCC) needed
!               write(*,*)'3E Entering Qing-Xiongmagnetic ',addphasetypedef(jt)
! in TDB files ALWAYS average bohr magenton numbers
               phlista(lokph)%status1=ibset(phlista(lokph)%status1,PHBMAV)
               if(addphasetypedef(jt).eq.37) then
! BCC ........... very cryptic: space, " ", means not idividual IBM
                  call add_addrecord(lokph,'Y ',xiongmagnetic)
               elseif(addphasetypedef(jt).eq.25) then
! FCC and others
                  call add_addrecord(lokph,'N ',xiongmagnetic)
               else
                  write(*,13)lokph,addphasetypedef(jt)
13                format(78('*')/'3E unknown addition: ',2i7/78('*'))
               endif
            endif
            if(gx%bmperr.ne.0) goto 1000
         enddo phasetypes
!         write(*,607)trim(name1),iph
607      format('3E Entered phase ',a,i5)
      endif condis2
!      write(*,*)'readtdb 9B:',name1,nsl,phtype
!-------------------------------------------------------------------
   case(6) ! PARAMETER --------------------------------------------
!   elseif(line(4:13).eq.'PARAMETER ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
!   PARAMETER G(LIQUID,CR;0)  2.98150E+02  +24339.955-11.420225*T
      if(eolch(longline,nextc)) then
         if(.not.silent) write(kou,*)'Empty line after PARAMETER'
         gx%bmperr=4311; goto 1000
      endif
!      if(dodis.eq.1) write(*,*)'Reading disordered parameters'
!      write(*,*)'3E found parameter: ',totalpar,dodis,nl
! count parameter only when dodis=0
!      if(dodis.eq.1) totalpar=totalpar+1
      if(dodis.eq.0) totalpar=totalpar+1
      ip=nextc
      funname=longline(ip:)
! problem with default low T limit, can be ,, directly after parameter )
      kp=index(funname,' ')
      cbug=index(funname,'),')
! save position after parameter name in nextc
      if(cbug.gt.0 .and. cbug.lt.kp) then
         nextc=ip+cbug+1
         kp=cbug+1
!         write(*,*)'3E ,,2: ',trim(longline),ip,kp
      else
         cbug=index(funname,')')
         if(cbug.lt.kp) then
            nextc=ip+kp
         else
! We have spaces inside constituent arrays !!!
            kp=cbug+1
            nextc=ip+cbug
!            write(*,*)'3E spaces inside constituent array ',&
!                 trim(funname(1:kp)),kp,nextc
            funname(kp:)=' '
612         continue
            cbug=index(funname(1:kp),' ')
            if(cbug.gt.0) then
               funname(cbug:)=funname(cbug+1:)
               kp=kp-1
               goto 612
            endif
!            write(*,*)'3E spaces removed in constituent array? ',&
!                 trim(funname(1:kp)),kp,nextc
            kp=kp+1
         endif
      endif
      funname(kp:)=' '
! extract symbol, normally G or L but TC, BMAGN and others can occur
      lp1=index(funname,'(')
      name1=funname(1:lp1-1)
! strange bu when V0 is interpreted as LPZ !!!
!      write(*,*)'3E mpi: ',name1
      typty=0
! this "L " is kept for compatibility with old TDB files
      if(name1(1:2).eq.'G ' .or. name1(1:2).eq.'L ') then
         typty=1
      elseif(name1(1:3).eq.'TC ') then
         typty=2
! NOTE this is actually too long only 4 letters mpi should be allowed
      elseif(name1(1:6).eq.'BMAGN ') then
         typty=3
!      elseif(name1(1:3).eq.'V0 ') then
! Wow ... these not corrected when changing model_parameter_id !!!
!         typty=8
!      elseif(name1(1:3).eq.'VA ') then
!         typty=9
      endif
! we should handle also other parameter types
      if(typty.eq.0) then
! find the property associated with this symbol
!         write(*,*)'psym1: ',trim(name1)
         call get_parameter_typty(name1,lokph,typty,fractyp)
         if(gx%bmperr.ne.0) then
            lp2=len_trim(name1)
            do lp1=1,nundefmpi
               if(undefmpi(lp1)(1:lp2).eq.trim(name1)) goto 618
            enddo
            if(nundefmpi.lt.mundefmpi) then
               nundefmpi=nundefmpi+1
               undefmpi(nundefmpi)=trim(name1)
            else
               write(*,*)'3B too many model parameter identifier errors',&
                    mundefmpi
            endif
            if(.not.silent) write(kou,*) &
                 ' *** WARNING unknown parameter identifier, "',&
                 trim(name1),'" on line: ',nl
618         continue
            gx%bmperr=0; typty=0
            warning=.TRUE.
         endif
!         write(*,*)'psym2: ',typty,fractyp
      endif
! fractyp 1 is normal or ordered part if there is a disordered part
      fractyp=1
!       write(*,*)'readtdb: PAR',name1,typty
! extract phase name and constituent array
      lp1=index(funname,'(')
      lp2=index(funname,',')
      name2=funname(lp1+1:lp2-1)
      dispar: if(dodis.eq.1) then
! first check if phase name is a disordered part, if not skip
! then change phase name to ordered phase and set fractyp=2
! and add a suffix D to parameter symbol
         do jl=1,disparttc
            if(name2.eq.dispartph(jl)) goto 710
         enddo
!         notusedpar=notusedpar+1
! not disordered phase, skip this parameter
!         goto 100
         reason=1
         goto 888
!-----------------------
! This parameter was added to notusedpar at first run, correct that now
710      continue
         notusedpar=notusedpar-1
!         write(*,*)'Entering disordered parameter to: ',thisdis,jl
         thisdis=jl
         if(dbcheck) write(*,887)notusedpar,longline(ip:ip+55)
887      format('3E restored: ',i5,': ',a)
!         write(*,*)'Entering disordered parameter to: ',ordpartph(thisdis)
!         write(*,*)'3E ',longline(1:len_trim(longline))
         name2=ordpartph(jl)
         fractyp=2
      endif dispar
!---------------------- check phase is entered,
! the database may contain many phases that are not selected
!      if(name2(1:2).eq.'ZR') write(*,*)'3E parameter for phase: ',trim(name2)
      call find_phase_by_name_exact(name2,jph,kkk)
!      write(*,*)'readtdb 19: ',jph,gx%bmperr,name2
      if(gx%bmperr.ne.0) then
! Why is ZRTE not accepted?? ... exact match with first phase was not OK! suck
         if(gx%bmperr.eq.4121) &
              write(*,*)'3E WARNING parameter with ambiguous phase name',&
              ' ignored: ',trim(name2)
! this parameter is not entered as phase not entered
!         notusedpar=notusedpar+1
!         gx%bmperr=0; goto 100
         gx%bmperr=0; reason=2; goto 888
!         goto 1000
      endif
! extract constituent array, remove final ) and decode
! constituent names can be very long ....
      lokph=phases(jph)
      if(btest(phlista(lokph)%status1,PHIONLIQ)) then
! check if ionic liquid for handling neutrals ... with or without *:
         ionliq=.TRUE.
      else
         ionliq=.FALSE.
      endif
      name4=funname(lp2+1:)
! find terminating )
      lp1=index(name4,')')
!      if(name2(1:7).eq.'FCC_L12') then
!         write(*,*)'3E constituent array: ',trim(name4)
!      endif
      if(lp1.le.0) then
         if(.not.silent) then
! problem with space in constituent array ...
            write(kou,*) &
                 '3E WARNING missing ")" in parameter constituent array "',&
                 trim(name2),',',trim(name4),'", line:',nl
            write(*,*)'3E funname: ',trim(funname(lp2+1:))
            write(*,*)'3E longline: ',trim(longline)
         endif
         warning=.TRUE.
!         notusedpar=notusedpar+1
!         goto 100
         reason=3
         goto 888
      else
         name4(lp1:)=' '
      endif
! Handling of ionic liquid parameters for neutrals
      if(ionliq) then
         nsl=index(name4,':')
!         write(*,*)'3E ionic liquid parameter: ',trim(name4),nsl
         if(nsl.le.0) then
            name4(3:)=name4
            name4(1:2)='*:'
!            write(*,*)'3E Added wildcard to parameter: ',trim(name4)
         endif
      endif
297   continue
!
      call decode_constarr(lokph,name4,nsl,endm,nint,lint,ideg)
      if(ocv()) write(*,303)'readtdb 303: ',name4(1:len_trim(name4)),&
           nsl,endm(1),endm(2),nint,((lint(ip,jp),ip=1,2),jp=1,nint)
303   format(a,a,2i4,2x,2i3,' : ',3(2i3,2x))
      if(gx%bmperr.ne.0) then
! error here can mean parameter with un-selected constituent, i.e. no error
!         write(*,*)'3E: decode',ionliq,tdbv,nsl,gx%bmperr
         if(ionliq .and. tdbv.eq.1 .and. nsl.eq.1) then
! handle parameters in ionic liquids with only neutrals in second sublattice
! in TC one can have no constituent there or an arbitrary constituent,
! in OC the constituent in sublattice 1 must be a *
            nsl=2
            endm(2)=endm(1)
            endm(1)=-99
! shift any interaction from sublattice 1 to 2
            do ip=1,nint
!               write(*,*)'3E lint: ',lint(1,ip),lint(2,ip)
               lint(2,ip)=2
            enddo
            if(ocv()) write(*,303)'modif endmem: ',name4(1:len_trim(name4)),&
                 nsl,endm(1),endm(2),nint,((lint(ip,jp),ip=1,2),jp=1,nint)
            gx%bmperr=0
         else
            if(ocv()) write(*,*)'Skipping parameter: ',name4(1:len_trim(name4))
!            notusedpar=notusedpar+1
!            gx%bmperr=0; goto 100
            gx%bmperr=0; reason=4; goto 888
         endif
      endif
!      if(nint.gt.1) then
! lint(1,1) is species of first, lint(1,2) in second interaction
!          write(*,305)'readtdb 305: ',endm(1),nint,lint(2,1),lint(2,2)
!      endif
305    format(a,5i4)
!---------------- encode function
!      if(dodis.eq.1) write(*,*)'We are here 1'
      ip=0
      jp=0
400    continue
      ip=ip+1
405    continue
      ch1=funname(ip:ip)
! accept the first 8 letters and numbers of phase name
      if((ch1.ge.'A' .and. ch1.le.'Z') .or. &
         (ch1.ge.'0' .and. ch1.le.'9')) goto 400
      if(ch1.ne.' ') then
         funname(ip:)=funname(ip+1:)
         jp=jp+1
         if(jp.lt.8) goto 405
         funname(ip+1:)=' '
      endif
      funname='_'//funname
!-------------------------------------------------
! now read the function, start from position nextc
!      write(*,398)'3E ,,: ',trim(longline),nextc
      longline=longline(nextc:)
!410    continue
      jp=len_trim(longline)
      if(longline(jp:jp).ne.'!') then
         if(.not.silent) write(kou,410)nl,ip,longline(1:ip)
410      format('Error, parameter line not ending with !',2i5/a)
         gx%bmperr=4312; goto 1000
      endif
! extract bibliographic reference if any
! NOTE: a legal ending is ;,,,!
      refx='none'
      kp=jp-1
      do while(longline(kp:kp).ne.';')
         kp=kp-1
         if(kp.lt.1) then
! illegal termination of function in TDB file
            if(.not.silent) write(kou,417)nl
417 format('No final ; of function in TDB file, around line: ',i5)
            gx%bmperr=4013; goto 1000
         endif
      enddo
      kp=kp+2
! longline(kp:kp) is character after "; " or ";," 
! next is upper temperature limit or , meaning default.  We have a "!" at end
430   continue
      if(eolch(longline,kp)) continue
      if(longline(kp:kp).eq.',') then
         kp=kp+1
      elseif(longline(kp:kp).eq.'!') then
         goto 433
      else
!    ; 6000 N 91DIN !
!   kp=^                 => index(...,' ')=5; kp=kp+4
         kp=kp+index(longline(kp:),' ')-1
      endif
! next is N or ,
      if(eolch(longline,kp)) continue
      if(longline(kp:kp).ne.'!') then
         kp=kp+1
      endif
      if(eolch(longline,kp)) continue
      if(kp.lt.jp) then
         refx=longline(kp:jp-1)
         call capson(refx)
      else
         refx=' '
      endif
! ------------------- we found the reference, continue with the expression
433   continue
! replace any # by ' '
412   continue
      jss=index(longline(1:jp),'#')
      if(jss.gt.0) then
         longline(jss:jss)=' '
         goto 412
      endif
!      write(*,*)'3E Entering function 2: ',funname,trim(longline)
      lrot=0
!      call store_tpfun(funname,longline,lrot,.TRUE.)
      call store_tpfun(funname,longline,lrot,rewindx)
!          write(*,17)lokph,typty,nsl,lrot,(endm(i),i=1,nsl)
17 format('readtdb 17: ',4i3,5x,10i3)
!         write(*,404)'readtdb entpar: ',refx,fractyp,nint,ideg
404   format(a,a,i3,2x,10i3)
      if(gx%bmperr.ne.0) then
         if(.not.silent) write(kou,*)'Error set: ',gx%bmperr,lrot,' ',&
              funname(1:len_trim(funname)),' around line: ',nl
         goto 1000
      else
!         if(dodis.eq.1) write(*,*)'We are here 2'
         call enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,&
              lrot,refx)
         if(ocv()) write(*,407)'3E Entered parameter: ',lokph,typty,gx%bmperr
!         write(*,407)'Entered parameter: ',lokph,typty,gx%bmperr
         if(gx%bmperr.ne.0) then
! error entering parameter, not fatal
!            if(dodis.eq.1 .and. .not.silent) &
!                 write(*,408)'3E parameter warning:',gx%bmperr,nl,&
!                 funname(1:40)
!408         format(a,i6,' line ',i5,': ',a)
!            if(.not.(gx%bmperr.ne.4096 .or. gx%bmperr.ne.4066)) then
!               goto 1000
! Error 4096 means "no such constituent" and 4066 "... in a sublattice"
! Error 4154 means no reference but the parameter has been entered
            if(gx%bmperr.eq.4096 .or. gx%bmperr.eq.4066 .or. &
                 gx%bmperr.eq.4154) then
! this means the user has not selected this component or forgot reference
!            write(*,*)'readtdb entparerr: ',gx%bmperr,' >',&
!                 funname(1:len_trim(funname))
! error 4154 means missing reference but the parameter is entered
               if(gx%bmperr.eq.4154 .and. .not.silent) then
                  write(*,409)gx%bmperr,nl
409               format('3E Parameter reference missing ',i6,&
                       ', around line: ',i7,', continuing')
                  warning=.TRUE.
               endif
            else
! Other errors than 4096, 4066 and 4154 are fatal
               goto 1000
            endif
            gx%bmperr=0
!         else
         endif
         enteredpar=enteredpar+1
!         write(*,407)'3E Entered parameter: ',lokph,typty,gx%bmperr,enteredpar
407      format(a,4i5)
      endif
! there cannot be any error when we come here ...
!      if(gx%bmperr.ne.0 .and. .not.silent) &
!           write(*,*)'3E parameter function error: ',gx%bmperr
      goto 100
!------------------------------------------------------------------
! this is end of PARAMETER keyword
888   continue
! TAFID with 9000 parameters have about 100 unused when all selected
! reason 1= parameter not part of disordered fraction set after rewind
! reason 2= phase not entered
! reason 3= constituent array error
! reason 4= constituent array not selected
      if(reason.ne.1) then
         notusedpar=notusedpar+1
         if(dbcheck) write(*,889)reason,notusedpar,longline(ip:ip+55)
889      format('3E unused: ',i2,i5,': ',a)
!      else
! parameters in disordered part read after rewinding
!         notusedpar=notusedpar-1
      endif
      goto 100
!------------------------------------------------------------------
!   elseif(line(2:17).eq.'TYPE_DEFINITION ') then
   case(7) !TYPE_DEFINITION 
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! TYPE_DEFINITION & GES A_P_D BCC_A2 MAGNETIC  -1.0    4.00000E-01 !
      nytypedef=nytypedef+1
      typedefchar(nytypedef)=longline(nextc+1:nextc+1)
! in TC the same typedef "letter" can be used several times
      do ip=1,nooftypedefs
         if(typedefchar(nytypedef).eq.typedefchar(ip)) then
            write(*,*)'3E Same typedef again, "',&
                 typedefchar(nytypedef),'", ignoring second or later occurance'
            nytypedef=nytypedef-1
            goto 88
         endif
      enddo
      nooftypedefs=nytypedef
      if(nooftypedefs.gt.maxtypedefs) then
         write(*,*)'3E Too many TYPE_DEFINITION, modify in readtdb'
         gx%bmperr=4399; goto 1000
      endif
      ip=nextc+3
!      newtypedef: if(index(longline(ip:),' SEQ').gt.0) then
      newtypedef: if(index(longline,' SEQ').gt.0) then
         typedefaction(nytypedef)=100
      else
         km=index(longline,' EINSTEIN ')
         einstein: if(km.gt.0) then
            typedefaction(nytypedef)=1905
            exit newtypedef
         endif einstein
         km=index(longline,' MAGNETIC ')
!         write(*,*)'3E typedef: ',trim(longline),km
         magnetic: if(km.gt.0) then
            ip=km+9
!73           format(a,i3,' "',a,'"')
            call getrel(longline,ip,xxx)
            if(buperr.ne.0) then
               gx%bmperr=buperr; goto 1000
            endif
            if(xxx.eq.zero) then
! this is Qing-Xion magnetic model, next number is 0.37 for BCC or 0.25
               call getrel(longline,ip,xxx)
               if(buperr.ne.0) then
                  gx%bmperr=buperr; goto 1000
               endif
               typedefaction(nytypedef)=int(1.0D2*xxx)
!               write(*,*)'3E Qing-Xiong magnetic model',nytypedef,&
!                    typedefaction(nytypedef)
            else
! this can be -1 for BCC or -3 for FCC, HCP and other phases
               typedefaction(nytypedef)=int(xxx)
            endif
         else
            km=index(longline,' DIS_PART ')
            never=1
!            write(*,*)'3E sigma1: ',trim(longline),km,never
            if(km.eq.0) then
! Allow for NEVER_DIS ...
               km=index(longline,' NEVER')
! this is for disordered SIGMA etc.
               if(never.gt.0) then
                  never=-1
               endif
            endif
!            write(*,*)'3E sigma2: ',trim(longline),km,never
            if(km.gt.0) then
! disordered part, either DIS_PART or NEVER_DIS several checks
               disparttc=disparttc+1
! find the ordered phase name, we have to go backwards from km
               ip=km-1
81             continue
               if(longline(ip:ip).eq.' ') then
                  ordpartph(disparttc)=' '
! The ordpartph is not correct
                  ordpartph(disparttc)=longline(ip+1:km)
! if the ordered part rejected skip this TYPE_DEF
               else
                  ip=ip-1
                  goto 81
               endif
               orddistyp(disparttc)=never
! extract the disordered part phase name
               ip=index(longline(km+2:),' ')
               dispartph(disparttc)=longline(km+2+ip:)
! find the end of phase name, a space or a , there is always a space after ,
               ip=index(dispartph(disparttc),' ')
               km=index(dispartph(disparttc),',')
               if(km.gt.0 .and. km.lt.ip) ip=km
!               if(ip.le.0) ip=1
               dispartph(disparttc)(ip:)=' '
! if ordered part rejected all OK
               do jt=1,nphrej
                  if(ordpartph(disparttc).eq.phreject(jt)) then
                     write(*,*)'3E ordered part rejected, keep disordered'
                     goto 84
                  endif
               enddo
               if(.not.silent) write(kou,82)disparttc, &
                    trim(ordpartph(disparttc)),&
                    trim(dispartph(disparttc)),orddistyp(disparttc)
82             format('3E Found a type_def DIS_PART:',i2,&
                    ' with ',a,' and ',a,' type:',i2)
! if the disordered part phase already entered give warning and advice
               call find_phase_by_name(dispartph(disparttc),iph,ics)
               if(gx%bmperr.ne.0) then
                  gx%bmperr=0
               else
                  if(.not.silent) write(kou,83)dispartph(disparttc)
83                format('3E *** Warning, the disordered phase is already',&
                       ' entered ***'/' Please rearrange the TDB file so',&
                       ' this TYPE_DEF comes before'/&
                       ' the PHASE keyword for the disordered phase: ',a/&
                       ' *** The disordered part ignored ***')
                  disparttc=disparttc-1
                  warning=.TRUE.
               endif
84             continue
            else
               km=index(longline,' NEVER ')
!               write(*,*)'3E sigma3: ',trim(longline),km
               if(km.gt.0) then
! this is for disordered SIGMA etc.
                  write(*,*)'3E Not yet implemented NEVER'
               else
                  typedefaction(nytypedef)=99
                  if(.not.silent) &
                       write(kou,87)nl,longline(1:min(78,len_trim(longline)))
87                format('3E WARNING ignoring TYPE_DEF on line ',i5,':'/a)
                  warning=.TRUE.
!               write(*,*)' WARNING SET TRUE <<<<<<<<<<<<<<<<<<<<<<<<<<<'
               endif
            endif
         endif magnetic
      endif newtypedef
88    continue
!---------------------------------------------------------------------
!   elseif(line(2:20).eq.'LIST_OF_REFERENCES ' .or. &
!          line(2:16).eq.'ADD_REFERENCES ') then
   case(8,9) ! LIST_OF_REFERENCES and ADD_REFERENCES bibliography
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! LIST_OF_REFERENCES
! NUMBER  SOURCE
!   REF283  'Alan Dinsdale, SGTE Data for Pure Elements,
!          Calphad Vol 15(1991) p 317-425,
!          also in NPL Report DMA(A)195 Rev. August 1990'
!       write(kou,*)'Does not handle REFERENCES'
! skip the line with "NUMBER  SOURCE"
! position ip after "NUMBER  SOURCE"
      ip=index(longline,'NUMBER  SOURCE')+14
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'Empty reference line',nl
         gx%bmperr=4313; goto 1000
      endif
      if(longline(ip:ip).eq.'!') then
!         write(*,*)'No references at all'
         goto 100
      endif
!      write(*,*)'list_of_references text length: ',len_trim(longline),ip
! some reference lists like those from SSUB has no single quotes
      kp=index(longline(ip:),"'")
      citationmarks: if(kp.gt.0) then
775      continue
! reference symbol is refx; reference text in reftext
         refx=longline(ip:ip+kp-2)
         if(longline(ip+kp:ip+kp).eq."'") then
! two ' after each other, a dummy reference
            reftext=' '
            ip=ip+kp+1
            kkk=1
!            write(*,*)'dummy: ',refx,' next >',longline(ip:ip+20),'<'
         else
            jp=ip+kp+1+index(longline(ip+kp+1:),"'")
            reftext=longline(ip+kp:jp-2)
            ip=jp
! when all works replace multiple spaces by a single one in reftext
            kkk=len_trim(reftext)
            kp=index(reftext(1:kkk),'  ')
            do while(kp.gt.0)
               reftext(kp:)=reftext(kp+1:)
               kkk=kkk-1
               kp=index(reftext(1:kkk),'  ')
            enddo
         endif
!         write(*,776)refx,nrefs,ip,jp,reftext(1:kkk)
776      format('Reference: ',a,3i5/a)
! this will not create bibliographic references that has not been referenced
         call tdbrefs(refx,reftext(1:kkk),1,ix)
         nrefs=nrefs+1
!         write(*,*)'added biblio ',refx,'>',longline(ip-5:ip+5),'<'
         if(eolch(longline,ip)) then
            gx%bmperr=4313; goto 1000
         endif
         if(longline(ip:ip).ne.'!') then
            kp=index(longline(ip:),"'")
            goto 775
         endif
      else
! references without citation marks
! ip is at the start of the reference id, look for space
         if(.not.silent) write(kou,*) &
              'Cannot handle references without citation marks',nl
         gx%bmperr=4314; goto 1000
      endif citationmarks
777   continue
!      write(*,*)'Read ',nrefs,' references, ending at',nl
!----------------------------------------------------------------
   case(10) ! ASSESSED_SYSTEMS
      if(.not.silent) write(kou,*) &
           '3E cannot handle ASSESSED_SYSTEMS ending at ',nl
!      warning=.TRUE.
! skip lines until !
      do while(index(line,'!').le.0)
         read(21,110)line
         nl=nl+1
!         call replacetab(line,nl)
      enddo
!------------------------------------------------------------------
   case(11) ! DATABASE_INFORMATION
! skip this as checktdb2 has already presented the information
!     if(.not.silent) write(kou,*)'3E Cannot handle DATABASE_INFORMATION at ',nl
!      warning=.TRUE.
! skip lines until !
!      write(*,*)'3E reading database information'
!      write(*,*)'3E ',trim(line)
!      ll=index(line,'!')
!      write(*,*)'3E value of ll: ',ll
! this loop probably meaningless as we have read up to ! already ...
!      write(*,*)'3E found this line: ',nl
      do while(index(line,'!').le.0)
         read(21,110)line
         nl=nl+1
!         ll=index(line,'!')
!         write(*,*)'3E value of ll: ',ll
!         call replacetab(line,nl)
      enddo
!------------------------------------------------------------------
   case(12) ! VERSION, recognize OC1
780   continue
      if(eolch(line,ip)) then
         read(21,110)line
         nl=nl+1
         call replacetab(line,nl)
         goto 780
      else
         if(line(ip:ip).eq.'!' .and. .not.silent) then
            write(kou,*)'Found VERSION keyword but no specification'
         else
            if(line(ip:ip+3).eq.'OC1 ') tdbv=2
         endif
      endif
! skip lines until !
      do while(index(line,'!').le.0)
         read(21,110)line
         nl=nl+1
         call replacetab(line,nl)
      enddo
!------------------------------------------------------------------
   case(13) ! DEFAULT_COMMAND, handle REJECT only
! skip lines until !
      do while(index(line,'!').le.0)
         read(21,110)line
         nl=nl+1
         call replacetab(line,nl)
      enddo
! replace - by _  ... can be dangerous for electrons /-
790   continue
      ip=index(line,'-')
      if(ip.gt.0) then
         line(ip:ip)='_'
         goto 790
      endif
! here I handle only reject phase
791   continue
      call getext(line,nextc,1,name1,' ',ix)
      if(name1(1:ix).eq.'REJECT_PHASE') then
793      continue
! save phase names to be rejected in a structure            
         call getext(line,nextc,1,name1,' ',ix)
         if(name1(1:1).eq.' ' .or. name1(1:1).eq.'!') then
            goto 794
         else
            nphrej=nphrej+1
            if(nphrej.gt.maxrejph) then
               write(*,*)'3E Too many phases to reject, increase maxrejph'
            else
               write(*,*)'3E rejected phase: ',name1
               phreject(nphrej)=name1
            endif
         endif
         goto 793
      elseif(name1(1:7).eq.'DEF_SYS' .or. &
           name1(1:13).eq.'DEFINE_SYSTEM') then
! ignore default define_system... as, Va and /- are always entered by default
         continue
      else
         write(*,*)'3E WARNING: ignoring default command: ',trim(name1)
      endif
794   continue
! rejected phases OK
!      do zz=1,nphrej
!         write(*,*)'3E rejected phase: ',phreject(zz)
!      enddo
!--------------------------------- DEFINE
      case(14) !ignore without warning
         write(*,*)'3E ignoring DEFINE keyword'
         continue
   end select
!-------------------------------------------------------- end select
   if(gx%bmperr.ne.0 .and. .not.silent) then
      write(kou,711)gx%bmperr,nl,trim(line)
711   format('3E error: ',i5,' around line ',i7,': '/a)
! this error means reference error
      if(gx%bmperr.eq.4154) gx%bmperr=0
   endif
! look for next KEYWORD
   goto 100
!--------------------------------------------------------
!----- reading functions at the end from a TDB file
800 continue
   if(eolch(line,nextc)) then
      if(.not.silent) write(kou,*) &
           'Function name must be on same line as FUNCTION'
      gx%bmperr=4315; goto 1000
   endif
   ipp=nextc+index(line(nextc:),' ')
   name1=line(nextc:ipp-1)
!         write(*,18)'function >',name1,'< ',nextc,ipp
!18       format(a,a,a,2i4)
! old code
   longline=' '
   longline=line(ipp:)
810 continue
   jp=max(len_trim(longline),1)
!   write(*,811)jp,longline(jp:jp),longline(1:jp)
811 format('3E ll: ',i3,' "',a1,'" ',a)
   if(longline(jp:jp).eq.'!') then
! replace # by ' '
820   continue
      jss=index(longline(1:jp),'#')
      if(jss.gt.0) then
         longline(jss:jss)=' '
         goto 820
      endif
! check if function is entered as undefined, exact match of name required
      call find_tpfun_by_name_exact(name1,nr,notent)
      if(gx%bmperr.eq.0) then
         if(notent) then
!            write(*,*)'Entering function: ',name1
! entering a function may add new unentered functions ... last argument TRUE
!            write(*,*)'3E Entering function 3: ',name1,len_trim(longline)
!            lrot=0
!            call store_tpfun(name1,longline,lrot,.TRUE.)
            call store_tpfun(name1,longline,lrot,rewindx)
            if(gx%bmperr.ne.0) then
! one may have error here
               if(.not.silent) write(kou,*)'Failed entering function: ',name1
               goto 1000
            endif
            if(ocv()) write(*,*)'Entered function: ',name1
            nofunent=nofunent+1
         else
!            write(*,*)'3E referenced: ',trim(name1),nr,&
!                 tpfuns(nr)%rewind,rewindx
            if(tpfuns(nr)%rewind.eq.rewindx) then
! Function entered and referenced, check if duplicate!
               write(*,828)trim(name1),nl,rewindx
828            format('3E WARNING duplicate function ',a,' at line: ',2i5)
               warning=.TRUE.
            endif
         endif
      else
! igore the function as it is not referenced.  Reset error code
         gx%bmperr=0
      endif
   else
830   continue
      nl=nl+1
      read(21,110)line
!            write(kou,101)'readtdb 2: ',nl,line(1:40)
! skip lines with a $ in first position
      if(line(1:1).eq.'$')goto 830
      call replacetab(line,nl)
      longline=longline(1:jp)//line
      goto 810
   endif
   goto 100
!   endif barafun
!---------------------------------------------------------
! We have now read all
!--------------------------------------------------------
1000 continue
!   write(*,1111)totalpar,totalpar-notusedpar
!   write(*,1111)totalpar,enteredpar,notusedpar
   write(*,1111)totalpar,enteredpar
1111 format(/'Out of ',i5,' model parameters ',i5,' have been entered'/)
   if(warning) then
1001  continue
      write(*,*)
! if silent set ignore warnings
      if(.not.silent) then
         do jss=1,nundefmpi
            write(*,1008)undefmpi(jss)
1008        format('3E *** WARNING unused model parameter identifier ',a,&
                 ' in some phases')
         enddo
         write(kou,1003)
1003     format(/'There were warnings, check them carefully'/&
              'and press RETURN if you wish to continue.')
         read(kiu,1004)ch1
1004     format(a)
!         if(ch1.eq.'N') stop 'warnings reading database'
!         if(ch1.ne.'Y') then
!            write(kou,*)'Please answer Y or N'
!            goto 1001
!         endif
      endif
   endif
!   write(*,*)'3E At label 1000'
   if(buperr.ne.0 .or. gx%bmperr.ne.0) then
      if(gx%bmperr.eq.0) gx%bmperr=buperr
      if(.not.silent) write(kou,1002)gx%bmperr,nl
1002   format('Error ',i5,', occured at TDB file line ',i7)
!      write(*,*)'Do you want to continue at your own risk anyway?'
!      read(*,1008)ch1
!1008  format(a)
!      if(ch1.eq.'Y') then
!         write(*,*)'Now any kind of error may occur .... '
!         buperr=0
!         gx%bmperr=0
!         goto 100
!      endif
   endif
   close(21)
! read numbers, value after / is maximum
! endmember, interactions, property,
! tpfuns, composition sets, equilibria
! state variable functions, references, additions
   if(ocv()) write(*,1007)noofel,maxel,noofsp,maxsp,noofph,maxph,&
        noofem,100000,noofint,100000,noofprop,100000,&
        notpf(),maxtpf,highcs,2*maxph,eqfree-1,maxeq,&
        nsvfun,maxsvfun,reffree-1,maxrefs,addrecs,csfree-1
1007 format('Created records for elements, species, phases: ',2x,&
          3(i4,'/',i4,1x)/&
          'end members, interactions, properties: ',10x,&
          3(i4,'/',i4,1x)/&
          'TP-funs, max and free composition sets, equilibria: ',10x,&
          3(i4,'/',i4,1x)/&
          'state variable functions, references, additions: ',&
          3(i4,'/',i4,1x)/)
! a special warning message as it may be scrolled away by all references
!   write(*,*)'Any warnings?',warning
! nonzero multiuse will prompt a warning in the monitor
   firsteq%multiuse=0
   if(gx%bmperr.eq.0 .and. warning) firsteq%multiuse=-1
   return
!--------------------------------------------------------------------------
! errors and rewinds
1010 continue
   if(.not.silent) write(kou,*)'I/O error opening file: ',gx%bmperr
   return
!-----------------------------------------------------
! end of file found, act differently if reading functions
2000 continue
   rewindx=rewindx+1
   rewindfile: if(dodis.eq.0 .and. disparttc.gt.0) then
! rewind to read disordred parts
      if(.not.silent) &
           write(kou,*)'3E Rewinding to read disordered part of phases'
      rewind(21)
      dodis=1
      nl=0
      goto 100
   elseif(.not.onlyfun) then
! rewind to read referenced functions and references !!
      dodis=2
      rewind(21)
      onlyfun=.TRUE.
      nofunent=0
!      write(*,2002)gx%bmperr
2002 format('Found end-of-file, rewind to find functions',i5)
      nl=0
      goto 100
   elseif(nofunent.gt.0) then
! rewind if there were functions entered last time
      rewind(21)
      norew=norew+1
!     write(*,*)'Found functions: ',nofunent,' rewinding again',norew,gx%bmperr
!      if(newfun.gt.0) then
!          write(*,*)'Read ',newfun+nfail,' functions, entered ',newfun,&
!               ' rewinding ',norew
!         newfun=0
      nofunent=0
      nl=0
      goto 100
   else
! check if there are any unentered functions
      call list_unentered_funs(kou,nr)
      if(nr.gt.0) then
         if(.not.silent) write(kou,*)'Number of missing function: ',nr
         gx%bmperr=4186
      endif
! check if any function not entered
      onlyfun=.FALSE.
   endif rewindfile
   goto 1000
! end of file while looking for ! terminating a keyword
2200 continue
   if(.not.silent) write(kou,2210)nl,longline(1:72)
2210 format('End of file at ',i5,' looking for end of keyword:'/a)
   gx%bmperr=4316
   goto 1000
 end subroutine readtdb

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine readtdbsilent
!\begin{verbatim}
 subroutine readtdbsilent
!\end{verbatim} %+
   globaldata%status=ibset(globaldata%status,GSSILENT)
   return
 end subroutine readtdbsilent

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine any_disordered_part
!\begin{verbatim}
 subroutine any_disordered_part(lin,ndisph,dispartph,ordpartph,orddistyp)
! reading data from a TDB/PDB file with selection of elements
! extract only TYPE_DEFS for order/disorder
!-------------------------------------------------------
! Not all TYPE_DEFS implemented
!-------------------------------------------------------
   implicit none
   integer lin,ndisph
   character dispartph(*)*(*),ordpartph(*)*(*)
   integer orddistyp(*)
!\end{verbatim}
   character line*128,longline*1024,phase*24
   integer ip,jp,zp
   ndisph=0
   loop1: do while(.true.)
      read(lin,100,end=900)line
100   format(a)
      ip=1
      if(eolch(line,ip)) cycle loop1
      if(line(ip:ip).eq.'$') cycle loop1
      typedef: if(line(ip:ip+7).eq.'TYPE_DEF') then
! search for ! meaning end of keyword
         longline=line(ip:)
         ip=len_trim(longline)
         loop2: do while(longline(ip:ip).ne.'!')
            read(lin,100,end=900)line
            longline(ip+1:)=line
            ip=len_trim(longline)
         enddo loop2
!         write(*,*)'3E type_def 1: ',longline(1:ip)
! the important part is "GES" followed by "A_P_D" (or AMEND_PHASE_DEFINITION)
! followed by phase name and followed by "DIS_PART" or "NEVer" and a phase name
         jp=index(longline,' GES ')
         if(jp.le.0) exit typedef
! below is a clumsy way to extract phase names for ordered/disordered parts
! skip the first item after "GES" (should be AMEND_PHASE_DESCRIPTION or abbrev
         zp=jp+4
         if(eolch(longline,zp)) then
            exit typedef
         endif
! skip to next item, if none loop, else extract next item
         jp=index(longline(zp:),' ')
         if(zp.le.0) exit typedef
         zp=zp+jp
!         write(*,*)'3E zp1: ',trim(longline),zp
         if(eolch(longline,zp)) then
            exit typedef
         endif
         phase=longline(zp:)
         jp=index(phase,' ')
         phase(jp:)=' '
! this should be a phase name, maybe the ordered part
!         write(*,*)'3E phase name? ',trim(phase)
! check what comes after phase name, 
! the important is "DISORDERED_PART" or "NEVER_DISORDERED" or abbreviations
         ip=zp+jp
!         write(*,*)'3E after phase: "',trim(longline(ip:)),'"'
         if(eolch(longline,ip)) exit typedef
         if(longline(ip:ip+2).eq.'NEV') then
            ndisph=ndisph+1
            orddistyp(ndisph)=1
         elseif(longline(ip:ip+2).eq.'DIS') then
            ndisph=ndisph+1
            orddistyp(ndisph)=-1
         else
            exit typedef
         endif
!         write(*,'(a,a,2i3)')'3E zp2: ',trim(longline),orddistyp(ndisph),ip
! find space after NEVER or DIS_PART
         jp=ip+index(longline(ip:),' ')
! disordered phase name should be now
!         write(*,*)'3E after dis/nev: "',trim(longline(ip:)),'"',jp
         if(eolch(longline,jp)) then
            write(*,'(a/a)')'3E no disordered phase! ',trim(longline)
            exit typedef
         endif
         dispartph(ndisph)=longline(jp:)
! name terminated by space, comma or !
         ip=index(dispartph(ndisph),',')
         if(ip.gt.0) dispartph(ndisph)(ip:)=' '
         ip=index(dispartph(ndisph),'!')
         if(ip.gt.0) dispartph(ndisph)(ip:)=' '
!         ip=index(dispartph,' ')
!         dispartph(ndisph)(ip:)=' '
         ordpartph(ndisph)=phase
!         write(*,'(a,i2,5a,i3)')'3E ord/dis: ',ndisph,' "',&
!              trim(ordpartph(ndisph)),'" + "',trim(dispartph(ndisph)),'"',&
!              orddistyp(ndisph)
      endif typedef
   enddo loop1
! eof
900 continue
   rewind(lin)
1000 continue
   return
 end subroutine any_disordered_part

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine readpdb
!\begin{verbatim}
 subroutine readpdb(filename,nel,selel,options)
! reading data from a PDB file with selection of elements
   implicit none
   integer nel
   character filename*(*),selel(*)*2,options*(*)
!\end{verbatim}
   character line*150,elsym*2,name1*24,name2*80,elsyms(10)*2
   character longline*10000,reftext*512
   character phtype*1,ch1*1,const(maxsp)*24,name3*24,funname*60,name4*60
   character refx*16
   character (len=1), dimension(20) :: typedefchar
   integer, dimension(20) :: typedefaction
   integer, dimension(5) :: addphasetypedef
   double precision mass,h298,s298
   integer, dimension(10) :: knr,endm
! lint(1,*) is sublattice, lint(2,*) is species
   double precision stoik(10),xsl,xxx
   integer lint(2,3),nytypedef,nextc,keyw,tdbv
   integer typty,fractyp,lp1,lp2,ix,jph,kkk,lcs,nint,noelx
   logical onlyfun,nophase,ionliq,notent
   integer norew,newfun,nfail,nooftypedefs,nl,ipp,jp,jss,lrot,ip,jt
   integer nsl,ll,kp,nr,nrr,mode,lokph,lokcs,km,nrefs,ideg,iph,ics
! disparttc and dispartph to handle phases with disordered parts
   integer nofunent,disparttc,dodis,jl,nd1,thisdis
   integer excessmodel,modelcode,noofadds,noofdet,permut,enteredpar
   character*24 dispartph(5),ordpartph(5)
   integer add(10),allel,cbug,rewindx,orddistyp(5),emodel
   character modelname*72
   integer, parameter :: nadditions=6
   character*128, dimension(10) :: detail
   character*12, dimension(8), parameter :: addition=&
        ['IMAGF       ','IMAGB       ','IWMAGF      ','IWMAGB      ',&
         'DEBYE1      ','EINSTEIN1   ','            ','            ']
   logical warning,verbose,subord
! set to TRUE if element present in database
   logical, allocatable :: present(:)
! to prevent any output
   logical silent,havedisorder
!  mmyfr
! if warning is true at the end pause before listing bibliography
!
   emodel=0
   warning=.FALSE.
   silent=.FALSE.
   verbose=.TRUE.
!   if(btest(globaldata%status,GSSILENT)) then
!      silent=.TRUE.
!      write(*,*)'3E reading database silent'
!   endif
!   write(*,*)'reading a PDB file'
   if(.not.(index(filename,'.pdb').gt.0 &
       .or. index(filename,'.PDB').gt.0)) then
! no extention provided
      filename(len_trim(filename)+1:)='.PDB'
   endif
   if(nel.gt.0) then
      allocate(present(nel))
      present=.FALSE.
   else
      allel=0
      allocate(present(maxel))
      present=.FALSE.
   endif
! disparttc counts the number of disordered phases to read, the
! phase names are in dispartph(1..disparttc)
! dodis is nonzero only when reading the disordered part of phases.
   disparttc=0
   dodis=0
   open(21,file=filename,access='sequential',form='formatted',&
        err=1010,iostat=gx%bmperr,status='old')
   onlyfun=.FALSE.
   tdbv=1
   norew=0
   newfun=0
   nfail=0
   nrefs=0
   nooftypedefs=0
   enteredpar=0
! nophase set false after reading a PHASE keyword, 
! expecting next keyword to be CONSTITUENT
   nophase=.TRUE.
! return here after rewind
90  continue
   nl=0
! return here to look for a new keyword, end-of-file OK here
100 continue
   read(21,110,end=2000)line
110 format(a)
   nl=nl+1
   if(len_trim(line).gt.128) then
      write(*,121)nl
121   format(' *** Warning: line ',i5,' has characters beyond position 128,'/&
           'some information may be lost')
   endif
! look for first nonblank character
   ipp=1
   if(eolch(line,ipp)) goto 100
   if(line(ipp:ipp).eq.'$') goto 100
! replace TAB by space
   call replacetab(line,nl)
! debug
!  write(*,*)'input: ',trim(line),nl
!---------------------------------------------------------
! handle all PDB keywords except function, 5
120 continue
   keyw=ispdbkeyword(line,nextc)
!   write(*,*)'Keyword? ',trim(line),keyw,nophase
   if(.not.onlyfun) then
      if(keyw.eq.0) then
         write(*,123)trim(line)
123      format(//' *** Warning, not a PDB keyword: ',a)
      endif
   endif
   if(keyw.eq.0) then
      ip=1
      if(.not.eolch(line,ip)) then
         if(ocv()) write(*,*)'Ignoring line: ',nl,ip,line(ip:ip+20)
      endif
      if(.not.onlyfun) write(*,*)' *** Ignoring unknown keyword: ',line(1:24)
!      read(*,122)ch1
122   format(a)
      goto 100
   endif
   if(.not.nophase .and. keyw.ne.4) then
! after a PHASE keyword one should have a CONSTITUENT
      if(.not.silent) write(kou,*)'expeciting CONSTITUENT: ',trim(line)
      warning=.TRUE.
   endif
! check there is a ! in line, otherwise read until we find an exclamation mark
!   ip=1
   longline=line
   ip=len_trim(longline)+1
   do while(index(longline,'!').le.0)
      read(21,110,err=2200)line
      nl=nl+1
      if(line(1:1).ne.'$') then
         call replacetab(line,nl)
         longline(ip:)=line
         ip=len_trim(longline)+1
         if(ip.gt.len(longline)-80) then
            if(.not.silent) write(kou,69)nl,ip,longline(1:72)
69          format('Overflow in longline ',2i5,' for line starting:'/a)
            gx%bmperr=4304; goto 1000
         endif
      endif
   enddo
!
   if(onlyfun) then
      if(keyw.eq.5) then
!         write(*,*)'At 100: ',trim(longline)
         goto 800
      endif
      goto 100
   endif
!
  select case(keyw)
! case default means keyw not understood
   case default
      write(*,*)'no such keyword? ',trim(longline)
   case(1) !element ------------------------------------------------
! ELEMENT CR   BCC_A2                    5.1996E+01  4.0500E+03  2.3560E+01!
! SAME AS TDB
      ip=nextc
      if(eolch(longline,ip)) then
         write(kou,*)'No element name after ELEMENT keyword on line ',nl
         gx%bmperr=4305; goto 1000
      endif
      elsym=longline(ip:ip+1)
      if(elsym.eq.'/-' .or. elsym.eq.'VA') goto 100
! allow lower case in TDB file ...
      call capson(elsym)
      if(nel.gt.0) then
! check if element among selected, if nel=0 accept all
         do jt=1,nel
            if(elsym.eq.selel(jt)) goto 76
         enddo
! ignore this element as not selected
         goto 100
      else
         allel=allel+1
         jt=allel
      endif
! mark we found a selected element
76    continue
      if(allocated(present)) then
         present(jt)=.TRUE.
      endif
! we seem to miss the first letter of the reference state below ??
      ip=ip+len_trim(elsym)-1
      if(eolch(longline,ip)) then
         name1='DUMMY'
         mass=one
         h298=zero
         s298=zero
      else
! extract the reference phase, third argument is 1 meaning until next space
! ix is the length of the reference phase (irrelevant here)
! ip is updated to character after the name extracted
         call getext(longline,ip,1,name1,' ',ix)
! after the name should be mass, H298-H0 and S298, ignore errors
         call getrel(longline,ip,mass)
         if(buperr.ne.0) then
            mass=one; buperr=0
         endif
         call getrel(longline,ip,h298)
         if(buperr.ne.0) then
            h298=zero; buperr=0
         endif
         call getrel(longline,ip,s298)
         if(buperr.ne.0) then
            s298=zero; buperr=0
         endif
! in OC one can have the element name in addition to its symbol
         name2=elsym
      endif
      call store_element(elsym,name2,name1,mass,h298,s298)
      if(gx%bmperr.ne.0) goto 1000
   case(2) !SPECIES -------------------------------------------------
! SPECIES O3PU2                       O3PU2!
! SAME AS TDB
      ip=nextc
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'Line after SPECIES keyword empty'
         gx%bmperr=4306; goto 1000
      endif
      name1=longline(ip:)
! find first space after non-space
      jp=index(name1,' ')
      name1(jp:)=' '
      ip=ip+jp
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'No stoichiometry for species: ',name1
         warning=.TRUE.
         goto 100
      endif
      name2=longline(ip:)
      jp=index(name2,' ')
      name2(jp:)=' '
      call decode_stoik(name2,noelx,elsyms,stoik)
      if(gx%bmperr.ne.0) goto 1000
! check elements exist
      call enter_species(name1,noelx,elsyms,stoik)
!      write(*,*)'3E: entering species error: ',gx%bmperr
      if(gx%bmperr.ne.0) then
! if element not selected just skip the species
         if(gx%bmperr.eq.4046) then
            gx%bmperr=0; goto 100
         else
            if(.not.silent) write(kou,*)'Error entering species: ',name1,name2
            goto 1000
         endif
      endif
!-----------------------------------------------------------------------
   case(5) ! function
! all functions entered at the end at label 800
! SAME AS TDB
      if(index(longline,'!').le.0) then
         write(*,*)' Error, terminating ! not found for function!!',nl
         gx%bmperr=4307; goto 1000
      endif
!-------------------------------------------------------------------------
!   elseif(line(2:7).eq.'PHASE ') then
   case(3) ! PHASE
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! PHASE LIQUID MODEL_CODE ADDITION_CODE DETAIL_CODE !
      if(nophase) then
         nophase=.false.
! give a warning if any selected element is not present
         if(allocated(present)) then
            funname=' '
            kkk=1
            do jt=1,max(nel,allel)
               if(.not.present(jt)) then
                  funname(kkk:)=selel(jt)
                  kkk=len_trim(funname)+2
               endif
            enddo
            if(kkk.gt.1) then
               if(.not.silent) write(kou,68)funname(1:kkk)
68             format(/' *** Warning, elements not present in database: ',a/)
            endif
            deallocate(present)
         endif
      else
         if(.not.silent) write(kou,*) &
              'Error, a PHASE keyword must be followed by its CONSTIT'
         gx%bmperr=4308; goto 1000
      endif
      ip=nextc
!      if(eolch(longline,ip)) then
!         if(.not.silent) write(kou,*)'line after PHASE empty'
!         goto 100
!      endif
! this extracts the phase name from position ip terminated by a space
! ix is the length of the phase name, ip is updated to position after the name
! note getext increments ip by 1 before extracting the name!
      call getext(longline,ip,1,name1,' ',ix)
      if(name1(1:1).eq.' ') then
         write(kou,*)'Missing phase name after PHASE keyword'
         gx%bmperr=4399; goto 1000
      endif
! convert phase name to upper case
      call capson(name1)
! phytype was used in TDB to specify gas and liquid, not used here
      phtype=' '
!-----------------------------------------------
! extract phase model codes like CEF3D(10.0 4.0 16.0)
! additions and details ...
! NEED to modify getext to read up to )
      call getext(longline,ip,1,name2,' ',ix)
      if(name2(1:1).eq.' ') then
         write(kou,*)'Missing phase model'
         gx%bmperr=4399; goto 1000
      endif
      call capson(name2)
      ix=index(name2,'(')
      if(ix.gt.0) then
         modelname=name2
       else
         modelname=name2
       endif
! if name2 contains a ( but no ) append that part ...
      if(index(name2,'(').gt.0) then
         if(index(name2,')').eq.0) then
            jp=index(longline(ip:),')')
!      write(*,178)'CEFx: ',ip,jp,longline(ip:ip+jp)
            name2(len_trim(name2)+1:)=longline(ip:ip+jp)
!      write(*,178)'CEFy: ',trim(name2)
            ip=ip+jp
         endif
      endif
! known models: SUBxxx CEFxxx I2SL 
! check model is known
      thisdis=0
      permut=0
      ionliq=.FALSE.
      havedisorder=.FALSE.
      if(name2(1:6).eq.'IDEAL ') then
! ideal model
         stoik(1)=one
         modelcode=0
         nsl=1
      elseif(name2(1:3).eq.'SUB') then
! substitutional model
         modelcode=1
         nsl=1
         stoik(1)=one
         if(name2(4:7).eq.'RKM ') then
            excessmodel=1
         elseif(name2(4:6).eq.'PK ') then
            excessmodel=2
         elseif(name2(4:7).eq.'MIX ') then
            excessmodel=3
         else
            write(*,*)'Unknown model for ',trim(name1),': ',trim(name2)
            gx%bmperr=4399; goto 1000
         endif
      elseif(name2(1:3).eq.'CEF') then
! CEF model, the 4th letter is the number of sublattices
         modelcode=2
         nsl=ichar(name2(4:4))-ichar('0')
         if(nsl.lt.1 .or. nsl.gt.9) then
            write(*,*)'Wrong number of sublattices for phase: ',trim(name1),nsl
            gx%bmperr=4399; goto 1000
         endif
         jp=5
         if(name2(5:5).eq.'D') then
! there can be a disordered part and ordered part should be subtracted
            thisdis=1
            jp=6
            havedisorder=.TRUE.
            subord=.TRUE.
         elseif(name2(5:5).eq.'N') then
! there can be a disordered part and the ordered part should NOT be subtracted
! like for SIGMA etc
            thisdis=1
            jp=6
            havedisorder=.TRUE.
            subord=.FALSE.
         elseif(nsl.ge.4) then
!            write(*,*)'3E checking permut: ',name2(5:5)
! >=4 sublattice models can have permutations
            if(name2(5:5).eq.'F') then
               permut=1
            elseif(name2(5:5).eq.'B') then
               permut=2
            endif
            if(permut.gt.0) then
               if(name2(6:6).eq.'D') then
                  havedisorder=.TRUE.
                  subord=.TRUE.
                  jp=7
               else
                  jp=6
               endif
            else
               jp=5
            endif
         endif
! here jp should be at the position of (
         if(name2(jp:jp).ne.'(') then
            write(*,*)'Missing ( after CEFj ',trim(name1),': ',trim(name2),jp
            gx%bmperr=4399; goto 1000
         else
            jp=jp+1
!           write(*,181)'CEF: ',nsl,jp,trim(name2)
181         format(a,2i3,1x,a)
            stoik=zero
            do ll=1,nsl
               call getrel(name2,jp,stoik(ll))
               if(buperr.ne.0) then
                   write(*,*)'Site ratio wrong: ',trim(name1),ll,buperr
                   gx%bmperr=4399; goto 1000
               endif
            enddo
! maybe there should be a check that ) is the next character?? Can be a space?
!            write(*,182)'Sites: ',nsl,(stoik(ll),ll=1,nsl)
182         format(a,i2,9(1pE12.4))
         endif
      elseif(name2(1:5).eq.'I2SL ') then
! ionic liquid model
         modelcode=3
         ionliq=.TRUE.
         stoik(1)=one
         stoik(2)=one
      elseif(name2(1:5).eq.'KFGL ') then
! Kaphor-Frohberg-Gye-Lehman IRSID slag model
         modelcode=4
         stoik(1)=one
      else
         write(*,*)'Unknown model for phase: ',trim(name1)
         gx%bmperr=4399; goto 1000
      endif
! additions (or details), there can be none or several
      noofadds=0
      noofdet=0
180   continue
!      write(*,*)'phase: ',ip,trim(line)
      call getext(longline,ip,1,name2,' ',ix)
! if name2 contains a ( but no ) append that part ...
! here detail can be complicated with internal ( ... ) ....
      if(index(name2,'(').gt.0) then
         if(index(name2,')').eq.0) then
            jp=index(longline(ip:),')')
            name2(len_trim(name2)+1:)=longline(ip:ip+jp)
            ip=ip+jp
         endif
      endif
      if(name2(1:1).ne.' ') then
! additions are IMAGF, IMAGB, IWMAGF, IWMAGB, DEBYE1, EINSTEIN1
         add1: do jp=1,nadditions
            if(name2(1:12).eq.addition(jp)) then
               noofadds=noofadds+1
               add(noofadds)=jp
               goto 180
            endif
         enddo add1
! details will be handled later, save the detail text
         if(name2(1:6).eq.'DETAIL') then
            noofdet=noofdet+1
            detail(noofdet)=name2
         endif
      endif
!---------------------------------------------------------------------
! The constituent line must follow PHASE before any new phase
   case(4) !    CONSTITUENT LIQUID:L :CR,FE,MO :  !
! the phase must have been defined
      if(nophase) then
         if(.not.silent) write(kou,*) &
              'A CONSTITUENT keyword not directly preceeded by PHASE!'
         gx%bmperr=4308; goto 1000
      endif
      nophase=.true.
      condis1: if(dodis.eq.1) then
         if(thisdis.eq.0) goto 100
! we skip the constituent line and go directly to create disordered fractions
         goto 395
      elseif(disparttc.gt.0 .and. thisdis.lt.0) then
! this is a disordered part, skip
         goto 100
      endif condis1
!360    continue
      jp=len_trim(longline)
!      write(*,*)'readpdb gas1: ',nl,jp,longline(1:jp)
! eliminate all after the exclamation mark
!      longline(jp+1:)=' '
! 
      ip=index(longline,' :')+2
!      write(*,*)'readpdb gas2: ',jp,longline(1:jp)
      ll=0
      nr=0
      nrr=0
!      write(*,*)'readpdb 3C: ',ll,nr,nsl,longline(ip:jp)
! mode=1 indicates to getname that / + - are allowed in species names
      mode=1
370   continue
      if(ll.ge.1) then
         knr(ll)=nr
         if(nr.le.0) then
            if(ocv()) then
               write(*,*)'Skipping phase due to missing constituents: ',name1
!              write(*,378)name1,ll
378            format('Phase ',a,' has no constituents in sublattice ',i2)
! Not a fatal error when elements have been selected but skip this phase
            endif
            goto 100
         endif
      endif
      ll=ll+1
!      write(*,*)'start sublat ',ll,nsl,nr,ip
      if(ll.gt.nsl) goto 390
      nr=0
380   continue
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'Error extracting constituents 1'
         gx%bmperr=4309; goto 1000
      endif
      nr=nr+1
      nrr=nrr+1
!      write(*,379)'readpdb 3CXX: ',ip,nr,longline(ip:ip+10)
379   format(a,2i4,' >',a,'< >',a,'< >',a,'<')
      call getname(longline,ip,name3,mode,ch1)
!      write(*,379)'readpdb 3CY: ',ip,nr,longline(ip:ip+10),name3,ch1
      if(buperr.ne.0) then
!         write(*,381)'readpdb 3E: ',ll,nr,longline(1:ip+5),ip,name3
381      format(a,2i4,' "',a,'" ',i5,1x,a,'"',a)
         gx%bmperr=buperr; goto 1000
      endif
!      write(*,381)'readpdb 3E: ',ll,nr,longline(1:ip+5),ip,name3,ch1
      const(nrr)=name3
! bypass any "major" indicator %
      if(ch1.eq.'%') ip=ip+1
      if(eolch(longline,ip)) then
!         if(.not.silent) write(kou,*)'Error extracting constituents 2'
         gx%bmperr=4309; goto 1000
      endif
! check that const(nrr) among the selected elements ...
!      write(*,*)'Testing constituent: ',name3,nr
      call find_species_record_noabbr(name3,lp1)
      if(gx%bmperr.ne.0) then
! this species is not present, not a fatal error, skip it and continue
!         write(*,*)'Skipping constituent: ',name3
         gx%bmperr=0; nrr=nrr-1; nr=nr-1
      endif
      ch1=longline(ip:ip)
      if(ch1.eq.',') then
         ip=ip+1; goto 380
      elseif(ch1.eq.':') then
         ip=ip+1; goto 370
      endif
      if(ch1.ne.'!') goto 380
! when an ! found the list of constutents is finished.  But we
! should have found a : before the !
      if(.not.silent) write(kou,*)'Found "!" before terminating ":"'
      gx%bmperr=4310; goto 1000
!      write(*,*)'Species terminator error: ',ch1,nl
!      gx%bmperr=4157; goto 1000
390    continue
! name2 is model, ignored on reading TDB
!      ionliq=.FALSE.
!      if(phtype(1:1).eq.'Y') then
!         name2='IONIC_LIQUID '
!         ionliq=.TRUE.
!      else
!         name2='CEF-TDB-RKM? '
!      endif
      if(ocv()) write(*,*)'readpdb 9: ',name1,nsl,knr(1),knr(2),phtype
395   continue
! disordered phase is part of model
      condis2: if(dodis.eq.1) then
! if we have a disordered part do not enter the phase, add disordered fracs!
! the ordered phase name is ordpart(thisdis)
!         call find_phase_by_name(ordpartph(thisdis),iph,ics)
!         if(gx%bmperr.ne.0) then
! NOTE THE ORDERED PHASE MAY NOT BE ENTERED DUE TO COMPONENTS!!
!            if(.not.silent) write(kou,396)thisdis,ordpartph(thisdis)
!396         format('Disordered phase skipped as no ordered: ',i3,' "',a,'"')
!            warning=.TRUE.
!            gx%bmperr=0
!            goto 100
!         else
!            if(.not.silent) write(kou,*) &
!                 '3E Adding disordered fraction set: ',ordpartph(thisdis)
!         endif
! we are creating the phase, there is only one composition set
!         call get_phase_compset(iph,1,lokph,lokcs)
!         if(gx%bmperr.ne.0) goto 1000
! ch1 is suffix for parameters, always D
!         ch1='D'
! jl=0 if NDM (sigma)
! jl=1 if phase can be totally disordered (but can have interstitials)
! nd1 is the number of sublattices to sum into disordered set
         if(dispartph(thisdis)(1:7).eq.'FCC_A1 ' .or. &
              dispartph(thisdis)(1:7).eq.'BCC_A2 ' .or. &
              dispartph(thisdis)(1:7).eq.'HCP_A3 ') then
! if disordred phase is FCC, BCC or HCP then set jl=1 and nd1 to 2 or 4
            if(phlista(lokph)%noofsubl.le.5) nd1=4
            if(phlista(lokph)%noofsubl.le.3) nd1=2
            if(.not.silent) write(kou,397) &
                 ordpartph(thisdis)(1:len_trim(ordpartph(thisdis))),nd1
!         write(*,399)
!399      format('Phase names for disordered parts of FCC, BCC and HCP must',&
!              ' start with:'/'  A1_ , A2_ and A3_ respectivly!'/&
!              ' and have an interstitial sublattice')
            jl=1
         elseif(dispartph(thisdis)(1:3).eq.'A1_' .or. &
              dispartph(thisdis)(1:3).eq.'A2_' .or. &
              dispartph(thisdis)(1:3).eq.'A3_') then
! if disordred phase is FCC, BCC or HCP then set jl=1 and nd1 to 2 or 4
            if(phlista(lokph)%noofsubl.le.5) nd1=4
            if(phlista(lokph)%noofsubl.le.3) nd1=2
            if(.not.silent) write(kou,397) &
                 ordpartph(thisdis)(1:len_trim(ordpartph(thisdis))),nd1
397         format('Phase ',a,&
                 ' has an order/disorder partition model summing first ',i2)
            jl=1
         elseif(dispartph(thisdis)(1:4).eq.'DIS_') then
! disordered part of sigma, mu etc.
            jl=0; nd1=phlista(lokph)%noofsubl
         else
! probably disordered part of sigma, mu etc.
            write(kou,495)trim(ordpartph(thisdis))
495         format(' *** WARNING: Non-standard name of disordered phase: ',a/&
                 ' Assuming ordered phase will never be disordered',&
                 ' like sigma, mu etc')
            jl=0; nd1=phlista(lokph)%noofsubl
         endif
!         if(jl.eq.0 .and. .not.silent) write(kou,398)trim(ordpartph(thisdis))
398      format(' Assuming that phase ',a,' cannot completely disorder')
! add DIS_PART from TDB
         write(*,*)'3E adding fraction set 2:',csfree,highcs
         call add_fraction_set(iph,ch1,nd1,jl)
         if(gx%bmperr.ne.0) then
            if(.not.silent) write(kou,*) &
                 '3E Error entering disordered fraction set: ',gx%bmperr
            goto 1000
         endif
         write(*,*)'3E added fraction set 2:',csfree,highcs
         if(jl.eq.0) then
! we must set the correct formula unit of the disordered phase, on the
! TDB file it is unity.  Sum up the sites for the ordered phase in lokcs
            xxx=zero
            do ll=1,nd1
               xxx=xxx+firsteq%phase_varres(lokcs)%sites(ll)
            enddo
            firsteq%phase_varres(lokcs)%disfra%fsites=xxx
         else
            xxx=one
         endif
!         if(.not.silent) write(kou,601) &
!              dispartph(thisdis)(1:len_trim(dispartph(thisdis))),ch1,nd1,jl,xxx
601      format('3E Add parameters from disordered part: ',a,5x,a,2x,2i3,F12.4)
      else
!--- ENTER PHASE with constituents and model.  Modelcode is CEF etc
! in OC phtype G means it is first, L that it is first after G
! remaining phases ordered alphabetically
         if(name1(1:4).eq.'GAS ') then
            phtype='G'
         elseif(name1(1:7).eq.'LIQUID ') then
            phtype='L'
         endif
         if(modelname(1:5).eq.'I2SL ') phtype='Y'
!
! THIS IS READPDB subroutine .... OBSOLETE
!
!         write(*,*)'3E enter phase: ',trim(name1),' ',subord,havedisorder
!         call enter_phase(name1,nsl,knr,const,stoik,modelname,phtype)
         call enter_phase(name1,nsl,knr,const,stoik,modelname,phtype,warning,&
              emodel)
!      write(*,*)'readpdb 9A: ',gx%bmperr
         if(gx%bmperr.ne.0) then
            if(gx%bmperr.eq.4121) then
               if(.not.silent) write(kou,*) &
                    'Phase ',trim(name1),&
                    ' is ambiguous or short for another phase'
            endif
            goto 1000
         endif
! we may need to add things to the phase
         call find_phase_by_name(name1,iph,ics)
         if(gx%bmperr.ne.0) goto 1000
         call get_phase_compset(iph,1,lokph,lokcs)
         if(gx%bmperr.ne.0) goto 1000
! add disordered fraction set if any
         if(havedisorder) then
! this is the suffix for the disordered parameters
            ch1='D'
            if(subord) then
! we should subract the ordered G as disordered
! if 4 or 5 subl in ordered the first 4 are summed in the disordered
! if 2 or 3 subl in ordered the first 2 are summed in the disordered
! jl is stored in fracset%totdis and 1 means subtract ordered as disordered
               jl=1
               if(phlista(lokph)%noofsubl.le.5) nd1=4
               if(phlista(lokph)%noofsubl.le.3) nd1=2
! THIS IS CONFUSED!! never set PHSUBO
! The idea was that FCC/BCC/HCP ordering may be calculated without subtracting
! the ordered set but then
!               write(*,*)'3E Setting bit: ', trim(name1),PHSUBO,permut
!               phlista(lokph)%status1=ibset(phlista(lokph)%status1,PHSUBO)
!               firsteq%phase_varres(lokcs)%status2=&
!                    ibset(firsteq%phase_varres(lokcs)%status2,PHSUBO)
            else
! we should just add the ordered and disordered G (without the config.entropy)
! all ordered sublattices merged into a single one for the disordered phase
               jl=0
               nd1=phlista(lokph)%noofsubl
            endif
            call add_fraction_set(iph,ch1,nd1,jl)
            if(gx%bmperr.ne.0) then
               if(.not.silent) write(kou,*) &
                 '3E Error entering disordered fraction set: ',gx%bmperr
               goto 1000
            endif
! we must set the correct formula unit of the disordered phase
            xxx=zero
            do ll=1,nd1
               xxx=xxx+firsteq%phase_varres(lokcs)%sites(ll)
            enddo
! here to number of sites are stored in the disordered fraction set!!
            firsteq%phase_varres(lokcs)%disfra%fsites=xxx
         endif
! set permutations
         if(permut.eq.1) then
            phlista(lokph)%status1=ibset(phlista(lokph)%status1,PHFORD)
         elseif(permut.eq.2) then
            write(*,*)'3E bcc permutations not implemented'
            gx%bmperr=4399; goto 1000
!            phlista(lokph)%status1=ibset(phlista(lokph)%status1,PHBORD)
         endif
! any typedefs? only magnetic handelled at present
!         call find_phase_by_name(name1,iph,lcs)
!         write(*,*)'readpdb 9X: ',gx%bmperr
!         if(gx%bmperr.ne.0) then
!            if(.not.silent) write(kou,*)'Phase ',name1,' is ambiguous'
!            goto 1000
!         endif
!         lokph=phases(iph)
!      write(*,*)'3E typedefs for ',name1(1:20),lokph,xx
!         phasetypes: do jt=1,xx
!          write(*,*)'3E typedef ',jt,addphasetypedef(jt)
!   character*12, dimension(8), parameter :: addition=&
!        ['IMAGF       ','IMAGB       ','IWMAGF      ','IWMAGB      ',&
!         'DEBYE1      ','EINSTEIN1   ','            ','            ']
!
! THIS READPDB ... not readtdb
         do jt=1,noofadds
            if(add(jt).eq.1) then
! Inden magnetic for FCC
               call add_addrecord(lokph,'N',indenmagnetic)
            elseif(add(jt).eq.2) then
! Inden magnetic for BCC
               call add_addrecord(lokph,'Y',indenmagnetic)
            elseif(add(jt).eq.3) then
! Xiong-Inden magnetic for non-BCC
               call add_addrecord(lokph,'N',xiongmagnetic)
            elseif(add(jt).eq.4) then
! Xiong-Inden magnetic for BCC
               call add_addrecord(lokph,'Y',xiongmagnetic)
!            elseif(add(jt).eq.7) then
! Volume model volmod1 always entered at call enter_phase(...
!               write(*,*)'3E adding volume model ',jt
!               call add_addrecord(lokph,'Y',volmod1)
            else
! no other implemented
               write(*,*)'Addition not implemented: ',addition(add(jt))
               gx%bmperr=4399; goto 1000
            endif
         enddo
      endif condis2
!      write(*,*)'readpdb 9B:',name1,nsl,phtype
!-------------------------------------------------------------------
   case(6) ! PARAMETER --------------------------------------------
! for PDB databases
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
!   PARAMETER G(LIQUID,CR;0)  2.98150E+02  +24339.955-11.420225*T
      if(eolch(longline,nextc)) then
         if(.not.silent) write(kou,*)'Empty line after PARAMETER'
         gx%bmperr=4311; goto 1000
      endif
      ip=nextc
      funname=longline(ip:)
      kp=index(funname,' ')
      cbug=index(funname,'),')
! save position after parameter name in nextc
      if(cbug.gt.0 .and. cbug.lt.kp) then
         nextc=ip+cbug+1
         kp=cbug+1
      else
         nextc=ip+kp
      endif
      funname(kp:)=' '
! extract symbol, normally G or L but TC, BMAGN and others can occur
      lp1=index(funname,'(')
      name1=funname(1:lp1-1)
! handle parameters for disordered part
! NOTE WE SHOULD NEVER HAVE A PARAMETER IDENTIFIER ENDING WITH D ...
      if(funname(lp1-1:lp1-1).eq.'D') then
         fractyp=2
         name1(lp1-1:lp1-1)=' '
      else
         fractyp=1
      endif
      typty=0
! this is kept for compatibility with TDB files generated by TC
      if(name1(1:2).eq.'G ' .or. name1(1:2).eq.'L ') then
         typty=1
      elseif(name1(1:3).eq.'TC ') then
         typty=2
      elseif(name1(1:5).eq.'BMAG ') then
         typty=3
      elseif(name1(1:6).eq.'BMAGN ') then
         typty=3
      endif
! we should handle also other parameter types
      if(typty.eq.0) then
! find the property associated with this symbol
!         write(*,*)'psym1: ',name1(1:len_trim(name1))
         call get_parameter_typty(name1,lokph,typty,fractyp)
         if(gx%bmperr.ne.0) then
            if(.not.silent) write(kou,*) &
                 ' *** Illegal parameter identifier on line: ',nl
            gx%bmperr=0; typty=0
            warning=.TRUE.
         endif
!         write(*,*)'psym2: ',typty,fractyp
      endif
! only fractyp 1 on TDB files until I implemented disordered part
!      fractyp=1
!       write(*,*)'readpdb: PAR',name1,typty
! extract phase name and constituent array
      lp1=index(funname,'(')
      lp2=index(funname,',')
      name2=funname(lp1+1:lp2-1)
!      dispar: if(dodis.eq.1) then
! first check if phase name is a disordered part, if not skip
! then change phase name to ordered phase and set fractyp=2
! and add a suffix D to parameter symbol
!         do jl=1,disparttc
!            if(name2.eq.dispartph(jl)) goto 710
!         enddo
! not disordered phase, skip this parameter
!         goto 100
!-----------------------
!710      continue
!         write(*,*)'Entering disordered parameter to: ',thisdis,jl
!         thisdis=jl
!         write(*,*)'Entering disordered parameter to: ',ordpartph(thisdis)
!         write(*,*)'3E ',longline(1:len_trim(longline))
!         name2=ordpartph(jl)
!         fractyp=2
!      endif dispar
      call find_phase_by_name_exact(name2,jph,kkk)
!      write(*,*)'readpdb 19: ',jph,gx%bmperr,name2
      if(gx%bmperr.ne.0) then
!         write(*,*)'Skipping parameter due to phase: ',name2
         gx%bmperr=0; goto 100
!         goto 1000
      endif
! extract constituent array, remove final ) and decode
! constituent names can be very long ....
      lokph=phases(jph)
      name4=funname(lp2+1:)
! find terminating )
      lp1=index(name4,')')
      if(lp1.le.0) then
         if(.not.silent) write(kou,*) &
              'Possible error in constituent array? ',name4,', line:',nl
         warning=.TRUE.
         goto 100
      else
         name4(lp1:)=' '
      endif
297   continue
!
      call decode_constarr(lokph,name4,nsl,endm,nint,lint,ideg)
      if(ocv()) write(*,303)'readpdb 303: ',name4(1:len_trim(name4)),&
           nsl,endm(1),endm(2),nint,((lint(ip,jp),ip=1,2),jp=1,nint)
303   format(a,a,2i4,2x,2i3,' : ',3(2i3,2x))
      if(gx%bmperr.ne.0) then
! error here can mean parameter with un-selected constituent, i.e. no error
!         write(*,*)'3E: decode',ionliq,tdbv,nsl,gx%bmperr
         if(ionliq .and. tdbv.eq.1 .and. nsl.eq.1) then
! handle parameters in ionic liquids with only neutrals in second sublattice
! in TC one can have no constituent there or an arbitrary constituent,
! in OC the constituent in sublattice 1 must be a *
            nsl=2
            endm(2)=endm(1)
            endm(1)=-99
! shift any interaction from sublattice 1 to 2
            do ip=1,nint
!               write(*,*)'3E lint: ',lint(1,ip),lint(2,ip)
               lint(2,ip)=2
            enddo
            if(ocv()) write(*,303)'modif endmem: ',name4(1:len_trim(name4)),&
                 nsl,endm(1),endm(2),nint,((lint(ip,jp),ip=1,2),jp=1,nint)
            gx%bmperr=0
         else
            if(ocv()) write(*,*)'Skipping parameter: ',name4(1:len_trim(name4))
            gx%bmperr=0; goto 100
!         write(*,*)'readpdb error: ',gx%bmperr,name4
!         goto 1000
         endif
      endif
! check number of sublattices correct
      if(.not.ionliq) then
         lokcs=phlista(lokph)%linktocs(1)
         if(fractyp.eq.1) then
            if(phlista(lokph)%noofsubl.ne.nsl) then
               write(*,*)'Wrong number of sublattices 1',&
                    phlista(lokph)%noofsubl,nsl
               gx%bmperr=4399; goto 1000
            endif
         else
            if(firsteq%phase_varres(lokcs)%disfra%ndd.ne.nsl) then
               write(*,*)'Wrong number of sublattices 2',&
                    firsteq%phase_varres(lokcs)%disfra%ndd,nsl
               gx%bmperr=4399; goto 1000
            endif
         endif
      endif
!      if(nint.gt.1) then
! lint(1,1) is species of first, lint(1,2) in second interaction
!          write(*,305)'readpdb 305: ',endm(1),nint,lint(2,1),lint(2,2)
!      endif
305    format(a,5i4)
!---------------- encode function
!      if(dodis.eq.1) write(*,*)'We are here 1'
      ip=0
      jp=0
400    continue
      ip=ip+1
405    continue
      ch1=funname(ip:ip)
! accept the first 8 letters and numbers of phase name
      if((ch1.ge.'A' .and. ch1.le.'Z') .or. &
         (ch1.ge.'0' .and. ch1.le.'9')) goto 400
      if(ch1.ne.' ') then
         funname(ip:)=funname(ip+1:)
         jp=jp+1
         if(jp.lt.8) goto 405
         funname(ip+1:)=' '
      endif
      funname='_'//funname
!-------------------------------------------------
! now read the function, start from position nextc
      longline=longline(nextc:)
!410    continue
      jp=len_trim(longline)
      if(longline(jp:jp).ne.'!') then
         if(.not.silent) write(kou,410)nl,ip,longline(1:ip)
410      format('Error, parameter line not ending with !',2i5/a)
         gx%bmperr=4312; goto 1000
      endif
! extract reference if any
! NOTE: a legal ending is ;,,,!
      refx='none'
      kp=jp-1
      do while(longline(kp:kp).ne.';')
         kp=kp-1
         if(kp.lt.1) then
! illegal termination of function in TDB file
            if(.not.silent) write(kou,417)nl
417 format('No final ; of function in TDB file, around line: ',i5)
            gx%bmperr=4013; goto 1000
         endif
      enddo
      kp=kp+2
! longline(kp:kp) is character after "; " or ";," 
! next is upper temperature limit or , meaning default.  We have a "!" at end
430   continue
      if(eolch(longline,kp)) continue
      if(longline(kp:kp).eq.',') then
         kp=kp+1
      elseif(longline(kp:kp).eq.'!') then
         goto 433
      else
!    ; 6000 N 91DIN !
!   kp=^                 => index(...,' ')=5; kp=kp+4
         kp=kp+index(longline(kp:),' ')-1
      endif
! next is N or ,
      if(eolch(longline,kp)) continue
      if(longline(kp:kp).ne.'!') then
         kp=kp+1
      endif
      if(eolch(longline,kp)) continue
      if(kp.lt.jp) then
         refx=longline(kp:jp-1)
         call capson(refx)
      else
         refx=' '
      endif
! ------------------- we found the reference, continue with the expression
433   continue
! replace any # by ' '
412   continue
      jss=index(longline(1:jp),'#')
      if(jss.gt.0) then
         longline(jss:jss)=' '
         goto 412
      endif
!      write(*,*)'3E Entering function 2: ',funname,len_trim(longline)
!      lrot=0
!      call store_tpfun(funname,longline,lrot,.TRUE.)
      call store_tpfun(funname,longline,lrot,rewindx)
!          write(*,17)lokph,typty,nsl,lrot,(endm(i),i=1,nsl)
17 format('readpdb 17: ',4i3,5x,10i3)
!         write(*,404)'readpdb entpar: ',refx,fractyp,nint,ideg
404   format(a,a,i3,2x,10i3)
      if(gx%bmperr.ne.0) then
         if(.not.silent) write(kou,*)'Error set: ',gx%bmperr,lrot,' ',&
              funname(1:len_trim(funname)),' around line: ',nl
         goto 1000
      else
! THIS IS IN READPDB
!         if(dodis.eq.1) write(*,*)'We are here 2'
         call enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,&
              lrot,refx)
         if(ocv()) write(*,407)'Entered parameter: ',lokph,typty,gx%bmperr
407      format(a,3i5)
         if(gx%bmperr.ne.0) then
! error entering parameter, not fatal
            if(dodis.eq.1 .and. .not.silent) &
                 write(*,408)'3E parameter warning:',gx%bmperr,nl,&
                 funname(1:40)
408         format(a,i6,' line ',i5,': ',a)
            if(.not.(gx%bmperr.ne.4096 .or. gx%bmperr.ne.4066)) goto 1000
! ignore error 4096 meaning "no such constituent" or "... in a sublattice"
!            write(*,*)'readpdb entparerr: ',gx%bmperr,' >',&
!                 funname(1:len_trim(funname))
!            if(gx%bmperr.eq.7778 .and. .not.silent) &
            if(gx%bmperr.eq.4153 .and. .not.silent) &
                 write(*,*)'3E Error 7778 at line: ',nl
            gx%bmperr=0
!         elseif(dodis.eq.1) then
!            write(*,*)'Disordered parameter should be entered ok'
         else
            enteredpar=enteredpar+1
         endif
      endif
      if(gx%bmperr.ne.0 .and. .not.silent) write(*,*)'3E error 1: ',gx%bmperr
!      if(verbose) write(*,*)'Entered parameter: ',trim(name1)
!------------------------------------------------------------------
   case(7) !SPECIAL
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
      write(*,*)'SPECIAL not yet implemented'
      goto 100
! 
! TYPE_DEFINITION & GES A_P_D BCC_A2 MAGNETIC  -1.0    4.00000E-01 !
      nytypedef=nooftypedefs+1
      nooftypedefs=nytypedef
      typedefchar(nytypedef)=longline(nextc+1:nextc+1)
      ip=nextc+3
      newtypedef: if(longline(ip:ip+2).eq.'SEQ') then
         typedefaction(nytypedef)=100
      else
         km=index(longline,' MAGNETIC ')
         magnetic: if(km.gt.0) then
            ip=km+9
!73           format(a,i3,' "',a,'"')
            call getrel(longline,ip,xxx)
            if(buperr.ne.0) then
               gx%bmperr=buperr; goto 1000
            endif
! this can be -1 for BCC or -3 for FCC, HCP and other phases
            typedefaction(nytypedef)=int(xxx)
         else
            km=index(longline,' DIS_PART ')
            if(km.gt.0) then
! disordered part, several checks
               disparttc=disparttc+1
! find the ordered phase name, we have to go backwrds from km
               ip=km-1
81             continue
               if(longline(ip:ip).eq.' ') then
!                  ordpartph(disparttc)=' '
                  ordpartph(disparttc)=longline(ip+1:km)
               else
                  ip=ip-1
                  goto 81
               endif
! extract the disordered part phase name
               ip=index(longline(km+2:),' ')
               dispartph(disparttc)=longline(km+2+ip:)
! find the end of phase name, a space or a , there is always a space ...
               ip=index(dispartph(disparttc),' ')
               km=index(dispartph(disparttc),',')
               if(km.gt.0 .and. km.lt.ip) ip=km
!               if(ip.le.0) ip=1
               dispartph(disparttc)(ip:)=' '
! ths is reading PDB file ...
               if(.not.silent) write(kou,82) &
                    disparttc,trim(ordpartph(disparttc)),&
                    trim(dispartph(disparttc)),orddistyp(disparttc)
!                    longline(1:len_trim(longline))
82             format('Found a type_def DIS_PART:',i2,1x,a,1x,a)
! if the disordered part phase already entered give advice
               call find_phase_by_name(dispartph(disparttc),iph,ics)
               if(gx%bmperr.ne.0) then
                  gx%bmperr=0
               else
                  if(.not.silent) write(kou,83)dispartph(disparttc)
83                format(' *** Warning, the disordered phase is already',&
                       ' entered ***'/' Please rearrange the TDB file so',&
                       ' this TYPE_DEF comes before'/&
                       ' the PHASE keyword for the disordered phase: ',a/&
                       ' *** The disordordered part ignored ***')
                  disparttc=disparttc-1
                  warning=.TRUE.
               endif
            else
               typedefaction(nytypedef)=99
               if(.not.silent) &
                    write(kou,87)nl,longline(1:min(78,len_trim(longline)))
87             format('Skipping this TYPE_DEFINITION on line ',i5,':'/a)
               warning=.TRUE.
!               write(*,*)' WARNING SET TRUE <<<<<<<<<<<<<<<<<<<<<<<<<<<'
            endif
         endif magnetic
      endif newtypedef
!---------------------------------------------------------------------
! BIBLIOGRAPHY
   case(8) 
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
!   REF283  'Alan Dinsdale, SGTE Data for Pure Elements,
!          Calphad Vol 15(1991) p 317-425,
!          also in NPL Report DMA(A)195 Rev. August 1990'
      ip=1
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'Empty reference line',nl
         gx%bmperr=4313; goto 1000
      endif
      if(longline(ip:ip).eq.'!') then
!         write(*,*)'No references at all'
         goto 100
      endif
!      write(*,*)'list_of_references text length: ',len_trim(longline),ip
! some reference lists like those from SSUB has no single quotes
      kp=index(longline(ip:),"'")
      citationmarks: if(kp.gt.0) then
775      continue
! reference symbol is refx; reference text in reftext
         refx=longline(ip:ip+kp-2)
         if(longline(ip+kp:ip+kp).eq."'") then
! two ' after each other, a dummy reference
            reftext=' '
            ip=ip+kp+1
            kkk=1
!            write(*,*)'dummy: ',refx,' next >',longline(ip:ip+20),'<'
         else
            jp=ip+kp+1+index(longline(ip+kp+1:),"'")
            reftext=longline(ip+kp:jp-2)
            ip=jp
! when all works replace multiple spaces by a single one in reftext
            kkk=len_trim(reftext)
            kp=index(reftext(1:kkk),'  ')
            do while(kp.gt.0)
               reftext(kp:)=reftext(kp+1:)
               kkk=kkk-1
               kp=index(reftext(1:kkk),'  ')
            enddo
         endif
!         write(*,776)refx,nrefs,ip,jp,reftext(1:kkk)
776      format('Reference: ',a,3i5/a)
! this will not create bibliographic references that has not been referenced
         call tdbrefs(refx,reftext(1:kkk),1,ix)
         nrefs=nrefs+1
!         write(*,*)'added biblio ',refx,'>',longline(ip-5:ip+5),'<'
         if(eolch(longline,ip)) then
            gx%bmperr=4313; goto 1000
         endif
         if(longline(ip:ip).ne.'!') then
            kp=index(longline(ip:),"'")
            goto 775
         endif
      else
! references without citation marks
! ip is at the start of the reference id, look for space
         if(.not.silent) write(kou,*) &
              'Cannot handle references without citation marks',nl
         gx%bmperr=4314; goto 1000
      endif citationmarks
777   continue
!      write(*,*)'Read ',nrefs,' references, ending at',nl
!------------------------------------------------------------------
   case(9) ! TABLE_OF_MODELS
      write(*,*)'Table of models not implemented'
!------------------------------------------------------------------
   case(10) ! TABLE_OF_IDENTIFIERS
      write(*,*)'Table of idenifiers not implemented'
!------------------------------------------------------------------
   case(11) ! DATABASE_INFORMATION
      if(.not.silent) write(kou,*)'Cannot handle DATABASE_INFORMATION at ',nl
!      warning=.TRUE.
! skip lines until !
      do while(index(line,'!').le.0)
         read(21,110)line
         nl=nl+1
         call replacetab(line,nl)
      enddo
!----------------------------------------------------------------
   case(12) ! ASSESSED_SYSTEMS
      if(.not.silent) write(kou,*) &
           'Cannot handle ASSESSED_SYSTEMS ending at ',nl
!      warning=.TRUE.
! skip lines until !
      do while(index(line,'!').le.0)
         read(21,110)line
         nl=nl+1
         call replacetab(line,nl)
      enddo
!------------------------------------------------------------------
   case(13) ! DEFAULTS
      write(*,*)'Cannot handle defaults'
!------------------------------------------------------------------
   case(14) ! VERSION
780   continue
      write(kou,*)'Found VERSION keyword: ',trim(longline)
!----------------------------------------
   case(15) ! INCLUDE_FILE
      write(*,*)'INCLUDE_FILE not implemented yet'
!----------------------------------------
   case(16) ! CHECKSUM
      write(*,*)'CHECKSUM not implemented yet'
   end select
!=====================================================
   if(gx%bmperr.ne.0 .and. .not.silent) then
      write(kou,711)gx%bmperr,nl,trim(line)
711   format('3E error: ',i5,' around line ',i7,': '/a)
! this error means reference error
      if(gx%bmperr.eq.4154) gx%bmperr=0
   endif
! look for next KEYWORD
   goto 100
!--------------------------------------------------------
!----- reading functions at the end from a PDB file
800 continue
!   if(eolch(line,nextc)) then
!      if(.not.silent) write(kou,*) &
!           'Function name must be on same line as FUNCTION'
!      gx%bmperr=4315; goto 1000
!   endif
!   ipp=nextc+index(line(nextc:),' ')
   ip=nextc
   call getext(longline,ip,1,name1,' ',ix)
!         write(*,18)'function >',name1,'< ',nextc,ipp
!18       format(a,a,a,2i4)
! old code
!   longline=' '
!   longline=line(ipp:)
!810 continue
!   jp=max(len_trim(longline),1)
!   write(*,811)jp,longline(jp:jp),longline(1:jp)
!811 format('3E ll: ',i3,' "',a1,'" ',a)
!   if(longline(jp:jp).eq.'!') then
! replace # by ' '
   jp=max(len_trim(longline),1)
820 continue
   jss=index(longline(1:jp),'#')
   if(jss.gt.0) then
      longline(jss:jss)=' '
      goto 820
   endif
! check if function is entered as undefined, exact match of name required
!   write(*,*)'TPFUN: ',trim(name1),':  ',notent,ip,longline(ip:ip+30)
   call find_tpfun_by_name_exact(name1,nr,notent)
   if(gx%bmperr.eq.0) then
!      write(*,*)'Found TPfunction: ',trim(name1),notent
      if(notent) then
! entering a function may add new unentered functions ... last argument TRUE
!            write(*,*)'3E Entering function 3: ',name1,len_trim(longline)
!            lrot=0
!         call store_tpfun(name1,longline(ip:),lrot,.TRUE.)
         call store_tpfun(name1,longline(ip:),lrot,rewindx)
         if(gx%bmperr.ne.0) then
! one may have error here
            if(.not.silent) write(kou,*)'Failed entering function: ',name1
            goto 1000
         endif
!         write(*,*)'Entered missing function: ',name1
         nofunent=nofunent+1
!      elseif
! this is the PDB reader
! function and its expression entered, check if at the same rewind!
       endif
   else
! reset error code
      gx%bmperr=0
   endif
!   else
!830   continue
!      nl=nl+1
!      read(21,110)line
!            write(kou,101)'readpdb 2: ',nl,line(1:40)
! skip lines with a $ in first position
!      if(line(1:1).eq.'$')goto 830
!      call replacetab(line,nl)
!      longline=longline(1:jp)//line
!      goto 810
!   endif
   goto 100
!   endif barafun
!---------------------------------------------------------
! We have now read all
!--------------------------------------------------------
1000 continue
   if(warning) then
1001  continue
! if silent set ignore warnings
      if(.not.silent) then
         write(kou,1003)
1003     format(/'There were warnings, check them carefully'/&
              'and press RETURN if you wish to continue.')
         read(kiu,1004)ch1
1004     format(a)
!         if(ch1.eq.'N') stop 'warnings reading database'
!         if(ch1.ne.'Y') then
!            write(kou,*)'Please answer Y or N'
!            goto 1001
!         endif
      endif
   endif
! 4057
   if(gx%bmperr.ge.4000 .and. gx%bmperr.le.4399) then
      write(*,*)trim(bmperrmess(gx%bmperr))
   endif
!   write(*,*)'3E At label 1000'
   if(buperr.ne.0 .or. gx%bmperr.ne.0) then
      if(gx%bmperr.eq.0) gx%bmperr=buperr
      if(.not.silent) write(kou,1002)gx%bmperr,nl
1002   format('Error ',i5,', occured at PDB file line ',i7)
!      write(*,*)'Do you want to continue at your own risk anyway?'
!      read(*,1008)ch1
!1008  format(a)
!      if(ch1.eq.'Y') then
!         write(*,*)'Now any kind of error may occur .... '
!         buperr=0
!         gx%bmperr=0
!         goto 100
!      endif
   endif
   close(21)
! read numbers, value after / is maximum
! endmember, interactions, property,
! tpfuns, composition sets, equilibria
! state variable functions, references, additions
   if(ocv()) write(*,1007)noofel,maxel,noofsp,maxsp,noofph,maxph,&
        noofem,100000,noofint,100000,noofprop,100000,&
        notpf(),maxtpf,highcs,2*maxph,eqfree-1,maxeq,&
        nsvfun,maxsvfun,reffree-1,maxrefs,addrecs,csfree-1
1007 format('Created records for elements, species, phases: ',2x,&
          3(i4,'/',i4,1x)/&
          'end members, interactions, properties: ',10x,&
          3(i4,'/',i4,1x)/&
          'TP-funs, max and free composition sets, equilibria: ',10x,&
          3(i4,'/',i4,1x)/&
          'state variable functions, references, additions: ',&
          3(i4,'/',i4,1x)/)
   return
1010 continue
   if(.not.silent) write(kou,*)'I/O error opening file: ',gx%bmperr
   return
!-----------------------------------------------------
! end of file found, act differently if reading functions
2000 continue
!   rewind: if(dodis.eq.0 .and. disparttc.gt.0) then
! rewind to read disordred parts
!      if(.not.silent) &
!           write(kou,*)'Rewind to read disordered parts of phases: ',disparttc
!      rewind(21)
!      dodis=1
!      nl=0
!      goto 100
!   elseif(.not.onlyfun) then
   rewindx=rewindx+1
   rewind: if(.not.onlyfun) then
! rewind to read referenced functions
!      dodis=2
!      write(*,*)'Rewinding to read functions first time'
      rewind(21)
      onlyfun=.TRUE.
      nofunent=0
!      write(*,2002)gx%bmperr
2002 format('Found end-of-file, rewind to find functions',i5)
      nl=0
      goto 100
!   elseif(norew.lt.3) then
! rewind if there were functions entered last time
!      rewind(21)
!      norew=norew+1
!    write(*,*)'Found functions: ',nofunent,' rewinding again',norew,gx%bmperr
!      if(newfun.gt.0) then
!          write(*,*)'Read ',newfun+nfail,' functions, entered ',newfun,&
!               ' rewinding ',norew
!         newfun=0
!      nofunent=0
!      nl=0
!      goto 100
   elseif(nofunent.gt.0) then
! rewind if there were functions entered last time
      rewind(21)
      norew=norew+1
!    write(*,*)'Found functions: ',nofunent,' rewinding again',norew,gx%bmperr
!      if(newfun.gt.0) then
!          write(*,*)'Read ',newfun+nfail,' functions, entered ',newfun,&
!               ' rewinding ',norew
!         newfun=0
      nofunent=0
      nl=0
      goto 100
   else
! check if there are any unentered functions
      call list_unentered_funs(kou,nr)
      if(nr.gt.0) then
         if(.not.silent) write(kou,*)'Number of missing function: ',nr
         gx%bmperr=4186
      endif
! check if any function not entered
      onlyfun=.FALSE.
   endif rewind
   goto 1000
! end of file while looking for ! terminating a keyword
2200 continue
   if(.not.silent) write(kou,2210)nl,longline(1:72)
2210 format('End of file at ',i5,' looking for end of keyword:'/a)
   gx%bmperr=4316
   goto 1000
 end subroutine readpdb

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine write_pdbformat
!\begin{verbatim}
 subroutine write_pdbformat(unit)
! write a PDB database
   implicit none
   integer unit
!\end{verbatim}
   write(*,*)'PDB format output not yet implemented'
1000 continue
   return
 end subroutine write_pdbformat

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine checkdb
!\begin{verbatim}
 subroutine checkdb(filename,ext,nel,selel)
! checking a TDB/PDB file exists and return the elements
   implicit none
   integer nel
   character filename*(*),ext*4,selel(*)*2
!\end{verbatim} %+
   character line*256,ext2*4
   integer ipp,nl,kk
!
   ext2=ext
   call capson(ext2)
   if(.not.(index(filename,ext).gt.0 &
       .or. index(filename,ext2).gt.0)) then
! no extention provided
      filename(len_trim(filename)+1:)=ext2
   endif
   open(21,file=filename,access='sequential',form='formatted',&
        err=1010,iostat=gx%bmperr,status='old')
! if first line of file is "$OCVERSION ..." the text is displayed once
   read(21,110)line
   if(line(1:11).eq.'$OCVERSION ') then
      write(kou,117)trim(line(12:))
117   format(/'TDB file id: ',a/)
   endif
   rewind(21)
! just check for ELEMENT keywords
! return here to look for a new keyword, end-of-file OK here
   nl=0
   nel=0
100 continue
   read(21,110,end=2000)line
110 format(a)
   nl=nl+1
! One should remove TAB characters !! ??
   call replacetab(line,ipp)
   ipp=1
   if(eolch(line,ipp)) goto 100
   if(line(ipp:ipp).eq.'$') goto 100
! look for ELEMENT keyword, ipp=1
   ipp=istdbkeyword(line,kk)
   if(ipp.ne.1) goto 100
!
! ignore /- and VA
   if(eolch(line,kk)) goto 100
   if(line(kk:kk+1).eq.'/-' .or. line(kk:kk+1).eq.'VA') goto 100
   nel=nel+1
   selel(nel)=line(kk:kk+1)
!      write(*,111)nl,line(1:20)
!111   format('Read line ',i5,': ',a)
   goto 100
!---------
1000 continue
   return
! error
1010 continue
   goto 1000
! end of file
2000 continue
   close(21)
   goto 1000
   return
 end subroutine checkdb

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine checkdb2
!\begin{verbatim}
 subroutine checkdb2(filename,ext,nel,selel)
! checking a TDB/PDB file exists and return the elements
! It also writes 15 lines from any "DATABASE_INFO" in the file
   implicit none
   integer nel
   character filename*(*),ext*4,selel(*)*2
!\end{verbatim}
   character line*256,ext2*4
   integer ipp,nl,kk,dbinfo
!
   ext2=ext
   dbinfo=0
   call capson(ext2)
   if(.not.(index(filename,ext).gt.0 &
       .or. index(filename,ext2).gt.0)) then
! no extention provided
      filename(len_trim(filename)+1:)=ext2
   endif
   open(21,file=filename,access='sequential',form='formatted',&
        err=1010,iostat=gx%bmperr,status='old')
! if first line of file is "$OCVERSION ..." the text is displayed once
   read(21,110)line
   if(line(1:11).eq.'$OCVERSION ') then
      write(kou,117)trim(line(12:))
117   format(/'TDB file id: ',a/)
   endif
   rewind(21)
! just check for ELEMENT and DATABASE_INFO keywords
! return here to look for a new keyword, end-of-file OK here
   nl=0
   nel=0
100 continue
   read(21,110,end=2000)line
110 format(a)
   nl=nl+1
! One should remove TAB characters !! ??
   call replacetab(line,ipp)
   ipp=1
   if(eolch(line,ipp)) goto 100
   if(line(ipp:ipp).eq.'$') goto 100
! look for ELEMENT keyword, ipp=1
   ipp=istdbkeyword(line,kk)
   if(ipp.eq.11 .and. dbinfo.eq.0) then
! DATABASE_INFORMATION keyword, ipp=11
!      write(*,*)'3E at line ',nl,': ',trim(line)
      dbinfo=1
      write(kou,200)trim(line)
200   format(/'This database has information to users, please read carefully'/a)
      do while(index(line,'!').le.0)
         read(21,110)line
         write(kou,110)trim(line)
      enddo
      write(kou,*)
   endif
   if(ipp.ne.1) goto 100
!
! ignore /- and VA
   if(eolch(line,kk)) goto 100
   if(line(kk:kk+1).eq.'/-' .or. line(kk:kk+1).eq.'VA') goto 100
   nel=nel+1
   selel(nel)=line(kk:kk+1)
!      write(*,111)nl,line(1:20)
!111   format('Read line ',i5,': ',a)
   goto 100
!---------
1000 continue
   return
! error
1010 continue
   goto 1000
! end of file
2000 continue
   close(21)
   goto 1000
   return
 end subroutine checkdb2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
! moved from gtp3C

!\addtotable subroutine save_datformat
!\begin{verbatim}
 subroutine save_datformat(filename,version,kod,ceq)
! writes a SOLGASMIX DAT format file. not (ever?) finished
   implicit none
   integer kod
   character filename*(*),version*(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer ntpf,last,i1,i2,i3,npows,lut,ip,jp,nstoi,lokph,isp,f1,nphstoi,nphmix
   integer, dimension(:), allocatable :: ncon,phmix,phstoi,estoi,endx,endy
   integer nelectrons,lokcs,nsubl,isubl,mphstoi,k1,lcase,mult(10),check
   integer cation,anion,firstcation,ilevel,intconst(9),intconstx(9),ideg
   integer lokdis,wildp,havemag,offset,nn,paratyp,maxideg,wildloop,intconsty(9)
! iset is sets of interaction constituents generated for wildcards
   integer, allocatable, dimension(:,:) :: iset
   logical logok,nogas,ionliq,wildcard,iliqwild,excessparam,skipfc,magloop
   character ch1*1,line*16,powers*80,model*24,constext*80,elsym*2
! overflow in text line before label 210
   character text*2048
   character date*12,hour*12,phunique*4,phdummy*4
! this is for mixture phases that have names with are not unique first 4 chars
   character phcharged(50)*24
   integer phchargedx(50),nnn,noelel,donotincrement
! this is to check we have correct number of endmembers
   integer end1mem,end2mem,lineb
   type(gtp_tpfun2dat), dimension(:), allocatable :: tpfc
   type(gtp_endmember), pointer :: endmember,nextcation,samecation
   double precision, allocatable, dimension(:) :: constcomp,constcompiliq
   double precision valency(9),ccc,cationval,factor,disfactor,aff,partc,parbm
   double precision extcpar(0:7),exbmpar(0:7),xxx
   double precision, parameter :: maxcc=1.0D2
   TYPE(gtp_phase_add), pointer :: addrec
! These are to handle problems ....
   integer warnings,decimals,missend(9),thisend(9),www,xnooffr(0:9)
! indices for excess parameters
   integer exix(500),lastix,firstix
! we must probably create a stack for excess parameters
   type intstack
      type(gtp_interaction), pointer :: intlink
   end type intstack
   type(intstack), dimension(5) :: saveint
   type(gtp_interaction), pointer :: intparam
   type(gtp_phase_varres), pointer :: varres,disvarres,fedup
   type(gtp_property), pointer :: property,nextprop,savedproperty
!
   inquire(file=filename,exist=logok)
   date=' '
   hour=' '
   if(logok) then
      line=' '
      last=len(line)
      call gparcdx('File exists, overwrite?',line,last,1,ch1,'N','SAVE solgas')
      if(ch1.ne.'Y') then
         write(*,*)'Better luck next time!'
         goto 1000
      endif
   endif
   ntpf=freetpfun-1
! allocate coefficient arrays for all TP functions (incl parameters) and 5 extra
   write(*,*)'3C Allocating for TP functions: ',ntpf+5
   allocate(tpfc(ntpf+5))
!   write(*,*)'TPfuns and parameters: ',ntpf
! in this call all tpfuns are converted to arrays of coefficients
! each tpfc(i) represent TPfunction i (note parameters are also TP functions!)
! text returns the powers of T used
   call tpfun2coef(tpfc,ntpf,npows,powers)
   if(gx%bmperr.ne.0) goto 1000
! listing of all tpfuns
!   write(*,*)trim(text)
!   do i1=3,ntpf
!      write(text,76)i1
!76    format(' ranges, TP function number s ',i5,' *****************')
!      call list_tpascoef(kou,text,i1,npows,tpfc)
!   enddo
   if(noofph.le.0) then
      write(*,*)'3C No data so nothing to save'
      goto 1000
   endif
   warnings=0
   disfactor=one
   lut=21
   open(lut,file=filename,access='sequential',status='unknown')
   write(*,*)'Writing on file: ',trim(filename)
   text='System'
   ip=len_trim(text)+2
   do i1=1,noofel
      text(ip:)=trim(ellista(elements(i1))%symbol)//'-'
      ip=len_trim(text)+1
   enddo
   call date_and_time(date,hour)
   text(ip-1:)=' generated from TDB file by OC '//version//' '//date(1:4)//&
        '.'//date(5:6)//'.'//date(7:8)//' : '//hour(1:2)//'.'//hour(3:4)
   write(lut,100)trim(text)
99 format(a)
100 format(1x,a)
!------------------- we have to sort the phases as SOLGASMIX wants
! and list constituents in gas, mixtures, stoichiometric 
   allocate(ncon(noofph))
   allocate(endx(noofph))
   allocate(endy(noofph))
   allocate(phmix(noofph))
   allocate(phstoi(noofph))
   allocate(estoi(noofph))
! SOLGASMIX phase names must start with 4 unique letters, when TDB files
! has phases with same first 4 characters add a prefix
   phunique='P000'
   ncon=0
   endx=0
   endy=0
   phmix=0
   phstoi=0
   estoi=0
   nelectrons=0
! check for gas phase, it must be the first phase and name start with GAS
   lokph=phases(1)
   nogas=.true.
   if(phlista(lokph)%name(1:3).eq.'GAS') nogas=.false.
   nphmix=0
   nphstoi=0
   phchargedx=0
   nnn=0
   phloop1: do i1=1,noofph
      lokph=phasetuple(i1)%lokph
      if(ceq%phase_varres(phlista(lokph)%linktocs(1))%phstate.eq.PHSUS) then
! skip phases with suspended default composition set
         write(*,*)'3C skipping phase loop 1: ',phlista(lokph)%name
         cycle phloop1
      endif
      skipfc=.false.
      if(phlista(lokph)%nooffs.gt.1) then
         lokcs=phlista(lokph)%linktocs(1)
!         write(*,105)trim(phlista(lokph)%name),&
!              ceq%phase_varres(lokcs)%disfra%ndd
105      format('Phase ',a,' only disordered saved ',i3)
         skipfc=.true.
      endif
      ncon(i1)=phlista(lokph)%tnooffr-phlista(lokph)%noofsubl
      if(ncon(i1).eq.0) then
! this phase has fixed composition
         nphstoi=nphstoi+1
         phstoi(nphstoi)=i1
! a stoichiometric phase cannot have a charge ...
      else
         nphmix=nphmix+1
         phmix(nphmix)=i1
         if(btest(phlista(lokph)%status1,PHEXCB)) then
! phases with electrones must have the same name for the e(...) as below ...
            nelectrons=nelectrons+1
            estoi(i1)=-noofel-nelectrons
         endif
! prepare a dummy prefix
         phdummy=phlista(lokph)%name(1:4)
         jp=0
         dupname: do i3=1,noofph
            if(i3.ne.lokph .and. phdummy.eq.phlista(i3)%name(1:4)) then
               jp=1; exit dupname
            endif
         enddo dupname
         if(jp.gt.0) then
! we must increment phunique even if phase is not ionic !!
            call incunique(phunique)
! save in phasecharged only if estoi(i1) nonzero
            if(estoi(i1).ne.0) then
               nnn=nnn+1
               if(nnn.gt.50) then
                  write(*,*)'3C too many phases that has name change'
                  gx%bmperr=4399; goto 1000
               endif
               phchargedx(nnn)=i1
               phcharged(nnn)=phunique//'_'//phlista(lokph)%name
               write(*,*)'3C modified charged phase name: ',trim(phcharged(nnn))
            endif
         endif
! should ncon be the number of endmembers?? YES
! NOTE for ionic liquid with neutrals the DAT format requires that the neutrals
! are repeated for each cation, thus the same equation here!!
! if skipfc is TRUE only for disordered fraction set
         i3=1
         if(skipfc) then
            lokcs=phlista(lokph)%linktocs(1)
            varres=>ceq%phase_varres(lokcs)
! this is to check how the ordered phase constituents
!            ip=0
!            do i2=1,phlista(lokph)%noofsubl
!               write(*,*)'3CA: ',lokph,(phlista(lokph)%constitlist(ip+nn),&
!                    nn=1,phlista(lokph)%nooffr(i2))
!               ip=ip+phlista(lokph)%nooffr(i2)
!               i3=i3*phlista(lokph)%nooffr(i2)
!            enddo
!            write(*,*)'3C number of endmembers 1: ',lokph,i3
            disvarres=>ceq%phase_varres(varres%disfra%varreslink)
!            write(*,*)'3ZZ: ',disvarres%sites(1),disvarres%sites(2)
! there must be a disfra record, 
! the number of sublattices and constituents in each sublattice found there
            ip=0
            i3=1
            nsubl=varres%disfra%ndd
            do i2=1,varres%disfra%ndd
!               write(*,*)'3ZB: ',varres%disfra%nooffr(2),&
!                    (varres%disfra%splink(nn+ip),&
!                    nn=1,varres%disfra%nooffr(i2))
!               ip=ip+varres%disfra%nooffr(i2)
               i3=i3*varres%disfra%nooffr(i2)
            enddo
!            write(*,*)'3C number of endmembers 2: ',lokph,i3
            disfactor=varres%disfra%fsites
!            write(*,*)'3ZC factor: ',disfactor,varres%disfra%latd
         else
            nsubl=phlista(lokph)%noofsubl
            do i2=1,phlista(lokph)%noofsubl
               i3=i3*phlista(lokph)%nooffr(i2)
            enddo
            disfactor=one
         endif
!         write(*,*)'3C nonsuspended phase constituents: ',i1,i3
         ncon(i1)=i3
! for check at the end
         endx(i1)=ncon(i1)
      endif
   enddo phloop1
! now can we write the line with overall phase information ... suck
   ip=1
   noelel=noofel+nelectrons
   write(text(ip:),110)noelel
   ip=len_trim(text)+1
! number of mixture phases and for each mixture the number of endmembers
! if nogas is TRUE add a phase with zero endmembers first
   if(nogas) then
      write(text(ip:),109)nphmix+1,0
      ip=len_trim(text)+1
109   format(2i4)
   else
      write(text(ip:),110)nphmix
      ip=len_trim(text)+1
110   format(i4)
112   format(i5)
   endif
   ph1: do i1=1,noofph
      lokph=phasetuple(i1)%lokph
      if(ceq%phase_varres(phlista(lokph)%linktocs(1))%phstate.eq.PHSUS) then
!         write(*,*)'3C skipping phase loop 2: ',phlista(lokph)%name
         cycle ph1
      endif
! Write the number of constituents in mixures (including gas if present)      
!      write(*,*)'3C mixture constituents: ',i1,ncon(i1)
      if(ncon(i1).gt.0) then
         write(text(ip:),112)ncon(i1)
         ip=len_trim(text)+1
         if(ip.gt.72) then
!            write(lut,100)trim(text)
! According to Ted
            write(lut,99)trim(text)
            ip=1
         endif
      endif
   enddo ph1
! finally the number of stoichiometric phases using i5
   write(text(ip:),112)nphstoi
! NOTE format 100 adds an initial space on the line
!   write(lut,100)trim(text)
! According to Ted
   write(lut,99)trim(text)
!   write(*,*)'3C elements mm: ',trim(text)
!------------------ system components including electrons for charged phases
   ip=1
   text=' '
   lcase=ichar('a')-ichar('A')
   do i1=1,noofel
! second letter lower case
      elsym=ellista(elements(i1))%symbol
      if(elsym(2:2).ne.' ') then
         elsym(2:2)=char(ichar(elsym(2:2))+lcase)
      endif
      text(ip:)=elsym
      ip=ip+25
      if(ip.gt.51) then
         write(lut,100)trim(text)
         ip=1
         text=' '
      endif
   enddo
! electrons
   nnn=1
   do i1=1,noofph
      if(estoi(i1).lt.0) then
         if(phchargedx(nnn).eq.i1) then
            text(ip:)='e('//trim(phcharged(nnn))//')'
            nnn=nnn+1
         else
            lokph=phasetuple(i1)%lokph
            text(ip:)='e('//trim(phlista(lokph)%name)//')'
         endif
         ip=ip+25
         if(ip.gt.51) then
            write(lut,100)trim(text)
            ip=1
            text=' '
         endif
      endif
   enddo
   if(ip.gt.1) then
      write(lut,100)trim(text)
   endif
! allocate an array for constituent stoichiometry
!   if(noofel+nelectrons.gt.50) &
   if(noelel.gt.50) &
        write(*,*)'Allocating large constituent array: ',noelel
   allocate(constcomp(noelel))
!----------------------------- system component mass, electrons 0.00054858???
   ip=1
   text=' '
   do i1=1,noofel
      write(text(ip:),130)ellista(elements(i1))%mass
130   format(F25.8)
      ip=ip+25
      if(ip.gt.51) then
         write(lut,100)trim(text)
         ip=1
         text=' '
      endif
   enddo
! electrons
   do i1=1,nelectrons
      write(text(ip:),130)5.4858D-4
      ip=ip+25
      if(ip.gt.51) then
         write(lut,100)trim(text)
         ip=1
         text=' '
      endif
   enddo
   if(ip.gt.1) then
      write(lut,100)trim(text)
   endif
!---------------------------------T powers, always the same line 
!   if(npows.eq.9) then
! 10 here are the allowed powers: 0 1 100 2 3 -1 ; 7 -9 -2  any any
!                                 1 2   3 4 5  6   7  8  9  10  11
! Those after the ; are special. 100 means T*ln(T)
   if(npows.le.15) then
! the first 7 digits should be 9 1..6
!      write(lut,140)trim(powers(36:))
!      write(lut,140)trim(powers(36:))
!140   format('   9   1   2   3   4   5   6',a)
! it does not seem to matter what is on these lines ...
      write(lut,140)
      write(lut,140)
!140   format('   6   1   2   3   4   5   6')
! According to Ted
140   format('6    1  2  3  4  5  6  ')
   else
      write(*,*)'3C too many different T powers: ',npows
      stop
   endif
!-------------------------------------- end of header section
! SOLGASMIX phase names must start with 4 unique letters, when TDB files
! has phases with same first 4 characters add a prefix
   phunique='P000'
! data for mixtures
! First the endmembers
   mphstoi=1
   phases1: do i1=1,noofph
      lokph=phasetuple(i1)%lokph
      skipfc=.false.
      if(ceq%phase_varres(phlista(lokph)%linktocs(1))%phstate.eq.PHSUS) then
! skip phases with suspended default composition set
!         write(*,*)'3C skipping phase loop 3: ',phlista(lokph)%name
         cycle phases1
      endif
! havemag nonzero if there are magnetic parameters
! magloop set to TRUE to list magnetic excess parameters
      havemag=0
      magloop=.FALSE.
      if(phlista(lokph)%nooffs.gt.1) then
! skip first ordered fraction set
         skipfc=.true.
      endif
      if(i1.eq.phstoi(mphstoi)) then
!         write(*,*)'3C skipping stoichiometric ',trim(phlista(lokph)%name)
         mphstoi=mphstoi+1
         cycle phases1
!      else
!         write(*,*)'3C parameters for mixture ',trim(phlista(lokph)%name)
      endif
      lokcs=phlista(lokph)%linktocs(1)
      varres=>ceq%phase_varres(lokcs)
! if disordered fraction set, set varres to point to disordered phase_varres
      if(skipfc) then
         varres=>ceq%phase_varres(lokcs)
         fedup=>varres
!         write(*,*)'3C disordered part: ',varres%disfra%ndd
         varres=>ceq%phase_varres(varres%disfra%varreslink)
      endif
      nsubl=1
      ionliq=.false.
! phase model ane expected endmembers
! we calculate the number of endmembers, end1mem is needed for DAT file
! end2mem is actual.  Error is not the same
      end1mem=0
      if(btest(phlista(lokph)%status1,PHIONLIQ)) then
         model='SUBI'
         nsubl=2
         ionliq=.true.
! there can just be one ionic liquid ... ??
         allocate(constcompiliq(noelel))
      elseif(btest(phlista(lokph)%status1,PHID)) then
         model='IDMX'         
      else
! there are phases with other bits which will not work but they are rarely set
! now for sublattices ...
         nsubl=phlista(lokph)%noofsubl
         offset=nsubl
         if(btest(phlista(lokph)%status1,PHFORD)) then
! NOTE varres is the disordered fraction set
            nsubl=size(varres%sites)
!            write(*,141)trim(phlista(lokph)%name),nsubl
141         format('Phase ',a,' has FCC permutated parameters, ignore ordered',&
                 i3)
         elseif(btest(phlista(lokph)%status1,PHBORD)) then
            nsubl=size(varres%sites)
!            write(*,142)trim(phlista(lokph)%name),nsubl
142         format('Phase ',a,' has BCC permutated parameters, ignore ordered',&
                 i3)
         elseif(btest(phlista(lokph)%status1,PHMFS)) then
            nsubl=size(varres%sites)
!            write(*,143)trim(phlista(lokph)%name),nsubl
143         format('Phase ',a,' has disorded fraction sets, ignore ordered',i3)
         endif
         if(nsubl.gt.1) then
            model='SUBL'
         else
            model='RKMP'
         endif
! fill values in xnooffr
         if(skipfc) then
            xnooffr=0
            do i2=1,nsubl
               xnooffr(i2)=xnooffr(i2-1)+fedup%disfra%nooffr(i2)
            enddo
            xnooffr(0)=1
         else
            xnooffr=0
            do i2=1,nsubl
               xnooffr(i2)=xnooffr(i2-1)+phlista(lokph)%nooffr(i2)
            enddo
            xnooffr(0)=1
         endif
!         write(*,*)'3C xnooffr: ',(xnooffr(i2),i2=0,nsubl)
! magnetism?
         addrec=>phlista(lokph)%additions
         lastadd: do while(associated(addrec))
! no need to increment CHTD except for magnetism
!            write(*,*)'3C additions?: ',phlista(lokph)%name,addrec%type
            if(addrec%type.eq.1) then
               aff=addrec%aff
               havemag=3
               model(5:5)='M'
               write(*,*)'3C magnetic phase 2: ',phlista(lokph)%name
            elseif(addrec%type.ne.7) then
! ignore addrec%type=7 which is volume model               
               write(*,*)'3C WARNING addition type: ',addrec%type,' ignored'
            endif
            addrec=>addrec%nextadd
         enddo lastadd
      endif
! prepare a dummy prefix
      phdummy=phlista(lokph)%name(1:4)
      jp=0
      name2: do i3=1,noofph
         if(i3.ne.lokph .and. phdummy.eq.phlista(i3)%name(1:4)) then
!            write(*,*)'Duplicate name',i3,lokph,phdummy,' ? ',&
!                 phlista(i3)%name(1:4)
            jp=1; exit name2
         endif
      enddo name2
      if(jp.gt.0) then
         call incunique(phunique)
         phdummy=phunique
!         write(*,*)'3C prefixing TDB phase name ',&
!              phdummy//'_'//trim(phlista(lokph)%name),i1,lokph
      else
         phdummy=' '
      endif
! According to Ted
      if(phdummy(1:1).eq.' ') then
         write(*,180)trim(phlista(lokph)%name),trim(model),&
              nsubl,ncon(i1),disfactor
180      format('3C mixture: ',a,' with model ',a,2i4,F12.4,a)
         write(lut,201)phlista(lokph)%name,nsubl,trim(model)
      else
         warnings=warnings+1
         write(*,180)phdummy//'_'//trim(phlista(lokph)%name),trim(model),&
              nsubl,ncon(i1),disfactor,' with name change'
         write(lut,201)phdummy//'_'//phlista(lokph)%name,nsubl,trim(model)
      endif
201   format(a,5x,'= MIXTURE PHASE =',i3/a)
      if(havemag.ne.0) then
         if(aff.eq.one) then
! Inden BCC magnetism
            write(lut,202)-aff,0.4
         else
! Inden FCC, HCP and other structures
            write(lut,202)-one/aff,0.28
         endif
202      format(F8.6,2x,F10.6)
      endif
!-------------------- we must repeat the endmember loop below for interactions
205   continue
      missend(1)=1
      do ip=2,nsubl
         missend(ip)=missend(ip-1)+phlista(lokph)%nooffr(ip-1)
      enddo
      endmember=>phlista(lokph)%ordered
      if(associated(phlista(lokph)%disordered)) then
! skip writing ordered part, nsubl set above!!
!         if(.not.skipfc) then
!            write(*,*)'3C We have disorderd fraction set but skipfc not set!'
!         else
!            write(*,*)'3C Skipfc set correctly',nsubl
!         endif
         write(*,*)'BEWARE skipping ordered part of :',&
              trim(phlista(lokph)%name),nsubl,offset
         endmember=>phlista(lokph)%disordered
      endif
!      write(*,*)'3C first the endmembers',nsubl
! endmember parameters, when they are done loop for excess parameters
      excessparam=.FALSE.
!===========================================================================
!================================== big loop for endmembers and interactions
! when all endmembers written then set excesspara=.true. and jump back here
207   continue
      if(ionliq) then
         nextcation=>endmember%nextem
         cation=endmember%fraclinks(1,1)
         if(.not.excessparam) then
! check if there is a missing endmember, skip wildcard parameters
            if(cation.ne.-99 .and. anion.ne.-99) then
               if(cation.ne.missend(1) .or. &
                    endmember%fraclinks(2,1).ne.missend(2)) then
                  write(*,*)'3C first endmember missing for liquid: ',&
                       missend(1),missend(2)
                  stop 'Check if inonic liquid has all endmember parameters'
               endif
            endif
         endif
! NOTE there can be missing endmembers!!
!         write(*,*)'3C firstcation: ',cation
         firstcation=cation
         iliqwild=.false.
         if(firstcation.eq.-99) then
            iliqwild=.true.
         else
            ccc=one
         endif
      endif
      lokcs=phlista(lokph)%linktocs(1)
      varres=>ceq%phase_varres(lokcs)
!--------------------------------------------------------------------
! here starts the loop for all parameters
! i1 is the index of this phase in the SOLGASMIX order
      allend: do while(associated(endmember))
! we have to generate two lines by extracting the endmember and constituents
! we may have to do this loop several times for the same phase to list
! the endmembers in correct order, at least for the ionic liquid phase
! For the ionic liquid all endmembers with the same cation must come together
         constcomp=zero
         if(ionliq) constcompiliq=zero
         constext=' '
         text=' '
         ip=1
         valency=zero
         wildcard=.false.
         if(.not.ionliq) then
!-----------------------------------------------
! for all other mixtures except ionic liquid ... note there are some tests
! of ionliq here as this loop originally was also for ionic liquids ...
            sloop1: do isubl=1,nsubl
! this is the loop for the constituents in sublattices
               if(skipfc) then
! We should skip the ordered sublattices
! for isubl=2 we should use the constituents in the last sublattce
                  isp=endmember%fraclinks(isubl,1)
!                  write(*,*)'3C constituent 1: ',isp,offset
               else
                  isp=endmember%fraclinks(isubl,1)
               endif
               intconst(isubl)=isp
               if(isp.eq.-99) then
! this means wildcard in this sublattice
                  wildcard=.true.
                  constext(ip:)='*:'
                  ip=ip+2
!                  if(ionliq .and. wildcard .and. isubl.eq.1) then
!                     iliqwild=.true.
!                     wildcard=.false.
!                  endif
! Hm we should add stoichiometric factors for all constituents in this subl
! For ionliq this means neutrals on sublattice 2
!>> QUESTION >> the DAT format repeats all neutrals for all cations
!>>>>>>>>>>>>>> with the stoichiometry of the cation !!!
!                  if(ionliq) valency(1)=one
                  cycle sloop1
               endif
               if(skipfc) then
! which index should be used to find the constituent in last sublattice
!                  write(*,*)'3C disordered species: ',isp
!                  nn=phlista(lokph)%constitlist(isp)
                  i3=firsteq%phase_varres(lokcs)%disfra%splink(isp)
!                  write(*,*)'3C disordered species: ',isp,nn,i3
                  isp=i3
!                  write(*,*)'3C species: ',splista(isp)%symbol
               else
                  isp=phlista(lokph)%constitlist(isp)
               endif
               if(btest(splista(isp)%status,SPVA)) then
                  valency(isubl)=zero
! according to the example I have the stoichiometry should be 1 for (cation:VA)
                  if(ionliq) valency(2)=-one
               else
                  valency(isubl)=splista(isp)%charge
                  if(abs(valency(isubl)).lt.1.0D-6) valency(isubl)=zero
               endif
! here we cannot have ionic liquid here!
!               if(ionliq .and. isubl.eq.2) then
!                  write(*,*)'3C we cannot have an ionic liquid here!'
!                  do i3=1,noofel
!                     constcomp(i3)=-constcomp(i3)*valency(2)
!                  enddo
!               elseif(estoi(i1).lt.0) then
               if(estoi(i1).lt.0) then
! charged sublattice phase.  Electronic stoichiometry should be positive!
!                  constcomp(-estoi(i1))=constcomp(-estoi(i1))+&
                  constcomp(-estoi(i1))=constcomp(-estoi(i1))-&
                       valency(isubl)*varres%sites(isubl)
!               write(*,901)'3C e-stoik:',isubl,-estoi(i1),&
!                    valency(isubl),varres%sites(isubl),constcomp(-estoi(i1))
901               format(a,2i3,3F10.2)
               endif
               call lower_case_species_name(constext,ip,isp)
               constext(ip:ip+1)=':'
               ip=ip+1
               do i2=1,splista(isp)%noofel
! this is a loop for the components of the endmember constituents
                  i3=ellista(splista(isp)%ellinks(i2))%alphaindex
                  if(i3.eq.0) then
! skip vacancies
                     continue
                  elseif(ionliq) then
                     write(*,*)'#C we should never be here if ionic liquid 2'
                     if(isubl.eq.1) then
                        constcomp(i3)=constcomp(i3)+&
                             splista(isp)%stoichiometry(i2)
                     else
                        constcomp(i3)=constcomp(i3)+&
                             splista(isp)%stoichiometry(i2)*valency(1)
                     endif
                  else
!
! here the stoichiometry of the endmember is added together
!
                     if(skipfc) then
                        constcomp(i3)=constcomp(i3)+&
                             splista(isp)%stoichiometry(i2)*&
                             varres%disfra%dsites(isubl)
                     else
                        constcomp(i3)=constcomp(i3)+&
                             splista(isp)%stoichiometry(i2)*varres%sites(isubl)
                     endif
                  endif
               enddo
            enddo sloop1
! for endmembers check that there is no missing endmember
            missend1: if(.not.excessparam) then
               www=0
               donotincrement=0
               miss7: do i2=1,nsubl
                  if(intconst(i2).eq.-99) then
! if we find a wildcard endmember do not increment missend !!!
                     www=0; goto 1814
                  endif
                  thisend(i2)=missend(i2)
                  if(intconst(i2).ne.missend(i2)) then
                     www=77
! this endmember is not the expected one.  There can be several missing
! but we should expect the one following.  That means we should reset
! constituents expected in higher sublattices ....                      
                     missend(i2)=intconst(i2)+1
                     if(donotincrement.eq.0) donotincrement=i2
                  endif
               enddo miss7
               if(www.ne.0) then
                  write(*,48)'3C *** Phase ',trim(phlista(lokph)%name),&
                       ' missing endmember: ',(thisend(i2),i2=1,nsubl)
!                  write(*,49)'3C found endmember: ',(intconst(i2),i2=1,nsubl)
                  warnings=warnings+1
48                format(a,a,a,9(1x,i3,':'))
49                format(a,19x,9(1x,i3,':'))
               endif
! increment constituents from the end for next test
! To handle also disordered fraction sets use varres pointer
! xnooffr(0) initially 1, xnooffr(j) is sum of constituents to and including j
               if(donotincrement.ne.nsubl) missend(nsubl)=missend(nsubl)+1
               do i2=nsubl,2,-1
                  if(missend(i2).gt.xnooffr(i2)) then
                     missend(i2)=xnooffr(i2-1)+1
                     missend(i2-1)=missend(i2-1)+1
                  endif
               enddo
1814           continue
!               write(*,49)'3C expecting:       ',(missend(i2),i2=1,nsubl)
            endif missend1
         else
!--------------------------------------------------------------------
! This is exclusivly for inonic liquids, loop second sublattice first ...
! this is the loop for the constituents in sublattices
! Hm we should add stoichiometric factors for all constituents in this subl
!            write(*,*)'3C we are here 1 ',excessparam,firstcation
            if(.not.iliqwild) then
               isp=phlista(lokph)%constitlist(cation)
               intconst(1)=cation
               valency(1)=splista(isp)%charge
               cationval=valency(1)
               do i2=1,splista(isp)%noofel
! this is a loop for the components of the endmember constituents
                  i3=ellista(splista(isp)%ellinks(i2))%alphaindex
                  if(i3.eq.0) then
! skip vacancies
                     continue
                  else
                     constcomp(i3)=constcomp(i3)+&
                          splista(isp)%stoichiometry(i2)
                  endif
               enddo
               call lower_case_species_name(constext,ip,isp)
               constext(ip:ip+1)=':'
               ip=ip+1
            else
               valency(1)=one
            endif
! what about neutrals?
            anion=endmember%fraclinks(2,1)
            intconst(2)=anion
            isp=phlista(lokph)%constitlist(anion)
            missend2: if(.not.excessparam) then
               if(cation.ne.missend(1) .or. anion.ne.missend(2)) then
                  write(*,47)'3C **** liquid missing endmember: ',&
                       missend(1),missend(2)
47                format(a,2i5,5x,2i5)
                  stop 'Missing endmember in ionic liquid'
                  warnings=warnings+1
! avoid having several errors due to a missing cation:anion pair
                  missend(1)=cation
               endif
! Hm, cation should not be incremented here ...
!               missend(1)=cation+1
               missend(2)=anion+1
               if(anion.eq.phlista(lokph)%tnooffr) then
                  missend(1)=missend(1)+1
                  missend(2)=phlista(lokph)%nooffr(1)+1
               endif
            endif missend2
            if(btest(splista(isp)%status,SPVA)) then
! according to the example I have the stoichiometry should be 1 for (cation:VA)
               valency(2)=-one
            else
               valency(2)=splista(isp)%charge
               if(abs(valency(2)).lt.1.0D-6) valency(2)=zero
            endif
! This is values in the stoichiometry line ....
            do i3=1,noofel
               constcomp(i3)=-constcomp(i3)*valency(2)
            enddo
            call lower_case_species_name(constext,ip,isp)
            constext(ip:ip+1)=':'
            ip=ip+1
            do i2=1,splista(isp)%noofel
! this is a loop for the components of the endmember constituents
               i3=ellista(splista(isp)%ellinks(i2))%alphaindex
               if(i3.eq.0) then
! skip vacancies
                  continue
               elseif(ionliq .and. iliqwild) then
! For neutrals in ionic liquid we must multiply with ccc (the cation valency)
                  constcomp(i3)=constcomp(i3)+&
                       splista(isp)%stoichiometry(i2)*valency(1)*ccc
               else
                  constcomp(i3)=constcomp(i3)+&
                       splista(isp)%stoichiometry(i2)*valency(1)
               endif
            enddo
!            write(*,917)'3C Ionliq endmember: ',constext(1:ip-2),iliqwild,ccc,&
!                 valency(2),(constcomp(i3),i3=1,noofel)
917         format(a,a,L3,2F10.2/10F7.3)
!------------------ end special ionic liquid
         endif
!         write(*,*)'3C we are here 2 '
         endorexcess: if(excessparam) then
! we can have several excess parameters for each endmember
            intparam=>endmember%intpointer
            ilevel=0
!            write(*,*)'3C we are here 3 '
            intree: do while(associated(intparam))
! we must save intparam%nextlink to be able to follow the parameter tree
               ilevel=ilevel+1
               saveint(ilevel)%intlink=>intparam%nextlink
               isp=intparam%fraclink(1)
               intconst(nsubl+ilevel)=isp
               isp=phlista(lokph)%constitlist(isp)
               property=>intparam%propointer
! Check if endmember contains wildcard
               if(wildcard .and. associated(property)) then
                  write(*,903)'3C Expanding wildcard interaction: ',&
                       trim(phlista(lokph)%name),trim(constext),&
                       (intconst(k1),k1=1,nsubl+ilevel)
903               format(a,a,',',a,2x,6i4)
! we should make a loop fof all constituents in sublattice with wildcard
! and write the same parameter for all.  There can be several wildcards!!
! like G(C1_MO2,Zr+2:*:*), where *=(O-2,Va) in both cases
! wildloop expanded constituent sets returned in iset, allocated inside
                  call expand_wildcards(intconst,nsubl+ilevel,&
                       wildloop,iset,lokph)
!                  wildloop=1
! replace current intconst with values in iset and loop below back to 310
                  do k1=1,nsubl+ilevel
                     intconsty(k1)=intconst(k1)
                     intconst(k1)=iset(k1,wildloop)
                  enddo
!                  write(*,324)'3C wildloop1: ',wildloop,&
!                       (intconst(k1),k1=1,nsubl+ilevel)
                  savedproperty=>property
               else
                  wildloop=0
               endif
! return here with new set of constituents in intconst if wildloop not zero
310            continue
               maxideg=-1
               extcpar=zero; exbmpar=zero
               intproploop: do while(associated(property))
! Check type of excess parameter and what kind to be listed ....
                  if(magloop) then
                     if(property%proptype.eq.2) then
! this is Curie/Neel temperature
                        do ideg=0,property%degree
                           f1=property%degreelink(ideg)
                           if(f1.gt.0) then
                              extcpar(ideg)=tpfc(f1)%cfun%coefs(1,1)
                           else
                              write(*,315)' 3C zero excess TC parameter: ',&
                  trim(tpfuns(property%degreelink(property%degree))%symbol)
!                                   trim(phlista(lokph)%name),ideg,ilevel
315                           format(a,a,5i5)
                              extcpar(ideg)=zero
                           endif
                        enddo
!                        write(*,*)'3C excess TC: ',f1,partc
                        paratyp=17
                        if(ideg.gt.maxideg) maxideg=ideg
                     elseif(property%proptype.eq.3) then
! This is BMAGN
                        do ideg=0,property%degree
                           f1=property%degreelink(ideg)
                           if(f1.gt.0) then
                              exbmpar(ideg)=tpfc(f1)%cfun%coefs(1,1)
                           else
                              write(*,315)' 3C zero excess BM parameter: ',&
                  trim(tpfuns(property%degreelink(property%degree))%symbol)
                              exbmpar(ideg)=zero
                           endif
                        enddo
                        paratyp=17
                        if(ideg.gt.maxideg) maxideg=ideg
                     endif
                     property=>property%nextpr
                     cycle intproploop
                  elseif(property%proptype.ne.1) then
! we should have a loop here also as G not always first parameter
                     continue
                  endif
! write the identification of the excess parameter ....
! The list of constituents (in intconst) arranged in ascending order
                  call intsort(intconst,nsubl+ilevel,intconstx)
! write interaction level (2=binary, 3=ternary ...)
! Then constituent indices in acending order (maybe rearrange intconst)
! finally the degree (number of Redlich-Kister parameters)
!                  write(*,907)'3C solgasorder: ',nsubl+ilevel,&
!                       (intconstx(k1),k1=1,nsubl+ilevel),property%degree+1
! write an excess parameter
                  write(lut,208)nsubl+ilevel,&
                       (intconstx(k1),k1=1,nsubl+ilevel),property%degree+1
907               format(a,10i5)
208               format(i5/10i5)
! write the expression of the excess parameter .... (Redlich-Kister ??)
                  alldegs: do ideg=0,property%degree
                     f1=property%degreelink(ideg)
! excess parameters has just the coefficients
!                     call list_tpascoef(lut,text,f1,npows,tpfc)
                     if(f1.lt.1) then
! This means one RK parameter is zero!! L(FCC,NB:C,Va,1) is zero !!1
!                        write(*,*)'3C No function?: ',f1,ideg,property%degree
                        write(*,*)'3C zero RK paramameter: ',&
                             tpfuns(property%degreelink(property%degree))%symbol
                        write(lut,307)0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
307                     format(6(1x,G15.8))
                        cycle alldegs
                     endif
                     if(tpfc(f1)%nranges.gt.1) then
                        write(*,*)'3C excess parameter with T-ranges!'
                        stop
                     endif
! This gave compiler error on MacOS 10.13 ??? GNU Fortran 5.2 ...
!                     write(lut,311)(tpfc(f1)%cfun%coefs(jj,1),jj=1,6)
! write another excess parameter.  What about magnetism and paratype???
                     write(lut,311)tpfc(f1)%cfun%coefs(1,1),&
                          tpfc(f1)%cfun%coefs(2,1),tpfc(f1)%cfun%coefs(3,1),&
                          tpfc(f1)%cfun%coefs(4,1),tpfc(f1)%cfun%coefs(5,1),&
                          tpfc(f1)%cfun%coefs(6,1)
311                  format(6(1x,G15.8))
                  enddo alldegs
                  property=>property%nextpr
               enddo intproploop
!               if(wildloop.gt.0) write(*,*)'3C wildloop2: ',wildloop
! magnetic excess parameter not written above but here
!               write(*,*)'3C exit intproploop',magloop,paratyp,&
!                    associated(property)
               if(magloop .and. paratyp.eq.17) then
                  paratyp=4
                  call intsort(intconst,nsubl+ilevel,intconstx)
!                  write(*,907)'3C solgasorder: ',nsubl+ilevel,&
!                       (intconstx(k1),k1=1,nsubl+ilevel),1
                  write(lut,208)nsubl+ilevel,&
                       (intconstx(k1),k1=1,nsubl+ilevel),maxideg
                  write(lut,323)(extcpar(ideg),exbmpar(ideg),ideg=0,maxideg-1)
323               format(2F12.3)
! end of output of magnetic excess parameter
               endif
! If this is a wildcard parameter maybe it should be written several times
               if(wildcard) then
                  if(wildloop.gt.1) then
                     wildloop=wildloop-1
                     do k1=1,nsubl+ilevel
                        intconst(k1)=iset(k1,wildloop)
                     enddo
!                     write(*,324)'3C next expanded: ',wildloop,&
!                          (intconst(k1),k1=1,nsubl+ilevel)
324                  format(a,i3,2x,10i4)
                     property=>savedproperty
                     goto 310
                  else
! deallocate iset and restore intconst as we may have higher interactions ...
                     deallocate(iset)
                     do k1=1,nsubl+ilevel
                        intconst(k1)=intconsty(k1)
                     enddo
                  endif
               endif
! Take link to higher
               intparam=>intparam%highlink
               do while(ilevel.gt.0 .and. .not.associated(intparam))
! go down the saved links
                  intparam=>saveint(ilevel)%intlink
                  ilevel=ilevel-1
               enddo
            enddo intree
         else
! here we are writing endmembers, we have generated the endmember symbol,
! for the parameters follow the property link
!            write(*,*)'3C We are here 4'
            property=>endmember%propointer
            if(wildcard .and. associated(property)) then
               write(*,*)'3C ERROR! Endmember parameter with wildcard: ',&
                    trim(phlista(lokph)%name),',',trim(constext)
            endif
            paratyp=4
            partc=zero; parbm=zero
! return here if we find a magnetic property first
333         continue
            propem: if(associated(property)) then
! some endmembers may not have a property record!!
               if(property%proptype.ne.1) then
! for magnetism we can have proptype 1 and 2 (TC and BMAGN)
! They can be before the G parameter in the TDB file.
!                  write(*,*)'3C magnetic 1: ',trim(phlista(lokph)%name),&
!                       havemag,property%proptype,paratyp
                  if(havemag.ne.0) then
                     if(property%proptype.eq.2) then
! this is Curie/Neel temperature
                        f1=property%degreelink(0)
                        partc=tpfc(f1)%cfun%coefs(1,1)
!                        write(*,*)'3C endmember TC: ',f1,partc
                        paratyp=16
                     elseif(property%proptype.eq.3) then
! This is BMAGN
                        f1=property%degreelink(0)
                        parbm=tpfc(f1)%cfun%coefs(1,1)
!                        write(*,*)'3C endmember BMAGN: ',f1,parbm
                        paratyp=16
                     else
                        write(*,*)'3C skipping magnetic endmember property: ',&
                             property%proptype
                        exit propem
                     endif
                  else
                     write(*,*)'3C unknown endmember property: ',&
                          property%proptype
                     exit propem
                  endif
                  if(associated(property%nextpr)) then
                     property=>property%nextpr
                     goto 333
                  endif
!               else
!                  paratyp=4
               endif
! this line with the stoichiometry of the endmember should be written
! together with the type of coefficients and number of ranges
! it may require several lines
               write(text,210)constcomp
! THIS IS THE STOICHIMETRY OF THE ENDMEMBER, with 6 decimal digits
! If this format is changed the output routine list_tpascoef must be changed!
!210            format(60(1x,F11.6))
! ERNESTO GEIGER complained it did not work ... this is stoichiometry format
210            format(60(1x,F7.2))
! Check if any value in contcomp is greated than 1000, could give overflow
! Check also if two decimals not enough
               do i3=1,noofel
                  if(constcomp(i3).gt.maxcc) then
                     warnings=warnings+1
                     write(*,206)trim(phlista(lokph)%name),i3,constcomp(i3)
206                  format('3C *** Warning stoichiometry factor >100: ',&
                          a,i4,F10.2)
                  endif
                  decimals=int(1.0D5*constcomp(i3))
                  xxx=1.0D-5*dble(decimals)
                  if(abs(xxx-constcomp(i3)).gt.1.0D-6) then
                     warnings=warnings+1
                     write(*,203)trim(constext),i3,constcomp(i3)
203                  format('3C *** Warning stoichiometry with >5 decimals: ',&
                          a,i4,2F10.6)
                  endif
               enddo
! property record has property=1 it is G; take care of magnetic properties
!               write(*,*)'3C havemag: ',trim(phlista(lokph)%name),havemag
               magprop: if(havemag.gt.0) then
                  nextprop=>property%nextpr
334               continue
!                  write(*,*)'3C magnetic 2: ',trim(phlista(lokph)%name),&
!                       property%proptype,associated(nextprop)
                  if(associated(nextprop)) then
                     if(nextprop%proptype.eq.2) then
! this is Curie/Neel temperature
                        f1=nextprop%degreelink(0)
                        partc=tpfc(f1)%cfun%coefs(1,1)
!                        write(*,*)'3C endmember TC2: ',f1,partc
                        paratyp=16
                     elseif(nextprop%proptype.eq.3) then
! This is BMAGN
                        f1=nextprop%degreelink(0)
                        parbm=tpfc(f1)%cfun%coefs(1,1)
!                        write(*,*)'3C endmember BMAGN2: ',f1,parbm
                        paratyp=16
                     else
                        write(*,*)'3C ignoring endmember property: ',&
                             nextprop%proptype
                     endif
                  else
                     exit magprop
                  endif
                  nextprop=>nextprop%nextpr
                  goto 334
               endif magprop
! property record has still property=1 it is G
               f1=property%degreelink(0)
               if(f1.gt.0) then
                  factor=one
                  if(ionliq .and. iliqwild) then
                     write(lut,211)constext(1:ip-2),ccc
! According to Ted
211                  format(a,40x,' * ',F12.2)
! We must multiply tpfc(f1) with ccc, store in tpfc(jp) coefficient function!
                     jp=ntpf+1
                     call tpmult(f1,jp,ccc,tpfc)
                     call list_tpascoef(lut,text,paratyp,jp,npows,factor,tpfc)
                     if(paratyp.eq.16) write(lut,222)partc,parbm
222                  format(2G15.8)
                  else
! according to Ted: endmember symbol 
!                     write(*,99)constext(1:ip-2)
                     write(lut,99)constext(1:ip-2)
                     call list_tpascoef(lut,text,paratyp,f1,npows,factor,tpfc)
                     if(paratyp.eq.16) write(lut,222)partc,parbm
                  endif
               else
                  write(*,*)'3 C missing function for endmember property',&
                       constext(3:ip-2)
               endif
               endy(i1)=endy(i1)+1
            endif propem
         endif endorexcess
! take next endmember
!         write(*,*)'3C We are here 5'
         if(.not.ionliq) then
            endmember=>endmember%nextem
         else
! find next endmember with the same cation, liquids without cations? !!
! if none set endmember=>nextcation
! if nextcation has same cation as firstcation we have finished!
240         continue
            iliqwild=.false.
!            write(*,241)'ionliq done:   ',firstcation,cation,&
!                 endmember%fraclinks(1,1),endmember%fraclinks(2,1)
            endmember=>endmember%nextem
            if(associated(endmember)) then
!               write(*,241)'ionliq ass:    ',firstcation,cation,&
!                    endmember%fraclinks(1,1),endmember%fraclinks(2,1)
241            format(a,2i3,2x,2i3)
               if(endmember%fraclinks(1,1).eq.-99) then
                  iliqwild=.true.
! ccc  is the valency of the cation used to multiply the neutral parameter
                  ccc=cationval
               elseif(endmember%fraclinks(1,1).ne.cation) then
                  goto 240
               endif
            else
!               write(*,*)'3C we are here 6: ',associated(nextcation)
               endmember=>nextcation
               if(associated(endmember)) then
                  nextcation=>nextcation%nextem
                  cation=endmember%fraclinks(1,1)
!               write(*,241)'ionliq notaass: ',firstcation,cation,&
!                    endmember%fraclinks(1,1),endmember%fraclinks(2,1)
! we have looped through all cations
                  if(cation.eq.firstcation) exit allend
! there were just one cation but some neutrals (already listed)
                  if(endmember%fraclinks(1,1).eq.-99) exit allend
!               else
! no more cations, finished!
!                  write(*,*)'3C no nextcation!'
               endif
            endif
         endif
      enddo allend
! ------------------- end of endmembers, constituents and excess parameters ??
      if(model(1:4).eq.'IDMX') cycle phases1
      if(excessparam) goto 297
! After endmembers for sublattice phases write number of sublattices and sites
      if(model(1:4).eq.'SUBL') then
         write(lut,250)nsubl
         if(skipfc) then
            write(lut,260)(varres%disfra%dsites(isubl),isubl=1,nsubl)
         else
            write(lut,260)(ceq%phase_varres(lokcs)%sites(isubl),isubl=1,nsubl)
250      format(1x,i4)
260      format(1x,8F9.5)
         endif
      endif
!      write(*,*)'3C here 8: ',phlista(lokph)%name,model
      if(model(1:4).eq.'SUBL' .or. model(1:4).eq.'SUBI') then
! number of constituents in each sublattice
         if(skipfc) then
            write(lut,270)(varres%disfra%nooffr(isubl),isubl=1,nsubl)
         else
            write(lut,270)(phlista(lokph)%nooffr(isubl),isubl=1,nsubl)
270         format(9i5)
         endif
      endif
! For all phases with sublattices we should write the constituents of each
! problem here for UC2_C11A, constituent in first sublattice ignored
      if(nsubl.eq.1) goto 280
      i3=0
!      do isubl=1,phlista(lokph)%noofsubl
      do isubl=1,nsubl
         constext=' '
         ip=1
         if(skipfc) then
            nn=varres%disfra%nooffr(isubl)
         else
            nn=phlista(lokph)%nooffr(isubl)
         endif
         do i2=1,nn
            i3=i3+1
            if(skipfc) then
               isp=firsteq%phase_varres(lokcs)%disfra%splink(i3)
            else
               isp=phlista(lokph)%constitlist(i3)
            endif
            jp=ip
            call lower_case_species_name(constext,ip,isp)
            ip=jp+25
            if(ip.ge.75) then
               write(lut,100)trim(constext)
               constext=' '
               ip=1
            endif
         enddo
!         write(*,271)'3C constext: ',trim(constext),isubl,i2,i3,ip
!271      format(a,a,4i4)
         ip=len_trim(constext)
! for a single component names ip=1 here ...
         if(ip.gt.1 .or. constext(1:1).ne.' ') then
            write(lut,100)trim(constext)
         endif
      enddo
280   continue
      if(model(1:4).eq.'SUBI') then
! There should be a line with just a "2" ???
         write(lut,272)
272      format('   2')
! for ionic liquid list abs(valencies) of constituents, one line per sublattice
         ip=1
         isp=1
         constext=' '
         do i2=1,phlista(lokph)%nooffr(1)
            ccc=splista(phlista(lokph)%constitlist(isp))%charge
            write(constext(ip:),274)ccc
274         format(F10.5)
            ip=len_trim(constext)
            if(ip.gt.69) then
               write(lut,99)trim(constext)
               ip=1
               constext=' '
            endif
            isp=isp+1
         enddo
         if(ip.gt.1) then
            write(lut,99)trim(constext)
         endif
         ip=1
         constext=' '
         do i2=1,phlista(lokph)%nooffr(2)
! For anions the charge as a positive value, for Va unity, for neutrals zero
            if(btest(splista(phlista(lokph)%constitlist(isp))%status,SPVA)) then
               ccc=one
            else
               ccc=abs(splista(phlista(lokph)%constitlist(isp))%charge)
            endif
            write(constext(ip:),274)ccc
            ip=len_trim(constext)
            if(ip.gt.69) then
               write(lut,99)trim(constext)
               ip=1
               constext=' '
            endif
            isp=isp+1
         enddo
         if(ip.gt.1) then
            write(lut,99)trim(constext)
         endif
      endif
!      if(phlista(lokph)%noofsubl.gt.1) then
!      lastix=0
      if(nsubl.gt.1) then
! A very strange output of integers representing endmembers
         jp=1
         mult=1
!         do isubl=phlista(lokph)%noofsubl,1,-1
         do isubl=nsubl,1,-1
            mult(isubl)=jp
            if(skipfc) then
               jp=jp*varres%disfra%nooffr(isubl)
            else
               jp=jp*phlista(lokph)%nooffr(isubl)
            endif
         enddo
!         write(*,278)'3C mult2: ',jp,(mult(ip),ip=1,phlista(lokph)%noofsubl)
278      format(a,10i4)
         do isubl=1,nsubl
            lastix=0
            text=' '
            ip=3
            k1=0
            i2=0
290         continue
               k1=k1+1
               i3=0
292            continue
                  lastix=lastix+1
                  exix(lastix)=k1
! the use of text here will be made redundant
                  call wriint(text,ip,k1)
                  ip=ip+3
                  i2=i2+1
                  i3=i3+1
                  if(i3.lt.mult(isubl)) goto 292
               if(skipfc) then
                  if(k1.gt.varres%disfra%nooffr(isubl)) k1=0
                  if(k1.eq.varres%disfra%nooffr(isubl) .and. isubl.gt.1) k1=0
               else
                  if(k1.gt.phlista(lokph)%nooffr(isubl)) k1=0
                  if(k1.eq.phlista(lokph)%nooffr(isubl) .and. isubl.gt.1) k1=0
               endif
            if(i2.lt.jp) goto 290
! According to Markus Piro one should have 19 values per line, 18*4+3=75
! New code using i4 format
            lineb=1
            firstix=1
            do while(lastix.gt.lineb)
               lineb=min(firstix+18,lastix)
!               write(*,*)'3E firstix: ',firstix,lineb
               write(lut,'(19i4)')(exix(isp),isp=firstix,lineb)
               firstix=lineb+1
            enddo
! output below is wrong and removed redunant
!            isp=1
!            do while(len_trim(text(isp:))-76.gt.0)
! Corrected 2020-11-12 with the help from Max Poschmann and Markus Piro
!            do while(len_trim(text(isp:))-76.gt.0)
!               write(lut,99)trim(text(isp:isp+74))
!               isp=isp+75
!               lineb=75
!               do while(text(isp+lineb:isp+lineb).ne.' ')
! increment lineb until we find a space
!                  lineb=lineb+1
!               enddo
!               write(*,*)'3E linebreak: "',text(isp+lineb-1:isp+lineb-1),&
!                    '" and "',text(isp+lineb:isp+lineb),'"',lineb
!               write(lut,99)trim(text(isp:isp+lineb-1))
!               isp=isp+lineb
!            enddo
!            if(len_trim(text(isp:)).gt.0) write(lut,99)trim(text(isp:))
         enddo
      endif
!...................... repeat loop for excess parameters
297   continue
      if(.not.excessparam) then
! repeat the endmember loop again for interaction parameters (and magnetism??)
!         write(*,*)'3C Now the excess parameters',nsubl
         excessparam=.true.
! if magnetic we have FIRST loop all excess parameters for magnetic parameters
         if(havemag.ne.0) magloop=.TRUE.
! and then again for the G parameters .... SUCK
         endmember=>phlista(lokph)%ordered
         if(associated(phlista(lokph)%disordered)) then
            endmember=>phlista(lokph)%disordered
         endif
!         if(magloop) write(*,*)'3C First magnetic excess parameters'
         goto 207
      elseif(magloop) then
! First finish the magetic excess parameter parameters with a zero
         write(lut,555)
555      format(' 0',30x,' = end of magnetic excess parameters')
! here we write the Gibbs energy excess parameters
         magloop=.FALSE.
         endmember=>phlista(lokph)%ordered
         if(associated(phlista(lokph)%disordered)) then
            endmember=>phlista(lokph)%disordered
         endif
!         write(*,*)'3C Gibbs energy excess parameters after magnetic'
         goto 207
      endif
! terminate the excess parameters for this phase with a line starting with 0
      write(lut,300)
300   format(' 0')
   enddo phases1
!-------------------------------------------------------
! now data for stoichiometric phases
   mphstoi=1
!   write(*,*)
!   write(*,*)'3C loop for compounds ',nphstoi
!
   phases2: do i1=1,noofph
      lokph=phasetuple(i1)%lokph
      if(ceq%phase_varres(phlista(lokph)%linktocs(1))%phstate.eq.PHSUS) then
! skip phases with suspended default composition set
!         write(*,*)'3C skipping phase loop 4: ',phlista(lokph)%name
         cycle phases2
      endif
      if(i1.ne.phstoi(mphstoi)) then
!         write(*,*)'3C skipping mixture ',trim(phlista(lokph)%name),&
!              i1,mphstoi,phstoi(mphstoi)
         cycle phases2
      endif
      mphstoi=mphstoi+1
      skipfc=.FALSE.
      factor=one
      if(phlista(lokph)%nooffs.gt.1) then
! skip first composition set
         skipfc=.true.
      endif
! magnetism?
      havemag=0
      addrec=>phlista(lokph)%additions
      lastadd2: do while(associated(addrec))
! no need to increment CHTD except for magnetism
!         write(*,*)'3C additions?: ',phlista(lokph)%name,addrec%type
         if(addrec%type.eq.1) then
            havemag=3
            write(*,*)'3C magnetic phase 1: ',phlista(lokph)%name
            aff=addrec%aff
         elseif(addrec%type.ne.7) then
! type 7 is volume
            write(*,*)'3C WARNING addition type: ',addrec%type,' ignored'
         endif
         addrec=>addrec%nextadd
      enddo lastadd2
      lokcs=phlista(lokph)%linktocs(1)
      varres=>ceq%phase_varres(lokcs)
      nsubl=1
      ionliq=.false.
      nsubl=phlista(lokph)%noofsubl
      if(skipfc) then
         factor=varres%disfra%fsites
         varres=>ceq%phase_varres(varres%disfra%varreslink)
         if(btest(phlista(lokph)%status1,PHMFS)) then
            nsubl=size(varres%sites)
         endif
         endmember=>phlista(lokph)%disordered
      else
! there is just one endmember!!
         endmember=>phlista(lokph)%ordered
      endif
! prepare a dummy prefix for compounds ... NOT NECESSARY
!      phdummy=phlista(lokph)%name(1:4)
!      jp=0
!      do i3=1,noofph
!         if(i3.ne.lokph .and. phdummy.eq.phlista(i3)%name(1:4)) jp=1
!      enddo
!      if(jp.gt.0) then
!         warnings=warnings+1
!         call incunique(phunique)
!         phdummy=phunique
!         write(*,*)'3C prefixing TDB phase name ',&
!              phdummy//'_'//trim(phlista(lokph)%name),i1
!      else
!         phdummy=' '
!      endif
      phdummy=' '
      if(phdummy(1:1).eq.' ') then
         write(*,477)trim(phlista(lokph)%name),nsubl,factor
477      format('3C Compound: ',a,i3,F12.3,a)
! write on file
         write(lut,500)phlista(lokph)%name,factor
500      format(1x,a,5x,'= COMPOUND PHASE = ',F12.4)
      else
         write(*,477)phdummy//'_'//trim(phlista(lokph)%name),nsubl,factor,&
              ' with name change'
         write(lut,500)phdummy//'_'//phlista(lokph)%name,factor
      endif
      if(havemag.ne.0) then
         if(aff.eq.one) then
            write(lut,202)-aff,0.4
         else
            write(lut,202)-one/aff,0.28
         endif
      endif
      constext=' '
      ip=1
      constcomp=zero
      sloop2: do isubl=1,nsubl
! this is the loop for the constituents in sublattices
         if(.not.associated(endmember)) then
            write(*,*)'3C no parameter!! ',phlista(lokph)%name
            cycle sloop2
         endif
         isp=endmember%fraclinks(isubl,1)
         if(isp.eq.-99) then
! this means wildcard in this sublattice
            write(*,*)'3C *** ERROR! Wildcard in a stoichiometric compound!!!'
            constext(ip:)='*:'
            ip=ip+2
            cycle sloop2
         endif
! Hm we should add stoichiometric factors for all constituents in this subl
         isp=phlista(lokph)%constitlist(isp)
         if(btest(splista(isp)%status,SPVA)) then
            write(*,*)'3C Warning: vacancy in stoichiometric compound!!'
         endif
         write(constext(ip:),99)trim(splista(isp)%symbol)//':'
         ip=len_trim(constext)+1
         do i2=1,splista(isp)%noofel
! this is a loop for the components of the endmember constituents
            i3=ellista(splista(isp)%ellinks(i2))%alphaindex
            if(i3.eq.0) then
! skip vacancies
               continue
            else
               constcomp(i3)=constcomp(i3)+splista(isp)%stoichiometry(i2)*&
                    varres%sites(isubl)
            endif
         enddo
      enddo sloop2
! we may come here if there are no endmembers!
      if(.not.associated(endmember)) then
         write(*,*)'3C skipping this phase'
         cycle phases2
      endif
! for the parameters follow the property link
      property=>endmember%propointer
      if(associated(property)) then
! For a compound do not write any constituent array
!         write(lut,100)constext(1:ip-2)
! this line should be written together with the type of coefficients and ranges
! it may require several lines
         write(text,210)constcomp
! Check if any value in contcomp is greated than 1000, could give overflow
         do i3=1,noofel
            if(constcomp(i3).gt.maxcc) then
               warnings=warnings+1
               write(*,206)trim(phlista(lokph)%name)
            endif
         enddo
         paratyp=4
         partc=zero; parbm=zero
         if(havemag.ne.0) paratyp=16
! what about several properties?? YES
575      continue
         if(property%proptype.eq.1) then
            f1=property%degreelink(0)
            if(f1.gt.0) then
               call list_tpascoef(lut,text,paratyp,f1,npows,factor,tpfc)
            else
               write(*,*)'missing endmember parameter'
            endif
         elseif(property%proptype.eq.2) then
            f1=property%degreelink(0)
            partc=tpfc(f1)%cfun%coefs(1,1)
         elseif(property%proptype.eq.3) then
            f1=property%degreelink(0)
            parbm=tpfc(f1)%cfun%coefs(1,1)
         else
            write(*,*)'3C ignoring compound property ',property%proptype
         endif
         property=>property%nextpr
         if(associated(property)) goto 575
      endif
      if(paratyp.eq.16) write(lut,222)partc,parbm
   enddo phases2
! At the end some dummy line for the pure elements??
   write(lut,602)
602 format('###################################################')
   goto 900
!----------------------- ???
! At the end some dummy line for the pure elements??
   do i1=1,noofel
      write(lut,605)ellista(elements(i1))%symbol
605   format(1x,a2,22x,'#')
      constcomp=zero
      constcomp(i1)=one
      write(lut,610)constcomp
610   format('   4  1',12F7.1)
      write(lut,620)
620   format('  6001.0000     0.00000000     0.00000000     0.00000000',&
           '     0.00000000 '/' 0.00000000     0.00000000 '/&
           ' 1 0.00000000       0.00')
   enddo
!
900 continue
   do i1=1,noofph
      if(endx(i1).ne.endy(i1)) then
         lokph=phases(i1)
         write(*,911)trim(phlista(lokph)%name),endx(i1),endy(i1)
911      format('3C Endmembers missing for ',a,&
              ', should have ',i3,' endmembers, has ',i3)
      endif
   enddo
   write(*,700)noofph,nphmix,nphstoi
700 format('3C written data for ',i4,' phases: ',i3,' mixtures and ',&
         i4,' compounds')
   if(warnings.gt.0) write(*,701)warnings
701 format(' *** Attention: there were ',i3,' warnings!')
! 
1000 continue
! Finished SOLGASMIX outpur
   if(allocated(tpfc)) deallocate(tpfc)
   if(gx%bmperr.ne.0) then
      write(*,1009)trim(filename),gx%bmperr
1009  format(/' *** Output terminated on ',a,' due to error ',i5/)
   elseif(date(1:4).ne.'    ') then
      write(*,1010)trim(filename)
1010  format('3C Output finished on ',a/)
   else
      write(*,1020)trim(filename)
1020  format('3C no output on ',a/)
   endif
   close(lut)
   return
 end subroutine save_datformat

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine incunique
!\begin{verbatim}
 subroutine incunique(text)
   character text*(*)
!\end{verbatim} %+
   integer j1,j2,j3
   j1=len(text)
!   write(*,*)'3C phunique 1: ',text
   loop: do while(j1.ge.1)
      j2=ichar(text(j1:j1))-ichar('0')
! this position is not a number, exit
      if(j2.lt.0) exit loop
      if(j2.lt.9) then
! increment the number and exit
         text(j1:j1)=char(j2+1+ichar('0'))
         exit loop
      elseif(j2.eq.9) then
         text(j1:j1)='0'
         j1=j1-1
      else
! this position is not a number, exit
         exit loop
      endif
   enddo loop
!   write(*,*)'3C phunique 2: ',text
   return
 end subroutine incunique

 !/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine expand_wildcards
!\begin{verbatim}
 subroutine expand_wildcards(intconst,nconst,wildloop,iset,lokph)
! Expand a wildcard constituent with all constituents it replaces
! There can be several wildcards
! intconst is the original set of constuents including the wildcards (-99)
! nconst is the number of constituents
! wildloop is set to the number of times the interaction is repeated
! iset is a matrix with the expanded constituents
! phrecord is the phase record where one can find the phase structure
   implicit none
   integer intconst(*)
   integer, allocatable, dimension(:,:) :: iset
   integer nconst,wildloop,lokph
!\end{verbatim} %+
   integer la,lb,lc,lz,ja,jb,jc,jz,ka,kb,nexp
   integer, allocatable, dimension(:) :: multi
!   write(*,10)'3C in expand_wildcard: ',nconst,(intconst(la),la=1,nconst)
10 format(a,i3,2x,10i4)
   nexp=1
   allocate(multi(phlista(lokph)%noofsubl))
   multi=1
   do la=1,phlista(lokph)%noofsubl
      if(intconst(la).eq.-99) then
         multi(la)=nexp
         nexp=nexp*phlista(lokph)%nooffr(la)
      endif
   enddo
!   write(*,*)'3C expand: ',nconst,nexp
   allocate(iset(nconst,nexp))
! initiate iset to original constituents (with wildcards)
   do la=1,nexp
      do ja=1,nconst
         iset(ja,la)=intconst(ja)
      enddo
   enddo
!   do ja=1,nexp
!      write(*,10)'3C before expanded: ',ja,(iset(la,ja),la=1,nconst)
!   enddo
! loop several times expanding one sublattice with wildcard each time
   ja=1
   lat1: do la=1,phlista(lokph)%noofsubl
      if(iset(la,1).eq.-99) then
         ka=1
         do while(ka.lt.nexp)
            jc=ja
            do jb=1,phlista(lokph)%nooffr(la)
               do jz=1,multi(la)
                  iset(la,ka)=jc
                  ka=ka+1
               enddo
               jc=jc+1
            enddo
         enddo
      endif
      ja=ja+phlista(lokph)%nooffr(la)
   enddo lat1
   wildloop=nexp
!   do ja=1,wildloop
!      write(*,10)'3C after expanded:  ',ja,(iset(la,ja),la=1,nconst)
!   enddo
1000 continue
   return
 end subroutine expand_wildcards
 
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine intsort
!\begin{verbatim}
 subroutine intsort(intc,nint,intx)
! This is just another stupid sorting subroutine   
! intc is not changed
   implicit none
   integer intc(*),intx(*),nint
!\end{verbatim} %+
   integer byte,jj
   if(nint.lt.2) then
      write(*,*)'*** ERROR: intsort called with too few constituents',nint
      stop
   endif
   do byte=1,nint
      intx(byte)=intc(byte)
   enddo
   do while(byte.gt.0)
! values in intx are never zero
      byte=0
      do jj=2,nint
         if(intx(jj-1).gt.intx(jj)) then
            byte=intx(jj)
            intx(jj)=intx(jj-1)
            intx(jj-1)=byte
         endif
      enddo
   enddo
1000 continue
   return
 end subroutine intsort

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine notallowlisting
!\begin{verbatim}
 logical function notallowlisting(privil)
! check if user is allowed to list data
   double precision privil
!\end{verbatim}
   logical ok
   ok=.TRUE.
   if(proda.ne.zero) then
      write(*,*)'3E Testing if listing allowed'
      if(privil.ne.proda) then
         ok=.FALSE.
      else
         write(*,*)'Sorry you are not allowed to list data'
      endif
   else
      ok=.FALSE.
   endif
   notallowlisting=ok
   return
 end function notallowlisting

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine lower_case_species_name
!\begin{verbatim}
 subroutine lower_case_species_name(constext,ip,isp)
! writes a species name using lower case for second letter of element
   implicit none
   character constext*(*)
   integer ip,isp
!\end{verbatim}
   integer iel,jp,lcase,kp
   character elsym*2,name*24
   jp=1
   name=' '
   lcase=ichar('a')-ichar('A')
   do iel=1,splista(isp)%noofel
      elsym=ellista(splista(isp)%ellinks(iel))%symbol
      kp=0
      if(elsym(2:2).ne.' ') then
         elsym(2:2)=char(ichar(elsym(2:2))+lcase)
         name(jp:)=elsym
         jp=jp+2
      else
         name(jp:)=elsym
         jp=jp+1
         kp=1
      endif
! 3rd argument 0 means no sign
      if(abs(splista(isp)%stoichiometry(iel)-one).gt.1.0D-6 .or. &
           (iel.lt.splista(isp)%noofel .and. kp.eq.1)) then
         call wrinum(name,jp,6,0,splista(isp)%stoichiometry(iel))
         if(buperr.ne.0) then
            write(*,*)'3C buperr 2: ',trim(name),buperr
            buperr=0
         endif
      endif
   enddo
! species may have a charge
   if(splista(isp)%charge.eq.one) then
      name(jp:jp)='+'
      jp=jp+1
   elseif(splista(isp)%charge.eq.-one) then
      name(jp:jp)='-'
      jp=jp+1
   elseif(abs(splista(isp)%charge).gt.1-0D-6) then
      call wrinum(name,jp,6,1,splista(isp)%charge)
   endif
!   write(*,*)'3C suck: lower case name: ',trim(name)
   constext(ip:)=name
   ip=len_trim(constext)+1
1000 continue
   return
 end subroutine lower_case_species_name

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

