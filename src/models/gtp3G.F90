!
! gtp3G included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     11. Status for things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine change_element_status
!\begin{verbatim}
 subroutine change_element_status(elname,nystat,ceq)
! change the status of an element, can affect species and phase status
! nystat:0=entered, 1=suspended, -1 special (exclude from sum of mole fraction)
!
! suspending elements for each equilibrium separately not yet implemented
!
   implicit none
   character elname*(*)
   integer nystat
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer iel,lokel
   call find_element_by_name(elname,iel)
   if(gx%bmperr.ne.0) goto 1000
   lokel=elements(iel)
   write(*,*)'3G Changing element status not yet implemented'
   goto 1000
   if(btest(ellista(iel)%status,elsus)) then
! element already suspended, quit it should be suspended again ....
      if(nystat.eq.1) goto 1000
! element status should be changed from suspended to entered
      ellista(iel)%status=ibclr(ellista(iel)%status,elsus)
      call restore_species_implicitly_suspended
      call restore_phases_implicitly_suspended
   elseif(nystat.eq.1) then
! element should be changed from entered to suspended
      ellista(iel)%status=ibset(ellista(iel)%status,elsus)
      call suspend_species_implicitly(ceq)
      call suspend_phases_implicitly(ceq)
   endif
1000 continue
   return
 end subroutine change_element_status

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function testelstat
!\begin{verbatim}
 logical function testelstat(iel,status)
! return value of element status bit
   implicit none
   integer iel,status
!\end{verbatim}
   integer lokel
   if(iel.gt.0 .and. iel.lt.noofel) then
      lokel=elements(iel)
      if(btest(ellista(lokel)%status,status)) then
! btest(iword,bit) .true. if bit set in iword
! iword=ibclr(iword,bit) to clear bit bit in iword
! iword=ibset(iword,bit) to set bit bit in iword
         testelstat=.true.
      else
         testelstat=.false.
      endif
   else
      gx%bmperr=4042
   endif
 end function testelstat

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine change_species_status
!\begin{verbatim}
 subroutine change_species_status(spname,nystat,ceq)
! change the status of a species, can affect phase status
! nystat:0=entered, 1=suspended
   implicit none
   integer nystat
   character spname*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer loksp
   call find_species_record(spname,loksp)
   if(gx%bmperr.ne.0) goto 1000
   write(*,*)'3G Changing species status not yet implemented'
   goto 1000
   if(btest(splista(loksp)%status,spsus)) then
! species already suspended, quit if it should be suspended again ....
      if(nystat.eq.1) goto 1000
! restore the species (and phases) unless implicitly suspended
      if(btest(splista(loksp)%status,spimsus)) then
! species cannot be entered as it is implicitly suspended (some element susp)
         gx%bmperr=4085; goto 1000
      endif
      splista(loksp)%status=ibclr(splista(loksp)%status,spsus)
      call restore_phases_implicitly_suspended
   elseif(nystat.eq.1) then
! suspend the species and possibly some phases
      splista(loksp)%status=ibset(splista(loksp)%status,spsus)
      call suspend_phases_implicitly(ceq)
   endif
1000 continue
   return
 end subroutine change_species_status

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function testspstat
!\begin{verbatim}
 logical function testspstat(isp,status)
! return value of species status bit
   implicit none
   integer isp,status
!\end{verbatim}
   integer loksp
   if(isp.gt.0 .and. isp.le.noofsp) then
      loksp=species(isp)
      if(btest(splista(loksp)%status,status)) then
! btest(iword,bit) .true. if bit set in iword
! iword=ibclr(iword,bit) to clear bit bit in iword
! iword=ibset(iword,bit) to set bit bit in iword
         testspstat=.true.
      else
         testspstat=.false.
      endif
   else
      gx%bmperr=4051
   endif
 end function testspstat

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function get_phase_status
!\begin{verbatim}
 integer function get_phase_status(iph,ics,text,ip,val,ceq)
! return phase status as text and amount formula units in val
! for entered and fix phases also phase amounts.
! OLD Function value: 1=entered, 2=fix, 3=dormant, 4=suspended, 5=hidden
   implicit none
   character text*(*)
   integer iph,ics,ip
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision val
!\end{verbatim} %+
   integer ists,lokph,lokcs,j
! write current status
   ists=0
   val=-one
   if(iph.gt.0 .and. iph.le.noph()) then
      call get_phase_compset(iph,ics,lokph,lokcs)
!old      if(btest(phlista(lokph)%status1,phhid)) then
!old         text='HIDDEN'; ip=6
!old         ists=5
!old      elseif(btest(ceq%phase_varres(lokcs)%status2,CSSUS)) then
!              entered,   fix,   suspended,   dormant
! bit setting: 00         01   , 10           11
!old           if(btest(ceq%phase_varres(lokcs)%status2,CSFIXDORM)) then
!old              text='DORMANT'; ip=7
!old              ists=3
!old           else
!old              text='SUSPENDED'; ip=9
!old              ists=4
!old           endif
!old      elseif(btest(ceq%phase_varres(lokcs)%status2,CSFIXDORM)) then
!old         text='FIXED'; ip=5
!         val=ceq%phase_varres(lokcs)%amount(1)
!old         val=ceq%phase_varres(lokcs)%amfu
!old         ists=2
!old      else
!old         text='ENTERED'; ip=7
!old         val=ceq%phase_varres(lokcs)%amfu
!old         ists=1
!old      endif
! new way, test PHSTATE
      j=ceq%phase_varres(lokcs)%phstate
!z      if(j.lt.-4 .or. j.gt.2) then
! I had an erroor here when plotting map2 macro because after the second
! map command I had 2 liquid compsets and during the first mapping I had
! only one liquid so I think
!z         ip=j
!z         j=0
!z         if(btest(ceq%phase_varres(lokcs)%status2,CSSUS)) then
!z            if(btest(ceq%phase_varres(lokcs)%status2,CSFIXDORM)) then
!z               j=-2
!z            else ! suspended
!z               j=3
!z            endif
!z         elseif(btest(ceq%phase_varres(lokcs)%status2,CSFIXDORM)) then
! fix
!z            j=2
!z         else ! entered
!z            j=0
!z         endif
! save this status .... ???
!z         write(*,16)'3G PHSTATE wrong, fixing ...',iph,ics,j,ip,&
!z              ceq%phase_varres(lokcs)%status2
!z         ceq%phase_varres(lokcs)%phstate=j
!z      endif
      select case(j)
      case default
         write(*,16)'3G: PHSTATE not correct: ',iph,ics,j,ip,&
              ceq%phase_varres(lokcs)%status2
16       format(a,4i3,2x,z16)
         gx%bmperr=4324
      case(phfixed) ! fix 2
         text='FIXED'
         ip=5
         val=ceq%phase_varres(lokcs)%amfu
         ists=phfixed
      case(-1,0,1) ! entered (unstable, unknown, stable)
         text='ENTERED'
         ip=7
         val=ceq%phase_varres(lokcs)%amfu
         ists=phentered
      case(phdorm) ! dormant -2
         text='DORMANT'
         ip=7
         ists=phdorm
      case(phsus) ! suspended -3
         text='SUSPENDED'
         ip=9
         ists=phsus
      case(phhidden) ! hidden -4
         text='HIDDEN'
         ip=6
         ists=phhidden
      end select
   else
!      write(*,*)'No such phase'
      gx%bmperr=4050; goto 1000
   endif
   get_phase_status=ists
!   write(*,*)'3G: PHSTAT value: ',ists
!   write(*,*)'3G: gps: ',ip
1000 continue
   return
 end function get_phase_status

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function test_phase_status
!\begin{verbatim}
 integer function test_phase_status(iph,ics,val,ceq)
! Almost same as get_..., returns phase status as function value but no text
! value is amfu
! nystat:-4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fix
! this is different from in change_phase .... one has to make up one's mind
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer iph,ics
   double precision val
!\end{verbatim}
   integer ists,lokph,lokcs,j,ip
   character text*24
! new code
   ists=0
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 900
   ists=ceq%phase_varres(lokcs)%phstate
   val=ceq%phase_varres(lokcs)%amfu
   goto 900
!============================================= code below redundant?
   ists=0
   ip=1
   val=-one
   ists=get_phase_status(iph,ics,text,ip,val,ceq)
   goto 900
!------------------
   if(iph.gt.0 .and. iph.le.noph()) then
      call get_phase_compset(iph,ics,lokph,lokcs)
! biet set means false ....
!z      if(btest(phlista(lokph)%status1,phhid)) then
! hidden
!z         ists=5
!z      elseif(btest(ceq%phase_varres(lokcs)%status2,CSSUS)) then
!              entered,   fix,   suspended,   dormant
! bit setting: 00         01   , 10           11
!z           if(btest(ceq%phase_varres(lokcs)%status2,CSFIXDORM)) then
!z              ists=3
!z           else
!z              ists=4
!z           endif
!z      elseif(btest(ceq%phase_varres(lokcs)%status2,CSFIXDORM)) then
!z         val=ceq%phase_varres(lokcs)%amfu
!z         ists=2
!z      else
!z         ists=1
!z         val=ceq%phase_varres(lokcs)%amfu
!z      endif
! new way, test PHSTATE
      j=ceq%phase_varres(lokcs)%phstate
      select case(ceq%phase_varres(lokcs)%phstate)
      case default
         write(*,*)'PHSTAT outside range -4:2: ',j
      case(phfixed) ! fix +2
         if(ists.ne.2) write(*,*)'wrong PHSTAT',ists,j
      case(-1,0,1) ! entered (unstable, unknown, stable)
         if(ists.ne.1) write(*,*)'wrong PHSTAT',ists,j
      case(phdorm) ! dormant -2
         if(ists.ne.3) write(*,*)'wrong PHSTAT',ists,j
      case(phsus) ! suspended -3
         if(ists.ne.4) write(*,*)'wrong PHSTAT',ists,j
      case(phhidden) ! hidden -4
         if(ists.ne.5) write(*,*)'wrong PHSTAT',ists,j
      end select
   else
!      write(*,*)'No such phase'
      gx%bmperr=4050; goto 1000
   endif
900 continue
   test_phase_status=ists
1000 continue
   return
 end function test_phase_status

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_phase_status_bit
!\begin{verbatim}
 subroutine set_phase_status_bit(lokph,bit)
! set the status bit "bit" in status1, cannot be done outside this module
! as the phlista is private
! These bits do not depend on the composition set
   implicit none
   integer lokph,bit
!\end{verbatim} %+
   integer lokcs,j
   if(bit.lt.0 .or. bit.gt.31) then
      write(*,*)'Illegal phase bit number'
      gx%bmperr=4325; goto 1000
   elseif(lokph.le.0 .or. lokph.gt.noofph) then
      write(*,*)'Illegal phase in call to set_phase_status_bit'
      gx%bmperr=4326; goto 1000
   endif
!   write(*,99)'sphs1bit: ',lokph,bit,phlista(lokph)%status1
99 format(a,2i3,z8)
   phlista(lokph)%status1=ibset(phlista(lokph)%status1,bit)
   if(bit.eq.PHHID) then
! if bit is PHHID, i.e. hidden, set PHSTATE in all phase_varres record to -4
      do j=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(j)
! eventually, this must be set in all equilibrium records now just firsteq ??
         firsteq%phase_varres(lokcs)%phstate=-4
      enddo
   endif
1000 continue
   return
 end subroutine set_phase_status_bit

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine clear_phase_status_bit
!\begin{verbatim} %-
 subroutine clear_phase_status_bit(lokph,bit)
! clear the status bit "bit" in status1, cannot be done outside this module
! as the phlista is private
   implicit none
   integer lokph,bit
!\end{verbatim} %+
   integer lokcs,j
   if(bit.lt.0 .or. bit.gt.31) then
      write(*,*)'Illegal phase bit number'
      gx%bmperr=4325; goto 1000
   endif
   phlista(lokph)%status1=ibclr(phlista(lokph)%status1,bit)
   if(bit.eq.PHHID) then
      write(*,*)'clear_bit: Not implemented to change PHSTATE'
! if bit is PHHID, i.e. hidden, set PHSTATE in all phase_varres record to 0
      do j=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(j)
! eventually, this must be set in all equilibrium records now just firsteq ??
         firsteq%phase_varres(lokcs)%phstate=phentered
      enddo
   endif
1000 continue
   return
 end subroutine clear_phase_status_bit

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function test_phase_status_bit
!\begin{verbatim} %-
 logical function test_phase_status_bit(iph,ibit)
! return TRUE is status bit ibit for  phase iph, is set
! because phlista is private.  Needed to test for gas, ideal etc, 
! DOES NOT TEST STATUS like entered/fixed/dormant/suspended
   implicit none
   integer iph,ibit
!\end{verbatim}
   integer lokph
   if(iph.gt.0 .and. iph.le.noofph) then
      lokph=phases(iph)
   else
      gx%bmperr=4050; goto 1000
   endif
   if(btest(phlista(lokph)%status1,ibit)) then
      test_phase_status_bit=.true.
   else
      test_phase_status_bit=.false.
   endif
1000 continue
   return
 end function test_phase_status_bit

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine change_many_phase_status
!\begin{verbatim}
 subroutine change_many_phase_status(phnames,nystat,val,ceq)
! change the status of many phases. 
! nystat:-4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fix
! phnames is a list of phase names or *S (all suspeded) *D (all dormant) or
! *E (all entered (stable, unknown, unstable), *U all unstable
! If just * then change_phase_status is called directly
! It calls change_phase_status for each phase
   implicit none
   character phnames*(*)
   integer nystat
   double precision val
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer qph,ics,oldstat,ipos,slen,lokph,lokcs,phcsx
   character name*24
! CCI correted size of phnames
!   write(*,*)'3G phnames: ',trim(phnames),' >',phnames(1:1),'<'
   if(phnames(1:1).eq.'*') then
!      write(*,*)'3G star'
      if(phnames(2:2).eq.'S') then
         oldstat=-3
      elseif(phnames(2:2).eq.'D') then
         oldstat=-2
      elseif(phnames(2:2).eq.'E') then
! all entered (stable, unstable, unknown)
         oldstat=0
      elseif(phnames(2:2).eq.'U') then
! all unstable phases (not those which ar efix!)
         oldstat=1
      elseif(phnames(2:2).eq.' ') then
         qph=-1
!         write(*,*)'3G star',qph,ics,nystat,val
         call change_phase_status(qph,ics,nystat,val,ceq)
         goto 1000
      else
         write(*,*)'Illegal selection of old phase status after *'
         gx%bmperr=4327; goto 1000
      endif
! loop for all phases to find those with correct old status
      do qph=1,noofph
! we cannot loop for ics as we do not know lokph
         ics=1
         call get_phase_compset(qph,ics,lokph,lokcs)
200      continue
! stable phases has ceq%phase_varres(lokcs)%phstate = 1
! fix phases =2
         ipos=oldstat-ceq%phase_varres(lokcs)%phstate
!         write(*,*)'3G entered: ',qph,ics,oldstat,ipos
         if((oldstat.ne.1 .and. ipos.eq.0) .or. &
              (oldstat.eq.0 .and. abs(ipos).eq.1) .or.&
              (oldstat.eq.1 .and. ipos.gt.0)) then
! *U=nystat means all all with phstate <=0 that means ipos=1-0; 1-(-1)=2 etc
!              (oldstat.eq.1 .and. abs(ipos).gt.0)) then
! this comp.set has correct old phase status
            call change_phase_status(qph,ics,nystat,val,ceq)
            if(gx%bmperr.ne.0) goto 1000
         endif
! take next composition set if any, else next phase
         if(ics.lt.phlista(lokph)%noofcs) then
            ics=ics+1
            lokcs=phlista(lokph)%linktocs(ics)
            goto 200
         endif
      enddo
   else
! we have one or more specific phase names separated by space or comma
! ipos is updated inside getext, The 3rd argument of getext is JTYP
!      JTYP DEFINES THE TERMINATION OF A STRING
!      1 TEXT TERMINATED BY SPACE OR ","
!      2 TEXT TERMINATED BY SPACE
!      3 TEXT TERMINATED BY ";" OR "."
!      4 TEXT TERMINATED BY ";"
!      5 TEXT UP TO END-OF-LINE
!      6 TEXT UP TO AND INCLUDING ";"
!    >31, THE CHAR(JTYP) IS USED AS TERMINATING CHARACTER
      ipos=0
500   continue
         call getext(phnames,ipos,1,name,' ',slen)
!         write(*,*)'3G phase1: ',slen,' ',name
         if(name(1:1).eq.' ') goto 1000
!         write(*,*)'3G phase2: ',name
!         call find_phase_by_name(name,qph,ics)
! phcsx=-1 means that all composition sets should have new status
         phcsx=-1
         call find_phasex_by_name(name,phcsx,qph,ics)
         if(gx%bmperr.ne.0) then
!            write(*,*)'No phase called "',name(1:len_trim(name)),'"'
            gx%bmperr=0
         else
! we may have to make a loop for all composition sets
! A phase without composition set specification but with several composition 
! sets should have all composition sets changed to the new status
! UNLESS the status is FIX
            if(ics.lt.0) then
! we should never loop to set all composition sets to FIXED
! if another composition set than 1 was to be set fixed ics is not negative
               if(nystat.eq.PHFIXED) then
                  slen=1
               else
                  slen=-ics
               endif
!             write(*,*)'3G Status changed for several composition sets: ',slen
               do ics=1,slen
                  call change_phase_status(qph,ics,nystat,val,ceq)
                  if(gx%bmperr.ne.0) goto 1000
               enddo
            else
!               write(*,*)'3G changing status for a single phase',nystat
               call change_phase_status(qph,ics,nystat,val,ceq)
               if(gx%bmperr.ne.0) goto 1000
            endif
         endif
         goto 500
   endif
1000 continue
   return
 end subroutine change_many_phase_status

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine change_phtup_status
!\begin{verbatim} %-
 subroutine change_phtup_status(phtupx,nystat,val,ceq)
! change the status of a phase tuple. Also used when setting phase fix etc.
! nystat:-4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fix
! qph can be -1 meaning all or a specifix phase index. ics compset
! 
   implicit none
   integer phtupx,nystat
   double precision val
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,iph,ics
   if(phtupx.lt.0) then
! change status for all phases to nystat
      call change_many_phase_status('* ',nystat,val,ceq)
   else
!      lokph=phasetuple(phtupx)%phaseix
      lokph=phasetuple(phtupx)%lokph
      iph=phlista(lokph)%alphaindex
      ics=phasetuple(phtupx)%compset
!   write(*,77)'3G Test: ',phlista(lokph)%name,phtupx,lokph,iph,phases(iph)
!77 format(a,a,10i5)
      call change_phase_status(iph,ics,nystat,val,ceq)
   endif
1000 continue
   return
 end subroutine change_phtup_status

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine change_phase_status
!\begin{verbatim} %-
 subroutine change_phase_status(qph,ics,nystat,val,ceq)
! change the status of a phase. Also used when setting phase fix etc.
! nystat:-4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fix
! qph can be -1 meaning all or a specifix phase index. ics compset
! 
   implicit none
   integer qph,ics,nystat
   double precision val
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,lokcs,iph,ip,mcs
   character line*80,phname*32
!   write(*,11)'3G in change_phase_status: ',qph,ics,nystat,val
11 format(a,3i5,1pe14.6)
   if(qph.eq.-1) then
! this means all phases. All phases cannot be set fix
      if(nystat.eq.3) then
         gx%bmperr=4152; goto 1000
      endif
      iph=1
      ics=1
   else
! a specific phase
      iph=qph
   endif
! return here for next phase
100 continue
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3G: Phase and status: ',iph,ceq%phase_varres(lokcs)%phstate
   if(ceq%phase_varres(lokcs)%phstate.eq.phfixed) then
! this phase and composition set is fix, remove condition
! unless the new status is also FIX
      if(nystat.ne.phfixed) then
         call get_phase_name(iph,ics,phname)
         line=' NOFIX='//phname(1:len_trim(phname))
         ip=1
!         write(*,*)'Remove fix phase: ',line(1:len_trim(line))
         call set_condition(line,ip,ceq)
         if(gx%bmperr.ne.0) then
!            write(*,*)'Failed to remove fix phase as condition'
            goto 1000
         endif
      endif
   endif
   bigif: if(ceq%phase_varres(lokcs)%phstate.eq.phhidden) then
! phase is hidden, quit if it should be hidden again
!   bigif: if(btest(phlista(lokph)%status1,phhid)) then
      if(nystat.eq.phhidden) goto 900
!      phlista(lokph)%status1=ibclr(phlista(lokph)%status1,phhid)
!??? this phase must be added in phlista ??? no it is already there ???
      write(*,*)'Unifished handling of hide/not hide ...'
      gx%bmperr=4095; goto 900
   elseif(nystat.eq.phhidden) then
! phase is not hidden but should be set as hidden,
! Always applies to all composition sets
! clear all entered/suspended/dormant/fix for all composition sets
      phlista(lokph)%status1=ibset(phlista(lokph)%status1,phhid)
      do mcs=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(mcs)
         ceq%phase_varres(lokcs)%phstate=phhidden         
! also set amounts and dgm to zero
         ceq%phase_varres(lokcs)%amfu=zero
         ceq%phase_varres(lokcs)%netcharge=zero
         ceq%phase_varres(lokcs)%dgm=zero
      enddo
   else !bigif
      lokcs=phlista(lokph)%linktocs(ics)
! changing FIX/ENTERED/SUSPENDED/DORMANT for a composition set
! input nystat:0=entered, 3=fix, 1=suspended, 2=dormant
! bit setting: 00         01   , 10           11  !! BITS NO LONGER USED
!      write(*,71)'3G new status: ',iph,ics,lokph,lokcs,nystat,phentered,val
71    format(a,6i5,1pe14.6)
      if(nystat.eq.phentered .or. nystat.eq.phentunst .or. &
           nystat.eq.phentstab) then
! set enterered with amount val and dgm zero
!         write(*,*)'Setting phase as entered',nystat
!         ceq%phase_varres(lokcs)%phstate=phentered
         ceq%phase_varres(lokcs)%phstate=nystat
         ceq%phase_varres(lokcs)%amfu=val
         ceq%phase_varres(lokcs)%netcharge=zero
         ceq%phase_varres(lokcs)%dgm=zero
      elseif(nystat.eq.phsus) then
! set suspended with amount and dgm zero
         ceq%phase_varres(lokcs)%phstate=phsus
         ceq%phase_varres(lokcs)%amfu=zero
         ceq%phase_varres(lokcs)%netcharge=zero
         ceq%phase_varres(lokcs)%dgm=zero
      elseif(nystat.eq.phdorm) then
! set dormant with amount and dgm zero
         ceq%phase_varres(lokcs)%phstate=phdorm
         ceq%phase_varres(lokcs)%amfu=zero
         ceq%phase_varres(lokcs)%netcharge=zero
         ceq%phase_varres(lokcs)%dgm=zero
      elseif(nystat.eq.phfixed) then
! to allow MAPPHASEFIX=3
         ceq%phase_varres(lokcs)%phstate=phfixed
         ceq%phase_varres(lokcs)%amfu=val
         ceq%phase_varres(lokcs)%netcharge=zero
         ceq%phase_varres(lokcs)%dgm=zero
! also set as condition
         call get_phase_name(iph,ics,phname)
         line=' FIX='//phname(1:len_trim(phname))//' =='
         ip=len_trim(line)+2
         call wrinum(line,ip,6,0,val)
         if(buperr.ne.0) goto 1000
         ip=1
!         write(*,*)'phase fix condition: ',line(1:40)
         call set_condition(line,ip,ceq)
      endif
   endif bigif
900 continue
! check if loop
   if(qph.eq.-1) then
      lokph=phases(iph)
      if(ics.lt.phlista(lokph)%noofcs) then
         ics=ics+1
      elseif(iph.lt.noofph) then
         iph=iph+1
         ics=1
      else
         goto 1000
      endif
      goto 100
   endif
1000 continue
!   write(*,*)'error code: ',gx%bmperr
   return
 end subroutine change_phase_status

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine mark_stable_phase
!\begin{verbatim} &-
 subroutine mark_stable_phase(iph,ics,ceq)
! change the status of a phase. Does not change fix status
! called from meq_sameset to indicate stable phases (nystat=1)
! nystat:-4 hidden, -3 suspended, -2 dormant, -1,0,1 entered, 2 fix
! 
   implicit none
   integer iph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokph,lokcs
!   write(*,11)'3G mark as stable: ',iph,ics,phentstab
11 format(a,3i5,1pe14.6)
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3G: Phase and status: ',iph,ceq%phase_varres(lokcs)%phstate
   if(ceq%phase_varres(lokcs)%phstate.eq.phhidden) then
      write(*,*)'Error calling mark_stable for hidden phase'
      gx%bmperr=4095; goto 1000
   elseif(ceq%phase_varres(lokcs)%phstate.le.phdorm) then
      write(*,*)'Cannot make suspended or doremant phases as stable'
      gx%bmperr=4095; goto 1000
   elseif(ceq%phase_varres(lokcs)%phstate.eq.phfixed) then
! do nothing
      goto 1000
   else
      ceq%phase_varres(lokcs)%phstate=phentstab
   endif
1000 continue
!   write(*,*)'error code: ',gx%bmperr
   return
 end subroutine mark_stable_phase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     12. Unfinished things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_unit
!\begin{verbatim}
 subroutine set_unit(property,unit)
! set the unit for a property, like K, F or C for temperature
! >>>> unfinished
   implicit none
   character*(*) property,unit
!\end{verbatim}
   write(*,*)'Not implemented yet'
   gx%bmperr=4078
1000 continue
   return
 end subroutine set_unit

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine save_results
!\begin{verbatim}
 subroutine save_results(lut,iph,ics,long)
! write calculated results for a phase for later use in POST
   implicit none
   integer lut,iph,ics,long
!\end{verbatim}
   write(*,*)'Not implemented yet'
   gx%bmperr=4078
! header with abbreviations
!    call list_abbrev(lut)
! first conditions ...
!    call list_conditions(lut)
! Global values of G, N, V etc
!    call list_global_results(lut)
! Element data
!    call list_components_results(lut)
! Phases and composition sets
!   call dump_phase_results(lut)
1000 return
 end subroutine save_results

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_constituent_reference_state
!\begin{verbatim}
 subroutine set_constituent_reference_state(iph,icon,asum)
! determine the end member to calculate as reference state for this constituent
! Used when giving a chemical potential for a constituent like MU(GAS,H2O)
   implicit none
   integer iph,icon
   double precision asum
!\end{verbatim}
   type(gtp_endmember), pointer :: endmemrec
   integer lokph,nsl,ll,jcon,loksp,loksp2,lokcs
!
   lokph=phases(iph)
   loksp=phlista(lokph)%constitlist(icon)
   nsl=phlista(lokph)%noofsubl
   endmemrec=>phlista(lokph)%ordered
   asum=one
   lokcs=phlista(lokph)%linktocs(1)
   if(nsl.eq.1) then
      asum=firsteq%phase_varres(lokcs)%sites(1)
      emlist1: do while(associated(endmemrec))
         if(endmemrec%fraclinks(nsl,1).eq.icon) goto 300
         endmemrec=>endmemrec%nextem
      enddo emlist1
   else
! several sublattices OK if same species or vacancies in other sublattices
      asum=zero
      emlist2: do while(associated(endmemrec))
         do ll=1,nsl
            jcon=endmemrec%fraclinks(ll,1)
            if(jcon.ne.icon) then
               loksp2=phlista(lokph)%constitlist(jcon)
               if(loksp2.eq.loksp) then
! same species in this sublattice, add sites to asum
                  asum=asum+firsteq%phase_varres(lokcs)%sites(ll)
               elseif(.not.btest(splista(loksp2)%status,spva)) then
! other species (not vacancies) in this sublattice, skip this end member
                  goto 200
               endif
            else
               asum=asum+firsteq%phase_varres(lokcs)%sites(ll)
            endif
         enddo
! this endmember OK
         goto 300
! not this end member
200       continue
         endmemrec=>endmemrec%nextem
      enddo emlist2
   endif
! this phase cannot exist for species icon as pure
   gx%bmperr=4112; goto 1000
300 continue
1000 continue
   return
 end subroutine set_constituent_reference_state

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine elements2components1
!\begin{verbatim}
 subroutine elements2components1(nspel,dum,ncmp,cmpstoi,ceq)
! converts a stoichiometry array for a species from elements to components
! This subroutine, is it used to get activity for a constituent in gtp3F
! dum is is no longer used
   implicit none
   integer nspel,ncmp
! cmpstoi is stoichiometry as element, changed to be as components
   double precision cmpstoi(*),dum(*)
   double precision, allocatable :: stoi(:)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   double precision, parameter :: small=1.0d-30
   integer ic,jc,ns
   allocate(stoi(noofel))
   do ic=1,noofel
      stoi(ic)=cmpstoi(ic)
      cmpstoi(ic)=zero
   enddo
! use the ceq%complist(ic)%invcompstoi
!   do ic=1,noofel
!      cmpstoi(ic)=zero
!   enddo
! not sure about the indices here .... ????
!   write(*,*)'e2c: ',noofel,nspel,stoi(1),ceq%invcompstoi(1,1)
   do ic=1,noofel
! convert elements to components, if the elements are components no problem
      do jc=1,noofel
         cmpstoi(ic)=cmpstoi(ic)+ceq%invcompstoi(ic,jc)*stoi(ic)
      enddo
   enddo
!   write(*,7)'3G 1: ',(stoi(ic),ic=1,noofel)
!   write(*,7)'3G 2: ',(stoi(ic),ic=1,noofel)
7  format(a,10(1pe12.4))
! MODIFIED HERE 190710/BoS, return stoichiometry for ALL components
   ncmp=noofel
   goto 1000
!---------------------
!  skip code below ...   
   ncmp=0
   ic=0
   ns=0
200 continue
   ic=ic+1
   if(ic.lt.noofel) then
      if(abs(cmpstoi(ic)).lt.small) then
         do jc=ic,noofel
            cmpstoi(jc)=cmpstoi(jc+1)
         enddo
      else
         ncmp=ncmp+1
!         write(*,*)'c2c1: ',ic,ncmp
      endif
      goto 200
   elseif(abs(cmpstoi(ic)).gt.small) then
!      write(*,*)'c2c2: ',ic,ncmp,cmpstoi(ic)
      ncmp=ncmp+1
   endif
!   write(*,190)ic,(cmpstoi(i),i=1,ncmp)
!190 format('e2c3: ',i3,10F7.3)
1000 continue
   return
 end subroutine elements2components1

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\!
!>     13. Internal stuff
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
 
!\addtotable subroutine termterm
!\begin{verbatim}
 subroutine termterm(string,ich,kpos,lpos,value)
! search for first occurance of + - = > or <
! if + or - then also extract possible value after sign
! value is coefficient for NEXT term (if any)
! IF WE FIND A ( accept all characters up to ), constitunets can have + or -
! kpos is last character in THIS state variable, lpos where NEXT may start
   implicit none
   character string*(*)
   integer kpos,ich,lpos
   double precision value
!\end{verbatim}
   integer ipos,jpos,i1
   logical afterlp
   character ch1*1
   character (len=1), dimension(6), parameter :: chterm=&
        ['+','-','=','<','>',':']
!
   afterlp=.FALSE.
   ich=0
   sloop: do ipos=1,len_trim(string)
      ch1=string(ipos:ipos)
! I do not check for nested ( ) or ...
      if(ch1.eq.'(') then
         afterlp=.TRUE.
      elseif(ch1.eq.')') then
         afterlp=.FALSE.
      endif
! accept all characters between ( )
      if(afterlp) cycle sloop
      do i1=1,6
         if(ch1.eq.chterm(i1)) then
            kpos=ipos; ich=i1; exit sloop
         endif
      enddo
   enddo sloop
! different actions depending on ich
!   write(*,17)'3G termterm: ',trim(string),string(1:kpos),ich,kpos
17 format(a,' "',a,'" >',a,'< ',2i3)
   select case(ich)
   case default
      write(*,*)'3G wrong ich case: ',ich
   case(0)
! no terminator, just return with position pointer after the text
      continue
      kpos=len_trim(string)+1
   case(1,2)
! there is a - or + sign, collect value in front of next term
      lpos=kpos+1
      call getrel(string,lpos,value)
      if(buperr.ne.0) then
! a sign not followed by number means unity
         buperr=0; value=one
      else
! lpos first character after number, a number must be followed by a "*"
         if(string(lpos:lpos).ne.'*') then
            write(*,*)'3G syntax error missing *: ',string(1:lpos+5),lpos
            gx%bmperr=4130
         else
            lpos=lpos+1
         endif
      endif
      if(ich.eq.2) value=-value
   case(3,4,5)
! there is an = sign, or > or <, just set back the pointer
      kpos=ipos
      lpos=0
   case(6)
! there is an : sign, meaning a condition number, must be followed by =
      if(string(kpos+1:kpos+1).ne.'=') then
         gx%bmperr=4328; goto 1000
      endif
      kpos=ipos+1
      lpos=0
   end select
1000 continue
!   write(*,17)'3G termterm: ',trim(string),string(1:lpos),ich,lpos
   return
 end subroutine termterm

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine alphaelorder
!\begin{verbatim}
 subroutine alphaelorder
! arrange new element in alphabetical order
! also make alphaindex give alphabetical order
   implicit none
!\end{verbatim} %+
   character symb1*2
   integer i,j
   symb1=ellista(noofel)%symbol
!  write(6,*)'alphaelorder 1: ',symb1,noofel
   loop1: do i=1,noofel-1
      if(symb1.lt.ellista(elements(i))%symbol) then
         loop2: do j=noofel,i+1,-1
            elements(j)=elements(j-1)
            ellista(elements(j))%alphaindex=j
         enddo loop2
!        write(6,*)'alphaelorder 3: ',i
         elements(i)=noofel
         ellista(elements(i))%alphaindex=i
         exit
      endif
   enddo loop1
!  write(6,*)'alphaelorder 4: ',(elements(k),k=1,noofel)
 END subroutine alphaelorder

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine alphasporder
!\begin{verbatim}
 subroutine alphasporder
! arrange new species in alphabetical order
! also make alphaindex give alphabetical order
   implicit none
!\end{verbatim} %+
   character symb1*24
   integer i,j
   symb1=splista(noofsp)%symbol
!  write(6,*)'alphasporder 1: ',symb1(1:6),noofsp
   loop1: do i=1,noofsp-1
      if(symb1.lt.splista(species(i))%symbol) then
!        write(6,*)'alphasporder 2; ',symb1,splista(species(i))%symbol
         loop2: do j=noofsp,i+1,-1
            species(j)=species(j-1)
            splista(species(j))%alphaindex=j
         enddo loop2
         species(i)=noofsp
!        write(6,*)'alphasporder 3:',i
         splista(species(i))%alphaindex=i
         exit
      endif
   enddo loop1
!  write(6,*)'alphasporder 4: ',(species(k),k=1,noofsp)
 END subroutine alphasporder

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine alphaphorder
!\begin{verbatim}
 subroutine alphaphorder(tuple)
! arrange last added phase in alphabetical order
! also make alphaindex give alphabetical order
! phletter G and L and I have priority
! tuple is returned as position in phase tuple
   implicit none
   integer tuple
!\end{verbatim}
   character symb1*24,ch1*1,ch2*1
   integer iph,lokph,j,lokcs
!
   symb1=phlista(noofph)%name
   ch1=phlista(noofph)%phletter
! one more phase in "phases" array
   phases(noofph)=noofph
!  write(6,75)'alphaphorder 1: ',noofph,ch1,symb1(1:6)
!75 format(A,I3,1x,A,1x,A)
   loop1: do iph=1,noofph-1
      lokph=phases(iph)
      ch2=phlista(lokph)%phletter
!      write(6,76)'alphaphorder 2A: ',iph,lokph,ch1,ch2
76 format(A,2I3,1x,A,1x,A)
! phaseletter different, if ch1=G insert it here
      if(ch1.eq.'G') goto 300
      if(ch2.eq.'G') goto 200
      liquid: if(ch1.eq.'L') then
         if(ch2.eq.'G') goto 200
         if(ch2.eq.'L') goto 100
         goto 300
      endif liquid
      if(ch2.eq.'L') goto 200
      solution: if(ch1.eq.'S') then
         if(ch2.eq.'G' .or. ch2.eq.'L') goto 200
         if(ch2.eq.'S') goto 100
         goto 300
      endif solution
      if(ch2.eq.'S') goto 200
      compound: if(ch1.eq.'C') then
         if(ch2.eq.'C') goto 100
         goto 200
      endif compound
! here phletter of lokph and the new phase are the same
100   continue
!     write(6,*)'alphaphorder 2B: ',symb1,phlista(lokph)%name
      if(symb1.lt.phlista(lokph)%name) goto 300
200    continue
   enddo loop1
! exit loop, add new phase last
!   lokph=phases(noofph)
   iph=phases(noofph)
300 continue
!   write(*,*)'3G new phase position: ',iph
!  write(6,77)'alphaphorder 2C: ',iph,lokph,phlista(lokph)%name
!77 format(A,2I3,1X,A)
! insert phase here at iph, shift down trailing phase indices
! also OK if new phase should be last
   loop2: do j=noofph,iph+1,-1
! update index of trailing phases, loop from the end not to overwrite
      phases(j)=phases(j-1)
      phlista(phases(j))%alphaindex=j
   enddo loop2
! index of new phase
!  write(6,*)'alphaphorder 4: ',lokph,iph,noofph
   phases(iph)=noofph
   phlista(noofph)%alphaindex=iph
   nooftuples=nooftuples+1
   tuple=iph
!   write(*,771)iph,phasetuple(iph),phlista(noofph)%name
771 format('3G tuple: ',i5,': ',4(i8,1x),2x,a)
! link to first compset set when phase_varres record connected
!   write(*,777)'3G phase tuple position: ',iph,noofph,lokph,lokcs,tuple
777 format(a,10i5)
   return
 END subroutine alphaphorder

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine check_alphaindex
!\begin{verbatim}
 subroutine check_alphaindex
! just for debugging, check that ellist(i)%alphaindex etc is  correct
   implicit none
!\end{verbatim}
   integer i,j,k,l
   write(kou,*)
   write(6,77)(ellista(elements(i))%symbol,i=1,noofel)
77  format(20(1x,A2))
   write(6,78)(splista(species(i))%symbol,i=1,noofsp)
78  format(20(1x,a6))
   write(6,*)'element alphaindex'
   check1:  do i=1,noofel
      j=ellista(elements(i))%alphaindex
      write(6,*)i,j,elements(i),ellista(i)%symbol
   enddo check1
   write(6,*)'species alphaindex'
   check2: do i=1,noofsp
      j=species(i)
      k=splista(j)%alphaindex
      l=splista(species(j))%alphaindex
      write(6,79)i,k,j,l,splista(j)%symbol
   enddo check2
79  format(4i4,1x,A)
   check3: do i=1,noofsp
      write(6,*)i,splista(i)%alphaindex,splista(i)%symbol
   enddo check3
 END subroutine check_alphaindex

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_constitlist
!\begin{verbatim}
 subroutine create_constitlist(constitlist,nc,klist)
! creates a constituent list ...
   implicit none
   integer, dimension(*) :: klist
   integer, dimension(:), allocatable :: constitlist
   integer nc
!\end{verbatim}
   integer ic
   ALLOCATE(constitlist(nc))
   DO ic=1,nc
      constitlist(ic)=klist(ic)
   enddo
   return
 END subroutine create_constitlist

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_parrecords
!\begin{verbatim}
 subroutine create_parrecords(lokph,lokcs,nsl,nc,nprop,iva,ceq)
! fractions and results arrays for a phase for parallel calculations
! location is returned in lokcs
! nsl is sublattices, nc number of constituents, nprop max number if propert,
! iva is an array which is set as constituent status word (to indicate VA)
! ceq is always firsteq ???
!
! BEWARE not adopted for threads
!
! >>> changed all firsteq below to ceq????
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer, dimension(*) :: iva
   integer lokph,lokcs, nsl, nc, nprop
!\end{verbatim}
   integer ic,nnc
! find free record, free list csfree maintained in FIRSTEQ only!
!   write(*,*)'3G maxcalcprop: ',nprop
   lokcs=csfree
   if(csfree.le.0) then
! This means no free phase_varres records.
! csfree is set to -1 by the statement csfree=phase_varres(lokcs)%next below
! when reserving the last free record.  The same for the other free lists
      gx%bmperr=4094; goto 1000
   endif
! the free list of phase_varres record is only maintained in firsteq
! but all equilibria have identical allocation of phase_varres records
! the free list is created when starting OC, each record points to the next
! After composition sets has been entered and deleted it may be different
! highcs should always be the index of the highest used record
   csfree=firsteq%phase_varres(lokcs)%nextfree
!   write(*,*)'3G looking for free varres record 1:',lokcs,csfree
! wrong ...   if(csfree.gt.highcs) highcs=csfree
! The varres record used will be, csfree is updated
!   write(*,*)'3G looking for free varres record 2:',lokcs,csfree
   firsteq%phase_varres(lokcs)%nextfree=0
   firsteq%phase_varres(lokcs)%status2=0
   ic=newhighcs(.true.)
   if(lokcs.gt.highcs) highcs=lokcs
! added integer status array constat. Set CONVA bit from iva array
!   write(*,*)'Allocate constat 2: ',nc,lokcs
   if(.not.allocated(ceq%phase_varres(lokcs)%constat)) then
! remove CSDEL if set
   firsteq%phase_varres(lokcs)%status2=&
        ibclr(firsteq%phase_varres(lokcs)%status2,CSDEL)
! already allocated error for the Al-Ni case, why?
! Maybe if composition set has been deleted without releasing allocated arrays?
      allocate(ceq%phase_varres(lokcs)%constat(nc))
   endif
!   write(*,*)'3G compset: ',trim(phlista(lokph)%name),nc,lokcs,&
!        size(ceq%phase_varres(lokcs)%constat)
!    write(*,33)nc,(iva(i),i=1,nc)
   do ic=1,nc
      ceq%phase_varres(lokcs)%constat(ic)=iva(ic)
   enddo
! allocate fraction and default fraction arrays
   allocate(ceq%phase_varres(lokcs)%yfr(nc))
   allocate(ceq%phase_varres(lokcs)%mmyfr(nc))
   do ic=1,nc
      ceq%phase_varres(lokcs)%yfr(ic)=one
      ceq%phase_varres(lokcs)%mmyfr(ic)=zero
   enddo
!   write(*,*)'Allocated mmyfr: ',lokcs,nc,nprop
! abnorm initiated to unity to avoid trouble at first calculation
   ceq%phase_varres(lokcs)%abnorm=one
   allocate(ceq%phase_varres(lokcs)%sites(nsl))
!
   if(btest(phlista(lokph)%status1,PHIONLIQ)) then
! for ionic liquid the sites may depend on composition
! I get error these already allocated. Why ??
      if(.not.allocated(ceq%phase_varres(lokcs)%dpqdy)) then
         allocate(ceq%phase_varres(lokcs)%dpqdy(nc))
         allocate(ceq%phase_varres(lokcs)%d2pqdvay(nc))
      endif
   endif
!
! result arrays for a phase for use in parallel processing
   ceq%phase_varres(lokcs)%nprop=nprop
   allocate(ceq%phase_varres(lokcs)%listprop(nprop))
   allocate(ceq%phase_varres(lokcs)%gval(6,nprop))
!   write(*,*)'Allocated gval: ',nprop,nc
   allocate(ceq%phase_varres(lokcs)%dgval(3,nc,nprop))
   nnc=nc*(nc+1)/2
!   write(*,*)'Allocated dgval: ',nprop,nc,nnc
   allocate(ceq%phase_varres(lokcs)%d2gval(nnc,nprop))
!   write(*,*)'Allocated d2gval: ',nprop,nc,nnc
! zero everything
   ceq%phase_varres(lokcs)%listprop=0
!   ceq%phase_varres(lokcs)%amount=zero
   ceq%phase_varres(lokcs)%amfu=zero
   ceq%phase_varres(lokcs)%netcharge=zero
   ceq%phase_varres(lokcs)%dgm=zero
   ceq%phase_varres(lokcs)%gval=zero
   ceq%phase_varres(lokcs)%dgval=zero
   ceq%phase_varres(lokcs)%d2gval=zero
! Mark there is no disordered phase_varres record
   ceq%phase_varres(lokcs)%disfra%varreslink=0
!   write(*,*)'parrecords: ',lokcs,nsl,nc
1000 continue
   return
 end subroutine create_parrecords

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_interaction
!\begin{verbatim}
 subroutine create_interaction(intrec,mint,lint,intperm,intlinks)
! creates a parameter interaction record 
! with permutations if intperm(1)>0
   implicit none
   type(gtp_interaction), pointer :: intrec
   integer, dimension(2,*) :: lint,intlinks
   integer, dimension(*) :: intperm
   integer mint
!\end{verbatim}
   integer permut,emperm,nz,nq,lqq,ii,ll
!
!   write(*,5)'create interaction:',mint,lint(1,mint),lint(2,mint),&
!        (intperm(ii),ii=1,6)
5  format(a,i5,2x,2i3,2x,6i3)
   allocate(intrec)
! note that the order of values in intperm here is not the same in 
! fccpermuts or bccpermuts ??  Intlinks is the same
   permut=intperm(1)
   if(permut.le.0) then
! This is a default for no permutations, store 1's
      permut=0
      allocate(intrec%noofip(2))
      intrec%noofip(1)=1
      intrec%noofip(2)=1
      allocate(intrec%sublattice(1))
      allocate(intrec%fraclink(1))
   elseif(mint.eq.1) then
! Intperm contains information as created by fccpermut or bccpermut
! intperm(1) and 2 are related to mint=1 (level 1 interaction),
! intperm(3) to mint=2
! The values are stored in noofip(1) and intperm(2..) in noofip(2..)
! For mint=1 intperm(1..2) are stored in noofipermt(1..2)
!   intperm(1) is the number of interaction permutations for each
!    endmember permutation.
!   intperm(2) are the number total number of permutations on level 1
!   The number of endmember permutations is thus intperm(2)/intperm(1)
!      write(*,17)'intrec: ',mint,intperm(1),intperm(2)
      permut=intperm(2)
      nz=intperm(2)
      allocate(intrec%noofip(2))
      intrec%noofip(1)=intperm(1)
      intrec%noofip(2)=intperm(2)
      allocate(intrec%sublattice(nz))
      allocate(intrec%fraclink(nz))
      nq=0
   elseif(mint.eq.2) then
! For mint=2 intperm(3) is stored in noofip(1) and intperm(4..) after that
!   if intperm(3)>1 then there are intperm(3) number of limits in
!   intperm(2..) for each lower order interaction.
! Example endmember A:A:A:A; no permutations
! 1st level intperm(1)=1, intperm(2)=4; permutations AX:A:A:A, A:AX:A:A etc
! 2nd level intperm(1)=4, inteprm(2..4)=(3, 2, 1, 0)
!   3 permutations for AX:A:A:A: AX:AX:A:A; AX:A:AX:A; AX:A:A:AX
!   2 permutations for A:AX:A:A: A:AX:AX:A A:AX:A:AX;
!   1 permutation  for A:A:AX:A: A:A:AX:AX;
!   0 permutations for A:A:A:AX: none
! If noofpermut>1 the index selected of noofip is by the permutation of 
! the lower order interaction
! the value in intpermut(4+intperm(3)) is total number of permutations
      lqq=intperm(4+intperm(3))
!      write(*,17)'intrec: ',mint,intperm(3),(intperm(3+ii),ii=1,intperm(3))
17    format(a,2i4,2x,10i4)
      permut=intperm(3)
      emperm=intperm(2)/intperm(1)
      allocate(intrec%noofip(permut+2))
      nz=0
      intrec%noofip(1)=intperm(3)
      do ii=1,permut
         intrec%noofip(1+ii)=intperm(3+ii)
         nz=nz+intperm(3+ii)
      enddo
!      write(*,19)'ci: ',nz,emperm,permut,(intrec%noofip(j),j=1,permut+2)
19    format(a,10i4)
! AX:AX:A:A; 1 endmember permutation, 4 1st level permutations; 6 2nd level
! emperm=1; intperm(3)=4, intparm(4..6)=(3,2,1,0), nz=1*6=6
! AX:AX:B:B; 6 endmember permutation, 6 1st level permutations; 6 2nd level
! emperm=6; nz=1; nz=1*6=6
! number of permutations is related to the previous level
!      nz=nz*emperm
      nz=lqq
!      write(*,*)'Level 2 permutations: ',nz
      allocate(intrec%sublattice(nz))
      allocate(intrec%fraclink(nz))
! Save at the end the total number of permutations stored
      intrec%noofip(permut+2)=nz
      nq=intperm(2)
!      write(*,19)'c2: ',nz,emperm,permut,(intrec%noofip(j),j=1,permut+2)
!      write(*,17)'level 2 permutations: ',nz,emperm,nq,lqq
   else
      write(*,*)'Create_interaction called with too many permutations'
      gx%bmperr=4329; goto 1000
   endif
   if(permut.eq.0) then
! this is again a default when there are no permutations
      intrec%sublattice(1)=lint(1,mint)
      intrec%fraclink(1)=lint(2,mint)
   else
! We can have cases like noofiperumt(1)=1; noofip(2)=4 or
! noofip(1)=4; noofip(2..5)=(4, 3, 2, 1)
! nq is 0 for first level, intperm(2) for second level
      do ll=1,nz
         intrec%sublattice(ll)=intlinks(1,nq+ll)
         intrec%fraclink(ll)=intlinks(2,nq+ll)
      enddo
!      write(*,99)'isp: ',mint,&
!           (intrec%sublattice(ll),intrec%fraclink(ll),ll=1,nz)
99    format(a,i2,8(2x,2i3))
   endif
   nullify(intrec%propointer)
   nullify(intrec%nextlink)
   nullify(intrec%highlink)
! nullify Kohler-Toop link
!   write(*,*)'3G nullifying tooprec pointer'
   nullify(intrec%tooprec)
   intrec%status=0
   noofint=noofint+1
   intrec%antalint=noofint
1000 continue
   return
 end subroutine create_interaction

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_endmember
!\begin{verbatim}
 subroutine create_endmember(lokph,newem,noperm,nsl,endm,elinks)
! create endmember record with nsl sublattices with endm as constituents
! noperm is number of permutations
! endm is the basic endmember (if there are permutations)
! elinks are the links to constituents for all permutations
   implicit none
   integer endm(*)
   integer lokph,noperm,nsl
   type(gtp_endmember), pointer :: newem
   integer, dimension(nsl,noperm) :: elinks
!\end{verbatim}
   integer is,ndemr,noemr,nn
   allocate(newem)
   nullify(newem%nextem)
   allocate(newem%fraclinks(nsl,noperm))
   if(noperm.eq.1) then
      do is=1,nsl 
         newem%fraclinks(is,1)=endm(is)
      enddo
   else
!      write(*,*)'3G permutations: ',noperm,nsl
!      write(*,7)((elinks(is,nn),is=1,4),nn=1,noperm)
7     format('3G ce1: ',4(4i3,2x))
      newem%fraclinks=elinks
   endif
! zero or set values
   newem%noofpermut=noperm
   newem%phaselink=lokph
   noofem=noofem+1
   newem%antalem=noofem
   nullify(newem%propointer)
   nullify(newem%intpointer)
! indicate that oendmemarr and denmemarr must be renewed ???
   noemr=0
   ndemr=0
1000 continue
   return
 end subroutine create_endmember

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_proprec
!\begin{verbatim}
 subroutine create_proprec(proprec,proptype,degree,lfun,refx)
! reservs a property record from free list and insert data
   implicit none
   TYPE(gtp_property), pointer :: proprec
   integer proptype,degree,lfun
   character refx*(*)
!\end{verbatim} %+
   integer j,iref
   character notext*32
   if(degree.lt.0 .or. degree.gt.9) then
      gx%bmperr=4063; goto 1000
   endif
   allocate(proprec)
! enter data in reserved record
   allocate(proprec%degreelink(0:degree))
   nullify(proprec%nextpr)
!   if(proptype.ge.100) write(*,*)'property type: ',proptype
   proprec%proptype=proptype
! also save %modelparamid for unformatted files ...
   if(proptype.gt.100) then
! this is a property with a constituent suffix like MQ&FE
      proprec%modelparamid=propid(proptype/100)%symbol
!      write(*,*)'3G proptype ',propid(proptype/100)%symbol,proptype
   else
      proprec%modelparamid=propid(proptype)%symbol
!      write(*,*)'3G proptype ',propid(proptype)%symbol,proptype
   endif
   proprec%degree=degree
   do j=0,degree
      proprec%degreelink(j)=0
   enddo
   proprec%degreelink(degree)=lfun
   proprec%reference=adjustl(refx)
! create reference record if new, can be amended later
   call capson(refx)
   notext='*** Not set by database or user '
!------counter
   noofprop=noofprop+1
   proprec%antalprop=noofprop
!   write(*,11)refx,notext
!11 format('create proprec: ',a,a)
   call tdbrefs(refx,notext,0,iref)
   proprec%extra=0
1000 continue
   return
 end subroutine create_proprec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine extend_proprec
!\begin{verbatim}
 subroutine extend_proprec(current,degree,lfun)
! extends a property record  and insert new data
   implicit none
   integer degree,lfun
   type(gtp_property), pointer :: current
!\end{verbatim}
   integer oldeg,j
   integer :: savedegs(0:9)
! save degreelinks ... maybe not necessary ....
   oldeg=current%degree
!    write(*,*)'extend_proprec 1: ',current,degree,lfun,oldeg
   do j=0,9
      savedegs(j)=0
   enddo
   do j=0,oldeg
      savedegs(j)=current%degreelink(j)
   enddo
! important to get it correct here
   deallocate(current%degreelink)
   allocate(current%degreelink(0:degree))
   current%degree=degree
   do j=0,current%degree
      current%degreelink(j)=0
   enddo
   do j=0,oldeg
      current%degreelink(j)=savedegs(j)
   enddo
   current%degreelink(degree)=lfun
1000 continue
   return
 end subroutine extend_proprec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine add_fraction_set
!\begin{verbatim}
 subroutine add_fraction_set(iph,id,ndl,totdis)
! add a new set of fractions to a phase, usually to describe a disordered state
! like the "partitioning" in old TC
!
! BEWARE this is only done for firsteq, illegal when having more equilibria
!
! id is a letter used as suffix to identify the parameters of this set
! ndl is the last original sublattice included in the (first) disordered set
! ndl can be 1 meaning sublattice 2..nsl are disordered, or nsl meaning all are
!     disordered
! totdis=0 if phase never disorder totally (like sigma)
!
! For a phase like (Al,Fe,Ni)3(Al,Fe,Ni)1(C,Va)4 to add (Al,Fe,Ni)4(C,Va)4
! icon=1 2 3 1 2 3 4 5 with ndl=2
! For a phase like (Fe,Ni)10(Cr,Mo)4(Cr,Fe,Mo,Ni)16 then
! icon=2 4 1 3 1 2 3 4 with ndl=3
! This subroutine will create the necessary data to calculate the
! disordered fraction set from the site fractions.
!
! IMPORTANT (done): for each composition set this must be repeated
! if new composition sets are created it must be repeated for these
!
! IMPORTANT (not done): order the constituents alphabetically in each disorderd
! sublattice otherwise it will not be possible to enter parameters correctly
!
   implicit none
   integer iph,ndl,totdis
   character id*1
!\end{verbatim}
! ceq probably not needed as firsteq is declared as pointer
!   TYPE(gtp_equilibrium_data), target :: ceq
   TYPE(gtp_fraction_set), target :: fsdata
! jsp(i) contains species locations of disordered constituent i
! jy2x(i) is the disordered fraction to which site fraction i should be added
! y2x(i) is the site ration factor for multiplying sitefraction i when added
! ispord and ispold are needed to sort disordered constituents
   integer jsp(maxconst,2),jy2x(maxconst),iva(maxconst)
   integer ispord(maxconst),ispold(maxconst),nrj3(2),nrj4(2)
   integer lokph,lokcs,nsl,ii,nrj1,nrj2,nlat,lokx,l2
   integer ll,kk,jall,nnn,mmm,ioff,koff,jl,j1,j2,ix,is,jj,k,ijcs,nydis,nyttcs
   double precision sum,div
!
   if(.not.allowenter(2)) then
      gx%bmperr=4125
      goto 1000
   endif
! this subroutine can only be called when there is only one equilibrium
   lokph=phases(iph)
! phase must not have any suspended constituents nor any composition sets
   if(phlista(lokph)%noofcs.gt.1) then
      gx%bmperr=4029; goto 1000
   else
      lokcs=phlista(lokph)%linktocs(1)
      if(btest(firsteq%phase_varres(lokcs)%status2,CSCONSUS)) then
         gx%bmperr=4030; goto 1000
      endif
   endif
   nsl=phlista(lokph)%noofsubl
   if(ndl.le.1 .or. ndl.gt.nsl) then
! ndl must be larger than 2 and lesser or equal to nsl
      gx%bmperr=4076; goto 1000
   endif
! location of first composition set, there may be more
   if(btest(phlista(lokph)%status1,phmfs)) then
! disordered fractions already set
      gx%bmperr=4077; goto 1000
   endif
!   write(*,*)'3G in add_fr: ',iph,id,ndl,totdis
! we must organise a constituent list for the disordered fractions by
! scanning the constituents in the current phlista(lokph)%constitlist
! we must also contruct the way site fractions should be added
   ii=0
   nrj1=1
   nrj2=0
   nlat=0
   lokx=0
   l2=1
   iva=0
   subloop: do ll=1,nsl
      constloop: do kk=1,phlista(lokph)%nooffr(ll)
         ii=ii+1
         if(nrj2.lt.nrj1) then
            nrj2=nrj2+1
            lokx=lokx+1
            jy2x(ii)=lokx
            jsp(nrj2,l2)=phlista(lokph)%constitlist(ii)
!            write(*,46)'new 1: ',nrj2,l2,ii,nlat,jsp(nrj2,l2),jy2x(ii)
         else
            do jall=nrj1,nrj2
               if(phlista(lokph)%constitlist(ii).eq.jsp(jall,l2)) then
! this constituent already found in another sublattice to be merged
!                  write(*,*)'same: ',jall,nlat,jall+nlat,ii,jy2x(jall+nlat)
                  jy2x(ii)=jy2x(jall+nlat);  goto 50
               endif
            enddo
! new constituent
            nrj2=nrj2+1
            lokx=lokx+1
            jy2x(ii)=lokx
            jsp(nrj2,l2)=phlista(lokph)%constitlist(ii)
!            write(*,46)'new 2: ',nrj2,l2,ii,nlat,jsp(nrj2,l2),jy2x(ii)
46          format(a,10i3)
! if vacancy set that bit in iva
            if(btest(firsteq%phase_varres(lokcs)%constat(ii),conva)) then
               iva(nrj2)=ibset(iva(nrj2),conva)
            endif
!             write(*,*)'addfs 7B: ',ll,ii,nrj2
50           continue
         endif
      enddo constloop
      if(ll.eq.ndl) then
! next sublattices (if any) will be summed to second disordered sublattice
         nrj3(1)=nrj2
         nrj3(2)=0
! bug??
         nlat=ii
         nrj1=1
         nrj2=0
! nrj4 is the number of constituents in ordered phase thst is summed
! to first disordered sublattice.  Needed below to rearrange jy2x
         nrj4(1)=ii
         nrj4(2)=0
         if(ndl.lt.nsl) l2=2
!          write(*,*)'addfs 7C: ',ll,ndl,nrj1,nrj2,nrj3
      elseif(ll.eq.nsl) then
! this may never be executed if ndl=nsl but we set nrj3(2)=0 above
         nrj3(2)=nrj2
         nrj4(2)=ii-nrj4(1)
      endif
   enddo subloop
!   write(*,53)'add_fraction_set 2: ',(jy2x(i),i=1,ii)
53 format(a,20i3)
! added fsites to handle the case when reading sigma etc from a TDB file
! as the TDB file format assumes 1 site.  Default is 1.0, changed externally
   fsdata%fsites=one
!   write(*,*)'3G Set fsites: ',fsdata%fsites,ndl,totdis,nnn
!
!   write(*,53)'add_fraction_set 3: ',nrj1,nrj2,nrj3,nrj4
   fsdata%latd=ndl
   fsdata%tnoofyfr=phlista(lokph)%tnooffr
   fsdata%varreslink=lokcs
! totdis=1 means disordered fcc, bcc, ncp. totdis=0 means sigma
   fsdata%totdis=totdis
   fsdata%id=id
! one or 2 disordered sublattices
   nnn=1
   if(ndl.lt.nsl) nnn=2
! try to allow more than one interstitial sublattice ... NO
   if(nsl-ndl.gt.1) then
      write(*,*)'3G *** Error max one sublattices outside the disordered set'
      gx%bmperr=4399
      goto 1000
   endif
!   write(*,'(a,10i5)')'3G disordered sublattices:',nsl,ndl,nnn,fsdata%latd
   allocate(fsdata%dsites(nnn))
   fsdata%ndd=nnn
   allocate(fsdata%nooffr(nnn))
   fsdata%nooffr(1)=nrj3(1)
   if(nnn.eq.2) fsdata%nooffr(2)=nrj3(2)
! nrj3(1) are the number of constituents on first sublattice, nrj3(2) on 2nd
   mmm=nrj3(1)+nrj3(2)
   fsdata%tnoofxfr=mmm
   allocate(fsdata%splink(mmm))
   allocate(fsdata%y2x(phlista(lokph)%tnooffr))
   allocate(fsdata%dxidyj(phlista(lokph)%tnooffr))
!    write(*,*)'add_fs dxidyj: ',phlista(lokph)%tnooffr
! the constituents in jsp(i..n,subl) must be ordered alphabetically!!!
! get the species number in alphadetical order
   ioff=0
   koff=0
   do l2=1,nnn
      do jl=1,nrj3(l2)
!         write(*,*)'l2 loop: ',jsp(i,l2)
         ispord(jl)=splista(jsp(jl,l2))%alphaindex
      enddo
!       write(*,47)1,(ispord(i),i=1,nrj3(l2))
47     format('add_fs ',i1,': ',20i3)
!                 species, noofsp, origonal order
      call sortin(ispord,nrj3(l2),ispold)
      if(buperr.ne.0) then
         gx%bmperr=buperr; goto 1000
      endif
! when rearranging jsp(1..n,l2) we must also rearrange y2x
! for 2nd sublattice add nrj3(1) to ispold
      if(l2.eq.2) then
         ioff=nrj4(1)
         koff=nrj3(1)
      endif
!       write(*,47)2,(jy2x(ioff+i),i=1,nrj4(l2))
! this must be possible to do smarter .....
      do j2=1,nrj4(l2)
         do j1=1,nrj3(l2)
            if(jy2x(ioff+j2).eq.ispold(j1)+koff) then
               jy2x(ioff+j2)=j1+koff; goto 77
            endif
         enddo
77        continue
      enddo
      do j1=1,nrj3(l2)
         ispord(j1)=jsp(ispold(j1),l2)
      enddo
      do j1=1,nrj3(l2)
         jsp(j1,l2)=ispord(j1)
      enddo
!       write(*,47)5,(jsp(i,l2),i=1,nrj3(l2))
   enddo
   fsdata%splink=0
!
   do jl=1,phlista(lokph)%tnooffr
      fsdata%y2x(jl)=jy2x(jl)
   enddo
   ix=0
   do l2=1,nnn
      do jl=1,nrj3(l2)
         ix=ix+1
         fsdata%splink(ix)=jsp(jl,l2)
      enddo
   enddo
!    write(*,*)'addfs splink: ',fsdata%splink
!
   is=0
   sum=zero
   do ll=1,ndl
!      sum=sum+phlista(lokph)%sites(ll)
      sum=sum+firsteq%phase_varres(lokcs)%sites(ll)
   enddo
   fsdata%dsites(1)=sum
   if(ndl.lt.nsl) then
      sum=zero
      do ll=ndl+1,nsl
!         sum=sum+phlista(lokph)%sites(ll)
         sum=sum+firsteq%phase_varres(lokcs)%sites(ll)
      enddo
      fsdata%dsites(2)=sum
   endif
!
   jj=0
   sum=fsdata%dsites(1)
!   write(*,*)'3G sum: ',ndl,sum,fsdata%dsites
   do ll=1,nsl
      if(ll.gt.ndl) sum=fsdata%dsites(2)
!      div=phlista(lokph)%sites(ll)/sum
      div=firsteq%phase_varres(lokcs)%sites(ll)/sum
!       write(*,78)'add_fs 5A ',div,phlista(lokph)%sites(ll),sum
!78     format(a,6F10.7)
      do k=1,phlista(lokph)%nooffr(ll)
         jj=jj+1
         fsdata%dxidyj(jj)=div
      enddo
   enddo
!   write(*,99)'add_fs 5B ',fsdata%dxidyj
99 format(a,6(F10.7))
   firsteq%phase_varres(lokcs)%disfra=fsdata
   firsteq%phase_varres(lokcs)%status2=&
        ibset(firsteq%phase_varres(lokcs)%status2,CSDLNK)
! we have to reserve a phase_varres record for calculations
!  ... det gäller att hålla tungan rätt i mun ...
!   nprop=10
!   call create_parrecords(nyttcs,nnn,mmm,nprop,iva,firsteq)
   call create_parrecords(lokph,nyttcs,nnn,mmm,maxcalcprop,iva,firsteq)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3G created disordered phase_varres: ',csfree,highcs,nyttcs
   fsdata%varreslink=nyttcs
! note ceq is firsteq but declared target
!   write(*,*)'3G disordered fraction set',nyttcs
!*?   fsdata%phdapointer=>ceq%phase_varres(nyttcs)
   firsteq%phase_varres(nyttcs)%phlink=lokph
   firsteq%phase_varres(nyttcs)%prefix=' '
   firsteq%phase_varres(nyttcs)%suffix=' '
   do ll=1,nnn
      firsteq%phase_varres(nyttcs)%sites(ll)=fsdata%dsites(ll)
   enddo
   firsteq%phase_varres(nyttcs)%status2=0
   firsteq%phase_varres(nyttcs)%status2=&
        ibset(firsteq%phase_varres(nyttcs)%status2,CSDFS)
! finally copy fsdata to the link in lokcs
   call copy_fracset_record(lokcs,fsdata,firsteq)
   if(gx%bmperr.ne.0) goto 1000
! if there are several composition sets create fracset records for each
200 continue
!   if(firsteq%phase_varres(lokcs)%next.gt.0) then
!      lokcs=firsteq%phase_varres(lokcs)%next
   do ijcs=2,phlista(lokph)%noofcs
      lokcs=phlista(lokph)%linktocs(ijcs)
! one must also create parrecords for these !!!
!      call create_parrecords(nydis,nnn,mmm,nprop,iva,firsteq)
      call create_parrecords(lokph,nydis,nnn,mmm,maxcalcprop,iva,firsteq)
      if(gx%bmperr.ne.0) goto 1000
      fsdata%varreslink=nydis
! set pointer also
!*?      fsdata%phdapointer=firsteq%phase_varres(nydis)
      firsteq%phase_varres(nydis)%phlink=lokph
      firsteq%phase_varres(nydis)%prefix=' '
      firsteq%phase_varres(nydis)%suffix=' '
      do ll=1,nnn
         firsteq%phase_varres(nydis)%sites(ll)=fsdata%dsites(ll)
      enddo
      firsteq%phase_varres(nydis)%status2=0
      firsteq%phase_varres(nydis)%status2=&
           ibset(firsteq%phase_varres(nyttcs)%status2,CSDFS)
! This does not create a new record
!       firsteq%phase_varres(lokcs)%disfra=fsdata
! but this seems to work
      call copy_fracset_record(lokcs,fsdata,firsteq)
      if(gx%bmperr.ne.0) goto 1000
      firsteq%phase_varres(lokcs)%status2=&
           ibset(firsteq%phase_varres(lokcs)%status2,CSDLNK)
      goto 200
   enddo
! set status bit for multiple/disordered fraction sets and no of fraction sets
   phlista(lokph)%status1=ibset(phlista(lokph)%status1,PHMFS)
   phlista(lokph)%nooffs=2
1000 continue
!   write(*,*)'3G exit add_fraction_set: ',fsdata%fsites,nnn
! NOTE fsdata&fsites updated in calling routine.  A bit strange but ...
   return
! nydis
 end subroutine add_fraction_set  ! no ceq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine copy_fracset_record
!\begin{verbatim}
 subroutine copy_fracset_record(lokcs,disrec,ceq)
! attempt to create a new disordered record  ??? this can probably be done
! with just one statement .. but as it works I am not changing right now
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_fraction_set) :: disrec
   integer lokcs
!\end{verbatim}
   TYPE(gtp_fraction_set) :: discopy
! the hard way ??
   discopy%fsites=disrec%fsites
   discopy%latd=disrec%latd
   discopy%ndd=disrec%ndd
   discopy%tnoofxfr=disrec%tnoofxfr
   discopy%tnoofyfr=disrec%tnoofyfr
   discopy%varreslink=disrec%varreslink
!*?   discopy%phdapointer=>disrec%phdapointer
   discopy%totdis=disrec%totdis
   discopy%id=disrec%id
   allocate(discopy%dsites(disrec%ndd))
   allocate(discopy%nooffr(disrec%ndd))
   allocate(discopy%splink(disrec%tnoofxfr))
   allocate(discopy%y2x(disrec%tnoofyfr))
   allocate(discopy%dxidyj(disrec%tnoofyfr))
!
   discopy%dsites=disrec%dsites
   discopy%nooffr=disrec%nooffr
   discopy%splink=disrec%splink
   discopy%y2x=disrec%y2x
   discopy%dxidyj=disrec%dxidyj
!
!    write(*,*)'copyfs 1: ',lokcs,discopy%varreslink,disrec%varreslink
   ceq%phase_varres(lokcs)%disfra=discopy
!    write(*,*)'copyfs 2: ',phase_varres(lokcs)%disfra%varreslink
1000 continue
   return
 end subroutine copy_fracset_record

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine suspend_species_implicitly
!\begin{verbatim}
 subroutine suspend_species_implicitly(ceq)
! loop through all entered species and suspend those with an element suspended
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer loksp,iel,lokel
   sploop: do loksp=1,noofsp
      if(.not.btest(splista(loksp)%status,spsus)) then
         elloop: do iel=1,splista(loksp)%noofel
            lokel=splista(loksp)%ellinks(iel)
            if(btest(ellista(lokel)%status,elsus)) then
! an element is suspended, suspend this species implicitly
               splista(loksp)%status=ibset(splista(loksp)%status,spsus)
               splista(loksp)%status=ibset(splista(loksp)%status,spimsus)
               goto 200
            endif
         enddo elloop
      endif
200    continue
   enddo sploop
1000 continue
   return
 end subroutine suspend_species_implicitly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine suspend_phases_implicitly
!\begin{verbatim}
 subroutine suspend_phases_implicitly(ceq)
! loop through all entered phases and suspend constituents and
! SUSPEND phases with all constituents in a sublattice suspended
!   dimension lokcs(9)
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,lokcs,ncc,kk,kkl,nek,icon,ll,loksp,jl
!
! BEWARE not adopted fro parallel processing
!
   phloop: do lokph=1,noofph
      if(.not.btest(phlista(lokph)%status1,phhid)) then
! locate all composition sets and store indices in lokcs
         ncc=phlista(lokph)%noofcs
         kk=0
         sublloop: do ll=1,phlista(lokph)%noofsubl
            kkl=kk
            nek=0
            constloop: do icon=1,phlista(lokph)%nooffr(ll)
               kk=kk+1
               loksp=phlista(lokph)%constitlist(kk)
               if(btest(splista(loksp)%status,spsus)) then
! a constituent is suspended, mark this also in constat for all comp.sets
                  compsets: do jl=1,ncc
                     lokcs=phlista(lokph)%linktocs(jl)
                     ceq%phase_varres(lokcs)%constat(kk)=&
                          ibset(ceq%phase_varres(lokcs)%constat(kk),consus)
                     ceq%phase_varres(lokcs)%constat(kk)=&
                         ibset(ceq%phase_varres(lokcs)%constat(kk),conimsus)
! mark that some constituents are suspended in this composition set
                     ceq%phase_varres(lokcs)%status2=&
                          ibset(ceq%phase_varres(lokcs)%status2,CSCONSUS)
                  enddo compsets
                  goto 200
               else
                  nek=nek+1
               endif
            enddo constloop
            if(nek.eq.0) then
! this sublattice has all constituents suspended, hide/suspend the phase
               phlista(lokph)%status1=ibset(phlista(lokph)%status1,phhid)
               phlista(lokph)%status1=ibset(phlista(lokph)%status1,phimhid)
! also set amount to zero ??
               compsets2: do jl=1,ncc
                  lokcs=phlista(lokph)%linktocs(jl)
!                  ceq%phase_varres(lokcs)%amount=zero
                  ceq%phase_varres(lokcs)%amfu=zero
                  ceq%phase_varres(lokcs)%netcharge=zero
               enddo compsets2
            endif
            goto 300
200         continue
            kk=kkl+phlista(lokph)%nooffr(ll)
            kkl=kk-1
         enddo sublloop
300      continue
      endif
   enddo phloop
1000 continue
   return
 end subroutine suspend_phases_implicitly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine restore_species_implicitly_suspended
!\begin{verbatim}
 subroutine restore_species_implicitly_suspended
! loop through all implicitly suspended species and restore those with
! all elements enteded
   implicit none
!\end{verbatim} %+
   integer loksp,lokel
   sploop: do loksp=1,noofsp
      if(btest(splista(loksp)%status,spimsus)) then
         elloop: do lokel=1,splista(loksp)%noofel
! an element is suspended, keep species suspended
            if(btest(ellista(lokel)%status,elsus)) goto 200
         enddo elloop
! all elements entered, restore species as entered
         splista(loksp)%status=ibclr(splista(loksp)%status,spsus)
         splista(loksp)%status=ibclr(splista(loksp)%status,spimsus)
      endif
200    continue
   enddo sploop
1000 continue
   return
 end subroutine restore_species_implicitly_suspended

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine restore_phases_implicitly_suspended
!\begin{verbatim}
 subroutine restore_phases_implicitly_suspended
! loop through all implicitly suspended phases and restore those with
! at least one constituent entered in each sublattice
   implicit none
!\end{verbatim}
   integer lokph,ll,kk,kkl,icon,loksp
   phloop: do lokph=1,noofph
      if(btest(phlista(lokph)%status1,phimhid)) then
         kk=0
         sublloop: do ll=1,phlista(lokph)%noofsubl
            kkl=kk
            constloop: do icon=1,phlista(lokph)%nooffr(ll)
               kk=kk+1
               loksp=phlista(lokph)%constitlist(kk)
               if(.not.btest(splista(loksp)%status,spsus)) goto 200
            enddo constloop
! all constituents in this sublattice are suspended, keep the phase hidden
            goto 300
200          continue
            kk=kkl+phlista(lokph)%nooffr(ll)
            kkl=kk-1
         enddo sublloop
! all sublattices have at least one constituent entered, restore it
         phlista(lokph)%status1=ibclr(phlista(lokph)%status1,phhid)
         phlista(lokph)%status1=ibclr(phlista(lokph)%status1,phimhid)
300       continue
      endif
   enddo phloop
1000 continue
   return
 end subroutine restore_phases_implicitly_suspended

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine add_to_reference_phase
!\begin{verbatim}
 subroutine add_to_reference_phase(loksp)
! add this element to the reference phase
! loksp: species index of new element
   implicit none
   integer loksp
!\end{verbatim}
! one must extend all arrays in phlista, phase_varres and phase_varres
   integer lokph,noc,i,nprop,mc2,lokcs
   integer, dimension(maxel) :: isave
   lokph=0
   lokcs=phlista(lokph)%linktocs(1)
! constitlist
   noc=phlista(lokph)%tnooffr
   do i=1,noc
      isave(i)=phlista(lokph)%constitlist(i)
   enddo
   deallocate(phlista(lokph)%constitlist)
   noc=noc+1
   allocate(phlista(lokph)%constitlist(noc))
   isave(noc)=loksp
   do i=1,noc
      phlista(lokph)%constitlist(i)=isave(i)
   enddo
   phlista(lokph)%tnooffr=noc
   phlista(lokph)%nooffr(1)=noc
! phase_varres, no data need saving
!  write(*,*)'Deallocate constat 5: ',size(firsteq%phase_varres(lokcs)%constat)
   deallocate(firsteq%phase_varres(lokcs)%constat)
   deallocate(firsteq%phase_varres(lokcs)%yfr)
   deallocate(firsteq%phase_varres(lokcs)%mmyfr)
!   write(*,*)'Allocate constat 5: ',noc
   allocate(firsteq%phase_varres(lokcs)%constat(noc))
   firsteq%phase_varres(lokcs)%constat(noc)=0
   allocate(firsteq%phase_varres(lokcs)%yfr(noc))
   allocate(firsteq%phase_varres(lokcs)%mmyfr(noc))
   firsteq%phase_varres(lokcs)%yfr=one
   firsteq%phase_varres(lokcs)%mmyfr=zero
   nprop=firsteq%phase_varres(lokcs)%nprop
   deallocate(firsteq%phase_varres(lokcs)%dgval)
   deallocate(firsteq%phase_varres(lokcs)%d2gval)
   allocate(firsteq%phase_varres(lokcs)%dgval(3,noc,nprop))
   mc2=noc*(noc+1)/2
   allocate(firsteq%phase_varres(lokcs)%d2gval(mc2,nprop))
! ready!!
1000 continue
   return
 end subroutine add_to_reference_phase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
