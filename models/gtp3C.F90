!
! gtp3C included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     7. List data
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_all_elements(unit)
! lists elements
   implicit none
   integer unit
!\end{verbatim} %+
   integer jl,ipos
   character line*80
   line=' '
   write(unit,10)noofel
10  format(/'List of ',i2,' elements'/ &
        ' No Sym Name',10X,'Reference state',12X,&
        'Mass  H298-H0   S298    Status')
   loop1: do jl=-1,noofel
      ipos=1
      call list_element_data(line,ipos,elements(jl))
      if(gx%bmperr.ne.0) goto 1000
      write(unit,100)jl,line(1:ipos)
   enddo loop1
100 format(i3,2x,A)
1000 continue
   return
 end subroutine list_all_elements

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine list_all_elements2(unit)
! lists elements
   implicit none
   integer unit
!\end{verbatim}
   integer jl
   character line*80
   line=' '
   loop1: do jl=-1,noofel
      write(unit,100) ellista(jl)%symbol,ellista(jl)%ref_state,&
           ellista(jl)%mass,ellista(jl)%h298_h0,ellista(jl)%s298
   enddo loop1
100 format('ELEMENT ',A,'  ',A,3(1pe12.4),' !')
1000 continue
   return
 END subroutine list_all_elements2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_all_components(unit,ceq)
! lists the components for an equilibrium
   implicit none
   integer unit
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer jl,loksp
   character symbol*24
   double precision moles,masspercent,chempot
   moles=zero
   masspercent=zero
   chempot=zero
   write(unit,10)
10  format('List of components'/ &
        'No Symbol',19X,'Moles',6x,'Mass %',5x,'Chem pot',3x,'Ref. state')
   loop1: do jl=1,noofel
      loksp=ceq%complist(jl)%splink
      symbol=splista(loksp)%symbol
      write(unit,100)jl,symbol,moles,masspercent,chempot,&
           ceq%complist(jl)%refstate
   enddo loop1
100 format(i2,1x,A,3(1PE11.3),1X,A)
1000 continue
   return
 end subroutine list_all_components

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_element_data(text,ipos,elno)
   implicit none
   character text*(*)
   integer ipos,elno
!\end{verbatim}
   if(elno.lt.-1 .or. elno.gt.noofel) then
      gx%bmperr=4042
      goto 1000
   endif
   if(ipos.lt.1 .or. ipos.ge.len(text)) then
      gx%bmperr=4043
      goto 1000
   endif
   text(ipos:ipos+2)=ellista(elno)%symbol
   text(ipos+3:ipos+16)=ellista(elno)%name
   text(ipos+17:ipos+40)=ellista(elno)%ref_state
   write(text(ipos+41:ipos+73),100)ellista(elno)%mass,&
        ellista(elno)%h298_h0,ellista(elno)%s298,ellista(elno)%status
100 format(1x,f7.3,1x,f7.2,1x,f7.3,1x,z8)
   ipos=len_trim(text)
!   write(*,*)'x:',text(1:79)
1000 continue
   return
 END subroutine list_element_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_species_data(text,ipos,spno)
   implicit none
   character text*(*)
   integer ipos,spno
!\end{verbatim} %+
   character dummy*48
   integer jpos
   if(spno.lt.1 .or. spno.gt.noofsp) then
!       write(*,*)'in list_species_data'
      gx%bmperr=4051
      goto 1000
   endif
   if(ipos.lt.1 .or. ipos.ge.len(text)) then
      gx%bmperr=4043
      goto 1000
   endif
   text(ipos:ipos+24)=splista(spno)%symbol
   text(ipos+25:ipos+25)=' '
   dummy=' '
   call encode_stoik(dummy,jpos,spno)
   text(ipos+26:ipos+48)=dummy(1:min(23,jpos))
   if(jpos.gt.23) text(ipos+46:ipos+48)='<.>'
   text(ipos+49:ipos+49)=' '
   write(text(ipos+50:ipos+59),100)splista(spno)%mass
   write(text(ipos+60:ipos+65),105)splista(spno)%charge
100 format(F10.3)
105 format(F6.1)
   text(ipos+66:)=' '
!    write(*,120)splista(spno)%status
   write(text(ipos+66:ipos+73),120)splista(spno)%status
120 format(Z8)
   ipos=ipos+73
1000 continue
   return
 END subroutine list_species_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine list_species_data2(text,ipos,loksp)
! loksp is species record ...
   implicit none
   character text*(*)
   integer ipos,loksp
!\end{verbatim}
   character dummy*24
   integer jpos
   if(loksp.lt.1 .or. loksp.gt.noofsp) then
!       write(*,*)'in list_species_data'
      gx%bmperr=4051
      goto 1000
   endif
   if(ipos.lt.1 .or. ipos.ge.len(text)) then
      gx%bmperr=4043
      goto 1000
   endif
   text(ipos:ipos+24)=splista(loksp)%symbol
   text(ipos+25:ipos+25)=' '
   dummy=' '
   call encode_stoik(dummy,jpos,loksp)
   text(ipos+26:ipos+48)=dummy(1:jpos)
!   text(ipos+49:ipos+49)=' '
!   write(text(ipos+50:ipos+59),100)splista(loksp)%mass
!   write(text(ipos+60:ipos+65),105)splista(loksp)%charge
100 format(F10.3)
105 format(F6.1)
!   text(ipos+66:)=' '
!    write(*,120)splista(loksp)%status
!   write(text(ipos+66:ipos+73),120)splista(loksp)%status
120 format(Z8)
!   ipos=ipos+73
1000 continue
   return
 END subroutine list_species_data2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_all_species(unit)
   implicit none
   integer unit
!\end{verbatim}
   integer jl,ipos
   character line*80
   write(unit,10)noofsp
10  format(/'List of ',i3,' species'/ &
        '  No Symbol',20X,'Stoichiometry',12X,'Mass      Charge Status')
   loop1: do jl=1,noofsp
      ipos=1
      call list_species_data(line,ipos,species(jl))
      if(gx%bmperr.ne.0) goto 1000
      write(unit,100)jl,line(1:ipos)
   enddo loop1
100 format(i4,1x,A)
1000 continue
   return
 END subroutine list_all_species

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_sorted_phases(unit,ceq)
! short list with one line for each phase
! suspended phases merged into one line
! stable first, then entered ordered in driving force order, then dormant
! also in driving force order.  Only 10 of each, the others lumped together
   implicit none
   integer unit
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer jl,jk,ics,lokph,lokcs,kp,ndorm,nsusp,nent,nstab,iph,jph
   character line*80,phname*24,trailer*28,chs*1,csname*36,susph*4096,ch1*1
   integer, dimension(:), allocatable :: entph,dorph
   TYPE(gtp_phase_varres), pointer :: csrec
   double precision am1,am2
!
   allocate(entph(nooftuples))
   allocate(dorph(nooftuples))
   nstab=0; nent=0; ndorm=0; nsusp=1
   susph=' '
   phloop: do jk=1,noofph
      lokph=phases(jk)
      csloop: do ics=1,phlista(lokph)%noofcs
!         write(*,17)'3C sort1: ',nent,(entph(iph),iph=1,nent)
17       format(a,i3,2x,16(i4))
         lokcs=phlista(lokph)%linktocs(ics)
         csrec=>ceq%phase_varres(lokcs)
         if(csrec%phstate.ge.PHENTSTAB) then
            if(nent.eq.0) then
               nent=1;
               entph(nent)=lokcs
!               write(*,*)'3C first phase stable: ',nent,nent,lokcs
            else
! FIX and STABLE phases first in order of amount
               do iph=1,nent
                  am1=csrec%amfu*csrec%abnorm(1)
                  am2=ceq%phase_varres(entph(iph))%amfu*&
                       ceq%phase_varres(entph(iph))%abnorm(1)
                  if(am1.lt.am2) cycle
!                  if(csrec%amfu.lt.ceq%phase_varres(entph(iph))%amfu) cycle
! this is the place for this phase, shift later down
                  do jph=nent,iph,-1
                     entph(jph+1)=entph(jph)
                  enddo
                  exit
               enddo
! according to new fortran standard loop variable at exit is high limit+1
!               write(*,18)'3C inserted stable phase ',iph,lokcs,csrec%amfu
18             format(a,2i4,1pe12.4)
               entph(iph)=lokcs
               nent=nent+1
!               write(*,*)'3C stable phase: ',nent,iph,lokcs
            endif
         elseif(csrec%phstate.eq.PHENTERED .or. &
              csrec%phstate.eq.PHENTUNST) then
            if(nent.eq.0) then
               nent=1
               entph(nent)=lokcs
!               write(*,69)'3C first phase unstable: ',nent,nent,lokcs,csrec%dgm
            else
! ENTERED, not stable, sort after all stable phase and with smallest DGM first
               do iph=1,nent
! bypass all stable phases
                  if(ceq%phase_varres(entph(iph))%amfu.gt.zero) cycle
                  if(csrec%dgm.lt.ceq%phase_varres(entph(iph))%dgm) cycle
! this is the place for this phase, shift later phases down
                  do jph=nent,iph,-1
                     entph(jph+1)=entph(jph)
                  enddo
                  exit
               enddo
! according to new fortran standard loop variable at exit is high limit+1
!               write(*,18)'3C inserted ustable phase ',iph,lokcs,csrec%dgm
               entph(iph)=lokcs
               nent=nent+1
!               write(*,69)'3C unstable phase: ',iph,nent,lokcs,csrec%dgm
69             format(a,3i4,1pe12.4)
            endif
         elseif(csrec%phstate.eq.PHDORM) then
            if(ndorm.eq.0) then
               ndorm=ndorm+1
               dorph(ndorm)=lokcs
!               write(*,*)'3C first dormant phase: ',ndorm,ndorm,lokcs
            else
! DORMANT sort after with smallest (least nagative) DGM first
               do iph=1,ndorm
                  if(csrec%dgm.lt.ceq%phase_varres(dorph(iph))%dgm) cycle
! this is the place for this phase, shift later down
                  do jph=ndorm,iph,-1
                     dorph(jph+1)=dorph(jph)
                  enddo
                  exit
               enddo
! according to new fortran standard loop variable at exit is high limit+1
               dorph(iph)=lokcs
               ndorm=ndorm+1
!               write(*,*)'3C dormant phase: ',iph,ndorm,lokcs
            endif
         elseif(csrec%phstate.eq.PHSUS) then
! skip composition set number and pre/suffixes at present ....
            susph(nsusp:)=trim(phlista(lokph)%name)//', '
            nsusp=len_trim(susph)+2
         endif
      enddo csloop
   enddo phloop
! we have now sorted stable, entered and dormant phases
   write(unit,10)
10 format(/'List of stable and entered phases'/ &
        '  No tup Name',22x,'Mol.comp. Comp/FU   dGm/RT  Status1  Status2')
!        '  No tup Name',22x,'Mol.comp. At/F.U.   dGm/RT  Status1  Status2')
! come back here for dormant phases
200 continue
   jph=0
   entlist: do iph=1,nent
      trailer=' '
      lokcs=entph(iph)
      csrec=>ceq%phase_varres(lokcs)
      lokph=csrec%phlink
      phname=phlista(lokph)%name
! how do I to know composition set number???  
! Aha!! in phasetuple(phase_varres(lokcs)%phtupx)%compset
      if(phlista(lokph)%noofcs.gt.1) then
         ics=phasetuple(ceq%phase_varres(lokcs)%phtupx)%compset
         chs=char(ichar('0')+ics)
         kp=len_trim(csrec%prefix)
         if(kp.gt.0) then
            csname=csrec%prefix(1:kp)//'_'//phname
         else
            csname=phname
         endif
         kp=len_trim(csrec%suffix)
         if(kp.gt.0) csname=csname(1:len_trim(csname))//'_'//csrec%suffix(1:kp)
         csname=csname(1:len_trim(csname))//'#'//chs//trailer
      else
         csname=phname
      endif
! phase names for composition sets can be larger than 24, remove middle part
      jl=len_trim(csname)
      if(jl.gt.24) then
         csname=csname(1:12)//'..'//csname(jl-9:jl)
      endif
      write(unit,112)phlista(lokph)%alphaindex,csrec%phtupx,csname, &
           csrec%amfu*csrec%abnorm(1),&
           csrec%abnorm(1),csrec%dgm,phlista(lokph)%status1,&
           ceq%phase_varres(lokcs)%status2,ch1
112   format(2i4,1x,a24,1PE10.2,1x,0PF8.2,1PE10.2,2(0p,z8),a1)
!            write(unit,111)jk,csrec%phtupx,csname, &
!                 csrec%abnorm(1),csrec%dgm,&
!                 phlista(lokph)%status1,ceq%phase_varres(lokcs)%status2,ch1
!111         format(2i4,1x,a24,'       0.0',1x0PF8.2,1PE10.2,2(0p,z8),a1)
      if(csrec%dgm.lt.zero) then
         jph=jph+1
         if(jph.gt.10) then
            write(unit,*)' ... remaining phases further from stability'
            exit entlist
         endif
      endif
   enddo entlist
!
   if(ndorm.eq.0) goto 400
   write(unit,210)
210 format(/'List of dormant phases'/ &
         '  No tup Name',22x,'Mol.comp.  Comp/FU  dGm/RT   Status1 Status2')
!         '  No tup Name',22x,'Mol.comp.  At/F.U.  dGm/RT   Status1 Status2')
   jph=0
   dorlist: do iph=1,ndorm
      trailer=' '
      lokcs=dorph(iph)
      csrec=>ceq%phase_varres(lokcs)
      lokph=csrec%phlink
      phname=phlista(lokph)%name
! how do I to know composition set number???  
! Aha!! in phasetuple(phase_varres(lokcs)%phtupx)%compset
      if(phlista(lokph)%noofcs.gt.1) then
         ics=phasetuple(ceq%phase_varres(lokcs)%phtupx)%compset
         chs=char(ichar('0')+ics)
         kp=len_trim(csrec%prefix)
         if(kp.gt.0) then
            csname=csrec%prefix(1:kp)//'_'//phname
         else
            csname=phname
         endif
         kp=len_trim(csrec%suffix)
         if(kp.gt.0) csname=csname(1:len_trim(csname))//'_'//csrec%suffix(1:kp)
         csname=csname(1:len_trim(csname))//'#'//chs//trailer
      else
         csname=phname
      endif
! phase names for composition sets can be larger than 24, remove middle part
      jl=len_trim(csname)
      if(jl.gt.24) then
         csname=csname(1:12)//'..'//csname(jl-9:jl)
      endif
      write(unit,112)phlista(lokph)%alphaindex,csrec%phtupx,csname, &
           csrec%amfu*csrec%abnorm(1),&
           csrec%abnorm(1),csrec%dgm,phlista(lokph)%status1,&
           ceq%phase_varres(lokcs)%status2,ch1
      jph=jph+1
      if(jph.gt.10) then
         write(unit,*)' ... other phases further from stability'
         exit dorlist
      endif
   enddo dorlist
! list suspended phases without composition set numbers
400 continue
   if(nsusp.gt.1) then
      write(unit,300)
300   format(/'List of suspended phases:')
! First indentation 4, for 2nd and later lines 4 also
      call wrice2(unit,2,4,78,1,susph(1:nsusp-3))
   endif
1000 continue
   return
 end subroutine list_sorted_phases

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_all_phases(unit,ceq)
! short list with one line for each phase
! suspended phases merged into one line
   implicit none
   integer unit
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! separate entered/fixed form suspended/dormant
   integer jl,jk,ics,lokph,lokcs,kp,ndorm,nsusp
   character line*80,phname*24,trailer*28,chs*1,csname*36,susph*4096,ch1*1
!   type(gtp_phasetuple), allocatable :: dormant
   TYPE(gtp_phase_varres), pointer :: csrec
   susph=' '
   nsusp=1
   write(unit,10)nooftuples
10  format(/'List of ',i3,' phases'/ &
         '  No tup Name',22x,'Mol.comp. Comp/FU   dGm/RT  Status1  Status2')
!         '  No tup Name',22x,'Mol.comp. At/F.U.   dGm/RT  Status1  Status2')
   jl=0
   trailer=' '
!   write(*,*)'In list_all_phases',noofph
!   allocate(dormant(noofph))
!   dormant=0
   ndorm=0
! come back here for listing dormant phases
20 continue
!
   phloop: do jk=1,noofph
      line=' '
! list in alphabetical order except gas and liquid(s) first
      lokph=phases(jk)
      csloop: do ics=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(ics)
         csrec=>ceq%phase_varres(lokcs)
!         write(*,*)'lpd: 69: ',jk,ics,lokph,lokcs
         if(ndorm.ge.0) then
            if(csrec%phstate.eq.PHDORM) then
               ndorm=ndorm+1
               cycle
            elseif(csrec%phstate.eq.PHSUS) then
! skip composition set number and pre/suffixes at present ....
               susph(nsusp:)=phlista(lokph)%name(1:&
                    len_trim(phlista(lokph)%name))//', '
               nsusp=len_trim(susph)+2
               cycle
            endif
         elseif(csrec%phstate.ne.PHDORM) then
! when ndorm<0 skip all ohases that are suspended, entered or fix
            cycle
         endif
         phname=phlista(lokph)%name
         jl=jl+1
!         write(*,70)'lpd: 70:',phname,phlista(lokph)%noofcs
!70       format(a,a24,5i6)
         if(phlista(lokph)%noofcs.gt.1) then
            chs=char(ichar('0')+ics)
            kp=len_trim(csrec%prefix)
            if(kp.gt.0) then
               csname=csrec%prefix(1:kp)//'_'//phname
            else
               csname=phname
            endif
            kp=len_trim(csrec%suffix)
            if(kp.gt.0) &
                 csname=csname(1:len_trim(csname))//'_'//csrec%suffix(1:kp)
            csname=csname(1:len_trim(csname))//'#'//chs//trailer
         else
            csname=phname
         endif
! phase names for composition sets can be larger than 24, remove middle part
         jl=len_trim(csname)
         if(jl.gt.24) then
            csname=csname(1:12)//'..'//csname(jl-9:jl)
         endif
         if(ceq%phase_varres(lokcs)%phstate.eq.phfixed) then
            ch1='F'
         elseif(ceq%phase_varres(lokcs)%phstate.eq.phentstab) then
            ch1='S'
         elseif(ceq%phase_varres(lokcs)%phstate.eq.phentered) then
            ch1='E'
         elseif(ceq%phase_varres(lokcs)%phstate.eq.phentunst) then
            ch1='U'
         elseif(ceq%phase_varres(lokcs)%phstate.eq.phdorm) then
            ch1='D'
         elseif(ceq%phase_varres(lokcs)%phstate.eq.phsus) then
            ch1='X'
         else
            write(*,*)'3C unknown state: ',ceq%phase_varres(lokcs)%phstate
         endif
!
         if(csrec%amfu.ne.zero) then
            if(csrec%dgm.eq.zero) then
!               write(unit,110)jk,ics,csname, &
               write(unit,110)jk,csrec%phtupx,csname, &
                    csrec%amfu*csrec%abnorm(1),&
                    csrec%abnorm(1),phlista(lokph)%status1,&
                    ceq%phase_varres(lokcs)%status2,ch1
110            format(2i4,1x,a24,1PE10.2,1x,0PF8.2,'       0.0',2(0p,z8),a1)
!110            format(2i4,1x,a24,1PE10.2,1x,0PF9.2,'       0.0',2(0p,z8))
            else
!               write(unit,112)jk,ics,csname, &
               write(unit,112)jk,csrec%phtupx,csname, &
                    csrec%amfu*csrec%abnorm(1),&
                    csrec%abnorm(1),csrec%dgm,&
                    phlista(lokph)%status1,ceq%phase_varres(lokcs)%status2,ch1
112            format(2i4,1x,a24,1PE10.2,1x,0PF8.2,1PE10.2,2(0p,z8),a1)
!112            format(2i4,1x,a24,1PE10.2,1x,0PF9.2,1PE10.2,2(0p,z8))
            endif
         else
!            write(unit,111)jk,ics,csname, &
            write(unit,111)jk,csrec%phtupx,csname, &
                 csrec%abnorm(1),csrec%dgm,&
                 phlista(lokph)%status1,ceq%phase_varres(lokcs)%status2,ch1
111         format(2i4,1x,a24,'       0.0',1x0PF8.2,1PE10.2,2(0p,z8),a1)
!111         format(2i4,1x,a24,'       0.0',1x0PF9.2,1PE10.2,2(0p,z8))
         endif
      enddo csloop
   enddo phloop
   if(ndorm.gt.0) then
      write(unit,200)
200   format(/'List of dormant phases'/ &
           '  No tup Name',22x,'Mol.comp.  Comp/FU  dGm/RT   Status1 Status2')
!           '  No tup Name',22x,'Mol.comp.  At/F.U.  dGm/RT   Status1 Status2')
      ndorm=-1
      goto 20
   endif
! list suspended phases without composition set numbers
   if(nsusp.gt.1) then
      write(unit,300)
300   format(/'List of phases that are suspended:')
! First indentation 4, for 2nd and later lines 4 also
      call  wrice2(unit,2,4,78,1,susph(1:nsusp-3))
   endif
1000 continue
! temporary list all phase tuples
!   do jl=1,nooftuples
!      lokph=phases(phasetuple(jl)%phase)
!      lokcs=phlista(lokph)%linktocs(phasetuple(jl)%compset)
!      write(*,600)jl,phasetuple(jl)%phase,phasetuple(jl)%compset,lokcs,&
!           firsteq%phase_varres(lokcs)%phtupx
!600   format('Phase tuple: ',3i4,' backlink: ',5i4)
!   enddo
   return
 END subroutine list_all_phases

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_global_results(lut,ceq)
! list G, T, P, V and some other things
   implicit none
   integer lut
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character encoded*64
   double precision x1,x2,x3,xn,rtn
!
!   write(kou,*)'gtp3C: output unit: ',lut
   encoded=' '
   call get_state_var_value('T ',x1,encoded,ceq)
   call get_state_var_value('P ',x2,encoded,ceq)
! We must use VS to get SER reference 
   call get_state_var_value('VS ',x3,encoded,ceq)
! this will write error message if any and reset the code
   if(.not.gtp_error_message(0)) then
! no error, list the data
      write(lut,10)x1,x1-273.15,x2,x3
10    format('T= ',F9.2,' K (',F9.2,' C), P= ',1pe11.4,&
        ' Pa, V= ',1pe11.4,' m3')
      rtn=globaldata%rgas*x1
   else
      rtn=one
   endif
! problem with N, should not take into account the atoms/formula units?
   call get_state_var_value('N ',xn,encoded,ceq)
   call get_state_var_value('B ',x2,encoded,ceq)
   if(.not.gtp_error_message(0)) then
      write(lut,11)xn,x2,rtn
11    format('N= ',1pe12.4,' moles, B= ',1pe12.4,' g, RT= ',1pe12.4,' J/mol')
   endif
! we must use suffix S to have values referred to SER
   call get_state_var_value('GS ',x1,encoded,ceq)
   call get_state_var_value('HS ',x2,encoded,ceq)
   call get_state_var_value('SS ',x3,encoded,ceq)
   if(.not.gtp_error_message(0)) then
      write(lut,12)x1,x1/xn,x2,x3
12    format('G= ',1pe12.5,' J, G/N= ',1pe11.4,' J/mol, H= ',1pe11.4,&
           ' J, S= ',1pe11.4,' J/K')
   endif
1000 continue
   return
 end subroutine list_global_results

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_components_result(lut,mode,ceq)
! list one line per component (name, moles, x/w-frac, chem.pot. reference state
! mode 1=mole fractions, 2=mass fractions
   implicit none
   integer lut,mode
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character svtext*64,encoded*64,name*24
   integer ie,kl
   double precision x1,x2,x3,x4,rtn
   encoded=' '
   if(mode.eq.1) then
      write(lut,7)
!7     format('Component name',11x,'Moles',7x,'Mole-fracs  Chem.potent. ',&
7     format('Component name',4x,'Moles',6x,'Mole-fr  Chem.pot/RT  ',&
           'Activities  Ref.state')
   elseif(mode.eq.2) then
      write(lut,9)
9     format('Component name',4x,'Moles',6x,'Mass-fr  Chem.pot/RT  ',&
           'Activities  Ref.state')
   endif
   call get_state_var_value('T ',x1,encoded,ceq)
   rtn=globaldata%rgas*x1
   do ie=1,noofel
      call get_component_name(ie,name,ceq)
      kl=len_trim(name)
      svtext='N('//name(1:kl)//') '
!      write(*,*)'state variable :',svtext
      call get_state_var_value(svtext,x1,encoded,ceq)
      if(gx%bmperr.ne.0) goto 1000
!
      if(mode.eq.1) then
         svtext='X('//name(1:kl)//') '
      elseif(mode.eq.2) then
         svtext='W('//name(1:kl)//') '
      endif
      call get_state_var_value(svtext,x2,encoded,ceq)
      if(gx%bmperr.ne.0) goto 1000
! This should be read from component record .... ???? YES
      svtext='MU('//name(1:kl)//') '
!      write(*,*)'state variable :',svtext
      call get_state_var_value(svtext,x3,encoded,ceq)
      if(gx%bmperr.ne.0) then
         write(*,*)'3C Error line 659: ',trim(svtext),gx%bmperr
         gx%bmperr=0; x3=1.0D2*rtn
      endif
! divide mu with RT, lnac
      if(abs(x3).gt.1.0D-30) then
         x3=x3/rtn
      else
         x3=zero
      endif
      x4=exp(x3)
! reference state, by default "SER (default)" take from component record
!      if(ceq%complist(ie)%phlink.gt.0) then
      encoded=ceq%complist(ie)%refstate
!      else
! default name of reference state
!         encoded='SER (default)'
!      endif
      write(lut,10)name(1:16),x1,x2,x3,x4,encoded(1:16)
!10    format(a,3(1pe12.4),2x,a)
10    format(a,1pe12.4,0pf9.5,2(1pe12.4),2x,a)
   enddo
1000 continue
   return
 end subroutine list_components_result

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_phases_with_positive_dgm(mode,lut,ceq)
! list one line for each phase+comp.set with positive dgm on device lut
! The phases must be dormant or the result is in error.  mode is not used
   implicit none
   integer mode,lut
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character name*24
!   character*10, dimension(-3:2) :: status=&
!        ['SuspendedEntered   ','Fix       ','Dormant   ','Suspended ']
   integer once,iph,lokph,ics,lokcs,kkz,jd
   integer, dimension(:), allocatable :: phtupx
   integer, dimension(:), allocatable :: isort
   double precision xxx
!   write(*,*)'In list_phases_with_positive_dgm'
   once=0
   do iph=1,noofph
      lokph=phases(iph)
      csloop: do ics=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(ics)
         if(ceq%phase_varres(lokcs)%phstate.lt.PHDORM) cycle csloop
!         if(abs(ceq%phase_varres(lokcs)%netcharge).gt.1.0d-6) then
!            write(*,*)'ignoring phase with net charge: ',iph,ics
!            cycle csloop
!         endif
         if(ceq%phase_varres(lokcs)%dgm.gt.1.0D-3) then
            if(once.eq.0) then
               allocate(phtupx(nooftuples))
            endif
            once=once+1
            if(once.eq.1) write(lut,109)
109         format(/' *** There are phases which would like to be stable')
            phtupx(once)=ceq%phase_varres(lokcs)%phtupx
!            write(*,98)once,phtupx(once),phasetuple(phtupx(once))%phase,iph,&
!                 lokcs,ceq%phase_varres(lokcs)%dgm,&
!                 ceq%phase_varres(lokcs)%netcharge
98          format('3C dgm: ',5i4,2(1pe12.4))
         endif
      enddo csloop
   enddo
   goto 1000
   if(once.gt.0) then
      write(lut,110)once
110   format(/' *** ',i3,' Phases which would like to be stable in order')
      allocate(isort(once))
!      call sortrdd(pdgm,once,isort)
!      if(buperr.ne.0) then
!         write(*,*)'Error sorting fractions',buperr
!         goto 1000
!      endif
      do jd=1,once
! add next line when we have sorted
!         isort(jd)=jd
         isort(jd)=phtupx(jd)
! This is getting messy again, the phase tuple index is at present
! the index to phase_varres +1 (as index 1 is the stable reference phase)
!         iph=phasetuple(phtupx(isort(jd)))%phaseix
! removing redundant call to get_phase_compset
!         iph=phasetuple(phtupx(isort(jd)))%ixphase
!         ics=phasetuple(phtupx(isort(jd)))%compset
!         call get_phase_compset(iph,ics,lokph,lokcs)
!         if(gx%bmperr.ne.0) goto 1000
         lokph=phasetuple(phtupx(isort(jd)))%lokph
         lokcs=phasetuple(phtupx(isort(jd)))%lokvares
!         write(*,117)jd,isort(jd),iph,ics,lokcs,phtupx(isort(jd)),&
!              ceq%phase_varres(lokcs)%dgm
117      format('3C Phase: ',2i3,2i5,2i7,1pe10.2)
!         call get_phasetup_name(phasetuple(isort(jd)),name)
!         kkz=test_phase_status(iph,ics,xxx,ceq)
!            write(*,*)'3C: error: ',name,lokcs,kkz
! old kkz.le.2 means entered or fixed
!            if(kkz.le.3) then
! now: kkz= -3,    -2,         -1,         0,           1,         2 
! means SUSPEND, DORMANT, ENTENTED/UNST, ENTERED, ENTERD/STABLE, FIXED
         kkz=ceq%phase_varres(lokcs+1)%phstate
!         if(kkz.ge.PHDORM) then
!            write(lut,120)name,phstate(kkz),ceq%phase_varres(lokcs+1)%dgm
120         format('Phase: ',a,' Status: ',a,' Driving force:',1pe12.4)
!         endif
      enddo
   endif
1000 continue
   return
 end subroutine list_phases_with_positive_dgm


!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_phase_results(iph,jcs,mode,lut,once,ceq)
! list results for a phase+comp.set on lut
! mode specifies the type and amount of results,
! unit digit:   0=mole fraction,      othewise mass fractions
! 10th digit:   0=only composition,   10=also constitution
! 100th digit:  0=value order,        100=alphabetical order
! 1000th digit: 0=all phases,         1000=only stable phases
   implicit none
   integer iph,jcs,mode,lut
   logical once
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   character text*256,phname*24,status*10
   character (len=24), dimension(:), allocatable :: consts
!    character*24, allocatable (:) :: consts
   double precision xmol(maxel),wmass(maxel),totmol,totmass,amount,abv,mindgm
   double precision, dimension(:), allocatable :: ymol
   integer lokph,lokcs,kode,nz,jl,nk,ll,ip,kstat
   mindgm=1.0D-10
   if(ocv()) write(*,*)'mode: ',mode
   if(iph.lt.1 .or. iph.gt.noofph) then
!       write(*,*)'lpr ',iph,jcs,mode
      gx%bmperr=4050; goto 1000
   endif
   lokph=phases(iph)
   if(btest(phlista(lokph)%status1,phhid)) then
! phase is hidden
      gx%bmperr=4119; goto 1000
   endif
!
! .gt.9
!
   if(jcs.lt.0 .or. jcs.gt.phlista(lokph)%noofcs) then
      gx%bmperr=4072; goto 1000
   elseif(jcs.eq.0) then
      jcs=1
   endif
   lokcs=phlista(lokph)%linktocs(jcs)
!   write(*,*)'3C lpr 2: ',jcs,phlista(lokph)%noofcs,lokcs
! get name with pre- and suffix
   call get_phase_name(iph,jcs,phname)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,11)'Phase name: ',iph,jcs,phname
!11 format(a,2i3,'"',a,'"')
   if(mode.ge.1000) then
! if mode>=1000 list stable phases only (dgm<0 )
!      if(ceq%phase_varres(lokcs)%amount(1).eq.zero) then
      if(abs(ceq%phase_varres(lokcs)%netcharge).gt.1.0d-6) then
         if(ceq%phase_varres(lokcs)%phstate.gt.phentered) then
            write(lut,18)phname(1:len_trim(phname)),&
                 ceq%phase_varres(lokcs)%netcharge
18          format('Phase: ',a,' has stable status with net charge: ',F6.3)
            goto 1000
         endif
      endif
      if(ceq%phase_varres(lokcs)%amfu.eq.zero) then
! skip phases with zero amount unless expcitly stable or positive dgm
         if(ceq%phase_varres(lokcs)%dgm.eq.zero) then
!            if(ceq%phase_varres(lokcs)%phstate.ne.PHFIXED) goto 1000
            if(ceq%phase_varres(lokcs)%phstate.lt.PHENTSTAB) goto 1000
         elseif(ceq%phase_varres(lokcs)%dgm.lt.mindgm) then
            goto 1000
         endif
      endif
   endif
! phase status (except hidden) .... use get_phase_status instead ???
!   if(btest(ceq%phase_varres(lokcs)%status2,cssus)) then
!      if(btest(ceq%phase_varres(lokcs)%status2,csfixdorm)) then
   if(ceq%phase_varres(lokcs)%phstate.eq.PHDORM) then
      status='Dormant'
      kstat=4
! skip dormant phases unless positive driving force
!      if(ceq%phase_varres(lokcs)%dgm.le.mindgm) goto 1000
   elseif(ceq%phase_varres(lokcs)%phstate.eq.PHSUS) then
! skip suspended phases
      status='Suspended'
      goto 1000
!      if(btest(ceq%phase_varres(lokcs)%status2,csfixdorm)) then
   elseif(ceq%phase_varres(lokcs)%phstate.eq.PHFIXED) then
      status='Fixed'
      kstat=2
   else
      status='Entered'
      kstat=1
! skip phase with net charge
!      if(abs(ceq%phase_varres(lokcs)%netcharge).gt.1.0D-6) goto 1000
! skip entered phases that have positive driving force, why??
!      if(ceq%phase_varres(lokcs)%dgm.gt.zero) goto 1000
   endif
   if(phname(1:1).lt.'A' .or. phname(1:1).gt.'Z') then
! in some cases unprintable phase names appears!!
      write(lut,19)iph,jcs,lokph,lokcs
19    format('Illegal  phase name: ',10i5)
   endif
!X   write(lut,20)phname,status,ceq%phase_varres(lokcs)%dgm
20  format(/'Phase: ',A,' Status: 'A,' Driving force: ',1PE12.4)
!------------------------
!   xmol=zero
!   wmass=zero
   call calc_phase_molmass(iph,jcs,xmol,wmass,totmol,totmass,amount,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'Error: ',gx%bmperr; goto 1000
   endif
!   write(*,99)'xmol: ',xmol
!99 format(a,6(1pe12.4))
   kode=mod(mode,10)
!   write(*,*)'lpr: ',mode,kode
   abv=ceq%phase_varres(lokcs)%abnorm(1)
! a shorter output
!   write(lut,700)phname,status(1:1),totmol,totmass*0.001, &
!        amount*ceq%rtn*ceq%phase_varres(lokcs)%gval(3,1),&
!        ceq%phase_varres(lokcs)%amfu,abv,ceq%phase_varres(lokcs)%dgm
! try to fill upp the phase name with '.'
   nz=len_trim(phname)
   phname(nz+1:)='.........................'
   if(kode.eq.0) then
! The volume value here is WRONG: ceq%phase_varres(lokcs)%gval(3,1) !!! ???
      if(once) write(lut,699)'Moles     '
      once=.FALSE.
699   format(/'Name                Status ',a,' Volume',&
           '    Form.Units Cmp/FU dGm/RT  Comp:')
!          '    Form.U      At/FU DGM    Fracs:')
      write(lut,700)phname,status(1:1),totmol,&
           amount*ceq%rtn*ceq%phase_varres(lokcs)%gval(3,1),&
           ceq%phase_varres(lokcs)%amfu,abv,ceq%phase_varres(lokcs)%dgm,'X:'
!      phase status moles/mass   (volume FU)  Atomes/FU DGM Content
700 format(a,1x,a,  1pe11.3,     2(1pe10.2),1x,0pF7.2,1pe10.2,2x,a)
!X      if(abs(ceq%phase_varres(lokcs)%netcharge).gt.1.0D-6) then
!X         write(lut,28)totmol,totmass*0.001, &
!X              amount*ceq%rtn*ceq%phase_varres(lokcs)%gval(3,1),&
!X              ceq%phase_varres(lokcs)%netcharge
!X      else
!X         write(lut,25)totmol,totmass*0.001, &
!X              amount*ceq%rtn*ceq%phase_varres(lokcs)%gval(3,1)
!X      endif
!X      write(lut,21)ceq%phase_varres(lokcs)%amfu,abv
21    format('Formula Units: ',1pe12.4,', Moles of atoms/FU: ',1pe12.4,&
           ', Molar content:')
   else
      if(once) write(lut,699)'Mass      '
      once=.FALSE.
      write(lut,700)phname,status(1:1),totmass*0.001, &
           amount*ceq%rtn*ceq%phase_varres(lokcs)%gval(3,1),&
           ceq%phase_varres(lokcs)%amfu,abv,ceq%phase_varres(lokcs)%dgm,'W:'
!X      write(lut,25)totmol,totmass*0.001,&
!X           amount*ceq%rtn*ceq%phase_varres(lokcs)%gval(3,1)
!X      write(lut,22)ceq%phase_varres(lokcs)%amfu,abv
22    format('Formula Units: ',1pe12.4,', Moles of atoms/FU: ',1pe12.4,&
           ', Mass fractions:')
   endif
25  format('Moles',1PE12.4,', Mass',1PE12.4,' kg, Volume',1PE12.4,' m3')
28  format('Moles',1PE11.3,' Mass',1PE11.3,' kg, Volume',1PE11.3,' m3,',&
         ' Charge: ',1pe11.3)
! composition
   nz=noofel
   allocate(consts(nz))
   consts=' '
   do jl=1,nz
      consts(jl)=splista(ceq%complist(jl)%splink)%symbol
   enddo
!    write(*,187)'lpr: ',consts
!187 format(a,20(1x,a2))
   if(kode.eq.0) then
      call format_phase_composition(mode,nz,consts,xmol,lut)
   else
      call format_phase_composition(mode,nz,consts,wmass,lut)
   endif
   deallocate(consts)
   if(gx%bmperr.ne.0) goto 1000
!-------------------------------------
! constitution only if nonzero tenth-digit of mode or if GAS
300 continue
   if(.not.btest(phlista(lokph)%status1,PHGAS)) then
      if(mod(mode/10,10).le.0) goto 900
   endif
   write(lut,310,advance='no')
310  format(' Constitution: ')
!---------------
   nk=0
   sublatloop: do ll=1,phlista(lokph)%noofsubl
      nz=phlista(lokph)%nooffr(ll)
!      if(phlista(lokph)%noofsubl.gt.1) then
      if(size(ceq%phase_varres(lokcs)%sites).gt.1) then
!         write(lut,320)ll,nz,phlista(lokph)%sites(ll)
         if(ll.gt.1) then
            write(lut,319)ll,nz,ceq%phase_varres(lokcs)%sites(ll)
         else
            write(lut,320)ll,nz,ceq%phase_varres(lokcs)%sites(ll)
         endif
319      format(15x,'Sublattice ',i2,' with ',i5,' constituents and ',&
              F12.6,' sites')
320      format('Sublattice ',i2,' with ',i5,' constituents and ',&
              F12.6,' sites')
!      elseif(phlista(lokph)%sites(ll).eq.one) then
      elseif(ceq%phase_varres(lokcs)%sites(ll).eq.one) then
         write(lut,321)nz
321      format('There are ',i5,' constituents:')
      else
!         write(lut,322)nz,phlista(lokph)%sites(ll)
         write(lut,322)nz,ceq%phase_varres(lokcs)%sites(ll)
322      format('Single lattice with ',i5,' constituents and ',&
              F12.6,' sites')
      endif
      text=' '; ip=1
      allocate(consts(nz))
      allocate(ymol(nz))
      consts=' '
      do jl=1,nz
!         jcons=splista(phlista(lokph)%constitlist(nk+jl))%alphaindex
         consts(jl)=' '
         if(phlista(lokph)%constitlist(nk+jl).gt.0) then
            consts(jl)=splista(phlista(lokph)%constitlist(nk+jl))%symbol
         else
            consts(jl)='*'
         endif
         ymol(jl)=ceq%phase_varres(lokcs)%yfr(nk+jl)
      enddo
      call format_phase_composition(mode,nz,consts,ymol,lut)
      deallocate(consts)
      deallocate(ymol)
      if(gx%bmperr.ne.0) goto 1000
      nk=nk+nz
   enddo sublatloop
900 continue
! write an empty line after each phase ...
   write(lut,*)
1000 continue
   return
 end subroutine list_phase_results

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_short_results(lut,ceq)
! list short results for all stable phases (for debugging) lut
   implicit none
   integer lut
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer iph,ics,lokph,lokcs,i1,i2
   phaseloop: do iph=1,noofph
      lokph=phases(ics)
      compsets: do ics=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(ics)
         if(ceq%phase_varres(lokcs)%phstate.ge.PHENTSTAB) then
            write(lut,110)phlista(lokph)%name,ics,ceq%phase_varres(lokcs)%amfu
110         format(a,i2,4(1pe12.4))
         endif
      enddo compsets
   enddo phaseloop
1000 continue
   return
 end subroutine list_short_results

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine format_phase_composition(mode,nv,consts,vals,lut)
! list composition/constitution in alphabetical or value order
! entalsiffra 0 mole fraction, 1 mass fraction, 3 mole percent, 4 mass percent
! tiotalsiffra alphabetical order ... ??
! mode >100 else alphabetical order
! nv is number of components/constitunents (in alphabetical order in consts)
! components/constituents in consts, fractions in vals
   implicit none
   integer nv,mode,lut
   character consts(nv)*(*)
   double precision vals(nv)
!\end{verbatim}
   integer maxl,jl,kp,ncol,nrow2,nvrest,n1,nempty,n3r,n4r
   character names(4)*12
   integer, dimension(:), allocatable :: isort
! 3-13 position name, 12 positions value (1pe12.5), 2 positions separator
! NOTE components can have negative fractions but not constituents
! so leave one blank after component names
! Constituents with names longer than 13 will be written A23456..12345
! with 6 initial characters, two dots and then the 5 last characters
! Max 4 columns with 18 positions(=72) plus 3*2=6 position separator,
! min 3 columns with 24 positions(=72) plus 2*2=4 position separator
!
! max length of names and number of columns
   maxl=0
   do jl=1,nv
      kp=len_trim(consts(jl))
      if(kp.gt.maxl) then
         maxl=kp
      endif
   enddo
   if(maxl.le.4) then
! use 4 columns if names are short
      ncol=4
   else
      ncol=3
   endif
! number of rows is needed to have valuses in columns decending like:
!  FE  0.75 SI 0.05 Ti 0.02 C 0.01
!  CR  0.20 Mn 0.04 V  0.01
!-----------------------------------
   nrow2=(nv+ncol-1)/ncol
! always use isort for the order, if alphabetical isort(i)=i
   allocate(isort(nv+4))
!   if(mode.ge.100) then
!   write(*,*)'3C mode: ',mode,mod(mode,100),mode-100*mod(mode,100)
!   if(mod(mode,10).eq.0) then
   if(mod(mode/100,10).eq.0) then
! value order
      call sortrdd(vals,nv,isort)
      if(buperr.ne.0) then
         write(*,*)'Error sorting fractions',buperr
         gx%bmperr=buperr; goto 1000
      endif
   else
! if alphabetical order just set isort(i)=i, same index as for vals
      do jl=1,nv
         isort(jl)=jl
      enddo
   endif
! list constituents in the order of isort
   if(ncol.eq.4) then
! All names max 4 characters, 4 columns: 1 + 4+1+13+2 +20 +20 +18 =  79
      nvrest=nv
      n1=1
! number of empty colums in last row is 4*nrow2-nv
      nempty=4*nrow2-nv
! 3rd and 4th column may start from one or two indices less
      n3r=2*nrow2
      n4r=3*nrow2
      if(nempty.eq.3) then
         n3r=n3r-1
         n4r=n4r-2
      elseif(nempty.eq.2) then
         n4r=n4r-1
      endif
100   continue
! this can be quite complicated as last row may be partially empty as 
      if(nvrest.ge.4) then
         names(1)=consts(isort(n1))
         names(2)=consts(isort(n1+nrow2))
         names(3)=consts(isort(n1+n3r))
! 4th column may be empty after first row
         if(n1+n4r.le.nv) then
            names(4)=consts(isort(n1+n4r))
            write(lut,110)names(1)(1:4),vals(n1),&
                 names(2)(1:4),vals(n1+nrow2),names(3)(1:4),vals(n1+n3r),&
                 names(4)(1:4),vals(n1+n4r)
110         format(1x,a,1x,1pe13.5,3(2x,a,1x,1pe13.5))
            nvrest=nvrest-4
         else
            write(lut,110)names(1)(1:4),vals(n1),&
                 names(2)(1:4),vals(n1+nrow2),names(3)(1:4),vals(n1+n3r)
            nvrest=nvrest-3
         endif
         n1=n1+1
      else
! last row can be 1 to 3 columns
         names(1)=consts(isort(n1))
         if(nvrest.gt.1) then
            names(2)=consts(isort(n1+nrow2))
            if(nvrest.gt.2) then
               names(3)=consts(isort(n1+n3r))
               write(lut,110)names(1)(1:4),vals(n1),&
                    names(2)(1:4),vals(n1+nrow2),names(3)(1:4),vals(n1+n3r)
            else
               write(lut,110)names(1)(1:4),vals(n1),&
                    names(2)(1:4),vals(n1+nrow2)
            endif
         else
            write(lut,110)names(1)(1:4),vals(n1)
         endif
         nvrest=0
      endif
      if(nvrest.gt.0) goto 100
   else
! All listed names have max 13 characters, longer names are truncated
      nvrest=nv
      n1=1
! number of empty columns in last row
      nempty=3*nrow2-nv
! 3rd column may start from an indices less
      n3r=2*nrow2
      if(nempty.eq.2) then
         n3r=n3r-1
      endif
200   continue
      if(nvrest.ge.3) then
         names(1)=consts(isort(n1))
         names(2)=consts(isort(n1+nrow2))
         if(n1+2*nrow2.le.nv) then
            names(3)=consts(isort(n1+n3r))
            write(lut,210)names(1),vals(n1),names(2),vals(n1+nrow2),&
                 names(3),vals(n1+n3r)
210         format(1x,a,1pe12.5,2(2x,a,1pe12.5))
            nvrest=nvrest-3
         else
            write(lut,210)names(1),vals(n1),names(2),vals(n1+nrow2)
            nvrest=nvrest-2
         endif
         n1=n1+1
      else
! last row can be 1 or 2 columns
         names(1)=consts(isort(n1))
         if(nvrest.gt.1) then
            names(2)=consts(isort(n1+nrow2))
            write(lut,210)names(1),vals(n1),names(2),vals(n1+nrow2)
         else
            write(lut,210)names(1),vals(n1)
         endif
         nvrest=0
      endif
      if(nvrest.gt.0) goto 200
   endif
!
1000 continue
   return
 end subroutine format_phase_composition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_many_formats(cline,last,ftyp,unit1)
! lists all data in different formats: SCREEN/TDB/MACRO/LaTeX/ODB
!                                         1    2    3     4    5
! unfinished
   implicit none
   character cline*(*)
   integer last,unit1,ftyp
!\end{verbatim}
   integer iph,ipos,kousave,unit,isp
   character text*64, text2*2000,fil*64
   character date*8,CHTD*1
! if not screen then ask for file name
! for screen outpot of file use /option= ...
   if(ftyp.ne.1) then
      call gparcd('Output file: ',cline,last,1,fil,'database',q1help)
      ipos=len_trim(fil)
! it is impossible to have a blank name here, check if there is an extension
      iph=index(fil,'.')
      if(iph.gt.0) then
! do not allow empty extensions
         if(fil(iph+1:iph+1).ne.' ') ipos=0
      endif
      if(ipos.gt.0) then
         if(ftyp.eq.2) then
! TDB file a la TC
            fil(ipos+1:)='.TDB'
         elseif(ftyp.eq.3) then
            fil(ipos:)='.OCM'
         elseif(ftyp.eq.4) then
            fil(ipos:)='.tex'
         elseif(ftyp.eq.5) then
! PDB file
            fil(ipos:)='.PDB'
         endif
      endif
! check if file exists ... overwriting not allowed ...
      open(unit=31,file=fil,access='sequential',status='new',err=900)
      kousave=unit
      unit=31
   endif
   call date_and_time(date)
   select case(ftyp) 
   case default
      write(kou,*)'No such format'
!----------------------------------------------------------
! This can be written to file using the /output option
   case(1) ! ftyp=1 SCREEN format
      call list_all_elements(kou)
      if(gx%bmperr.ne.0) goto 1000
      call list_all_species(kou)
      if(gx%bmperr.ne.0) goto 1000
      call list_all_funs(kou)
      if(gx%bmperr.ne.0) goto 1000
      do iph=1,noph()
         call list_phase_data(iph,' ',kou)
         if(gx%bmperr.ne.0) goto 1000
      enddo
! list reference phase last
      iph=0
      call list_phase_data(0,' ',kou)
! finally list the data bibliography
      write(kou,*)
      call list_bibliography(' ',kou)
!--------------------------------------------------------------
! write on unit
   case(2) ! ftyp=2 TDB format
! CHTD1 keeps track of type definitions, note: incremented before use
      CHTD='0'
      write(unit,106)date(1:4),date(5:6),date(7:8)
106   format('$ Database file written by Open Calphad ',a,'-',a,'-',a/)
      call list_all_elements2(unit)
      write(unit,107)
107   format(/'$ =================',/)
      text=' '
      sploop: do isp=1, nosp()
! skip vacancy species and species that are elements
         iph=species(isp)
         ipos=1
         call list_species_data2(text,ipos,iph)
! not very logical, using species index below and location above ... suck
         if(testspstat(isp,SPEL) .or. testspstat(isp,SPVA)) then
            cycle sploop
         endif
         write(unit,110)text(1:len_trim(text))
110      format('SPECIES ',A,' !')     
      end do sploop
      write(unit,107)
      text2=' '
! skip the first two functions which are R and RTLNP (using R)
! write RTLNP in correct TDB form here
      text2='FUNCTION RTLNP 10 8.31451*R*LN(1.0D-5*P); 20000 N !'
      write(unit,112)text2(1:len_trim(text2))
112   format(a)
!
      tpfuns: do iph=3, notpf()! freetpfun-1
         text2='FUNCTION '
         call list_tpfun(iph,0,text2(10:))
! skip functions with names staring with _ as they are parameters
         if(text2(10:10).eq.'_') cycle tpfuns
! for the remaining functions OC writes them with = T_low ...
! and for TC one must remove the = sign
         ipos=index(text2,'=')
         text2(ipos:ipos)=' '
! then add a ! at the end
         ipos=len_trim(text2)
         text2(ipos+1:)=' !'
         call  wrice2(unit,0,8,78,1,text2)
      end do tpfuns
      write(unit,107)
      write(unit,130)
130   format(/'TYPE_DEFINITION % SEQ * !'/ &
          'DEFINE_SYSTEM_DEFAULT ELEMENT 2 !'/ &
          'DEFAULT_COMMAND DEF_SYS_ELEMENT  VA /- !'/)
      write(unit,107)
      do iph=1, noph()
         call list_phase_data2(iph,ftyp,CHTD,unit)
      enddo
      write(unit,107)
      write(unit,140)
140   format(/' LIST_OF_REFERENCES'/ ' NUMBER  SOURCE')
      call list_bibliography(' ',unit)
      write(unit,141)
141   format('!')
      close(unit)
!--------------------------------------------------------------
   case(3) ! ftyp=3 MACRO format
      write(kou,*)'MACRO not implemented yet'
!--------------------------------------------------------------
   case(4) ! ftyp=4 LATEX format
      write(kou,*)'LaTeX not implemented yet'
!--------------------------------------------------------------
   case(5) ! ftyp=5 Open Calphad TDB format
!      write(kou,*)'PDB not implemented yet'
      call write_pdbformat(unit)
   end select
!--------------------------------------------------------------
   goto 1000
! error
900 continue
!  write(kou,*)'File already exist, overwriting not allowed'
   close(31)
   gx%bmperr=4190
1000 continue
!   unit=kousave
   return
 end subroutine list_many_formats

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_phase_model(iph,ics,lut,CHTD,ceq)
! list model (no parameters) for a phase on lut
   implicit none
   integer iph,ics,lut
   TYPE(gtp_equilibrium_data), pointer :: ceq
   character CHTD*1
!\end{verbatim}
   character phname*24,l78*78
!   integer, dimension(maxsubl) :: endm,ilist
   integer lokcs,knr,kmr,ll,ip,lokph,ftyp
   TYPE(gtp_fraction_set) :: disfra
   type(gtp_phase_add), pointer :: addrec
   double precision rl
! screen
   ftyp=1
! if ics=0 list fractions for all composition sets
   lokph=phases(iph)
! name, model name
! sublattices, status,
! additions
! sites, constituents and fractions in each disordered constituents
! number of disordered sublattices
! sites, constituents and fractions in each disordered constituents
   if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
!      write(*,*)'No subch composition set'
      gx%bmperr=4072; goto 1000
   elseif(ics.eq.0) then
      ics=1
   endif
   lokcs=phlista(lokph)%linktocs(ics)
   call get_phase_name(iph,ics,phname)
   write(lut,110)phname,phlista(lokph)%models(1:40),&
        phlista(lokph)%noofsubl,phlista(lokph)%status1,&
        ceq%phase_varres(lokcs)%status2
110 format(a,' model: ',a/'  Number of sublattices: ',i2,&
        ', status: ',z8,1x,z8,5x)
   addrec=>phlista(lokph)%additions
   lastadd: do while(associated(addrec))
      call list_addition(lut,CHTD,phname,ftyp,addrec)
      addrec=>addrec%nextadd
   enddo lastadd
! return here if more composition sets
200 continue
   rl=zero
   knr=0
   kmr=0
! return here for each sublattice
   do ll=1,phlista(lokph)%noofsubl
      rl=rl+one
      kmr=kmr+phlista(lokph)%nooffr(ll)
      l78='Subl. '; ip=7
      call wrinum(l78,ip,2,0,rl)
      l78(ip:)=', sites: '; ip=ip+9
!      call wrinum(l78,ip,6,0,phlista(lokph)%sites(ll))
      call wrinum(l78,ip,6,0,ceq%phase_varres(lokcs)%sites(ll))
      l78(ip:)=', const.: '; ip=ip+10
! return here for each new constituent in this sublattice
320    continue
      knr=knr+1
      if(phlista(lokph)%constitlist(knr).gt.0) then
         l78(ip:)=splista(phlista(lokph)%constitlist(knr))%symbol
      else
         l78(ip:)='*'
      endif
      ip=len_trim(l78)+2
      l78(ip-1:ip-1)='='
! The fractions for normal sublattice done by list result or list phase-const
      call wrinum(l78,ip,6,0,ceq%phase_varres(lokcs)%yfr(knr))
      l78(ip:ip+1)=', '
      ip=ip+2
      if(ip.gt.60) then
         write(lut,330)l78(1:ip-3)
330       format(2x,a)
         l78=' '
         ip=4
      endif
      if(knr.lt.kmr) goto 320
      if(ip.gt.4) write(lut,330)l78(1:ip-3)
   enddo
   if(btest(phlista(lokph)%status1,PHMFS)) then
! the phase has disordered fractions
! ?? does the = here make a copy?  I just want a pointer ...
      disfra=ceq%phase_varres(lokcs)%disfra
      lokcs=disfra%varreslink
      if(disfra%ndd.eq.1) then
         write(lut,410)disfra%latd
410      format(4x,'Disordred fractions adding all fractions from all ',&
              i2,' sublattices together')
      else
         write(lut,420)disfra%latd
420      format(4x,'Disordred fractions adding fractions from first ',i2,&
              ' sublattices together in'/&
              4x,'the first disordered sublattice',&
              ' and the remaining fractions in the second.')
      endif
! write the disordered constituents and fractions
      ll=0
      rl=zero
      knr=0
      kmr=0
! return here for second sublattice (if any)
430   continue
      ll=ll+1
      rl=rl+one
      kmr=kmr+disfra%nooffr(ll)
      l78='Subl. '; ip=7
      call wrinum(l78,ip,2,0,rl)
      l78(ip:)=', sites: '; ip=ip+9
      call wrinum(l78,ip,6,0,disfra%dsites(ll))
      l78(ip:)=', const.: '; ip=ip+10
! return here for each new constituent in this sublattice
440   continue
      knr=knr+1
      l78(ip:)=splista(disfra%splink(knr))%symbol
! list fractions in disordered sublattice as this is the only place for that
      ip=len_trim(l78)+2
      l78(ip-1:ip-1)='='
      call wrinum(l78,ip,6,0,ceq%phase_varres(lokcs)%yfr(knr))
      l78(ip:)=','
      ip=ip+2
      if(ip.gt.60) then
         write(lut,450)l78(1:ip-3)
450       format(4x,a)
         l78=' '
         ip=4
      endif
      if(knr.lt.kmr) goto 440
      if(ip.gt.4) write(lut,330)l78(1:ip-3)
      if(ll.lt.disfra%ndd) goto 430
   endif
1000 continue
   return
 end subroutine list_phase_model

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_phase_data(iph,CHTD,lut)
! list parameter data for a phase on unit lut
   implicit none
   integer iph,lut
   character CHTD*1
!\end{verbatim} %+
   integer typty,parlist,typspec,lokph,nsl,nk,ip,ll,jnr,ics,lokcs
   integer nint,ideg,ij,kk,iel,ncsum,kkx,kkk,jdeg,iqnext,iqhigh,lqq,nz,ik
   integer intpq,linkcon,ftyp
   character text*2048,phname*24,prop*32,funexpr*1024
   character special*8
!   integer, dimension(2,3) :: lint
! ?? increased dimension of lint ??
   integer, dimension(2,5) :: lint
   integer, dimension(maxsubl) :: endm,ilist
   logical subref,noelin1
   type(gtp_fraction_set), pointer :: disfrap
! a smart way to have an array of pointers
   TYPE intrecarray 
      type(gtp_interaction), pointer :: p1
   end TYPE intrecarray
   type(intrecarray), dimension(20) :: intrecstack
   type(gtp_property), pointer :: proprec
   type(gtp_interaction), pointer :: intrec
   type(gtp_endmember), pointer :: endmemrec
   TYPE(gtp_fraction_set) :: disfra
   TYPE(gtp_phase_add), pointer :: addrec
!
! output on screen
   ftyp=1
   if(iph.lt.0 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   elseif(noofel.eq.0) then
! this needed as there is a reference phase with iph=0 when there are elements
      goto 1000
   endif
!   write(*,*)'lpd 1:',iph,phases(iph)
   if(iph.gt.0) then
      lokph=phases(iph)
   else
      lokph=0
   endif
   ics=1
   phname=phlista(lokph)%name
   nsl=phlista(lokph)%noofsubl
   special=' '
! indicate some status bit specially
   if(btest(phlista(lokph)%status1,PHFORD)) special(1:1)='F'
   if(btest(phlista(lokph)%status1,PHBORD)) special(1:1)='B'
   if(btest(phlista(lokph)%status1,PHSORD)) special(1:1)='S'
   if(btest(phlista(lokph)%status1,PHIONLIQ)) special(1:1)='I'
   if(btest(phlista(lokph)%status1,PHMFS)) special(2:2)='D'
! This subroutine is independent of current equilibrium, use firsteq
!   write(lut,10)phname,phlista(lokph)%status1,special,&
!        nsl,(phlista(lokph)%sites(ll),ll=1,nsl)
   lokcs=phlista(lokph)%linktocs(ics)
   write(lut,10)phname,phlista(lokph)%status1,special,&
        nsl,(firsteq%phase_varres(lokcs)%sites(ll),ll=1,nsl)
10  format(/'Phase: ',A,', Status: ',Z8,2x,a/'  Subl:',I3,10(1x,F7.3))
   nk=0
   text='Constituents: '
   ip=15
   sublatloop: do ll=1,nsl
      constloop: do ik=1,phlista(lokph)%nooffr(ll)
         nk=nk+1
         jnr=phlista(lokph)%constitlist(nk)
         if(jnr.gt.0) then
            text(ip:)=splista(jnr)%symbol
         else
            text(ip:)='*'
         endif
         ip=len_trim(text)+1
!         text(ip:ip)=','
         text(ip:ip)=' '
         ip=ip+1
      enddo constloop
      text(ip-1:ip)=':  '
      ip=ip+1
   enddo sublatloop
   call wrice2(lut,2,4,78,-1,text)
!    write(lut,17)text(1:ip)
!17  format(A)
! additions
   addrec=>phlista(lokph)%additions
   lastadd: do while(associated(addrec))
      call list_addition(lut,CHTD,phname,ftyp,addrec)
      addrec=>addrec%nextadd
   enddo lastadd
60 continue
! parameters for end members using site fractions
   if(btest(phlista(lokph)%status1,PHMFS)) then
      subref=.FALSE.
   else
      subref=.TRUE.
   endif
   parlist=1
!--------------------------------------------------
! return here to list disordered parameters
100 continue
! parlist changed below for disordered fraction set
   if(parlist.eq.1) then
      endmemrec=>phlista(lokph)%ordered
   else
      if(ocv()) write(*,*)'Listing disordred parameters ',nsl
      endmemrec=>phlista(lokph)%disordered
      disfrap=>firsteq%phase_varres(lokcs)%disfra
   endif
   endmemberlist: do while(associated(endmemrec))
      do ll=1,nsl
!         ilist(ll)=emlista(lokem)%fraclinks(ll,1)
         ilist(ll)=endmemrec%fraclinks(ll,1)
         if(ilist(ll).gt.0) then
            if(parlist.eq.2) then
! what is disfra here??!!
               endm(ll)=disfra%splink(ilist(ll))
            else
               endm(ll)=phlista(lokph)%constitlist(ilist(ll))
            endif
         else
! wildcard, write '*'
            endm(ll)=-99
         endif
      enddo
      nint=0
      ideg=0
      call encode_constarr(text,nsl,endm,nint,lint,ideg)
      if(gx%bmperr.ne.0) goto 1000
      proprec=>endmemrec%propointer
      ptyloop: do while(associated(proprec))
         ij=proprec%proptype
         if(ij.ge.100) then
            typty=ij/100
            typspec=mod(ij,100)
         else
            typty=ij
         endif
         if(typty.gt.0 .and. typty.le.ndefprop) then
            prop=propid(typty)%symbol
            if(parlist.eq.2) then
! disordered endmember parameter
               kk=len_trim(prop)+1
               prop(kk:kk)='D'
            endif
            if(btest(propid(typty)%status,IDELSUFFIX)) then
! property like ZZ&<element>(phase,constituent array)
! the element index should be in typsepc
               iel=typspec
               if(iel.ge.0 .and. iel.le.noofel) then
!                  prop=propid(typty)%symbol
                  prop=prop(1:len_trim(prop))//'&'&
                       //ellista(elements(iel))%symbol
               else
                  gx%bmperr=4082; goto 1000
               endif
            elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
! property like mobility, MQ&<constituent#sublat>(phase,constituent array)
! the suffix is a constituent
               iel=typspec
               if(iel.gt.0 .and. iel.le.phlista(lokph)%tnooffr) then
                  if(parlist.eq.2) then
! we must consider parlist, take disordered constituent list
! we have no current equilibrium record but can use firsteq!!
!                     lokcs=phlista(lokph)%linktocs(1)
!                     write(*,*)'3C: endmember typspec 1: ',iel
!                     write(*,*)'3C splink: ',disfrap%splink
                     linkcon=disfrap%splink(iel)
!                     write(*,*)'3C: endmember typspec 2: ',linkcon
                     ll=0
!                     ll=1
! linkcon has nothing to do with which sublattice, ignore ll
!                     if(linkcon.gt.disfrap%nooffr(1)) ll=2
                     prop=prop(1:len_trim(prop))//'&'&
                          //splista(linkcon)%symbol
!                     write(*,*)'3C We are here',linkcon,disfrap%nooffr(1),ll
                     prop=prop(1:len_trim(prop))
!                     goto 120
                     goto 121
                  else
                     linkcon=phlista(lokph)%constitlist(iel)
                     if(linkcon.le.0) then
                        write(*,*)'Illegal use of wildcard 1'
                        gx%bmperr=4286; goto 1000
                     endif
                     prop=prop(1:len_trim(prop))//'&'&
                          //splista(linkcon)%symbol
! also add the sublattice number ...
                     ncsum=0
                     do ll=1,phlista(lokph)%noofsubl
                        ncsum=ncsum+phlista(lokph)%nooffr(ll)
                        if(iel.le.ncsum) goto 120
                     enddo
                  endif
! error if sublattice not found
                  write(kou,*)'Error in constituent depended parameter id'
                  gx%bmperr=4287; goto 1000
! jump here to append sublattice
120               continue
!                  write(*,*)'property 1: ',prop(1:10),ll
                  if(ll.gt.1) then
                     prop=prop(1:len_trim(prop))//'#'//char(ll+ichar('0'))
!                  else
!                     prop=prop(1:len_trim(prop))//'#'//char(ll+ichar('0'))
                  endif
121               continue
               else
                  write(kou,*)'lpd 7B: ',iel,typty
                  gx%bmperr=4082; goto 1000
               endif
            endif
         else
! unknown property ...
            write(*,*)'unknown property type xx: ',ij,typty,typspec
            prop='ZZ'
         endif
! note changes here must be repeated for interaction parameters below
         write(funexpr,200)prop(1:len_trim(prop)),&
              phname(1:len_trim(phname)),text(1:len_trim(text))
200      format(A,'(',A,',',A,') ')
         ip=len_trim(funexpr)+1
! subtract reference states
         if(subref .and. typty.eq.1) then
            call subrefstates(funexpr,ip,lokph,parlist,endm,noelin1)
            if(noelin1) then
! this can happen for ionic liquids with just neutrals in sublattice 2
! replace the constituent in sublattice 1 with "*" !!!
!               write(*,*)'before: ',funexpr(1:ip)
               kk=index(funexpr,',')
               ik=index(funexpr,':')
               funexpr(kk+1:)='*'//funexpr(ik:)
               ip=len_trim(funexpr)+2
!               write(*,*)'after:  ',funexpr(1:ip)
            endif
         endif
! this writes the expression
         call list_tpfun(proprec%degreelink(0),1,funexpr(ip:))
         ip=len_trim(funexpr)
         funexpr(ip+1:)=' '//proprec%reference
         ip=len_trim(funexpr)
! nice output over several lines if needed with indentation 12 spaces
         call wrice2(lut,2,12,78,1,funexpr(1:ip))
         proprec=>proprec%nextpr
      enddo ptyloop
      if(endmemrec%noofpermut.gt.1) then
         intpq=0
         if(associated(endmemrec%intpointer)) then
            intpq=endmemrec%intpointer%antalint
         endif
!         write(kou,207)endmemrec%antalem,endmemrec%noofpermut,intpq
207      format('@$ Endmember, permutations, interaction: ',3i5)
      endif
      endmemrec=>endmemrec%nextem
   enddo endmemberlist
!-----------------------------------------------------------------------
! parameters for interactions using site fractions
   if(parlist.eq.1) then
      endmemrec=>phlista(lokph)%ordered
   else
      endmemrec=>phlista(lokph)%disordered
   endif
   intlist1: do while(associated(endmemrec))
      intrec=>endmemrec%intpointer
      if(associated(intrec)) then
!         write(*,*)'intlist 1B: ',intrec%status
         do ll=1,nsl
            kkx=endmemrec%fraclinks(ll,1)
            if(kkx.eq.-99) then
! wildcard
               endm(ll)=-99
            elseif(parlist.eq.2) then
               endm(ll)=disfra%splink(kkx)
            else
               endm(ll)=phlista(lokph)%constitlist(kkx)
            endif
         enddo
      endif
      nint=0
      intlist2: do while(associated(intrec))
         nint=nint+1
         intrecstack(nint)%p1=>intrec
         lint(1,nint)=intrec%sublattice(1)
         kkk=intrec%fraclink(1)
         if(parlist.eq.2) then
            lint(2,nint)=disfra%splink(kkk)
         else
            lint(2,nint)=phlista(lokph)%constitlist(kkk)
         endif
         proprec=>intrec%propointer
         ptyloop2: do while(associated(proprec))
!            typty=proprec%proptype
            ij=proprec%proptype
            if(ij.ge.100) then
               typty=ij/100
               typspec=mod(ij,100)
            else
               typty=ij
            endif
!            typspec=proprec%proptype
!            if(typspec.gt.100) then
!               typty=typspec/100
!               typspec=mod(typty,100)
!            else
!               typty=typspec
!            endif
            if(typty.gt.0 .and. typty.le.ndefprop) then
               prop=propid(typty)%symbol
               if(parlist.eq.2) then
! disordered interaction parameter
                  kk=len_trim(prop)+1
                  prop(kk:kk)='D'
               endif
               if(btest(propid(typty)%status,IDELSUFFIX)) then
! property like ZZ&<element>(phase,constituent array)
! the element index should be in typsepc
                  iel=typspec
                  if(iel.ge.0 .and. iel.le.noofel) then
                     prop=prop(1:len_trim(prop))//'&'&
                          //ellista(elements(iel))%symbol
                  else
!                          write(*,*)'lpd 7: ',iel,typty
                     gx%bmperr=4082; goto 1000
                  endif
               elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
! property like mobility MQ&<constiutent#sublatt>(phase,constituent array)
! the suffix is a constituent
                  iel=typspec
                  if(iel.gt.0 .and. iel.le.phlista(lokph)%tnooffr) then
                     if(parlist.eq.2) then
! we must consider parlist, take disordered constituent list
! we have no current equilibrium record but can use firsteq!!
!                        write(*,*)'3C: typspec: 3 ',typty,iel,prop(1:10)
                        linkcon=disfrap%splink(iel)
!                        write(*,*)'3C: typspec: 4 ',typty,linkcon,prop(1:10)
                        ll=1
                        if(iel.gt.disfrap%nooffr(1)) ll=2
                        prop=prop(1:len_trim(prop))//'&'&
                             //splista(linkcon)%symbol
                        goto 220
                     else
                        linkcon=phlista(lokph)%constitlist(iel)
                        if(linkcon.le.0) then
!                           write(*,*)'Illegal use of wildcard 2'
                           gx%bmperr=4286; goto 1000
                        endif
                        prop=prop(1:len_trim(prop))//'&'&
                             //splista(linkcon)%symbol
! also add the sublattice number ...
                        ncsum=0
                        do ll=1,phlista(lokph)%noofsubl
                           ncsum=ncsum+phlista(lokph)%nooffr(ll)
                           if(iel.le.ncsum) goto 220
                        enddo
                     endif
! there cannot be any errors here ....
!                     write(*,*)'Never never error 2'
                     gx%bmperr=4288; goto 1000
220                  continue
!                     write(*,*)'property 2: ',prop(1:10),ll
! add sublattice index only if not unity
                     if(ll.gt.1) then
                        prop=prop(1:len_trim(prop))//'#'//char(ll+ichar('0'))
                     endif
                  else
!                          write(*,*)'lpd 7: ',iel,typty
                     gx%bmperr=4082; goto 1000
                  endif
               endif
            else
! unknown property ...
               write(*,*)'unknown property type yy: ',typty
               prop='ZZ'
            endif
! if disordered fraction set add D, already set above ??!!
!         if(parlist.eq.2) then
!            prop=prop(1:len_trim(prop))//'D'
!         endif
! note changes here must be repeated for endmember parameters above
            degree: do jdeg=0,proprec%degree
               if(proprec%degreelink(jdeg).eq.0) then
!                  write(*,*)'Ignoring function link'
                  cycle degree
               endif
               call encode_constarr(text,nsl,endm,nint,lint,jdeg)
               write(funexpr,300)prop(1:len_trim(prop)), &
                    phname(1:len_trim(phname)),text(1:len_trim(text))
300            format(A,'(',A,',',A,') ')
               ip=len_trim(funexpr)+1
               call list_tpfun(proprec%degreelink(jdeg),1,funexpr(ip:))
               ip=len_trim(funexpr)
               funexpr(ip+1:)=' '//proprec%reference
               ip=len_trim(funexpr)
               call wrice2(lut,4,12,78,1,funexpr(1:ip))
            enddo degree
            proprec=>proprec%nextpr
         enddo ptyloop2
! list temporarily the number of permutations
         if(intrec%noofip(1).gt.1 .or. intrec%noofip(2).gt.1) then
            if(nint.eq.1) then
               nz=intrec%noofip(2)
            else
               nz=size(intrec%sublattice)
               lqq=intrec%noofip(size(intrec%noofip))
               if(lqq.ne.nz) then
                  write(*,*)'Not same: ',intrec%antalint,nz,lqq
               endif
!               write(*,301)nz,intrec%noofip
301            format('noofip: ',10i3)
!               nz=intrec%noofip(intrec%noofip(1)+2)
            endif
            iqnext=0
            iqhigh=0
            if(associated(intrec%highlink)) then
               iqhigh=intrec%highlink%antalint
            endif
            if(associated(intrec%nextlink)) then
               iqnext=intrec%nextlink%antalint
            endif
            write(*,302)intrec%antalint,nz,nint,iqhigh,iqnext
302         format('@$ Interaction, permutations, level, high, next: ',5i5)
         endif
         intrec=>intrec%highlink
         empty: do while(.not.associated(intrec))
            if(nint.gt.0) then
! restore pointers in same clumsy way
               intrec=>intrecstack(nint)%p1
               intrec=>intrec%nextlink
!               write(*,*)'poping a pointer from intrecstack',ninit
               nint=nint-1
            else
               exit intlist2
            endif
         enddo empty
      enddo intlist2
      endmemrec=>endmemrec%nextem
   enddo intlist1
! check if there are other fraction lists
!   parlist=parlist+1, hm parlist can only be 1 or 2
!   write(*,*)'checking for disordered parameters'
   if(parlist.eq.1 .and. associated(phlista(lokph)%disordered)) then
      subref=.TRUE.
!      lokcs=phlista(lokph)%cslink
      lokcs=phlista(lokph)%linktocs(ics)
! does this make a copy?  Maybe it should be a pointer
      disfra=firsteq%phase_varres(lokcs)%disfra
      write(lut,810)disfra%fsites
810    format('Disordered fraction set parameters, factor: ',F10.4,2x,10('-'))
      nsl=disfra%ndd
      parlist=2
      if(ocv()) write(*,*)'Jump back to list disordered',nsl,parlist
      goto 100
   endif
1000 continue
   return
 END subroutine list_phase_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine list_phase_data2(iph,ftyp,CHTD,lut)
! list parameter data for a phase on unit lut in ftyp format, ftyp=2 is TDB
   implicit none
   integer iph,lut,ftyp
   character CHTD*1
!\end{verbatim}
   integer typty,parlist,typspec,lokph,nsl,nk,ip,ll,jnr,ics,lokcs,isp
   integer nint,ideg,ij,kk,iel,ncsum,kkx,kkk,jdeg,iqnext,iqhigh,lqq,nz,ik
   integer intpq,linkcon
   character text*1024,phname*24,prop*32,funexpr*1024
   character special*8
   integer, dimension(2,3) :: lint
   integer, dimension(maxsubl) :: endm,ilist
   logical subref,noelin1
   type(gtp_fraction_set), pointer :: disfrap
! a smart way to have an array of pointers
   TYPE intrecarray 
      type(gtp_interaction), pointer :: p1
   end TYPE intrecarray
   type(intrecarray), dimension(20) :: intrecstack
   type(gtp_property), pointer :: proprec
   type(gtp_interaction), pointer :: intrec
   type(gtp_endmember), pointer :: endmemrec
   TYPE(gtp_fraction_set) :: disfra
   TYPE(gtp_phase_add), pointer :: addrec
! an empty line first
   write(lut,*)
! for type definitions
   if(iph.lt.0 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   elseif(noofel.eq.0) then
! this needed as there is a reference phase with iph=0 when there are elements
      goto 1000
   endif
!   write(*,*)'lpd 1:',iph,phases(iph)
   if(iph.gt.0) then
      lokph=phases(iph)
   else
      lokph=0
   endif
   ics=1
   phname=phlista(lokph)%name
   nsl=phlista(lokph)%noofsubl
   special=' '
   special(1:1)='%'
   isp=1
! indicate some status bit specially, not useful for TDB files ...
!   if(btest(phlista(lokph)%status1,PHFORD)) then
!      special(2:2)='F'
!      isp=2
!   elseif(btest(phlista(lokph)%status1,PHBORD)) then
!      special(2:2)='B'
!      isp=2
!   elseif(btest(phlista(lokph)%status1,PHSORD)) then
!      special(2:2)='S'
!      isp=2
!   elseif(btest(phlista(lokph)%status1,PHIONLIQ)) then
!      special(2:2)='I'
!      isp=2
!   endif
! here isp can be 1 or 2
!   if(btest(phlista(lokph)%status1,PHMFS)) then
!      isp=isp+1
!      special(isp:isp)='D'
!   endif
   if(btest(phlista(lokph)%status1,PHIONLIQ)) then
      lokcs=len_trim(phname)+1
      phname(lokcs:)=':Y'
   elseif(btest(phlista(lokph)%status1,PHGAS)) then
      phname='GAS:G'
   elseif(btest(phlista(lokph)%status1,PHLIQ)) then
      phname='LIQUID:L'
   endif
   if(btest(phlista(lokph)%status1,PHMFS)) then
!      write(*,*)'3C typedef character 1 ',ichar(CHTD),' "',chtd,'"'
      CHTD=char(ichar(CHTD)+1)
      isp=isp+1
      special(isp:isp)=CHTD
      if(.not.btest(globaldata%status,GSSILENT)) then
         write(kou,53)
53       format(' *** Warning: disordered fraction sets need manual editing!')
      endif
      write(lut,55)CHTD,phname(1:len_trim(phname)),phname(1:len_trim(phname))
55    format('$ *** Warning: disordered fraction sets need manual editing!'/&
           ' TYPE_DEFINITION ',a,' GES A_P_D ',a,' DIS_PART DIS_',a,' !')
   endif
! additions
   addrec=>phlista(lokph)%additions
   lastadd: do while(associated(addrec))
! no need to increment CHTD except for magnetism
      if(addrec%type.eq.1) then
!         write(*,*)'3C typedef character 2 ',ichar(CHTD),' "',chtd,'"'
         CHTD=char(ichar(CHTD)+1)
         isp=isp+1
         special(isp:isp)=CHTD
      endif
      call list_addition(lut,CHTD,phname,ftyp,addrec)
      addrec=>addrec%nextadd
   enddo lastadd
60 continue
! This subroutine is independent of current equilibrium, use firsteq
!   write(lut,10)phname,phlista(lokph)%status1,special,&
!        nsl,(phlista(lokph)%sites(ll),ll=1,nsl)
!   write(*,*)'3C phase: ',phname,special
   lokcs=phlista(lokph)%linktocs(ics)
   write(lut,10,advance='no')phname(1:len_trim(phname)),special(1:isp),&
        nsl,(firsteq%phase_varres(lokcs)%sites(ll),ll=1,nsl)
10  format(' PHASE ',A,1x,a,1x,I2,10(1x,F7.3))
   write(lut,11)
11 format('!')
   nk=0
   text='   CONSTITUENT '//phname(1:len_trim(phname))//' :'
   ip=len_trim(text)+1
   sublatloop: do ll=1,nsl
      constloop: do ik=1,phlista(lokph)%nooffr(ll)
         nk=nk+1
         jnr=phlista(lokph)%constitlist(nk)
         if(jnr.gt.0) then
            text(ip:)=splista(jnr)%symbol
         else
            text(ip:)='*'
         endif
         ip=len_trim(text)+1
!         text(ip:ip)=','
         text(ip:ip)=' '
         ip=ip+1
      enddo constloop
      text(ip-1:ip)=':  '
      ip=ip+1
   enddo sublatloop
   text(ip-2:)=':!'
   call wrice2(lut,2,4,78,-1,text)
!    write(lut,17)text(1:ip)
!17  format(A)
! remove any :Y, :L or :G
   ip=index(phname,':')
   if(ip.gt.0) phname(ip:)=' '
! parameters for end members using site fractions
   if(btest(phlista(lokph)%status1,PHMFS)) then
      subref=.FALSE.
   else
      subref=.TRUE.
   endif
   parlist=1
!--------------------------------------------------
! return here to list disordered parameters
100 continue
! parlist changed below for disordered fraction set
   if(parlist.eq.1) then
      endmemrec=>phlista(lokph)%ordered
   else
      if(ocv()) write(*,*)'Listing disordred parameters ',nsl
      endmemrec=>phlista(lokph)%disordered
      disfrap=>firsteq%phase_varres(lokcs)%disfra
   endif
   endmemberlist: do while(associated(endmemrec))
      do ll=1,nsl
!         ilist(ll)=emlista(lokem)%fraclinks(ll,1)
         ilist(ll)=endmemrec%fraclinks(ll,1)
         if(ilist(ll).gt.0) then
            if(parlist.eq.2) then
! what is disfra here??!!
               endm(ll)=disfra%splink(ilist(ll))
            else
               endm(ll)=phlista(lokph)%constitlist(ilist(ll))
            endif
         else
! wildcard, write '*'
            endm(ll)=-99
         endif
      enddo
      nint=0
      ideg=0
      call encode_constarr(text,nsl,endm,nint,lint,ideg)
      if(gx%bmperr.ne.0) goto 1000
      proprec=>endmemrec%propointer
      ptyloop: do while(associated(proprec))
         ij=proprec%proptype
         if(ij.ge.100) then
            typty=ij/100
            typspec=mod(ij,100)
         else
            typty=ij
         endif
         if(typty.gt.0 .and. typty.le.ndefprop) then
            prop=propid(typty)%symbol
            if(parlist.eq.2) then
! disordered endmember parameter
               kk=len_trim(prop)+1
               prop(kk:kk)='D'
            endif
            if(btest(propid(typty)%status,IDELSUFFIX)) then
! property like ZZ&<element>(phase,constituent array)
! the element index should be in typsepc
               iel=typspec
               if(iel.ge.0 .and. iel.le.noofel) then
!                  prop=propid(typty)%symbol
                  prop=prop(1:len_trim(prop))//'&'&
                       //ellista(elements(iel))%symbol
               else
                  gx%bmperr=4082; goto 1000
               endif
            elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
! property like mobility, MQ&<constituent#sublat>(phase,constituent array)
! the suffix is a constituent
               iel=typspec
               if(iel.gt.0 .and. iel.le.phlista(lokph)%tnooffr) then
                  if(parlist.eq.2) then
! we must consider parlist, take disordered constituent list
! we have no current equilibrium record but can use firsteq!!
!                     lokcs=phlista(lokph)%linktocs(1)
!                     write(*,*)'3C: endmember typspec 1: ',iel
                     linkcon=disfrap%splink(iel)
!                     write(*,*)'3C: endmember typspec 2: ',linkcon
                     ll=1
                     if(linkcon.gt.disfrap%nooffr(1)) ll=2
                     prop=prop(1:len_trim(prop))//'&'&
                          //splista(linkcon)%symbol
                     goto 120
                  else
                     linkcon=phlista(lokph)%constitlist(iel)
                     if(linkcon.le.0) then
!                        write(*,*)'Illegal use of wildcard 1'
                        gx%bmperr=4286; goto 1000
                     endif
                     prop=prop(1:len_trim(prop))//'&'&
                          //splista(linkcon)%symbol
! also add the sublattice number ...
                     ncsum=0
                     do ll=1,phlista(lokph)%noofsubl
                        ncsum=ncsum+phlista(lokph)%nooffr(ll)
                        if(iel.le.ncsum) goto 120
                     enddo
                  endif
! error if sublattice not found
                  write(kou,*)'Error in constituent depended parameter id'
                  gx%bmperr=4287; goto 1000
! jump here to append sublattice
120               continue
!                  write(*,*)'property 1: ',prop(1:10),ll
                  prop=prop(1:len_trim(prop))//'#'//char(ll+ichar('0'))
               else
                  write(kou,*)'lpd 7B: ',iel,typty
                  gx%bmperr=4082; goto 1000
               endif
            endif
         else
! unknown property ...
            write(*,*)'unknown property type xx: ',ij,typty,typspec
            prop='ZZ'
         endif
! if disordered fraction set add D, already done above
!         if(parlist.eq.2) then
!            prop=prop(1:len_trim(prop))//'D'
!         endif
! note changes here must be repeated for interaction parameters below
         write(funexpr,200)prop(1:len_trim(prop)),&
              phname(1:len_trim(phname)),text(1:len_trim(text))
200      format('   PARAMETER ',A,'(',A,',',A,') ')
         ip=len_trim(funexpr)+1
!-------------------------------- this is not done for TDB files
! subtract reference states
!         if(subref .and. typty.eq.1) then
!            call subrefstates(funexpr,ip,lokph,parlist,endm,noelin1)
!            if(noelin1) then
! this can happen for ionic liquids with just neutrals in sublattice 2
! replace the constituent in sublattice 1 with "*" !!!
!               write(*,*)'before: ',funexpr(1:ip)
!               kk=index(funexpr,',')
!               ik=index(funexpr,':')
!               funexpr(kk+1:)='*'//funexpr(ik:)
!               ip=len_trim(funexpr)+2
!               write(*,*)'after:  ',funexpr(1:ip)
!            endif
!         endif
! this writes the expression, problem if function is zero
         call list_tpfun(proprec%degreelink(0),1,funexpr(ip:))
! remove = sign
         ip=index(funexpr,'=')
         funexpr(ip:ip)=' '
         ip=len_trim(funexpr)
         funexpr(ip+1:)=' '//proprec%reference
         ip=len_trim(funexpr)
         funexpr(ip+1:)=' !'
! nice output over several lines if needed with indentation 12 spaces
         call wrice2(lut,2,12,78,1,funexpr(1:ip+2))
         proprec=>proprec%nextpr
      enddo ptyloop
      if(endmemrec%noofpermut.gt.1) then
         intpq=0
         if(associated(endmemrec%intpointer)) then
            intpq=endmemrec%intpointer%antalint
         endif
!         write(kou,207)endmemrec%antalem,endmemrec%noofpermut,intpq
207      format('@$ Endmember, permutations, interaction: ',3i5)
      endif
      endmemrec=>endmemrec%nextem
   enddo endmemberlist
!-----------------------------------------------------------------------
! parameters for interactions using site fractions
   if(parlist.eq.1) then
      endmemrec=>phlista(lokph)%ordered
   else
      endmemrec=>phlista(lokph)%disordered
   endif
   intlist1: do while(associated(endmemrec))
      intrec=>endmemrec%intpointer
      if(associated(intrec)) then
!         write(*,*)'intlist 1B: ',intrec%status
         do ll=1,nsl
            kkx=endmemrec%fraclinks(ll,1)
            if(kkx.eq.-99) then
! wildcard
               endm(ll)=-99
            elseif(parlist.eq.2) then
               endm(ll)=disfra%splink(kkx)
            else
               endm(ll)=phlista(lokph)%constitlist(kkx)
            endif
         enddo
      endif
      nint=0
      intlist2: do while(associated(intrec))
         nint=nint+1
         intrecstack(nint)%p1=>intrec
         lint(1,nint)=intrec%sublattice(1)
         kkk=intrec%fraclink(1)
         if(parlist.eq.2) then
            lint(2,nint)=disfra%splink(kkk)
         else
            lint(2,nint)=phlista(lokph)%constitlist(kkk)
         endif
         proprec=>intrec%propointer
         ptyloop2: do while(associated(proprec))
!            typty=proprec%proptype
            ij=proprec%proptype
            if(ij.ge.100) then
               typty=ij/100
               typspec=mod(ij,100)
            else
               typty=ij
            endif
!            typspec=proprec%proptype
!            if(typspec.gt.100) then
!               typty=typspec/100
!               typspec=mod(typty,100)
!            else
!               typty=typspec
!            endif
            if(typty.gt.0 .and. typty.le.ndefprop) then
               prop=propid(typty)%symbol
               if(parlist.eq.2) then
! disordered interaction parameter
                  kk=len_trim(prop)+1
                  prop(kk:kk)='D'
               endif
               if(btest(propid(typty)%status,IDELSUFFIX)) then
! property like ZZ&<element>(phase,constituent array)
! the element index should be in typsepc
                  iel=typspec
                  if(iel.ge.0 .and. iel.le.noofel) then
                     prop=prop(1:len_trim(prop))//'&'&
                          //ellista(elements(iel))%symbol
                  else
!                          write(*,*)'lpd 7: ',iel,typty
                     gx%bmperr=4082; goto 1000
                  endif
               elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
! property like mobility MQ&<constiutent#sublatt>(phase,constituent array)
! the suffix is a constituent
                  iel=typspec
                  if(iel.gt.0 .and. iel.le.phlista(lokph)%tnooffr) then
                     if(parlist.eq.2) then
! we must consider parlist, take disordered constituent list
! we have no current equilibrium record but can use firsteq!!
!                        write(*,*)'3C: typspec: 3 ',typty,iel,prop(1:10)
                        linkcon=disfrap%splink(iel)
!                        write(*,*)'3C: typspec: 4 ',typty,linkcon,prop(1:10)
                        ll=1
                        if(iel.gt.disfrap%nooffr(1)) ll=2
                        prop=prop(1:len_trim(prop))//'&'&
                             //splista(linkcon)%symbol
                        goto 220
                     else
                        linkcon=phlista(lokph)%constitlist(iel)
                        if(linkcon.le.0) then
!                           write(*,*)'Illegal use of wildcard 2'
                           gx%bmperr=4286; goto 1000
                        endif
                        prop=prop(1:len_trim(prop))//'&'&
                             //splista(linkcon)%symbol
! also add the sublattice number ...
                        ncsum=0
                        do ll=1,phlista(lokph)%noofsubl
                           ncsum=ncsum+phlista(lokph)%nooffr(ll)
                           if(iel.le.ncsum) goto 220
                        enddo
                     endif
! there cannot be any errors here ....
!                     write(*,*)'Never never error 2'
                     gx%bmperr=4288; goto 1000
220                  continue
!                     write(*,*)'property 2: ',prop(1:10),ll
                     prop=prop(1:len_trim(prop))//'#'//char(ll+ichar('0'))
                  else
!                          write(*,*)'lpd 7: ',iel,typty
                     gx%bmperr=4082; goto 1000
                  endif
               endif
            else
! unknown property ...
               write(*,*)'unknown property type yy: ',typty
               prop='ZZ'
            endif
! if disordered fraction set add D, already set above ??!!
!         if(parlist.eq.2) then
!            prop=prop(1:len_trim(prop))//'D'
!         endif
! note changes here must be repeated for endmember parameters above
            degree: do jdeg=0,proprec%degree
               if(proprec%degreelink(jdeg).eq.0) then
!                  write(*,*)'Ignoring function link'
                  cycle degree
               endif
               call encode_constarr(text,nsl,endm,nint,lint,jdeg)
               write(funexpr,300)prop(1:len_trim(prop)), &
                    phname(1:len_trim(phname)),text(1:len_trim(text))
300            format('PARAMETER ',A,'(',A,',',A,') ')
               ip=len_trim(funexpr)+1
               call list_tpfun(proprec%degreelink(jdeg),1,funexpr(ip:))
! remove = sign
               ip=index(funexpr,'=')
               funexpr(ip:ip)=' '
               ip=len_trim(funexpr)
               funexpr(ip+1:)=' '//proprec%reference
               ip=len_trim(funexpr)
               funexpr(ip+1:)=' !'
               call wrice2(lut,4,12,78,1,funexpr(1:ip+2))
            enddo degree
            proprec=>proprec%nextpr
         enddo ptyloop2
! list temporarily the number of permutations
         if(intrec%noofip(1).gt.1 .or. intrec%noofip(2).gt.1) then
            if(nint.eq.1) then
               nz=intrec%noofip(2)
            else
               nz=size(intrec%sublattice)
               lqq=intrec%noofip(size(intrec%noofip))
               if(lqq.ne.nz) then
                  write(*,*)'Not same: ',intrec%antalint,nz,lqq
               endif
!               write(*,301)nz,intrec%noofip
301            format('noofip: ',10i3)
!               nz=intrec%noofip(intrec%noofip(1)+2)
            endif
            iqnext=0
            iqhigh=0
            if(associated(intrec%highlink)) then
               iqhigh=intrec%highlink%antalint
            endif
            if(associated(intrec%nextlink)) then
               iqnext=intrec%nextlink%antalint
            endif
            write(*,302)intrec%antalint,nz,nint,iqhigh,iqnext
302         format('@$ Interaction, permutations, level, high, next: ',5i5)
         endif
         intrec=>intrec%highlink
         empty: do while(.not.associated(intrec))
            if(nint.gt.0) then
! restore pointers in same clumsy way
               intrec=>intrecstack(nint)%p1
               intrec=>intrec%nextlink
!               write(*,*)'poping a pointer from intrecstack',ninit
               nint=nint-1
            else
               exit intlist2
            endif
         enddo empty
      enddo intlist2
      endmemrec=>endmemrec%nextem
   enddo intlist1
! check if there are other fraction lists
!   parlist=parlist+1, hm parlist can only be 1 or 2
!   write(*,*)'checking for disordered parameters'
   if(parlist.eq.1 .and. associated(phlista(lokph)%disordered)) then
      write(lut,810)
810    format('$ Disordered fraction parameters:',20('-'))
      subref=.TRUE.
!      lokcs=phlista(lokph)%cslink
      lokcs=phlista(lokph)%linktocs(ics)
! does this make a copy?  Maybe it should be a pointer
      disfra=firsteq%phase_varres(lokcs)%disfra
      nsl=disfra%ndd
      parlist=2
      if(ocv()) write(*,*)'Jump back to list disordered',nsl,parlist
      goto 100
   endif
1000 continue
   return
 END subroutine list_phase_data2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine subrefstates(funexpr,jp,lokph,parlist,endm,noelin1)
! list a sum of reference states for a G parameter
! like "-H298(BCC_A2,FE)-3*H298(GRAPITE,C)"
   implicit none
   integer jp,lokph,parlist,endm(*)
   character funexpr*(*)
   logical noelin1
!\end{verbatim}
! special care for ionic liquid as sites varies ...
   character text*80,els*2
   integer element(maxel),lokel
   double precision coef(maxel),xx,pqval(2)
!   TYPE(gtp_fraction_set) :: disfra
   TYPE(gtp_fraction_set), pointer :: disfra
   integer nsl,lokcs,ie,ll,jsp,nrel,ik,je,more,is,ip
!
   noelin1=.FALSE.
   lokcs=phlista(lokph)%linktocs(1)
   if(btest(phlista(lokph)%status1,PHIONLIQ)) goto 210
   if(parlist.eq.1) then
      nsl=phlista(lokph)%noofsubl
   else
! for disordered fraction set always use 1 as factor ??
! How about bcc with C? the second sublattice should count ...
! CONCLUSION: If disordered fraction set has 2 sublattices calculate
! should disfra be a pointer?? It seems to work like this ....
!      disfra=firsteq%phase_varres(lokcs)%disfra      
      disfra=>firsteq%phase_varres(lokcs)%disfra
      nsl=disfra%ndd
   endif
   ie=0
   sublat: do ll=1,nsl
      jsp=endm(ll)
      if(jsp.gt.0) then
         nrel=splista(jsp)%noofel
         elem: do ik=1,nrel
            do je=1,ie
               if(splista(jsp)%ellinks(ik).eq.element(je)) then
                  if(parlist.eq.1) then
                     coef(je)=coef(je)+&
                     firsteq%phase_varres(lokcs)%sites(ll)*&
                     splista(jsp)%stoichiometry(ik)
!                     phlista(lokph)%sites(ll)*splista(jsp)%stoichiometry(ik)
                  else
                     coef(je)=coef(je)+&
                     disfra%dsites(ll)*splista(jsp)%stoichiometry(ik)
                  endif
                  goto 200
               endif
            enddo
! new element, increment ie and initiate coef
! ignore the element VA with element index 0
            if(splista(jsp)%ellinks(ik).eq.0) goto 200
            ie=ie+1
            element(ie)=splista(jsp)%ellinks(ik)
            if(parlist.eq.1) then
               coef(ie)=&
                    firsteq%phase_varres(lokcs)%sites(ll)*&
                    splista(jsp)%stoichiometry(ik)
!                    phlista(lokph)%sites(ll)*splista(jsp)%stoichiometry(ik)
            else
! if a single disordered sublattice ignore the number of sites !!!
               if(nsl.eq.1) then
                  coef(ie)=splista(jsp)%stoichiometry(ik)
               else
                  coef(ie)=disfra%dsites(ll)*splista(jsp)%stoichiometry(ik)
               endif
            endif
200          continue
         enddo elem
      else
! wildcard, ignore references
         continue
      endif
   enddo sublat
   goto 300
!------------------------------------------------------------
! ionic liquid special, 2 sublattices but sites varies with charges
210 continue
   ie=0
   jsp=endm(1)
   if(jsp.gt.0) then
      pqval(2)=splista(jsp)%charge
   else
      pqval(2)=one
   endif
   jsp=endm(2)
   if(jsp.gt.0) then
      if(btest(splista(jsp)%status,SPVA)) then
         pqval(1)=one
      else
         pqval(1)=-splista(jsp)%charge
         if(pqval(1).eq.zero) then
            noelin1=.TRUE.
            pqval(2)=one
         endif
      endif
   else
!      write(*,*)'Illegal with wildcards in 2nd sublattice'
      gx%bmperr=4262; goto 1000
   endif
   ionsl: do ll=1,2
      jsp=endm(ll)
      if(jsp.lt.0) cycle
      nrel=splista(jsp)%noofel
      ionel: do ik=1,nrel
         do je=1,ie
            if(splista(jsp)%ellinks(ik).eq.element(je)) then
               coef(je)=coef(je)+&
                    pqval(ll)*splista(jsp)%stoichiometry(ik)
               cycle ionel
            endif
         enddo
! new element, increment ie and initiate coef
! ignore the element VA with element index 0
         if(splista(jsp)%ellinks(ik).ne.0) then
            ie=ie+1
            element(ie)=splista(jsp)%ellinks(ik)
            coef(ie)=&
                 pqval(ll)*splista(jsp)%stoichiometry(ik)
            endif
      enddo ionel
   enddo ionsl
!------------------------------------------------------------
! sort the elements
300 continue
   more=0
   do je=1,ie-1
      if(element(je).gt.element(je+1)) then
         is=element(je)
         element(je)=element(je+1)
         element(je+1)=is
         xx=coef(je)
         coef(je)=coef(je+1)
         coef(je+1)=xx
         more=1
      endif
   enddo
   if(more.gt.0) goto 300
! list the elements as -10*H298(SER,element)
!    write(*,*)'subrefstate 2:',ie,(element(i),i=1,ie)
   ip=1
   text=' '
   do je=1,ie
      if(coef(je).ne.one) then
         call wrinum(text,ip,10,6,-coef(je))
         text(ip:ip)='*'
      else
         text(ip:ip)='-'
      endif
      ip=ip+1
      lokel=element(je)
      els=ellista(lokel)%symbol
      if(ellista(lokel)%refstatesymbol.eq.0) then
         text(ip:)='H298(SER,'//els(1:len_trim(els))//')'
      else
         text(ip:)='G(SER,'//els(1:len_trim(els))//')'
      endif
      ip=len_trim(text)+1
   enddo
!    write(*,*)'subrefstate 9: ',ip,text(1:ip)
   funexpr(jp:)=text
   jp=jp+ip
1000 continue
   return
 end subroutine subrefstates

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine encode_stoik(text,ipos,spno)
! generate a stoichiometric formula of species from element list
   implicit none
   integer ipos,spno
   character text*(*)
!\end{verbatim}
   character elnam*2,ltext*60
   integer eli,noelx,iel,isto,jpos,ich,nlen
   double precision stoi,charge
   if(spno.lt.1 .or. spno.gt.noofsp) then
!       write(*,*)'in encode_stoik'
      gx%bmperr=4051
      goto 1000
   endif
   ipos=1
   noelx=splista(spno)%noofel
!  write(6,*)'encode_stoik 1: ',spno,noelx
   loop1: do iel=1,noelx
      eli=splista(spno)%ellinks(iel)
      elnam=ellista(eli)%symbol
!     write(6,*)'encode_stoik 2: ',eli,elnam
      if(elnam(2:2).ne.' ') then
         ltext(ipos:ipos+1)=elnam
         nlen=2
      else
         ltext(ipos:ipos)=elnam
         nlen=1
      endif
      ipos=ipos+nlen
      stoi=splista(spno)%stoichiometry(iel)
      isto=int(stoi)
      if(abs(dble(isto)-stoi).lt.1.0D-3) then
! handle integer stoichiometries nicely
         if(isto.gt.99) then
            write(ltext(ipos:ipos+2),200)isto
200         format(I3)
            ipos=ipos+3
         elseif(isto.gt.9) then
            write(ltext(ipos:ipos+1),205)isto
205         format(I2)
            ipos=ipos+2
         elseif(isto.gt.1) then
            write(ltext(ipos:ipos),210)isto
210          format(i1)
            ipos=ipos+1
!           write(6,*)'encode_stoik 4B: ',ltext(ipos-3:ipos)
         elseif(nlen.eq.1 .and. iel.ne.noelx) then
            ltext(ipos:ipos)='1'
            ipos=ipos+1
         endif
      else
! stoichiometry is a non-integer value
         jpos=ipos
         call wrinum(ltext,ipos,8,0,stoi)
         if(buperr.ne.0) then
            gx%bmperr=buperr; goto 1000
         endif
! remove trailing zeroes
300       continue
         if(ltext(ipos:ipos).eq.'0') then
            ipos=ipos-1; goto 300
         endif
      endif
   enddo loop1
   charge=splista(spno)%charge
   ich=int(charge)
!  write(6,*)'encode_stoik 5: ',ich,charge
   if(ich.lt.zero) then
! limit output to integer charges <10
      ltext(ipos:ipos+3)='/-'//char(ichar('0')-ich)
      ipos=ipos+3
   elseif(charge.gt.zero) then
      ltext(ipos:ipos+3)='/+'//char(ichar('0')+ich)
      ipos=ipos+3
   endif
   text=ltext
   ipos=ipos-1
!  write(6,*)'encode_stoik 6: ',ipos,ltext(1:ipos)
1000 continue
   return
 END subroutine encode_stoik

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine decode_stoik(name,noelx,elsyms,stoik)
! decode a species stoichiometry in name to element index and stoichiometry
! all in upper case
   implicit none
   character name*(*),elsyms(*)*2
   double precision stoik(*)
   integer noelx
!\end{verbatim}
   character lname*72,ch2*2
   double precision xx
   integer ip,jp
   lname=name
   call capson(lname)
   noelx=0
   ip=1
! expect element symbol
   if(eolch(lname,ip)) then
! empty line, expected species stoichiometry
      gx%bmperr=4083; goto 1000
   endif
!   write(*,*)'3C decode_stoik 1: ',name
100 continue
   ch2=lname(ip:ip+1)
!   write(*,*)'Looking for element: ',ip,ch2
   if(ch2(2:2).ge.'A' .and. ch2(2:2).le.'Z') then
      noelx=noelx+1
      elsyms(noelx)=ch2
      ip=ip+2
   elseif(ch2(1:1).ge.'A' .and. ch2(1:1).le.'Z') then
      noelx=noelx+1
      elsyms(noelx)=ch2(1:1)
      ip=ip+1
   elseif(ch2(1:1).eq.'/') then
! electron is always /-, if /+ is given change sign in lname
      noelx=noelx+1
      elsyms(noelx)='/-'
      if(ch2(2:2).eq.'+') then
         lname(ip+1:ip+1)='-'
         ip=ip+1
      elseif(ch2(2:2).eq.'-') then
         ip=ip+2
      else
! do not accept Fe/2 for Fe/+2, always require + or -
!         write(*,*)'Charge must always be given as /+ or /-'
         gx%bmperr=4289; goto 1000
      endif
!      write(*,*)'Found charge: ',ip,noelx,'>',lname(ip:ip+5),'<'
   else
      goto 900
   endif
! an element found, no stoichiometry number means stoik=1
!   write(*,17)'decode_stoik 2: ',ip,ch2,lname(ip:ip+5)
17 format(a,i3,'>',a,'<>',a,'<')
   if(lname(ip:ip).eq.' ') then
      stoik(noelx)=one
   else
      jp=ip
      call getrel(lname,ip,xx)
!      write(*,*)'decode_stoik 3: ',jp,ip,buperr,xx
      if(buperr.eq.0) then
         stoik(noelx)=xx
      else
! Strange error entering stoichiometry U1EUO3.83, ip=4, jp=2 and buperr=1937
! getrel evidently did not find the "1".  Check explictly if lname(jp:jp)
! is a number!
         if(lname(jp:jp).ge.'1' .and. lname(jp:jp).le.'9') then
            stoik(noelx)=dble(ichar(lname(jp:jp))-ichar('0'))
            ip=jp+1
            buperr=0
            goto 100
         else
! accept missing stoichiometry value as 1, it is accepted to write cao as cao
            stoik(noelx)=one
!         buperr=0
! the error can be due to another element follows directly, restore ip an check
!         ip=jp
!         goto 100
         endif
      endif
! in one case of missing stoichiometry ip exceeded length of lname
!      write(*,*)'decode_stoik 4: ',stoik(noelx),buperr
      fraction: if(buperr.eq.0 .and. lname(ip:ip).eq.'/') then
! a stoichiometric factor followed by / without sign will be interpreted
! as a fraction like AL2/3O.  Note AL2/+3 means AL2 with charge +3
         jp=ip+1
         if(.not.(lname(jp:jp).eq.'+' .or. lname(jp:jp).eq.'-')) then
            call getrel(lname,jp,xx)
!            write(*,*)'decode_stoik 5: ',ip,jp,buperr,xx
            if(buperr.eq.0) then
               stoik(noelx)=stoik(noelx)/xx
               ip=jp
            else
               buperr=0
            endif
!         else
!            write(*,*)'Interpret / as charge!'
         endif
      else
!         write(*,*)'3C decode: ',ip,trim(lname)
         buperr=0
      endif fraction
      if(ip.lt.len(lname)) goto 100
   endif
900 continue
   if(noelx.eq.0) then
      write(*,*)'3C error in species stoichiometry: ',trim(name),ip
      gx%bmperr=4084
   endif
!    write(*,19)(stoik(i),i=1,noelx)
!19 format('decode_stoik 5: ',5(1PE12.3))
1000 continue
   return
 end subroutine decode_stoik

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine encode_constarr(constarr,nsl,endm,nint,lint,ideg)
! creates a constituent array
   implicit none
   character constarr*(*)
   integer, dimension(*) :: endm
   integer nsl,nint,ideg
   integer, dimension(2,*) :: lint
!\end{verbatim}
   integer ip,mint,ll,l2
   ip=1
   constarr=' '
   mint=1
!  if(nint.gt.0) then
!     write(*,*)'encode_contarr ',lint(1,1),lint(2,1)
!  endif
   do ll=1,nsl
      if(endm(ll).gt.0) then
         constarr(ip:)=splista(endm(ll))%symbol
      else
         constarr(ip:)='*'
      endif
      ip=len_trim(constarr)
      if(mint.le.nint) then
!        write(*,*)'encode_contarr ',lint(1,1),lint(2,1)
         do l2=mint,nint
            if(lint(1,mint).eq.ll) then
               constarr(ip+1:ip+1)=','
               ip=ip+2
               constarr(ip:)=splista(lint(2,mint))%symbol
               ip=len_trim(constarr)
               mint=mint+1
            endif
         enddo
      endif
      constarr(ip+1:ip+1)=':'
      ip=ip+2
   enddo
   constarr(ip-1:ip-1)=';'
   constarr(ip:ip)=char(ideg+ichar('0'))
   return
 end subroutine encode_constarr

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine decode_constarr(lokph,constarr,nsl,endm,nint,lint,ideg)
! deconde a text string with a constituent array
! a constituent array has <species> separated by , or : and ; before degree
   implicit none
   character constarr*(*)
   integer endm(*),lint(2,*)
   integer nsl,nint,ideg,lokph,lord
!\end{verbatim}
   character const*24,ch1*1
   integer ll,ip,jp,isep,loksp,mord,isp,jsp,nord
   integer constlist(5),klok(5),knr(2)
!
   nint=0; ideg=0; ll=1
   endm(ll)=0
   ip=1
! write(*,*)'decode_constarr 1: ',ip,constarr
   if(eolch(constarr,ip)) then
      gx%bmperr=4061; goto 1000
   endif
   jp=ip-1
!  write(*,*)'decode_constarr 2: ',ip,jp
   loop: do while(.true.)
! find separators between constituents, no spaces allowed
      jp=jp+1
      ch1=biglet(constarr(jp:jp))
!       write(*,*)'decode_constarr 3: ',jp,ch1
      letter: if(ch1.eq.',') then
         isep=1
      elseif(ch1.eq.':') then
         isep=2
      elseif(ch1.eq.';') then
         isep=3
      elseif(ch1.eq.' ') then
         isep=4
      elseif(.not.(ch1.ge.'A' .and. ch1.le.'Z')) then
!          write(*,*)'decode_constarr 3B: ',jp,ip,ch1
         if(jp.gt.ip) then
! accept 0-9 and _ and . and / and + and - 
! after the first character of a constituent
!             write(*,24)'decode constarr 24A: "',ch1
            if(.not.((ch1.ge.'0' .and. ch1.le.'9') .or. &
                 ch1.eq.'_' .or. ch1.eq.'.' .or. &
                 ch1.eq.'/' .or. ch1.eq.'+' .or. ch1.eq.'-')) then
!               write(*,24)'3C: decode constarr 24B: "',ch1
24             format(a,a,'"')
               gx%bmperr=4062; goto 1000
            endif
         elseif(ch1.ne.'*') then
! last possibility: wildcard
!             write(*,24)'decode constarr 24C: "',ch1
            gx%bmperr=4062; goto 1000
         endif
!          write(*,24)'decode constarr 24D: "',ch1
         cycle
      else
         cycle
      endif letter
! we have a species name between ip and jp
      const=constarr(ip:jp-1)
      call find_species_record_exact(const,loksp)
      if(gx%bmperr.ne.0) then
         if(const(1:2).eq.'* ') then
! wildcard, the parameter is independent of the fraction in this sublattice
            loksp=-99; gx%bmperr=0
         else
            goto 1000
         endif
      endif
!       write(*,11)'decode constarr 11: ',ip,jp,loksp,const
!11     format(a,3i4,'"',a,'"')
      place: if(endm(ll).eq.0) then
! first constituent of sublattice ll independent of separator
         endm(ll)=loksp
      else
         lint(1,nint)=ll
         lint(2,nint)=loksp
      endif place
      next: if(isep.eq.1) then
! separator was a , next constituent an interaction
         nint=nint+1
      elseif(isep.eq.2) then
!  separator was a ":" meaning new sublattice
         ll=ll+1
         endm(ll)=0
      elseif(isep.eq.3) then
! this is end of constituent array, followed ba a degree 0-9
         ideg=ichar(constarr(jp+1:jp+1))-ichar('0')
         if(ideg.lt.0 .or. ideg.gt.9) then
! a degree must be between 0 and 9
            gx%bmperr=4063; goto 1000
         endif
         exit loop
      elseif(isep.eq.4) then
         exit loop
      endif next
! beginning of next constituent
      ip=jp+1
   enddo loop
! number of sublattices
   nsl=ll
! make sure the constituents are in alphabetcal order for each sublattice.
!--------------------------------------------------------
! Special order of constituents for ionic liquid ....
   if(btest(phlista(lokph)%status1,PHIONLIQ)) then
      constlist(1)=endm(1)
      if(nsl.ne.2) then
         if(nsl.eq.1) then
! when ionic liquid parameters entered from TDB-TC files parameters
! with just neutrals may have only one sublattice.  Error cleared by
! the readtdb subroutine.  
! BUT we must sort constituents on the sublattice, must be only neutrals ...
! I hope that will be chacked later ...
            do jsp=1,nint
               constlist(1+jsp)=lint(2,jsp)
            enddo
! simple bubble sort of constlist
44          continue
            do jsp=1,nint
               if(constlist(jsp+1).lt.constlist(jsp)) then
                  lord=constlist(jsp)
                  constlist(jsp)=constlist(jsp+1)
                  constlist(jsp+1)=lord
                  goto 44
               endif
            enddo
         endif
         endm(1)=constlist(1)
         do jsp=1,nint
            lint(2,jsp)=constlist(1+jsp)
         enddo
!         if(ocv()) write(*,*)'Ionic liquid has always 2 sublattices'
         gx%bmperr=4255; goto 1000
      endif
      lord=1
      do jsp=1,nint
         if(lint(1,jsp).eq.1) then
            lord=lord+1
            constlist(lord)=lint(2,jsp)
         endif
      enddo
      knr(1)=lord
      lord=lord+1
      constlist(lord)=endm(2)
      do jsp=1,nint
         if(lint(1,jsp).eq.2) then
            lord=lord+1
            constlist(lord)=lint(2,jsp)
         endif
      enddo
      knr(2)=lord-knr(1)
      call sort_ionliqconst(lokph,1,knr,constlist,klok)
      if(gx%bmperr.ne.0) then
         write(*,*)'3C Error return from sort_ionliqconst ',gx%bmperr
         goto 1000
      endif
!      write(*,65)lord,(klok(ll),ll=1,lord)
65    format('3C from sort: ',i5,5x,5i3)
      lord=0
      endm(1)=klok(1)
      do jsp=2,knr(1)
         lord=lord+1
         lint(1,lord)=1
         lint(2,lord)=klok(lord+1)
      enddo
      endm(2)=klok(lord+2)
      do jsp=2,knr(2)
         lord=lord+1
         lint(1,lord)=2
         lint(2,lord)=klok(lord+2)
      enddo
!      write(*,66)endm(1),endm(2),(lint(1,ll),lint(2,ll),ll=1,nint)
66    format('decode: ',2i5,5x,3(2i3,2x))
      goto 1000
   endif
!--------------------------------------------------------
! first the endmember must be in order of the constituents, except wildcard
   order1: do mord=1,nint
      ll=lint(1,mord)
      isp=lint(2,mord)
      jsp=endm(ll)
! we can have isp or jsp or both negative if wildcard, WILDCARD ALWAYS IN ENDM
      if(isp.lt.0 .and. jsp.lt.0) then
! only one wildcard in each sublattice
         gx%bmperr=4032; goto 1000
      elseif(isp.lt.0 .and. jsp.gt.0) then
         endm(ll)=isp
         lint(2,mord)=jsp
      elseif(isp.gt.0 .and. jsp.lt.0) then
         endm(ll)=jsp
         lint(2,mord)=isp
      elseif(splista(isp)%alphaindex.lt.splista(jsp)%alphaindex) then
         endm(ll)=isp
         lint(2,mord)=jsp
      endif
   enddo order1
! then order if there are two interacting in same sublattice
! There are almost never more than 3 constituents interacting in one sublattice
   order2: do mord=1,nint
      ll=lint(1,mord)
      order3: do nord=mord+1,nint
         if(lint(1,nord).eq.ll) then
            isp=lint(2,nord)
            jsp=lint(2,mord)
            if(isp.lt.0 .or. jsp.lt.0) then
               gx%bmperr=4032; goto 1000
            endif
            if(splista(isp)%alphaindex.lt.splista(jsp)%alphaindex) then
               lint(2,mord)=isp
               lint(2,nord)=jsp
            endif
         endif
      enddo order3
   enddo order2
!  write(*,77)(splista(endm(i))%alphaindex,i=1,nsl), &
!       (lint(1,j),lint(2,j),j=1,nint)
!77 format('decode_contarr 7: ',3I3,5x,2i2)
1000 continue
   return
 end subroutine decode_constarr

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_bibliography(bibid,lut)
! list bibliographic references
   implicit none
   integer lut
   character bibid*(*)
!\end{verbatim}
   character longline*2048
   integer ir,jp,nl,ll,maxl
   if(lut.eq.kou) then
      write(lut,10)reffree-1
!   else
!      write(lut,11)reffree-1
   endif
10  format('There are ',i5,' bibliographic references')
11  format('$ There are ',i5,' bibliographic references')
   maxl=0
   do ir=1,reffree-1
      if(bibid(1:1).ne.' ' .and. &
           .not.compare_abbrev(bibid,bibrefs(ir)%reference)) cycle
      longline=bibrefs(ir)%reference
      longline(17:17)="'"
      jp=18
!      nl=size(bibrefs(ir)%refspec)
!      do ll=1,nl
!         longline(jp:)=bibrefs(ir)%refspec(ll)
!         jp=jp+64
!      enddo
! this require Fortran standard 2003/2008      
!      longline(jp:)=bibrefs(ir)%nyrefspec
      ll=bibrefs(ir)%wprefspec(1)
! loadc/storc are WPACK routines to store/load characters in integer arrays
      call loadc(2,bibrefs(ir)%wprefspec,longline(jp:jp+ll-1))
      jp=len_trim(longline)+1
      longline(jp:jp)="'"
      call wrice(lut,0,17,78,longline(1:jp))
      maxl=maxl+1
      if(lut.ne.kou .and. maxl.gt.50) then
! Thermo-Calc limit is 150 lines for each LIST_OF_REFERENCES on a TDB file
         write(lut,17)
17       format(' !'//' ADD_REFERENCES'/'  NUMBER  SOURCE'/" dummy ' '")
         maxl=0
      endif
   enddo
!   write(*,*)'3C refs: ',reffree,maxl
1000 continue
   return
 end subroutine list_bibliography

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_conditions(lut,ceq)
! lists conditions on lut
   implicit none
   integer lut
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character*1024 text
   integer kl
   text=' '
   call get_all_conditions(text,0,ceq)
   if(gx%bmperr.ne.0) goto 1000
   kl=index(text,'CRLF')
   if(kl.gt.1) then
      call wrice2(lut,2,4,78,1,text(1:kl-1))
   endif
   write(lut,50)text(kl+4:len_trim(text))
50 format(a)
1000 continue
   return
 end subroutine list_conditions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine get_one_condition(ip,text,seqz,ceq)
! list the condition with the index seqz into text
! It lists also fix phases and conditions that are not active
   implicit none
   integer ip,seqz
   character text*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer jl,iterm,indx(4)
   TYPE(gtp_condition), pointer :: last,current
   type(gtp_state_variable), pointer :: svrrec
   double precision wone
!
   if(ip.le.0) ip=1
   text(ip:)=' '
   if(.not.associated(ceq%lastcondition)) then
!      write(*,*)'3C No conditions at all'
      gx%bmperr=4143; goto 1000
   endif
   last=>ceq%lastcondition
   current=>last
70 continue
      if(current%seqz.eq.seqz) goto 100
      current=>current%next
      if(.not.associated(current,last)) goto 70
! no condition with this index found
      gx%bmperr=4131; goto 1000
!
100 continue
   iterm=1
! return here for each term if several
150 continue
   do jl=1,4
      indx(jl)=current%indices(jl,iterm)
   enddo
!   write(*,*)'3C g1c: ',indx
   if(abs(current%condcoeff(iterm)-one).gt.1.0D-10) then
      wone=current%condcoeff(iterm)+one
      if(abs(wone).lt.1.0D-10) then
         text(ip:ip)='-'
         ip=ip+1
      else
! not +1 or -1, write number
! if iterm=1 no not write a positive sign
         if(iterm.eq.1) then
            call wrinum(text,ip,8,1,current%condcoeff(iterm))
         else
            call wrinum(text,ip,8,0,current%condcoeff(iterm))
         endif
         text(ip:ip)='*'
         ip=ip+1
      endif
   elseif(iterm.gt.1) then
! must be a + in front of second and later terms even if coeff is +1
      text(ip:ip)='+'
      ip=ip+1
   endif
! why is ceq needed?? BECAUSE COMPONENTS CAN BE DIFFERENT   ... hm?? !! 
!   call encode_state_variable2(text,ip,current%statev,indx,&
!        current%iunit,current%iref,ceq)
   svrrec=>current%statvar(1)
   call encode_state_variable(text,ip,svrrec,ceq)
   if(iterm.lt.current%noofterms) then
      iterm=iterm+1; goto 150
   endif
! write = followed by the value 
   if(text(ip:ip).ne.' ') ip=ip+1
   text(ip:)='='
   ip=ip+1
!   write(*,*)'3C symlink: ',current%symlink1,current%prescribed
   if(current%symlink1.gt.0) then
! the value is a symbol
      text(ip:)=svflista(current%symlink1)%name
      ip=len_trim(text)+1
   else
      call wrinum(text,ip,10,0,current%prescribed)
   endif
1000 continue
   return
 end subroutine get_one_condition

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine list_experiments(lut,ceq)
! list all experiments into text
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
      call get_one_experiment(ip,text,seqz,.TRUE.,ceq)
!      write(*,*)'3C Back from get_one'
      if(gx%bmperr.ne.0) then
! error code for no more experiments or inactive experiment
!         write(*,*)'3C error line 3117: ',gx%bmperr,seqz,text(1:ip)
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
 end subroutine list_experiments

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine get_one_experiment(ip,text,seqz,eval,ceq)
! list the experiment with the index seqz into text
! It lists also experiments that are not active ??
! UNFINISHED current value should be appended
   implicit none
   integer ip,seqz
   character text*(*)
   logical eval
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer jl,iterm,indx(4),symsym
   TYPE(gtp_condition), pointer :: last,current
   type(gtp_state_variable), pointer :: svrrec
   double precision wone,xxx
   character actual_arg*16
!
   if(ip.le.0) ip=1
   text(ip:)=' '
   if(.not.associated(ceq%lastexperiment)) then
!      write(*,*)'3C No experiments'
      gx%bmperr=4249; goto 1000
   endif
   last=>ceq%lastexperiment
   current=>last
!   write(*,*)'3C index of last experiment: ',current%seqz
70 continue
!   write(*,*)'3C experiment number: ',seqz,current%seqz
   if(current%seqz.eq.seqz) goto 100
   current=>current%next
   if(.not.associated(current,last)) goto 70
! no experiment with this index found or it is inactivated
   gx%bmperr=4131; goto 1000
!
100 continue
   if(current%active.eq.1) then
!      write(*,*)'3C Experiment not active '
      gx%bmperr=4218; goto 1000
   endif
   iterm=1
150 continue
!   write(*,*)'3C Testing is symbol or state variable record',&
!        allocated(current%statvar)
   nostv: if(.not.allocated(current%statvar)) then
! an experiment is a symbol!!! Then statvar is not allocated
      symsym=current%statev
!      write(*,*)'3C A symbol, not a state variable for this experiment',symsym
! get the symbol name
      text=svflista(symsym)%name
      ip=len_trim(text)+1
!      text(ip-1:ip-1)='='
!      write(*,*)'3C experiment: ',text(1:ip),ip
   else
!      write(*,*)'3C This experiment has a state variable record',&
!           allocated(current%statvar),allocated(current%indices),iterm
      symsym=0
! these are not needed??
!      do jl=1,4
!         indx(jl)=current%indices(jl,iterm)
!      enddo
!      if(abs(current%condcoeff(iterm)-one).gt.1.0D-10) then
!         wone=current%condcoeff(iterm)+one
!         if(abs(wone).lt.1.0D-10) then
!            text(ip:ip)='-'
!            ip=ip+1
!         else
! not +1 or -1, write number
!            call wrinum(text,ip,8,1,current%condcoeff(iterm))
!            text(ip:ip)='*'
!            ip=ip+1
!         endif
!      elseif(iterm.gt.1) then
! must be a + in front of second and later terms
!         text(ip:ip)='+'
!         ip=ip+1
!      endif
! why is ceq needed?? BECAUSE COMPONENTS CAN BE DIFFERENT   ... hm?? !! 
!   call encode_state_variable2(text,ip,current%statev,indx,&
!        current%iunit,current%iref,ceq)
      svrrec=>current%statvar(1)
      call encode_state_variable(text,ip,svrrec,ceq)
      if(iterm.lt.current%noofterms) then
         iterm=iterm+1; goto 150
      endif
   endif nostv
!   write(*,*)'3C ok here',symsym
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
!   write(*,*)'3C experiment line 2: ',text(1:ip),ip
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
!   write(*,*)'3C experiment line 3: ',text(1:ip),ip
   if(current%symlink2.gt.0) then
! the value is a symbol
      text(ip:)=svflista(current%symlink1)%name
      ip=len_trim(text)+1
   else
!      call wrinum(text,ip,10,0,current%uncertainty)
      call wrinum(text,ip,8,0,current%uncertainty)
   endif
!   write(*,*)'3C ok here 2',symsym,text(1:ip)
!   write(*,*)'3C experiment line 2: ',text(1:ip),ip
   if(current%experimenttype.eq.100) then
      text(ip:ip)='%'
      ip=ip+1
   endif
!   write(*,*)'3C ok here 3',symsym
! if eval TRUE add the current value of the experiment after a $ sign
! TROUBLE GETTING WRONG VALUE HERE WHEN USER DEFINED REFERENCE STATES
   if(.not.eval) then
      text(ip+2:)='$ ?? '
      goto 1000
   endif
   if(symsym.eq.0) then
      call state_variable_val(svrrec,xxx,ceq)
   else
!      write(*,*)'3C ok here 4',symsym
      actual_arg=' '
      xxx=evaluate_svfun_old(symsym,actual_arg,1,ceq)
   endif
   if(gx%bmperr.ne.0) then
! it is maybe a derivative ... 
      write(*,*)'3C we cannot evaluate a derivative here ...',gx%bmperr
! but meq_evaluate_svfun not available here ... it is part of the minimizer
!      gx%bmperr=0
!      xxx=meq_evaluate_svfun(symsym,actual_arg,0,ceq)
!   endif
!   if(gx%bmperr.ne.0) then
      write(*,*)'3C Error evaluating symbol: ',gx%bmperr
      text(ip:)=' $ ?? '
      ip=ip+5
      gx%bmperr=0
   else
!   write(*,*)'3C experimental state variable current value: ',xxx
      text(ip:)=' $'
      ip=ip+3
!      call wrinum(text,ip,12,0,xxx)
      call wrinum(text,ip,8,0,xxx)
!      write(*,*)'3C experiment line 3: ',text(1:ip),ip
   endif
!   write(*,*)'3C ok here 5'
1000 continue
!   write(*,*)'3C experiment line 4: ',text(1:ip),ip,gx%bmperr
   return
 end subroutine get_one_experiment

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine get_all_conditions(text,mode,ceq)
! list all conditions if mode=0, experiments if mode=1, -1 if no numbers
   implicit none
   integer mode
   character text*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   TYPE(gtp_condition), pointer :: last,current,first
   type(gtp_state_variable), pointer :: svrrec
   character phname*32
   integer ntot,nc,ip,iterm,iph,ics,jl
   double precision value,wone
   integer indx(4)
   ntot=0
   text=' '
   if(mode.eq.1) then
! cannot enter experiments yet
      goto 1000
   endif
   if(noofel.eq.0) then
      text='CRLF No elements'
      goto 1000
   endif
   last=>ceq%lastcondition
   if(.not.associated(last)) then
      if(mode.eq.-1) then
         text=' '
      else
! The CRLF indicates CR+LF at output
         write(text,50)noofel+2
50       format('CRLF Degrees of freedom are ',i3)
      endif
      goto 1000
   endif
   current=>last%next
   first=>current
   nc=1
   ip=1
100 continue
! conditions can also be fixed phases !!!
   ntot=ntot+1
   if(current%active.ne.0) then
! if active is nonzero the condition is not active
      goto 200
   endif
   if(mode.ne.-1) then
! no condition numbers for mode=-1
      call wriint(text,ip,nc)
! number the conditions
      text(ip:)=':'
!   ip=ip+2
! No space after :
      ip=ip+1
   endif
   iterm=1
   if(current%statev.lt.0) then
! handle FIX phases
      iph=-current%statev
      ics=current%iref
      call get_phase_name(iph,ics,phname)
      if(gx%bmperr.ne.0) then
         write(*,*)'3C list condition error for phase ',iph,ics
         gx%bmperr=4178; goto 1000
      endif
      text(ip:)='<'//phname
      ip=len_trim(text)+3
      text(ip-2:ip-1)='>='
      value=current%prescribed
      if(value.lt.1.0d-8) then
         value=zero
      endif
      call wrinum(text,ip,4,0,value)
      goto 190
   endif
! return here for each term if several
150 continue
   do jl=1,4
      indx(jl)=current%indices(jl,iterm)
   enddo
!   if(iterm.gt.1) write(*,152)'3C 150: ',iterm,indx,current%condcoeff(iterm)
152 format(a,5i4,1pe12.4)
   if(abs(current%condcoeff(iterm)-one).gt.1.0D-10) then
      wone=current%condcoeff(iterm)+one
      if(abs(wone).lt.1.0D-10) then
         text(ip:ip)='-'
         ip=ip+1
      else
! not +1 or -1, write number
!         write(*,*)'3C list cond: ',current%condcoeff(iterm),one,wone
         if(iterm.eq.1) then
! do not write a + in front of first term
            call wrinum(text,ip,8,0,current%condcoeff(iterm))
         else
            call wrinum(text,ip,8,1,current%condcoeff(iterm))
         endif
         text(ip:ip)='*'
         ip=ip+1
      endif
   elseif(iterm.gt.1) then
! must be a + or - in front of second and later terms
      text(ip:ip)='+'
      ip=ip+1
   endif
! why is ceq needed?? BECAUSE COMPONENTS CAN BE DIFFERENT   ... hm?? !! 
!   write(*,*)'3C encode: ',current%statev,indx
!   call encode_state_variable2(text,ip,current%statev,indx,&
!        current%iunit,current%iref,ceq)
!   svrrec=>current%statvar(1)
   svrrec=>current%statvar(iterm)
   if(svrrec%argtyp.eq.3) then
!      write(*,153)svrrec%argtyp,svrrec%phase,svrrec%compset,svrrec%component
153 format('3C gac 2: ',4i4)
   endif
   call encode_state_variable(text,ip,svrrec,ceq)
   if(iterm.lt.current%noofterms) then
      iterm=iterm+1; goto 150
   endif
! problem with current position ... LNAC(CR) had the last ) overwritten ...
!   write(*,157)ip,text(1:ip)
157 format('3C gc: ',i2,'"',a,'"')
   if(text(ip:ip).ne.' ') ip=ip+1
   text(ip:)='='
   ip=ip+1
   if(current%symlink1.gt.0) then
! the value is a symbol
!      write(*,*)'3C value is a symbol: ',current%symlink1
      text(ip:)=svflista(current%symlink1)%name
      ip=len_trim(text)+1
   else
      call wrinum(text,ip,10,0,current%prescribed)
   endif
190 continue
   if(ip.ge.len(text)) then
      write(*,*)'3C text: "',text,'" ',ip,len(text)
   endif
   text(ip:ip)=', '
   ip=ip+2
   nc=nc+1
200 continue
   current=>current%next
   if(.not.associated(current,first)) goto 100
! there can be non-active conditions only
   if(nc.gt.1) then
! write without the last ,
      text(ip-2:)=' '
!      write(kou,99)text(1:ip-3)
!99    format(a)
   endif
   if(mode.eq.0) then
! the degrees of freedoms   
      write(text(ip:),50)noofel+3-nc
   endif
1000 return
 end subroutine get_all_conditions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine list_defined_properties(lut)
! lists all parameter identifiers allowed
   implicit none
   integer lut
!\end{verbatim}
   character special*32,tdep*1,pdep*1
   integer typty,kk
   write(lut,10)
10 format('Index Ident  T P Specification',23x,' Status Note')
!10 format('Index  Symbol Specification',26x,' Status Note')
   do typty=1,ndefprop
      special=' '
      if(btest(propid(typty)%status,IDELSUFFIX)) then
         special='&<element>'
      elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
         special='&<constituent#sublattice>'
      endif
      kk=len_trim(special)
      if(kk.gt.0) then
         special(kk+1:)=';'
         kk=kk+2
      else
         kk=1
      endif
      tdep='T'
      pdep='P'
      if(btest(propid(typty)%status,IDNOTP)) then
!         special(kk:)='Not T- and P-dependent'
         tdep='-'
         pdep='-'
      elseif(btest(propid(typty)%status,IDONLYP)) then
!         special(kk:)='Not T-dependant'
         tdep='-'
      elseif(btest(propid(typty)%status,IDONLYT)) then
!         special(kk:)='Not P-dependant'
         pdep='-'
      endif
      write(lut,50)typty,propid(typty)%symbol,tdep,pdep,special,&
           propid(typty)%status,propid(typty)%note
50    format(i5,2x,a,2x,a,1x,a,2x,a,2x,z8,1x,a)
   enddo
1000 continue
   return
 end subroutine list_defined_properties

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine find_defined_property(symbol,mode,typty,iph,ics)
! searches the propid list for one with symbol or identifiction typty
! if mode=0 then symbol given, if mode=1 then typty given
! symbol can be TC(BCC), BM(FCC), MQ&FE(HCP) etc, the phase must be 
! given in symbol as otherwise it is impossible to find the consititent!!!
! A constituent may have a sublattice specifier, MQ&FE#3(SIGMA)
   implicit none
   integer mode,typty,iph,ics
   character symbol*(*)
!\end{verbatim}
   character phsym*24,specid*24,nude*4
   integer splink,k1,k2,lattice,lokph,ityp,iel,kk,ll,jj
   integer jtyp
!   write(*,7)'3C fdp 1: ',symbol(1:5),mode,typty,iph,ics
7  format(a,a,5i5)
   if(mode.eq.0) then
! parameter identifier given, can include & # and ( ) like MQ&FE#3(SIGMA)
      lattice=0
      nude=' '
      specid=' '
      k1=index(symbol,'&')
      if(k1.gt.0) then
         nude=symbol(1:k1-1)
         k2=index(symbol,'#')
         if(k2.eq.0) then
            k2=index(symbol,'(')
            if(k2.eq.0) then
!               write(*,*)'3C: Missing phase specifier in property symbol 1'
!               write(*,*)'Error in symbol: ',symbol
               gx%bmperr=4290; goto 1000
            endif
         else
            lattice=ichar(symbol(k2+1:k2+1))-ichar('0')
            if(lattice.le.0 .or. lattice.gt.9) then
!               write(*,*)'3C Sublattice outside range in property symbol'
               gx%bmperr=4290; goto 1000
            endif
         endif
         specid=symbol(k1+1:k2-1)
         call capson(specid)
      endif
! there must be a phase name within ( )
      k1=index(symbol,'(')
      if(k1.gt.0) then
         k2=index(symbol,')')
         if(k2.lt.k1) then
!            write(*,*)'3C Missing phase specifier in property symbol 2'
!            write(*,*)'Symbol: ',symbol
            gx%bmperr=4291; goto 1000
         endif
         phsym=symbol(k1+1:k2-1)
         call find_phase_by_name(phsym,iph,ics)
         if(gx%bmperr.ne.0) goto 1000
         lokph=phases(iph)
         if(nude(1:1).eq.' ') nude=symbol(1:k1-1)
      elseif(mode.ne.0) then
         write(*,*)'3C Missing phase specifier in property symbol 3'
         write(*,*)'Symbol: ',symbol,mode
         gx%bmperr=4291; goto 1000
!      else
! mode=0 means just ignore
!         write(*,*)'3C mode: ',mode,iph,ics
!         goto 1000
      endif
! now nude is the property id, lokph is phase location, specid is element or
! constituent symbol, lattice is sublattice number
! skip index 1 as G is a state variable
      call capson(nude)
!      write(*,*)'3C fdp 2: ',iph,ics,nude
      do ityp=2,ndefprop
!         write(*,*)'3C fdp 3: ',ityp,nude,propid(ityp)%symbol
         if(propid(ityp)%symbol.ne.nude) cycle
         if(btest(propid(ityp)%status,IDELSUFFIX)) then
! element specifier, IBM&CR(BCC) (when we have element specific Bohr magnetons)
!            write(*,*)'3C fdp 4: element: ',specid
            call find_element_by_name(specid,iel)
            if(gx%bmperr.ne.0) goto 1000
            typty=100*ityp+iel
            goto 200
         elseif(btest(propid(ityp)%status,IDCONSUFFIX)) then
! constituent specifier, for example: MQ&FE#3(SIGMA)
!            write(*,*)'3C fdp 5: constituent: ',specid
            kk=0
            do ll=1,phlista(lokph)%noofsubl
               do jj=1,phlista(lokph)%nooffr(ll)
                  kk=kk+1
                  splink=phlista(lokph)%constitlist(kk)
                  if(splink.le.0) then
!                     write(*,*)'3C Illegal use of woildcard 3'
                     gx%bmperr=4286; goto 1000
                  endif
                  if(specid.eq.splista(splink)%symbol .and. &
                       (lattice.eq.0 .or. lattice.eq.ll)) then
                     typty=100*ityp+kk
                     goto 200
                  endif
               enddo
            enddo
         else
! property without specifier like TC(FCC)
            typty=ityp
            goto 200
         endif
      enddo
! if we come here we have not found the constituent or element or property
! it may be OK anyway if this is a call to test if symbol exists ??
!      write(*,*)'3C Illegal property symbol'
      gx%bmperr=4290; goto 1000
! we must return property number, phase location, element
! the value TYPTY stored in property records is "idprop" or
! if IDELSUFFIX set then 100*"idprop"+ellista index of element
! if IDCONSUFFIX set then 100*"idprop"+constituent index
200   continue
   else
! indices given, typty, iph and ics, construct the symbol
! if typty>100 there is also an element or constituent specifier
      lokph=phases(iph)
!      write(*,*)'3C fdp 10: ',typty,iph,ics,lokph
      ityp=typty
      jtyp=-1
      if(ityp.gt.100) then
         ityp=typty/100
         jtyp=typty-100*ityp
      endif
      if(ityp.le.1 .or. ityp.gt.ndefprop) then
!         write(*,*)'3C Property number outside range ',ityp,typty
         gx%bmperr=4292; goto 1000
      endif
      symbol=propid(ityp)%symbol
      if(btest(propid(ityp)%status,IDELSUFFIX)) then
! could one have /- as specifier??? NO !! But maye Va
         if(jtyp.lt.0) then
!            write(*,*)'3C Missing element index in property symbol'
            gx%bmperr=4290; goto 1000
         endif
         if(jtyp.lt.0 .or. jtyp.gt.noofel) then
!            write(*,*)'3C Too high element index in property symbol'
            gx%bmperr=4290; goto 1000
         endif
         symbol=symbol(1:len_trim(symbol))//'&'//ellista(jtyp)%symbol
      elseif(btest(propid(ityp)%status,IDCONSUFFIX)) then
         if(jtyp.lt.0) then
!            write(*,*)'3C Missing constituent index in property symbol'
            gx%bmperr=4290; goto 1000
         endif
         if(iph.le.0 .or. iph.gt.noofph) then
!            write(*,*)'3C Illegal phase location in property symbol'
            gx%bmperr=4290; goto 1000
         endif
         kk=0
         do ll=1,phlista(lokph)%noofsubl
            do jj=1,phlista(lokph)%nooffr(ll)
               kk=kk+1
               if(kk.eq.jtyp) then
                  splink=phlista(lokph)%constitlist(kk)
                  if(splink.le.0) then
!                     write(*,*)'3C Illegal use of woildcard 4'
                     gx%bmperr=4286; goto 1000
                  endif
                  specid=splista(splink)%symbol
                  if(ll.gt.1) then
                     specid=specid(1:len_trim(specid))//&
                          '#'//char(ichar('0')+ll)
                  endif
                  goto 400
               endif
            enddo
         enddo
! we come here is we failed to find the constituent
         write(*,*)'3C Illegal constituent index in property symbol'
         gx%bmperr=4290; goto 1000
400      continue
         symbol=symbol(1:len_trim(symbol))//'&'//specid
      elseif(jtyp.gt.0) then
         write(*,*)'3C This property has no specifier'
         gx%bmperr=4290; goto 1000
      endif
! add the phase
!      write(*,*)'3C fdp 11: ',lokph,ics
      symbol=symbol(1:len_trim(symbol))//'('//phlista(lokph)%name
      if(ics.lt.0 .or. ics.gt.phlista(lokph)%noofcs) then
!         write(*,*)'3C No such composition set'
         gx%bmperr=4072; goto 1000
      endif
      if(ics.gt.1) symbol=symbol(1:len_trim(symbol))//'#'//char(ichar('0')+ics)
      symbol=symbol(1:len_trim(symbol))//')'
!      write(*,*)'3C fdp 12: ',symbol(1:20)
   endif
1000 continue
   return
 end subroutine find_defined_property

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine line_with_phases_withdgm0(line,ceq)
! used in amend lines with stored STEP/MAP results
! enter first 6, two .. and last 2 characters of phase names with abs(dgm)<1-8
! line LIQUID#2, PHYRRO..#2 
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   character line*(*)
!\end{verbatim}
   integer iph,ik,jk,tup
   character name*32
   ik=1
! number of phases is equal to number of phase tuples?? no
   do iph=1,noofph
!         iph=phasetuple(phtupx(isort(jd)))%phaseix
!         ics=phasetuple(phtupx(isort(jd)))%compset
!         call get_phase_compset(iph,ics,lokph,lokcs)
      tup=iph
100   continue
      if(abs(ceq%phase_varres(phasetuple(tup)%lokvares)%dgm).lt.1.0D-9) then
         call get_phasetup_name(tup,name)
         if(gx%bmperr.ne.0) goto 1000
         jk=len_trim(name)
         if(ik+10.gt.len(line)) then
            line(ik:)=' ...'
         elseif(jk.gt.8) then
            line(ik:)=name(1:6)//'..'
            line(ik+8:)=name(jk-1:jk)
            ik=ik+11
         else
            line(ik:)=name
            ik=len_trim(line)+2
         endif
!      else
!         continue
      endif
! find higher composition sets of this phase
      tup=phasetuple(tup)%nextcs
      if(tup.gt.0) goto 100
   enddo
!   write(*,*)'3C phaseline: ',line
1000 continue
   return
 end subroutine line_with_phases_withdgm0

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine list_equilibria_details(mode,teq)
! not used yet ...
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: teq
   integer mode
!\end{verbatim}
   TYPE(gtp_equilibrium_data), pointer :: ceq
!   TYPE(gtp_phase_varres) :: varres
   integer ieq,noofeq,iph
   noofeq=noeq()
   select case(mode)
   case default
      write(*,*)'3C No such mode: ',mode
!--------------------------------------------------
   case(1) ! list equilibria and some general data
      write(*,10)noofeq
10    format('3C Number of equilibria: ',i3)
      do ieq=1,noofeq
         ceq=>eqlista(ieq)
         write(*,11)ceq%eqno,ceq%eqname
11       format('3C Equilibrium ',i3,', ',a)
      enddo
!--------------------------------------------------
   case(100:199) ! list phase varres data for phase mod(mode,100)
      iph=mod(mode,100)
      if(iph.eq.0) then
         write(*,*)'3C all phases'
      else
         write(*,*)'3C phase ',iph
      endif
   end select
1000 continue
   return
 end subroutine list_equilibria_details

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine save_datformat(filename,kod,ceq)
! writes a DAT format file. kod is not used (yet)
   implicit none
   integer kod
   character filename*(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ntpf,last,i1,i2,i3,npows,lut,ip,jp,nstoi,lokph,isp,f1,nphstoi,nphmix
   integer, dimension(:), allocatable :: ncon,phmix,phstoi,estoi
   integer nelectrons,lokcs,nsubl,isubl,mphstoi,k1,lcase,mult(10)
   logical logok,nogas,ionliq,wildcard
   character ch1*1,line*16,text*512,powers*80,model*24,constext*80,elsym*2
   type(gtp_tpfun2dat), dimension(:), allocatable :: tpfc
   type(gtp_endmember), pointer :: endmember
   double precision, allocatable, dimension(:) :: constcomp,constcompiliq
   double precision valency(9),ccc
! we must probably create a stack
!   type(gtp_interaction), dimension(5), pointer :: nint
   type(gtp_phase_varres), pointer :: varres
   type(gtp_interaction), pointer :: hint
   type(gtp_property), pointer :: property
   inquire(file=filename,exist=logok)
   if(logok) then
      line=' '
      last=len(line)
      call gparcd('File exists, overwrite?',line,last,1,ch1,'N',q1help)
      if(ch1.ne.'Y') then
         write(*,*)'Better luck next time!'
         goto 1000
      endif
   endif
   ntpf=freetpfun-1
! allocate coefficient arrays for all TP functions (incl parameters) and 5 extra
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
   open(lut,file=filename,access='sequential',status='unknown')
   write(*,*)'Writing on file: ',trim(filename)
   text='System'
   ip=len_trim(text)+2
   do i1=1,noofel
      text(ip:)=trim(ellista(elements(i1))%symbol)//'-'
      ip=len_trim(text)+1
   enddo
   text(ip-1:)=' generated from TDB file by OC'
   write(lut,100)trim(text)
99 format(a)
100 format(1x,a)
!------------------- we have to sort the phases as SOLGASMIX wants
! and list constituents in gas, mixtures, stoichiometric 
   allocate(ncon(noofph))
   allocate(phmix(noofph))
   allocate(phstoi(noofph))
   allocate(estoi(noofph))
   ncon=0
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
   do i1=1,noofph
      lokph=phasetuple(i1)%lokph
      if(phlista(lokph)%nooffs.gt.1) then
         write(*,105)trim(phlista(lokph)%name)
105      format('Cannot handle phase ',a,' with several fraction sets')
         gx%bmperr=4399; goto 1000
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
            nelectrons=nelectrons+1
            estoi(i1)=-noofel-nelectrons
         endif
! should ncon be the number of endmembers?? YES
! NOTE for ionic liquid with neutrals the DAT format requires that the neutrals
! are repeated for each cation, thus the same equation here!!
         i3=1
         do i2=1,phlista(lokph)%noofsubl
            i3=i3*phlista(lokph)%nooffr(i2)
         enddo
         ncon(i1)=i3
      endif
   enddo
! now can we write the line with overall phase information ... suck
   ip=1
   write(text(ip:),110)noofel+nelectrons
   ip=len_trim(text)+1
! number of mixture phases and for each mixture the number of endmembers
! if nogas is TRUE add 1 with zero endmembers first
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
      if(ncon(i1).gt.0) then
         write(text(ip:),112)ncon(i1)
         ip=len_trim(text)+1
         if(ip.gt.72) then
!            write(lut,100)trim(text)
! Accordung to Ted
            write(lut,99)trim(text)
            ip=1
         endif
      endif
   enddo ph1
! number of stoichiometric phases using i5
   write(text(ip:),112)nphstoi
! NOTE format 100 adds an initial space on the line
!   write(lut,100)trim(text)
! Accordung to Ted
   write(lut,99)trim(text)
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
   do i1=1,noofph
      if(estoi(i1).lt.0) then
         lokph=phasetuple(i1)%lokph
         text(ip:)='e('//trim(phlista(lokph)%name)//')'
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
   allocate(constcomp(noofel+nelectrons))
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
   if(npows.eq.9) then
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
      write(*,*)'Please fix the T powers manually'
      stop
   endif
!-------------------------------------- end of header section
! data for mixtures
! First the endmembers
   mphstoi=1
   phases1: do i1=1,noofph
!      if(i1.gt.1) exit phases
      lokph=phasetuple(i1)%lokph
      if(i1.eq.phstoi(mphstoi)) then
!         write(*,*)'3C skipping stoichiometric ',trim(phlista(lokph)%name)
         mphstoi=mphstoi+1
         cycle phases1
      else
         write(*,*)'3C parameters for mixture ',trim(phlista(lokph)%name)
      endif
      lokcs=phlista(lokph)%linktocs(1)
      varres=>ceq%phase_varres(lokcs)
      nsubl=1
      ionliq=.false.
! phase model
      if(btest(phlista(lokph)%status1,PHFORD)) then
         write(*,141)trim(phlista(lokph)%name)
141      format('Phase ',a,' has FCC permutated parameters, ignored')
         cycle phases1
      elseif(btest(phlista(lokph)%status1,PHBORD)) then
         write(*,142)trim(phlista(lokph)%name)
142      format('Phase ',a,' has BCC permutated parameters, ignored')
         cycle phases1
      elseif(btest(phlista(lokph)%status1,PHMFS)) then
         write(*,143)trim(phlista(lokph)%name)
143      format('Phase ',a,' has disorded fraction sets, ignored')
         cycle phases1
      elseif(btest(phlista(lokph)%status1,PHIONLIQ)) then
         model='SUBI'
         nsubl=2
         ionliq=.true.
! there can just be one ionic liquid ... ??
!         if(.not.allocated(constcompiliq)) then
            allocate(constcompiliq(noofel+nelectrons))
!         endif
      elseif(btest(phlista(lokph)%status1,PHID)) then
         model='IDMX'         
      else
! there are phases with other bits which will not work but they are rarely set
! now for sublattices ...
         nsubl=phlista(lokph)%noofsubl
         if(nsubl.gt.1) then
            model='SUBL'
         else
            model='RKMP'
         endif
      endif
! UNFINISHED: Magnetism? A suffix M ...
      write(lut,200)phlista(lokph)%name,trim(model)
!200   format(1x,a,5x,40('=')/1x,a)
! According to Ted
200   format(a,5x,40('=')/a)
!-------------------- we must repeat this endmember loop for interactions
      endmember=>phlista(lokph)%ordered
! endmember parameters
      do while(associated(endmember))
! we have to generate two lines by extracting the endmember and constituents
         constcomp=zero
         if(ionliq) constcompiliq=zero
         constext=' '
         text=' '
         ip=1
!         jp=1
         valency=zero
         sloop: do isubl=1,nsubl
! this is the loop for the constituents in sublattices
            isp=endmember%fraclinks(isubl,1)
            if(isp.eq.-99) then
! this means wildcard in this sublattice
               write(*,*)'3C Beware! Parameter with wildcard constituent'
               wildcard=.true.
               constext(ip:)='*:'
               ip=ip+2
! Hm we should add stoichiometric factors for all constituents in this subl
! For ionliq this means neutrals on sublattice 2
!>> QUESTION >> the DAT format repeats all neutrals for all cations
!>>>>>>>>>>>>>> with the stoichiometry of the cation !!!
               if(ionliq) valency(1)=one
               cycle sloop
            endif
            isp=phlista(lokph)%constitlist(isp)
            if(btest(splista(isp)%status,SPVA)) then
!               write(*,*)'3C vacancy'
               valency(isubl)=zero
! according to the example I have the stoichiometry should be 1 for (cation:VA)
!               if(ionliq) valency(2)=-valency(1)
               if(ionliq) valency(2)=-one
            else
               valency(isubl)=splista(isp)%charge
               if(abs(valency(isubl)).lt.1.0D-6) valency(isubl)=zero
            endif
            if(ionliq .and. isubl.eq.2) then
               do i3=1,noofel
                  constcomp(i3)=-constcomp(i3)*valency(2)
               enddo
            elseif(estoi(i1).lt.0) then
! charged sublattice phase
               constcomp(-estoi(i1))=constcomp(-estoi(i1))+&
                    valency(isubl)*varres%sites(isubl)
!               write(*,901)'3C e-stoik:',isubl,-estoi(i1),&
!                    valency(isubl),varres%sites(isubl),constcomp(-estoi(i1))
901            format(a,2i3,3F10.2)
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
                  if(isubl.eq.1) then
                     constcomp(i3)=constcomp(i3)+splista(isp)%stoichiometry(i2)
                  else
                     constcomp(i3)=constcomp(i3)+&
                          splista(isp)%stoichiometry(i2)*valency(1)
                  endif
               else
                  constcomp(i3)=constcomp(i3)+splista(isp)%stoichiometry(i2)*&
                       varres%sites(isubl)
!                  if(estoi(i1).lt.0) then
!                     write(*,901)'3C c-stoik:',isubl,i3,&
!                          splista(isp)%stoichiometry(i2),varres%sites(isubl),&
!                          constcomp(i3)
!                  endif
               endif
            enddo
         enddo sloop
! for the parameters follow the property link
         property=>endmember%propointer
         if(associated(property)) then
! some endmembers may not have a property record!!
! this line should be written together with the type of coefficients and ranges
! it may require several lines
            write(text,210)constcomp
!210         format(50F7.1)
! according to Ted
210         format(50F6.1)
! what about several properties??
            f1=property%degreelink(0)
!            write(*,*)'TP function pointer is ',f1
            if(f1.gt.0) then
               if(ionliq .and. wildcard) then
!                  write(lut,207)
!207               format(' *** A neutral should be repeated for each cation',&
!                       ' multipled with its valency')
                  do k1=1,phlista(lokph)%nooffr(1)
! find the charge of all cations
                     isp=phlista(lokph)%constitlist(k1)
                     ccc=splista(isp)%charge
! write the endmember without the first *: and final : with factor ccc
                     write(lut,211)constext(3:ip-2),ccc
!211                  format(1x,a,40x,' * ',F12.2)
! According to Ted
211                  format(a,40x,' * ',F12.2)
! use one of the "extra" coefficient function!
                     jp=ntpf+1
                     call tpmult(f1,jp,ccc,tpfc)
                     call list_tpascoef(lut,text,jp,npows,tpfc)
                  enddo
               else
! write the endmember without the final :
!                  write(lut,100)constext(1:ip-2)
! according to Ted
                  write(lut,99)constext(1:ip-2)
                  call list_tpascoef(lut,text,f1,npows,tpfc)
               endif
            else
               write(*,*)'3 C missing function for endmember property',&
                    constext(3:ip-2)
            endif
         endif
         endmember=>endmember%nextem
      enddo
! ------------------- end of endmembers, constituents and excess parameters ??
      if(model(1:4).eq.'IDMX') cycle phases1
! For sublattice phases we write number of sublattices and sites
      if(model(1:4).eq.'SUBL') then
         write(lut,250)nsubl
         write(lut,260)(ceq%phase_varres(lokcs)%sites(isubl),isubl=1,nsubl)
250      format(1x,i4)
260      format(1x,8F9.5)
      endif
      if(model(1:4).eq.'SUBL' .or. model(1:4).eq.'SUBI') then
! number of constituents in each sublattice
         write(lut,270)(phlista(lokph)%nooffr(isubl),isubl=1,nsubl)
270      format(9i5)
      endif
! For all mixtures we should write the constituents of all sublattices
      i3=0
      do isubl=1,phlista(lokph)%noofsubl
         constext=' '
         ip=1
         do i2=1,phlista(lokph)%nooffr(isubl)
            i3=i3+1
            isp=phlista(lokph)%constitlist(i3)
            jp=ip
            call lower_case_species_name(constext,ip,isp)
            ip=jp+25
            if(ip.ge.75) then
               write(lut,100)trim(constext)
               constext=' '
               ip=1
            endif
         enddo
         ip=len_trim(constext)
         if(ip.gt.1) then
            write(lut,100)trim(constext)
         endif
      enddo
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
      if(phlista(lokph)%noofsubl.gt.1) then
! A very strange output of integers representing endmembers
         jp=1
         mult=1
         do isubl=phlista(lokph)%noofsubl,1,-1
            mult(isubl)=jp
            jp=jp*phlista(lokph)%nooffr(isubl)
         enddo
!         write(*,278)'3C mult2: ',jp,(mult(ip),ip=1,phlista(lokph)%noofsubl)
278      format(a,10i4)
         do isubl=1,phlista(lokph)%noofsubl
            text=' '
            ip=3
            k1=0
            i2=0
290         continue
               k1=k1+1
               i3=0
292            continue
                  call wriint(text,ip,k1)
                  ip=ip+3
                  i2=i2+1
                  i3=i3+1
                  if(i3.lt.mult(isubl)) goto 292
               if(k1.gt.phlista(lokph)%nooffr(isubl)) k1=0
               if(k1.eq.phlista(lokph)%nooffr(isubl) .and. isubl.gt.1) k1=0
            if(i2.lt.jp) goto 290
!            if(i2.lt.phlista(lokph)%nooffr(isubl)) goto 290
! According to Markus Piro one should have 19 values per line, 18*4+3=75
            isp=1
            do while(len_trim(text(isp:))-76.gt.0)
               write(lut,99)trim(text(isp:isp+74))
               isp=isp+75
            enddo
            if(len_trim(text(isp:)).gt.0) write(lut,99)trim(text(isp:))
         enddo
      endif
! repeat the endmember loop again for interaction parameters (and magnetism??)
! I still have to figure out how to reference interacting constituents
! for the moment just terminate with a line starting with 0
      write(lut,300)
300   format(' 0')
   enddo phases1
!-------------------------------------------------------
! now data for stoichiometric phases
   mphstoi=1
   phases2: do i1=1,noofph
      lokph=phasetuple(i1)%lokph
      if(i1.ne.phstoi(mphstoi)) then
!         write(*,*)'3C skipping mixture ',trim(phlista(lokph)%name)
         cycle phases2
      endif
      mphstoi=mphstoi+1
      write(*,*)'3C parameters for compound ',trim(phlista(lokph)%name)
      lokcs=phlista(lokph)%linktocs(1)
      varres=>ceq%phase_varres(lokcs)
      nsubl=1
      ionliq=.false.
      nsubl=phlista(lokph)%noofsubl
      write(lut,500)phlista(lokph)%name
500   format(1x,a,5x,20('='),' compound')
! there is just one endmember!!
      endmember=>phlista(lokph)%ordered
      constext=' '
      ip=1
      constcomp=zero
      sloop2: do isubl=1,nsubl
! this is the loop for the constituents in sublattices
         isp=endmember%fraclinks(isubl,1)
         if(isp.eq.-99) then
! this means wildcard in this sublattice
            write(*,*)'3C Beware! Wildcard in a stoichiometric compound!!!'
            constext(ip:)='*:'
            ip=ip+2
            cycle sloop2
         endif
! Hm we should add stoichiometric factors for all constituents in this subl
         isp=phlista(lokph)%constitlist(isp)
         if(btest(splista(isp)%status,SPVA)) then
            write(*,*)'3C vacancy in stoichiometric compound!!'
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
! for the parameters follow the property link
      property=>endmember%propointer
      if(associated(property)) then
! For a compound do not write any constituent array
!         write(lut,100)constext(1:ip-2)
! this line should be written together with the type of coefficients and ranges
! it may require several lines
         write(text,210)constcomp
! what about several properties??
         f1=property%degreelink(0)
         if(f1.gt.0) then
            call list_tpascoef(lut,text,f1,npows,tpfc)
         else
            write(*,*)'missing endmember parameter'
         endif
      endif
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
   write(*,*)' *** WARNING: DAT format not fully implemented'
! 
1000 continue
   if(allocated(tpfc)) deallocate(tpfc)
   write(*,1010)trim(filename)
1010 format(/'Output finished on ',a)
   close(lut)
   return
 end subroutine save_datformat

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

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
   if(abs(splista(isp)%charge).gt.1-0D-6) then
      call wrinum(name,jp,6,1,splista(isp)%charge)
   endif
!   write(*,*)'3C suck: lower case name: ',trim(name)
   constext(ip:)=name
   ip=len_trim(constext)+1
1000 continue
   return
 end subroutine lower_case_species_name

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 logical function gtp_error_message(reset)
! tests the error code and writes the error message (if any) 
! and reset error code if reset=0
! if reset >0 that is set as new error message
! if reset <0 the error code is not changed
! return TRUE if error code set, FALSE if error code is zero
   implicit none
   integer reset
!\end{verbatim}
   if(gx%bmperr.ne.0) then
      if(gx%bmperr.ge.4000 .and. gx%bmperr.le.nooferm) then
         write(kou,10)gx%bmperr,bmperrmess(gx%bmperr)
10       format(' *** Error ',i5/a)
      elseif(gx%bmperr.ne.0) then
         write(*,20)gx%bmperr
20       format('3C Error without message: ',i7)
      endif
      if(reset.eq.0) then
! if reset zero reset error code
         gx%bmperr=0
      elseif(gx%bmperr.gt.0) then
! if reset positive set this as error code
         gx%bmperr=reset
      endif
! if reset negative do not change error code.  Set function to TRUE
      gtp_error_message=.TRUE.
   else
! no error, return false
      gtp_error_message=.FALSE.
   endif
1000 continue
   return
 end function gtp_error_message
   
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/
! subroutines extracted from the user interface

!\begin{verbatim}
  subroutine enterphase(cline,last)
! interactive entering of phase
    character cline*(*)
    integer last
!    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    character name1*24,text*80,name3*24,model*72,phtype*1,ch1*1
    integer nsl,defnsl,icon,ll,jp
    double precision sites(9)
    character (len=34) :: quest1='Number of sites on sublattice xx: '
! constituent indices in a phase
    integer, dimension(maxconst) :: knr
! array with constituents in sublattices when entering a phase
    character, dimension(maxconst) :: const*24
    logical once
!
    call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
! ionic liquid require special sorting of constituents on anion sublattice
    call capson(name1)
    defnsl=1
    if(name1(1:4).eq.'GAS ') then
       phtype='G'
       model='CEF-RKM'
    elseif(name1(1:7).eq.'LIQUID ') then
       phtype='L'
       model='CEF-RKM'
    elseif(name1(1:9).eq.'IONIC_LIQ') then
       phtype='L'
       model='IONIC_LIQUID'
       defnsl=2
    else
       phtype='S'
       model='CEF-RKM'
    endif
    call gparid('Number of sublattices: ',cline,last,nsl,defnsl,q1help)
    if(buperr.ne.0) goto 900
    if(nsl.le.0) then
       write(kou,*)'At least one configurational space!!!'
       goto 1000
    elseif(nsl.ge.10) then
       write(kou,*)'Maximum 9 sublattices'
       goto 1000
    endif
    icon=0
    sloop: do ll=1,nsl
! 'Number of sites on sublattice xx: '
!  123456789.123456789.123456789.123
       once=.true.
4042   continue
       write(quest1(31:32),4043)ll
4043   format(i2)
       call gparrd(quest1,cline,last,sites(ll),one,q1help)
       if(buperr.ne.0) goto 900
       if(sites(ll).le.1.0D-6) then
          write(kou,*)'Number of sites must be larger than 1.0D-6'
          if(once) then
             once=.false.
             goto 4042
          else
             goto 1000
          endif
       endif
! This should be extended to allow several lines of input
! 4 means up to ;
       once=.true.
4045   continue
       call gparc('All Constituents: ',cline,last,4,text,';',q1help)
       if(buperr.ne.0) goto 900
       knr(ll)=0
       jp=1
4047   continue
       if(eolch(text,jp)) goto 4049
       if(model(1:13).eq.'IONIC_LIQUID ' .and. ll.eq.1 &
            .and. knr(1).eq.0) then
! a very special case: a single "*" is allowed on 1st sublattice for ionic liq
          if(text(jp:jp).eq.'*') then
             icon=icon+1
             const(icon)='*'
             knr(1)=1
             cycle sloop
          endif
       endif
       call getname(text,jp,name3,1,ch1)
       if(buperr.eq.0) then
          icon=icon+1
          const(icon)=name3
          knr(ll)=knr(ll)+1
!          write(*,66)'constituent: ',knr(ll),icon,jp,const(icon)
66        format(a,3i3,a)
! increment jp to bypass a separating , 
          jp=jp+1
          goto 4047
       elseif(once) then
!          write(kou,*)'Input error ',buperr,', at ',jp,', please reenter'
          buperr=0; once=.false.; goto 4045
       else
          goto 1000
       endif
       buperr=0
4049   continue
    enddo sloop
    call enter_phase(name1,nsl,knr,const,sites,model,phtype)
    if(gx%bmperr.ne.0) goto 1000
900 continue
    if(buperr.ne.0) gx%bmperr=buperr
1000 continue
  end subroutine enterphase

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

!\begin{verbatim}
  subroutine listoptcoeff(lut)
! listing of optimizing coefficients
    integer lut
!    integer lut,mexp
!    double precision errs(*)
!    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    type(gtp_equilibrium_data), pointer :: neweq
    integer i1,i2,j1,j2,j3
    character name1*24,line*80
    double precision xxx
!
    write(lut,610)
610 format(/'List of coefficents with non-zero values'/&
         'Name  Current value   Start value    Scaling factor',&
         ' RSD')
    name1=' '
    do i1=0,size(firstash%coeffstate)-1
       coeffstate: if(firstash%coeffstate(i1).ge.10) then
! optimized variable, read from TP constant array
          call get_value_of_constant_index(firstash%coeffindex(i1),xxx)
          call makeoptvname(name1,i1)
          write(lut,615)name1(1:3),xxx,&
               firstash%coeffstart(i1),firstash%coeffscale(i1),zero
615       format(a,2x,4(1pe15.6))
          if(firstash%coeffstate(i1).eq.11) then
! there is a prescribed minimum
             write(lut,616)' minimum ',firstash%coeffmin(i1)
616          format(6x,'Prescribed ',a,': ',1pe12.4)
          elseif(firstash%coeffstate(i1).eq.12) then
! there is a prescribed maximum
             write(lut,616)' maximum ',firstash%coeffmax(i1)
          elseif(firstash%coeffstate(i1).eq.13) then
! there are prescribed minimum and maximum
             write(lut,617)firstash%coeffmin(i1),firstash%coeffmax(i1)
617          format(6x,'Prescribed min and max: ',2(1pe12.4))
          elseif(firstash%coeffstate(i1).gt.13) then
             write(lut,*)'Wrong coefficent state, set to 10'
             firstash%coeffstate(i2)=10
          endif
       elseif(firstash%coeffstate(i1).gt.0) then
! fix variable status
          call get_value_of_constant_index(firstash%coeffindex(i1),xxx)
          call makeoptvname(name1,i1)
          write(lut,615)name1(1:3),xxx
       elseif(firstash%coeffscale(i1).ne.0) then
! coefficient with negative status, status set to 1
          call get_value_of_constant_index(firstash%coeffindex(i1),xxx)
          write(lut,619)i1,firstash%coeffscale(i1),xxx,zero
619       format('Wrong state for coefficient ',i3,4(1pe12.4))
          firstash%coeffstate(i1)=1
       endif coeffstate
    enddo
1000 continue
    return
  end subroutine listoptcoeff

!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/

