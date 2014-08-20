!
! included in pmod25.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     11. Save and read things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine gtpsave(filename,str)
! save all data on file, unformatted, TDB or macro NOT EVEN STARTED
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
! separate UNFORMATTED, TDB and MACRO
   if(str(1:1).eq.'U') then
      call gtpsaveu(filename)
   else
      call gtpsavetm(filename,str)
   endif
1000 continue
   return
 end subroutine gtpsave

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!-\begin{verbatim}
 subroutine gtpsaveu(filename)
! MUST BE COMPLETELY REWRITTEN AGAIN ....
! save all data on file, formatted or unformatted
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
! equilibrium record(s) with conditions, componenets, phase_varres records etc
   implicit none
   character*(*) filename
!-\end{verbatim}
!
   character id*55,comment*72
   logical form
   integer freetpfun,i,isp,jph,kontroll,lenc,lokph
10  format(i8)
! only unformatted output
   form=.FALSE.
   comment='no comments'
   lenc=11
   if(index(filename,'.').eq.0) then
      if(form) then
         filename(len_trim(filename)+1:)='.ocf'
      else
         filename(len_trim(filename)+1:)='.ocu'
      endif
   endif
   write(*,*)'Unformatted save not implemented yet'
   goto 1000
!
   if(form) then
      open(21,file=filename,access='sequential',status='unknown',&
           form='formatted')
   else
      open(21,file=filename,access='sequential',status='unknown',&
           form='unformatted')
   endif
!       123456789.123456789.123456789.123456789.123456789.12345
   id='This is a save file for the free thermodynamic system, version: '
! this contol number will be written regularly on the file and checked on read
   kontroll=175638
!------ write first some id, version etc.
   if(form) then
      write(21,100)id,gtpversion,lenc
      write(21,101)comment(1:lenc)
      write(21,102)noofel,noofsp,noofph
100    format('HEADER: ',a,', 'a/'&',a)
101    format('TEXT: ',a)
102    format('NOOFEL, NOOFSP, NOOFPH: ',3i10)
   else
      write(21)id,gtpversion,lenc
      write(21)comment(1:lenc)
      write(21)noofel,noofsp,noofph
      write(21)1,kontroll
   endif
!----------------------------------------------------------------------
!
! it is extremely important to keep the order of the records as they
! are linked using indices
!
!---------- elementlist
   if(form) then
      write(21,120)(i,ellista(i),i=1,noofel)
! note double (( because format repeated for each element))
120 format(('ELEMENT ',i3,': ',a2,', ',a24,', ',a24/'&',3(1PE16.8)/'&',4i10))
   else
      write(21)(ellista(i),i=1,noofel)
      write(21)2,kontroll
   endif
!---------- specieslist
   do isp=1,noofsp
      if(form) then
         write(21,130)isp,splista(isp)%symbol,splista(isp)%mass,&
              splista(isp)%charge
         write(21,131)splista(isp)%noofel,splista(isp)%status, &
              splista(isp)%alphaindex
         write(21,132)(splista(isp)%ellinks(i),i=1,splista(isp)%noofel)
         write(21,133)(splista(isp)%stoichiometry(i),i=1,splista(isp)%noofel)
130       format('SPECIES ',i4,': ',A,2(1PE16.8))
131       format('&',i10,Z10,i10)
132       format('&',10i7)
133       format(('&',5(1PE16.8)))
      else
         write(21)splista(isp)%symbol,splista(isp)%mass,splista(isp)%charge
         write(21)splista(isp)%noofel,splista(isp)%status, &
              splista(isp)%alphaindex
         write(21)(splista(isp)%ellinks(i),i=1,splista(isp)%noofel)
         write(21)(splista(isp)%stoichiometry(i),i=1,splista(isp)%noofel)
      endif
   enddo
   if(.not.form) write(21)3,kontroll
!----------- phaselist, start from 0 (reference phase)
! including sublattces, endmembers, interactions, properties etc
   do jph=0,noofph
      lokph=phases(jph)
!      write(*,*)'calling savephase: ',jph,lokph
      call savephase(21,form,lokph)
      if(gx%bmperr.ne.0) goto 1000
   enddo
   if(.not.form) write(21)4,kontroll
!------------- tpfuns
   call tpfunsave(21,form)
   if(.not.form) write(21)5,kontroll
!------------- state variable functions
   call svfunsave(21,form,firsteq)
   if(.not.form) write(21)6,kontroll
!------------- references
   call refsave(21,form)
   if(.not.form) write(21)7,kontroll
!-------------------------------------------------------
! write the equilibrium records, at present for FIRSTEQ only
! conditions, components, phase_varres for all composition sets etc
   call saveequil(21,form,firsteq)
   if(.not.form) write(21)8,kontroll
!-------------------------------------------------------
   if(form) then
      write(21,*)'- END OF DATA - '
   else
      write(21)'- END OF DATA - '
      write(21)'- END OF DATA - '
   endif
900 continue
   close(21)
1000 continue
   return
 end subroutine gtpsaveu

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine savephase(lut,form,jph)
! save data for phase jph except data in the equilibrium record
! For phases with disordered set of parameters we must access the number of
! sublattices via firsteq ... UNFINISHED
   implicit none
   integer lut,jph
   logical form
!\end{verbatim}
   integer doneord,i,j,level,lokcs,nem,noi,nop,nox,nsl,nup
   type(gtp_endmember), pointer :: emrec
   type(gtp_interaction), pointer :: intrec
   type(gtp_property), pointer :: proprec
   type saveint
      type(gtp_interaction), pointer :: p1
   end type saveint
   type(saveint), dimension(:), pointer :: stack
   type(gtp_phase_add), pointer :: addlink
!
   allocate(stack(5))
!   write(*,*)'in savephase ',jph,form
   if(form) then
! as phlista contain pointers one has to write each item separately
      write(21,140)jph,phlista(jph)%name,&
           phlista(jph)%models,phlista(jph)%phletter,&
           phlista(jph)%status1,phlista(jph)%alphaindex,phlista(jph)%noofcs,&
           phlista(jph)%nooffs
140    format(('PHASE ',i3,': ',a24/'&',a,1x,a/'&',z10,2i10/'&',5i10))
   else
      write(21)jph,phlista(jph)%name,&
           phlista(jph)%models,phlista(jph)%phletter,&
           phlista(jph)%status1,phlista(jph)%alphaindex,phlista(jph)%noofcs,&
           phlista(jph)%nooffs
   endif
!----------- sublatticerecord
   nsl=phlista(jph)%noofsubl
   if(form) then
!      write(21,150)jph,nsl,phlista(jph)%cslink,&
      write(21,150)jph,nsl,phlista(jph)%linktocs,&
           phlista(jph)%tnooffr
      write(21,151)(phlista(jph)%nooffr(i),i=1,nsl)
      write(21,152)(phlista(jph)%sites(i),i=1,nsl)
      write(21,153)(phlista(jph)%constitlist(i),i=1,phlista(jph)%tnooffr)
150   format('SUBLATTICE RECORD',i3,': ',i3,2i10)
151   format('&',10i4)
152   format(('&',5(1PE15.7)))
153   format(('&',15i5))
   else
!      write(21)nsl,phlista(jph)%cslink,phlista(jph)%tnooffr
      write(21)nsl,phlista(jph)%linktocs,phlista(jph)%tnooffr
      write(21)(phlista(jph)%nooffr(i),i=1,nsl)
      write(21)(phlista(jph)%sites(i),i=1,nsl)
      write(21)(phlista(jph)%constitlist(i),i=1,phlista(jph)%tnooffr)
   endif
!--------- endmember list, interaction tree and property records
! save all parameter data starting from the endmember list
   doneord=0
   emrec=>phlista(jph)%ordered
! there can be phases without any ordered parameters ...
   if(.not.associated(emrec)) goto 400
! we come back here if there are disordered parameters
200 continue
! if doneord=1 then we have listed the ordered parameters
   if(doneord.eq.1) then
      emrec=>phlista(jph)%disordered
   endif
!   write(*,*)'savephase label 200: ',jph,doneord,nsl
   emlista: do while(associated(emrec))
      proprec=>emrec%propointer
      intrec=>emrec%intpointer
      nop=0
      noi=0
      nem=0
      if(associated(proprec)) nop=1
      if(associated(intrec)) noi=1
      if(associated(emrec%nextem)) nem=1
      if(form) then
         write(lut,210)emrec%noofpermut,emrec%phaselink,emrec%antalem,nop,noi
210      format(3x,'ENDMEMBER: ',3i7,2x,2i2)
         do j=1,emrec%noofpermut
            write(lut,211)(emrec%fraclinks(i,j),i=1,nsl)
         enddo
211      format('& ',10i4)
      else
         write(lut)emrec%noofpermut,emrec%phaselink,emrec%antalem,nop,noi,nem
         do j=1,emrec%noofpermut
            write(lut)(emrec%fraclinks(i,j),i=1,nsl)
         enddo
      endif
! endmember property
      emproplista: do while(associated(proprec))
         nox=0
         if(associated(proprec%nextpr)) nox=1
         if(form) then
            write(lut,220)proprec%reference,proprec%proptype,&
                 proprec%degree,proprec%extra,proprec%antalprop,nox
220         format(9x,'PROPERTY: ',a,2x,2i3,2i7,i2)
            do i=0,proprec%degree
               call save1tpfun(21,form,proprec%degreelink(i))
            enddo
         else
            write(lut)proprec%reference,proprec%proptype,&
                 proprec%degree,proprec%extra,proprec%antalprop,nox
            do i=0,proprec%degree
               call save1tpfun(21,form,proprec%degreelink(i))
            enddo
         endif
         proprec=>proprec%nextpr
      enddo emproplista
! interaction tree
      level=0
300   continue
      intlista: do while(associated(intrec))
! noi is next, nup is higher, nop is property
         noi=0
         nup=0
         nop=0
         if(associated(intrec%nextlink)) noi=1
         if(associated(intrec%highlink)) nup=1
         if(associated(intrec%propointer)) nop=1
!         write(*,*)'saving interaction: ',intrec%status
310      continue
         if(form) then
!            write(lut,320)intrec%noofpermut,intrec%status,noi,nup,nop
! code saving intrec need update here
            write(lut,320)intrec%noofip(1),intrec%status,noi,nup,nop
            do i=1,intrec%noofip(1)
               write(lut,321)intrec%sublattice(i),intrec%fraclink(i)
            enddo
320         format(6x,'INTERACTION: ',2i7,2x,3i2)
321         format('& ',2i7)
         else
            write(lut)intrec%noofip,intrec%status,noi,nup,nop
            do i=1,intrec%noofip(1)
               write(lut)intrec%sublattice(i),intrec%fraclink(i)
            enddo
         endif
! interaction property
         proprec=>intrec%propointer
         intproplista: do while(associated(proprec))
            nox=0
            if(associated(proprec%nextpr)) nox=1
            if(form) then
               write(lut,220)proprec%reference,proprec%proptype,&
                    proprec%degree,proprec%extra,proprec%antalprop,nox
               do i=0,proprec%degree
                  call save1tpfun(21,form,proprec%degreelink(i))
               enddo
            else
               write(lut)proprec%reference,proprec%proptype,&
                    proprec%degree,proprec%extra,proprec%antalprop,nox
               do i=0,proprec%degree
                  call save1tpfun(21,form,proprec%degreelink(i))
               enddo
            endif
            proprec=>proprec%nextpr
         enddo intproplista
! take link to higher higher interaction
         level=level+1
         if(level.gt.5) then
!            write(*,*)'Too many interaction levels'
            gx%bmperr=4164; goto 1000
         endif
!         write(*,*)'pushing level; ',level,intrec%status
         stack(level)%p1=>intrec
         intrec=>intrec%highlink
      enddo intlista
! pop previous intrec and take link to next interaction
!      write(*,*)'pop level: ',level
      if(level.gt.0) then
         intrec=>stack(level)%p1
!         write(*,*)'popping ',intrec%status
         intrec=>intrec%nextlink
!         if(associated(intrec)) then
!            write(*,*)'next of poped: ',intrec%status
!         endif
         level=level-1
         goto 300
      endif
!---- next endmember
      emrec=>emrec%nextem
   enddo emlista
! no more endmembers, check if the disordered (if any) has been written
400 continue
   if(doneord.eq.0) then
      if(associated(phlista(jph)%disordered)) then
! there are some disordered parameters
! we should maybe write the disfra record .... or is that done in saveequil??
! we have to change nsl ...three % vojvoj
         doneord=1
!         lokcs=phlista(jph)%cslink
         lokcs=phlista(jph)%linktocs(1)
         nsl=firsteq%phase_varres(lokcs)%disfra%ndd
         if(form) then
            write(21,490)'DISORDERED PARAMETERS: ',nsl
490         format(a,i3)
         else
! write that there are disordered parameters and the number of sublattices
            write(21)2
            write(21)nsl
         endif
! emrec should already be null but for security ....
         nullify(emrec)
         goto 200
      elseif(.not.form) then
! we must write that there are no disordered parameters
         write(21)0
      endif
   endif
!------ additions list
500 continue
   addlink=>phlista(jph)%additions
   addition: do while(associated(addlink))
      if(form) then
         nox=0
         if(associated(addlink%nextadd)) nox=1
         write(21,520)addlink%type,addlink%addrecno,addlink%aff,nox
         if(addlink%type.eq.1) then
! hm, these are expressions, not TPFUNS
!            call save1tpfun(21,form,addlink%explink(1))
!            call save1tpfun(21,form,addlink%explink(2))
         else
            write(*,*)'Saving unknown addition record type ',addlink%type
         endif
520       format('ADDRECORD ',3i7,2x,i2)
521       format('&',2i10)
      else
         write(21)addlink%type,addlink%addrecno,addlink%aff,nox
         if(addlink%type.eq.1) then
! magnetic, 2 links, hm AFF above is just for magnetic too ....
!            call save1tpfun(21,form,addlink%explink(1))
!            call save1tpfun(21,form,addlink%explink(2))
         else
            write(*,*)'Not saving unknown addition record type ',addlink%type
         endif
      endif
      addlink=>addlink%nextadd
   enddo addition
   if(.not.form) write(21)0,0,0,0
! ---- all done !! ?? !!
1000 continue
   return
 end subroutine savephase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine saveequil(lut,form,ceq)
! save data for an equilibrium record
   implicit none
   logical form
   integer lut
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character text*256
   type(gtp_phase_varres), pointer :: firstvarres
   TYPE(gtp_fraction_set), pointer :: fslink
   TYPE(gtp_condition), pointer :: condrec
   integer i,isp,j,k,kl,lokcs,lokph,mc,mc2,nsl
   write(*,17)ceq%eqname,ceq%eqno,ceq%status,ceq%next
   if(form) then
      write(lut,17)ceq%eqname,ceq%eqno,ceq%status,ceq%next
17    format('EQUILIBRIUM: ',a,i5,z8,i5)
   else
      write(lut)ceq%eqname,ceq%eqno,ceq%status,ceq%next
   endif
! ignore svfunres and eq_tpres
!---- components
   do i=1,noofel
      isp=ceq%complist(i)%splink
      if(form) then
         write(lut,21)splista(isp)%symbol,ceq%complist(i)%phlink,&
              ceq%complist(i)%status,&
              ceq%complist(i)%refstate,ceq%complist(i)%tpref
21       format(a,2i7,z8,2x/'& ',a,2(1pe12.4))
      else
         write(lut)ceq%complist(i)%splink
         write(lut)ceq%complist(i)%phlink,&
              ceq%complist(i)%status,&
              ceq%complist(i)%refstate,ceq%complist(i)%tpref
      endif
   enddo
   do i=1,noofel
      if(form) then
         write(lut,30)(ceq%compstoi(j,i),j=1,noofel)
30       format(6(1pe12.4))
      else
         write(lut)(ceq%compstoi(j,i),j=1,noofel)
      endif
   enddo
!---- conditions
   call get_all_conditions(text,0,ceq)
   if(gx%bmperr.ne.0) goto 1000
   kl=index(text,'CRLF')
   if(form) then
      write(lut,50)text(1:kl-1)
50    format(a)
   else
      write(lut)kl-1
      if(kl.gt.1) then
         write(lut)text(1:kl-1)
      endif
   endif
!---- experiments
   call get_all_conditions(text,1,ceq)
   if(gx%bmperr.ne.0) goto 1000
   kl=len_trim(text)
   if(form) then
      if(kl.gt.1) write(lut,50)text(1:kl)
   else
      write(lut)kl-1
      if(kl.gt.1) then
         write(lut)text(1:kl-1)
      endif
   endif
!---- varres records, one for each composition set
   if(form) then
      write(21,70)highcs
70    format('Number of composition sets: ',i7)
   else
      write(21)highcs
   endif
   compset: do j=1,highcs-1
! loop for all composition sets
!      write(*,*)'composition set index: ',j
      firstvarres=>ceq%phase_varres(j)
      if(btest(firstvarres%status2,CSDFS)) then
! this phase_varres/parres records belong to disordered fraction_set
! A big tricky to find the number of sublattices and constituents ....
         lokph=firstvarres%phlink
!         lokcs=phlista(lokph)%cslink
         lokcs=phlista(lokph)%linktocs(1)
         nsl=ceq%phase_varres(lokcs)%disfra%ndd
         mc=ceq%phase_varres(lokcs)%disfra%tnoofxfr
      else
         lokph=0; lokcs=0
         nsl=phlista(firstvarres%phlink)%noofsubl
         mc=phlista(firstvarres%phlink)%tnooffr
      endif
      mc2=mc*(mc+1)/2
      formatted: if(form) then
         write(21,170)j,firstvarres%nextfree,firstvarres%phlink,&
              firstvarres%status2
         write(21,171)firstvarres%prefix,firstvarres%suffix
         write(21,175)firstvarres%abnorm
         write(21,174)(firstvarres%constat(i),i=1,mc)
         write(21,172)(firstvarres%yfr(i),i=1,mc)
         write(21,172)(firstvarres%mmyfr(i),i=1,mc)
         write(21,173)(firstvarres%sites(i),i=1,nsl)
! derivatives of site fractions for ionic liquids ignored at present
!         write(21)(firstvarres%dsitesdy(i),i=1,mc)
!         write(21)(firstvarres%d2sitesdy2(i),i=1,mc2)
170      format('PHASE_VARRES: ',i5,2i10,z10)
171      format('&',a,1x,a)
172      format(('&',5(1PE15.7)))
173      format(('&',5(1PE15.7)))
174      format(('&',8(1X,Z8)))
175      format(('&',2(1PE15.7)))
!
         write(21,180)j,firstvarres%amfu,firstvarres%netcharge,&
              firstvarres%dgm,firstvarres%nprop
         write(21,181)(firstvarres%gval(i,1),i=1,6)
         do k=1,mc
            write(21,182)(firstvarres%dgval(i,k,1),i=1,3)
         enddo
         write(21,183)(firstvarres%d2gval(i,1),i=1,mc2)
180      format('PARRES: ',i3,3(1PE13.5),1x,i2)
181      format('&',6(1PE13.5))
182      format('&',3(1PE16.8))
183      format(('&',5(1PE15.7)))
         if(btest(firstvarres%status2,CSDLNK)) then
            write(*,*)'fraction set linked from: ',j,fslink%varreslink
            fslink=>firstvarres%disfra
            write(21,300)fslink%latd,fslink%ndd,fslink%tnoofxfr,&
                 fslink%tnoofyfr,fslink%totdis,fslink%varreslink,fslink%id
            write(21,302)fslink%nooffr,fslink%splink
            write(21,301)fslink%dsites
            write(21,302)fslink%y2x
            write(21,301)fslink%dxidyj
300         format('DISORDRED FRACTION SET: ',6i7,' ',a)
301         format(('&',5(1PE15.7)))
302         format(('&',10i7))
         endif
      else
! unformatted
         write(21)firstvarres%nextfree,firstvarres%phlink,&
              firstvarres%status2
         write(21)firstvarres%prefix,firstvarres%suffix
         write(21)firstvarres%abnorm
         write(21)(firstvarres%constat(i),i=1,mc)
         write(21)(firstvarres%yfr(i),i=1,mc)
!         write(*,*)'mmyfr: ',mc,size(firstvarres%mmyfr)
         write(21)(firstvarres%mmyfr(i),i=1,mc)
         write(21)(firstvarres%sites(i),i=1,nsl)
! These should only be interesting for ionic liquids and in that case
! only the dimension, not the values
!          write(21)(firstvarres%dsitesdy(i),i=1,mc)
!          write(21)(firstvarres%d2sitesdy2(i),i=1,mc2)
         lokph=firstvarres%phlink
!         write(*,*)'checking disordered fraction set'
         fsrec: if(btest(firstvarres%status2,CSDLNK)) then
! we must indicate on the file a fraction_set record follows!
            write(*,*)'fraction set linked from: ',j,fslink%varreslink
            write(21)1
            fslink=>firstvarres%disfra
            write(21)fslink%latd,fslink%ndd,fslink%tnoofxfr,&
                 fslink%tnoofyfr,fslink%totdis,fslink%varreslink,fslink%id
            write(21)fslink%nooffr,fslink%splink
            write(21)fslink%dsites
            write(21)fslink%y2x
            write(21)fslink%dxidyj
         else
! no disordered fraction set record
            write(21)0
         endif fsrec
         write(21)firstvarres%amfu,firstvarres%netcharge,firstvarres%dgm,&
              firstvarres%nprop
! only G values saved ???? well maybe not even those ...
         write(21)(firstvarres%gval(i,1),i=1,6)
         do k=1,mc
            write(21)(firstvarres%dgval(i,k,1),i=1,3)
         enddo
         write(21)(firstvarres%d2gval(i,1),i=1,mc2)
      endif formatted
   enddo compset
1000 continue
   return
 end subroutine saveequil

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine gtpread(filename,str)
! being updated
! read all data in the following order
! header
! element list
! species list
! phase list with sublattices, endmembers, interactions and parameters etc
! tpfuns
! state variable functions
! references
! equilibrium record(s) with conditions, componenets, phase_varres records etc
!
   implicit none
   character*(*) filename,str
!\end{verbatim}
   character id*55,start*8,endoffile*16,version*6
   integer highexpr,freexpr,freetpfun,i,i1,i2,isp,jph,kontroll,lenc,nel
   TYPE(gtp_phase_varres) :: firstvarres
! strange error with len()
10  format(i8)
   if(index(filename,'.').eq.0) then
      filename(len_trim(filename)+1:)='.gtp'
   endif
   kontroll=175638
   open(21,file=filename,access='sequential',status='unknown',&
        form='unformatted')
!------ read some identification etc
   read(21)id,version,lenc
   if(version.ne.gtpversion) then
      write(*,11)id,version,gtpversion
11     format('File not same version as program: ',A,2x,a)
      gx%bmperr=2901; goto 1000
   endif
   read(21)str(1:lenc)
   str(lenc+1:)=' '
!   write(*,*)'text on saved file: ',str(1:lenc+1)
   read(21)noofel,noofsp,noofph
!-------
   read(21)i1,i2
   if(i1.ne.1 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 1'
      gx%bmperr=4165; goto 1000
!   else
!      write(*,*)'read header OK'
   endif
!---------- elementlist
   read(21)(ellista(i),i=1,noofel)
!-------
   read(21)i1,i2
   if(i1.ne.2 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 2'
      gx%bmperr=4165; goto 1000
!   else
!      write(*,*)'read elements OK'
   endif
!---------- specieslist
   do isp=1,noofsp
      read(21)splista(isp)%symbol,splista(isp)%mass,splista(isp)%charge
      read(21)splista(isp)%noofel,splista(isp)%status, &
           splista(isp)%alphaindex
!       write(*,*)'bmpread 1F: ',isp,splista(isp)%noofel,splista(isp)%symbol
      if(isp.gt.1) then
         nel=splista(isp)%noofel
         allocate(splista(isp)%ellinks(nel))
         allocate(splista(isp)%stoichiometry(nel))
      endif
      read(21)(splista(isp)%ellinks(i),i=1,splista(isp)%noofel)
      read(21)(splista(isp)%stoichiometry(i),i=1,splista(isp)%noofel)
   enddo
!-------
   read(21)i1,i2
   if(i1.ne.3 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 3'
      gx%bmperr=4165; goto 1000
!   else
!      write(*,*)'read species OK'
   endif
!----------- phaselist, starting from 0, the reference phase
   noofem=0
   noofint=0
   noofprop=0
   do jph=0,noofph
      call readphase(21)
      if(gx%bmperr.ne.0) goto 1000
   enddo
!--------
   read(21)i1,i2
   if(i1.ne.4 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 4'
      gx%bmperr=4165; goto 1000
!   else
!      write(*,*)'read phases OK'
   endif
!---------- tpfuns
   call tpfunread(21)
   if(gx%bmperr.ne.0) goto 1000
!--------
   read(21)i1,i2
   if(i1.ne.5 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 5'
      gx%bmperr=4165; goto 1000
!   else
!      write(*,*)'read TPFUNS OK'
   endif
!---------- state variable functions
   call svfunread(21)
   if(gx%bmperr.ne.0) goto 1000
!--------
   read(21)i1,i2
   if(i1.ne.6 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 6'
      gx%bmperr=4165; goto 1000
!   else
!      write(*,*)'read state variable functions OK'
   endif
!---------- references
   call refread(21)
   if(gx%bmperr.ne.0) goto 1000
!--------
   read(21)i1,i2
   if(i1.ne.7 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 7'
      gx%bmperr=4165; goto 1000
!   else
!      write(*,*)'read references OK'
   endif
!---------- equilibrium record
   call readequil(21,firsteq)
   if(gx%bmperr.ne.0) goto 900
!   if(gx%bmperr.ne.0) goto 1000
!--------
   read(21)i1,i2
   if(i1.ne.8 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 8'
      gx%bmperr=4165; goto 1000
!   else
!      write(*,*)'read equilibrium records OK'
   endif
!------ read all ??
   endoffile=' '
   read(21,end=800,err=800)endoffile
800 continue
   if(endoffile.ne.'- END OF DATA - ') then
      write(kou,811)endoffile
811    format('End of file mark not was expted:'/'>',A,'<')
       gx%bmperr=4166; goto 1000
   endif
! emergency exit
900 continue
   close(21)
! this is just setting the link from elements list to ellista etc.
   do i=1,noofel
      elements(ellista(i)%alphaindex)=i
   enddo
   do i=1,noofsp
      species(splista(i)%alphaindex)=i
   enddo
   do i=1,noofph
      phases(phlista(i)%alphaindex)=i
   enddo
!
! program dies between the return here and line after the call in ocmon ...
! because I used = instead of => in readequil ....
!
   write(*,1007)noofel,noofsp,noofph,&
        notpf(),csfree-1,eqfree-1,nsvfun,reffree-1,addrecs
1007 format('Created records for elements, species, phases: ',2x,3i5/&
          'end members, interactions, properties: ',10x,3i5/&
          'TP-funs, composition sets, equilibria: ',10x,3i5/&
          'state variable functions, references, additions: ',4i5)
1000 continue
   return
 end subroutine gtpread

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readphase(lut)
! read data for phlista and all endmembers etc
! works for test case without disordered fraction test
   implicit none
   integer lut
!\end{verbatim}
   integer firstendmem,i,i1,i2,i3,i4,jph,level,nem,noi,nop,nox,nup,nsl
   type(gtp_endmember), pointer :: emrec
   type(gtp_interaction), pointer :: intrec
   type(gtp_property), pointer :: proprec
   type saveint
      type(gtp_interaction), pointer :: p1
      integer noi
   end type saveint
   type(saveint), dimension(:), pointer :: stack
   type(gtp_phase_add), pointer :: addlink
!
   allocate(stack(5))
!   write(*,*)'in readphase:'
! as phlista record contain pointers each item must be read separately
!----------- phlista
   read(21)jph,phlista(jph)%name,&
        phlista(jph)%models,phlista(jph)%phletter,&
        phlista(jph)%status1,phlista(jph)%alphaindex,phlista(jph)%noofcs,&
        phlista(jph)%nooffs
!   write(*,*)'read some phase data',jph,phlista(jph)%name
!----------- sublatticelist
!   read(lut)phlista(jph)%noofsubl,phlista(jph)%cslink,phlista(jph)%tnooffr
   read(lut)phlista(jph)%noofsubl,phlista(jph)%linktocs,phlista(jph)%tnooffr
   nsl=phlista(jph)%noofsubl
   allocate(phlista(jph)%nooffr(nsl))
   allocate(phlista(jph)%sites(nsl))
   read(lut)(phlista(jph)%nooffr(i),i=1,nsl)
   read(lut)(phlista(jph)%sites(i),i=1,nsl)
   allocate(phlista(jph)%constitlist(phlista(jph)%tnooffr))
   read(lut)(phlista(jph)%constitlist(i),i=1,phlista(jph)%tnooffr)
!------ endmember records, these must be allocated and linked
   nullify(phlista(jph)%ordered)
   nullify(phlista(jph)%disordered)
!   write(*,*)'read some more phase data',nsl
! return here when an endmember list empty. Read also disordered endmembers
200 continue
   read(lut)nem
!   write(*,*)'read even more phase data',nem
   if(nem.eq.0) then
! if morendmem=0 there are no more endmembers
      goto 500
   endif
   nem=1
   firstendmem=nem
!   write(*,202)'reading parameters: ',phlista(jph)%name,nsl,firstendmem
202 format(a,a,10i4)
   newendmem: do while(nem.eq.1)
      if(firstendmem.eq.1) then
! the first ordered endmember
         call readendmem(21,nsl,phlista(jph)%ordered,nop,noi,nem)
         emrec=>phlista(jph)%ordered
         firstendmem=0
      elseif(firstendmem.eq.2) then
! we must also read the number of sublattices in the disordered part
! and later more data may be needed here.
         read(21)nsl
! the first disordered endmember
         call readendmem(21,nsl,phlista(jph)%disordered,nop,noi,nem)
         emrec=>phlista(jph)%disordered 
         firstendmem=0
      else
! this the second or later endmember in the same list
         call readendmem(21,nsl,emrec%nextem,nop,noi,nem)
         emrec=>emrec%nextem
      endif
      if(nop.eq.1) then
! there are some property records !!
         call readproprec(21,emrec%propointer,nox)
         proprec=>emrec%propointer
         do while(nox.eq.1)
            call readproprec(21,proprec%nextpr,nox)
            proprec=>proprec%nextpr
         enddo
      endif
      inttree: if(noi.eq.1) then
! first interaction parameter record
         level=0
         call readintrec(21,emrec%intpointer,noi,nup,nop)
         intrec=>emrec%intpointer
!         write(*,13)'read interaction: ',intrec%status,noi,nup,nop
13       format(a,10i4)
300      continue
         if(nop.eq.1) then
            call readproprec(21,intrec%propointer,nox)
            proprec=>intrec%propointer
            do while(nox.eq.1)
               call readproprec(21,proprec%nextpr,nox)
               proprec=>proprec%nextpr
            enddo
         endif
! push before going to higher
330      continue
         level=level+1
         stack(level)%p1=>intrec
         stack(level)%noi=noi
!         write(*,13)'pushed interaction: ',intrec%status,0,0,0,level
         higher: if(nup.eq.1) then
! go to higher level and save intrec
            call readintrec(21,intrec%highlink,noi,nup,nop)
            intrec=>intrec%highlink
!            write(*,13)'read higher interaction: ',intrec%status,noi,nup,nop
            if(nop.eq.1) then
! there are some property records !!
               call readproprec(21,intrec%propointer,nox)
               proprec=>intrec%propointer
               do while(nox.eq.1)
                  call readproprec(21,proprec%nextpr,nox)
                  proprec=>proprec%nextpr
               enddo
            endif
            goto 330
         endif higher
! we come here when no higher records, pop records from stack
350      continue
         pop: if(level.gt.0) then
            intrec=>stack(level)%p1
            noi=stack(level)%noi
            level=level-1
!            write(*,13)'poped interaction: ',intrec%status,0,0,0,level
            if(noi.eq.1) then
               call readintrec(21,intrec%nextlink,noi,nup,nop)
               intrec=>intrec%nextlink
!               write(*,13)'read interaction: ',intrec%status,noi,nup,nop
               goto 300
            else
               goto 350
            endif
         endif pop
      endif inttree
   enddo newendmem
! we come nere when no more endmembers 
! have we checked the disordered list?
   goto 200
!------ additions list
500 continue
   nullify(phlista(jph)%additions)
510 continue
   read(21)i1,i2,i3,i4
!   write(*,*)'Addition; ',i1
   if(i1.eq.1) then
! here addition record should be created but as I have not managed to
! save the functions I just skip this for the moment and use create_magrec
      call create_magrec_inden(addlink,i3)
      if(gx%bmperr.ne.0) goto 1000
      if(.not.associated(phlista(jph)%additions)) then
         phlista(jph)%additions=>addlink
      else
         addlink%nextadd=>addlink
         addlink=>addlink%nextadd
      endif
      goto 510
   endif
1000 continue
   return
 end subroutine readphase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readendmem(lut,nsl,emrec,nop,noi,nem)
! allocates and reads an endmember record
   implicit none
   integer lut,nsl,nop,noi,nem
   type(gtp_endmember), pointer :: emrec
!\end{verbatim}
   integer i,j,i1,i2,i3,i4,jph
   integer level,nox,nup
   allocate(emrec)
   read(lut)emrec%noofpermut,emrec%phaselink,emrec%antalem,nop,noi,nem
!   write(*,17)'readendmem: ',nsl,emrec%noofpermut,emrec%antalem,nop,noi,nem
!17 format(a,6i5)
   allocate(emrec%fraclinks(nsl,emrec%noofpermut))
   do j=1,emrec%noofpermut
      read(lut)(emrec%fraclinks(i,j),i=1,nsl)
   enddo
   nullify(emrec%nextem)
   nullify(emrec%propointer)
   nullify(emrec%intpointer)
   noofem=noofem+1
   return
 end subroutine readendmem

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readproprec(lut,proprec,nox)
! allocates and reads a property record
   implicit none
   integer lut,nox
   type(gtp_property), pointer :: proprec
!\end{verbatim}
   integer i
   allocate(proprec)
   read(lut)proprec%reference,proprec%proptype,&
        proprec%degree,proprec%extra,proprec%antalprop,nox
   allocate(proprec%degreelink(0:proprec%degree))
!   write(*,17)'readprop: ',proprec%proptype,proprec%degree,&
!        proprec%antalprop,nox
!17 format(a,6i5)
   do i=0,proprec%degree
      call read1tpfun(21,proprec%degreelink(i))
   enddo
   nullify(proprec%nextpr)
   noofprop=noofprop+1
   return
 end subroutine readproprec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readintrec(lut,intrec,noi,nup,nop)
! allocates and reads an interaction record UNFINISHED
   integer none
   integer lut,noi,nup,nop
   type(gtp_interaction), pointer :: intrec
!\end{verbatim}
   allocate(intrec)
!   read(lut)intrec%noofip,intrec%status,noi,nup,nop
!   write(*,17)'readint: ',intrec%noofip,intrec%status,noi,nup,nop
!17 format(a,6i5)
!   allocate(intrec%sublattice(intrec%noofip))
!   allocate(intrec%fraclink(intrec%noofip))
!   do i=1,intrec%noofip(1)
!      read(lut)intrec%sublattice(i),intrec%fraclink(i)
!   enddo
   nullify(intrec%nextlink)
   nullify(intrec%highlink)
   nullify(intrec%propointer)
   return
 end subroutine readintrec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readequil(lut,ceq)
! Read equilibria records from a file
   implicit none
   integer lut
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character text*512
   type(gtp_phase_varres), pointer :: firstvarres
   TYPE(gtp_fraction_set) :: fslink
   integer i,ic,ierr,ip,isp,ivar,j,jp,k,lokcs,lokph,mc,mc2,nprop,nsl
! containing conditions, components and phase varres records for wach compset
   read(lut)ceq%eqname,ceq%eqno,ceq%status,ceq%next
!   write(*,*)'Reading equilibrium: ',ceq%eqname
!----- components
!   allocate(ceq%complist(noofel)) already allocated for 20
!   write(*,*)'Size of component array: ',size(ceq%complist)
   do i=1,noofel
      read(lut)isp
      ceq%complist(i)%splink=isp
      read(lut)ceq%complist(i)%phlink,&
           ceq%complist(i)%status,&
           ceq%complist(i)%refstate,ceq%complist(i)%tpref
   enddo
! stoichiometry conversion matrix, already allocated??
!   allocate(ceq%compstoi(noofel,noofel))
   ceq%compstoi=zero
   do i=1,noofel
      read(lut)(ceq%compstoi(j,i),j=1,noofel)
   enddo
! calculate the inverse stoichiometry matrix ... later
   call mdinv(maxel,maxel+1,ceq%compstoi,ceq%invcompstoi,ic,ierr)
!   write(*,*)'Calculatied inverse stoichiometry matrix'
!   do j=1,noofel
!      write(*,*)(ceq%compstoi(i,j),i=1,noofel)
!   enddo
!   do j=1,noofel
!      write(*,*)(ceq%invcompstoi(i,j),i=1,noofel)
!   enddo
!----- conditions, can be empty
   read(lut)ip
   if(ip.gt.0) then
      read(lut)text(1:ip)
!      write(*,*)'Conditions: ',text(1:ip)
! set the conditions, ip will be incremented by 1 in enter_condintion
      jp=0
      call set_condition(text(1:ip),jp,firsteq)
      if(gx%bmperr.ne.0) goto 1000
   endif
!---- experiments
   read(lut)ip
   if(ip.gt.0) then
      read(lut)text(1:ip)
   endif
!----------- phase_varres record
   read(lut)highcs
!   write(*,*)'Number of phase_varres records: ',highcs-1
!
!   write(*,*)'phase_varres size: ',size(ceq%phase_varres)
   do j=1,highcs-1
!      write(*,*)'reading phase_varres ',j
!------------------------------------------
! DEBUGPROBLEM BEWARE, using = instead of => below took 2 days to find
!------------------------------------------
! >>>      firstvarres=ceq%phase_varres(j)    <<< error
      firstvarres=>ceq%phase_varres(j)
      read(lut)firstvarres%nextfree,firstvarres%phlink,&
           firstvarres%status2
      read(lut)firstvarres%prefix,firstvarres%suffix
      read(lut)firstvarres%abnorm
! check interconecctions, firstvarres%phlink is phase record index
! from phlista(firstvarres%phlink)%clink one should find this record (j)
!      jxph=firstvarres%phlink
      if(btest(firstvarres%status2,CSDFS)) then
! this phase_varres records belong to a disordered fraction_set
         lokph=firstvarres%phlink
!         lokcs=phlista(lokph)%cslink
         lokcs=phlista(lokph)%linktocs(1)
         nsl=ceq%phase_varres(lokcs)%disfra%ndd
         mc=ceq%phase_varres(lokcs)%disfra%tnoofxfr
      else
         nsl=phlista(firstvarres%phlink)%noofsubl
         mc=phlista(firstvarres%phlink)%tnooffr
      endif
!       write(*,*)'bmpread 78: ',j,nsl,mc
      mc2=mc*(mc+1)/2
! added integer status array constat
!      write(*,*)'Allocate constat 1: ',nc
      allocate(firstvarres%constat(mc))
      allocate(firstvarres%yfr(mc))
      allocate(firstvarres%sites(nsl))
      read(lut)(firstvarres%constat(i),i=1,mc)
      read(lut)(firstvarres%yfr(i),i=1,mc)
      allocate(firstvarres%mmyfr(mc))
      read(lut)(firstvarres%mmyfr(i),i=1,mc)
      read(lut)(firstvarres%sites(i),i=1,nsl)
! these are ignored, values not important but must be allocated
! for ionic liquid as (2,mc) and (2,mc2)
!       read(lut)(firstvarres%dsitesdy(2,i),i=1,mc)
!       read(lut)(firstvarres%d2sitesdy2(2,i),i=1,mc2)
      read(lut)ivar
      if(ivar.eq.1) then
! extra fraction set
         write(*,*)'reading extra fraction set for ',j
         read(lut)fslink%latd,fslink%ndd,fslink%tnoofxfr,&
              fslink%tnoofyfr,fslink%totdis,fslink%varreslink,fslink%id
         allocate(fslink%nooffr(fslink%ndd))
         allocate(fslink%dsites(fslink%ndd))
         allocate(fslink%splink(fslink%tnoofxfr))
         allocate(fslink%y2x(fslink%tnoofyfr))
         allocate(fslink%dxidyj(fslink%tnoofyfr))
         read(lut)fslink%nooffr,fslink%splink
         read(lut)fslink%dsites
         read(lut)fslink%y2x
         read(lut)fslink%dxidyj
! now copy fslink to the correct record and then deallocate fslink arrays
         call copy_fracset_record(j,fslink,ceq)
         deallocate(fslink%nooffr)
         deallocate(fslink%dsites)
         deallocate(fslink%splink)
         deallocate(fslink%y2x)
         deallocate(fslink%dxidyj)
      endif
! The result data
      read(lut)firstvarres%amfu,firstvarres%netcharge,firstvarres%dgm, &
           firstvarres%nprop
      nprop=firstvarres%nprop
      allocate(firstvarres%listprop(nprop))
      allocate(firstvarres%gval(6,nprop))
      allocate(firstvarres%dgval(3,mc,nprop))
      allocate(firstvarres%d2gval(mc2,nprop))
      read(lut)(firstvarres%gval(i,1),i=1,6)
      do k=1,mc
         read(21)(firstvarres%dgval(i,k,1),i=1,3)
      enddo
      read(lut)(firstvarres%d2gval(i,1),i=1,mc2)
!      write(*,*)'phase_varres size: ',j,size(ceq%phase_varres)
   enddo
!
   csfree=highcs
1000 continue
   return
 end subroutine readequil

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
 subroutine new_gtp
! deletes all data
!
! MUST BE COMPLETELY REWRITTEN
!
! this is also needed before reading a new file (or same file ...)
! we must go through all records and delete and deallocate each
! separately.  Copied from "bmpread" so some comments irrelevant.
!   implicit double precision (a-h,o-z)
   implicit none 
!\end{verbatim}
   integer isp,j,lokcs,lokph,mc,mc2,nel,nsl
   TYPE(gtp_equilibrium_data) :: ceq
   TYPE(gtp_fraction_set) :: fslink
!   write(*,*)'Removing current data is not yet implemented'
   gx%bmperr=4172
   goto 1000
!---------- elementlist, just restore free list, done at the end
!---------- specieslist, we have to deallocate
   ceq=firsteq
   do isp=1,noofsp
      nel=splista(isp)%noofel
      deallocate(splista(isp)%ellinks)
      deallocate(splista(isp)%stoichiometry)
   enddo
!---------- components, just restore free list
! phases, many records
!----------- sublatticelist
   do j=0,noofph
      deallocate(phlista(j)%nooffr)
      deallocate(phlista(j)%sites)
      deallocate(phlista(j)%constitlist)
!      phlista(j)%addlink=0
      phlista(j)%noofcs=0
      phlista(j)%nooffs=0
!      phlista(j)%endmember(1)=0
!      phlista(j)%endmember(2)=0
!      phlista(j)%endmember(3)=0
   enddo
! some records can have holes in the list
!    write(*,18)'bmpread: ',highem,highexpr,highadd,highcs
!----------- phase_varres list
   do j=1,highcs-1
      if(btest(ceq%phase_varres(j)%status2,CSDFS)) then
! this phase_varres records belong to disordered fraction_set
         lokph=ceq%phase_varres(j)%phlink
!         lokcs=phlista(lokph)%cslink
         lokcs=phlista(lokph)%linktocs(j)
         nsl=ceq%phase_varres(lokcs)%disfra%ndd
         mc=ceq%phase_varres(lokcs)%disfra%tnoofxfr
      else
         nsl=phlista(ceq%phase_varres(j)%phlink)%noofsubl
         mc=phlista(ceq%phase_varres(j)%phlink)%tnooffr
      endif
      mc2=mc*(mc+1)/2
! added integer status array constat
      deallocate(ceq%phase_varres(j)%constat)
      deallocate(ceq%phase_varres(j)%yfr)
      deallocate(ceq%phase_varres(j)%mmyfr)
      deallocate(ceq%phase_varres(j)%sites)
      deallocate(ceq%phase_varres(j)%listprop)
      deallocate(ceq%phase_varres(j)%gval)
      deallocate(ceq%phase_varres(j)%dgval)
      deallocate(ceq%phase_varres(j)%d2gval)
   enddo
!------ endmember list
! endmembers etc are stored only up to highem etc.  Read these values now
! for all free lists
!   do j=1,highem-1
!      if(emlista(j)%phaselink.gt.0) then
!         nsl=phlista(emlista(j)%phaselink)%noofsubl
!      else
!         nsl=0
!      endif
! for unused emlista records the fraclinks have not been allocated.
!      if(nsl.gt.0) then
!         deallocate(emlista(j)%fraclinks)
!      endif
!   enddo
!------ additions list
!   do j=1,highadd-1
!      if(addlink%type.eq.1) then
! magnetic, 2 links
!         deallocate(addlink%link)
!      endif
!   enddo
!------ tpfunction expressions and other lists
! TP functions, normally 500 tpfuns are created
   call tpfun_deallocate
!----------- tpfunction result list
! reference records
   do j=1,reffree-1
      deallocate(reflista(j)%refspec)
   enddo
!----------------
   deallocate(ellista)
   deallocate(elements)
!    deallocate(splista)
!    deallocate(species)
!    deallocate(
! now initiate all lists and a little more
   call init_gtp
!   write(*,*)'new 5: ',highem,highexpr,freetpfun
!
1000 continue
   return
 end subroutine new_gtp

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\begin{verbatim}
 subroutine svfunsave(lut,form,ceq)
! saves all state variable functions on a file
   implicit none
   integer lut
   logical form
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character text*512,symbols(20)*32,afterdot*32
   integer indstv(4),indices(4)
   integer ii,ip,ipos,istv,js,jt,kl,ks,lrot
   if(form) then
      write(lut,17)nsvfun
17    format('STATE VARIABLE FUNCTIONS: ',i7)
   else
      write(lut)nsvfun
   endif
   do lrot=1,nsvfun
      ipos=1
      if(svflista(lrot)%narg.eq.0) goto 500
      js=0
      jt=0
100 continue
      jt=jt+1
      js=js+1
      ip=1
      symbols(js)=' '
      istv=svflista(lrot)%formal_arguments(1,jt)
      if(istv.lt.0) then
! function refer to another function
         symbols(js)=svflista(-istv)%name
      else
         do ii=1,4
            indices(ii)=svflista(lrot)%formal_arguments(1+ii,jt)
         enddo
         call encode_state_variable(symbols(js),ip,istv,indices,&
              svflista(lrot)%formal_arguments(6,jt), &
              svflista(lrot)%formal_arguments(7,jt),ceq)
         if(svflista(lrot)%formal_arguments(8,jt).ne.0) then
! a derivative!!!
            jt=jt+1
            afterdot=' '
            ip=1
            call encode_state_variable(afterdot,ip,&
                 svflista(lrot)%formal_arguments(1,jt),indices,&
                 svflista(lrot)%formal_arguments(6,jt), &
                 svflista(lrot)%formal_arguments(7,jt),ceq)
            symbols(js)=symbols(js)(1:len_trim(symbols(js)))//'.'//afterdot
         endif
      endif
      if(jt.lt.svflista(lrot)%narg) goto 100
500   continue
      kl=len_trim(svflista(lrot)%name)
      text(ipos:ipos+kl+1)=svflista(lrot)%name(1:kl)//'= '
      ipos=ipos+kl+2
      call wrtfun(text,ipos,svflista(lrot)%linkpnode,symbols)
      if(pfnerr.ne.0) then
         write(kou,*)'Putfun error listing funtion ',ks,pfnerr
         gx%bmperr=4142; goto 1000
      endif
      if(form) then
         call wrice2(lut,2,4,78,1,text(1:ipos-1))
      else
         write(lut)ipos-1
         write(lut)text(1:ipos-1)
      endif
   enddo
1000 continue
   return
 end subroutine svfunsave

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine refsave(lut,form)
! saves references on a file
   implicit none
   integer lut
   logical form
!\end{verbatim}
   character longline*2048
   integer ir,jp,ll,nl
   if(form) then
      write(lut,10)reffree-1
   else
      write(lut)reffree-1
   endif
10  format('REFERENCES: ',i5)
   do ir=1,reffree-1
      longline=reflista(ir)%reference
      jp=17
      nl=size(reflista(ir)%refspec)
      do ll=1,nl
         longline(jp:)=reflista(ir)%refspec(ll)
         jp=jp+64
      enddo
      jp=len_trim(longline)
      if(form) then
         call wrice(lut,0,17,78,longline(1:jp))
      else
         write(lut)jp
         write(lut)longline(1:jp)
      endif
   enddo
1000 continue
   return
 end subroutine refsave

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine svfunread(lut)
! read a state variable function from save file and store it.
! by default there are some state variable functions, make sure
! they are deleted.  Done here just by setting nsvfun=0
   implicit none
   integer lut
!\end{verbatim}
   integer nsvfun,i,ip,nsvfunfil
   character*512 text
   nsvfun=0
   read(lut)nsvfunfil
!   write(*,*)'Number of state variable functions: ',nsvfunfil
   do i=1,nsvfunfil
      read(lut)ip
!      write(*,*)'Number of characters: ',ip
      text=' '
      read(lut)text(2:ip)
!      write(*,*)text(2:ip)
      ip=1
      call enter_svfun(text,ip,firsteq)
      if(gx%bmperr.ne.0) then
!         write(*,*)'Error entering svf from file',gx%bmperr
         goto 1000
      endif
   enddo
1000 continue
   return
 end subroutine svfunread

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine refread(lut)
! read references from save file
   implicit none
   integer lut
!\end{verbatim}
   character text*512
   integer i,iref,jp,nrefs
   read(lut)nrefs
!   write(*,*)'Number of references: ',nrefs
   do i=1,nrefs
      read(lut)jp
      read(lut)text(1:jp)
!      write(*,*)text(1:jp)
      call tdbrefs(text(1:16),text(17:jp),0,iref)
      if(gx%bmperr.ne.0) goto 1000
   enddo
1000 continue
   return
 end subroutine refread

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine readtdb(filename)
! reading data from a TDB file.  No selection of elements.
! Several TYPE_DEFS not handelled
! first read only FUNCTION (maybe several times), then the rest
! Only minimally edited iutput files from TC accepted
   implicit none
   character filename*(*)
!\end{verbatim}
   character line*100,elsym*2,name1*24,name2*24,elsyms(10)*2
   character longline*2048,phtype*1,ch1*1,const(maxsp)*24,name3*24,funname*60
   character refx*16
   character (len=1), dimension(20) :: typedefchar
   integer, dimension(20) :: typedefaction
   integer, dimension(5) :: addphasetypedef
   double precision mass,h298,s298
   integer, dimension(10) :: elidx,knr,endm
! lint(1,*) is sublattice, lint(2,*) is species
   double precision stoik(10),xsl,xxx
   integer lint(2,3),noofphasetype,nytypedef
   integer typty,fractyp,lp1,lp2,ipnewref,ix,jph,kkk,lcs,nint,noelx
   logical onlyfun
   integer norew,newfun,nfail,nooftypedefs,nl,ipp,jp,jss,lrot,ip,jt
   integer nsl,ll,kp,nr,nrr,mode,lokph,ip1,km,nrefs,ideg,ierr,iph,ics
!
   write(*,*)'reading a TDB file'
   if(.not.(index(filename,'.tdb').gt.0 &
       .or. index(filename,'.TDB').gt.0)) then
! no extention provided
      filename(len_trim(filename)+1:)='.TDB'
   endif
   open(21,file=filename,access='sequential',form='formatted',&
        err=1010,iostat=ierr,status='old')
   onlyfun=.TRUE.
   norew=0
   newfun=0
   nfail=0
   nooftypedefs=0
! return here after rewind
90  continue
   nl=0
! return here for new keyword
100 continue
   read(21,110,end=2000)line
110 format(a)
   nl=nl+1
   barafun: if(onlyfun) then
! first read only functions until all has been read
      if(line(2:10).eq.'FUNCTION ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! FUNCTION GHSERCR    2.98150E+02  -8856.94+157.48*T-26.908*T*LN(T)
         name1=line(11:18)
! special case, error in TDB file, UN_ASS is only 6 characters
         if(name1(1:6).eq.'UN_ASS') then
            name1=line(11:16); ipp=18
         else
            ipp=20
         endif
         longline=' '
         longline=line(ipp:)
111       continue
         jp=len_trim(longline)
         if(longline(jp:jp).eq.'!') then
! replace # by ' '
112          continue
            jss=index(longline(1:jp),'#')
            if(jss.gt.0) then
               longline(jss:jss)=' '
               goto 112
            endif
! here search functions already enterd, there are freetpfun-1 of them
!+             do jss=1,freetpfun-1
!                write(*,*)'readtdb 111: ',jss,symbol,tpfuns(jss)%symbol
!+                if(name1.eq.tpfuns(jss)%symbol) then
!+                   goto 100
!+                endif
!+             enddo
!             write(*,*)'Entering function: ',name1
            call enter_tpfun(name1,longline,lrot)
            if(gx%bmperr.ne.0) then
! one may have error here if function calls other functions not entered, 4002
! or if the function is already entered, 4026
               if(gx%bmperr.eq.4002 .or. gx%bmperr.eq.4026) then
                  nfail=nfail+1
                  gx%bmperr=0; goto 100
               endif
               goto 1000
            endif
!             write(*,*)'Entered function: ',name1
            newfun=newfun+1
         else
            nl=nl+1
            read(21,110)line
!            write(kou,101)'readtdb 2: ',nl,line(1:40)
            longline=longline(1:jp)//line
            goto 111
         endif
      endif
      goto 100
   endif barafun
!    write(kou,*)'readtdb 100: ',nl,line(1:40)
!101 format('readtdb 1: ',i3,'>',a,'<')
   if(line(2:9).eq.'ELEMENT ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! ELEMENT CR   BCC_A2                    5.1996E+01  4.0500E+03  2.3560E+01!
      elsym=line(10:11)
      if(elsym.eq.'/-' .or. elsym.eq.'VA') goto 100
      name1=line(15:39)
      ip=40
      call getrel(line,ip,mass)
      if(buperr.ne.0) then
         mass=one; buperr=0
      endif
      call getrel(line,ip,h298)
      if(buperr.ne.0) then
         h298=zero; buperr=0
      endif
      call getrel(line,ip,s298)
      if(buperr.ne.0) then
         s298=zero; buperr=0
      endif
      name2=elsym
!       write(*,*)'readtdb 7B: ',elsym,mass
      call new_element(elsym,name2,name1,mass,h298,s298)
      if(gx%bmperr.ne.0) goto 1000
   elseif(line(2:9).eq.'SPECIES ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! SPECIES O3PU2                       O3PU2!
      name1=line(10:34)
      name2=line(38:)
      call decode_stoik(name2,noelx,elsyms,stoik)
      if(gx%bmperr.ne.0) goto 1000
! check elements exist
      call new_species(name1,noelx,elsyms,stoik)
      if(gx%bmperr.ne.0) goto 1000
   elseif(line(2:10).eq.'FUNCTION ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! FUNCTION GHSERCR    2.98150E+02  -8856.94+157.48*T-26.908*T*LN(T)
      name1=line(11:18)
      longline=' '
      longline=line(20:)
300    continue
      jp=len_trim(longline)
      if(longline(jp:jp).eq.'!') then
! all functions should already be entered
!          write(*,*)'Skipping function: ',name1
         continue
!          call enter_tpfun(name1,longline,lrot)
!          if(gx%bmperr.ne.0) goto 1000
      else
         nl=nl+1
         read(21,110)line
!          write(kou,101)'readtdb 2: ',nl,line(1:40)
         longline=longline(1:jp)//line
         goto 300
      endif
   elseif(line(2:7).eq.'PHASE ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! PHASE LIQUID:L %  1  1.0  !
      name1=line(8:)
      jp=index(name1,' ')
      ip=7+jp
      if(jp.gt.0) then
         name1(jp:)=' '
      endif
      jp=index(name1,':')
!      write(*,*)'readtdb 11: ',name1,ip,jp
! phytype
      if(jp.gt.0) then
         phtype=name1(jp+1:jp+1)
         name1(jp:)=' '
      else
         phtype=' '
      endif
! phase type code
      noofphasetype=0
      ip=ip+1
      jp=ip
      name2=line(ip:jp)
310    jp=jp+1
      typedefcheck: if(line(jp:jp).ne.' ') then
         ch1=line(jp:jp)
         do jt=1,nooftypedefs
            if(ch1.eq.typedefchar(jt)) goto 320
         enddo
         goto 310
320       continue
         if(typedefaction(jt).eq.99) then
! ignore TYPE_DEF SEQ
            continue
         elseif(typedefaction(jt).eq.-1 .or. &
              typedefaction(jt).eq.-3) then
! magnetic addition, save for after phase created
            noofphasetype=noofphasetype+1
            addphasetypedef(noofphasetype)=typedefaction(jt)
         else
            continue
         endif
         goto 310
      endif typedefcheck
      name2='TDB file model: '//name2
! sublattices
      call getrel(line,jp,xsl)
      if(buperr.ne.0) then
         gx%bmperr=buperr; goto 1000
      endif
      nsl=int(xsl)
      do ll=1,nsl
         call getrel(line,jp,stoik(ll))
      if(buperr.ne.0) then
         gx%bmperr=buperr; goto 1000
      endif
      enddo
!      write(kou,*)'readtdb 3A: ',nsl,(stoik(i),i=1,nsl)
!    CONSTITUENT LIQUID:L :CR,FE,MO :  !
      read(21,110)line
      nl=nl+1
!      write(kou,*)'readtdb 3B: ',nl,line(1:40)
      kp=index(line,' :')
      longline=line(kp:)
360    continue
      jp=len_trim(longline)
!      write(*,*)'readtdb gas1: ',nl,jp,longline(1:jp)
      if(longline(jp:jp).ne.'!') then
! continue reading tdbfile
         nl=nl+1
         read(21,110)line
! modified 2012.08.08/BoS: skip  3 spaces at the beginning of each line
! modified 2013.01.06/BoS: skip  4 spaces at the beginning of each line !!??
         longline=longline(1:jp)//line(5:)
!         write(*,*)'readtdb gas1: ',jp,longline(jp-5:jp+10)
         goto 360
      endif
! eliminate all after the exclamation mark
      longline(jp:)=' '
! 
      ip=index(longline,' :')+2
!      write(*,*)'readtdb gas2: ',jp,longline(1:jp)
      ll=0
      nr=0
      nrr=0
!      write(*,*)'readtdb 3C: ',ll,nr,nsl,longline(ip:jp)
! mode=1 indicates to getname that / + - are allowed in species names
      mode=1
370   continue
      if(ll.ge.1) knr(ll)=nr
      ll=ll+1
!      write(*,*)'start sublat ',ll,nsl,nr,ip
      if(ll.gt.nsl) goto 390
      nr=0
380   continue
      nr=nr+1
      nrr=nrr+1
!      write(*,379)'readtdb 3CXX: ',ip,nr,longline(ip:ip+10)
379   format(a,2i4,' >',a,'< >',a,'< >',a,'<')
      call getname(longline,ip,name3,mode,ch1)
!      write(*,379)'readtdb 3CY: ',ip,nr,longline(ip:ip+10),name3,ch1
      if(buperr.ne.0) then
         write(*,381)'readtdb 3E: ',ll,nr,longline(1:ip+5),ip,name3
381      format(a,2i4,' "',a,'" ',i5,1x,a,'"',a)
         gx%bmperr=buperr
         goto 1000
      endif
!      write(*,381)'readtdb 3E: ',ll,nr,longline(1:ip+5),ip,name3,ch1
      const(nrr)=name3
      ip=ip+1
385    continue
      if(ch1.eq.',') goto 380
      if(ch1.eq.':') goto 370
      if(ch1.eq.' ') then
! on ouput from TC a sublattice constituent list ends by " : "
         ip=ip+2
         goto 370
      endif
! this is the major constituent indicator, ignore
      if(ch1.eq.'%') then
         ch1=longline(ip:ip)
         ip=ip+1
!          write(*,*)'readtdb 6A: ',ip,ch1
         goto 385
      endif
      write(*,*)'Species terminator error: ',ch1,nl
      gx%bmperr=4157; goto 1000
390    continue
! name2 is model, ignored on reading TDB
      name2='CEF-TDB-RKM? '
!       write(*,*)'readtdb 9: ',name1,nsl,knr(1),knr(2),phtype
      call new_phase(name1,nsl,knr,const,stoik,name2,phtype)
!      write(*,*)'readtdb 9A: ',gx%bmperr
      if(gx%bmperr.ne.0) goto 1000
! any typedefs? only magnetic handelled at present
      call find_phase_by_name(name1,iph,lcs)
!       write(*,*)'readtdb 9X: ',gx%bmperr
      if(gx%bmperr.ne.0) goto 1000
      lokph=phases(iph)
!      write(*,*)'typedefs for ',name1(1:20),lokph,noofphasetype
      phasetypes: do jt=1,noofphasetype
!          write(*,*)'typedef ',jt,addphasetypedef(jt)
         if(addphasetypedef(jt).eq.-1) then
            call add_magrec_inden(lokph,1,-1)
         elseif(addphasetypedef(jt).eq.-3) then
            call add_magrec_inden(lokph,1,-3)
         endif
         if(gx%bmperr.ne.0) goto 1000
      enddo phasetypes
!      write(*,*)'readtdb 9B:',name1,nsl,phtype
   elseif(line(4:13).eq.'PARAMETER ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
!   PARAMETER G(LIQUID,CR;0)  2.98150E+02  +24339.955-11.420225*T
      funname=line(14:)
      kp=index(funname,' ')
      funname(kp:)=' '
! extract symbol, normally G or L but TC, BMAGN and others can occur
      lp1=index(funname,'(')
      name1=funname(1:lp1-1)
! this is kept for compatibility with TDB files generated by TC
      if(name1(1:2).eq.'G ' .or. name1(1:2).eq.'L ') then
         typty=1
      elseif(name1(1:3).eq.'TC ') then
         typty=2
      elseif(name1(1:6).eq.'BMAGN ') then
         typty=3
      endif
! only fractype 1 on TDB files unles I handle disordered part
      fractyp=1
!       write(*,*)'readtdb: PAR',name1,typty
! extract phase name and constituent array
      lp1=index(funname,'(')
      lp2=index(funname,',')
      name2=funname(lp1+1:lp2-1)
      call find_phase_by_name_exact(name2,jph,kkk)
!       write(*,*)'readtdb 19: ',jph,gx%bmperr,name2
      if(gx%bmperr.ne.0) goto 1000
! extract constituent array, remove final ) and decode
      lokph=phases(jph)
      name3=funname(lp2+1:)
      lp1=len_trim(name3)
      name3(lp1:)=' '
      call decode_constarr(name3,nsl,endm,nint,lint,ideg)
!      write(*,303)'readtdb 303: ',name3,nsl,endm(1),nint,&
!           ((lint(i,j),i=1,2),j=1,nint)
303   format(a,a,2i3,2x,i2,' : ',3(2i3,2x))
      if(gx%bmperr.ne.0) then
         write(*,*)'readtdb error: ',gx%bmperr,name3
         goto 1000
      endif
      if(nint.gt.1) then
! lint(1,1) is species of first, lint(1,2) in second interaction
!          write(*,305)'readtdb 305: ',endm(1),nint,lint(2,1),lint(2,2)
      endif
305    format(a,5i4)
!---------------- encode function
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
! now read the function
      jp=kp+13
      longline=' '
      longline=line(jp:)
410    continue
      jp=len_trim(longline)
      if(longline(jp:jp).eq.'!') then
! extract reference if any
         kp=jp-1
         do while(longline(kp:kp).eq.' ')
            kp=kp-1
         enddo
         lp1=kp
         do while(longline(kp:kp).ne.'N')
            kp=kp-1
            if(kp.lt.1) then
! illegal termination of function in TDB file
               write(*,413)nl
413 format('illegal termination of function in TDB file, around line: ',i5)
               gx%bmperr=4013; goto 1000
            endif
         enddo
         refx=longline(kp+1:lp1)
! replace # by ' '
412       continue
         jss=index(longline(1:jp),'#')
         if(jss.gt.0) then
            longline(jss:jss)=' '
            goto 412
         endif
         call enter_tpfun(funname,longline,lrot)
!          write(*,17)lokph,typty,nsl,lrot,(endm(i),i=1,nsl)
17 format('readtdb 17: '4i3,5x,10i3)
!         write(*,404)'readtdb entpar: ',refx,fractyp,nint,ideg
404      format(a,a,i3,2x,10i3)
         call enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,&
              lrot,refx)
         if(gx%bmperr.ne.0) then
            write(*,*)'readtdb entparerror: ',gx%bmperr
            goto 1000
         endif
      else
! continue reading tdbfile
         nl=nl+1
         read(21,110)line
!          write(kou,101)'readtdb 2: ',nl,line(1:40)
         longline=longline(1:jp)//line
         goto 410
      endif
   elseif(line(2:17).eq.'TYPE_DEFINITION ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! TYPE_DEFINITION & GES A_P_D BCC_A2 MAGNETIC  -1.0    4.00000E-01 !
      nytypedef=nooftypedefs+1
      nooftypedefs=nytypedef
      typedefchar(nytypedef)=line(18:18)
      newtypedef: if(line(20:22).eq.'SEQ') then
         typedefaction(nytypedef)=100
      else
         km=index(line,' MAGNETIC ')
         magnetic: if(km.gt.0) then
            ip=km+9
!             write(*,73)'decode typedef ',ip,line(ip:80)
!73           format(a,i3,' "',a,'"')
            call getrel(line,ip,xxx)
            if(buperr.ne.0) then
               gx%bmperr=buperr; goto 1000
            endif
! this can be -1 for BCC or -3 for FCC, HCP and other phases
            typedefaction(nytypedef)=int(xxx)
         else
            typedefaction(nytypedef)=99
            write(kou,*)'Cannot handle this TYPE_DEFINITION:'
            write(kou,*)line(17:)
         endif magnetic
      endif newtypedef
   elseif(line(2:20).eq.'LIST_OF_REFERENCES ' .or. &
          line(2:16).eq.'ADD_REFERENCES ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! LIST_OF_REFERENCES
! NUMBER  SOURCE
!   REF283  'Alan Dinsdale, SGTE Data for Pure Elements,
!          Calphad Vol 15(1991) p 317-425,
!          also in NPL Report DMA(A)195 Rev. August 1990'
!       write(kou,*)'Does not handle REFERENCES'
! skip the line with "NUMBER  SOURCE"
      nl=nl+1
      read(21,110)line
      nrefs=0
775    continue
      nl=nl+1
      read(21,110)line
! new reference, finish if empty line or line with just "!"
      ip=1
      if(eolch(line,ip)) goto 777
! some reference lists like those from SSUB has no sigle quotes
      kp=index(line,"'")
      citationmarks: if(kp.gt.0) then
         refx=line(ip:kp-1)
         longline=' '
         jp=1
         line=line(kp+1:)
         ip=1
776       continue
!          write(*,771)ip,kp,line(1:len_trim(line))
!771       format('tdbref 5: ',2i4,' "',a,'"')
         kp=index(line,"'")
         if(kp.gt.0) then
! second ' means end of reference text
!             write(*,772)'readtdb 6: ',jp,ip,kp,line(ip:kp)
            longline(jp:)=line(ip:kp-1)
            jp=jp+kp-ip
         else
            kp=len_trim(line)
!             write(*,772)'readtdb 7: ',jp,ip,kp,line(ip:kp)
!772          format(a,3i4,1x,a)
            longline(jp:)=line(ip:)
! add a space at newline
            jp=jp+kp-ip+2
778          continue
            nl=nl+1
            read(21,110)line
            ip=1
            if(eolch(line,ip)) goto 778
            goto 776
         endif
!          write(*,774)refx,longline(1:jp)
774       format('readtdb 8: ',a,'=',a)
         call tdbrefs(refx,longline(1:jp),1,ix)
         nrefs=nrefs+1
         goto 775
      else
! references without citation marks
! ip is at the start of the reference id, look for space
         refx=' '
         ipnewref=ip
780       continue
         line=line(ip:)
         kp=index(line,' ')
         refx=line(1:kp-1)
         longline=line(kp+1:)
         jp=len_trim(longline)+2
! whole first line belong to this reference, check if more lines
782       continue
         nl=nl+1
         read(21,110)line
         ip=1
         if(eolch(line,ip)) goto 777
         if(line(ip:ip).eq.'!') then
! end of all references
            if(refx(1:1).ne.' ') then
               call tdbrefs(refx,longline(1:jp),1,ix)
               nrefs=nrefs+1
            endif
            goto 777
         elseif(ip.eq.ipnewref) then
! new reference
            if(refx(1:1).ne.' ') then
!                write(*,*)'readtbd refx: ',refx,jp,ipnewref
               call tdbrefs(refx,longline(1:jp),1,ix)
               nrefs=nrefs+1
            endif
            goto 780
         else
! add this line to current reference
            longline(jp:)=line(ip:)
            jp=len_trim(longline)+2
            goto 782
         endif
      endif citationmarks
777    continue
!      write(*,*)'Read ',nrefs,' references'
! look for next KEYWORD
   endif
   goto 100
1000 continue
   if(buperr.ne.0 .or. gx%bmperr.ne.0) then
      write(*,1005)buperr,gx%bmperr,nl
1005   format('Error ',i5,' or ',i5', occured at TDB file line ',i7)
   endif
   close(21)
! endmember, interactions,property,tpfuns,composition sets,equilibria
! references,additions
   write(*,1007)noofel,noofsp,noofph,noofem,noofint,noofprop,&
        notpf(),csfree-1,eqfree-1,nsvfun,reffree-1,addrecs
1007 format('Created records for elements, species, phases: ',2x,3i5/&
          'end members, interactions, properties: ',10x,3i5/&
          'TP-funs, composition sets, equilibria: ',10x,3i5/&
          'state variable functions, references, additions: ',4i5)
   return
1010 continue
   write(*,*)'I/O error opening file: ',ierr
   gx%bmperr=7000+ierr; 
   return
! end of file, check if reading functions
2000 continue
   if(onlyfun) then
      rewind(21)
      norew=norew+1
      if(newfun.gt.0) then
!          write(*,*)'Read ',newfun+nfail,' functions, entered ',newfun,&
!               ' rewinding ',norew
         newfun=0
         nfail=0
      else
!         write(*,*)'All functions read, rewinding to read the rest'
         onlyfun=.FALSE.
      endif
      goto 90
   else
      goto 1000
   endif
 end subroutine readtdb

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!-\begin{verbatim}
 subroutine gtpsavetm(filename,str)
! MUST BE COMPLETELY REWRITTEN AGAIN ....
! save all data on file, formatted or unformatted
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
!-\end{verbatim}
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
! unfinished ....
! open file and write (either as TDB or MACRO):
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
