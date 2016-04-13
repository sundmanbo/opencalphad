!
! gtp3E included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     9. Save and read things from files
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
      call gtpsaveu(filename,str(3:))
   elseif(str(1:1).eq.'D') then
      call gtpsavedir(filename,str(3:))
   elseif(str(1:1).eq.'T') then
      call gtpsavetdb(filename,str(3:))
   elseif(str(1:1).eq.'L') then
      call gtpsavelatex(filename,str(3:))
   else
      call gtpsavetm(filename,str)
   endif
1000 continue
   return
 end subroutine gtpsave

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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

!\begin{verbatim} %-
 subroutine gtpsavetdb(filename,specification)
! save all data in TDB format on an file
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

!\begin{verbatim}
 subroutine gtpsaveu(filename,specification)
! save all data unformatted on an file
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
!
   character id*40,comment*72,endoffile*16,mark*8
   integer i,isp,jph,kontroll,lokph,lut
!
   if(index(filename,'.').eq.0) then
      filename(len_trim(filename)+1:)='.ocu'
   endif
   lut=21
   open(lut,file=filename,access='sequential',status='unknown',&
           form='unformatted',iostat=gx%bmperr,err=1000)
   id='This is a save file for OC version:    '
   comment=specification
! this control number will be written regularly on the file and checked on read
   kontroll=175638
   mark=' MARK '//char(13)//char(10)
!>>>>> 1: write first some id, version etc.
   write(lut)id,savefile,comment,globaldata
   write(lut)noofel,noofsp,noofph,nooftuples
   write(lut)1,kontroll,mark
!----------------------------------------------------------------------
! note the use of gtp_xxx_version to handle versions
!----------------------------------------------------------------------
!
! it is extremely important to keep the order of the records as they
! are linked using indices
!
!>>>>> 2: elementlist
   if(ocv()) write(*,*)'Writing elements'
   write(lut)gtp_element_version
   do i=1,noofel
      write(lut)ellista(i)
   enddo
!-----------
   write(lut)2,kontroll,mark
!>>>>> 3: specieslist
   if(ocv()) write(*,*)'Writing species'
   write(lut)gtp_species_version
   do isp=1,noofsp
      write(lut)splista(isp)%symbol,splista(isp)%mass,splista(isp)%charge
      write(lut)splista(isp)%noofel,splista(isp)%status, &
           splista(isp)%alphaindex
      write(lut)(splista(isp)%ellinks(i),i=1,splista(isp)%noofel)
      write(lut)(splista(isp)%stoichiometry(i),i=1,splista(isp)%noofel)
   enddo
   write(lut)3,kontroll,mark
!>>>>> 4: phaselist, start from 0 (reference phase)
! including sublattces, endmembers, interactions, properties etc
! save version of various records
   if(ocv()) write(*,*)'Writing phases'
   write(lut)gtp_phase_version,gtp_endmember_version,gtp_interaction_version,&
        gtp_property_version
   if(noofph.gt.0) then
      do jph=0,noofph
         lokph=phases(jph)
         call savephase(lut,lokph)
         if(gx%bmperr.ne.0) goto 1000
         if(ocv()) write(*,*)'Saved phase: ',jph
      enddo
   endif
   write(lut)(phasetuple(i),i=1,nooftuples)
   write(lut)4,kontroll,mark
!------------- tpfuns
!>>>>> 20: tpfuns
   if(ocv()) write(*,*)'Writing tpfuns'
   call tpfunsave(lut,.FALSE.)
   write(lut)5,kontroll,mark
!------------- state variable functions
!>>>>> 30: svfuns
   if(ocv()) write(*,*)'Writing state variable functions'
   call svfunsave(lut,firsteq)
   write(lut)6,kontroll,mark
!   write(*,*)'Writing mark: ',6,kontroll,mark
!------------- references
!>>>>> 40: bibliographic references
   if(ocv()) write(*,*)'Writing references'
   call bibliosave(lut)
   write(lut)7,kontroll,mark
!-------------------------------------------------------
! write the equilibrium records, at present for FIRSTEQ only
! conditions, components, phase_varres for all composition sets etc
!>>>>> 50: equilibria
   if(ocv()) write(*,*)'Writing equilibria'
   write(lut)gtp_equilibrium_data_version,gtp_component_version,&
        gtp_phase_varres_version
   call saveequil(lut,firsteq)
   write(lut)8,kontroll,mark
!-------------------------------------------------------
   endoffile='- END OF DATA - '
   write(lut)endoffile
900 continue
   close(lut)
1000 continue
   return
 end subroutine gtpsaveu

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine savephase(lut,lokph)
! save data for phase at location lokph (except data in the equilibrium record)
! For phases with disordered set of parameters we must access the number of
! sublattices via firsteq
   implicit none
   integer lut,lokph
!\end{verbatim}
   integer doneord,i,j,level,lokcs,nem,noi,nop,nox,nsl,nup,noendm,fipsize
   type(gtp_endmember), pointer :: emrec
   type(gtp_interaction), pointer :: intrec
   type(gtp_property), pointer :: proprec
! to keep track of interaction records
   type saveint
      type(gtp_interaction), pointer :: p1
   end type saveint
   type(saveint), dimension(:), pointer :: stack
   type(gtp_phase_add), pointer :: addlink
   if(ocv()) write(*,*)'In savephase'
   allocate(stack(5))
!>>>>> 5: phase header
   write(lut)lokph,phlista(lokph)%name,&
        phlista(lokph)%models,phlista(lokph)%phletter,&
        phlista(lokph)%status1,&
        phlista(lokph)%alphaindex,phlista(lokph)%noofcs,phlista(lokph)%nooffs
   nsl=phlista(lokph)%noofsubl
   emrec=>phlista(lokph)%ordered
   if(.not.associated(emrec)) then
      noendm=0
   else
      noendm=1
   endif
!>>>>> 6: sublattice info
   j=phlista(lokph)%tnooffr
   if(ocv()) write(*,10)j,lokph,size(phlista(lokph)%constitlist)
10 format('3E: ',3i20)
   if(ocv()) write(*,11)(phlista(lokph)%constitlist(i),i=1,j)
11 format('3E: ',20i3)
   write(lut)nsl,phlista(lokph)%linktocs,phlista(lokph)%tnooffr
   write(lut)(phlista(lokph)%nooffr(i),i=1,nsl),&
        (phlista(lokph)%constitlist(i),i=1,j),noendm
!--------- endmember list, interaction tree and property records
! save all parameter data starting from the endmember list
   doneord=0
   if(ocv()) write(*,*)'listing endmembers',doneord,nsl,noendm
! there can be phases without any ordered parameters ...
   if(.not.associated(emrec)) goto 400
! we come back here if there are disordered parameters
200 continue
! if doneord=1 then we have listed the ordered parameters
   if(doneord.eq.1) then
      emrec=>phlista(lokph)%disordered
      if(ocv()) write(*,*)'Saving disordered parameters'
   endif
   if(ocv()) write(*,*)'any endmember: ',doneord
   emlista: do while(associated(emrec))
      proprec=>emrec%propointer
      intrec=>emrec%intpointer
      nop=0
      noi=0
      nem=0
      if(associated(proprec)) nop=1
      if(associated(intrec)) noi=1
      if(associated(emrec%nextem)) nem=1
      if(ocv()) write(*,55)'writing endmember: ',nsl,emrec%noofpermut,&
           emrec%phaselink,emrec%antalem,nop,noi,nem
55    format(a,7i5)
!>>>>> 7: endmember record (basic or disordered)
      write(lut)emrec%noofpermut,emrec%phaselink,emrec%antalem,nop,noi,nem
      do j=1,emrec%noofpermut
         write(lut)(emrec%fraclinks(i,j),i=1,nsl)
      enddo
      emproplista: do while(associated(proprec))
         nox=0
         if(associated(proprec%nextpr)) nox=1
!>>>>> 8: endmember property record (loop)
         write(lut)proprec%reference,proprec%proptype,&
              proprec%degree,proprec%extra,proprec%antalprop,nox
         do i=0,proprec%degree
            call save1tpfun(lut,.FALSE.,proprec%degreelink(i))
         enddo
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
310      continue
!>>>>> 9: interaction record
! look in gtp3H, create_interaction for use of intec%noofip
         fipsize=size(intrec%noofip)
         write(lut)fipsize
         write(lut)intrec%noofip,intrec%status,noi,nup,nop
         do i=1,intrec%noofip(2)
            write(lut)intrec%sublattice(i),intrec%fraclink(i)
         enddo
! interaction property
         proprec=>intrec%propointer
         intproplista: do while(associated(proprec))
            nox=0
            if(associated(proprec%nextpr)) nox=1
!>>>>> 10: interaction property record (loop)
            write(lut)proprec%reference,proprec%proptype,&
                 proprec%degree,proprec%extra,proprec%antalprop,nox
            do i=0,proprec%degree
               call save1tpfun(lut,.FALSE.,proprec%degreelink(i))
            enddo
            proprec=>proprec%nextpr
         enddo intproplista
! take link to higher higher interaction
         level=level+1
         if(level.gt.5) then
!            write(*,*)'Too many interaction levels'
            gx%bmperr=4164; goto 1000
         endif
         stack(level)%p1=>intrec
         intrec=>intrec%highlink
      enddo intlista
! pop previous intrec and take link to next interaction
      if(level.gt.0) then
         intrec=>stack(level)%p1
         intrec=>intrec%nextlink
         level=level-1
         goto 300
      endif
!---- next endmember
      emrec=>emrec%nextem
   enddo emlista
! no more endmembers, check if the disordered (if any) has been written
400 continue
   if(doneord.eq.0) then
      if(ocv()) write(*,*)'any disordered endmembers?'
      if(associated(phlista(lokph)%disordered)) then
! there are some disordered parameters
! the disfra record is written in saveequil??
! we have to change nsl ...three % vojvoj
         doneord=1
         lokcs=phlista(lokph)%linktocs(1)
         nsl=firsteq%phase_varres(lokcs)%disfra%ndd
!>>>>> 11A: write disordered endmemebers
         write(lut)2,nsl
! emrec should already be null but for security ....
         nullify(emrec)
         goto 200
      else
! we must mark that there are no disordered parameters
!>>>>> 11B: no moe endmemebers
         write(lut)0,0
      endif
   endif
!------ additions list
500 continue
   addlink=>phlista(lokph)%additions
   addition: do while(associated(addlink))
      if(addlink%type.eq.1) then
!>>>>> 12A: additions id
         write(lut)addlink%type,addlink%addrecno,addlink%aff
      else
         write(*,*)'Not saving unknown addition record type ',addlink%type
      endif
      addlink=>addlink%nextadd
   enddo addition
!>>>>> 12B: mark end of data for phase
   write(lut)-1,-1,-1
1000 continue
   return
 end subroutine savephase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine saveequil(lut,ceq)
! save data for an equilibrium record
   implicit none
   integer lut
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character text*512
   type(gtp_phase_varres), pointer :: firstvarres
   TYPE(gtp_fraction_set), pointer :: fslink
!   TYPE(gtp_condition), pointer :: condrec
   integer i,isp,j,k,kl,lokcs,lokph,mc,mc2,nsl
!>>>>> 50:
   write(lut)ceq%eqname,ceq%eqno,ceq%status,ceq%next
! ignore svfunres and eq_tpres
!---- components
!>>>>> 51:
   do i=1,noofel
      isp=ceq%complist(i)%splink
      write(lut)isp
      write(lut)ceq%complist(i)%phlink,ceq%complist(i)%status,&
           ceq%complist(i)%refstate,ceq%complist(i)%tpref,&
           ceq%complist(i)%mass
   enddo
   do i=1,noofel
      if(ocv()) write(*,99)'comp.matrix: ',(ceq%invcompstoi(j,i),j=1,noofel)
   enddo
99 format(a,7e11.3)
   do i=1,noofel
      write(lut)(ceq%compstoi(j,i),j=1,noofel)
   enddo
!---- varres records, one for each composition set
!>>>>> 54:
   write(lut)highcs
   compset: do j=1,highcs
! loop for all composition sets
      firstvarres=>ceq%phase_varres(j)
      if(btest(firstvarres%status2,CSDFS)) then
! this phase_varres/parres records belong to disordered fraction_set
! A big tricky to find the number of sublattices and constituents ....
         lokph=firstvarres%phlink
         lokcs=phlista(lokph)%linktocs(1)
         nsl=ceq%phase_varres(lokcs)%disfra%ndd
         mc=ceq%phase_varres(lokcs)%disfra%tnoofxfr
      else
         lokph=0; lokcs=0
         nsl=phlista(firstvarres%phlink)%noofsubl
         mc=phlista(firstvarres%phlink)%tnooffr
      endif
      mc2=mc*(mc+1)/2
!>>>>> 55:
      write(lut)firstvarres%nextfree,firstvarres%phlink,&
           firstvarres%status2,firstvarres%phstate
      write(lut)firstvarres%prefix,firstvarres%suffix
      write(lut)firstvarres%abnorm
!>>>>> 56:
      write(lut)(firstvarres%constat(i),i=1,mc)
      write(lut)(firstvarres%yfr(i),i=1,mc)
      write(lut)(firstvarres%mmyfr(i),i=1,mc)
      write(lut)(firstvarres%sites(i),i=1,nsl)
! We do not save the cmuval array
! These should only be interesting for ionic liquids and in that case
! only the dimension, not the values
!          write(lut)(firstvarres%dsitesdy(i),i=1,mc)
!          write(lut)(firstvarres%d2sitesdy2(i),i=1,mc2)
      lokph=firstvarres%phlink
      fsrec: if(btest(firstvarres%status2,CSDLNK)) then
! we must indicate on the file a disordered fraction_set record follows!
         fslink=>firstvarres%disfra
         if(ocv()) write(*,*)'Disordered fraction set linked from: ',&
              j,fslink%varreslink
!>>>>> 57A: write disordered record, is is inside the phase_varres record
         write(lut)1
!>>>>> 58:
         write(lut)fslink%latd,fslink%ndd,fslink%tnoofxfr,&
              fslink%tnoofyfr,fslink%totdis,fslink%varreslink,fslink%id
         write(lut)fslink%nooffr,fslink%splink
         write(lut)fslink%dsites
         write(lut)fslink%y2x
         write(lut)fslink%dxidyj
      else
! no disordered fraction set record
!>>>>> 57B:
         write(lut)0
      endif fsrec
!>>>>> 59:
      write(lut)firstvarres%amfu,firstvarres%netcharge,firstvarres%dgm,&
           firstvarres%nprop
! only G values saved ???? well maybe not even those ...
!>>>>> 60:
      write(lut)(firstvarres%gval(i,1),i=1,6)
      do k=1,mc
         write(lut)(firstvarres%dgval(i,k,1),i=1,3)
      enddo
!>>>>> 61:
      write(lut)(firstvarres%d2gval(i,1),i=1,mc2)
   enddo compset
!---- conditions, write as text and recreate when reading file
   call get_all_conditions(text,0,ceq)
   if(gx%bmperr.ne.0) goto 1000
   kl=index(text,'CRLF')
!>>>>> 62:
   write(lut)kl-1
   if(kl.gt.1) then
      write(lut)' SET CONDITIONS ',text(1:kl-1)
   endif
!---- experiments
   call get_all_conditions(text,1,ceq)
   if(gx%bmperr.ne.0) goto 1000
   kl=len_trim(text)
!>>>>> 63:
   write(lut)kl-1
   if(kl.gt.1) then
      write(lut)' EXPERIMENTS    ',text(1:kl-1)
   endif
!>>>>>> 64: savesysmat
! NOTE:: ceq%sysmatdim negative, not initiallized??
! NOTE:: phasetuples not saved !!!
   write(lut)ceq%sysmatdim,ceq%nfixmu,ceq%nfixph
   if(ceq%nfixmu.gt.0) write(lut)(ceq%fixmu(kl),kl=1,ceq%nfixmu)
   if(ceq%nfixph.gt.0) write(lut)&
        (ceq%fixph(1,kl),ceq%fixph(2,kl),kl=1,ceq%nfixph)
   if(ceq%sysmatdim.gt.0) then
      do mc=1,ceq%sysmatdim
         write(lut)(ceq%savesysmat(mc,kl),kl=1,ceq%sysmatdim)
      enddo
   endif
1000 continue
   return
 end subroutine saveequil

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine svfunsave(lut,ceq)
! saves all state variable functions on a file
   implicit none
   integer lut
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character text*512,symbols(20)*32,afterdot*32
   integer ip,ipos,istv,js,jt,kl,ks,lrot
   type(gtp_state_variable), pointer :: svrrec
   write(lut)nsvfun
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
! the 1:10 was a new bug discovered in GNU fortran 4.7 and later
         call make_stvrec(svrrec,svflista(lrot)%formal_arguments(1:10,jt))
!         do ii=1,4
!            indices(ii)=svflista(lrot)%formal_arguments(1+ii,jt)
!         enddo
!         call encode_state_variable2(symbols(js),ip,istv,indices,&
!              svflista(lrot)%formal_arguments(6,jt), &
!              svflista(lrot)%formal_arguments(7,jt),ceq)
         call encode_state_variable(symbols(js),ip,svrrec,ceq)
         if(svflista(lrot)%formal_arguments(10,jt).ne.0) then
! a derivative!!!
            jt=jt+1
            afterdot=' '
            ip=1
            write(*,*)'What? Derivatives not implemented'
!            call encode_state_variable2(afterdot,ip,&
!                 svflista(lrot)%formal_arguments(1,jt),indices,&
!                 svflista(lrot)%formal_arguments(6,jt), &
!                 svflista(lrot)%formal_arguments(7,jt),ceq)
!            symbols(js)=symbols(js)(1:len_trim(symbols(js)))//'.'//afterdot
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
      write(lut)ipos-1
      write(lut)text(1:ipos-1)
   enddo
1000 continue
   return
 end subroutine svfunsave

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine bibliosave(lut)
! saves references on a file
   implicit none
   integer lut
!\end{verbatim}
   character longline*2048
   integer ir,jp,ll,nl
!>>>>> 40:
!   write(*,*)'Saving reference version and number of:',&
!        gtp_biblioref_version,reffree-1
   write(lut)gtp_biblioref_version,reffree-1
   do ir=1,reffree-1
      longline=bibrefs(ir)%reference
      jp=17
      nl=size(bibrefs(ir)%refspec)
      do ll=1,nl
         longline(jp:)=bibrefs(ir)%refspec(ll)
         jp=jp+64
      enddo
      jp=len_trim(longline)
!>>>>> 41:
      write(lut)jp
!>>>>> 42:
      write(lut)longline(1:jp)
   enddo
1000 continue
   return
 end subroutine bibliosave

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine gtpread(filename,str)
! read unformatted all data in the following order
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
   character id*40,endoffile*16,version*8,comment*72,mark*8
   integer i,i1,i2,i3,isp,jph,kontroll,nel,ivers,lin
10  format(i8)
   if(index(filename,'.').eq.0) then
      filename(len_trim(filename)+1:)='.ocu'
   endif
   kontroll=175638
   lin=21
   open(lin,file=filename,access='sequential',status='old',&
        form='unformatted',iostat=gx%bmperr,err=1100)
   if(ocv()) write(*,*)'Opening file: ',filename(1:len_trim(filename)),&
        ' for unformatted read'
!>>>>> 1: read some identification etc, SAVE VERBOSE option!
   if(ocv()) then
      i1=1
   else
      i1=0
   endif
   read(lin)id,version,comment,globaldata
   if(i1.eq.1) then
      globaldata%status=ibset(globaldata%status,GSVERBOSE)
   endif
   if(version.ne.savefile) then
      write(*,11)id,version,savefile
11     format('File not same version as program: ',A/a,' : ',a)
      gx%bmperr=2901; goto 1000
   endif
!   write(*,*)'comment: ',comment(1:len_trim(comment))
   str=comment
   read(lin)noofel,noofsp,noofph,nooftuples
   if(ocv()) write(*,*)'4 numbers: ',noofel,noofsp,noofph,nooftuples
!-------
   read(lin)i1,i2,mark
   if(i1.ne.1 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 1'
      gx%bmperr=4165; goto 1000
   elseif(ocv()) then
      write(*,*)'Control 1 OK'
   endif
!>>>>> 2: elementlist
   read(lin)ivers
   if(ivers.ne.gtp_element_version) then
      write(*,17)'Element',ivers,gtp_element_version
17    format(a,' record version error: ',2i4)
      gx%bmperr=7777; goto 1000
   endif
   do i=1,noofel
      read(lin)ellista(i)
   enddo
   do i=1,noofel
      elements(ellista(i)%alphaindex)=i
   enddo
   if(ocv()) write(*,19)(ellista(i)%alphaindex,i=1,noofel)
19 format('Ellista: ',100i3)
!-------
   read(lin)i1,i2,mark
   if(i1.ne.2 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 2'
      gx%bmperr=4165; goto 1000
   elseif(ocv()) then
      write(*,*)'Control 2 OK'
   endif
!>>>>> 3: specieslist
   read(lin)ivers
   if(ivers.ne.gtp_species_version) then
      write(*,17)'Species',ivers,gtp_species_version
      gx%bmperr=7777; goto 1000
   endif
   do isp=1,noofsp
      read(lin)splista(isp)%symbol,splista(isp)%mass,splista(isp)%charge
      read(lin)splista(isp)%noofel,splista(isp)%status, &
           splista(isp)%alphaindex
      if(isp.gt.1) then
         nel=splista(isp)%noofel
         allocate(splista(isp)%ellinks(nel))
         allocate(splista(isp)%stoichiometry(nel))
      endif
      read(lin)(splista(isp)%ellinks(i),i=1,splista(isp)%noofel)
      read(lin)(splista(isp)%stoichiometry(i),i=1,splista(isp)%noofel)
   enddo
   do i=1,noofsp
      species(splista(i)%alphaindex)=i
   enddo
!   write(*,22)(splista(i)%alphaindex,i=1,noofsp)
!22 format('3E splista: ',20i3)
!-------
   read(lin)i1,i2,mark
   if(i1.ne.3 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 3'
      gx%bmperr=4165; goto 1000
   elseif(ocv()) then
      write(*,*)'Control 3 OK'
   endif
!>>>>> 5: phaselist, starting from 0, the reference phase
   read(lin)ivers,i1,i2,i3
   if(ivers.ne.gtp_phase_version) then
      write(*,17)'Phase',ivers,gtp_phase_version
      gx%bmperr=7777; goto 1000
   endif
   if(i1.ne.gtp_endmember_version) then
      write(*,17)'Endmember',i1,gtp_endmember_version
      gx%bmperr=7777; goto 1000
   endif
   if(i2.ne.gtp_interaction_version) then
      write(*,17)'Interaction',i2,gtp_interaction_version
      gx%bmperr=7777; goto 1000
   endif
   if(i3.ne.gtp_property_version) then
      write(*,17)'Property',i3,gtp_property_version
      gx%bmperr=7777; goto 1000
   endif
   noofem=0
   noofint=0
   noofprop=0
   if(noofph.gt.0) then
      do jph=0,noofph
!>>>>> 5..12 inside readphase
         call readphase(lin,jph)
         if(gx%bmperr.ne.0) goto 1000
         if(ocv()) write(*,*)'Done reading phase: ',jph,' out of: ',noofph
      enddo
      do i=1,noofph
         phases(phlista(i)%alphaindex)=i
      enddo
   endif
   read(lin)(phasetuple(i),i=1,nooftuples)
!--------
   read(lin)i1,i2,mark
   if(i1.ne.4 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 4'
      gx%bmperr=4165; goto 1000
   elseif(ocv()) then
      write(*,*)'Control 4 OK'
   endif
!---------- tpfuns
!>>>>> 20.. inside tpfunread, skip functions already read
   call tpfunread(lin,.TRUE.)
!   write(*,*)'return with error code: ',gx%bmperr
   if(gx%bmperr.ne.0) then
! many functions already entered when reading parameters
      write(*,*)'Error reading TP functiona: ',gx%bmperr
      goto 1000
   endif
!--------
   read(lin)i1,i2,mark
   if(i1.ne.5 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 5'
      gx%bmperr=4165; goto 1000
   elseif(ocv()) then
      write(*,*)'read TPFUNS OK'
   endif
!---------- state variable functions
!>>>>> 30... inside svfunread
   call svfunread(lin)
   if(gx%bmperr.ne.0) goto 1000
!--------
   read(lin)i1,i2,mark
   if(i1.ne.6 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 6'
      gx%bmperr=4165; goto 1000
   elseif(ocv()) then
      write(*,*)'read state variable functions OK at mark ',i1,i2,mark
   endif
!---------- bibliographic references
!>>>>> 40.. inside refread
   call biblioread(lin)
   if(gx%bmperr.ne.0) goto 1000
!--------
   read(lin)i1,i2,mark
   if(i1.ne.7 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 7'
      gx%bmperr=4165; goto 1000
   elseif(ocv()) then
      write(*,*)'read references OK'
   endif
!---------- equilibrium record
!>>>>> 50.. inside readequil
   read(lin)i1,i2,i3
   if(i1.ne.gtp_equilibrium_data_version) then
      write(*,*)'Wrong version of equilibrium data record: ',i1,&
           gtp_equilibrium_data_version
      gx%bmperr=7777; goto 1000
   endif
   if(i2.ne.gtp_component_version) then
      write(*,*)'Wrong version of component record: ',i2,&
           gtp_component_version
      gx%bmperr=7777; goto 1000
   endif
   if(i3.ne.gtp_phase_varres_version) then
      write(*,*)'Wrong version of phase_varres record: ',i3,&
           gtp_phase_varres_version
      gx%bmperr=7777; goto 1000
   endif
   call readequil(lin,firsteq)
   if(gx%bmperr.ne.0) goto 900
!   if(gx%bmperr.ne.0) goto 1000
!--------
   read(lin)i1,i2,mark
   if(i1.ne.8 .and. i2.ne.kontroll) then
      write(*,*)'Read error at control 8'
      gx%bmperr=4165; goto 1000
   elseif(ocv()) then
      write(*,*)'read equilibrium records OK'
   endif
!------ read all ??
   endoffile=' '
   read(lin,end=800,err=800)endoffile
800 continue
   if(endoffile.ne.'- END OF DATA - ') then
      write(kou,811)endoffile
811   format('Unexpected end of file mark: '/'>',A,'<')
      gx%bmperr=4166; goto 1000
   elseif(ocv()) then
      write(kou,812)endoffile
812    format('Expected end of file mark found: '/'>',A,'<')
   endif
! emergency exit
900 continue
   close(lin)
!
1000 continue
   return
! error opening files
1100 continue
   write(*,1110)gx%bmperr,filename(1:len_trim(filename))
1110 format('I/O error: ',i5,', opening file; ',a)
   goto 1000
 end subroutine gtpread

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readphase(lin,jdum)
! read data for phlista and all endmembers etc
! works for test case without disordered fraction test
   implicit none
   integer lin,jdum
!\end{verbatim}
   integer firstendmem,i,i1,i2,i3,jph,level,nem,noi,nop,nox,nup,nsl,mult
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
   if(ocv()) write(*,*)'in readphase:'
! as the phlista record contain pointers each item must be read separately
!>>>>> 5: phase header
   read(lin)jph,phlista(jph)%name,&
        phlista(jph)%models,phlista(jph)%phletter,phlista(jph)%status1,&
        phlista(jph)%alphaindex,phlista(jph)%noofcs,phlista(jph)%nooffs
!>>>>> 6: sublattice info
   read(lin)phlista(jph)%noofsubl,phlista(jph)%linktocs,phlista(jph)%tnooffr
   nsl=phlista(jph)%noofsubl
   allocate(phlista(jph)%nooffr(nsl))
   allocate(phlista(jph)%constitlist(phlista(jph)%tnooffr))
   read(lin)(phlista(jph)%nooffr(i),i=1,nsl),&
        (phlista(jph)%constitlist(i),i=1,phlista(jph)%tnooffr),nem
!------ endmember records, these must be allocated and linked now
   nullify(phlista(jph)%ordered)
   nullify(phlista(jph)%disordered)
   nullify(emrec)
   if(associated(emrec)) then
      write(*,*)'nullify does not work'
      stop
   endif
   if(ocv()) write(*,*)'read endmember data',nsl,nem
! if nem=0 now there are no basic (ordered) endmember (can that happen?)
! return here when endmember list empty and there is a disordered list
   firstendmem=1
200 continue
   if(ocv()) write(*,202)'reading parameters: ',phlista(jph)%name,&
        nsl,firstendmem,nem
202 format(a,a,10i4)
!   newendmem: do while(nem.eq.1)
   newendmem: do while(nem.gt.0)
      if(associated(emrec)) then
!>>>>> 7C: the second or later endmember in the same list
         call readendmem(lin,nsl,emrec%nextem,nop,noi,nem)
         emrec=>emrec%nextem
      elseif(firstendmem.eq.1) then
!>>>>> 7A: the first (only or ordered) endmember
         call readendmem(lin,nsl,phlista(jph)%ordered,nop,noi,nem)
         emrec=>phlista(jph)%ordered
      elseif(firstendmem.eq.2) then
!>>>>> 7B: the first disordered endmember
         if(ocv()) write(*,*)'Reading isordered parameter list'
         call readendmem(lin,nsl,phlista(jph)%disordered,nop,noi,nem)
         emrec=>phlista(jph)%disordered 
         firstendmem=0
      endif
      if(nop.eq.1) then
!>>>>> 8A: endmember property (lookp)
         call readproprec(lin,emrec%propointer,nox)
         proprec=>emrec%propointer
         do while(nox.eq.1)
            call readproprec(lin,proprec%nextpr,nox)
            proprec=>proprec%nextpr
         enddo
      endif
      inttree: if(noi.eq.1) then
!>>>>> 9A: interaction record
         level=0
         call readintrec(lin,emrec%intpointer,mult,noi,nup,nop)
         intrec=>emrec%intpointer
         if(ocv()) write(*,13)'read interaction: ',intrec%status,noi,nup,nop
13       format(a,10i4)
300      continue
         if(nop.eq.1) then
!>>>>> 10A: interaction property record
            call readproprec(lin,intrec%propointer,nox)
            proprec=>intrec%propointer
            do while(nox.eq.1)
               call readproprec(lin,proprec%nextpr,nox)
               proprec=>proprec%nextpr
            enddo
         endif
! push before going to higher
330      continue
         level=level+1
         stack(level)%p1=>intrec
         stack(level)%noi=noi
         if(ocv()) write(*,13)'pushed interaction: ',intrec%status,0,0,0,level
         higher: if(nup.eq.1) then
!>>>>> 9B: go to higher level and save intrec
            call readintrec(lin,intrec%highlink,mult,noi,nup,nop)
            intrec=>intrec%highlink
!            write(*,13)'read higher interaction: ',intrec%status,noi,nup,nop
            if(nop.eq.1) then
!>>>>> 10B: there are some property records !!
               call readproprec(lin,intrec%propointer,nox)
               proprec=>intrec%propointer
               do while(nox.eq.1)
                  call readproprec(lin,proprec%nextpr,nox)
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
            if(ocv())write(*,13)'poped interaction: ',intrec%status,0,0,0,level
            if(noi.eq.1) then
!>>>>> 9C: 
               call readintrec(lin,intrec%nextlink,mult,noi,nup,nop)
               intrec=>intrec%nextlink
               if(ocv()) write(*,13)'read interaction: ',intrec%status,&
                    noi,nup,nop
               goto 300
            else
               goto 350
            endif
         endif pop
      endif inttree
   enddo newendmem
! we come nere when no more endmembers in this list
   if(firstendmem.eq.1) then
!>>>>> 11: if nem read here is zero there are no disordered endmembers
      if(ocv()) write(*,*)'checking for disordered endmembers'
      read(lin)nem,nsl
! we must nullify emrec to start a new list of endmembers
      nullify(emrec)
      if(nem.ne.0) then
         firstendmem=2
         if(ocv()) write(*,*)'Reading disordered parameters',nem,nsl
         goto 200
      endif
   endif
!------ additions list
!500 continue
   nullify(phlista(jph)%additions)
510 continue
   read(lin)i1,i2,i3
   if(ocv()) write(*,*)'Reading any addition; ',i1
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
   elseif(i1.eq.-1) then
! end of addition list
      continue
      if(i2.ne.i1 .and. i3.ne.i1) write(*,*)'end of phase error:',i2,i3
   endif
1000 continue
   return
 end subroutine readphase

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readendmem(lin,nsl,emrec,nop,noi,nem)
! allocates and reads an endmember record
   implicit none
   integer lin,nsl,nop,noi,nem
   type(gtp_endmember), pointer :: emrec
!\end{verbatim}
   integer i,j
   allocate(emrec)
!   write(*,*)'Going to read endmember record'
!>>>>> 7D: actually reading ....
   read(lin)emrec%noofpermut,emrec%phaselink,emrec%antalem,nop,noi,nem
   if(ocv()) write(*,17)'readendmem: ',nsl,emrec%noofpermut,emrec%phaselink,&
        emrec%antalem,nop,noi,nem
17 format(a,7i5)
   allocate(emrec%fraclinks(nsl,emrec%noofpermut))
   do j=1,emrec%noofpermut
      read(lin)(emrec%fraclinks(i,j),i=1,nsl)
   enddo
   nullify(emrec%nextem)
   nullify(emrec%propointer)
   nullify(emrec%intpointer)
   noofem=noofem+1
   return
 end subroutine readendmem

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readproprec(lin,proprec,nox)
! allocates and reads a property record
   implicit none
   integer lin,nox
   type(gtp_property), pointer :: proprec
!\end{verbatim}
   integer i
   allocate(proprec)
!>>>>> 8B: actually reading property record (endmember)
!>>>>> 10B: actually reading property record (interaction)
!   write(*,*)'Going to read property record'
   read(lin)proprec%reference,proprec%proptype,&
        proprec%degree,proprec%extra,proprec%antalprop,nox
   allocate(proprec%degreelink(0:proprec%degree))
   if(ocv()) write(*,17)'readprop: ',proprec%proptype,proprec%degree,&
        proprec%antalprop,nox
17 format(a,6i5)
!   write(*,*)'To read TP functions: ',proprec%degree,proprec%degreelink(0)
   do i=0,proprec%degree
      call read1tpfun(21,proprec%degreelink(i))
   enddo
   if(ocv()) write(*,*)'Read TP functions: ',proprec%degree,&
        proprec%degreelink(0)
   nullify(proprec%nextpr)
   noofprop=noofprop+1
   return
 end subroutine readproprec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readintrec(lin,intrec,mult,noi,nup,nop)
! allocates and reads an interaction record UNFINISHED
   implicit none
   integer lin,mult,noi,nup,nop
   type(gtp_interaction), pointer :: intrec
!\end{verbatim}
   integer fipsize,noofperm,i
! the storage of permutations in interaction records is complex ... one must
! take into account the number of permutations in lower order intecations ...
! for an fcc endmember A:A:A:B (4 perm) the binary interaction A:A:A,B:B has 
! 3; 3; 3 and 3 perms and the ternary A:A,B:A,B:B has 2; 2; 2; 2
! mult may not be needed ...
   allocate(intrec)
!>>>>> 9D: actually read the interaction record
   read(lin)fipsize
   allocate(intrec%noofip(fipsize))
   read(lin)intrec%noofip,intrec%status,noi,nup,nop
!   write(*,17)'3E readint: ',fipsize,intrec%status,noi,nup,nop,&
!        (intrec%noofip,i=1,fipsize)
17 format(a,5i4,2x,10i3)
   noofperm=intrec%noofip(2)
   allocate(intrec%sublattice(noofperm))
   allocate(intrec%fraclink(noofperm))
   do i=1,intrec%noofip(1)
      read(lin)intrec%sublattice(i),intrec%fraclink(i)
   enddo
   nullify(intrec%nextlink)
   nullify(intrec%highlink)
   nullify(intrec%propointer)
   return
 end subroutine readintrec

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine readequil(lin,ceq)
! Read equilibria records from a file
   implicit none
   integer lin
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character text*512,dum16*16
   type(gtp_phase_varres), pointer :: firstvarres
   TYPE(gtp_fraction_set) :: fslink
   integer i,ierr,ip,isp,ivar,j,jp,k,lokcs,lokph,mc,mc2,nprop,nsl,kp
   double precision, dimension(:,:), allocatable :: ca,ci
! containing conditions, components and phase varres records for wach compset
!>>>>> 50:
   read(lin)ceq%eqname,ceq%eqno,ceq%status,ceq%next
   if(ocv()) write(*,*)'Reading equilibrium: ',ceq%eqname
!----- components
!   allocate(ceq%complist(noofel)) already allocated for 20
!   write(*,*)'Size of component array: ',size(ceq%complist)
!>>>>> 51:
   do i=1,noofel
      read(lin)isp
      ceq%complist(i)%splink=isp
      read(lin)ceq%complist(i)%phlink,&
           ceq%complist(i)%status,&
           ceq%complist(i)%refstate,ceq%complist(i)%tpref,ceq%complist(1)%mass
! ============ never written ???
!      if(isp.gt.0) then
! user defined reference state only allocated if necessary
!         read(lin)j
!         allocate(ceq%complist(i)%endmember(j))
!         read(lin)ceq%complist(i)%endmember
!         read(lin)ceq%complist(i)%molat
!      endif
   enddo
!==================================
! stoichiometry conversion matrix, already allocated?? where??
   ceq%compstoi=zero
! calculate the inverse stoichiometry matrix
!   write(*,*)'Reading component stoichiometry matrix',noofel,maxel
   do j=1,noofel
      read(lin)(ceq%compstoi(j,i),i=1,noofel)
   enddo 
 ! this because mdinv did strange things inverting a larger matrix
   allocate(ca(noofel,noofel+1))
   allocate(ci(noofel,noofel))
   do i=1,noofel
      do j=1,noofel
         ca(i,j)=ceq%compstoi(i,j)
      enddo
   enddo
!   do j=1,noofel
!      write(*,99)'ca: ',(ca(j,i),i=1,noofel)
!   enddo 
99 format(a,7e11.3)
!   call mdinv(maxel-1,maxel,ceq%compstoi,ceq%invcompstoi,noofel,ierr)
   call mdinv(noofel,noofel+1,ca,ci,noofel,ierr)
!   write(*,*)'Inverting matrix',ierr
   do j=1,noofel
!      write(*,99)'ci: ',(ci(i,j),i=1,noofel)
      do i=1,noofel
         ceq%invcompstoi(i,j)=ci(i,j)
      enddo
   enddo
   deallocate(ca)
   deallocate(ci)
!----------- phase_varres record
!>>>>> 54:
   read(lin)highcs
   if(ocv()) then
      write(*,*)'Number of phase_varres records: ',highcs
      write(*,*)'phase_varres size: ',size(ceq%phase_varres)
   endif
   do j=1,highcs
!      write(*,*)'reading phase_varres ',j
!------------------------------------------
! DEBUGPROBLEM BEWARE, using = instead of => below took 2 days to find
!------------------------------------------
! >>>      firstvarres=ceq%phase_varres(j)    <<< error
      firstvarres=>ceq%phase_varres(j)
!>>>>> 55:
      read(lin)firstvarres%nextfree,firstvarres%phlink,&
           firstvarres%status2,firstvarres%phstate
      read(lin)firstvarres%prefix,firstvarres%suffix
      read(lin)firstvarres%abnorm
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
!      write(*,*)'Allocate constat 1: ',nsl,mc
      allocate(firstvarres%constat(mc))
      allocate(firstvarres%yfr(mc))
      allocate(firstvarres%sites(nsl))
! for ionic liquids allocate dpqdy
      if(btest(phlista(firstvarres%phlink)%status1,PHIONLIQ)) then
         if(ocv()) write(*,*)'Allocate dpqdy: ',mc
         allocate(firstvarres%dpqdy(mc))
      endif
!>>>>> 56:
      read(lin)(firstvarres%constat(i),i=1,mc)
      read(lin)(firstvarres%yfr(i),i=1,mc)
      allocate(firstvarres%mmyfr(mc))
      read(lin)(firstvarres%mmyfr(i),i=1,mc)
      read(lin)(firstvarres%sites(i),i=1,nsl)
! these are ignored, values not important but must be allocated
! for ionic liquid as (2,mc) and (2,mc2)
!       read(lin)(firstvarres%dsitesdy(2,i),i=1,mc)
!       read(lin)(firstvarres%d2sitesdy2(2,i),i=1,mc2)
!>>>>> 57:
      read(lin)ivar
      if(ivar.eq.1) then
! extra fraction set
         if(ocv()) write(*,*)'reading extra fraction set for ',j
!>>>>> 58:
         read(lin)fslink%latd,fslink%ndd,fslink%tnoofxfr,&
              fslink%tnoofyfr,fslink%totdis,fslink%varreslink,fslink%id
         allocate(fslink%nooffr(fslink%ndd))
         allocate(fslink%dsites(fslink%ndd))
         allocate(fslink%splink(fslink%tnoofxfr))
         allocate(fslink%y2x(fslink%tnoofyfr))
         allocate(fslink%dxidyj(fslink%tnoofyfr))
         read(lin)fslink%nooffr,fslink%splink
         read(lin)fslink%dsites
         read(lin)fslink%y2x
         read(lin)fslink%dxidyj
! now copy fslink to the correct record and then deallocate fslink arrays
         call copy_fracset_record(j,fslink,ceq)
         deallocate(fslink%nooffr)
         deallocate(fslink%dsites)
         deallocate(fslink%splink)
         deallocate(fslink%y2x)
         deallocate(fslink%dxidyj)
      endif
! The result data
!>>>>> 59:
      read(lin)firstvarres%amfu,firstvarres%netcharge,firstvarres%dgm, &
           firstvarres%nprop
      nprop=firstvarres%nprop
      allocate(firstvarres%listprop(nprop))
      allocate(firstvarres%gval(6,nprop))
      allocate(firstvarres%dgval(3,mc,nprop))
      allocate(firstvarres%d2gval(mc2,nprop))
!>>>>> 60:
      read(lin)(firstvarres%gval(i,1),i=1,6)
      do k=1,mc
         read(21)(firstvarres%dgval(i,k,1),i=1,3)
      enddo
!>>>>> 61:
      read(lin)(firstvarres%d2gval(i,1),i=1,mc2)
      if(ocv()) write(*,*)'phase_varres size: ',j,size(ceq%phase_varres)
   enddo
!----- conditions, can be empty, NOTE: entered after phase_varres
!>>>>> 62:
   read(lin)ip
   if(ip.gt.0) then
      read(lin)dum16,text(1:ip)
! set the conditions, ip will be incremented by 1 in enter_condition
! the text contains " number: variable=value, "
! we have to set each condition variable separately
      jp=1
      if(ocv()) write(*,*)'Conditions >',text(1:ip),'<',jp,ip
      cloop: do while(jp.lt.ip)
         k=index(text(jp:ip),':')
         if(k.le.0) exit cloop
         jp=jp+k
         kp=min(jp+index(text(jp+1:),' '),ip)
         if(kp.gt.jp) then
! remove , as that indicates more conditions on same line
            if(text(kp-1:kp-1).eq.',') then
               text(kp-1:kp-1)=' '
            endif
         else
            kp=ip
         endif
         jp=jp-1
!         write(*,*)'condition 1 >',text(jp:kp),'<',jp
         call set_condition(text(1:kp),jp,firsteq)
! jp automatically update inside set_condition
         if(gx%bmperr.ne.0) then
            write(*,*)'Error setting conditions'
            write(*,*)ip,' >',text(jp:kp),'<'
            goto 1000
         endif
      enddo cloop
   endif
!---- experiments
!>>>>> 63:
   read(lin)ip
   if(ip.gt.0) then
      read(lin)dum16,text(1:ip)
   endif
!>>>>>> 64: restore savesysmat
   read(lin)ceq%sysmatdim,ceq%nfixmu,ceq%nfixph
!   write(*,*)'savesysmat: ',ceq%sysmatdim,ceq%nfixmu,ceq%nfixph
   if(ceq%nfixmu.gt.0) then
      allocate(ceq%fixmu(ceq%nfixmu))
      read(lin)(ceq%fixmu(kp),kp=1,ceq%nfixmu)
   endif
   if(ceq%nfixph.gt.0) then
      allocate(ceq%fixph(2,ceq%nfixph))
      read(lin)(ceq%fixph(1,kp),ceq%fixph(2,kp),kp=1,ceq%nfixph)
   endif
   if(ceq%sysmatdim.gt.0) then
      allocate(ceq%savesysmat(ceq%sysmatdim,ceq%sysmatdim))
      do mc=1,ceq%sysmatdim
         read(lin)(ceq%savesysmat(mc,kp),kp=1,ceq%sysmatdim)
      enddo
   endif
! allocate and zero the array with current chemical potentials
   if(.not.allocated(ceq%cmuval)) allocate(ceq%cmuval(noofel))
   ceq%cmuval=zero
!
   write(*,*)'UNIFINSHED reading of equilibria ...',highcs
! what about the free list ....??? there can be free records below highcs ...
   csfree=highcs+1
1000 continue
   return
 end subroutine readequil

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine svfunread(lin)
! read a state variable function from save file and store it.
! by default there are some state variable functions, make sure
! they are deleted.  Done here just by setting nsvfun=0
   implicit none
   integer lin
!\end{verbatim}
   integer nsvfun,i,ip,nsvfunfil
   character*512 text
   nsvfun=0
   read(lin)nsvfunfil
   if(ocv()) write(*,*)'Number of state variable functions: ',nsvfunfil
   do i=1,nsvfunfil
      read(lin)ip
!      write(*,*)'Number of characters: ',ip
      text=' '
      read(lin)text(2:ip)
!      write(*,*)text(2:ip)
      ip=1
      call enter_svfun(text,ip,firsteq)
      if(gx%bmperr.ne.0) then
!         write(*,*)'Error entering svf from file',gx%bmperr
         if(gx%bmperr.ne.4136) goto 1000
         gx%bmperr=0
      endif
   enddo
1000 continue
   return
 end subroutine svfunread

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine biblioread(lin)
! read references from save file
   implicit none
   integer lin
!\end{verbatim}
   character text*512
   integer i,iref,jp,nrefs
!>>>>> 40: number of references
!   write(*,*)'Reading reference version and nummer of'
   read(lin)i,nrefs
!   write(*,*)i,nrefs
   if(i.ne.gtp_biblioref_version) then
      write(*,*)'Warning, the bibliographic references version not same'
   endif
   if(ocv()) write(*,*)'Reading bibligraphic references: ',nrefs,reffree
!   reffree=nrefs+1
   reffree=1
   do i=1,nrefs
!>>>>> 41: number characters to read
      read(lin)jp
!      write(*,*)'Length of text: ',i,jp
!>>>>> 42: text
      if(jp.gt.512) then
         write(*,*)'Too long bibliographic reference text',jp
         gx%bmperr=7777; goto 1000
      endif
      read(lin)text(1:jp)
!      write(*,*)text(1:jp)
      call tdbrefs(text(1:16),text(17:jp),0,iref)
      if(gx%bmperr.ne.0) goto 1000
   enddo
1000 continue
   return
 end subroutine biblioread

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

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
   if(ocv()) write(*,*)'Removing current data'
!---------- elementlist, no need to delete, just deallocate below
!>>>>> 2:
!---------- specieslist, we have to deallocate ?? maybe not ??
!>>>>> 3:
   if(btest(globaldata%status,GSNODATA)) then
      if(ocv()) write(*,*)'No thermodynamic data to delete'
      goto 600
   endif
   if(gtp_species_version.ne.1) then
      if(ocv()) write(*,17)'Species',1,gtp_species_version
17    format(a,' record version error: ',2i4)
      gx%bmperr=7777; goto 1000
   endif
   ceq=>firsteq
!   write(*,*)'3E No segmentation error A'
   do isp=1,noofsp
      nel=splista(isp)%noofel
      deallocate(splista(isp)%ellinks)
      deallocate(splista(isp)%stoichiometry)
   enddo
!---------- phases, many records, here we travese all endmembers etc
!>>>>> 4
!   write(*,*)'3E No segmentation error B'
   if(gtp_phase_version.ne.1) then
      if(ocv()) write(*,17)'Phase',1,gtp_phase_version
      gx%bmperr=7777; goto 1000
   endif
   if(gtp_endmember_version.ne.1) then
      if(ocv()) write(*,17)'Endmember',1,gtp_endmember_version
      gx%bmperr=7777; goto 1000
   endif
   if(gtp_interaction_version.ne.1) then
      if(ocv()) write(*,17)'Interaction',1,gtp_interaction_version
      gx%bmperr=7777; goto 1000
   endif
   if(gtp_property_version.ne.1) then
      if(ocv()) write(*,17)'Property',1,gtp_property_version
      gx%bmperr=7777; goto 1000
   endif
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
      do k=1,size(ceq%phase_varres)
         phdyn=>ceq%phase_varres(k)
         if(allocated(phdyn%gval)) then
            deallocate(phdyn%gval)
            deallocate(phdyn%dgval)
            deallocate(phdyn%d2gval)
!            write(*,*)'3E No segmentation error C7',j,k
         endif
      enddo
!      write(*,*)'3E No segmentation error C8',j
!      deallocate(ceq%phase_varres)
   enddo
!   write(*,*)'3E No segmentation error D1'
! I am not sure if this really releases all memory, how to check .... ???
   deallocate(eqlista)
!   write(*,*)'3E No segmentation error D2'
!------- deallocate elements, species and phases, will be allocated in init_gtp
   deallocate(ellista)
   deallocate(elements)
   deallocate(splista)
   deallocate(species)
   deallocate(phlista)
   deallocate(phases)
   deallocate(phasetuple)
!   write(*,*)'3E No segmentation error E'
!------ tpfunction expressions and other lists
!>>>>> 20: delete tpfuns
!   write(*,*)'Delete TP funs, just deallocate??'
   call delete_all_tpfuns
!   write(*,*)'Back from deleting all TP funs, this is fun!!'
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
!   write(*,*)'3E No segmentation error G'
!    deallocate( .... any more ???
!---------------------------
! now initiate all lists and a little more
   if(ocv()) write(*,*)'All data structures will be reinitiated'
! intv(1) negative means reinititate with same values as before
   intv(1)=-1
!   write(*,*)'3E No segmentation error H'
   call init_gtp(intv,dblv)
! after return firsteq must ve initiated ... maybe it should be done here ??
!
1000 continue
   return
 end subroutine new_gtp

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

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
!   write(*,*)'In delphase',lokph
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
!      write(*,*)'deallocate endmember record'
      deallocate(emrec)
! nextem do not need to be declared as target??
      emrec=>nextem
      emproplista: do while(associated(proprec))
         nextprop=>proprec%nextpr
!>>>>> 8: endmember property records
! functions and references deallocated separately
!         write(*,*)'deallocate endmember property record'
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
!         write(*,*)'Pushing ',level
         stack(level)%p1=>intrec%nextlink
         nextint=>intrec%highlink
         proprec=>intrec%propointer
!         write(*,*)'deallocate interaction record'
         deallocate(intrec)
         intproplista: do while(associated(proprec))
            nextprop=>proprec%nextpr
!>>>>> 10: interaction properties
!            write(*,*)'deallocate interaction property record'
            deallocate(proprec)
            proprec=>nextprop
         enddo intproplista
         intrec=>nextint
      enddo intlista
! pop the link to next interaction if any
      pop: if(level.gt.0) then
!         write(*,*)'popping interaction record',level
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
!      write(*,*)'disordered endmembers'
      emrec=>phlista(lokph)%disordered
      noendm=1
      goto 200
   endif
!   write(*,*)'finished parameter records'
!------ additions list
500 continue
   addlink=>phlista(lokph)%additions
   addition: do while(associated(addlink))
!>>>>> 12: additions
      nextadd=>addlink%nextadd
      if(addlink%type.eq.1) then
!>>>>> 12A: delete magnetic addition ...
         deallocate(addlink)
      else
         write(*,*)'Cannot delete unknown addition type ',addlink%type
      endif
      addlink=>nextadd
   enddo addition
!   write(*,*)'phase location: ',lokph,size(phlista(lokph)%nooffr),&
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
   return
 end subroutine delphase

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 logical function iskeyword(text,keyword,nextc)
! compare a text with a given keyword. Abbreviations allowed
! but the keyword and abbreviation must be surrounded by spaces
! nextc set to space character in text after the (abbreviated) keyword
   implicit none
   character text*(*),keyword*(*),key*64
   integer nextc
!\end{verbatim}
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

!\begin{verbatim}
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
   integer, parameter :: nkw=12
   character (len=kwl), dimension(nkw), parameter :: keyword=&
        ['ELEMENT             ','SPECIES             ',&
         'PHASE               ','CONSTITUENT         ',&
         'FUNCTION            ','PARAMETER           ',&
         'TYPE_DEFINITION     ','LIST_OF_REFERENCES  ',&
         'ADD_REFERENCES      ','ASSESSED_SYSTEMS    ',&
         'DATABASE_INFORMATION','VERSION             ']
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
   istdbkeyword=j
   return
 end function istdbkeyword

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

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

!\begin{verbatim}
 subroutine readtdb(filename,nel,selel)
! reading data from a TDB file with selection of elements
!-------------------------------------------------------
! Not all TYPE_DEFS implemented
!-------------------------------------------------------
   implicit none
   integer nel
   character filename*(*),selel(*)*2
!\end{verbatim}
   character line*100,elsym*2,name1*24,name2*24,elsyms(10)*2
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
   integer lint(2,3),noofphasetype,nytypedef,nextc,keyw,tdbv
   integer typty,fractyp,lp1,lp2,ix,jph,kkk,lcs,nint,noelx
   logical onlyfun,nophase,ionliq,notent
   integer norew,newfun,nfail,nooftypedefs,nl,ipp,jp,jss,lrot,ip,jt
   integer nsl,ll,kp,nr,nrr,mode,lokph,lokcs,km,nrefs,ideg,iph,ics
! disparttc and dispartph to handle phases with disordered parts
   integer nofunent,disparttc,dodis,jl,nd1,thisdis
   character*24 dispartph(5),ordpartph(5)
   logical warning
! set to TRUE if element present in database
   logical, allocatable :: present(:)
! to prevent any output
   logical silent
!  mmyfr
! if warning is true at the end pause before listing bibliography
   warning=.FALSE.
   silent=.FALSE.
   if(btest(globaldata%status,GSSILENT)) then
      silent=.TRUE.
!      write(*,*)'3E reading database silent'
   endif
   if(ocv()) write(*,*)'reading a TDB file'
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
! One should remove TAB characters !! ?? YES !!
!   if(line(1:1).eq.'$') goto 100
   ipp=1
   if(eolch(line,ipp)) goto 100
   if(line(ipp:ipp).eq.'$') goto 100
! replace TAB by space
   call replacetab(line,nl)
   goto 120
!----- this part moved to the end ....
   funfirst: if(onlyfun) then
! first read only functions until all has been read
!      if(line(2:10).eq.'FUNCTION ') then
!      write(*,*)'Input line >',line(1:20),'<'
      ipp=istdbkeyword(line,nextc)
      if(ipp.eq.5) then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! FUNCTION GHSERCR    2.98150E+02  -8856.94+157.48*T-26.908*T*LN(T)
!         name1=line(11:18)
! special case, error in TDB file, UN_ASS is only 6 characters
!         if(name1(1:6).eq.'UN_ASS') then
!            name1=line(11:16); ipp=18
!         else
!            ipp=20
!         endif
         if(eolch(line,nextc)) then
            if(.not.silent) &
                 write(kou,*)'Function name must be on same line as FUNCTION'
            gx%bmperr=4000; goto 1000
         endif
         ipp=nextc+index(line(nextc:),' ')
         name1=line(nextc:ipp-1)
!         write(*,18)'function >',name1,'< ',nextc,ipp
!18       format(a,a,a,2i4)
! old code
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
!            write(*,*)'3E Entering function 1: ',name1,len_trim(longline)
!            lrot=0
            call enter_tpfun(name1,longline,lrot,.TRUE.)
            if(gx%bmperr.ne.0) then
! one may have error here if function calls other functions not entered, 4002
! or if the function is already entered, 4026
               if(gx%bmperr.eq.4002.or. gx%bmperr.eq.4026) then
                  if(gx%bmperr.eq.4002) nfail=nfail+1
                  gx%bmperr=0; goto 100
               endif
               if(.not.silent) write(kou,*)'Failed entering function: ',name1
               gx%bmperr=4000
               goto 1000
            endif
            if(ocv()) write(*,*)'Entered function: ',name1
            newfun=newfun+1
         else
            nl=nl+1
            read(21,110)line
!            write(kou,101)'readtdb 2: ',nl,line(1:40)
            call replacetab(line,nl)
            longline=longline(1:jp)//line
            goto 111
         endif
      elseif(ipp.gt.0) then
! skip lines until !.  There can be a ! on the line with the keyword!
77       continue
         if(index(line,'!').le.0) then
            read(21,110,end=2000)line
            nl=nl+1
            call replacetab(line,nl)
            goto 77
         endif
      endif
      goto 100
   endif funfirst
!---------------------------------------------------------
! handle all TDB keywords except function
120 continue
   keyw=istdbkeyword(line,nextc)
   if(keyw.eq.0) then
      ip=1
      if(.not.eolch(line,ip)) then
         if(ocv()) write(*,*)'Ignoring line: ',nl,ip,line(ip:ip+20)
      endif
      goto 100
   elseif(onlyfun) then
      if(keyw.eq.5) goto 800
      goto 100
   endif
   if(.not.nophase .and. keyw.ne.4) then
! after a PHASE keyword one should have a CONSTITUENT
      if(.not.silent) write(kou,*)'expeciting CONSTITUENT: ',line(1:30)
      warning=.TRUE.
   endif
! check there is a ! in line, otherwise read until we find an exclamation mark
   ip=1
   longline(ip:)=line
   ip=len_trim(longline)+1
!   write(*,*)'new keyword ',ip,'>',longline(1:40)
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
            gx%bmperr=7777; goto 1000
         endif
      endif
   enddo
   if(dodis.eq.1) then
! if dodis=1 only read data for disordred phases
! PHASE=3, CONSTITUENT=4, PARAMETER=6 ... any more?      
      if(keyw.lt.3 .or. keyw.eq.5 .or. keyw.gt.6) goto 100
   endif
!
  select case(keyw)
   case default
      if(ocv()) write(*,*)'default case: ',keyw,line(1:30)
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
         gx%bmperr=7777; goto 1000
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
         if(ocv()) write(*,*)'Skipping database element: ',elsym
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
      call enter_element(elsym,name2,name1,mass,h298,s298)
      if(gx%bmperr.ne.0) goto 1000
   case(2) !SPECIES -------------------------------------------------
!   elseif(line(2:9).eq.'SPECIES ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! SPECIES O3PU2                       O3PU2!
      ip=nextc
      if(eolch(longline,ip)) then
         if(.not.silent) write(kou,*)'Line after SPECIES keyword empty'
         gx%bmperr=7777; goto 1000
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
!   elseif(line(2:10).eq.'FUNCTION ') then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! FUNCTION GHSERCR    2.98150E+02  -8856.94+157.48*T-26.908*T*LN(T)
!      name1=line(11:18)
!      longline=' '
!      longline=line(20:)
!300    continue
!      jp=len_trim(longline)
!      if(longline(jp:jp).eq.'!') then
!          write(*,*)'Skipping function: ',name1
! all functions entered at the end, skip until !
!      do while(index(longline,'!').le.0)
      if(index(longline,'!').le.0) then
         if(.not.silent) &
              write(*,*)' Error, terminating ! not found for funtion!!',nl
         gx%bmperr=7777; goto 1000
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
         gx%bmperr=7777; goto 1000
      endif
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
! phytype
      if(jp.gt.0) then
         phtype=name1(jp+1:jp+1)
         name1(jp:)=' '
      else
         phtype=' '
      endif
!      write(*,*)'nophase set to false, phase: ',name1
! phase type code
      noofphasetype=0
      ip=ip+1
      jp=ip
      name2=longline(ip:jp)
      thisdis=0
      phdis: if(dodis.eq.1) then
! special when reading disordered parts, check phase name equail
!         write(*,*)'Check if disordered part: ',dodis,name1
         do jt=1,disparttc
            if(name1.eq.dispartph(jt)) goto 307
         enddo
! not a disordered part
         goto 100
307      continue
         thisdis=jt
!         write(*,*)'Found disordered part: ',name1,thisdis
! we skip the rest of the phase line ...
         goto 100
      elseif(dodis.eq.0 .and. disparttc.gt.0) then
! we must not enter phases that are disordered parts
         do jt=1,disparttc
            if(name1.eq.dispartph(jt)) then
!               write(*,*)'Skip phase that is a disordered part: ',name1
               thisdis=-1
               goto 100
            endif
         enddo
      endif phdis
!      write(*,*)'Entering phase: ',name1
!      write(*,*)'Checking phase types for phase: ',name1,jp
! skip blanks, then read type code, finished by a blank
      if(eolch(longline,jp)) then
         if(.not.silent) &
              write(kou,*)'3E no phase typecode: ',name1(1:len_trim(name1))
         warning=.TRUE.
      endif
      jp=jp-1
310   jp=jp+1
! NOTE and FIX: type code expected to be after a single space: be flexible ??
      typedefcheck: if(longline(jp:jp).ne.' ') then
         ch1=longline(jp:jp)
!         write(*,*)'3E typedef: ',ch1,jp
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
!      write(*,*)'3E buperr: ',buperr ,jp
      call getrel(longline,jp,xsl)
      if(buperr.ne.0) then
         if(.not.silent) write(kou,*)'3E tdb: "',longline(1:jp),'"',buperr
         gx%bmperr=buperr; goto 1000
      endif
      nsl=int(xsl)
      do ll=1,nsl
         call getrel(longline,jp,stoik(ll))
         if(buperr.ne.0) then
            gx%bmperr=buperr; goto 1000
         endif
      enddo
!      write(*,*)'readtdb 3A: ',nsl,(stoik(ll),ll=1,nsl)
!---------------------------------------------------------------------
! The constituent line must follow PHASE before any new phase
   case(4) !    CONSTITUENT LIQUID:L :CR,FE,MO :  !
! the phase must have been defined
      if(nophase) then
         if(.not.silent) write(kou,*) &
              'A CONSTITUENT keyword not directly preceeded by PHASE!'
         gx%bmperr=7777; goto 1000
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
!      write(*,*)'readtdb gas1: ',nl,jp,longline(1:jp)
! eliminate all after the exclamation mark
!      longline(jp+1:)=' '
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
      if(ll.ge.1) then
         knr(ll)=nr
         if(nr.le.0) then
            if(ocv()) then
               write(*,*)'Skipping phase due to missing constituents: ',name1
!              write(*,378)name1,ll
378            format('Phase ',a,' has no constituents in sublattice ',i2)
! Not a fatal error when elements have been selected but skip this phase
!              gx%bmperr=7777; goto 1000
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
         gx%bmperr=7777; goto 1000
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
         if(.not.silent) write(kou,*)'Error extracting constituents 2'
         gx%bmperr=7777; goto 1000
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
      gx%bmperr=7777; goto 1000
!      write(*,*)'Species terminator error: ',ch1,nl
!      gx%bmperr=4157; goto 1000
390    continue
! name2 is model, ignored on reading TDB
      ionliq=.FALSE.
      if(phtype(1:1).eq.'Y') then
         name2='IONIC_LIQUID '
         ionliq=.TRUE.
      else
         name2='CEF-TDB-RKM? '
      endif
      if(ocv()) write(*,*)'readtdb 9: ',name1,nsl,knr(1),knr(2),phtype
395   continue
      condis2: if(dodis.eq.1) then
! if we have a disordered part do not enter the phase, add disordered fracs!
! the ordered phase name is ordpart(thisdis)
         call find_phase_by_name(ordpartph(thisdis),iph,ics)
         if(gx%bmperr.ne.0) then
! NOTE THE ORDERED PHASE MAY NOT BE ENTERED DUE TO COMPONENTS!!
            if(.not.silent) write(kou,396)thisdis,ordpartph(thisdis)
396         format('Disordered phase skipped as no ordered: ',i3,' "',a,'"')
            warning=.TRUE.
            gx%bmperr=0
            goto 100
         else
            if(.not.silent) write(kou,*) &
                 'Adding disordered fraction set: ',ordpartph(thisdis)
         endif
! we are creating the phase, there is only one composition set
         call get_phase_compset(iph,1,lokph,lokcs)
         if(gx%bmperr.ne.0) goto 1000
! ch1 is suffix for parameters, always D
         ch1='D'
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
         if(jl.eq.0 .and. .not.silent) write(kou,398)trim(ordpartph(thisdis))
398      format(' Assuming that phase ',a,' cannot completely disorder')
! add DIS_PART from TDB
         call add_fraction_set(iph,ch1,nd1,jl)
         if(gx%bmperr.ne.0) then
            if(.not.silent) write(kou,*) &
                 '3E Error entering disordered fraction set: ',gx%bmperr
            goto 1000
         endif
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
         if(.not.silent) write(kou,601) &
              dispartph(thisdis)(1:len_trim(dispartph(thisdis))),ch1,nd1,jl,xxx
601      format('Parameters from disordered part added: ',a,5x,a,2x,2i3,F12.4)
      else
         call enter_phase(name1,nsl,knr,const,stoik,name2,phtype)
!      write(*,*)'readtdb 9A: ',gx%bmperr
         if(gx%bmperr.ne.0) then
            if(gx%bmperr.eq.4121) then
               if(.not.silent) write(kou,*) &
                    'Phase ',name1(1:len_trim(name1)),&
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
!         write(*,607)name1,iph
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
         gx%bmperr=7777; goto 1000
      endif
!      if(dodis.eq.1) write(*,*)'Reading disordered parameters'
      ip=nextc
      funname=longline(ip:)
      kp=index(funname,' ')
! save position after parameter name in nextc
      nextc=ip+kp
      funname(kp:)=' '
! extract symbol, normally G or L but TC, BMAGN and others can occur
      lp1=index(funname,'(')
      name1=funname(1:lp1-1)
      typty=0
! this is kept for compatibility with TDB files generated by TC
      if(name1(1:2).eq.'G ' .or. name1(1:2).eq.'L ') then
         typty=1
      elseif(name1(1:3).eq.'TC ') then
         typty=2
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
! not disordered phase, skip this parameter
         goto 100
!-----------------------
710      continue
!         write(*,*)'Entering disordered parameter to: ',thisdis,jl
         thisdis=jl
!         write(*,*)'Entering disordered parameter to: ',ordpartph(thisdis)
!         write(*,*)'3E ',longline(1:len_trim(longline))
         name2=ordpartph(jl)
         fractyp=2
      endif dispar
      call find_phase_by_name_exact(name2,jph,kkk)
!      write(*,*)'readtdb 19: ',jph,gx%bmperr,name2
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
            gx%bmperr=0; goto 100
!         write(*,*)'readtdb error: ',gx%bmperr,name4
!         goto 1000
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
      longline=longline(nextc:)
!410    continue
      jp=len_trim(longline)
      if(longline(jp:jp).ne.'!') then
         if(.not.silent) write(kou,410)nl,ip,longline(1:ip)
410      format('Error, parameter line not ending with !',2i5/a)
         gx%bmperr=7777; goto 1000
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
      call enter_tpfun(funname,longline,lrot,.TRUE.)
!          write(*,17)lokph,typty,nsl,lrot,(endm(i),i=1,nsl)
17 format('readtdb 17: '4i3,5x,10i3)
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
!            write(*,*)'readtdb entparerr: ',gx%bmperr,' >',&
!                 funname(1:len_trim(funname))
            if(gx%bmperr.eq.7778 .and. .not.silent) &
                 write(*,*)'3E Error 7778 at line: ',nl
            gx%bmperr=0
!         elseif(dodis.eq.1) then
!            write(*,*)'Disordered parameter should be entered ok'
         endif
      endif
      if(gx%bmperr.ne.0 .and. .not.silent) write(*,*)'3E error 1: ',gx%bmperr
!------------------------------------------------------------------
!   elseif(line(2:17).eq.'TYPE_DEFINITION ') then
   case(7) !TYPE_DEFINITION 
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
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
               if(.not.silent) write(kou,82) &
                    disparttc,ordpartph(disparttc),dispartph(disparttc)
!                    longline(1:len_trim(longline))
!82             format('Found a type_def DIS_PART:',a,' : ',a)
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
            endif
         endif magnetic
      endif newtypedef
!---------------------------------------------------------------------
!   elseif(line(2:20).eq.'LIST_OF_REFERENCES ' .or. &
!          line(2:16).eq.'ADD_REFERENCES ') then
   case(8,9) ! LIST_OF_REFERENCES and ADD_REFERENCES 
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
         gx%bmperr=7777; goto 1000
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
            gx%bmperr=7777; goto 1000
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
         gx%bmperr=7777; goto 1000
      endif citationmarks
777   continue
!      write(*,*)'Read ',nrefs,' references, ending at',nl
!----------------------------------------------------------------
   case(10) ! ASSESSED_SYSTEMS
      if(.not.silent) write(kou,*) &
           'Cannot handle ASSESSED_SYSTEMS ending at ',nl
      warning=.TRUE.
! skip lines until !
      do while(index(line,'!').le.0)
         read(21,110)line
         nl=nl+1
         call replacetab(line,nl)
      enddo
!------------------------------------------------------------------
   case(11) ! DATABASE_INFORMATION
      if(.not.silent) write(kou,*)'Cannot handle DATABASE_INFORMATION at ',nl
      warning=.TRUE.
! skip lines until !
      do while(index(line,'!').le.0)
         read(21,110)line
         nl=nl+1
         call replacetab(line,nl)
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
   end select
   if(gx%bmperr.ne.0 .and. .not.silent) write(kou,*) &
        '3E errorcode 2: ',gx%bmperr
! look for next KEYWORD
   goto 100
!--------------------------------------------------------
!----- reading functions at the end
800 continue
!   barafun: if(onlyfun) then
! enter only functions that are undefined
!      if(line(2:10).eq.'FUNCTION ') then
!      write(*,*)'Input line >',line(1:20),'<'
!      ipp=istdbkeyword(line,nextc)
!      if(ipp.eq.5) then
!123456789.123456789.123456789.123456789.123456789.123456789.123456789.12345678
! FUNCTION GHSERCR    2.98150E+02  -8856.94+157.48*T-26.908*T*LN(T)
!         name1=line(11:18)
! special case, error in TDB file, UN_ASS is only 6 characters
!         if(name1(1:6).eq.'UN_ASS') then
!            name1=line(11:16); ipp=18
!         else
!            ipp=20
!         endif
   if(eolch(line,nextc)) then
      if(.not.silent) write(kou,*) &
           'Function name must be on same line as FUNCTION'
      gx%bmperr=4000; goto 1000
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
            call enter_tpfun(name1,longline,lrot,.TRUE.)
            if(gx%bmperr.ne.0) then
! one may have error here
               if(.not.silent) write(kou,*)'Failed entering function: ',name1
               goto 1000
            endif
            if(ocv()) write(*,*)'Entered function: ',name1
            nofunent=nofunent+1
         endif
      else
! reset error code
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
   if(warning) then
1001  continue
! if silent set ignore warnings
      if(.not.silent) then
         write(kou,1003)
1003     format(/'There were warnings, continue? Y/N')
         read(kiu,1004)ch1
1004     format(a)
         if(ch1.eq.'N') stop 'warnings reading database'
         if(ch1.ne.'Y') then
            write(kou,*)'Please answer Y or N'
            goto 1001
         endif
      endif
   endif
!   write(*,*)'3E At label 1000'
   if(buperr.ne.0 .or. gx%bmperr.ne.0) then
      if(gx%bmperr.eq.0) gx%bmperr=buperr
      if(.not.silent) write(kou,1002)gx%bmperr,nl
1002   format('Error ',i5', occured at TDB file line ',i7)
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
   rewind: if(dodis.eq.0 .and. disparttc.gt.0) then
! rewind to read disordred parts
      if(.not.silent) &
           write(kou,*)'Rewind to read disordered parts of phases: ',disparttc
      rewind(21)
      dodis=1
      nl=0
      goto 100
   elseif(.not.onlyfun) then
! rewind to read referenced functions
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
   endif rewind
   goto 1000
! end of file while looking for ! terminating a keyword
2200 continue
   if(.not.silent) write(kou,2210)nl,longline(1:72)
2210 format('End of file at ',i5,' looking for end of keyword:'/a)
   gx%bmperr=7777
   goto 1000
 end subroutine readtdb

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine checktdb(filename,nel,selel)
! checking a TDB file exists and return the elements
   implicit none
   integer nel
   character filename*(*),selel(*)*2
!\end{verbatim}
   character line*256
   integer ipp,nl,kk
!
   if(.not.(index(filename,'.tdb').gt.0 &
       .or. index(filename,'.TDB').gt.0)) then
! no extention provided
      filename(len_trim(filename)+1:)='.TDB'
   endif
   open(21,file=filename,access='sequential',form='formatted',&
        err=1010,iostat=gx%bmperr,status='old')
! just check for ELEMENT keywords
! return here to look for a new keyword, end-of-file OK here
   nl=0
   nel=0
100 continue
   read(21,110,end=2000)line
110 format(a)
   nl=nl+1
! One should remove TAB characters !! ??
   ipp=1
   if(eolch(line,ipp)) goto 100
   if(line(ipp:ipp).eq.'$') goto 100
! look for ELEMENT keyword, ipp=1
   ipp=istdbkeyword(line,kk)
   if(ipp.ne.1) goto 100
!
! ignore /- and VA
   if(line(kk+1:kk+2).eq.'/-' .or. line(kk+1:kk+2).eq.'VA') goto 100
   nel=nel+1
   selel(nel)=line(kk+1:kk+2)
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
 end subroutine checktdb

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!-\begin{verbatim}
 subroutine gtpsavetm(filename,str)
! save all data on file in a modified TDB format.  Also as macro and LaTeX
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
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!
 subroutine readtdbsilent
   globaldata%status=ibset(globaldata%status,GSSILENT)
   return
 end subroutine readtdbsilent
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!

