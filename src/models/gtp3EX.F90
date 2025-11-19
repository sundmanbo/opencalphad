!
! gtp3EX included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     9B. Section: read and save on files using XML bases XTDB
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine read_xtdb_dummy
!\begin{verbatim}
 subroutine read_xtdb_dummy(filename,nel,selel)
! reading data from an XTFB file and list data on screen
! if nel<=0 then only extract elements and retutn in selel (-nel dimension)
! if nel>0 extract data for the the nel elements in selel
!
   implicit none
   integer nel
   character filename*(*),selel(*)*2
!\end{verbatim}
! for extracting elements, phases and data from database
   integer warnings,ip,jp,kp,tagno,lastp,tag,fl,taglevel
!   character (len=2) :: present(78)
   character line*512,lline*2048
! when reading a tag with nested tag this is the end of the tad
! up to "nestedtags" levels of nested tags allowed
   integer, parameter :: nestedtags=10
   character(len=24), dimension(nestedtags) :: tagend
   character tagname*18
   type(gtp_xtdbcompatibility) :: pmpid
   logical onlyfun,comment

!   write(*,*)'3EX read_xtdb not implemented yet'
!   goto 1000
   if(len(xtdbtags(1)).gt.18) then
      write(*,*)'3EX xtdbtags longer than 18 characters, extend tagend!',&
           len(xtdbtags)
   endif
! initiating xtdbmpid
   call xtdbinitmpid(nxtdbmpids)
!
!
   if(index(filename,'.').eq.0) then
! no extention provided
      filename(len_trim(filename)+1:)='.XTDB'
   endif
!
   write(*,*)'Reading XTDB files not implemented'
   gx%bmperr=9999; goto 1000
!
   open(21,file=filename,access='sequential',form='formatted',&
        err=2010,iostat=gx%bmperr,status='old')
   onlyfun=.FALSE.
!
   warnings=0
! This is the current level of nexted tags 
   taglevel=0
   comment=.FALSE.
! if the tag and its attributes continue on next line lastp end of previous line
   lastp=0
! this is current number of lines read from file
   fl=1
!-------------------------
! Reading the XTDB file:
! 1. If no elements  provided just extract all the elements
! 2. If elements provides read first the Models tag and possibly modify the MPID
!    2.1 if an unknown model give a warning but no error unless used by a phase
! 3. If elements then read all species and phases that can form
!    3.1 There can be phases rejected or selected
!    3.2 If a phase have an unknown model skip it with a warning
! 4. Read the parameters for the phases and species entered
! 5. Maybe perform some conditional action (extra composition sets,
!       default constituents)
!-------------------------
   readfile: do while(.true.)
! return here to read next line from file
      read(21,110,end=900,err=2010)line
110   format(a)
      fl=fl+1
      ip=1
      if(taglevel.eq.0) then
! No current tag
! <tagname attributes /> or 
! <tag attribures until > (on several lines) with </tagname> on a later line
!    max nestedtags (=10) tags possible
! we expect to read the beginning of a tag or comment
!
!         call gettag(line,ip,tagname)
!
! extract tag and all its attributes (can be on several lines)
         jp=ip
!         call xtdbkey_old(line,jp,tagname)
         select case(tag)
! case default means keyw not understood
         case default
            write(*,*)'3EX no such XTDB tag: ',line(ip:jp)
! handle all tags here
         case(1)

         end select
      else
! we are reading a continuation line of a tag, can be a new tag
         continue
      end if
!
!
!
   end do readfile
! we have finished reading from the fil, we may need to read it agaon
900 continue
   close(21)
!
   if(warnings.ne.0) then
      write(*,*)'3EX warning ',warnings
   endif
!
1000 continue
   return
! error opening file
2000 continue
   write(*,'("3EX error opening file ",i7)')gx%bmperr
   goto 1000
2010 continue
   write(*,'("3EX error reading file ",i7,", line ",i7)')gx%bmperr,fl
   goto 1000
 end subroutine read_xtdb_dummy

!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

  subroutine xtdbinitmpid(ns)
    integer ns
! initiating xtdbmpid
! use the xtdbmodel records
    allocate(xtdbmodel(ns))
!
! maybe one should reconsider the indexing in the result data arrays ...
! In the new models all will have a GEIN model, some a magnetic model
! Note the liquid may be magnetic ...
! these are 9 models allowed in the AmendPhase tag
! for the ocix values: 1=G, 2=LNTH, 3=BMAGN, 4=TC, 5=NT, 6=G2
! Magnetic model 1, previously AFF=-1
!    write(*,*)'Initating xtdbmodel'
    xtdbmodel(1)%modelid='IHJBCC'
    xtdbmodel(1)%nmpid=2
    allocate(xtdbmodel(1)%mpid(2))
    allocate(xtdbmodel(1)%ocmpid(2))
    allocate(xtdbmodel(1)%ocix(2))
    xtdbmodel(1)%mpid(1)='TC'
    xtdbmodel(1)%mpid(2)='BMAG'
    xtdbmodel(1)%ocmpid(1)='TC'
    xtdbmodel(1)%ocmpid(2)='BMAG'
! The ocix is the location for the parameter in the result array in GTP ????
    xtdbmodel(1)%ocix(1)=4
    xtdbmodel(1)%ocix(2)=3
! Magnetic model 2, previously AFF=-3
    xtdbmodel(2)%modelid='IHJREST'
    xtdbmodel(2)%nmpid=2
    allocate(xtdbmodel(2)%mpid(2))
    allocate(xtdbmodel(2)%ocmpid(2))
    allocate(xtdbmodel(2)%ocix(2))
    xtdbmodel(2)%mpid(1)='TC'
    xtdbmodel(2)%mpid(2)='BMAG'
    xtdbmodel(2)%ocmpid(1)='TC'
    xtdbmodel(2)%ocmpid(2)='BMAG'
    xtdbmodel(2)%ocix(1)=4
    xtdbmodel(2)%ocix(2)=3
! Magnetic model 3, previously AFF=0
    xtdbmodel(3)%modelid='IHJQX'
    xtdbmodel(3)%nmpid=3
    allocate(xtdbmodel(3)%mpid(3))
    allocate(xtdbmodel(3)%ocmpid(3))
    allocate(xtdbmodel(3)%ocix(3))
    xtdbmodel(3)%mpid(1)='TC'
    xtdbmodel(3)%mpid(2)='NT'
    xtdbmodel(3)%mpid(3)='BMAG'
    xtdbmodel(3)%ocmpid(1)='TC'
    xtdbmodel(3)%ocmpid(2)='NT'
    xtdbmodel(3)%ocmpid(3)='BMAG'
    xtdbmodel(3)%ocix(1)=4
    xtdbmodel(3)%ocix(2)=5
    xtdbmodel(3)%ocix(3)=3
! Einstein
    xtdbmodel(4)%modelid='GEIN'
    xtdbmodel(4)%nmpid=1
    allocate(xtdbmodel(4)%mpid(1))
    allocate(xtdbmodel(4)%ocmpid(1))
    allocate(xtdbmodel(4)%ocix(1))
    xtdbmodel(4)%mpid(1)='LNTH'
    xtdbmodel(4)%ocmpid(1)='LNTH'
    xtdbmodel(4)%ocix(1)=2
! Liq2State
    xtdbmodel(5)%modelid='LIQ2STATE'
    xtdbmodel(5)%nmpid=2
    allocate(xtdbmodel(5)%mpid(2))
    allocate(xtdbmodel(5)%ocmpid(2))
    allocate(xtdbmodel(5)%ocix(2))
    xtdbmodel(5)%mpid(1)='LNTH'
    xtdbmodel(5)%mpid(2)='G2'
    xtdbmodel(5)%ocmpid(1)='LNTH'
    xtdbmodel(5)%ocmpid(2)='G2'
    xtdbmodel(5)%ocix(1)=2
    xtdbmodel(5)%ocix(2)=6
!
! These models have no model parameters
    xtdbmodel(6)%modelid='FCC4PERM'
    xtdbmodel(6)%nmpid=0
    xtdbmodel(7)%modelid='BCC4PERM'
    xtdbmodel(7)%nmpid=0
    xtdbmodel(8)%modelid='EEC'
    xtdbmodel(8)%nmpid=0
    xtdbmodel(9)%modelid='EBEF'
    xtdbmodel(9)%nmpid=0
!
! The DisorderedPart is a separate tag
!
    return
  end subroutine xtdbinitmpid
 
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine write_xtdbformat
!\begin{verbatim}
 subroutine write_xtdbformat(filename,ext)
! write an XTDB database
!--------------------------------------------------
! NOTE writing TDB files is in gtp3C.F90 by subroutine list_phase_data2
! The XTDB format is defined in gtp3_xml.F90
!--------------------------------------------------
   implicit none
   character*(*) filename,ext
!\end{verbatim}
   integer ip,jp,kp,lut,nl,ia,ib,ic,id,nk,isp
   integer lokph,ics,lokcs,proptyp,lokxx
!   integer, parameter :: maxint=20
   integer :: freemodels=0
   character date*20,tlim8*8,ch1*1,configmodel*32
!   character line*2000,text16*16,text80*80,text256*256,text512*512
! there was a TPfun in TAFID with more than 700 characters .... BAGAS?
   character line*2000,text16*16,text80*80,text256*256,text512*1024
! eventually list a short description of models used in this database
   character, dimension(50) :: usedmodels*24
   logical lrange
   TYPE(gtp_phase_varres), pointer :: phvarres
   TYPE(gtp_phase_add), pointer :: addrec
!
10 format(a)
!
!   write(*,*)'3E: in write_xtbformat: XTDB format output not yet finished
!   write(*,*)'3E: filename: "',trim(filename),'" "',ext,'"'
! make sure extention is .XTDB
   if(index(filename,ext).le.0) then
      ip=index(filename,' ')
      filename(ip:)='.'//ext
   endif
   write(*,*)'Output on: ',trim(filename)
! open the file
   lut=27
   open(lut,file=filename,access='sequential',form='formatted',&
        err=1100,iostat=gx%bmperr,status='unknown')
   nl=0
! write heading
! xtdbversion is in gtp3_xml.F90?
   call date_and_time(date)
   write(lut,80)trim(xtdbversion),version(2:7),date(1:4),date(5:6),date(7:8)
! There should probably be some intial XML stuff
80 format('<XTDB Version="',a,'" Software="OpenCalphad ',a,'" ',&
         'Date="',a,'-',a,'-',a,'" >')
   nl=nl+1
! values of lowtdef,hightdef,bibrefdef and eldef are set in gtp3_xml.F90
! but can be changed by user or when reading an XTDB database
   if(eldef(1:1).eq.' ') then
      write(lut,90)trim(lowtdef),trim(hightdef),trim(bibrefdef)
90    format('  <Defaults LowT="',a,'" HighT="',a,'" Bibref="',a,'" />')
   else
      write(lut,91)trim(lowtdef),trim(hightdef),trim(bibrefdef),trim(eldef)
91    format('  <Defaults LowT="',a,'" HighT="',a,'" Bibref="',a,&
           '" Elements="',a,'" />')
   endif
   nl=nl+1
   write(lut,95)ModelAppendXTDB
95 format('  <DatabaseInfo Info="Test version of XTDB format gtp3EX.FOR" />'/&
          '  <AppendXTDB Models="',a,'" />')
   nl=nl+2
! Writing in this order (option writing parameters by phase?)
! 1: Elements and species
! 2: Phases
! 3: Parameters
! 4: TPfuns
! 5: Bibliography
!-----------------------------
! 1: Elements and species
   do ia=1,noofel
! skip /- and VA
      write(lut,100)trim(ellista(ia)%symbol),trim(ellista(ia)%ref_state),&
           ellista(ia)%mass,ellista(ia)%h298_h0,ellista(ia)%s298
100   format(2x,'<Element Id="',a,'" Refstate="',a,'" Mass="',1PE12.6,&
           '" H298="',1PE12.6,'" S298="',1PE12.6'" />')
      nl=nl+1
   enddo
! list alse elements as species
   do ia=1,noofsp
      text80=' '
! ip is set to 1 inside encode_stoik, text is the stoichiometry with 8 digits
      call encode_stoik(text80,ip,8,ia)
! if MQMQA or UNIQUAC ???????
      write(lut,110)trim(splista(ia)%symbol),text80(1:ip)
110   format(2x,'<Species Id="',a,'" Stoichiometry="',a,'" />')
      nl=nl+1
   enddo
   write(*,*)'No check for MQMQA or UNIQUAC species'
!----------------------------------------------------------------------
! in gtp3C.F90 there is a subroutine list_phase_data2 for TDB files
! 2: Phases
   phaseloop: do ia=1,noofph
      lokph=phases(ia)
      ics=1
!      write(lut,*)'<!-- phase start -->'
!      write(*,*)'3E phase ',trim(phlista(lokph)%name),ia,lokph
! Default is CEF and solid phase, test some status bits
!-Bits in PHASE record STATUS1 there are also bits in each phase_varres record!
!  0 HID phase is hidden (not implemented)
!  1 IMHID phase is implictly hidden (not implemented)
!  2 ID phase is ideal, substitutional and no interaction
!  3 NOCV phase has no concentration variation
!  4 HASP phase has at least one parameter entered
!  5 FORD phase has 4 sublattice FCC ordering with parameter permutations
!  6 BORD phase has 4 sublattice BCC ordering with parameter permutations
!  7 SORD phase has TCP type ordering (not subract ordered as disordered, NEVER)
!  8 MFS phase has a disordered fraction set (DisorderedPart)
!  9 GAS this is the gas phase (first in phase list) 
! 10 LIQ phase is liquid (can be several but listed directly after gas)
! 11 IONLIQ phase has ionic liquid model (I2SL)
! 12 AQ1 phase has aqueous model (not implemented)
! 13 2STATE elemental liquid twostate model parameters (not same as I2SL!)
! 14 QCE phase has corrected quasichemical entropy (Hillerst-Selleby-Sundman)
! 15 CVMCE phase has some CVM ordering entropy (not implemented, SEE CVMTFL)
! 16 EXCB phase need explicit charge balance (has ions)
! 17 XGRID use extra dense grid in gridminimizer for this phase (not used ?)
! 18 MQMQA phase has FACT quasichem SRO model - implementation pending
! 19 NOCS not allowed to create composition sets for this phase
! 20 HELM parameters are for a Helmholz energy model (not implemented),
! 21 PHNODGDY2 phase model with no analytical 2nd derivatives (not implemented)
! 22 not used
! 23 EECLIQ this is the condensed phase (liquid) for highest entropy
! 24 PHSUBO special use testing models DO NOT USE
! 25 PALM interaction records numbered by PALMTREE NEEDED FOR PERMUTATIONS !!!
! 26 MULTI may be used with care
! 27 BMAV Xiong magnetic model with average Bohr magneton number
! 28 UNIQUAC The UNIQUAC fluid model
! 29 TISR phase has the TSIR entropy model (E Kremer)
! 30 PHSSRO phase has the tetrahedral FCC model for SRO (without LRO)
! 31 SROT phase has a tetrahedron quasichemical model -- NOT USED
! 32 CVMTFL phase has the tetrahedral FCC for LRO and SRO (not impl)
! some bits tested later for AmendPhase and DisorderedPart
      configmodel='CEF'; ch1='S'
      if(btest(phlista(lokph)%status1,PHGAS)) then
         configmodel='IDEAL'
         ch1='G'
      elseif(btest(phlista(lokph)%status1,PHIONLIQ)) then
         configmodel='I2SL'
         ch1='L'
      elseif(btest(phlista(lokph)%status1,PHMQMQA)) then
         configmodel='MQMQA'
         ch1='L'
      elseif(btest(phlista(lokph)%status1,PHUNIQUAC)) then
         configmodel='UNIQUAC'
         ch1='L'
      elseif(btest(phlista(lokph)%status1,PHLIQ)) then
         ch1='L'
      endif
      line='  <Phase Id="'//trim(phlista(lokph)%name)//&
           '" Configuration="'//trim(configmodel)//'" State="'//ch1//'" >'
      write(lut,10)trim(line)
      nl=nl+1
! sublattices/sublattices
      line='    <Sublattices NumberOf="'
      ip=len_trim(line)+1
      ib=phlista(lokph)%noofsubl
! wriint update the position
      call wriint(line,ip,ib)
      line(ip:)='" Multiplicities="'
      ip=len_trim(line)+1
! the sites are in phase_varres record ... index in %linktocs
      lokcs=phlista(lokph)%linktocs(1)
      phvarres=>firsteq%phase_varres(lokcs)
      do ic=1,ib
! ip postion (updates) 10 is max digits, 0 means no + sign (>0 means sign)
         call wrinum(line,ip,10,0,phvarres%sites(ic))
         ip=ip+1
      enddo
      line(ip-1:)='" >'
      write(lut,10)trim(line)
      nl=nl+1
! how are constiuents stored
!      write(*,70)phlista(lokph)%constitlist
70    format('3E constituents: ',20i3)
! constituents
! loop for all sublattices
      nk=0
      do ic=1,ib
         line='      <Constituents Sublattice="'
         ip=len_trim(line)+1
         call wriint(line,ip,ic)
         line(ip:)='" List="'
         ip=ip+8
! loop for all constituents
         do id=1,phlista(lokph)%nooffr(ic)
            nk=nk+1
            isp=phlista(lokph)%constitlist(nk)
            line(ip:)=splista(isp)%symbol
            ip=ip+len_trim(splista(isp)%symbol)+1
         enddo
         line(ip-1:)='" />'
!         write(*,10)'3E check; ',trim(line)
         write(lut,10)trim(line)
         nl=nl+1
      enddo
! end of sublattices
      write(lut,250)
250   format('    </Sublattices>')
      nl=nl+1
! DisorderedPart, in OC the ordered and disordered parameters in same phase
      if(btest(phlista(lokph)%status1,PHMFS)) then
         write(*,*)'3EX The phase '//trim(phlista(lokph)%name)//&
              ' has a DisorderedPart'
! in gtp_fraction_set the data for disordered phase can be found
!         write(*,*)'3E Disordere:',phvarres%disfra%latd,phvarres%disfra%totdis
         line='    <DisorderedPart Sum="'
         ip=len_trim(line)+1
         call wriint(line,ip,phvarres%disfra%latd)
         if(phvarres%disfra%totdis.eq.1) then
            line(ip:)='" Subtract="Y" />'
         else
            line(ip:)='" />'
         endif
!         if(btest(phlista(lokph)%status1,PHSORD)) then
! I am not sure how this is connected with Disordered_Part in OC
!            write(*,*)'3E do not subtract ordered as disordered'
!         endif
!         write(*,10)'3E check: ',trim(line)
         write(lut,10)trim(line)
         nl=nl+1
      endif
!------------------------------------
! AmendPhase: magnetism etc as additions
      text512=' '
      ip=1
      addrec=>phlista(lokph)%additions
      do while(associated(addrec))
         proptyp=addrec%type
!         write(*,*)'3E addrec%propval: ',proptyp
         select case(proptyp)
         case default
            write(*,*)'3E Unknown property: ',proptyp
         case(1) !  INDENMAGNETIC=1, BCC and other phases
            if(addrec%aff.eq.-1) then
               text512(ip:)='IHJBCC'; ip=ip+7
            elseif(addrec%aff.eq.-3) then
               text512(ip:)='IHJREST'; ip=ip+8
            endif
         case(2) !  XIONGMAGNETIC=2 same for all
            text512(ip:)='IHJQX'; ip=ip+6
         case(3)!   DEBYECP=3, not implemented
            continue
         case(4) !  EINSTEINCP=4 
            text512(ip:)='GEIN'; ip=ip+9
         case(5) !  TWOSTATEMODEL1=5
            text512(ip:)='LIQ2STATE'; ip=ip+10
         case(6) !  ELASTICMODEL1=6
         case(7) !  VOLMOD1=7
! OC by default set VOLMOL
            continue
         case(8) !   UNUSED_CRYSTBREAKDOWNMOD=8
         case(9) !   SECONDEINSTEIN=9
         case(10) !  SCHOTTKYANOMALY=10
         case(11) !  DIFFCOEFS=11
         end select
         addrec=>addrec%nextadd
      enddo
! these amendments/amendments are set by other bits
      if(btest(phlista(lokph)%status1,PHFORD)) then
         text512(ip:)='FCC4PERM'; ip=ip+9
      elseif(btest(phlista(lokph)%status1,PHBORD)) then
         text512(ip:)='BCC4PERM'; ip=ip+9
      endif
      if(ip.gt.1) then
!         write(*,*)'3E additions: ',trim(text512),ip
         if(len_trim(text512).gt.1) write(lut,270)trim(text512)
270      format('    <AmendPhase Models="',a,'" />')
         nl=nl+1
      endif
!------------------------------------------ suck ....
! CrystalStructure ?
      write(lut,290)
290   format('  </Phase>')
      nl=nl+1
!
!      write(*,*)'3E Finished this phase'
   enddo phaseloop
!   goto 900
   write(lut,*)' <!-- all phases written -->'
!-----------------------------
! 3: Parameters, phase by phase
! code here copied from gtp3C.F90 subroutine list_phase_data
   do ia=1,noph()
!      write(*,*)'3E calling write_parametrar',lut
      lokph=phases(ia)
      call write_parameters(lokph,lut,2,nl)
      if(gx%bmperr.ne.0) goto 1100
!      write(*,*)'3E Finished listing parameters for this phase'
   enddo
!----------------------------- CHANGE HERE if gtp3_xml.F90 changes
! 4: TPfuns  this is the tag xmlel(26) 
   write(lut,*)' <TPfun Id="R"     Expr="8.31451;" />'
   write(lut,*)' <TPfun Id="RTLNP" Expr="R*T*LN(1.0E-5)*P);" />'
   nl=nl+2
!   write(*,*)'3E tpfuns: ',notpf(),freetpfun
   tpfuns: do ia=3,notpf()
      text512=' '
      tlim8=' '
      lrange=.FALSE.
! A TPfun can be very long
      call list_tpfun(ia,0,text512)
      if(text512(1:1).eq.'_') cycle tpfuns
! check that there is a " N " in text512 to indicate end of expression
      if(index(text512,' N ').le.0) then
         ip=index(text512,'= ')-1
         write(*,*)text512(1:ip)
130      format('3E no end of TPfun ',a,' in text512')
         stop
      endif
! we have to format this using TPfun anda Trange tags
! we should use lowtdef and hightdef which are 8 characters
      ip=index(text512,'= ')+2
      tlim8=text512(ip:)
      jp=index(tlim8,' ')
      tlim8(jp:)=' '
      if(tlim8.ne.lowtdef) then
         line='  <TPfun Id="'//text512(1:ip-4)//'" LowT="'//trim(tlim8)//&
              '" Expr="'
      else
! Default LowT
         line='  <TPfun Id="'//text512(1:ip-4)//'" Expr="'
      endif
      ip=ip+index(text512(ip:),' ')
      kp=len_trim(line)+2
! ip is after lowT in text512 and kp is after Expr="
! there can be breakpoints in T
      tranges: do while(.TRUE.)
         jp=index(text512(ip:),';')
         if(jp.le.0) then
            write(*,*)'Missing ; at end of expression'
            stop
         endif
         line(kp:)=text512(ip:ip+jp)
         kp=kp+jp
         ip=ip+jp+1
! HighT limit or more ranges?
         jp=index(text512(ip:),' Y ')
         if(jp.gt.0) then
! more ranges, save end of this range do not check highT limit
            tlim8=text512(ip:ip+jp-2)
            if(tlim8.ne.hightdef) then
               if(lrange) then
                  line(kp:)='" HighT="'//text512(ip:ip+jp-2)//'" />'
               else
                  line(kp:)='" HighT="'//text512(ip:ip+jp-2)//'" >'
               endif
            else
               if(lrange) then
                  line(kp:)='" />'
               else
                  line(kp:)='" >'
               endif
            endif
            lrange=.TRUE.
            write(lut,10)trim(line)
            nl=nl+1
            ip=ip+jp+2
! This is the tag xmlel(26)
            line='    <Trange Expr="'
            kp=19
         else
! no more ranges, end of single TPfun or current Trange
            jp=index(text512(ip:),' ')
            tlim8=text512(ip:ip+jp-2)
            if(tlim8.ne.hightdef) then
               line(kp:)='" HighT="'//text512(ip:ip+jp-2)//'" />'
            else
! at the highT limit
               line(kp:)='" />'
            endif
            ip=ip+jp+1
            kp=len_trim(line)+1
            write(lut,10)line(1:kp)
            nl=nl+1
! if there has been Trange tags then end the TPfun, tag xmlel(26)
            if(lrange) write(lut,10)'  </TPfun>'
            cycle tpfuns
         endif
! there are more tranges
!         write(*,*)'3E next range: ',trim(text512(ip:))
      enddo tranges
   enddo tpfuns
!-----------------------------
! 5: Bibliography
!      write(*,*)'3EX reffree: ',reffree
      if(reffree.gt.1) then
         write(lut,400)
400      format('  <Bibliography>')
         nl=nl+1
         do ia=1,reffree-1
! Wow, reference texts are stored using storec/loadc ... max 1024 chars
            ip=bibrefs(ia)%wprefspec(1)
            line=' '
            call loadc(2,bibrefs(ia)%wprefspec,line(1:ip))
! check to < > and &
!            call check_illegal_xml(line,ip)
!            write(*,*)'Bibitem: ',trim(line),ip
            write(lut,410)trim(bibrefs(ia)%reference),trim(line)
410         format('    <Bibitem Id="',a,'" Text="',a,'" /> ')
            nl=nl+1
         enddo
         write(lut,420)trim(bibrefdef)
420      format(4x,'<Bibitem Id="Default" Text="',a,'" /> '/'  </Bibliography>')
         nl=nl+1
      else
      endif
!-----------------------------
! Models in AppendXTDB
! Finished !!
900 continue
   write(lut,990)
990 format('</XTDB>')
   nl=nl+1
!----------------------------
1000 continue
   close(lut)
   write(*,1010)nl,trim(filename)
!
1010 format('Written: ',i7,' lines on ',a)
   return
! error open or writing
1100 continue
   write(*,*)'Error during writing XTDB file',gx%bmperr
   goto 1000
!   
 end subroutine write_xtdbformat
 
!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

! \addtotable subroutine check_illegal_xml
 subroutine replace_illegal_xml(line,ip)
! replace < > & and " in references by
!         [ ] | and #
   character*(*) line
   integer ip
!
   integer jp
   lt: do while(.true.)
      jp=index(line(1:ip),'<')
      if(jp.eq.0) exit lt
      line(jp:jp)='['
   enddo lt
   gt: do while(.true.)
      jp=index(line(1:ip),'>')
      if(jp.eq.0) exit gt
      line(jp:jp)=']'
   enddo gt
   amp: do while(.true.)
      jp=index(line(1:ip),'&')
      if(jp.eq.0) exit amp
      line(jp:jp)='|'
   enddo amp
   quote: do while(.true.)
! this means " is replaced by a single '
      jp=index(line(1:ip),'''')
      if(jp.eq.0) exit quote
      line(jp:jp)='#'
   enddo quote
   return
 end subroutine replace_illegal_xml

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

! \addtotable subroutine exp2xml
 subroutine exp2xml(lut,expr,tag,nl)
! convert expr to TPfun/Parameter and Trange and write lines on unit lut
! tag is TPfun or Parameter and the first part of expr is the Id followd
! by the expression
   integer lut,nl
   character expr*(*),tag*(*)
!
   integer ip,jp,kp
   character tlim8*8,line*1000,bibref*16
   logical lrange
!   
! make sure there is a final ' N '
   ip=index(expr,' N ')
   if(ip.le.0) then
      write(*,*)'3E expression has no terminating N'
      stop
   elseif(tag(1:1).eq.'P') then
! extract reference if tag is Parameter
      bibref=trim(expr(ip+3:))
   else
      bibref=' '
   endif
! missing default bibref?? fixed
!   write(*,*)'3E exp2xml: "',trim(expr),'" "',trim(bibref)
! in tlow is after "= "
!   write(*,*)'3E exp2xml: ',len_trim(expr)
!   write(*,*)'3E exp2xml: ',trim(expr)
   lrange=.FALSE.
   ip=index(expr,'= ')+2
   tlim8=expr(ip:)
   jp=index(tlim8,' ')
   tlim8(jp:)=' '
   line='  <'//trim(tag)//' Id="'//expr(1:ip-3)//'"'
   kp=len_trim(line)
!   write(*,*)'3E line 1A: ',line(1:kp)
   if(tlim8.ne.lowtdef) then
      line(kp+2:)=' LowT="'//trim(tlim8)//'" Expr="'
   else
! Default LowT
      line(kp+2:)=' Expr="'
   endif
   kp=len_trim(line)+2
!   write(*,*)'3E line 1B: >',line(1:kp),'<'
! find space after lowT limit ...
!   write(*,*)'3E space1: ',expr(1:ip),ip
!   ip=ip+index(expr(ip:),' ')-1
   ip=ip+index(expr(ip:),' ')
!   write(*,*)'3E after space: ',expr(ip:ip+20)
! ip is after lowT in expr and kp is after Expr="
! there can be breakpoints in T
   tranges: do while(.TRUE.)
      jp=index(expr(ip:),';')
      if(jp.le.0) then
         write(*,*)'Missing ; at end of expression',ip,jp
         write(*,*)'3E expr :',trim(expr(ip:)),':'
         stop
      endif
!      write(*,*)'3E line 1C: >',line(1:kp),'<'
!      write(*,*)'3E exp 2: ',expr(ip:ip+20),ip
!      line(kp:)=expr(ip+1:ip+jp+1)
! problem here because initial sign disappeared!!!
      line(kp:)=expr(ip:ip+jp-1)
!      write(*,*)'3E line 2: >',trim(line),'<'
      kp=kp+jp
      ip=ip+jp+1
! HighT limit or more ranges?
      jp=index(expr(ip:),' Y ')
      if(jp.gt.0) then
! more ranges, save end of this range do not check highT limit
         tlim8=expr(ip:ip+jp-2)
         if(tlim8.ne.hightdef) then
            if(lrange) then
               line(kp:)='" HighT="'//expr(ip:ip+jp-2)//'" />'
            else
! there is a seconed or more ranges, the bibref should be part of TPfun tag
               line(kp:)='" HighT="'//expr(ip:ip+jp-2)//'" >'
               if(bibref(1:1).ne.' ') then
                  kp=len_trim(line)
                  line(kp:)=' Bibref="'//trim(bibref)//'" >'
                  bibref=' '
               endif
            endif
         else
            line(kp:)='" />'
         endif
         lrange=.TRUE.
!         write(*,*)'3E line: ',trim(line)
         write(lut,10)trim(line)
10       format(a)
         nl=nl+1
         nl=nl+1
         ip=ip+jp+2
         line='    <Trange Expr="'
         kp=19
      else
! no more ranges, end of single TPfun or current Trange
!         write(*,*)'3E no more ranges :',trim(expr(ip:)),': '
         jp=index(expr(ip:),' ')
         tlim8=expr(ip:ip+jp-2)
         if(tlim8.ne.hightdef) then
            line(kp:)='" HighT="'//expr(ip:ip+jp-2)//'" />'
         else
! at the highT limit
            line(kp:)='" />'
         endif
! add the Bibref if there has not been any Tranges
         if(bibref(1:1).ne.' ') then
            kp=len_trim(line)-1
            line(kp:)=' Bibref="'//trim(bibref)//'" />'
         endif
!         ip=ip+jp+1
         kp=len_trim(line)
!         write(*,*)'3E line 7: ',trim(line)
         write(lut,10)line(1:kp)
         nl=nl+1
! if there has been Trange tags then end the Parameter/Tpfun tag
!         write(*,*)trim(line)
         if(lrange) then
            write(lut,10)'  </'//trim(tag)//'>'
            nl=nl+1
         endif
         exit tranges
      endif
! loop if there are more tranges
   enddo tranges
1000 continue
   return
 end subroutine exp2xml

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\addtotable subroutine write_parameters
!\begin{verbatim}
 subroutine write_parameters(lokph,lut,typ,nl)
! code below same as in list_all_data ... maybe make it a subroutine ...
! lokph is phlista index, lut is output unit
! typ 1 for screen, 2 for XML, nl counts number of lines written
   implicit none
   integer lokph,lut,typ,nl
!\end{verbatim}
   integer parlist,ll,nsl,ij,ideg,typty,typspec,kk,kkx,iel,linkcon,nz,ics
   integer ncsum,ip,intpq,prplink,nint,lokcs,ik,jdeg,kkk,iqhigh,iqnext,lqq
!
   integer ilist(9), endm(15), lint(2,3)
   character text*2000,prop*32,funexpr*1000,phname*24,ch1*1,toopsp(3)*24
   logical mqmqa,noelin1,subref
! This is the index of PARAMETER in xtdbtags
   integer, parameter :: xmlpartag=17
!
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_endmember), pointer :: endmemrec
   TYPE(gtp_interaction), pointer :: intrec
! a smart way to have an array of pointers
   TYPE intrecarray 
      type(gtp_interaction), pointer :: p1
   end TYPE intrecarray
   integer, parameter :: maxstack=20
   TYPE(intrecarray), dimension(maxstack) :: intrecstack
   TYPE(gtp_fraction_set), pointer :: disfrap
   type(gtp_tooprec), pointer :: tooprec
!
!--------------------------------------------------
! return here to list disordered parameters
!   write(lut,*)'<!-- start write_parameters -->'
!   write(*,*)'3E In write_parameters',lokph,lut,typ
   phname=phlista(lokph)%name
   parlist=1
   ics=1
!   tooptop=0
   mqmqa=.FALSE.
100 continue
! parlist changed below for disordered fraction set
   if(parlist.eq.1) then
      endmemrec=>phlista(lokph)%ordered
      nsl=phlista(lokph)%noofsubl
   else
!      write(*,*)'3E Listing disordred parameters 1',nsl
      endmemrec=>phlista(lokph)%disordered
      disfrap=>firsteq%phase_varres(lokcs)%disfra
      nsl=disfrap%ndd
   endif
!   write(*,*)'3E Listing parameters 2, nsl=',nsl
   endmemberlist: do while(associated(endmemrec))
      do ll=1,nsl
!         ilist(ll)=emlista(lokem)%fraclinks(ll,1)
         ilist(ll)=endmemrec%fraclinks(ll,1)
         if(ilist(ll).gt.0) then
            if(parlist.eq.2) then
! what is disfra here??!!
!               write(*,*)'3E disfra?: ',disfra%splink(ilist(ll)),&
!                    disfrap%splink(ilist(ll))
!               endm(ll)=disfra%splink(ilist(ll))
               endm(ll)=disfrap%splink(ilist(ll))
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
      text=' '
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
!            if(parlist.eq.2) then
! disordered endmember parameter
! DO NOT ADD D on disordered parameter identifiers
!               kk=len_trim(prop)+1
!               prop(kk:kk)='D'
!            endif
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
!                     write(*,*)'3E: endmember typspec 1: ',iel
!                     write(*,*)'3E splink: ',disfrap%splink
                     linkcon=disfrap%splink(iel)
!                     write(*,*)'3E: endmember typspec 2: ',linkcon
                     ll=0
!                     ll=1
! linkcon has nothing to do with which sublattice, ignore ll
!                     if(linkcon.gt.disfrap%nooffr(1)) ll=2
                     prop=prop(1:len_trim(prop))//'&'&
                          //splista(linkcon)%symbol
!                     write(*,*)'3E We are here',linkcon,disfrap%nooffr(1),ll
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
                  write(*,*)'Error in constituent depended parameter id'
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
                  write(*,*)'lpd 7B: ',iel,typty
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
! check if FNN MQMQA parameter ...
         if(mqmqa) then
! ilist is index in fraction list, same as index in mqmqa_data%contyp
!            intpq=ilist(1)
!            write(*,*)'3E check if SNN parameter',intpq,&
!                 mqmqa_data%contyp(5,intpq)
            if(mqmqa_data%contyp(5,ilist(1)).le.0) goto 203
         endif
         if(typ.eq.1) then
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
         endif
203      continue
! this writes the expression as for a TDB file
         call list_tpfun(proprec%degreelink(0),1,funexpr(ip:))
!         write(*,*)' >>>> fun? ',trim(funexpr(ip:)),lut
         ip=len_trim(funexpr)
         if(len_trim(proprec%reference).le.0) then
            funexpr(ip+1:)=' Default '
         else
            funexpr(ip+1:)=' '//proprec%reference
         endif
         ip=len_trim(funexpr)
         if(typ.eq.1) then
! nice output over several lines if needed with indentation 12 spaces
            call wrice2(lun,2,12,78,1,funexpr(1:ip))
         else
!            write(*,10)funexpr(1:ip)
!            write(lut,10)funexpr(1:ip)
! this convert TDB expression to XTDB format, this is the Parameter tag
!            call exp2xml(lut,funexpr,xmlel(xmlpartag),nl)
            call exp2xml(lut,funexpr,xtdbtags(xmlpartag),nl)
         endif
         proprec=>proprec%nextpr
      enddo ptyloop
      if(typ.eq.1) then
         if(btest(phlista(lokph)%status1,PHFORD).or. &
              btest(phlista(lokph)%status1,PHBORD)) then
!      if(endmemrec%noofpermut.gt.1) then
            intpq=0
            if(associated(endmemrec%intpointer)) then
               intpq=endmemrec%intpointer%antalint
            endif
            prplink=0
            if(associated(endmemrec%propointer)) prplink=1
! keep this output for the moment
            if(parlist.eq.1) write(*,207)endmemrec%antalem,&
                 endmemrec%noofpermut,intpq,prplink
207         format('3E Endmember check: permut, interaction, pty: ',4i5)
         endif
      endif
      endmemrec=>endmemrec%nextem
   enddo endmemberlist
!   write(*,*)'3E Finished listing endmember parameters',parlist,nsl
!-----------------------------------------------------------------------
! parameters for interactions using site fractions
   if(parlist.eq.1) then
      endmemrec=>phlista(lokph)%ordered
   else
      if(.not.associated(phlista(lokph)%disordered)) &
         write(*,*)'3E Problems with disordered fraction set'
      endmemrec=>phlista(lokph)%disordered
      if(.not.associated(disfrap)) write(*,*)'3E disfrap problems!'
      nsl=disfrap%ndd
!      write(*,*)'new nsl',nsl
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
               endm(ll)=disfrap%splink(kkx)
            else
               endm(ll)=phlista(lokph)%constitlist(kkx)
            endif
         enddo
      endif
      nint=0
      intlist2: do while(associated(intrec))
         nint=nint+1
         if(nint.gt.maxstack) then
            write(*,*)'3E overflow in intrecstack 1'
            gx%bmperr=4399; goto 1000
         endif
         intrecstack(nint)%p1=>intrec
!--------------------------------------------------
         lint(1,nint)=intrec%sublattice(1)
         kkk=intrec%fraclink(1)
         if(parlist.eq.2) then
            lint(2,nint)=disfrap%splink(kkk)
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
            if(typty.gt.0 .and. typty.le.ndefprop) then
               prop=propid(typty)%symbol
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
!                        write(*,*)'3E: typspec: 3 ',typty,iel,prop(1:10)
                        linkcon=disfrap%splink(iel)
!                        write(*,*)'3E: typspec: 4 ',typty,linkcon,prop(1:10)
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
! note changes here must be repeated for endmember parameters above
            degree: do jdeg=0,proprec%degree
               if(proprec%degreelink(jdeg).eq.0) then
!                  write(*,*)'Ignoring function link'
                  cycle degree
               endif
               call encode_constarr(text,nsl,endm,nint,lint,jdeg)
               write(funexpr,300)trim(prop),trim(phname),trim(text)
300            format(A,'(',A,',',A,') ')
               ip=len_trim(funexpr)+1
               call list_tpfun(proprec%degreelink(jdeg),1,funexpr(ip:))
               ip=len_trim(funexpr)
! default reference missing for excess parameters
               if(len_trim(proprec%reference).le.0) then
                  funexpr(ip+1:)=' Default '
               else
                  funexpr(ip+1:)=' '//proprec%reference
               endif
!               funexpr(ip+1:)=' '//proprec%reference
               ip=len_trim(funexpr)
               if(typ.eq.1) then
                  call wrice2(lun,4,12,78,1,funexpr(1:ip))
               else
!                  write(*,10)funexpr(1:ip)
!                  write(lut,10)funexpr(1:ip)
! This is the Parameter xml tag
!                  call exp2xml(lut,funexpr,xmlel(xmlpartag),nl)
                  call exp2xml(lut,funexpr,xtdbtags(xmlpartag),nl)
10                format(a)
               endif
            enddo degree
            proprec=>proprec%nextpr
         enddo ptyloop2
! list temporarily the number of permutations for FCC and BCC ordering
         pdebug: if(typ.eq.1) then
            if(btest(phlista(lokph)%status1,PHFORD).or. &
                 btest(phlista(lokph)%status1,PHBORD)) then
               if(nint.eq.1) then
                  nz=intrec%noofip(2)
               else
                  nz=size(intrec%sublattice)
                  lqq=intrec%noofip(size(intrec%noofip))
                  if(lqq.ne.nz) then
                     write(*,*)'3E Not same 1: ',intrec%antalint,nz,lqq
                  endif
!               write(*,301)nz,intrec%noofip
301               format('noofip: ',10i3)
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
               prplink=0
               if(associated(intrec%propointer)) prplink=1
! keep this output for the moment
               if(parlist.eq.1) write(*,302)intrec%antalint,&
                    nz,nint,iqhigh,iqnext,prplink
302         format('3E Inter check 1: id, permut, level, high, next, pty: ',&
                 i5,i3,i3,i4,i4,i2)
            endif
         endif pdebug
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
!   write(*,*)'3E Finished listing interactions',parlist,nsl
! check if there are other fraction lists
!   parlist=parlist+1, hm parlist can only be 1 or 2
!   write(*,*)'checking for disordered parameters'
   if(parlist.eq.1 .and. associated(phlista(lokph)%disordered)) then
      subref=.TRUE.
!      lokcs=phlista(lokph)%cslink
      lokcs=phlista(lokph)%linktocs(ics)
! does this make a copy?  Maybe it should be a pointer. IT IS A POINTER!
      disfrap=>firsteq%phase_varres(lokcs)%disfra
      if(.not.associated(disfrap)) then
!         write(*,*)'disfrap OK'
!      else
         write(*,*)'3E disfrap not set, expect segmentation fault?'
      endif
      nsl=disfrap%ndd
!      write(*,810)disfrap%fsites,nsl
      write(lut,810)disfrap%fsites,nsl
810   format('<!-- Disordered fraction set factor: ',F10.4,&
           ' Sublattices: ',i3,' -->')
      parlist=2
!      write(*,*)'3E Jump back to list disordered parameters',nsl,parlist
      goto 100
   endif
! Check if there are toop/kohler ternaries
   tooprec=>phlista(lokph)%tooplast
   if(associated(tooprec)) then
      write(*,*)'3EX there are Toop/Kohler extrapolations'
      write(*,*)'There is some code in gtp3C.F90 to handle this'
   endif
!   write(lut,*)'<!-- exit write_parameters -->'
!   write(*,*)'3E listing by write_parameters
1000 continue
   return
 end subroutine write_parameters

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

!\addtotable subroutine debug_phaseparameters
!\begin{verbatim}
 subroutine debug_phaseparameters(lokph,lut,ceq)
! code to debug a phase structure of parameters
! It follows the data structure as in calcg_internal in gtp3X.F90
! lokph is phlista index, lut is output unit
   implicit none
   integer lokph,lut
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokres
! use composition set 1, phlista is protected inside GTP
   lokres=phlista(lokph)%linktocs(1)
! this call will eventually be replaced by call calcg_internal
! after setting a global debugpar variable.  The subroutines
! debug_endmember and debug_intrec may be integrated in calcg_internal
! to ensure that any new data structures are followed
   write(*,*)'3EX In debug_phaseparameters',lut
  call list_phaseparameters(lokph,lut,ceq%phase_varres(lokres),ceq)
   return
 end subroutine debug_phaseparameters
 
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_phaseparameters
!\begin{verbatim}
 subroutine list_phaseparameters(lokph,lut,cps,ceq)
! code to debug a phase structure of parameters
! It follows the data structure as in calcg_internal in gtp3X.F90
! lokph is phlista index, lut is output unit
   implicit none
   integer lokph,lut
   type(gtp_equilibrium_data), pointer :: ceq
   type(gtp_phase_varres), target :: cps
!\end{verbatim}
! Most declaration copied, some may not be needed
   integer fractype,epermut,maxprec,sameint,nprop,lokres,lokdiseq
   integer moded,nofc2,nsl,msl,qz,incffr(0:maxsubl),spmod,intlevel,ipermut
!
   TYPE(gtp_parcalc) :: gz
   TYPE(gtp_fraction_set), pointer :: fracset,dislink
   TYPE(gtp_phase_varres), pointer :: phres,phpart,phmain
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_endmember), pointer :: endmemrec
   TYPE(gtp_interaction), pointer :: intrec
! array with intrec pointers
   type intstackarray
      TYPE(gtp_interaction), pointer :: save
   end type intstackarray
   type(intstackarray), dimension(5) :: intstack
!   TYPE(gtp_pystack), pointer :: pystack
!   TYPE(gtp_phase_add), pointer :: addrec
! to handle parameters with wildcard constituent and other things
   logical wildc,nevertwice,first,chkperm,ionicliq,iliqsave,iliqva,iliqneut
   integer, parameter :: permstacklimit=150
   integer, dimension(permstacklimit) :: lastpmq,maxpmq
   integer pmq,nz
!
   write(*,*)'3EX in list_phaseparameters',lokph
   if(btest(phlista(lokph)%status1,PHMQMQA)) then
      write(*,*)'3EX phase has MQMQA model, no listing'
      goto 1000
   endif
   spmod=0
   if(btest(phlista(lokph)%status1,PHFORD) .or. &
        btest(phlista(lokph)%status1,PHBORD)) then
! PHPALM is needed for phases with permutations such as ordered FCC/BCC/HCP
      chkperm=.true.
! spmod tries to keep track of disordered/permutation of parameters? >10 permut
      spmod=10
      if(.not.btest(phlista(lokph)%status1,PHPALM)) then
!         write(*,*)'3X calling palmtree ',lokph,cps%phtupx
! This is needed only once unless parameters are changed.  It numbers the
! interaction records sequ+entially for the permutations
! the subroutine palmtree is in gtp3Y.F90 for some unknown reason ...
         call palmtree(lokph)
         if(gx%bmperr.ne.0) goto 1100
         write(lut,300)'Phase has parameter permutations'
300      format(a)
! this must be zeroed if a new interaction parameter is added
!         phlista(lokph)%status1=ibset(phlista(lokph)%status1,PHPALM)
      endif
   endif
!-----------------------------------------------------------------
50  continue
! local work arrays for products of Y and calculated parameters are allocated
   gz%nofc=phlista(lokph)%tnooffr
!
! dimension for number of parameter properties
   nprop=cps%nprop
! phres will point either to ordered or disordered results
! phmain will always point to record for ordered phase_varres
   phmain=>cps
   phres=>cps
!
   nsl=phlista(lokph)%noofsubl
   write(lut,10)trim(phlista(lokph)%name),phlista(lokph)%nooffs
10 format(a,' parameter structure with ',i3,' fraction sets')
   write(*,*)'3EX Inside list_phaseparameters ',phlista(lokph)%nooffs,nsl
!   
   fractype=0
! chkperm true if FCC/HCP or BCC permutation of ordered phases
   chkperm=.false.
   if(btest(phlista(lokph)%status1,PHFORD) .or. &
        btest(phlista(lokph)%status1,PHBORD)) then
! PHPALM is needed for phases with permutations such as ordered FCC/BCC/HCP
      chkperm=.true.
      if(.not.btest(phlista(lokph)%status1,PHPALM)) then
!         write(*,*)'3X calling palmtree ',lokph,cps%phtupx
! This is needed only once unless parameters are changed.
! interaction records sequentially for the permutations
         write(*,*)'3EX parameter permutations initiated'
         call palmtree(lokph)
         if(gx%bmperr.ne.0) goto 1000
      endif
   endif
!
! loop for different types of fractions: site fractions, mole fractions ...
   fractyp: do while(fractype.lt.phlista(lokph)%nooffs)
!
105 continue
      fractype=fractype+1
      write(*,*)'3EX fraction type: ',fractype
! Jump back here for parameters in disordered fraction set (fractype>1)
110   continue
      fracset=>phmain%disfra
! code to handle phases with two fraction sets (Disordered_2Part and _3Part)
      ftype: if(fractype.eq.1) then
         odtest: if(btest(phlista(lokph)%status1,PHMFS)) then
! there is a disordered fraction set            
            if(fracset%totdis.ne.0) then
               if(btest(phlista(lokph)%status1,PHSUBO)) then
                  write(lut,300)'Phase has Disordered_2Part model'
                  write(*,*)'Disordered_2Part model'
                  goto 106
               endif
               write(lut,300)'Phase has Disordered_3Part model'
               if(btest(phmain%status2,CSORDER)) then
! this is a Disordered_3Part model, the ordered part is calculated twice
                  nevertwice=.false.
!               else
! Disordered_3Part model and phase is disordered, skip ordered part (fractype=1)
!                  goto 105
! When listing include the listing of ordered parameters ...
               endif
            endif
         endif odtest
106      continue
! initiate variables for disordered part
         gz%nofc=phlista(lokph)%tnooffr
         msl=nsl
         incffr(0)=0
         do qz=1,nsl
            incffr(qz)=incffr(qz-1)+phlista(lokph)%nooffr(qz)
         enddo
      else
! parameters in disordered fraction set         
         msl=fracset%ndd
!         write(*,*)'3EX disordered fraction set',msl,associated(fracset)
         gz%nofc=fracset%tnoofxfr
         incffr(0)=0
         do qz=1,msl
            incffr(qz)=incffr(qz-1)+fracset%nooffr(qz)
         enddo
! no need to handle fractions and derivatives         
         dislink=>cps%disfra
         lokdiseq=dislink%varreslink
         phres=>ceq%phase_varres(lokdiseq)
! UNFINISHED: moded has no value here ....
         if(moded.gt.1) then
            nofc2=gz%nofc*(gz%nofc+1)/2
!               if(.not.allocated(phres%d2gval)) then
!                  allocate(phres%d2gval(nofc2,nprop))
!               endif
!               phres%d2gval=zero
         endif
      endif ftype
      first=.true.
!
! HERE WE START FOLLOWING THE LINKS BETWEEN ENDMEMBERS
! ordered fraction set listed first
      if(fractype.eq.1) then
         endmemrec=>phlista(lokph)%ordered
      else
         endmemrec=>phlista(lokph)%disordered
      endif
!      write(*,70)'3EX debug endmemberloop',nsl,msl,associated(endmemrec)
!70    format(a,2i3,l3)
!      write(*,*)'3EX we are here 1'
!      write(*,*)'3EX number of permutations: ',endmemrec%noofpermut
!      write(*,*)'3EX we are here 2'
      endmemloop: do while(associated(endmemrec))
! Note all interaction are calculated inside this loop!!!!
! The array maxpmq is used for interaction permutations.  It must be
! initialized to zero at the first endmember permutation.  It is set to
! limits for the interacton permutations for all interaction records.
         maxpmq=0
         maxprec=0
         epermut=0
         sameint=0
! not implemented for MQMQA
         empermut: do while(epermut.lt.endmemrec%noofpermut)
            epermut=epermut+1
!----------------------------------------------------------
! list endmembers
            call debug_endmemberpar(endmemrec,lut,lokph,msl,epermut,&
                 fractype,ceq)
!----------------------------------------------------------
            if(gx%bmperr.ne.0) goto 1100
!            write(*,*)'Listed excess parameters ',&
!                 epermut,endmemrec%noofpermut
! Excess parameters based on this endmember, also permutations ...
            pmq=1
            intlevel=0
            intrec=>endmemrec%intpointer
            interloop: do while(associated(intrec))
! return here if the interaction ?? maybe use cycle ? in gtp3X
200            continue
! each interaction has two pointer, to higher and to same level
! push link to same level on intstack and follow highlink
               intlevel=intlevel+1
               intstack(intlevel)%save=>intrec%nextlink
               pmq=intrec%order
!               write(*,*)'3EX interaction level ',intlevel
! this is label 220 also in gtp3X.F90
220            continue
               bford: if(chkperm) then
! complicated handling of permuted interaction parameter, see gtp3X,F90
! copied from gtp3X.F90
                  setipermut: if(maxpmq(pmq).eq.0) then
                     ipermut=1; lastpmq(pmq)=ipermut
! should I use gz%intlevel ??
                     maxpmq=intrec%noofip(intlevel)
                     plimit: if(ipermut.gt.maxpmq(pmq)) then
                        level: if(intlevel.eq.1) then
                           maxpmq(pmq)=maxpmq(pmq)+intrec%noofip(1)
                        elseif(intlevel.gt.2) then
                           write(*,*)'3EX permutation max intlevel 2'
                           gx%bmperr=4340; goto 1000
                        else
                           varying: if(intrec%noofip(1).eq.1) then
                              maxpmq(pmq)=maxpmq(pmq)+intrec%noofip(2)
                              if(ipermut.le.maxpmq(pmq)) goto 230
                           else
! complicated, see gtp3X at same place ....
                              nz=intrec%noofip(1)
                              if(maxpmq(pmq).gt.0) then
                                 if(intrec%noofip(1).eq.2) then
                                    maxpmq(pmq)=-maxpmq(pmq)
                                 else
                                    nz=mod(ipermut-1,intrec%noofip(1))
                                    if(nz.eq.0) then
                                       maxpmq(pmq)=maxpmq(pmq)
                                    else
                                       maxpmq(pmq)=maxpmq(pmq)+&
                                            intrec%noofip(nz+1)
                                    endif
                                 endif
                                 if(ipermut.le.maxpmq(pmq)) goto 230
                              else
                                 maxpmq(pmq)=intrec%noofip(2)-&
                                      maxpmq(pmq)
                                 if(ipermut.le.maxpmq(pmq)) goto 230
                              endif
                           endif varying
                        endif level
                        if(associated(intrec%highlink)) then
                           if(intlevel.eq.2) then
                              write(*,*)'3EX too high interaction'
                              gx%bmperr=4340; goto 1000
                           endif
                           goto 290
                        endif
                        if(intlevel.eq.0) exit interloop
                        pmq=intrec%order
                        nullify(intrec)
                        goto 295
                     endif plimit
230                  continue
                  endif setipermut
! this value of ipermut should vary with permutations 
                  write(*,*)'3EX ipermut: ',ipermut
                  lastpmq(pmq)=ipermut
               else
                  ipermut=1
               endif bford
!----------------------------------------------------------
! list excess parameters with permutationer
! subroutine debug_excesspar(intrec,lut,lokph,ipermut,intlevel,ceq)
               call debug_excesspar(intrec,lut,lokph,ipermut,intlevel,ceq)
               if(gx%bmperr.ne.0) goto 1000
!----------------------------------------------------------
! jump to 290 take higher excess several times for same lower order permutation
290            continue
! jump to 295 take higher excess ??
295            continue
! take highlink, if empty (not associated) pop intstack(level)
               intrec=>intrec%highlink
               nextint: do while(.not.associated(intrec))
!                  write(*,*)'3EX no higher link',intlevel
                  if(intlevel.le.0) then
                    exit nextint
                  endif
! intrec is set to the next interaction on same level, if empty decend furter
                  intrec=>intstack(intlevel)%save
                  intlevel=intlevel-1
               enddo nextint
            enddo interloop
! loop the whole data structure for each permutation of this endmember!!!
         enddo empermut
!-----------------------------------------
! Done all permutations and interaction of this endmember record
         endmemrec=>endmemrec%nextem
!         write(*,*)'Next endmember',associated(endmemrec)
      enddo endmemloop
! ?? check how this is used in calcg_internal in gtp3X.F90 ??
      write(*,*)'3EX Finished listing parameters for this phase'
      disord: if(fractype.eq.1 .and. btest(phlista(lokph)%status1,phmfs) &
           .and. btest(phmain%status2,csorder)) then
! we have traversing some of the parameter tree.
         returnoradd: if(first) then
! list (calculate) parameters in second (disordered) fraction set
            first=.false.
! in gtp3X we calculate for the first fraction set twice, now as disordered
! not needed for listing
!            goto 110
!         else
! here gtp3X returns to sums the disordered and ordered results and 
! maybe subtracs ordered as ordered.  Not needed when listing
         endif returnoradd
      else
! Here adding the two fraction sets, skip when listing 
      endif disord
! loop if more fraction types, fractype incremented when loop starts
   enddo fractyp
!----------------------------------------------
! original label ...
410 continue
!   fractionsets: if(btest(phlista(lokph)%status1,phmfs)) then
! both ordered and disordered listed above
!      write(*,*)'3EX Can be ignored?'
!   endif fractionsets
! finished this phase
1000 continue
   write(lut,1010)trim(phlista(lokph)%name)
1010 format(/'End of listing for phase ',a/60('-')/)
   return
! if error
1100 write(*,*)'Some error: ',gx%bmperr
   goto 1000
 end subroutine list_phaseparameters

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine debug_endmemberpar
!\begin{verbatim}
 subroutine debug_endmemberpar(endmemrec,lut,lokph,msl,epm,fractype,ceq)
! code to write a debug list of endmembers
! lokph is phlista index, lut is output unit
! cps is phase_varres record, ceq is equilibrium record
   implicit none
   integer lut,lokph,msl,epm,fractype
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_endmember), pointer :: endmemrec
   type(gtp_equilibrium_data), pointer :: ceq
!   type(gtp_phase_varres), target :: cps
!\end{verbatim}
   integer ll,idlist(9),id,ip,nsl
   character endmemconst*80
   character props*60
   TYPE(gtp_interaction), pointer :: intrec
   TYPE(gtp_pystack), pointer :: pystack
   TYPE(gtp_phase_add), pointer :: addrec
! to remember for fractset 2
   save nsl
!
!   write(*,*)'In debug_endmemberpar',msl
   endmemconst=':'
   ip=2
   if(fractype.eq.1) then
! >>>>>>>>> to be added: special for I2SL 
! save value of msl for use for fractype 2
      nsl=msl
      pyqloop: do ll=1,msl
! id is sequatial index of constituent over all sublattices ...
         id=endmemrec%fraclinks(ll,epm)
         idlist(ll)=id
         if(id.lt.0) then
            endmemconst(ip:)="*:"
            ip=ip+2
         else
! orderd or only fraction set has same constituents in sublattice 1 as ordered
            endmemconst(ip:)=&
                 trim(splista(phlista(lokph)%constitlist(id))%symbol)//':'
            ip=len_trim(endmemconst)+1
         endif
      enddo pyqloop
   else
! fractype 2, 
! If msl=2 then constituents in first and last sublattice, else just in first
! this is SPECIES index, alphabetical order.  Index in SPLISTA is SPECIES(..)
      idlist(1)=species(endmemrec%fraclinks(1,epm))
      endmemconst=':'//trim(splista(idlist(1))%symbol)//':'
      if(msl.eq.2) then
         ip=len_trim(endmemconst)+1
!  write(*,*)'3EX 2nd disordered sublattice: ',idlist(2),species(idlist(2))
         idlist(2)=species(endmemrec%fraclinks(2,epm))
         endmemconst(ip:)=trim(splista(idlist(2))%symbol)//':'
      endif
!      write(*,*)'Constituent: :',trim(endmemconst),idlist(1),idlist(2)
   endif
!
   proprec=>endmemrec%propointer
   write(*,*)'3EX Endmember: ',trim(endmemconst)
   write(lut,300,advance='no')trim(endmemconst)
300 format(5x,'Endmember ',a,i3)
! >>>>>>>> to be added: special for liquid2state   
   if(associated(proprec)) then
      props=' '
      id=1
      do while(associated(proprec))
         props(id:)=proprec%modelparamid
         id=len_trim(props)+2
         proprec=>proprec%nextpr
      enddo
! this is added to the line with the constituents
      write(lut,310)epm,trim(props)
310   format(i3,', MPIDs: ',a)
   else
      write(lut,320)epm
320   format(i3,', none')
   endif
1000 continue
   return
 end subroutine debug_endmemberpar

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable debug_excesspar
!\begin{verbatim}
 subroutine debug_excesspar(intrec,lut,lokph,ipermut,intlevel,ceq)
! code to debug list and endmember
! lokph is phlista index, lut is output unit
! cps is phase_varres record, ceq is equilibrium record
   implicit none
   integer lut,lokph,msl,epm,intlevel
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_interaction), pointer :: intrec
!   TYPE(gtp_endmember), pointer :: endmemrec
   type(gtp_equilibrium_data), pointer :: ceq
!   type(gtp_phase_varres), target :: cps
!\end{verbatim}
   integer nn,ll,isp,id
   character props*60
! from gtp3X:
   integer ipermut,intlat,ic,pmq
!
!   write(*,*)'3EX In debug_excesspar',intlevel,chkperm,ipermut
   intlat=intrec%sublattice(ipermut)
   ic=intrec%fraclink(ipermut)
   if(intlat.le.0 .or. ic.le.0) then
      write(*,*)'3EX illegal interaction constituent'
      gx%bmperr=4399; goto 1000
   endif
! ic is a sequential index to constituent fractions ... not species index
! extract the species number from phmain ...
! iliqsave and iliqva used in gtp3X but not here as we do not calculate anything
! but we can have wildcards ... maybe??
   isp=phlista(lokph)%constitlist(ic)
!   write(*,10)'3EX excess constituent/sublattice ',isp,ic,&
!        trim(splista(isp)%symbol),trim(splista(species(isp))%symbol)
!10 format(a,2i4,' symbol: ',a,' or ',a)
!
   proprec=>intrec%propointer
   props='MPIDs: '
   if(associated(proprec)) then
      id=8
      do while(associated(proprec))
         props(id:)=proprec%modelparamid
         id=len_trim(props)+2
         proprec=>proprec%nextpr
      enddo
   endif
   if(len_trim(props).eq.6) props='none'
   write(*,100)intlevel,trim(splista(isp)%symbol),intlat,trim(props)
   write(lut,100)intlevel,trim(splista(isp)%symbol),intlat,trim(props)
100 format(10x,'Excess level ',i3,'  ',a,'@',i1,', ',a)
!   bford: if(chkperm) then
! internal loop for for permutations for FCC/HCP and BCC with this endmember
!      write(*,*)'3EX this excess has permutations'
!   endif bford
1000 continue
   return
 end subroutine debug_excesspar

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
