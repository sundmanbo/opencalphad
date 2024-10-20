!
! gtp3EX included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     9B. Section: read and save on files using XML
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine read_xtdb
!\begin{verbatim}
 subroutine read_xtdb(filename,nel,selel)
! reading data from an XTFB file with selection of elements
!
! UNFINISHED REDUNDANT ... TO BE MODIFIED from reading PDB FORMAT to read XTDB
!   
   implicit none
   integer nel
   character filename*(*),selel(*)*2
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
   write(*,*)'3E in readxtdb! Not implemented yet'
   goto 1000
!
   emodel=0
   warning=.FALSE.
   silent=.FALSE.
   verbose=.TRUE.
!   if(btest(globaldata%status,GSSILENT)) then
!      silent=.TRUE.
!      write(*,*)'3E reading database silent'
!   endif
   write(*,*)'reading a XTDB file'
   if(.not.(index(filename,'.xtdb').gt.0 &
       .or. index(filename,'.XTDB').gt.0)) then
! no extention provided
      filename(len_trim(filename)+1:)='.XTDB'
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
! handle all XTDB tags 
120 continue
!   keyw=ispdbkeyword(line,nextc)
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
         if(ocv()) write(*,*)'3E Ignoring line: ',nl,ip,line(ip:ip+20)
      endif
      if(.not.onlyfun) write(*,*)'3E *** Ignoring unknown keyword: ',line(1:24)
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
      write(*,*)'3X read_xtb unfinished ....'
      stop
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
! THIS CODE NOT YET UPDATED
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
              '3E A CONSTITUENT keyword not directly preceeded by PHASE,',&
              ' line: ',nl
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
!      write(*,*)'readpdb 3E: ',ll,nr,nsl,longline(ip:jp)
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
!      write(*,379)'readpdb 3E: ',ip,nr,longline(ip:ip+10)
379   format(a,2i4,' >',a,'< >',a,'< >',a,'<')
      call getname(longline,ip,name3,mode,ch1)
!      write(*,379)'readpdb 3E: ',ip,nr,longline(ip:ip+10),name3,ch1
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
         if(modelname(1:5).eq.'I2SL ') then
            phtype='Y'
!            emodel=
         endif
!
! THIS IS READPDB subroutine .... to be changed to read_xtdb
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
         write(*,*)'3E this routine is wrong'
         stop
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
!      if(funname(lp1-1:lp1-1).eq.'D') then
!         fractyp=2
!         name1(lp1-1:lp1-1)=' '
!      else
!         fractyp=1
!      endif
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
!         call get_parameter_typty(name1,lokph,typty,fractyp)
         call get_parameter_typty(name1,lokph,typty)
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
! not disordered part, skip this parameter
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
         if(phlista(lokph)%noofsubl.eq.nsl) then
            fractyp=1
         elseif(phlista(lokph)%noofsubl.lt.nsl) then
! the disordered part has fewer sublattice ....
            fractyp=2
         else
            write(*,*)'Wrong number of sublattices 2',&
                 firsteq%phase_varres(lokcs)%disfra%ndd,nsl
            gx%bmperr=4399; goto 1000
         endif
! replaced this check by check above
!         if(fractyp.eq.1) then
!            if(phlista(lokph)%noofsubl.ne.nsl) then
!               write(*,*)'Wrong number of sublattices 1',&
!                    phlista(lokph)%noofsubl,nsl
!               gx%bmperr=4399; goto 1000
!            endif
!         else
!            if(firsteq%phase_varres(lokcs)%disfra%ndd.ne.nsl) then
!               write(*,*)'Wrong number of sublattices 2',&
!                    firsteq%phase_varres(lokcs)%disfra%ndd,nsl
!               gx%bmperr=4399; goto 1000
!            endif
!         endif
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
 end subroutine read_xtdb

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
   call date_and_time(date)
!   write(*,80)trim(XTDBversion),version(2:7),date(1:4),date(5:6),date(7:8)
   write(lut,80)trim(xtdbversion),version(2:7),date(1:4),date(5:6),date(7:8)
80 format('<?xml version="1.0"?>'/&
         '<?xml-model href="database.rng"',&
         ' schematypens="http://relaxng.org/ns/structure/1.0"',&
         ' type="application/xml"?>'/&
         '<Database version="',a,'">'/&
         '  <metadata>'/&
         '    <writer Software="OpenCalphad ',a,&
         '" Date="',a,'-',a,'-',a,'" />'/&
         '  </metadata>')
   nl=nl+6
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
! Writing in this order (option writing parameters by phase?)
! 1: Elements and species
! 2: TPfuns
! 3: Phases
! 4: Parameters
! 5: Bibliography
! 6: Models
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
!
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
!----------------------------- CHANGE HERE if gtp3_xml.F90 changes
! XTDBTPFUN=26
! 2: TPfuns  this is the tag xmlel(26) 
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
!----------------------------------------------------------------------
! in gtp3C.F90 there is a subroutine list_phase_data2 for TDB files
! 3: Phases
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
! DisorderedPart 
      if(btest(phlista(lokph)%status1,PHMFS)) then
         write(*,*)'3E The phase '//trim(phlista(lokph)%name)//&
              ' has a DisorderedPart'
! in gtp_fraction_set the data for disordered phase can be found
!         write(*,*)'3E Disordere:',phvarres%disfra%latd,phvarres%disfra%totdis
         if(phvarres%disfra%totdis.eq.1) then
            line='    <Disordered_3Part Sum="'
         else
            line='    <Disordered_2Part Sum="'
         endif
         if(btest(phlista(lokph)%status1,PHSORD)) then
! I am not sure how this is connected with Disordered_Part in OC
!            write(*,*)'3E do not subtract ordered as disordered'
         endif
         ip=len_trim(line)+1
         call wriint(line,ip,phvarres%disfra%latd)
         line(ip:)='" />'
!         write(*,10)'3E check: ',trim(line)
         write(lut,10)trim(line)
         nl=nl+1
!      else
!         write(*,*)'3E PHMFS bit not set'
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
            text512(ip:)='GLOWTEIN'; ip=ip+9
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
!      else
!         write(*,*)'3E this phase has no additions'
      endif
!------------------------------------------ suck ....
! CrystalStructure ?
      write(lut,290)
290   format('  </Phase>')
      nl=nl+1
!
!      write(*,*)'3E Finished this phase'
!      write(lut,*)'<!-- phase klar -->'
   enddo phaseloop
!   goto 900
!-----------------------------
! 4: Parameters, phase by phase
! code here copied from gtp3C.F90 subroutine list_phase_data
   do ia=1,noph()
!      write(*,*)'3E calling write_parametrar',lut
      lokph=phases(ia)
      call write_parameters(lokph,lut,2,nl)
      if(gx%bmperr.ne.0) goto 1100
!      write(*,*)'3E Finished listing parameters for this phase'
   enddo
!-----------------------------
! 5: Bibliography
!      write(*,*)'3E reffree: ',reffree
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
            call check_illegal_xml(line,ip)
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
! 6: Models
   if(includemodels) then
! if includemodels write all models, otherwise only those used
      call write_xtdbmodels1(lut,nl)
   endif
! finished
900 continue
   write(lut,990)
990 format('</Database>')
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
 subroutine check_illegal_xml(line,ip)
! replace < > & and " in references
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
      jp=index(line(1:ip),'"')
      if(jp.eq.0) exit quote
      line(jp:jp)='#'
   enddo quote
   return
 end subroutine check_illegal_xml

!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!!/!\!

! \addtotable subroutine exp2xml
 subroutine exp2xml(lut,expr,tag,nl)
! convert expr to TPfun/Parameter and Trange and write lines on unit lut
! tag is TPfun or Parameter and the first part of expr is the Id followd
! by the expression
   integer lut
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

!\addtotable subroutine write_xtdbmodels1
!\begin{verbatim}
 subroutine write_xtdbmodels1(lut,nl)
! write the model tags on unit lut with some explanations
   implicit none
   integer lut,nl
!\end{verbatim}
! Magnetic models
   write(lut,100)
100 format('  <Models>'/&
         4x,'This is a short explanation of XTDB model tags (or "elements") and their attributes,',&
         ' the models for the configurational entropy are not included.'/&
         4x,'The AmendPhase tag (nested inside a Phase tag) is used to specify some additional models for the phase'/&
         4x,'by using the attribute "Id" specified for most of the models below.'/&
         4x,'In these model tags there are model parameter identifiers (MPID) describing the dependence on composition, T and P.'/&
         4x,'A DisorderedPart tag must be nested inside the Phase tag as it has additional information.'/&
         4x,'The Toop and Kohler tags will normally appear together with model parameters for the binaries and ',&
         'has thus a phase attribute.'/&
         4x,'The EEC tag is global for the whole database if included.'/&
         4x,'Some model tags and MPIDs are tentative and some attributes of the tags are optional.')
   nl=nl+9
   write(lut,102)
102 format(4x,'<Magnetic Id="IHJBCC" MPID1="TC" MPID2="BMAGN" Aff=" -1.00" Bibref="82Her" > '/&
         6x,'f_below_TC= +1-0.905299383*TAO**(-1)-0.153008346*TAO**3-.00680037095*TAO**9-.00153008346*TAO**15; and'/&
         6x,'f_above_TC= -.0641731208*TAO**(-5)-.00203724193*TAO**(-15)-.000427820805*TAO**(-25); '/&
         6x,'in G=f(TAO)*LN(BMAGN+1) where TAO=T/TC.  Aff is the antiferromagnetic factor.'/&
         6x,'MPID1 is a combined Curie/Neel T and MPID2 the average Bohr magneton number.'/&
         4x,'</Magnetic>')
   nl=nl+6
!-----------------------------------------------
   write(lut,105)
105 format(4x,'<Magnetic Id="IHJREST"  MPID1="TC" MPID2="BMAGN" Aff=" -3.00" Bibref="82Her" > '/&
         6x,'f_below_TC= +1-0.860338755*TAO**(-1)-0.17449124*TAO**3-.00775516624*TAO**9-.0017449124*TAO**15; and '/&
         6x,'f_above_TC= -.0426902268*TAO**(-5)-.0013552453*TAO**(-15)-.000284601512*TAO**(-25); '/&
         6x,'in G=f(TAO)*LN(BMAGN+1) where TAO=T/TC.  Aff is the antiferromagnetic factor.'/&
         6x,'MPID1 is a combined Curie/Neel T and MPID2 the average Bohr magneton number.'/&
         4x,'</Magnetic>')
   nl=nl+6
!------------------------------------------ a lot of line-truncation errors
   write(lut,110)
110 format(4x,'<Magnetic Id="IHJQX" MPID1="CT" MPID2="NT" MPID3="BMAGN" Aff=" 0.00" Bibref="01Che 12Xio" > '/6x,&
         'f_below_TC= +1-0.842849633*TAO**(-1)-0.174242226*TAO**3-.00774409892*TAO**9-.00174242226*TAO**15-.000646538871*TAO**21;')
   write(lut,111)
111 format(6x,'f_above_TC= -.0261039233*TAO**(-7)-.000870130777*TAO**(-21)-.000184262988*TAO**(-35)-6.65916411E-05*TAO**(-49);'/&
         6x,'in G=f(TAO)*LN(BMAGN+1) where TAO=T/CT or T/NT.  Aff is the (redundant) antiferromagnetic factor.'/&
         6x,'MPID1 is the Curie T, MPID2 the Neel T and MPID3 the average Bohr magneton number.'/&
         4x,'</Magnetic>')
   nl=nl+5
!----------------------------------------------- "
   write(lut,120)
120 format(4x,'<Einstein Id="GLOWTEIN" MPID1="LNTH" Bibref="01Che" > '/&
         6x,'The Gibbs energy due to the Einstein low T vibrational model, G=1.5*R*THETA+3*R*T*LN(1-EXP(-THETA/T)).'/&
         6x,'MPID1 is the logarithm of the Einstein THETA.'/&
         4x,'</Einstein>')
   nl=nl+4
!-----------------------------------------------
   write(lut,125)
125 format(4x,'<Liquid2state Id="LIQ2STATE" MPID1="G2"  MPID2="LNTH" Bibref="88Agr 13Bec" > '/&
         6x,'Unified model for the liquid and the amorphous state which is treated as an Einstein solid.'/&
         6x,'MPID1 describes the stable liquid and the transition to the amorphous state and'/&
         6x,'MPID2 is the logarithm of the Einstein THETA for the amorphous phase.'/&
         4x,'</Liquid2state>')
   nl=nl+5
!-----------------------------------------------
   write(lut,130)
130 format(4x,'<Volume Id="VOLOWP" MPID1="V0"  MPID2="VA" MPID3="VB" Bibref="05Lu" > '/&
         6x,'The volume of a phase as function of T, moderate P and constitution via the model parameters:'/&
         6x,'MPID1 is the volume at the reference state, MPID2 is the integrated thermal expansion ',&
         'and MPID3 is the isothermal compressibilty at 1 bar.'/&
         4x,'</Volume>')
   nl=nl+4
!-----------------------------------------------
   write(lut,135)
135 format(4x,'<DisorderedPart Disordered=" " Subtract=" " Sum=" " Bibref="97Ans 07Hal" > '/&
         6x,'This tag is nested inside the ordered phase tag.  The disordered fractions are averaged over the ',&
         'number of ordered sublattices indicated by Sum.'/&
         6x,'The Gibbs energy is calculated separately for the ordered and disordered model parameters and added '/&
         6x,'but the configurational Gibbs energy is calculated only for the ordered phase.'/&
         6x,'If Subtract="Y" the Gibbs energy of the ordered phase is calculated a second time using the disordered ',&
         'fractions and subtracted.'/&
         6x,'Some software has no special disordered phase but all parameters are stored in the ordered one and'/&
         6x,'the parameters for the disordered phase has a suffix "D" (and different number of sublattices).'/&
         4x,'</DisorderedPart>')
   nl=nl+7
!-----------------------------------------------
   write(lut,140)
140 format(4x,'<Permutations Id="FCC4Perm" Bibref="09Sun" > '/&
         6x,'An FCC phase with 4 sublattices for the ordered tetrahedron use this model to indicate that parameters ',/&
         6x,'with permutations of the same set of constituents on identical sublattices are included only once.'/&
         4x,'</Permutations>')
   nl=nl+4
!-----------------------------------------------
   write(lut,142)
142 format(4x,'<Permutations Id="BCC4Perm" Bibref="09Sun" > '/&
         6x,'A BCC phase with 4 sublattices for the ordered asymmetric tetrahedron use this model to indicate that parameters '/&
         6x,'with permutations of the same set of constituents on identical sublattices are included only once.'/&
         4x,'</Permutations>')
   nl=nl+4
!-----------------------------------------------
   write(lut,145)
145 format(4x,'<EEC Id="EEC" Bibref="20Sun" > '/&
         6x,'The Equi-Entropy Criterion means that the software must ensure that solid phases with higher entropy than the ',&
         'liquid phase must not be stable. '/&
         4x,'</EEC>')
   nl=nl+3
!-----------------------------------------------
   write(lut,150)
150 format(4x,'<KohlerTernary Phase=" " Constituents=" " Bibref="01Pel" > '/&
         6x,'The symmetric Kohler model can be used for a specified ternary subsystem as described in the paper by Pelton.',/&
         6x,'The 3 constituents, separated by a space, can be in any order.'/&
         4x,'</KohlerTernary>') 
   nl=nl+4
!-----------------------------------------------
   write(lut,155)
155 format(4x,'<ToopTernary Phase=" " Constitutents=" " Bibref="01Pel" > '/&
         6x,'The asymmetric Toop model can be used for a specified ternary subsystem as described in the paper by Pelton.',/&
         6x,'The 3 constituents, separated by a space, must have the Toop constituent as the first one.'/&
         4x,'</ToopTernary>')
   nl=nl+4
!-----------------------------------------------
   write(lut,160)
160 format(4x,'<EBEF Id="EBEF" Bibref="18Dup" > '/&
         6x,'The Effective Bond Energy Formalism for phases with multiple sublattices using wildcards, "*", in the parameters',/&
         6x,'for sublattices with irrelevant constituents.  The parameters may also use the short form "constituent@sublattice" '/&
         6x,'in order to specify only the constituents in sublattices without wildcards.  It also requires ',&
         'the DisorderedPart model.'/&
         4x,'</EBEF>')
   nl=nl+5
!-----------------------------------------------
   write(lut,190)
190 format('  </Models>')
   nl=nl+1
!  Bibliography of models
   write(lut,200)
200 format('  <Bibliography> ')
   nl=nl+1
! the value of attribute "Id" in the Bibitem tag should
! be identical to the value of a "Bibref" attribute
   write(lut,210)
210 format('    <Bibitem Id="82Her" Text="S. Hertzman and B. Sundman,',&
         ' A Thermodynamic analysis of the Fe-Cr system,'&
         ' Calphad, Vol 6 (1982) pp 67-80" />'/&
         '    <Bibitem Id="88Agr" Text="J. Agren, Thermodynmaics of supercooled liquids and their glass transition,',&
         ' Phys Chem Liq, Vol 18 (1988) pp 123-139" />'/&
         '    <Bibitem Id="97Ans" Text="I. Ansara, N. Dupin, H. L. Lukas and B. Sundman, Thermodynamic assessment',&
         ' of the Al-Ni system, J All and Comp, Vol 247 (1997) pp 20-30" />'/&
         '    <Bibitem Id="01Che" Text="Q. Chen and B. Sundman, Modeling of Thermodynamic Properties for ',&
         'BCC, FCC, Liquid and Amorphous Iron, J Phase Eq, Vol 22 (2001) pp 631-644" />'/&
         '    <Bibitem Id="01Pel" Text="A. D. Pelton, A General Geometric Thermodynamic Model for Multicomponent ',&
         ' solutions, Calphad, Vol 25 (2001) pp 319-328" />'/&
         '    <Bibitem Id="05Lu" Text="X.-G. Lu, M. Selleby B. Sundman, Implementation of a new model for pressure',&
         ' dependence of condensed phases in Thermo-Calc, Calphad, Vol 29 (2005) pp 49-55" />'/&
         '    <Bibitem Id="07Hal" Text="B. Hallstedt, N. Dupin, M. Hillert, L. Hoglund, H. L. Lukas, J. C. Schuster',&
         ' and N. Solak, Calphad, Vol 31 (2007) pp 28-37" />'/&
         '    <Bibitem Id="09Sun" Text="B. Sundman, I. Ohnuma, N. Dupin, U. R. Kattner and S. G. Fries, An assessment',&
         ' of the entire Al-Fe system including D03 ordering, Acta Mater, Vol 57 (2009) pp 2896-2908" />'/&
         '    <Bibitem Id="12Xio" Text="W. Xiong, Q. Chen, P. A. Korzhavyi and M. Selleby, An improved magnetic',&
         ' model for thermodynamic modeling, Calphad, Vol 39 (2012) pp 11-20" />'/&
         '    <Bibitem Id="13Bec" Text="C. A. Becker, J. Agren, M. Baricco, Q Chen, S. A. Decterov, U. R. Kattner,',&
         ' J. H. Perepezko, G. R. Pottlacher and M. Selleby, Thermodynamic modelling of liquids: Calphad approaches',&
         ' and contributions from statistical physics, Phys Stat Sol B (2013) pp 1-20" />'/&
         '    <Bibitem Id="18Dup" Text="N. Dupin, U. R. Kattner, B. Sundman, M. Palumbo and S. G. Fries, Implementation of',&
         ' an Effective Bond Energy Formalism in the Multicomponent Calphad Approach, J Res NIST, Vol 123 (2018) 123020" />'/&
         '    <Bibitem Id="20Sun" Text="B. Sundman, U. R. Kattner, M. Hillert, M. Selleby, J. Agren, S. Bigdeli, Q. Chen,',&
         ' A. Dinsdale, B. Hallstedt, A. Khvan, H. Mao and R. Otis, A method for handling the extrapolation of solid crystalline',&
         ' phases to temperatures far above their melting point, Calphad, Vol 68 (2020) 101737" />'/&
         '  </Bibliography>')
   nl=nl+13
1000 continue
   return
 end subroutine write_xtdbmodels1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
! This is the index of PARAMETER in XMLEL
   integer, parameter :: xmlpartag=20
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
            call exp2xml(lut,funexpr,xmlel(xmlpartag),nl)
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
! redunant code
!         toop: if(associated(intrec%tooprec)) then
! Here we collect all Toop/Kohler tooprec and list them after the phase
!            write(*,*)'3E Found tooprecord',intrec%tooprec%uniqid
!            exit toop
!            tooprec=>intrec%tooprec
! tooplast is last tooprec that have had all links searched
!            tooplast=tooptop
!            tooplink=1
107         continue
!            write(*,210)'3E found toop record:',tooplink,tooptop,&
!                 tooprec%uniqid,&
!                 tooprec%toop,tooprec%const1,tooprec%const2,tooprec%const3,&
!                 tooprec%extra,associated(tooprec%next12),&
!                 associated(tooprec%next13),associated(tooprec%next23)
210         format(a,4i4,2x,3i3,2x,i5,' links: ',3l2)
! check if already found
!            do kk=1,tooptop
!               if(tooparray(kk)%p1%uniqid.eq.tooprec%uniqid) then
!                  write(*,*)'3E already found ',tooprec%uniqid
!                  tooplink=tooplink+1
!                  goto 108
!               endif
!            enddo
! a new tooprec, store it
!            tooptop=tooptop+1
!            tooparray(tooptop)%p1=>tooprec
!          write(*,*)'3E Adding a tooprec',tooptop,tooparray(tooptop)%p1%uniqid
!            tooplevel=0
! There are 3 links from each tooprecord, exact all records from these
! when empty go back one step of the stored tooprecord
! We do not have to go back further than tooplast.
108         continue
!            write(*,*)'3E tooprec links',tooplink,tooprec%uniqid
!            if(tooplink.eq.1) then
!               if(associated(tooprec%next12)) then
!                  tooprec=>tooprec%next12; goto 107
!               else
!                  tooplink=2
!               endif
!            endif
!            if(tooplink.eq.2) then
!               if(associated(tooprec%next13)) then
!                  tooprec=>tooprec%next13; goto 107
!               else
!                  tooplink=3
!               endif
!            endif
!            if(tooplink.eq.3 .and. associated(tooprec%next23)) then
!               tooprec=>tooprec%next23; goto 107
!            endif
! we have exhausted all links of this tooprec, go back to previous stored
! until we read tooplast
!            write(*,*)'3E no nextlinks: ',tooptop,tooplast,tooplevel
!            tooplevel=tooplevel+1
!            if(tooptop-tooplevel.gt.tooplast) then
!               write(*,*)'3E decreasing to ',tooptop-tooplevel,tooplast
!               tooprec=>tooparray(tooptop-tooplevel)%p1
!               tooplink=1
!               goto 108
!            endif
! We have checked all tooprecords linked from this binary
! Back to listing of parameters
!         endif toop
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
!            typspec=proprec%proptype
!            if(typspec.gt.100) then
!               typty=typspec/100
!               typspec=mod(typty,100)
!            else
!               typty=typspec
!            endif
            if(typty.gt.0 .and. typty.le.ndefprop) then
               prop=propid(typty)%symbol
!               if(parlist.eq.2) then
! disordered interaction parameter DO NOT ADD SUFFIX D any longer
!                  kk=len_trim(prop)+1
!                  prop(kk:kk)='D'
!               endif
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
                  call exp2xml(lut,funexpr,xmlel(xmlpartag),nl)
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


 
