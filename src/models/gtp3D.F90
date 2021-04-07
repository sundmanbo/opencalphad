!
! gtp3D included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     8. Interactive things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine ask_phase_constitution
!\begin{verbatim}
 subroutine ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
! interactive input of a constitution of phase iph
   implicit none
   integer last,iph,ics,lokcs
   character cline*(*)
!\end{verbatim} %+
! NOTE a strange bug when calculating a phase 
! the result is different if one sets the constitution explicitly the same!!
   character name1*24,quest*32
   double precision yarr(maxcons2),sites(maxsubl),qq(5),yyy,xxx,sss,ydef
   integer knl(maxsubl),knr(maxcons2)
   character line*64,ch1*1,crest*24
   character :: lastph*24='                        '
! changed default to N
!   character*1 :: chd='Y'
   character*1 :: chd='N'
   integer qph,lokph,nsl,kkk,loksp,ip,ll,nr,yrest
   TYPE(gtp_equilibrium_data), pointer :: ceq
   logical once
! save here to use the same default as last time
   save chd,lastph
   call gparcdx('Phase name: ',cline,last,1,name1,lastph,'?Amend phase constit')
   if(name1(1:2).eq.'* ') then
! this means all phases and composition sets
! If iph is -1 than this is not allowed!!
      if(iph.lt.0) then
         write(kou,*)'Wildcard not allowed in this case'
         goto 1000
      endif
      qph=-1
      iph=1
      ics=1
      call get_phase_name(iph,ics,name1)
      if(gx%bmperr.ne.0) goto 1000
   else
      qph=0
      call find_phase_by_name(name1,iph,ics)
      if(gx%bmperr.ne.0) goto 1000
! remember the phase name
      lastph=name1
   endif
100 continue
!   write(*,*)'3D spc 1',qph,iph,ics,name1
! skip hidden and suspended phases, test_phase_status return
! -4 hidden, -3 suspend, -2 dormant, -1,0, entered, 2 fixed
   if(qph.lt.0 .and. test_phase_status(iph,ics,xxx,ceq).le.PHDORM) goto 200
!   if(qph.lt.0 .and. (phase_status(iph,ics,PHHID,ceq) .or.&
!        phase_status(iph,ics,PHIMHID,ceq) .or.&
!        (phase_status(iph,ics,CSSUS,ceq) .and. &
!        .not.phase_status(iph,ics,CSFIXDORM,ceq)))) goto 200
!   lokph=phases(iph)
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   call get_phase_data(iph,ics,nsl,knl,knr,yarr,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! ask for amount of formula units, default is current amount
   yyy=ceq%phase_varres(lokcs)%amfu
   quest='Amount of '//name1
   call gparrdx(quest,cline,last,xxx,yyy,'?Amend phase constit')
! if input error quit asking more
   if(buperr.ne.0) then
      buperr=0; goto 1000
   endif
   ceq%phase_varres(lokcs)%amfu=xxx
! ask if we should set the current constitution, ignore default
!   write(*,*)'3D we are here!'
   call gparcdx('Current (Y), default (D) or new (N) constitution?',&
        cline,last,1,ch1,chd,'?Amend phase constit')
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      chd='Y'
! set the old constitution explicitly!!
! without this seemingly unnecessary call to set the same constution the
! calculate phase gives sometimes wrong values!
      call set_constitution(iph,ics,yarr,qq,ceq)
      goto 200
   elseif(ch1.eq.'d' .or. ch1.eq.'D') then
      chd='D'
      call set_default_constitution(iph,ics,ceq)
      goto 200
   else
! constitution entered interactivly
      chd='N'
   endif
! ask for constitution
   kkk=0
   nylat: do ll=1,nsl
      yrest=0
      sss=one
      ydef=one
      if(knl(ll).eq.1) then
         kkk=kkk+1; cycle nylat
      endif
      nycon: do nr=1,knl(ll)
         if(nr.eq.knl(ll) .and. yrest.eq.0) then
            cycle nycon
         endif
         kkk=kkk+1
         loksp=phlista(lokph)%constitlist(kkk)
         line='Fraction of '//splista(loksp)%symbol
         ip=len_trim(line)+1
         if(ll.gt.1) then
            line(ip:)='#'//char(ll+ichar('0'))
            ip=ip+2
         endif
         once=.true.
20        continue
         ydef=min(yarr(kkk),ydef)
         call gparrdx(line(1:ip+2),cline,last,xxx,ydef,'?Amend phase constit')
         if(buperr.ne.0) then
!            write(*,*)'3D Allow REST: ',trim(cline),last,buperr,yrest
            buperr=0
            if(yrest.eq.0) then
               crest=cline(last:last+3)
               call capson(crest)
               if(crest(1:4).eq.'REST') then
                  yrest=nr
                  last=len(cline)
                  crest=splista(loksp)%symbol
                  cycle nycon
               endif
            endif
         endif
         if(xxx.lt.zero) then
            if(once) then
               write(*,*)'A fraction must be greater than zero'
               yarr(kkk)=1.0D-12
               once=.false.
               goto 20
            else
               gx%bmperr=4146; goto 1000
            endif
         endif
         sss=sss-xxx
         if(sss.lt.zero) then
            xxx=max(sss+xxx,1.0D-12)
            sss=-1.0D12
            write(*,21)'Sum of fractions larger 1.0, fraction set to: ',xxx
21          format(a,1pe12.4)
            ydef=1.0D-12
         else
            ydef=sss
         endif
!         write(*,*)'ydef: ',ydef,sss
         yarr(kkk)=xxx
      enddo nycon
! if yrest is zero the last constituent is set to the rest, otherwise yrest
      if(yrest.eq.0) then
         kkk=kkk+1
         yarr(kkk)=max(sss,1.0D-12)
         write(*,21)'Last fraction set to: ',yarr(kkk)
      else
         yarr(yrest)=max(sss,1.0D-12)
         write(*,21)'Fraction of '//trim(crest)//' set to: ',yarr(yrest)
      endif
   enddo nylat
! set the new constitution
   call set_constitution(iph,ics,yarr,qq,ceq)
! if all phases loop
200 continue
   if(qph.lt.0) then
      if(gx%bmperr.eq.4050) then
! error no such phase, quit
         gx%bmperr=0; goto 1000
      elseif(gx%bmperr.eq.4072) then
! error no such composition set, take next phase
         gx%bmperr=0
         iph=iph+1
         ics=1
      else
         ics=ics+1
      endif
      call get_phase_name(iph,ics,name1)
      if(gx%bmperr.ne.0) goto 200
      goto 100
   endif
1000 continue
! return -1 as phase number of loop for all phases made
   if(qph.lt.0) iph=-1
   return
 end subroutine ask_phase_constitution

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine ask_phase_new_constitution
!\begin{verbatim} %-
 subroutine ask_phase_new_constitution(cline,last,iph,ics,lokcs,ceq)
! interactive input of a constitution of phase iph
   implicit none
   integer last,iph,ics,lokcs
   character cline*(*)
!\end{verbatim}
   character name1*24,quest*32
   double precision yarr(maxcons2),sites(maxsubl),qq(5),yyy,xxx,sss,ydef
   integer knl(maxsubl),knr(maxcons2)
   character line*64,ch1*1
   character*1 :: chd='Y'
   integer qph,lokph,nsl,kkk,loksp,ip,ll,nr
   TYPE(gtp_equilibrium_data), pointer :: ceq
   logical once
! save here to use the same default as last time
   save chd
   qph=0
!   call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
!   if(name1(1:2).eq.'* ') then
! this means all phases and composition sets
!      qph=-1
!      iph=1
!      ics=1
!      call get_phase_name(iph,ics,name1)
!      if(gx%bmperr.ne.0) goto 1000
!   else
!      qph=0
!      call find_phase_by_name(name1,iph,ics)
!      if(gx%bmperr.ne.0) goto 1000
!   endif
100 continue
!   write(*,*)'3D spc 1',qph,iph,ics,name1
! skip hidden and suspended phases, test_phase_status return
! -4 hidden, -3 suspend, -2 dormant, -1,0, entered, 2 fixed
!   if(qph.lt.0 .and. test_phase_status(iph,ics,xxx,ceq).le.PHDORM) goto 200
!   if(qph.lt.0 .and. (phase_status(iph,ics,PHHID,ceq) .or.&
!        phase_status(iph,ics,PHIMHID,ceq) .or.&
!        (phase_status(iph,ics,CSSUS,ceq) .and. &
!        .not.phase_status(iph,ics,CSFIXDORM,ceq)))) goto 200
!   lokph=phases(iph)
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   call get_phase_data(iph,ics,nsl,knl,knr,yarr,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! ask for amount of formula units, default is current amount
!   yyy=ceq%phase_varres(lokcs)%amfu
!   quest='Amount of '//name1
! NOTE name is not set!
!   call gparrd(quest,cline,last,xxx,yyy,q1help)
!   if input error quit asking more
!   if(buperr.ne.0) then
!      buperr=0; goto 1000
!   endif
!   ceq%phase_varres(lokcs)%amfu=abs(xxx)
! ask if we should set the default constitution
!   write(*,*)'3D we are really here?'
   call gparcdx('Default constitution?',cline,last,1,ch1,chd,&
        '?Amend phase constit')
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      call set_default_constitution(iph,ics,ceq)
      if(gx%bmperr.ne.0) goto 1000
      chd='Y'
      goto 200
   else
      chd='N'
   endif
! ask for constitution
   kkk=0
   nylat: do ll=1,nsl
      sss=one
      ydef=one
      nycon: do nr=1,knl(ll)-1
         kkk=kkk+1
         loksp=phlista(lokph)%constitlist(kkk)
         line='Fraction of '//splista(loksp)%symbol
         ip=len_trim(line)+1
         if(ll.gt.1) then
            line(ip:)='#'//char(ll+ichar('0'))
            ip=ip+2
         endif
         once=.true.
20        continue
         ydef=min(yarr(kkk),ydef)
         call gparrdx(line(1:ip+2),cline,last,xxx,ydef,'?Amend phase constit')
         if(xxx.lt.zero) then
            if(once) then
               write(*,*)'A fraction must be greater than zero'
               yarr(kkk)=1.0D-12
               once=.false.
               goto 20
            else
               gx%bmperr=4146; goto 1000
            endif
         endif
         sss=sss-xxx
         if(sss.lt.zero) then
            xxx=max(sss+xxx,1.0D-12)
            sss=-1.0D12
            write(*,21)'Sum of fractions larger 1.0, fraction set to: ',xxx
21          format(a,1pe12.4)
            ydef=1.0D-12
         else
            ydef=sss
         endif
!         write(*,*)'ydef: ',ydef,sss
         yarr(kkk)=xxx
      enddo nycon
! the last constituent is set to the rest
      kkk=kkk+1
      yarr(kkk)=max(sss,1.0D-12)
      write(*,21)'Last fraction set to: ',yarr(kkk)
   enddo nylat
! set the new constitution
   call set_constitution(iph,ics,yarr,qq,ceq)
! if all phases loop
200 continue
!   if(qph.lt.0) then
!      if(gx%bmperr.eq.4050) then
! error no such phase, quit
!         gx%bmperr=0; goto 1000
!      elseif(gx%bmperr.eq.4072) then
! error no such composition set, take next phase
!         gx%bmperr=0
!         iph=iph+1
!         ics=1
!      else
!         ics=ics+1
!      endif
!      call get_phase_name(iph,ics,name1)
!      if(gx%bmperr.ne.0) goto 200
!      goto 100
!   endif
1000 continue
! return -1 as phase number of loop for all phases made
   if(qph.lt.0) iph=-1
   return
 end subroutine ask_phase_new_constitution

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine enter_parameter_interactivly
!\begin{verbatim}
 subroutine enter_parameter_interactivly(cline,ip,mode)
! enter a parameter from terminal or macro
! NOTE both for ordered and disordered fraction set !!
! mode = 0 for entering
!        1 for listing on screen (kou)
   implicit none
   integer ip,mode
   character cline*(*)
!\end{verbatim}
   character name1*24,name2*24,longline*256,refx*16,elnam*24
! funame only 16 first characters will be used
   character name3*64,ch1*1,line*64,parname*64,funame*24
   integer typty,lint(2,5),fractyp,typty1,kp,lp1,kel,kq,iel,isp,lk3,lp2
   integer jph,ics,lokph,ll,k4,nint,jp,lsc,ideg,kk,lfun,nsl,loksp
   integer, dimension(maxsubl) :: endm(maxsubl)
   double precision xxx
!
   lk3=0
10  continue
   call gparcx('Parameter name: ',cline,ip,7,parname,' ','?Enter parameter')
! simple parameter names are like G(SIGMA,FE:CR:FE,CR;1)
   kp=index(parname,' ')
   parname(kp:)=' '
! extract symbol, normally G or L but TC and others can occur 
! for example a mobility like  MQ&FE+2#3 where FE+2#3 is a constinuent
! in sublattice 3
! NO ABBREVIATION IS ACCEPTED, for example not BM for BMAGN
   lp1=index(parname,'(')
!   write(*,*)'3D in parname: ',trim(parname),lp1
   if(lp1.le.1) then
      gx%bmperr=4027; goto 1000
   endif
! name1 is everything up to (
   name1=parname(1:lp1-1)
   call capson(name1)
! It can be a mobility with a & inside
   kel=index(name1,'&')
   if(kel.gt.0) then
! note that elnam may contain sublattice specification like Fe+2#2
      elnam=name1(kel+1:)
      name1=name1(1:kel-1)
!      write(*,*)'3D elnam: ',elnam
   endif
   kq=len_trim(name1)
!   write(*,*)'3D: fractyp: ',kq,name1(1:kq)
   if(name1(kq:kq).eq.'D') then
! A final "D" on the paramer symbol indicates fractyp=2
      name1(kq:kq)=' '
      fractyp=2
   else
      fractyp=1
   endif
! find the property associated with this symbol
   do typty=1,ndefprop
!      write(*,*)'Property symbol: "',propid(typty)%symbol,'"'
      if(name1(1:4).eq.propid(typty)%symbol) then
         goto 70
      endif
   enddo
! no matching symbol
   write(kou,*)'3D unknown parameter type, please reenter: ',&
        name1(1:len_trim(name1))
   parname=' '; goto 10
! typty is the parameter symbol index
70 continue
   typty1=typty
   iel=0; isp=0
! the beginning of the TP function name
!   funame='_'//propid(typty1)%symbol(1:1)
   funame='_'//propid(typty1)%symbol
   if(kel.gt.0) then
! there is a specifier, check if correct element or species
      kel=index(elnam,'#')
      if(kel.gt.0) then
! extract sublattice number 1-9 specification
         lk3=ichar(elnam(kel+1:kel+1))-ichar('0')
!         write(*,73)elnam(kel+1:kel+1),kel,elnam,lk3
!73       format('3D sublattice: "',a,'" position: ',i3,' in ',a,' : ',i3)
         elnam(kel:)=' '
      endif
      if(btest(propid(typty)%status,IDELSUFFIX)) then
!         write(*,*)'3D: elnam: ',kel,lk3,typty,elnam
         call find_element_by_name(elnam,iel)
         if(gx%bmperr.ne.0) then
            write(kou,*)'3D Unknown element ',elnam,&
                 ' in parameter type MQ, please reenter'
            goto 1000
!            parname=' '; gx%bmperr=0; goto 10
         endif
         typty=100*typty+iel
      elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
! to know the constituents we must know the phase but as we do not know 
! the phase name yet but check the species exists !!!
!         write(*,*)'3D: conname: ',kel,lk3,typty,elnam
         call find_species_by_name(elnam,isp)
         if(gx%bmperr.ne.0) then
! This is not an error, the species may simply not be selected !!!
            write(kou,*)'Unknown species ',trim(elnam),&
                 ' in parameter type MQ, please reenter',gx%bmperr
            goto 1000
!            parname=' '; gx%bmperr=0; goto 10
         endif
! convert from index to location, loksp
         loksp=species(isp)
         if(lk3.eq.0) then
! sublattice after # saved in lk3 above, if none (0) assume 1
            lk3=1
         endif
      else
         write(kou,*)'3D This model parameter identifier has no specifier'
         gx%bmperr=4168; goto 1000
      endif
! this is the property type stored in property record
   else
! check if there should be a specifier !!
      if(btest(propid(typty)%status,IDELSUFFIX) .or. &
           btest(propid(typty)%status,IDCONSUFFIX)) then
         write(*,*)'3D Parameter specifier missing'
         gx%bmperr=4169; goto 1000
      endif
   endif
! 4027?
! extract phase name and constituent array
   lp1=index(parname,'(')
   lp2=index(parname,',')
   if(lp2.lt.lp1) then
      gx%bmperr=4028; goto 1000
   endif
   name2=parname(lp1+1:lp2-1)
!    write(*,*)'enter_parameter_inter 1: ',lp1,lp2,name2
   call find_phase_by_name_exact(name2,jph,ics)
   if(gx%bmperr.ne.0) then
! special case for reference phase
      gx%bmperr=0;
      call capson(name2)
      if(name2.eq.'SELECT_ELEMENT_REFERENCE') then
         jph=0; ics=1
      else
         write(kou,*)'Unknown phase name, please reenter'
         kp=len(cline)
         goto 10
      endif
   endif
   lokph=phases(jph)
! add the full phase name to the function name.  Remove any _ or numbers ,,,
!   ll=len_trim(funame)+1
!   funame(3:)=phlista(lokph)%name
! only save first letter of parameter type, no problem with duplicate names
   ll=3
   funame(ll:)=phlista(lokph)%name
!   write(*,*)'3D funame 1: ',trim(funame),', ',name2
   ll=4
74 continue
   if(funame(ll:ll).eq.'_') then
      funame(ll:)=funame(ll+1:)
      goto 74
   elseif(ll.lt.9) then
      ll=ll+1
      goto 74
   endif
! eliminate anything from position 9
   funame(9:)=' '
!   write(*,*)'3D funame 2: ',trim(funame)
! if the parameter symbol has a constituent specification check that now
!   write(*,*)'3D lk3 and isp: ',lk3,isp
   if(lk3.gt.0 .and. isp.gt.0) then
! No check for elements ...
      k4=0
      do ll=1,phlista(lokph)%noofsubl
! careful ll is double letter l, not 11 (eleven)
         if(lk3.eq.0 .or. lk3.eq.ll) then
            do kk=1,phlista(lokph)%nooffr(ll)
               k4=k4+1
               if(phlista(lokph)%constitlist(k4).eq.loksp) goto 80
            enddo
         elseif(ll.lt.lk3) then
            k4=k4+phlista(lokph)%nooffr(ll)
         endif
      enddo
! constituent not found
      write(kou,*)'3D Parameter symbol contains unknown constituent'
      gx%bmperr=4066; goto 1000
! constituent found in right sublattice
80    continue
      typty=100*typty+k4
!      write(*,81)'3D: found: ',typty1,typty,lk3,k4,loksp
81    format(a,10i4)
   endif
!    write(*,*)'enter_parameter_inter 2: ',jph,lokph
! extract constituent array, remove final ) and decode
   name3=parname(lp2+1:)
   lp1=len_trim(name3)
! this removes the final )
   name3(lp1:)=' '
!
!   write(*,*)'3D decoding constituent array'
   call decode_constarr(lokph,name3,nsl,endm,nint,lint,ideg)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,83)'3D after d_c: ',name3(1:lp1),nint,(lint(2,kp),kp=1,nint)
83 format(a,a,i5,2x,5i4)
   kp=len_trim(funame)
   funame(kp+1:)=name3
   call capson(funame)
! finally remove all non-alphabetical characters in the function name by _
   kp=3
100 continue
   kp=kp+1
105 continue
!   ch1=parname(kp:kp)
   ch1=funame(kp:kp)
! should use ??
!   if(ucletter(ch1)) goto 100
   if(ch1.ge.'A' .and. ch1.le.'Z') goto 100
   if(ch1.ne.' ') then
!      parname(kp:)=parname(kp+1:)
      funame(kp:)=funame(kp+1:)
      if(kp.lt.16) goto 105
   endif
   funame(17:)=' '
   kp=len_trim(funame)
   if(kp.lt.16) funame(kp+1:kp+1)=char(ideg+ichar('0'))
!   write(*,*)'3D funame 3: ',trim(funame),', ',trim(name3)
!   parname='_'//parname
!-------------------------------------------------
! if mode=0 enter the parameter, 
! if mode=1 just list the parameter
! if mode=2 maybe amending (does FOOLED) work?
   if(mode.eq.1) then
      lfun=-1
!      write(*,*)'3D calling enter_parameter with lfun=',lfun
      call enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,&
           lfun,refx)
! this error means illegal reference ... 
! irrelevant but I am not sure where it is set ...
      if(gx%bmperr.eq.4154) gx%bmperr=0
      goto 1000
   endif
! continue here to enter the parameter
! If parameter has no T dependendence just ask for value
   if(btest(propid(typty1)%status,IDNOTP)) then
      write(kou,*)'This parameter can only be a constant'
      call gparrx('Value: ',cline,ip,xxx,zero,'?Enter parameter')
      if(buperr.ne.0) then
         xxx=zero; buperr=0
      endif
! the tpfun always want a low-T, expression; high-T N
      write(longline,110)xxx
110   format(' 1 ',1pe16.7,'; 20000 N ')
      jp=len_trim(longline)+2
      goto 200
   endif
   if(btest(propid(typty1)%status,IDONLYP)) then
      write(kou,*)'This parameter may not depend on T, only on P'
   endif
!-------------------------------------------------
! now read the function.
   call gparrx('Low  temperature limit /298.15/:',cline,ip,xxx,2.9815D2,&
        '?Enter parameter')
   if(buperr.ne.0) then
      buperr=0; longline=' 298.15 '
      jp=8
   else
      longline=' '
      jp=1
      call wrinum(longline,jp,8,0,xxx)
      if(buperr.ne.0) goto 1000
      jp=jp+1
   endif
!    write(*,152)-1,jp,longline(1:jp)
!-----------------------------------------------
! return here for new expression in another range
   lsc=1
115 continue
   call gparcx('Expression, end with ";":',cline,ip,6,line,';',&
        '?Enter parameter')
   if(buperr.ne.0) then
      buperr=0; line=';'
   endif
120 continue
   longline(jp:)=line
   jp=len_trim(longline)+1
!   write(*,152)0,jp,longline(1:jp)
   if(index(longline(lsc:),';').le.0) then
      call gparcx('&',cline,ip,6,line,';','?Enter parameter')
      if(buperr.ne.0) then
         buperr=0; line=';'
      endif
      goto 120
!   else
!      write(*,*)'Found ; at ',index(longline,';')
   endif
150 continue
   jp=jp+1
!   write(*,152)0,jp,longline(1:jp)
! lsc is positioned after the ; of previous ranges
   lsc=jp
!    write(*,152)1,ip,cline(1:ip)
   call gparrx('Upper temperature limit /6000/:',cline,ip,xxx,6.0D3,&
        '?Enter parameter')
   if(buperr.ne.0) then
      buperr=0; xxx=6.0D3
   endif
   call wrinum(longline,jp,8,0,xxx)
   if(buperr.ne.0) goto 1000
   call gparcdx('Any more ranges',cline,ip,1,ch1,'N','?Enter parameter')
   if(ch1.eq.'n' .or. ch1.eq.'N') then
      longline(jp:)=' N'
      jp=jp+3
   else
      longline(jp:)='Y'
      jp=jp+2
      goto 115
   endif
! jump here for parameters that are constants
200 continue
   call gparcdx('Reference symbol:',cline,ip,1,refx,'UNKNOWN',&
        '?Enter parameter')
   call capson(refx)
   longline(jp:)=refx
   jp=len_trim(longline)+1
!    write(*,252)2,jp,longline(1:jp)
252 format('3D ep: ',2i3,'>',a,'<')
!
   call capson(longline(1:jp))
!   write(*,*)'3D epi: ',longline(1:jp)
!   call enter_tpfun(parname,longline,lfun,.FALSE.)
!   write(*,*)'3D funame: ',trim(funame)
!   call store_tpfun(funame,longline,lfun,.FALSE.)
! last argumnent -1 means not reading from TDB file
   call store_tpfun(funame,longline,lfun,-1)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,290)'3D enter_par 7: ',lokph,nsl,nint,ideg,lfun,refx
290 format(a,5i4,1x,a)
!
   call enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,lfun,refx)
!
1000 continue
   return
 end subroutine enter_parameter_interactivly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine amend_global_data
!\begin{verbatim}
 subroutine amend_global_data(cline,ipos)
   implicit none
   character cline*(*)
   integer ipos
!\end{verbatim}
   character name*24,current*24,ch1*1,chd*1
   current=globaldata%name
!   write(*,*)'entering amend_global_data: ',cline(1:30)
   call gparcdx('System name: ',cline,ipos,1,name,current,'?Amend general')
   if(proper_symbol_name(name,0)) then
      globaldata%name=name
   else
      write(kou,*)'Illegal name ignored'
      goto 1000
   endif
100 continue
   chd='N'
   if(btest(globaldata%status,GSBEG)) then
      chd='B'
   elseif(btest(globaldata%status,GSADV)) then
      chd='E'
   else
      chd='F'
   endif
   call gparcdx('I am a beginner (B), freqent user (F) or expert (E): ',&
        cline,ipos,1,ch1,chd,'?Amend general')
   call capson(ch1)
   globaldata%status=ibclr(globaldata%status,GSBEG)
   globaldata%status=ibclr(globaldata%status,GSADV)
   globaldata%status=ibclr(globaldata%status,GSOCC)
   if(ch1.eq.'B') then
      globaldata%status=ibset(globaldata%status,GSBEG)
   elseif(ch1.eq.'E') then
      globaldata%status=ibset(globaldata%status,GSADV)
   else
      globaldata%status=ibset(globaldata%status,GSOCC)
   endif
120 continue
! is global minimization allowed?
   chd='Y'
   if(btest(globaldata%status,GSNOGLOB)) chd='N'
   call gparcdx('Global gridminimization allowed: ',&
        cline,ipos,1,ch1,chd,'?Amend general')
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      globaldata%status=ibclr(globaldata%status,GSNOGLOB)
   else
      globaldata%status=ibset(globaldata%status,GSNOGLOB)
   endif
! allow merging gridpoints after global?
   chd='Y'
   if(btest(globaldata%status,GSNOMERGE)) chd='N'
   call gparcdx('Merging gridpoints in same phase allowed: ',&
        cline,ipos,1,ch1,chd,'?Amend general')
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      globaldata%status=ibclr(globaldata%status,GSNOMERGE)
   else
      globaldata%status=ibset(globaldata%status,GSNOMERGE)
   endif
! GSNOACS can be changed interactivly, 0 means allowed
   chd='Y'
   if(btest(globaldata%status,GSNOACS)) chd='N'
   call gparcdx('Composition sets can be created automatically? ',&
        cline,ipos,1,ch1,chd,'?Amend general')
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      globaldata%status=ibclr(globaldata%status,GSNOACS)
   else
      globaldata%status=ibset(globaldata%status,GSNOACS)
   endif
! GSNOREMCS can be changed interactivly, 0 means not remove
   chd='Y'
   if(btest(globaldata%status,GSNOREMCS)) chd='N'
   call gparcdx('Delete unnecessary composition sets automatically? ',&
        cline,ipos,1,ch1,chd,'?Amend general')
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      globaldata%status=ibclr(globaldata%status,GSNOREMCS)
   else
      globaldata%status=ibset(globaldata%status,GSNOREMCS)
   endif
1000 continue
   return
 end subroutine amend_global_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine enter_bibliography_interactivly
!\begin{verbatim}
 subroutine enter_bibliography_interactivly(cline,last,mode,iref)
! enter a reference for a parameter interactivly
! mode=0 means enter, =1 amend
   implicit none
   character cline*(*)
   integer last,mode,iref
   logical twotries
!\end{verbatim}
! stupid with a variable called L80
   character line*256,refid*16,L80*80
   integer jl,ip
   call gparcx('Reference identifier:',cline,last,1,refid,' ',&
        '?Amend bibliography')
   if(buperr.ne.0 .or. refid(1:1).eq.' ') then
!      write(kou,*)'There must be an identifier'
      gx%bmperr=4155; goto 1000
   endif
   call capson(refid)
! check if unique, if mode=0 illegal
   do jl=1,reffree-1
      if(refid.eq.bibrefs(jl)%reference) then
         if(mode.eq.0) then
!            write(kou,*)'Reference identifier not unique'
            gx%bmperr=4156;goto 1000
         else
            goto 70
         endif
      endif
   enddo
! if mode=1 one should have found the reference
   if(mode.eq.1) then
      write(kou,*)'No such reference'
      goto 1000
   endif
70 continue
   ip=1
   line=' '
   twotries=.TRUE.
100 continue
!   call gparc('Reference text, end with ";":',cline,last,5,l80,';',q1help)
   call gparcx('Reference text, end with ";":',cline,last,5,l80,';',&
        '?Amend bibliography')
   line(ip:)=l80
   ip=len_trim(line)
   if(ip.le.1 .and. twotries) then
      twotries=.FALSE.
      write(kou,*)'There must be some bibilograpic text!'
      ip=1; goto 100
   elseif(line(ip:ip).ne.';') then
      twotries=.FALSE.
      write(*,*)'Terminate text with a ";"'
      ip=ip+1; goto 100
   elseif(.not.twotries) then
      if(mode.eq.1) then
         write(*,*)'Bibliogaphic reference unchanged'
      else
         write(*,*)'Bibliogaphic reference not entered'
      endif
      goto 1000
   else
      line(ip:)=' '
   endif
   call tdbrefs(refid,line,1,iref)
1000 continue
   return
 end subroutine enter_bibliography_interactivly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine enter_experiment
!\begin{verbatim}
 subroutine enter_experiment(cline,ip,ceq)
! enters an experiment, almost the same as set_condition   
   implicit none
   character cline*(*)
   integer ip
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
! New is set to the new condition or experiment
! in set_condition new is not used for anything.
! in this subroutine the new variable is removed from the condition list
! and instead added to the experimenal list
   integer kp,jc,istv,qp
   type(gtp_condition), pointer :: new,temp
!   integer nidlast,nidfirst,nidpre
   double precision xxx,yyy
   character usymbol*16,ch1*1
! do not allow experiments in first equilibrium!!
   if(ceq%eqno.eq.1) then
      write(kou,16)
16    format('Experiments are not allowed in the default equilibrium')
      goto 1000
   endif
!
! return here if more experiments
17 continue
! inside here things are done
!   write(*,*)'3D exp1: ',trim(cline),ip
   call set_cond_or_exp(cline,ip,new,1,ceq)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,'(a,a,2i4)')'3D exp2: ',trim(cline),ip,new%active
   if(new%active.ne.1) then
! the experiment is removed (inactivated) if activate is 1
! otherwise read the uncertainty to be set
!      write(*,*)'3D after set_c_or_e:',ip,': ',trim(cline),new%uncertainty
! set the default uncertainty to 10% of value 
      if(new%uncertainty.gt.zero) then
         yyy=1.0D-1*abs(new%prescribed)
      else
         yyy=new%uncertainty
      endif
      kp=ip
! bug reading the value after : ??
!      write(*,*)'3D uncertanity 2: ',ip,'"'//trim(cline)//'"'
! NOTE that gparcd increments ip before seaching for value
!      write(*,*)'3D Calling gparcd: ',trim(cline),ip
      call gparcdx('Uncertainty: ',cline,ip,1,usymbol,'1.0',&
           '?Enter experiment')
!      write(*,*)'3D extracted uncertainity: ',ip,buperr,'"'//usymbol//'"'
      jc=1
      call getrel(usymbol,jc,xxx)
!      read(*,'(a)')ch1
!      if(ch1.eq.'q') stop 'wrong place ...'
!      write(*,*)'3D even more: ',xxx,buperr
      if(buperr.eq.0) then
! usymbol is a numeric value !!
         if(xxx.le.zero) then
            write(*,*)'Uncertainty must not be zero, set to 0.1 of value'
            xxx=0.1*new%prescribed
         endif
         new%symlink2=0
         new%uncertainty=abs(xxx)
      else
! we should check that the symbol is not an expression ... how?
         buperr=0
         call capson(usymbol)
         call find_svfun(usymbol,istv)
!         write(*,*)'3D uncertainty symbol: ',usymbol,istv
         if(gx%bmperr.ne.0) then
            write(*,*)'3D No such symbol: ',usymbol,&
                 ' uncertainty set to 0.1 of value'
            xxx=0.1*new%prescribed
            new%symlink2=0
            new%uncertainty=abs(xxx)
         else
! check that the symbol is a constant
            if(.not.btest(svflista(istv)%status,SVCONST)) then
               write(*,*)'3D Experimental uncertainty symbol must be a value'
               gx%bmperr=4399; goto 1000
            endif
            new%symlink2=istv
         endif
      endif
! this is for relative errors, if last character is % it is a relative error!!
!      write(*,*)'3D relative errors: ',ip,len_trim(cline),'"',trim(cline),'"'
      if(ip.lt.len(cline)) then
         qp=len_trim(cline)
         if(cline(ip:ip).eq.'%' .or. cline(qp:qp).eq.'%') then
            if(new%experimenttype.eq.0) then
               new%experimenttype=100
!               write(*,*)'3D error is relative!'
            else
! the experiment is an inequality
               write(kou,*)'3D *** Inequalites must have absolute uncertainty'
!            new%experimenttype=101*new%experimenttype ???
            endif
         endif
      endif
! if weight is negative (meaning first experiment) set to unity
      if(ceq%weight.lt.zero) then
         ceq%weight=one
      endif
   endif
! any  more experiments?
!   write(*,*)'3D exp4: ',trim(cline),ip,kp,len(cline),len_trim(cline)
   if(kp.le.ip .and. len_trim(cline).gt.ip) then
!      write(*,*)'3D more experiments',trim(cline),kp,ip
      goto 17
   endif
1000 continue
 end subroutine enter_experiment

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function same_statevariable
!\begin{verbatim} %-
 logical function same_statevariable(svr1,svr2)
! returns TRUE if the state variable records are identical
   type(gtp_state_variable), pointer :: svr1,svr2
!\end{verbatim}
   logical same
   same=.FALSE.
   if(svr1%statevarid.ne.svr2%statevarid) goto 1000
   if(svr1%unit.ne.svr2%unit) goto 1000
   if(svr1%phref.ne.svr2%phref) goto 1000
   if(svr1%argtyp.ne.svr2%argtyp) goto 1000
   if(svr1%argtyp.gt.0) then
      if(svr1%phase.ne.svr2%phase) goto 1000
      if(svr1%argtyp.gt.1) then
         if(svr1%compset.ne.svr2%compset) goto 1000
         if(svr1%argtyp.gt.2) then
            if(svr1%component.ne.svr2%component) goto 1000
            if(svr1%argtyp.gt.3) then
               if(svr1%constituent.ne.svr2%constituent) goto 1000
            endif
         endif
      endif
   endif
! they are the same !!!
   same=.TRUE.
1000 continue
   same_statevariable=same
   return
 end function same_statevariable

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_condition
!\begin{verbatim}
 subroutine set_condition(cline,ip,ceq)
! to set a condition
   implicit none
   character cline*(*)
   integer ip
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
! New is set to the new condition or experiment
! in this subroutine new is not used for anything.
! in enter_experiment the new variable is removed from the condition list
! and instead added to the experimenal list
   type(gtp_condition), pointer :: new
!   write(*,*)'3D set_cond: ',cline(1:len_trim(cline)),ip
   call set_cond_or_exp(cline,ip,new,0,ceq)
1000 continue
! always mark that current equilibrium may not be consistent with conditions
   ceq%status=ibset(ceq%status,EQINCON)
   nullify(new)
 end subroutine set_condition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_cond_or_exp
!\begin{verbatim} %-
 subroutine set_cond_or_exp(cline,ip,new,notcond,ceq)
! decode an equilibrium condition, can be an expression with + and -
! the expression should be terminated with an = or value supplied on next line
! like "T=1000", "x(liq,s)-x(pyrrh,s)=0", "mu(cr)-1.5*mu(o)=muval"
! Illegal with number before first state variable !!!
! It can also be a "NOFIX=<phase>" or "FIX=<phase> value"
! The routine should also accept conditions identified with the "<number>:"
! where <number> is that preceeding each condition in a list_condition
! It should also accept changing conditions by <number>:=new_value
! The pointer to the (most recent) condition or experiment is returned in new
! notcond is 0 if a condition should be created, otherwise an experiment
   implicit none
   integer ip,notcond
   character cline*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_condition), pointer :: new
!\end{verbatim} %+
   integer nterm,kolon,iqz,krp,jp,istv,iref,iunit,jstv,jref,junit,jl,ks
   integer linkix,norem,ics,kstv,iph,nidfirst,nidlast,nidpre,qp,firstc,lpos
! a long line with conditions can create overflow and lost values ...
   character stvexp*500,stvrest*500,textval*32,c5*5,ch1
   character svtext*500,encoded*60,defval*18,actual_arg*24,svfuname*16
   integer indices(4),allterms(4,10),seqz,experimenttype
   integer ich,back,condvalsym,symsym,nextexp,colon
   double precision coeffs(10),xxx,value,ccc
   logical inactivate
! memory leak
   type(gtp_state_variable), target :: svrvar
   type(gtp_state_variable), pointer :: svr,svr2
   type(gtp_state_variable), dimension(10), target :: svrarr
   TYPE(gtp_condition), pointer :: temp
! safeguard: call old messy routine
!   call set_cond_or_exp_old(cline,ip,new,notcond,ceq)
!   return
!=========================================================================
   if(len_trim(cline).gt.400) then
      write(*,*)'3D *** Too long line with conditions:',len_trim(cline)
      gx%bmperr=4399; goto 1000
   endif
!
   nullify(temp)
   xxx=zero
   symsym=0
   iunit=0
   iref=0
   actual_arg=' '
!   write(*,*)'3D set cond or enter exper: ',trim(cline),ip
! return here to deconde another condition on the same line
50 continue
   nterm=0
   allterms=0
!==========================================================================
! return here to decode anther state variable term for condition
! step 1 extract the state variable termintade by + - = > < or :=
! NOT SUFFICIENT, a constituent can have a + or - !!!
55 continue
   experimenttype=0
   nullify(new)
   if(nterm.eq.0) then
! for second and later term coeffs already set below after call to termterm
      coeffs(1)=one
   endif
   indices=0
! the list of experiments changes ???
! NOTE we can have several conditions on the same line!!
! argument 4 equal to 5 of gpar* means extract the whole line
   stvexp=' '
   nextexp=ip
!   write(*,56)'3D scoe: ',nterm,ip,trim(cline)
56 format(a,2i3,' "',a,'" ')
   if(nterm.eq.0) then
! the whole line is read into stvexp, ip is increemented by 1
      call gparcdx('State variable: ',cline,ip,5,stvexp,'T','?Set condition')
   else
! the whole expression must have been entered on the same line
! note cline is updated below !!
      ip=ip-1
      call gparcdx(' ',cline,ip,5,stvexp,'!','?Set condition')
   endif
   if(stvexp(1:1).eq.' ') then
! if an expression is terminated with an empty line ask for value
      if(nterm.gt.0) goto 67
! if no terms and the line empty return error code for no condition
      gx%bmperr=4126; goto 1000
   elseif(stvexp(1:1).eq.'!') then
! this is an error while continuing reading an expression
      gx%bmperr=4126; goto 1000
   elseif(stvexp(1:3).eq.'FIX') then
! special case when called internally for setting phase fix
      inactivate=.FALSE.
      ip=5
      goto 299
   elseif(stvexp(1:5).eq.'NOFIX') then
      inactivate=.TRUE.
!      write(*,*)'3D Inactivate phase fix condition'
      ip=7
      goto 299
   endif
! this can be a condition or experiment ... and have several terms
! check for +, -, =, <, >, or :=  
! previous value of ip irrelevant, 
! ip points at terminator inside stvexp for current state variable
! if lpos>0 is where to start for next term
   call termterm(stvexp,ich,ip,lpos,ccc)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,48)'3D tt: ',ich,ip,lpos,trim(stvexp),stvexp(1:ip),ccc
48 format(a,3i4,' "',a,'" >',a,'< ',1pe12.4)
   if(ich.eq.6) then
! special case when condition number provided, extract the number, can be *
! meaning "all conditions", for example *:=NONE
      if(notcond.ne.0) then
!         write(*,*)'Experiments have no number'
         gx%bmperr=4131; goto 1000
      endif
      qp=1
      if(stvexp(qp:qp).eq.'*') then
!         write(*,*)'3D Special case of deleting all conditions'
! 0 means only conditions deleted, not the equilibrium
         call delete_all_conditions(0,ceq)
         goto 1000
      endif
      call getrel(stvexp,qp,xxx)
      if(buperr.ne.0) then
!         write(*,*)'No such condition number'
         gx%bmperr=4131; goto 1000
      endif
! the condition number must be an integer
      qp=-int(xxx)
! search for condition with number -qp
!      write(*,*)'3D looking for condition: ',-qp
! UNFINISHED: one should look for the qp:th ACTIVE condition ....
      temp=>ceq%lastcondition
!      write(*,*)'3D calling get_condition A'
      call get_condition(qp,svr,temp)
!      write(*,*)'3D Back from calling get_condition A'
      if(gx%bmperr.ne.0) goto 1000
! save link to old condition in new
      new=>temp
      xxx=new%prescribed
      nterm=1
!      write(*,*)'Found condition',-qp,xxx
! jump from here to 67 if condition specified as number:=value
      goto 67
   elseif(ich.lt.3 .and. nterm.eq.0) then
! first term of state variable which is an expression,
! this term is terminated by + or -
!      if(stvexp(1:1).ge.'0' .and. stvexp(1:1).le.'9') then
!      write(*,*)'3D extract coeff for first term, if none set to 1',ich
      firstc=1
      call getrel(stvexp,firstc,coeffs(1))
      if(buperr.ne.0) then
!         write(*,*)'3D error in coefficient for first condition term',buperr
! ignore error, means no coefficient
         buperr=0
         coeffs(1)=one
      else
! character after coefficient must be a *
         if(stvexp(firstc:firstc).ne.'*') then
            write(*,*)'3D coefficient is not terminated by *: ',&
                 stvexp(firstc:firstc)
            gx%bmperr=4130
         else
! update stvexp for next term, we must also update lpos ... lousy coding
            stvexp=stvexp(firstc+1:)
            lpos=lpos-firstc
         endif
      endif
   endif
!---------------------------------
! check it is a legal state variable, ignore terminator 
   svtext=stvexp(1:ip-1)
   symsym=0
!   write(*,*)'3D calling decode for "',trim(svtext),'" upt to:',ip-1
! memory leak
   svr=>svrvar
   call decode_state_variable(svtext,svr,ceq)
!   write(*,*)'3D state var: ',trim(svtext),gx%bmperr
   if(gx%bmperr.ne.0) then
! Experiments can be symbols
!      write(*,*)'3D not a state variable: ',svtext(1:5),gx%bmperr,notcond
      if(notcond.ne.0) then
         gx%bmperr=0
         svfuname=svtext
         call capson(svfuname)
!         write(*,*)'3D Searching for symbol: ',svfuname
!         call find_svfun(svfuname,symsym,ceq)
         call find_svfun(svfuname,symsym)
         if(gx%bmperr.ne.0) then
            write(*,*)'3D Experimental symbol neither state variable nor symbol'
            goto 1000
         endif
!         write(*,*)'3D experiment is a symbol ',symsym
! we do not have a state variable ... bypass some checks
         nullify(svr)
         goto 77
      else
         goto 1000
      endif
   endif
! convert to old state variable format
!   write(*,12)svr%argtyp,svr%phase,svr%compset,svr%component,svr%constituent
12 format('3D Decoded: ',5i5)
   indices=0
   if(svr%argtyp.eq.1) then
      indices(1)=svr%component
   elseif(svr%argtyp.eq.2) then
      indices(1)=svr%phase
      indices(2)=svr%compset
   elseif(svr%argtyp.eq.3) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%component
   elseif(svr%argtyp.eq.4) then
      indices(1)=svr%phase
      indices(2)=svr%compset
      indices(3)=svr%constituent
   endif
!   write(*,12)svr%argtyp,indices
!   do ks=1,4
!      allterms(nterm,ks)=indices(ks)
!   enddo
!   write(*,*)'3D newcond: ',svtext(1:20),ip
!----------------------------------------------------
77 continue
!   write(*,*)'3D error search ',notcond,ich,associated(temp)
! check that we we have a legal state variable for conditions
   if(notcond.eq.0) then
! it is a condition, check if allowed as condition
      istv=svr%oldstv
      if(istv.lt.0) then
! this means a symbol like TC or BMAGN, not allowed
         gx%bmperr=4127; goto 1000
      endif
      kstv=(svr%oldstv+1)/10+5
      if(kstv.eq.14 .or. kstv.eq.15) then
! this means state variables Q or DGM which cannot be used as condition
         gx%bmperr=4127; goto 1000
      endif
      if(istv.ge.3 .and. istv.le.5) then
! this is MU, AC and LNAC, do not allow with phase index (at present at least)
         if(indices(2).ne.0) then
            write(*,*)'Phase specific chemical potentials not allowed',&
                 ' as conditions'
            gx%bmperr=4127; goto 1000
         endif
      endif
! state variables with a single term will be prompted with current value
      encoded=' '
      call get_state_var_value(svtext,xxx,encoded,ceq)
      if(gx%bmperr.ne.0) then
! This error occurs when setting the first compositions before any calculations
         gx%bmperr=0; xxx=zero
      endif
   else
! It is an experiment, search for old value if experiment already entered
!      write(*,*)'3D experiment 1',notcond,ich
      if(ich.eq.4) then
         experimenttype=-1
      elseif(ich.eq.5) then
         experimenttype=1
      endif
      temp=>ceq%lastexperiment
      if(.not.associated(temp)) then
         xxx=zero
      else
! new is nullified, temp set above for search of conditions or experiments, 
88       continue
         temp=>temp%next
         if(symsym.eq.0) then
! new experiment is a state variable, what about temp?
!            write(*,*)'3D oldexp 1: ',symsym,temp%statev
            if(temp%statev.eq.0) then
               svr2=>temp%statvar(1)
               if(same_statevariable(svr,svr2) .and. &
                    experimenttype.eq.temp%experimenttype) then
                  xxx=temp%prescribed
! found experimental record, save link in new
                  new=>temp
               endif
            endif
            if(.not.associated(temp,ceq%lastexperiment)) goto 88
         else
! experiment is a symbol compare with other experiments for symbols
!            write(*,*)'3D oldexp 2: ',symsym,temp%statev
            if(symsym.eq.temp%statev .and. &
                 experimenttype.eq.temp%experimenttype) then
               xxx=temp%prescribed
! we have found a record for this experiment
               new=>temp
            else
               if(.not.associated(temp,ceq%lastexperiment)) goto 88
            endif
         endif
      endif
   endif
!   write(*,*)'3D Found old condition or experiment?',notcond
   if(notcond.eq.0) then
!----------------------------------------------------------------
! Only for conditions: save current term if several
      nterm=nterm+1
!      write(*,*)'3D several terms: ',nterm
      svrarr(nterm)=svr
!      write(*,*)'3D segfault search 3',nterm
! convert to old format, currently we need to store both formats ....
!   istv=svr%oldstv
      iref=svr%phref
      iunit=svr%unit
      do ks=1,4
         allterms(ks,nterm)=indices(ks)
      enddo
!      write(*,*)'3D segfault search 4',nterm,ich
!   write(*,*)'3D old indices:',nterm,indices
      if(ich.eq.1 .or. ich.eq.2) then
! terminator + or - means state variable expression with several terms
         if(nterm.gt.1) then
! UNFINISHED check the second or later state variable of same type as first
            continue
         endif
! multiterm expression, jump back to 55
!         write(*,*)'3D problems entering expression: ',trim(stvexp),lpos
         coeffs(nterm+1)=ccc
!         cline=stvexp(ip-1:); ip=1
         cline=stvexp(lpos:); ip=1
         goto 55
      endif
   else
! it is an existing experiment, we have only one term for experiments
      nterm=1
!      write(*,*)'3D segfault search 2C',nterm,associated(svr),notcond
      if(associated(svr)) then
         svrarr(nterm)=svr
!      else
!         write(*,*)'3D svr not associated, experiment is a symbol'
      endif
   endif
! jump here for qp:= or if several terms are terminated with empty line
67 continue
!==================================================================
! Step 2 ask for the numerical value or symbol, first insert a default value
   jp=1
   defval=' '
!   write(*,*)'3D default value: ',xxx,ip,' "',stvexp(ip:ip+5),'" '
   call wrinum(defval,jp,10,0,xxx)
   if(buperr.ne.0) then
      buperr=0; defval=' '
   endif
!157 continue
! stvexp is the whole line after the command
   colon=index(stvexp,':')
!   colon=index(cline,':')
!   write(*,*)'3D value: ',ip,' "',trim(stvexp),'" ',defval,colon
   call gparcdx('Value: ',stvexp,ip,1,textval,defval,'?Set condition')
!   write(*,*)'3D value: ',textval
   c5=textval(1:5)
   call capson(c5)
   none: if(c5.eq.'NONE ' .or. c5.eq.'<NONE' .or. c5.eq.'>NONE') then
      inactivate=.true.
!      write(*,158)'Inactivate condition: ',-qp,value,xxx
158   format(a,i5,2(1pe12.4))
      value=xxx
   else
      if(notcond.ne.0) then
         if(textval(1:1).eq.'<') then
            experimenttype=-1
            textval(1:1)=' '
         elseif(textval(1:1).eq.'<') then
            experimenttype=-1
            textval(1:1)=' '
         endif
      endif
      linkix=-1
      inactivate=.FALSE.
      jp=1
      call getrel(textval,jp,value)
      if(buperr.ne.0) then
! it can be a symbol
         buperr=0
         svfuname=textval
         call capson(svfuname)
!         call find_svfun(svfuname,condvalsym,ceq)
         call find_svfun(svfuname,condvalsym)
!         write(*,*)'3D Symbol link: ',textval(1:10),condvalsym,gx%bmperr
         if(gx%bmperr.ne.0) then
            write(*,*)'Condition value must be numeric or a symbol'; goto 1000
         endif
! only allowed if symbol is constant SVCONST or SVFVAL set 
! check we actually have correct symbol!!
!         write(*,*)'Symbol name: ',svflista(condvalsym)%name
         if(btest(svflista(condvalsym)%status,SVCONST) .or. &
              btest(svflista(condvalsym)%status,SVFVAL)) then
            linkix=condvalsym
            value=evaluate_svfun_old(linkix,actual_arg,1,ceq)
         else
            write(*,*)'3D Symbol must be constant or "evaluate explicit"'
            gx%bmperr=4293; goto 1000
         endif
      endif
! we must update ip in cline for uncertainty and another experiment
      if(colon.gt.0) then
         ip=colon
         istv=0
!         write(*,*)'3D changed experiment 1: ',ip,'"',trim(stvexp),'"',value 
      endif
   endif none
!
   findrecord: if(notcond.eq.0) then
! remove a condition
!      write(*,*)'3D avoiding creating expression:',associated(new)
      if(.not.associated(new)) then
! search if condition already exist
!         write(*,*)'3D searching for condition or experiment?'
         temp=>ceq%lastcondition
         if(nterm.eq.1) then
            call get_condition(nterm,svr,temp)
         else
            call get_condition_expression(nterm,svrarr,temp)
         endif
         if(gx%bmperr.ne.0 .and. inactivate) then
            write(kou,140)
140         format('Attempt to remove a non-existing condition')
            goto 1000
         endif
! the error code it will be tested below to create a condition record
      endif
   else
! remove or change an experiment
      if(.not.associated(new)) then
!         write(*,*)'3D First experiment: ',associated(temp)
! search for an experiment with state variable svr or symbol symsym
         temp=>ceq%lastexperiment
142      continue
!         write(*,*)'3D searching for experiment'
         if(symsym.eq.0) then
            istv=symsym
!            if(associated(temp)) then
!               write(*,*)'3D Calling get_condition C',allocated(temp%condcoeff)
!            else
!               write(*,*)'3D Calling get_condition C with temp null'
!            endif
            call get_condition(nterm,svr,temp)
!            write(*,*)'3D Back from calling get_condition C',gx%bmperr
            if(gx%bmperr.eq.0) then
! we must also test eperimenttype, if not same continue search
               if(temp%experimenttype.eq.experimenttype) goto 142
            else
! ensure temp is OK
               temp=>ceq%lastexperiment
!               write(*,*)'3D We are here',associated(temp),gx%bmperr
            endif
         else
!            write(*,*)'3D searching for experiment with symbol'
! temp is changed inside !!
            new=>temp
!            call get_experiment_with_symbol(symsym,experimenttype,temp)
            call get_experiment_with_symbol(symsym,experimenttype,new)
         endif
!      else
!         write(*,*)'We have a condition record in new', gx%bmperr
      endif
   endif findrecord
!======================================================
! step 3 create condition or experiment record, jump here from fix phase
199 continue
!   write(*,*)'3D at 199: ',associated(new),associated(temp),gx%bmperr,ip
   createrecord: if(gx%bmperr.eq.0) then
! no error code means we have found the condition/experiment
      if(.not.associated(new)) then
! the existing condition/experiment is in temp
         new=>temp
      endif
      if(inactivate) then
         new%active=1
!         write(*,*)'Inactivating condition',new%prescribed,new%active
      else
!         write(*,*)'3D looking for missing % ...'
! set the new value in the old condition/experiment remove any previous link!!
! linkix is link to symbol representing the value
         new%active=0
         new%prescribed=value
         new%symlink1=linkix
!         write(*,*)'3D Changing value of condition',istv,linkix,value
! special if istv=1 or 2 as ceq%tpfun should be updated
         if(istv.eq.1) then
! Save new T also locally in ceq
!            write(*,*)'3D we are here 1',ceq%tpval(1)
!            if(linkix.gt.0) then
!               write(*,*)'Cannot handle symbol as T value'
! it is allowed now
!               gx%bmperr=4293; goto 1000
!            endif
            ceq%tpval(1)=value
         elseif(istv.eq.2) then
!            if(linkix.gt.0) then
!               write(*,*)'Cannot handle symbol as P value'
! it is allowed now
!               gx%bmperr=4293; goto 1000
!            endif
            ceq%tpval(2)=value
         endif
! the uncertainty for experiments will be asked for later
! To avoid that valgrind compains uncertainty is not initiallized ...
!         write(*,*)'3D initiallizing uncertainty 2'
!         write(*,*)'3D changed experiment: ',ip,'"',trim(cline),'"',value
         new%uncertainty=zero
      endif
   else
! If we have an error from findrecord then create condition/experiment record
!      write(*,113)inactivate,gx%bmperr,nterm,istv,iunit,iref,linkix,&
!           allterms(1,1),value
113   format('Creating condition record',l2,2x,7i4,1pe12.4)
      gx%bmperr=0
      if(notcond.eq.0) then
         if(associated(ceq%lastcondition)) then
            seqz=ceq%lastcondition%seqz+1
         else
            seqz=1
         endif
         temp=>ceq%lastcondition
         allocate(ceq%lastcondition)
         new=>ceq%lastcondition
!         write(*,*)'3D new condition ',istv,symsym
      else
! it is an experiment ... ip OK here
!        write(*,13)'3D Creating new experiment record: ',symsym,value,linkix,ip
13       format(a,i5,1pe12.4,5i5)
         if(associated(ceq%lastexperiment)) then
            seqz=ceq%lastexperiment%seqz+1
         else
            seqz=1
         endif
         temp=>ceq%lastexperiment
         allocate(ceq%lastexperiment)
         new=>ceq%lastexperiment
!         write(*,*)'3D new experiment 2',istv,symsym,seqz
         istv=symsym
      endif
!      write(*,*)'3D we are here 3'
! for new conditions and experiments
      new%noofterms=nterm
      new%statev=istv
      new%iunit=iunit
      new%iref=iref
      new%active=0
      new%seqz=seqz
! To avoid that valgrind compains uncertainty is not initiallized ...
!      write(*,*)'3D initiallizing uncertainty 3',value
      new%uncertainty=zero
!      write(*,*)'3D experimenttype: ',experimenttype,symsym,nterm
      new%experimenttype=experimenttype
!      write(*,*)'3D symsym: ',symsym,nterm
      if(symsym.eq.0) then
! DO NOT allocate terms for condcoeff and indices if symbol
         allocate(new%condcoeff(nterm))
         allocate(new%indices(4,nterm))
!   write(*,*)'3D allocations ok',linkix,value
         do jl=1,nterm
            new%condcoeff(jl)=coeffs(jl)
!         write(*,111)'3D allterms:  ',istv,jl,(allterms(ks,jl),ks=1,4)
            do ks=1,4
               new%indices(ks,jl)=allterms(ks,jl)
            enddo
!         write(*,111)'3D in record: ',istv,jl,(new%indices(ks,jl),ks=1,4)
111         format(a,i3,i5,2x,4i4)
         enddo
! Only experiments can be symbols, what to do next?
      endif
!      write(*,*)'3D storing value: ',value,linkix
      if(linkix.lt.0) then
         new%prescribed=value
         new%symlink1=-1
      else
         new%symlink1=linkix
         value=evaluate_svfun_old(linkix,actual_arg,1,ceq)
!         write(*,*)'3D evaluating condition sysmbol ',linkix,value
         if(gx%bmperr.ne.0) then
            goto 1000
         endif
         new%prescribed=value
!         write(*,*)'3D  prescribed condition value: ',new%prescribed
      endif
! first test if condition on P or T is larger than 0.1
      if(istv.eq.1 .or. istv.eq.2) then
         if(value.lt.0.1D0) then
            gx%bmperr=4187; goto 1000
         endif
      endif
! special for T and P, change the local value and mark tpres
!      write(*,*)'3D set condition/enter experiment ',istv,value
      if(istv.eq.1) then
         ceq%tpval(1)=value
! Force recalculation of all TP functions but only in current equilibrium
         ceq%eq_tpres%tpused(1)=ceq%tpval(1)+one
!         write(*,*)'3D Changing tpused: ',ceq%eq_tpres%tpused(1)
      elseif(istv.eq.2) then
         ceq%tpval(2)=value
! Force recalculation of all TP functions.  Is there a better way?
         ceq%eq_tpres%tpused(2)=ceq%tpval(2)+one
      endif
! Another way to force recalculation of all TP functions by incrementing
! an integer in the tpfuns record of all TPFUN.  Used during assessments
!   do jl=1,freetpfun-1
!      tpfuns(jl)%forcenewcalc=tpfuns(jl)%forcenewcalc+1
!   enddo
!      write(*,*)'3D allocation of statvar ',symsym,istv,nterm
      if(symsym.eq.0) then
! store the state variable record in the condition, if symbol do not allocate
         allocate(new%statvar(nterm))
         do jl=1,nterm
            new%statvar(jl)=svrarr(jl)
         enddo
!      else
! experiment is a symbol, no statvar record !!
! The index to the state variable symbol is symsym stored in new%statev
!         write(*,*)'3D experiment is a symbol 2',istv,symsym,linkix
!         allocate(new%statvar(1))
!         new%statvar(1)%statevarid=0
!         new%statvar(1)%argtyp=-symsym
      endif
! link the new record into the condition list
!      write(*,*)'3D linking condition or experiment'
      if(associated(temp)) then
!       write(*,*)'Second or later condition'
         nidlast=temp%next%nid
         nidfirst=temp%nid
         nidpre=temp%previous%nid
         new%nid=nidlast+1
         temp%next%previous=>new
         new%next=>temp%next
         temp%next=>new
         new%previous=>temp
      else
! create the circular list
         new%nid=1
         new%next=>new
         new%previous=>new
      endif
      if(notcond.ne.0) then
! STRANGE ERRORS HERE
! we are actually entering an experiment, terminate here
! if textval(jp:jp) is ":" we have to step back ip one position
! increment ip with nextexp!
         ip=ip+nextexp
!         write(*,*)'3D exit? "',trim(cline),'" "',trim(textval),'"',&
!              ip,jp,nextexp
         if(cline(ip:ip).eq.':') then
! the colon should be the character at ip to extract uncertainty by getrel
            goto 1000
         endif
200      continue
         ip=ip+1
!         write(*,*)'3D looking for ":"',trim(cline),ip
         if(cline(ip:ip).eq.':' .or. ip.gt.len(cline)) goto 1000
         goto 200
!         if(ip.lt.len_trim(cline)) goto 200
! confusion of subtexts ... move ip forward to point at colon
!            write(*,*)'3D where is :? ',cline(ip+1:ip+1),ip
!            if(cline(ip+1:ip+1).ne.':') then
!               ip=ip+1; if(ip.lt.70) goto 200
!            endif
!            ip=ip+1  very confused !!
!            ip=ip-1
!         endif
! stop here so I can check
!         read(*,'(a)')ch1
!         if(ch1.eq.'q') stop 'cannot find ":"'
! allow for more experiments on the same line ...
!         cline=textval(jp:)
!         ip=1
!         goto 1000
      endif
   endif createrecord
!   write(*,*)'3D end of createrecord',ip,'"',trim(stvexp),'"'
!----------------------------------------------------------------
! if there is more in stvexp go back to label 50 ...
   if(.not.eolch(stvexp,ip)) then
!      write(*,*)'3d first character: "',stvexp(ip:ip),'" '      
! NOTE gparc skips the first character in cline
      if(stvexp(ip:ip).eq.':') then
! if experiment there can be an uncertainty ...
!         write(*,*)'3D where is the value?'
         cline=stvexp(ip:)
         cline(ip:ip)=' '
         ip=1
         goto 1000
      elseif(stvexp(ip:ip).eq.',') then
         cline=stvexp(ip:)
      else
         cline=stvexp(ip-1:)
      endif
!      write(*,*)'3d next condition: "',stvexp(ip:ip+20),'"'
      ip=1; goto 50
   endif
   goto 1000
!====================================================================
! Special below is for fix/unfix phases
299 continue
   if(notcond.ne.0) then
      write(kou,*)'3D Illegal to set a fix phase as experiment'
      gx%bmperr=4294; goto 1000
   endif
!   write(*,*)'3D fix phase 2: ',ip,stvexp(ip:40)
   call find_phase_by_name(stvexp(ip:),iph,ics)
   if(gx%bmperr.ne.0) then
      goto 1000
   endif
!   write(*,*)'3D Phase index: ',iph,ics
   nterm=1
   istv=-iph
   iref=ics
   iunit=0
   linkix=-1
   coeffs(1)=1.0D0
   do jl=1,4
      allterms(jl,1)=0
   enddo
! convert to state variable
!   write(*,*)'3D Setting svrarr(1) values'
   svrarr(1)%statevarid=istv
   svrarr(1)%oldstv=istv
   svrarr(1)%phase=ics
   svrarr(1)%unit=0
   svrarr(1)%argtyp=0
   svrarr(1)%phase=iph
   svrarr(1)%compset=ics
   svrarr(1)%component=0
   svrarr(1)%constituent=0
!
   temp=>ceq%lastcondition    
! if not inactivate get value
   if(inactivate) then
! bypass phase name
      ip=index(stvexp,' ')
   else
      ip=index(stvexp,'==')+2
      call getrel(stvexp,ip,value)
      if(buperr.ne.0) then
         write(*,*)'3D error setting fix amount ',ip,stvexp(1:40)
         gx%bmperr=4100; goto 1000
      endif
   endif
   svr=>svrarr(1)
!   write(*,*)'3D set_cond_or_exp for fix phase: ',svr%statevarid,svr%phase
! new must be unassociated, for inactivate temp will be set to condition.
   nullify(new)
   call get_condition(nterm,svr,temp)
!   write(*,*)'3D Back from get_condition ',gx%bmperr,ip
   goto 199
!
! finally, for conditions T or P copy value to ceq%tpval
! This may be a bit inconsistent .... but??
900 continue
   if(istv.eq.1 .and. iunit.eq.0 .and. iref.eq.0) then
      ceq%tpval(1)=value
   elseif(istv.eq.2 .and. iunit.eq.0 .and. iref.eq.0) then
      ceq%tpval(2)=value
   endif
! mark that any current results may be inconsistent with new conditions
!   globaldata%status=ibset(globaldata%status,GSINCON)
   ceq%status=ibset(ceq%status,EQINCON)
1000 continue
!   write(*,*)'exit set_condition, T= ',ceq%tpval(1)
! possible memory leaks
   nullify(svr)
   nullify(svr2)
!   nullify(svrarr)
   nullify(temp)
   return
 end subroutine set_cond_or_exp !svr ip

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine get_experiment_with_symbol
!\begin{verbatim} %-
 subroutine get_experiment_with_symbol(symsym,experimenttype,temp)
! finds an experiment with s symbol index symsym and exp.type
   implicit none
   integer symsym,experimenttype
   type(gtp_condition), pointer :: temp
! NOTE: temp must have been set to ceq%lastcondition before calling this
!\end{verbatim}
! pcond: pointer, to a gtp_condition record for this equilibrium
   type(gtp_condition), pointer :: pcond,last
   if(.not.associated(temp)) goto 900
   last=>temp
   pcond=>last
100 continue
   if(pcond%statev.eq.symsym .and. pcond%experimenttype.eq.experimenttype) then
! the index of the symbol is stored in statev, we have found the experiment
      goto 1000
   else
! Wow = here instead of => created a lot of problems!!!
      pcond=>pcond%next
! this is true unless we have circulated the whole list
      if(.not.associated(pcond,last)) goto 100
   endif
! we have not found this experiment
900 continue
   gx%bmperr=4131
1000 continue
   return
 end subroutine get_experiment_with_symbol

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine get_condition_expression
!\begin{verbatim}
 subroutine get_condition_expression(nterm,svrarr,pcond)
! I do not want to change get_condition ,,,,, suck
! finds a condition/experiment record with the given state variable expression
! If nterm<0 svr is irrelevant, the absolute value of nterm is the sequential
! number of the ACTIVE conditions
   implicit none
   integer nterm
   type(gtp_state_variable), dimension(*), target :: svrarr
! NOTE: pcond must have been set to ceq%lastcondition before calling this
! pcond: pointer, to a gtp_condition record for this equilibrium
   type(gtp_condition), pointer :: pcond
!\end{verbatim} %+
   type(gtp_condition), pointer :: last
   type(gtp_state_variable), pointer :: svr,condvar
   integer jj
!
!   write(*,*)'3D get_condition_expression: ',nterm
! start from first equilibrium in circular list
!   write(*,*)'3D one more line ...'
! this write statement caused crash if first condition had 2 terms ...
!   write(*,*)'3D size of pcond%statvar: ',size(pcond%statvar)
   if(nterm.gt.1) then
      write(*,*)'3D A condition with several terms sometimes causes crash'
      gx%bmperr=4399; goto 1000
   endif
   pcond=>pcond%next
   last=>pcond
100 continue
   write(*,*)'3D at label 100'
   ploop: do while(.true.)
      terms: do jj=1,nterm
         svr=>svrarr(jj)
!         write(*,*)'3D gce: ',jj,svr%component
         if(jj.gt.size(pcond%statvar)) then
            write(*,*)'3D too many terms in condition',jj,size(pcond%statvar)
            write(*,*)'3D more:',nterm,pcond%statvar(1)%oldstv
            gx%bmperr=4399; goto 1000
         endif
         condvar=>pcond%statvar(jj)
!         write(*,*)'3D get_condition: ',jj,condvar%oldstv,condvar%argtyp
! dissapointment, one cannot compare two structures ... unless pointers same
         if(condvar%oldstv.ne.svr%oldstv) goto 200
         if(condvar%argtyp.ne.svr%argtyp) goto 200
         if(condvar%phase.ne.svr%phase) goto 200
         if(condvar%compset.ne.svr%compset) goto 200
! skip fix phases
         if(condvar%statevarid.lt.0) goto 1000
! most conditions with 2 terms are x(a)-x(b) ore similar
!         write(*,*)'3D component: ',condvar%component,svr%component
         if(condvar%component.ne.svr%component) goto 200
         if(condvar%constituent.ne.svr%constituent) goto 200
         if(condvar%norm.ne.svr%norm) goto 200
         if(condvar%unit.ne.svr%unit) goto 200
      enddo terms
! we have found a condition with these state variables
!      write(*,*)'3D Found condition',pcond%active
      goto 1000
200   continue
      pcond=>pcond%next
      if(associated(pcond,last)) exit ploop
   enddo ploop
! we did not find this condition, maybe create it?
   gx%bmperr=4131; goto 1000
1000 continue
   return
 end subroutine get_condition_expression

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine get_condition
!\begin{verbatim}
 subroutine get_condition(nterm,svr,pcond)
! finds a condition/experiment record with the given state variable expression
! If nterm<0 svr is irrelevant, the absolute value of nterm is the sequential
! number of the ACTIVE conditions
   implicit none
   integer nterm
   type(gtp_state_variable), pointer :: svr
! NOTE: pcond must have been set to ceq%lastcondition before calling this
! pcond: pointer, to a gtp_condition record for this equilibrium
   type(gtp_condition), pointer :: pcond
!\end{verbatim} %+
   type(gtp_condition), pointer :: last
   type(gtp_state_variable), pointer :: condvar
   integer j1,num,iact
   if(.not.associated(pcond)) goto 900
!   write(*,*)'3D in get_condition: ',svr%statevarid,svr%oldstv,svr%argtyp,nterm
!   if(nterm.lt.0) write(*,*)'3D Condition number: ',-nterm
!   last=>pcond
! start from first equilibrium in circular list
   pcond=>pcond%next
   last=>pcond
   num=0
   iact=0
100 continue
      num=num+1
! iact is incremented with the active conditions
      if(pcond%active.eq.0) iact=iact+1
      if(nterm.lt.0) then
! we have found the active condition with number -nterm
!         write(*,102)'Cond #: ',pcond%active,nterm,iact,num,pcond%prescribed
!102      format(a,4i3,1pe12.4)
! pcond starts with the last equilibria, not the first ...
         if(iact+nterm.eq.0) goto 1000
      elseif(.not.allocated(pcond%condcoeff)) then
! problem when experiment is a symbol ...
!      elseif(.not.allocated(pcond%condcoeff) .and. istv.ne.0) then
! no coefficients allocated, it must be an experiment with a symbol as variable
!-         write(*,*)'3D experiment as symbol',pcond%statev,pcond%seqz,num
! we must transfer the symbol index ...
!         if(pcond%statev.eq. )then
!            goto 1000
!         endif
!         goto 200
         continue
      elseif(pcond%noofterms.eq.nterm) then
! we should never be here if nterm>1
         if(nterm.gt.1) then
            write(*,*)'3D call to get_contition with nterm: ',nterm
            gx%bmperr=4399; goto 1000
         endif
!         write(*,*)'3D nterm: ',nterm,pcond%noofterms
! experiments that are symbols have not allocated any coefficent record
         do j1=1,nterm
! if nterm>1 compare just nterm as this routine called for each term!!
            condvar=>pcond%statvar(j1)
!            write(*,*)'3D get_condition: ',j1,num,condvar%oldstv,condvar%argtyp
! dissapointment, one cannot compare two structures ... unless pointers same
!            if(condvar.ne.svr) goto 200
            if(condvar%oldstv.ne.svr%oldstv) goto 200
            if(condvar%argtyp.ne.svr%argtyp) goto 200
            if(condvar%phase.ne.svr%phase) goto 200
            if(condvar%compset.ne.svr%compset) goto 200
            if(condvar%statevarid.lt.0) goto 1000
! for fix phase the remaining have no importance
!            write(*,*)'3D component: ',condvar%component,svr%component
            if(condvar%component.ne.svr%component) goto 200
            if(condvar%constituent.ne.svr%constituent) goto 200
            if(condvar%norm.ne.svr%norm) goto 200
            if(condvar%unit.ne.svr%unit) goto 200
! we must have experimenttype=0 ??
!            if(condvar%experimenttype.ne.0) goto 200
         enddo
! we have found a condition with these state variables
!         write(*,*)'3D Found condition',pcond%active
         goto 1000
!      else
!         write(*,*)'3D ignoring condition with wrong number of terms',&
!              nterm,pcond%noofterms
      endif
200   continue
!      write(*,*)'Conditions not same'
      pcond=>pcond%next
      if(.not.associated(pcond,last)) goto 100
900 continue
!   write(*,*)'3D get_condition: No such condition or experiment'
   gx%bmperr=4131; goto 1000
1000 continue
 return
end subroutine get_condition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_condition2
!\begin{verbatim} %-
 subroutine get_condition2(nterm,coeffs,istv,indices,iref,iunit,pcond)
! finds a condition record with the given state variable expression
! nterm: integer, number of terms in the condition expression
! istv: integer, state variable used in the condition
! indices: 2D integer array, state variable indices used in the condition
! iref: integer, reference state of the condition (if applicable)
! iunit: integer, unit of the condition value
! NOTE: pcond must have been set to ceq%lastcond before calling this routine!!!
! pcond: pointer, to a gtp_condition record for this equilibrium
! NOTE: conditions like expressions x(mg)-2*x(si)=0 not implemeneted
! fix phases as conditions have negative condition variable
   implicit none
   TYPE(gtp_condition), pointer :: pcond
   integer, dimension(4,*) :: indices
   integer nterm,istv,iref,iunit
   double precision coeffs(*)
!\end{verbatim} %+
   TYPE(gtp_condition), pointer :: current,first
!   integer, dimension(4) :: indx
   integer ncc,nac,j1,j2
!   write(*,*)'looking for condition'
! pcond must have been set to ceq%lastcond before calling this routine!!!
   if(.not.associated(pcond)) goto 900
   first=>pcond%next
   current=>first
!   write(*,*)'get_condition start: ',current%statev,current%active
   ncc=1
   nac=0
   if(ocv()) write(*,98)'new:',0,nterm,istv,(indices(j1,1),j1=1,4),iref,iunit
98 format(a,2x,i2,5x,2i4,5x,4i4,5x,2i3)
100 continue
   if(ocv()) write(*,98)'old:' ,current%nid,current%noofterms,current%statev,&
        (current%indices(j1,1),j1=1,4),current%iref,current%iunit
   if(nterm.eq.0) then
! why nterm=0?  Check!!!
      if(ocv()) write(*,*)'3D get_condition: ',istv,ncc,nac
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
   if(ocv()) write(*,103)'Checking terms, istv, iref and unit ',&
        nac,ncc,nterm,current%noofterms
103 format(a,6i5)
   if(current%noofterms.ne.nterm .or. current%statev.ne.istv .or. &
        current%iref.ne.iref .or. current%iunit.ne.iunit) goto 200
   if(ocv()) write(*,*)'Checking indices'
   do j1=1,nterm
      do j2=1,4
         if(current%indices(j2,j1).ne.indices(j2,j1)) goto 200
      enddo
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
 end subroutine get_condition2

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine extract_stvr_of_condition
!\begin{verbatim} %-
 subroutine extract_stvr_of_condition(pcond,nterm)
! finds a condition record with the given state variable record
! returns it as a state variable record !!!
! nterm: integer, number of terms in the condition expression
! pcond: pointer, to a gtp_condition record
   implicit none
   TYPE(gtp_condition), pointer :: pcond
   integer nterm
!\end{verbatim}
   TYPE(gtp_condition), pointer :: current,first
!   integer, dimension(4) :: indx
   integer ncc,nac,j1,istv,iref,iunit
!
!   write(*,*)'not implemented!!'
   gx%bmperr=4078; goto 1000
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

!\addtotable subroutine locate_condition
!\begin{verbatim}
 subroutine locate_condition(seqz,pcond,ceq)
! locate a condition using a sequential number
   implicit none
   integer seqz
   type(gtp_condition), pointer :: pcond
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ij
!   write(*,*)'In locate_condition 1'
   pcond=>ceq%lastcondition
!   write(*,*)'In locate_condition 2',seqz
   do ij=1,seqz
      pcond=>pcond%next
! segmentation faults in this routine when locating ceq saved during step/map
!      write(*,*)'In locate_condition 3',ij
      if(seqz.gt.ij .and. associated(pcond,ceq%lastcondition)) then
!         write(*,*)'Locate condition called with illegal index: ',seqz
         gx%bmperr=4295; goto 1000
      endif
   enddo
1000 continue
   return
 end subroutine locate_condition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine apply_condition_value
!\begin{verbatim}
 subroutine apply_condition_value(current,what,value,cmix,ccf,ceq)
! This is called when calculating an equilibrium.
! It returns a condition at each call, at first call current must be nullified?
! When all conditions done the current is nullified again
! If what=-1 then return degrees of freedoms and maybe something more
! what=0 means calculate current values of conditions
! calculate the value of a condition, used in minimizing G
! ccf are the coefficients for conditions with several terms
   implicit none
   integer what,cmix(*)
   double precision value,ccf(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq    
   TYPE(gtp_condition), pointer :: current
!\end{verbatim} %+
! ceq is actually redundant as current is a pointer to condition list in ceq
   integer, dimension(4) :: indices
   integer iref,iunit,jx,istv,ip,linkix,nterms
   character encoded*60,actual_arg*60
   double precision xxx
!
100 continue
   if(current%active.ne.0) then 
! return 0 for inactive conditions
      cmix(1)=0;  goto 1000
   endif
   if(what.ge.0) goto 200
!----------------------------------------------------------
! Here we should return information about conditions on potentials (T, P, MU)
! and fix phases
   cmix(1)=0
   if(current%noofterms.gt.1) then
      if(current%statev.eq.111 .or.current%statev.eq.110) then
! allow several terms for mole fractions and \sum_i a_i*N(i)=0 !!
!         write(*,69)'3D in apply: ',current%statev,current%noofterms,&
!              ((current%indices(jx,nterms),jx=1,4),nterms=1,current%noofterms)
69       format(a,i4,i2,3(2x,4i5))
         nterms=current%noofterms
         do jx=1,nterms
            ccf(jx)=current%condcoeff(jx)
         enddo
!         write(*,68)nterms,(ccf(jx),jx=1,nterms)
!68       format('3D coeff: ',i2,6(1pe12.4))
! VERY CLUMSY but maybe good for the moment
      elseif(current%statev.eq.21 .or. &
           current%statev.eq.130) then
! implement S-S for EET calculations ...
! and y-y for constitutions ....
!         write(*,*)'3D Apply_condition with several terms 1',current%statev,&
!              current%iref,current%iunit
         nterms=current%noofterms
         do jx=1,current%noofterms
            ccf(jx)=current%condcoeff(jx)
!            write(*,*)'3D: ',jx,(current%indices(istv,jx),istv=1,4)
         enddo
!         gx%bmperr=4207; goto 900
      else
! cannot handle other conditions with several terms
         write(*,*)'3D Apply_condition with several terms 2',current%statev,&
              current%noofterms
         gx%bmperr=4207; goto 900
      endif
   else
! one term with coefficient one
      ccf(1)=one
      nterms=1
   endif
! for debugging
   istv=current%statev
   do jx=1,4
      indices(jx)=current%indices(jx,1)
   enddo
   iref=current%iref
   iunit=current%iunit
   ip=1
   encoded=' '
   actual_arg=' '
! fetch value of symbol link if any
   if(current%symlink1.gt.0) then
      linkix=current%symlink1
!      if(btest(svflista(linkix)%status,SVFVAL)) then
      if(btest(svflista(linkix)%status,SVCONST)) then
         xxx=svflista(linkix)%linkpnode%value
! wrong         xxx=svflista(linkix)%svfv
         ceq%svfunres(linkix)=xxx
!         write(*,*)'3D SVFVAL apply: ',linkix,xxx
      else
         actual_arg=' '
! no pointer to equilibrim record ... use firsteq ??
! NOTE if symbol has bit SVCONST set then do not evaluate, use
         xxx=evaluate_svfun_old(linkix,actual_arg,1,ceq)
         if(gx%bmperr.ne.0) then
            write(*,*)'3D error evaluate symbolic link as condition',linkix,xxx
            goto 1000
         endif
      endif
      current%prescribed=xxx
   endif
!------------------
   if(current%statev.lt.0) then
! a FIX PHASE condition has state variable equal to -iph, ics is stored in iref
      cmix(1)=4
      cmix(2)=-current%statev
      cmix(3)=current%iref
      value=current%prescribed
!      write(*,*)'3D Fix phase: ',-current%statev,current%iref,value
   elseif(current%statev.eq.1) then
! temperature
      cmix(1)=1
      value=current%prescribed
!      write(*,*)'3D conditon on T'
   elseif(current%statev.eq.2) then
! pressure
      cmix(1)=2
      value=current%prescribed
!      write(*,*)'3D conditon on P'
   elseif(current%statev.le.5) then
! potentials has statev=1..5 (T, P, MU, AC, LNAC)
      cmix(1)=3
      cmix(2)=current%statev
      cmix(3)=current%indices(1,1)
      value=current%prescribed
!      write(*,*)'3D condition on MU/AC/LNAC'
   elseif(current%statev.ge.10) then
! other condition must be on (normalized) extensive properties (N, X, H etc)
      cmix(1)=5
!      write(*,*)'3D Extensive condition: ',current%statev
! SPECIAL FOR CONDITIONS ON Y to inhibit grid minimizer
      if(current%statev.eq.130) cmix(1)=6
   else
!      write(*,*)'3D Illegal condition',current%statev
      gx%bmperr=4208; goto 1000
   endif
   goto 900
!--------------------------------------
! Here we should return extensive condition, maybe calculate value
200 if(what.ne.0) goto 300
   cmix(1)=0
! for debugging
   istv=current%statev
   do jx=1,4
      indices(jx)=current%indices(jx,1)
   enddo
   iref=current%iref
   iunit=current%iunit
   ip=1
   encoded=' '
   actual_arg=' '
!------------------
   if(current%statev.lt.10) goto 900
! condition must be on extensive properties (N, X, B, W, H etc)
   cmix(1)=5
   cmix(2)=current%statev
! indices are dimensioned (4,nterms)
   cmix(3)=current%indices(1,1)
   cmix(4)=current%indices(2,1)
   cmix(5)=current%indices(3,1)
   cmix(6)=current%indices(4,1)
! for one term set coefficient to one
   ccf(1)=one
! more than one term ... this is very clumy ...
   if(current%noofterms.gt.1) then
      do nterms=2,current%noofterms
! 7, 8, 9 10 for second term, 11, 12, 13 14 for third etc
         cmix(4*nterms-1)=current%indices(1,nterms)
         cmix(4*nterms)=current%indices(2,nterms)
         cmix(4*nterms+1)=current%indices(3,nterms)
         cmix(4*nterms+2)=current%indices(4,nterms)
      enddo
      ip=current%noofterms
      do jx=1,ip
         ccf(jx)=current%condcoeff(jx)
      enddo
!      write(*,211)'3D Many terms: ',(cmix(jx),jx=1,4*ip+2)
!      write(*,212)'3D more: ',(ccf(jx),jx=1,ip)
211   format(a,2i4,4(2x,44i3))
212   format(a,4(1pe12.4))
   endif
!   if(current%noofterms.gt.2) then
!      write(*,*)'3D Apply_condition with more than 2 terms',current%noofterms
!      gx%bmperr=4207; goto 1000
!   endif
   value=current%prescribed
   if(iunit.eq.100) then
! Prescribed value is in percent, divide value by 100
      value=1.0D-2*value
!      write(*,*)'3D iunit: ',iunit,value
   endif
   goto 900
!--------------------------------------
! this part is redundant ....
300   continue
!   write(*,*)'Calling apply_condition with illegal option'
   gx%bmperr=4296; goto 1000
!-----------------------------------------------------------
! maybe something common
900 continue
!
1000 continue
   return
 end subroutine apply_condition_value

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine condition_value
!\begin{verbatim} %-
 subroutine condition_value(mode,pcond,value,ceq)
! set (mode=0) or get (mode=1) a new value of a condition.  Used in mapping
   implicit none
   integer mode
   type(gtp_condition), pointer :: pcond
   type(gtp_equilibrium_data), pointer :: ceq
   double precision value
!\end{verbatim}
   if(mode.eq.0) then
! set the value
      pcond%prescribed=value
! special for T and P
      if(pcond%statev.eq.1) then
         ceq%tpval(1)=value
      elseif(pcond%statev.eq.2) then
         ceq%tpval(2)=value
      endif
   elseif(mode.eq.1) then
      value=pcond%prescribed
   else
      write(*,*)'Condition value called with illegal mode'
   endif
1000 continue
   return
 end subroutine condition_value

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine ask_default_constitution
!\begin{verbatim}
 subroutine ask_default_constitution(cline,last,iph,ics,ceq)
! set values of default constitution interactivly
! phase and composition set already given
   implicit none
   character cline*(*)
   integer last,iph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokph,lokcs,ky,ll,iy,jy,is,ip,abel,subl
   real mmyfr(maxconst)
   character quest*32,name*24,vdef*4,fdef*8
   double precision xxx,mass
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
! if PHNOCV set the composition is fixed
   if(btest(phlista(lokph)%status1,PHNOCV)) goto 1000
   write(*,10)ics
10 format('Give min or max fractions for composition set ',i2/&
        ' use < or negative value for max, > or positive for min',&
        ' or NONE for no default')
   name=' '
   ky=0
   do ll=1,phlista(lokph)%noofsubl
      if(phlista(lokph)%nooffr(ll).gt.1) then
! more than one constituent
         do iy=1,phlista(lokph)%nooffr(ll)
            ky=ky+1
!            call get_phase_constituent_name(iph,ky,name,subl)
            call get_constituent_name(iph,ky,name,mass)
            if(gx%bmperr.ne.0) then
               write(*,*)'3D default: ',iph,ky,iy
               goto 1000
            endif
            quest='Default for '//name(1:len_trim(name))//&
                 '#'//char(ichar('0')+ll)
! use current value as default if nonzero
            vdef=' '
            abel=10*abs(ceq%phase_varres(lokcs)%mmyfr(ky))
!            write(*,*)'3D abel:',ky,abel,ceq%phase_varres(lokcs)%mmyfr(ky)
            if(abel.ge.10) then
               vdef=' 1.0'
            elseif(abel.le.0) then
               vdef=' 0.1'
            else
               vdef=' 0.'//char(ichar('0')+abel)
            endif
            if(ceq%phase_varres(lokcs)%mmyfr(ky).lt.0.0) then
               vdef(1:1)='<'
            elseif(ceq%phase_varres(lokcs)%mmyfr(ky).gt.0.0) then
               vdef(1:1)='>'
            else
               vdef='NONE'
            endif
! modified for new online help
!            call gparcd(quest,cline,last,1,fdef,vdef,q1help)
            call gparcdx(quest,cline,last,1,fdef,vdef,&
                 '?Amend phase default constit')
            jy=1
            if(fdef(1:4).eq.'NONE') then
               xxx=0
               is=1
            elseif(eolch(fdef,jy)) then
               xxx=-1.0D-1
            else
               is=1
               if(fdef(jy:jy).eq.'<') then
                  is=-1
                  jy=jy+1
               elseif(fdef(jy:jy).eq.'>') then
                  jy=jy+1
               endif
!               write(*,*)'3D def1: ',fdef,jy
               call getrel(fdef,jy,xxx)
               if(buperr.ne.0) then
!                  write(*,*)'3D buperr ',buperr
                  buperr=0
               endif
               if(is.lt.0) xxx=-xxx
            endif
            if(abs(xxx).gt.one) xxx=sign(xxx,one)
!         write(*,*)'3D default: ',xxx
            mmyfr(ky)=real(xxx)
         enddo
      else
! a single constituent, we must increment ky as there may be more
         ky=ky+1
         mmyfr(ky)=1.0
      endif
   enddo
   call enter_default_constitution(iph,ics,mmyfr,ceq)
!   write(*,99)(mmyfr(jy),jy=1,ky)
99 format('3D defy: ',15(f5.1))
1000 continue
   return
 end subroutine ask_default_constitution

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine enter_default_constitution
!\begin{verbatim}
 subroutine enter_default_constitution(iph,ics,mmyfr,ceq)
! user specification of default constitution for a composition set
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer iph,ics
   real mmyfr(*)
!\end{verbatim}
   integer lokph,lokcs,jl,jk
!   write(*,*)'3D In enter_default_constitution ',iph,ics
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   jk=size(ceq%phase_varres(lokcs)%yfr)
!   write(*,909)lokph,lokcs,phlista(lokph)%tnooffr,ceq%eqno,&
!        size(ceq%phase_varres),size(ceq%phase_varres(lokcs)%mmyfr),jk
909 format('3D 2699: ',10i4)
!   write(*,46)'3D y: ',(ceq%phase_varres(lokcs)%yfr(jl),jl=1,jk)
46 format(a,10(F7.3))
   do jl=1,phlista(lokph)%tnooffr
      ceq%phase_varres(lokcs)%mmyfr(jl)=mmyfr(jl)
!      write(*,47)'3D jl: ',jl,mmyfr(jl),&
!           firsteq%phase_varres(lokcs)%mmyfr(jl),&
!           ceq%phase_varres(lokcs)%mmyfr(jl)
   enddo
47 format(a,i2,10F7.3)
! set bit indicating that this composition set has a default constitution
!   write(*,*)'3D enter_default_constitution?? ',lokcs,&
!        ceq%phase_varres(lokcs)%mmyfr(phlista(lokph)%tnooffr)
   ceq%phase_varres(lokcs)%status2=&
        ibset(ceq%phase_varres(lokcs)%status2,CSDEFCON)
1000 continue
   return
 end subroutine enter_default_constitution

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine set_input_amounts
!\begin{verbatim}
 subroutine set_input_amounts(cline,lpos,ceq)
! set amounts like n(specie)=value or b(specie)=value
! value can be negative removing amounts
! values are converted to moles and set or added to conditions
   implicit none
   integer lpos
   character cline*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   TYPE(gtp_state_variable), target :: svr1
   TYPE(gtp_state_variable), pointer :: svr
   TYPE(gtp_condition), pointer :: current,first,last
   character species*32,cval*16,statevar*4,condline*32
   integer ielno(10)
   double precision addval(maxel)
   integer k,loksp,istv,jel,ip
   double precision xval,sumstoi,xmols
! repeat reading until empty line
100 continue
   addval=zero
   call gparcx('Species and amount as N(..)= or B(...)= : ',&
        cline,lpos,1,species,' ','?Set input amounts')
   call capson(species)
   statevar=species(1:1)
   if(statevar.ne.'    ') then
      if(.not.(statevar(1:1).ne.'N' .or. statevar(1:1).ne.'B')) then
         write(*,*)'Illegal state variable for input amounts'
         goto 1000
      endif
      k=index(species,')')
      if(k.le.3) then
         write(*,*)'Species must be surrounded by ( )'
         gx%bmperr=4297; goto 1000
      endif
      cval=species(k+1:)
! this line gave compilation warning moving 32 bytes from a space of (max) 30
      species=species(3:k-1)//'  '
!      write(*,*)'3D Species: ',trim(species),' <',trim(cval),'> '
      if(index(species,',').gt.0 .or. index(species,'(').gt.0) then
         write(*,*)'Use only N(species) or B(species) in input amounts'
         goto 1000
      endif
   else
      goto 1000
   endif
   call find_species_record(species,loksp)
! not needed as we can access splista
!   call get_species_data(loksp,nspel,ielno,stoi,spmass,qsp)
   if(gx%bmperr.ne.0) goto 1000
! if user writes N(C)=2 the =2 will be in cval, if a space after = in cline
   if(cval(1:1).eq.'=') goto 300
!   goto 300
200 continue
! the user can also give values without = or with a space before =
! but no space allowed after =
!   write(*,*)'3D cline: ',trim(cline),lpos
   call gparcx('Amount: ',cline,lpos,1,cval,' ','?Set input amounts')
300 continue
   if(cval(1:1).eq.'=') cval(1:1)=' '
   ip=1
!   write(*,*)'3D cval: ',trim(cval),ip
   call getrel(cval,ip,xval)
   if(buperr.ne.0) then
      write(*,*)'Amount must be a real number'
      goto 1000
   endif
!   write(*,*)'3D xval: ',xval
! this return the internal code for N
! BUG here as svr is no longer allocated in decode_state_variable to avoid
! memory leaks
  svr=>svr1
   call decode_state_variable('N ',svr,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'Error decoding N in set_input_amounts'
      goto 1000
   endif
   istv=svr%oldstv
! if B convert to N: moles of species = input_mass/mass_of_species
! moles of element = stoiciometry_of_element/total_number_of_elements
   if(statevar(1:1).eq.'B') then
      write(kou,*)'Note: set input in mass converted to moles'
      xmols=xval/splista(loksp)%mass
   else
      xmols=xval
   endif
   sumstoi=zero
   do jel=1,splista(loksp)%noofel
      ielno(jel)=splista(loksp)%ellinks(jel)
      addval(ielno(jel))=splista(loksp)%stoichiometry(jel)*xmols
      sumstoi=sumstoi+splista(loksp)%stoichiometry(jel)
   enddo
! now create or att to existing conditions
   last=>ceq%lastcondition
   jel=1
   if(.not.associated(last)) goto 600
! return here to look for condition for another element
500 continue
!   write(*,*)'At 500',last%nid,last%next%nid   
   first=>last%next
   current=>first
! loop for all condition
510 continue
!   write(*,*)'loop: ',current%nid,current%indices(1,1),ielno(jel)
! check if this condition match amount of element jel
   if(current%noofterms.eq.1) then
      if(current%statev.eq.istv) then
         if(current%indices(1,1).eq.ielno(jel) .and. &
              current%indices(2,1).eq.0) then
! we have found an identical contition, add the new amount
! if condition not active (active=/=0) then activate and zero prescibed amount
            if(current%active.ne.0) then
               current%active=0
               current%prescribed=zero
            endif
            current%prescribed=current%prescribed+addval(ielno(jel))
            goto 700
         endif
      endif
   endif
   current=>current%next
!   write(*,*)'next: ',current%nid,first%nid
   if(.not.associated(current,first)) goto 510
600 continue
! new condition needed
   condline='N('//ellista(ielno(jel))%symbol&
      (1:len_trim(ellista(ielno(jel))%symbol))//')='
   ip=len_trim(condline)+1
   call wrinum(condline,ip,10,0,addval(ielno(jel)))
!   write(*,*)'new condition: ',condline
! set_condition starts by incementing ip
   ip=0
   call set_condition(condline,ip,ceq)
   if(gx%bmperr.ne.0) goto 1000
   if(.not.associated(last)) then
!  if ceq%lastcondition was not associated above the call to set_condition
!  will have set link in ceq%lastcondition
      last=>ceq%lastcondition
!      write(*,*)'condition id: ',last%nid
   endif
700 continue
   jel=jel+1
   if(jel.le.splista(loksp)%noofel) goto 500
! all elements for this species set as conditions, check if any more
   if(.not.eolch(cline,lpos)) goto 100
!
1000 continue
! possible memory leaks.  Maybe also current, last
   nullify(svr)
   return
 end subroutine set_input_amounts

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine get_parameter_typty
!\begin{verbatim}
 subroutine get_parameter_typty(name1,lokph,typty,fractyp)
! interpret parameter identifiers like MQ&C#2 in MQ&C#2(FCC_A1,FE:C) ...
! find the property associated with this symbol
   implicit none
   integer typty,fractyp,lokph
   character name1*(*)
!\end{verbatim}
   integer nr,typty1,iel,isp,kel,loksp,lk3,kq,k4,kk,ll
   character elnam*24
! It can be a mobility with a & inside
   kel=index(name1,'&')
   if(kel.gt.0) then
! note that elnam may contain sublattice specification like Fe+2#2
      elnam=name1(kel+1:)
      name1=name1(1:kel-1)
   endif
   kq=len_trim(name1)
!   write(*,*)'3D: fractyp: ',kq,name1(1:kq)
   if(name1(kq:kq).eq.'D') then
! A final "D" on the paramer symbol indicates fractyp=2
      name1(kq:kq)=' '
      fractyp=2
   else
      fractyp=1
   endif
!----------------------
!      write(*,*)'Property symbol: "',propid(nr)%symbol,'" >',name1(1:4),'<'
   do nr=1,ndefprop
      if(name1(1:4).eq.propid(nr)%symbol) then
         goto 70
      endif
   enddo
! no matching property identifier
   gx%bmperr=4292; goto 1000
!
70 continue
   typty=nr
   typty1=nr
   iel=0; isp=0
   if(kel.gt.0) then
! there is a specifier, check if correct element or species
      kel=index(elnam,'#')
      if(kel.gt.0) then
! extract sublattice number 1-9 specification
         lk3=ichar(elnam(kel+1:kel+1))-ichar('0')
!         write(*,73)elnam(kel+1:kel+1),kel,elnam,lk3
!73       format('3D sublattice: "',a,'" position: ',i3,' in ',a,' : ',i3)
         elnam(kel:)=' '
      else
         lk3=0
      endif
      if(btest(propid(typty)%status,IDELSUFFIX)) then
!         write(*,*)'3D: elnam: ',kel,lk3,typty,elnam
         call find_element_by_name(elnam,iel)
         if(gx%bmperr.ne.0) then
            write(kou,*)'3D Unknown element ',elnam,&
                 ' in parameter type MQ, please reenter'
            gx%bmperr=0; goto 1000
         endif
         typty=100*typty+iel
      elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
! to know the constituents we must know the phase but as we do not know 
! the phase name yet but check the species exists !!!
         call find_species_by_name(elnam,isp)
         if(gx%bmperr.ne.0) then
! This is not an error, the species may not be selected!!!
!            write(kou,*)'Unknown species ',elnam,&
!                 ' in parameter type MQ, please reenter',gx%bmperr
            gx%bmperr=0; goto 1000
         endif
! convert from index to location, loksp
         loksp=species(isp)
!         write(*,69)'3D: conname: ',kel,lk3,typty,isp,loksp,elnam
69       format(a,5i4,a)
! extract sublattice after #
      else
!         write(kou,*)'This property has no specifier'
         gx%bmperr=4168; goto 1000
      endif
! this is the property type stored in property record
   else
! check if there should be a specifier !!
      if(btest(propid(typty)%status,IDELSUFFIX) .or. &
           btest(propid(typty)%status,IDCONSUFFIX)) then
         write(*,77)propid(typty)%symbol
77       format('3D Missing specifier for model parameter idenifier ',a)
         gx%bmperr=4169; goto 1000
      endif
   endif
! if the parameter symbol has a constituent specification check that now
   if(isp.gt.0) then
      k4=0
      do ll=1,phlista(lokph)%noofsubl
         if(lk3.eq.0 .or. lk3.eq.ll) then
            do kk=1,phlista(lokph)%nooffr(ll)
               k4=k4+1
               if(phlista(lokph)%constitlist(k4).eq.loksp) goto 80
            enddo
         elseif(ll.lt.lk3) then
            k4=k4+phlista(lokph)%nooffr(ll)
         endif
      enddo
! constituent not found
      write(kou,*)'No such constituent'
      gx%bmperr=4066; goto 1000
! constituent found in right sublattice
80    continue
      typty=100*typty+k4
!      write(*,81)'3D: found: ',typty1,typty,lk3,k4,loksp
81    format(a,10i4)
   endif
1000 continue
   return
 end subroutine get_parameter_typty

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine enter_species_property
!\begin{verbatim}
 subroutine enter_species_property(loksp,nspx,value)
! enter an extra species property for species loksp
   implicit none
   integer loksp,nspx
   double precision value
!\end{verbatim} %+
! this is illegal for species that are elements ...
   if(btest(splista(loksp)%status,SPEL) .or. &
        btest(splista(loksp)%status,SPVA)) then
!      write(*,*)'Illegal to set this for element species'
      gx%bmperr=4298
   elseif(.not.allocated(splista(loksp)%spextra)) then
      write(*,*)'3D this species has no allocated extra data'
      gx%bmperr=4399; goto 1000
   elseif(nspx.gt.size(splista(loksp)%spextra)) then
      write(*,*)'3D species has not sufficient extra data allocated ',nspx
      gx%bmperr=4399; goto 1000
   else
      splista(loksp)%spextra(nspx)=value
   endif
1000 continue
   return
 end subroutine enter_species_property

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine suspend_unstable_sets
!\begin{verbatim}
 subroutine suspend_unstable_sets(mode,ceq)
! suspend extra composition sets that are not stable
   implicit none
   integer mode
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokph,ics,lokcs
!   loop for all phases
   phases: do lokph=1,noofph
      if(phlista(lokph)%noofcs.eq.1) cycle phases
      sets: do ics=2,phlista(lokph)%noofcs
! never change first composition set, even if not stable
         lokcs=phlista(lokph)%linktocs(ics)
         if(ceq%phase_varres(lokcs)%phstate.gt.0) cycle sets
         ceq%phase_varres(lokcs)%phstate=PHSUS
      enddo sets
   enddo phases
1000 continue
   return
 end subroutine suspend_unstable_sets

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine set_uniquac_species
!\begin{verbatim}
 subroutine set_uniquac_species(loksp)
! set the status bit and allocates spexttra array
   implicit none
   integer loksp
!\end{verbatim}
! this is illegal for species that are elements ...
   if(btest(splista(loksp)%status,SPEL) .or. &
        btest(splista(loksp)%status,SPVA)) then
      gx%bmperr=4298
   else
      splista(loksp)%status=ibset(splista(loksp)%status,SPUQC)
      if(.not.allocated(splista(loksp)%spextra)) then
         allocate(splista(loksp)%spextra(2))
         splista(loksp)%spextra=one
      endif
   endif
1000 continue
   return
 end subroutine set_uniquac_species
 
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine enter_material
!\begin{verbatim}
 subroutine enter_material(cline,last,nv,xknown,ceq)
! enter a material from a database
! called from user i/f
   implicit none
   integer last,nv
   character cline*(*)
   double precision xknown(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer nel,j1,j2,j3
   character material*72,database*72,selel(20)*2,ext*4,alloy(20)*2
   character majorel*2,ftype*1,bline*128,elnam*2
   double precision xalloy(20),rest,xxx,xxy
   logical byte
! these are saved for use in a subsequent call
   save selel,majorel,ftype,xalloy
!
   if(.not.btest(globaldata%status,GSNOPHASE)) then
! Ask for new alloy composition:
      if(ftype.eq.'Y') then
         rest=1.0D2
         bline='Mass % of '
      else
         rest=one
         bline='Mole fraction of '
      endif
      j2=len_trim(bline)+2
      do j1=1,noofel
         if(ellista(j1)%symbol.eq.majorel) cycle
         bline(j2:)=ellista(j1)%symbol
         xxy=xalloy(j1)
60       continue
         call gparrdx(bline,cline,last,xxx,xxy,'?Enter Material')
         if(buperr.ne.0 .or. xxx.le.zero) then
            write(*,*)'Illegal value for composition'
            goto 60
         endif
         xalloy(j1)=xxx
         rest=rest-xxx
      enddo
   else
      ext='.TDB'
      call gparcx('Database: ',cline,last,1,database,' ','?Enter matrial')
! this extracts all element symbols from database
      call checkdb(database,ext,nel,selel)
      if(gx%bmperr.ne.0) goto 1000
      write(kou,70)(selel(nv),nv=1,nel)
70    format('Elements: ',15(a2,', '))
! ask for major component
      call gparcx('Major element or material: ',cline,last,1,majorel,' ',&
           '?Enter material')
      call capson(majorel)
      do nv=1,nel
         if(majorel.eq.selel(nv)) goto 100
      enddo
      write(*,*)'No such element in the database'
      gx%bmperr=4399
      goto 1000
100   continue
      call gparcdx('Input in mass percent? ',cline,last,1,ftype,'Y',&
           '?Enter material')
      if(ftype.eq.'Y') then
         rest=1.0D2
         write(*,102)'mass percent'
      else
         rest=one
         write(*,102)'mole fractions'
      endif
102   format('Input expected in ',a/)
110   continue
      call gparcx('First alloying element:',cline,last,1,alloy(1),' ',&
           '?Enter matrial')
      nv=0
      call capson(alloy(1))
      do j1=1,nel
         if(alloy(1).eq.selel(j1)) goto 200
      enddo
      write(*,*)'No such element in database'
      goto 110
!-----
200   continue
      do j1=1,nv
         if(alloy(nv+1).eq.alloy(j1)) then
            write(*,*)'Alloying element already entered'
            goto 250
         endif
      enddo
      nv=nv+1
220   continue
      if(ftype.eq.'Y') then
         call gparrdx('Mass percent: ',cline,last,xalloy(nv),one,&
              '?Enter material')
         if(buperr.ne.0) then
            write(*,*)'Give a numeric value'; buperr=0
            goto 220
         endif
      else
         call gparrdx('Mole fraction: ',cline,last,xalloy(nv),1.0D-2,&
              '?Enter material')
         if(buperr.ne.0) then
            write(*,*)'Give a numeric value'; buperr=0
            goto 220
         endif
      endif
      if(xalloy(nv).le.zero) then
         write(*,*)'Composition must be positive!'
         goto 220
      endif
      rest=rest-xalloy(nv)
      if(rest.le.zero) then
         write(*,240)'zero!!'
240      format('Your major component composition is less than ')
         gx%bmperr=4399; goto 1000
      elseif(rest.le.5.0D-1) then
         write(*,240)'half the system!!'
      endif
250 continue
      if(nv.eq.1) then
         call gparcx('Second alloying element:',cline,last,1,alloy(2),' ',&
              '?Enter material')
         if(alloy(2).eq.'  ') goto 500
      elseif(nv.eq.2) then
         call gparcx('Third alloying element:',cline,last,1,alloy(3),' ',&
              '?Enter material')
         if(alloy(3).eq.'  ') goto 500
      else
         call gparcx('Next alloying element:',cline,last,1,&
              alloy(nv+1),' ','?Enter material')
         if(alloy(nv+1).eq.'  ') goto 500
      endif
      call capson(alloy(nv+1))
      do j1=1,nel
         if(alloy(nv+1).eq.selel(j1)) goto 200
      enddo
      write(*,*)'No such element in database'
      goto 250
!----------------------
! read the database including the major element
500   continue
!      write(*,505)'Comp: ',nv,(alloy(j1),xalloy(j1),j1=1,nv)
505   format(a,i2,2x,8(a2,F8.4,', '))
      nv=nv+1
      alloy(nv)=majorel
      xalloy(nv)=rest
!      write(*,505)'3D m1: ',nv,(alloy(j1),xalloy(j1),j1=1,nv)
      call readtdb(database,nv,alloy)
      if(gx%bmperr.ne.0) goto 1000
! order the amounts in xalloy in alphabetical order
      byte=.true.
      order: do while(byte)
         byte=.false.
         do j1=1,nv
            do j2=j1+1,nv
               if(alloy(j1).gt.alloy(j2)) then
                  byte=.true.
                  elnam=alloy(j1)
                  alloy(j1)=alloy(j2)
                  alloy(j2)=elnam
                  xxx=xalloy(j1)
                  xalloy(j1)=xalloy(j2)
                  xalloy(j2)=xxx
!                  write(*,505)'3D m1: ',nv,(alloy(j3),xalloy(j3),j3=1,nv)
                  cycle order
               endif
            enddo
         enddo
      enddo order
! these are saved until another enter material command
      do j1=1,nv
         xknown(j1)=xalloy(j1)
      enddo
!      write(*,505)'3D m2: ',nv,(alloy(j1),xknown(j1),j1=1,nv)
510   format('3D em: ',10(a2,F6.3,1x))
   endif
!----------------------------------
! set conditions for composition (replace major by N=1)
   bline=' '
   j2=len_trim(bline)+2
   do j1=1,nv
      if(alloy(j1).eq.majorel) cycle
      if(ftype.eq.'Y') then
         bline(j2:)='W%('//trim(alloy(j1))//')='
      else
         bline(j2:)='X('//trim(alloy(j1))//')='
      endif
      j2=len_trim(bline)+1
      call wrinum(bline,j2,10,0,xalloy(j1))
      j2=j2+2
   enddo
   bline(j2:)=' N=1 '
   j2=len_trim(bline)+2
   write(*,*)'3D em: ',trim(bline)
! set_condition will increment j1
   j1=1
   call set_condition(bline,j1,ceq)
1000 continue
   return
 end subroutine enter_material

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

