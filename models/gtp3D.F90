!
! gtp3D included in gtp3.F90
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     8. Interactive things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
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
   call gparc('Phase name: ',cline,last,1,name1,' ',q1help)
   if(name1(1:2).eq.'* ') then
! this means all phases and composition sets
      qph=-1
      iph=1
      ics=1
      call get_phase_name(iph,ics,name1)
      if(gx%bmperr.ne.0) goto 1000
   else
      qph=0
      call find_phase_by_name(name1,iph,ics)
      if(gx%bmperr.ne.0) goto 1000
   endif
100 continue
!   write(*,*)'spc 1',qph,iph,ics,name1
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
   call gparrd(quest,cline,last,xxx,yyy,q1help)
! if input error quit asking more
   if(buperr.ne.0) then
      buperr=0; goto 1000
   endif
   ceq%phase_varres(lokcs)%amfu=xxx
! ask if we should set the default constitution
   call gparcd('Default constitution?',cline,last,1,ch1,chd,q1help)
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
         call gparrd(line(1:ip+2),cline,last,xxx,ydef,q1help)
         if(xxx.lt.zero) then
            if(once) then
               write(*,*)'A Fraction must be greater than zero'
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
   character name3*64,ch1*1,line*64,parname*64
   integer typty,lint(2,5),fractyp,typty1,kp,lp1,kel,kq,iel,isp,lk3,lp2
   integer jph,ics,lokph,ll,k4,nint,jp,lsc,ideg,kk,lfun,nsl,loksp
   integer, dimension(maxsubl) :: endm(maxsubl)
   double precision xxx
!
10  continue
   call gparc('Parameter name: ',cline,ip,7,parname,' ',q1help)
! simple parameter names are like G(SIGMA,FE:CR:FE,CR;1)
   kp=index(parname,' ')
   parname(kp:)=' '
! extract symbol, normally G or L but TC and others can occur 
! for example a mobility like  MQ&FE+2#3 where FE+2#3 is a constinuent
! in sublattice 3
   lp1=index(parname,'(')
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
   write(kou,*)'unknown parameter type, please reenter: ',&
        name1(1:len_trim(name1))
   parname=' '; goto 10
!
70 continue
   typty1=typty
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
      endif
      if(btest(propid(typty)%status,IDELSUFFIX)) then
!         write(*,*)'3D: elnam: ',kel,lk3,typty,elnam
         call find_element_by_name(elnam,iel)
         if(gx%bmperr.ne.0) then
            write(kou,*)'Unknown element ',elnam,&
                 ' in parameter type MQ, please reenter'
            parname=' '; gx%bmperr=0; goto 10
         endif
         typty=100*typty+iel
      elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
! to know the constituents we must know the phase but as we do not know 
! the phase name yet but check the species exists !!!
!      write(*,*)'3D: conname: ',kel,lk3,typty,elnam
         call find_species_by_name(elnam,isp)
         if(gx%bmperr.ne.0) then
            write(kou,*)'Unknown species ',elnam,&
                 ' in parameter type MQ, please reenter',gx%bmperr
            parname=' '; gx%bmperr=0; goto 10
         endif
! convert from index to location, loksp
         loksp=species(isp)
! extract sublattice after #
      else
         write(kou,*)'This property has no specifier'
         gx%bmperr=4168; goto 1000
      endif
! this is the property type stored in property record
   else
! check if there should be a specifier !!
      if(btest(propid(typty)%status,IDELSUFFIX) .or. &
           btest(propid(typty)%status,IDCONSUFFIX)) then
         write(*,*)'Parameter specifier missing'
         gx%bmperr=4169; goto 1000
      endif
   endif
!
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
         write(kou,*)'unknown phase name, please reenter'
         kp=len(cline)
         goto 10
      endif
   endif
   lokph=phases(jph)
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
!    write(*,*)'enter_parameter_inter 2: ',jph,lokph
! extract constituent array, remove final ) and decode
   name3=parname(lp2+1:)
   lp1=len_trim(name3)
! this removes the final )
   name3(lp1:)=' '
!
   call decode_constarr(lokph,name3,nsl,endm,nint,lint,ideg)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,83)'3D after d_c: ',name3(1:lp1),nint,(lint(2,kp),kp=1,nint)
83 format(a,a,i5,2x,5i4)
! finally remove all non-alphabetical characters in the function name by _
   kp=0
100 continue
   kp=kp+1
105 continue
   ch1=parname(kp:kp)
! should use ??
!   if(ucletter(ch1)) goto 100
   if(ch1.ge.'A' .and. ch1.le.'Z') goto 100
   if(ch1.ne.' ') then
      parname(kp:)=parname(kp+1:)
      goto 105
   endif
   parname='_'//parname
!-------------------------------------------------
! if mode=0 enter the parameter, 
! if mode=1 just list the parameter
! if mode=2 maybe amending (does FOOLED) work?
   if(mode.eq.1) then
      lfun=-1
!      write(*,*)'3D calling enter_parameter with lfun=',lfun
      call enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,&
           lfun,refx)
      goto 1000
   endif
! continue here to enter the parameter
! If parameter has no T dependendence just ask for value
   if(btest(propid(typty1)%status,IDNOTP)) then
      write(kou,*)'This parameter can only be a constant'
      call gparr('Value: ',cline,ip,xxx,zero,q1help)
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
   call gparr('Low  temperature limit /298.15/:',cline,ip,xxx,2.9815D2,q1help)
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
   call gparc('Expression, end with ";":',cline,ip,6,line,';',q1help)
   if(buperr.ne.0) then
      buperr=0; line=';'
   endif
120 continue
   longline(jp:)=line
   jp=len_trim(longline)+1
!   write(*,152)0,jp,longline(1:jp)
   if(index(longline(lsc:),';').le.0) then
      call gparc('&',cline,ip,6,line,';',q1help)
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
   call gparr('Upper temperature limit /6000/:',cline,ip,xxx,6.0D3,q1help)
   if(buperr.ne.0) then
      buperr=0; xxx=6.0D3
   endif
   call wrinum(longline,jp,8,0,xxx)
   if(buperr.ne.0) goto 1000
   call gparcd('Any more ranges',cline,ip,1,ch1,'N',q1help)
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
   call gparcd('Reference symbol:',cline,ip,1,refx,'UNKNOWN',q1help)
   call capson(refx)
   longline(jp:)=refx
   jp=len_trim(longline)+1
!    write(*,252)2,jp,longline(1:jp)
252 format('ep: ',2i3,'>',a,'<')
!
   call capson(longline(1:jp))
!   write(*,*)'epi: ',longline(1:jp)
   call enter_tpfun(parname,longline,lfun,.FALSE.)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,290)'enter_par 7: ',lokph,nsl,nint,ideg,lfun,refx
290 format(a,5i4,1x,a)
!
   call enter_parameter(lokph,typty,fractyp,nsl,endm,nint,lint,ideg,lfun,refx)
!
1000 continue
   return
 end subroutine enter_parameter_interactivly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine amend_global_data(cline,ipos)
   implicit none
   character cline*(*)
   integer ipos
!\end{verbatim}
   character name*24,current*24,ch1*1,chd*1
   current=globaldata%name
!   write(*,*)'entering amend_global_data: ',cline(1:30)
   call gparcd('System name: ',cline,ipos,1,name,current,q1help)
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
   call gparcd('I am a beginner (B), freqent user (F) or expert (E): ',&
        cline,ipos,1,ch1,chd,q1help)
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
   call gparcd('Global gridminimization allowed: ',&
        cline,ipos,1,ch1,chd,q1help)
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      globaldata%status=ibclr(globaldata%status,GSNOGLOB)
   else
      globaldata%status=ibset(globaldata%status,GSNOGLOB)
   endif
! allow merging gridpoints after global?
   chd='Y'
   if(btest(globaldata%status,GSNOMERGE)) chd='N'
   call gparcd('Merging gridpoints in same phase allowed: ',&
        cline,ipos,1,ch1,chd,q1help)
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      globaldata%status=ibclr(globaldata%status,GSNOMERGE)
   else
      globaldata%status=ibset(globaldata%status,GSNOMERGE)
   endif
! GSNOACS can be changed interactivly, 0 means allowed
   chd='Y'
   if(btest(globaldata%status,GSNOACS)) chd='N'
   call gparcd('Composition sets can be created automatically? ',&
        cline,ipos,1,ch1,chd,q1help)
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      globaldata%status=ibclr(globaldata%status,GSNOACS)
   else
      globaldata%status=ibset(globaldata%status,GSNOACS)
   endif
! GSNOREMCS can be changed interactivly, 0 means not remove
   chd='Y'
   if(btest(globaldata%status,GSNOREMCS)) chd='N'
   call gparcd('Delete unnecessary composition sets automatically? ',&
        cline,ipos,1,ch1,chd,q1help)
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      globaldata%status=ibclr(globaldata%status,GSNOREMCS)
   else
      globaldata%status=ibset(globaldata%status,GSNOREMCS)
   endif
1000 continue
   return
 end subroutine amend_global_data

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_reference_interactivly(cline,last,mode,iref)
! enter a reference for a parameter interactivly
! this should be modified to allow amending an existing reference
   implicit none
   character cline*(*)
   integer last,mode,iref
!\end{verbatim}
! stupid with a variable called L80
   character line*256,refid*16,L80*80
   integer jl,ip
   call gparc('Reference identifier:',cline,last,1,refid,' ',q1help)
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
100 continue
   call gparc('Reference text, end with ";":',cline,last,5,l80,';',q1help)
   line(ip:)=l80
   ip=len_trim(line)
   if(ip.le.1) then
      write(kou,*)'There must be some reference text!'
      ip=1; goto 100
   elseif(line(ip:ip).ne.';') then
      ip=ip+1; goto 100
   else
      line(ip:)=' '
   endif
   call tdbrefs(refid,line,1,iref)
1000 continue
   return
 end subroutine enter_reference_interactivly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine enter_experiment(cline,ip,ceq)
! enters an experiment, almost the same as set_condition   
   implicit none
   character cline*(*)
   integer ip
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbation}
! New is set to the new condition or experiment
! in set_condition new is not used for anything.
! in this subroutine the new variable is removed from the condition list
! and instead added to the experimenal list
   type(gtp_condition), pointer :: new,temp
!   integer nidlast,nidfirst,nidpre
   double precision xxx,yyy
!   write(*,*)'experiments not iumplemented yet'
   call set_cond_or_exp(cline,ip,new,1,ceq)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3D Back in enter_experiment'
!   write(*,*)'3D Values: ',new%active,new%prescribed
!   write(*,*)'3D Segmentation fault before this line means a problem with new'
   if(new%active.ne.1) then
! the experiment is removed (inactivated) if activate is 1
! otherwise read the current uncertainty
      if(new%uncertainty.gt.zero) yyy=1.0D-1*abs(new%prescribed)
      call gparrd('Uncertainty: ',cline,ip,xxx,yyy,q1help)
      if(gx%bmperr.ne.0) then
         xxx=0.1*new%prescribed
      endif
      new%symlink2=0
      new%uncertainty=abs(xxx)
! if weight is sero set to unity
      if(ceq%weight.eq.zero) then
         ceq%weight=one
      endif
   endif
1000 continue
 end subroutine enter_experiment

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine set_condition(cline,ip,ceq)
! enters an experiment, almost the same as set_condition   
   implicit none
   character cline*(*)
   integer ip
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbation}
! New is set to the new condition or experiment
! in this subroutine new is not used for anything.
! in enter_experiment the new variable is removed from the condition list
! and instead added to the experimenal list
   type(gtp_condition), pointer :: new
   call set_cond_or_exp(cline,ip,new,0,ceq)
1000 continue
 end subroutine set_condition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 logical function same_statevariable(svr1,svr2)
! returns TRUE if the state variable records are identical
   type(gtp_state_variable), pointer :: svr1,svr2
!end{verbatim}
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

!\begin{verbatim}
 subroutine set_cond_or_exp(cline,ip,new,notcond,ceq)
! decode an equilibrium condition, can be an expression with + and -
! the expression should be terminated with an = or value supplied on next line
! like "T=1000", "x(liq,s)-x(pyrrh,s)=0", "2*mu(cr)-3*mu(o)=muval"
! It can also be a "NOFIX=<phase>" or "FIX=<phase> value"
! The routine should also accept conditions identified with the "<number>:"
! where <number> is that preceeding each condition in a list_condition
! It should also accept changing conditions by <number>:=new_value
! The pointer to the (most recent) condition is returned in new
! notcond is 0 if a condition should be created, otherwise an experiment
   implicit none
   integer ip,notcond
   character cline*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_condition), pointer :: new
!\end{verbatim}
   integer nterm,kolon,iqz,krp,jp,istv,iref,iunit,jstv,jref,junit,jl,ks
   integer linkix,norem,ics,kstv,iph,nidfirst,nidlast,nidpre,qp
   character stvexp*80,stvrest*80,textval*32,c4*4
   character svtext*60,encoded*60,defval*18
   integer indices(4),allterms(4,10),condno,seqz
   double precision coeffs(10),xxx,value
   logical inactivate
   type(gtp_state_variable), pointer :: svr,svr2
   type(gtp_state_variable), dimension(10), target :: svrarr
   TYPE(gtp_condition), pointer :: temp
!   TYPE(gtp_condition), pointer :: temp,new
50 continue
! extract symbol like T, X(FCC,CR), MU(C) etc up to space, +, - or = sign
   stvexp=' '
   call gparcd('State variable: ',cline,ip,ichar('='),stvexp,'T',q1help)
! can there be any error return??
   if(buperr.ne.0) goto 1000
!   write(*,*)'3D in set_condition: ',ip,stvexp(1:20)
! is the line empty?
   if(stvexp(1:1).eq.' ') then
      gx%bmperr=4126; goto 1000
   endif
!   write(*,*)'3D set_cond 0: ',stvexp(1:20)
   nterm=0
! check if there is a ":" in stvexp
!   write(*,*)'3D set_condition: ',stvexp(1:40),ip
   condno=0
   kolon=index(stvexp,':')
   if(kolon.gt.0) then
!      write(*,*)'3D: Found colon: ',kolon,stvexp(1:20)
!      if(len_trim(stvexp).gt.kolon) then
! is there any text after : ignore the text before ": "  WHY????
!         stvexp=stvexp(kolon+2:)
!      else
! the user specifies the condition by giving its number like "5:=none"
         iqz=1
         call getrel(stvexp,iqz,xxx)
         if(buperr.ne.0) then
            gx%bmperr=buperr; goto 1000
         endif
         condno=int(xxx)
!         write(*,*)'3D: Found condition number: ',condno,ip
! We cannot provide any default value
         defval=' '
         goto 157
!         goto 155
!      endif
   endif
   if(stvexp(1:1).eq.'*') then
! user can remove all conditions (except phase status) by the line
! *=NONE
      condno=-1
      goto 157
   endif
! check for phase status FIX or NOFIX, these are generated by change_status
! not given by user
   if(stvexp(1:3).eq.'FIX') then
      inactivate=.FALSE.
      goto 300
   elseif(stvexp(1:5).eq.'NOFIX') then
      inactivate=.TRUE.
!      write(*,*)'3D Inactivate phase fix condition'
      goto 300
   endif
! check if it is an expression with + or -
100 continue
! look for a ) followed by + or -
   krp=index(stvexp,')')
   if(krp.gt.0) then
      if(stvexp(krp+1:krp+1).eq.'+' .or. stvexp(krp+1:krp+1).eq.'-') then
         stvrest=stvexp(krp+1:)
         stvexp(krp+1:)=' '
      else
         stvrest=' '
      endif
   else
      stvrest=' '
   endif
!    write(*,*)'set_cond 2: ',krp,stvexp(1:20),':',stvrest(1:20)
! there can be a factor in front of the state variable
   jp=1
   call getrel(stvexp,jp,xxx)
   if(buperr.ne.0) then
!       write(*,*)'buperr ',buperr,jp
! 1035 and 1036 means a sign not followed by digits
      if(buperr.eq.1035) then
         xxx=one
         stvexp=stvexp(2:)
      elseif(buperr.eq.1036) then
         xxx=-one
         stvexp=stvexp(2:)
      else
         xxx=one
      endif
      buperr=0
   else
      if(stvexp(jp:jp).eq.'*') then
         stvexp=stvexp(jp+1:)
      else
         gx%bmperr=4130; goto 1000
      endif
   endif
! decode state variable expression
!    write(*,*)'3D, calling decode ',stvexp(1:20)
!   call decode_state_variable(stvexp,istv,indices,iref,iunit,svr,ceq)
   call decode_state_variable(stvexp,svr,ceq)
   if(gx%bmperr.ne.0) goto 1000
   svtext=stvexp
! convert to old format
!   if(svr%oldstv.ge.10) then
!      istv=10*(svr%oldstv-5)+svr%norm
!   else
   istv=svr%oldstv
!   endif
   iref=svr%phref
   iunit=svr%unit
! svr%argtyp specifies values in indices:
! svr%argtyp: 0=no arguments; 1=comp; 2=ph+cs; 3=ph+cs+comp; 4=ph+cs+const
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
!   else
!      write(*,*)'3D state variable has illegal argtyp: ',svr%argtyp
!      gx%bmperr=7775; goto 1000
   endif
!
   if(istv.lt.0) then
! istv < 0 means it is a parameter property symbol like TC, illegal as cond
      gx%bmperr=4127; goto 1000
   endif
! Use  kstv=(istv+1)/10+5 as check
   kstv=(istv+1)/10+5
!   write(*,*)'3D condition code: ',istv,kstv
!   if(kstv.eq.14 .or. kstv.eq.15 .or. kstv.ge.19) then
   if(kstv.eq.14 .or. kstv.eq.15) then
! the state variables Q and DG are illegal as conditions 
      gx%bmperr=4127; goto 1000
   endif
!   write(*,117)stvexp(1:10),istv,indices
!117 format('3D, sc: ',a,i7,5x,4i3)
   if(istv.ge.3 .and. istv.le.5) then
! this is MU, AC and LNAC, do not allow phase index (at present)
      if(indices(2).ne.0) then
      write(*,*)'Phase specific chemical potentials not allowed as conditions'
         gx%bmperr=7777; goto 1000
      endif
   endif
   if(nterm.gt.0) then
! it must be the same state variable in all terms
      if(istv.ne.jstv .or. iref.ne.jref .or. iunit.ne.junit) then
         gx%bmperr=4128; goto 1000
      endif
   else
      jstv=istv
      jref=iref
      junit=iunit
   endif
   nterm=nterm+1
! save svr record if several terms
   svrarr(nterm)=svr
   coeffs(nterm)=xxx
   do jl=1,4
      allterms(jl,nterm)=indices(jl)
   enddo
   stvexp=stvrest
!    write(*,*)'3D set_cond rest "',stvrest(1:20),'"'
   jp=1
   if(.not.eolch(stvexp,jp)) goto 100
!---------------
150 continue
!   write(*,*)'3D at 150',notcond
! the expression (or single variable) is decoded, get value after = sign
! for conditions provide the current value as default
! for experiments provide the old value if any ???
   if(notcond.eq.0) then
      encoded=' '
      call get_state_var_value(svtext,xxx,encoded,ceq)
      if(gx%bmperr.ne.0) then
! This error occurs when setting the first compositions before any calculations
!      write(*,152)gx%bmperr,svtext(1:len_trim(svtext))
!152   format('3D Cannot get current value of: ',a,', error: ',i5/&
!           'Setting default value to zero')
         gx%bmperr=0; xxx=zero
      endif
   else
! search for experimental value
!      write(*,*)'3D Check if old or current value',svr%statevarid
      nullify(new)
      temp=>ceq%lastexperiment
      if(.not.associated(temp)) then
         xxx=zero
      else
! loop for all conditions
151      continue
         temp=>temp%next
         svr2=>temp%statvar(1)
!         write(*,*)'3D we have an old experiment ',svr2%statevarid
         if(same_statevariable(svr,svr2)) then
            xxx=temp%prescribed
!            write(*,*)'3D We have found an old experiment',svr2%statevarid
            new=>temp
         else
            if(.not.associated(temp,ceq%lastexperiment)) goto 151
         endif
      endif
   endif
! jump here if we have found a : specifying a condition number
! Do not expect any = sign
155 continue
!   write(*,*)'3D we are here',notcond,xxx
   qp=1
   defval=' '
   call wrinum(defval,qp,10,0,xxx)
   if(buperr.ne.0) then
      buperr=0; defval=' '
   endif
157 continue
!   write(*,*)'3D we are here: ',defval
   call gparcd('Value: ',cline,ip,1,textval,defval,q1help)
!   write(*,*)'3D extract value: ',cline(1:ip),textval(1:len_trim(textval)),ip
!   if(buperr.ne.0) then
!      gx%bmperr=4129; buperr=0; goto 1000
!   endif
! value can be a symbol or a numeric value.  Symbols not allowed yet
!   write(*,*)'3D set cond/enter exper value: ',textval(1:20),notcond
   inactivate=.false.
   jp=1
   c4=textval(jp:jp+3)
   call capson(c4)
   inactive: if(c4.eq.'NONE') then
! value NONE means inactivate
! Problem when setting T=NONE, ceq%tpval(1) could have any value afterwards
      inactivate=.true.
!      write(*,*)'Inactivate condition: ',condno,value,xxx
      value=xxx
!      write(*,*)'3D inactivate'
   else
! It is not inactivate............
      if(notcond.ne.0) then
! If we have an experiment check if there is a : in textval
         call getrel(textval,jp,value)
         if(buperr.ne.0) then
! not a numeric value, error as symbols not allowed
            gx%bmperr=buperr; goto 1000
         endif
         linkix=-1
         ks=index(textval,':')
         if(ks.gt.0) then
! backspace to make ip point at the :
            ip=ip-(len_trim(textval)-ks)-1
!            write(*,*)'3D backspace: ',cline(1:ip),ip
         endif
!      endif
      else
! problem that the : for uncertainties is swollowed 
!         write(*,*)'gtp3D textval 1: ',textval(1:jp),jp
         call getrel(textval,jp,value)
!         write(*,*)'gtp3D textval 2: ',textval(1:jp),jp,buperr
         if(buperr.ne.0) then
! here one should look for a symbol
            buperr=0
            call capson(textval)
            do ks=1,nsvfun
               if(textval(1:16).eq.svflista(ks)%name) then
! found symbol, insert link
                  linkix=ks
                  goto 170
               endif
            enddo
! neither numeric value nor symbol, give error
! Some symbols are allowed for values .....
            gx%bmperr=4129; buperr=0; goto 1000
         else
            linkix=-1
         endif
      endif
170   continue
!      write(*,*)'gtp3D value: ',value
   endif inactive
! check if condition/experiment already exists. 
! If condno>0 it must exist or error
!   write(*,*)'3D we are here'
   if(notcond.eq.0) then
      temp=>ceq%lastcondition
   else
      if(associated(new)) then
!         write(*,*)'3D we have found the old experiment, just return',&
!              inactivate,new%active
! if experiment should be inactivated set that, else activate
         if(inactivate) then
            new%active=1
         else
            new%active=0
            new%prescribed=value
         endif
         goto 1000
      else
         temp=>ceq%lastexperiment
      endif
   endif
!   write(*,*)'3D condno: ',condno,notcond,value,buperr
   if(condno.eq.0) then
!      call get_condition(nterm,coeffs,jstv,allterms,iref,iunit,temp)
      call get_condition(nterm,svr,temp)
   elseif(condno.lt.0) then
      if(notcond.eq.0) then
         write(*,*)'3D Illegal to set all experiments to NONE'
         gx%bmperr=7766; goto 1000
      endif
      if(inactivate) then
! user has given *=NONE, remove all conditions except phase status FIX
! as inactive conditions are ignored just remove first acive until error
         norem=1
210      continue
         temp=>ceq%lastcondition
!         call get_condition(0,coeffs,norem,allterms,iref,iunit,temp)
         call get_condition(nterm,svr,temp)
         if(gx%bmperr.ne.0) then
            gx%bmperr=0; goto 1000
         endif
         if(temp%statev.gt.0) then
!            write(*,*)'3D inactivate ',norem,temp%statev
            temp%active=1
         endif
         goto 210
! no else path needed
      endif
   else
! user has given "5:= ..." and condition 5 must exist, otherwise error
!      call get_condition(0,coeffs,condno,allterms,iref,iunit,temp)
      qp=-condno
!      write(*,*)'3D calling get_condition ',qp
      call get_condition(qp,svr,temp)
      if(gx%bmperr.ne.0) then
         write(*,*)'3D Condition number error ',gx%bmperr
         goto 1000
      endif
   endif
277 continue
!   write(*,*)'3D At label 277: ',gx%bmperr,linkix
   newcond: if(gx%bmperr.ne.0) then
! condition does not exist, but if set equal to NONE just ignore it
      gx%bmperr=0
      if(inactivate) goto 900
   else !newcond
      if(inactivate) then
! remove condition
         temp%active=1
!         write(*,*)'inactivating: ',ceq%tpval(1)
      else
! condition already exist, just change value and activate
         temp%active=0
         temp%symlink1=-1
         if(linkix.lt.0) then
! special for T and P, change the local value also
!            write(*,*)'set condition ',istv,value
            if(istv.eq.1) then
               if(value.lt.1.0D-1) then
! illegal value of T
                  write(*,*)'Temparure must be larger than 0.1 K'
                  gx%bmperr=4187; goto 1000
               endif
               ceq%tpval(1)=value
               temp%prescribed=value
            elseif(istv.eq.2) then
               if(value.lt.1.0D-1) then
! illegal value of P
                  write(*,*)'Pressure must be larger than 0.1 Pa'
                  gx%bmperr=4187; goto 1000
               endif
               ceq%tpval(2)=value
               temp%prescribed=value
            else
!               write(*,*)'gtp3D: condition',value
               temp%prescribed=value
            endif
         else
            temp%symlink1=linkix
         endif
      endif
      goto 900
   endif newcond
   goto 500
!-----------------------------------------------------------------
! handle fix/nofix of a phase, a condition should be set inactive.
300 continue
   ip=ip+1
   if(notcond.ne.0) then
      write(kou,*)'Illegal to set a fix phase as experiment'
      gx%bmperr=7677; goto 1000
   endif
!   write(*,*)'3D fix phase: ',ip,cline(ip:40)
   call find_phase_by_name(cline(ip:),iph,ics)
   if(gx%bmperr.ne.0) goto 1000
   nterm=1
   jstv=-iph
   iref=ics
   iunit=0
   linkix=-1
   coeffs(1)=1.0D0
   do jl=1,4
      allterms(jl,1)=0
   enddo
! convert to state variable
!   write(*,*)'3D Setting svrarr(1) values'
   svrarr(1)%statevarid=jstv
   svrarr(1)%oldstv=jstv
   svrarr(1)%phase=ics
   svrarr(1)%unit=0
   svrarr(1)%argtyp=0
   svrarr(1)%phase=iph
   svrarr(1)%compset=ics
   svrarr(1)%component=0
   svrarr(1)%constituent=0
!
   if(notcond.eq.0) then
      temp=>ceq%lastcondition    
   else
      temp=>ceq%lastexperiment
   endif
!   write(*,*)'3D calling get_condition'
   svr=>svrarr(1)
   call get_condition(nterm,svr,temp)
!   write(*,*)'3D Back from get_condition ',gx%bmperr
   if(gx%bmperr.eq.0) then
      if(inactivate) then
! inactivate condition
         temp%active=1
      else
! set new value of prescribed amount, must be numerical, not symbol
         temp%active=0
         ip=index(cline,'==')+2
         call getrel(cline,ip,value)
         if(buperr.ne.0) then
!            write(*,*)'error setting fix amount old cond',ip,cline(1:40)
            gx%bmperr=4100; goto 1000
         endif
! value set explict,, not a symbol
         temp%prescribed=value
         temp%symlink1=-1
      endif
      goto 1000
   else
! if inactivate it is an error not finding the condition
!      write(*,*)'Finding condition error ',gx%bmperr
      if(inactivate) then
         write(*,*)'3D We should have found a condition ',gx%bmperr
         goto 1000
      endif
      gx%bmperr=0
   endif
! Here we create a new condition !!!
! get the value for the new condition
   ip=index(cline,'==')+2
   call getrel(cline,ip,value)
   if(buperr.ne.0) then
      write(*,*)'3D error setting fix amount ',ip,cline(1:40)
      gx%bmperr=4100; goto 1000
   endif
!   write(*,*)'3D Set fix phase amount: ',value
!-----------------------------------------------
! create a new condition record for this equilibrium (can be the first)
500 continue
! first test if condition on P or T is larger than 0.1
   if(jstv.eq.1 .or. jstv.eq.2) then
      if(value.lt.0.1D0) then
         gx%bmperr=4187; goto 1000
      endif
   endif
   if(notcond.eq.0) then
      if(associated(ceq%lastcondition)) then
         seqz=ceq%lastcondition%seqz+1
      else
         seqz=1
      endif
      temp=>ceq%lastcondition
      allocate(ceq%lastcondition)
      new=>ceq%lastcondition
   else
! it is an experiment
!      write(*,*)'3D We are after label 500',value,linkix
      if(associated(ceq%lastexperiment)) then
         seqz=ceq%lastexperiment%seqz+1
      else
         seqz=1
      endif
      temp=>ceq%lastexperiment
      allocate(ceq%lastexperiment)
      new=>ceq%lastexperiment
   endif
!   write(*,*)'3D allocating condition',seqz  
!   allocate(ceq%lastcondition)
!   new=>ceq%lastcondition
   new%noofterms=nterm
   new%statev=jstv
   new%iunit=iunit
   new%iref=iref
   new%active=0
   new%seqz=seqz
!    write(*,*)'allocating terms',nterm
   allocate(new%condcoeff(nterm))
   allocate(new%indices(4,nterm))
!   write(*,*)'3D allocations ok',linkix,value
   do jl=1,nterm
      new%condcoeff(jl)=coeffs(jl)
      do ks=1,4
         new%indices(ks,jl)=allterms(ks,jl)
      enddo
   enddo
   if(linkix.lt.0) then
      new%prescribed=value
      new%symlink1=-1
! special for T and P, change the local value
!      write(*,*)'set condition/enter experiment ',istv,jstv,value
      if(istv.eq.1) then
         ceq%tpval(1)=value
      elseif(istv.eq.2) then
         ceq%tpval(2)=value
      endif
   else
      new%symlink1=linkix
   endif
! store the state variable record in the condition
   allocate(new%statvar(nterm))
   do jl=1,nterm
      new%statvar(jl)=svrarr(jl)
   enddo
! link the new record into the condition list
!    write(*,*)'linking condition'
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
! we are actually entering an experiment, terminate here
!      write(*,*)'3D Finished creating condition record for an experiment ...'
      goto 1000
   endif
900 continue
   if(.not.eolch(cline,ip)) then
! look for more conditions.  Note that gparc increment ip by 1 at start
!       write(*,901)cline(ip-1:ip+20)
901   format(' >',a,"<")
      if(cline(ip:ip).ne.',') ip=ip-1
      goto 50
   endif
! finally, for conditions on T and P copy value to ceq%tpval
! This may be a bit inconsistent .... but??
   if(jstv.eq.1 .and. iunit.eq.0 .and. iref.eq.0) then
      ceq%tpval(1)=value
   elseif(jstv.eq.2 .and. iunit.eq.0 .and. iref.eq.0) then
      ceq%tpval(2)=value
   endif
! mark that any current results may be inconsistent with new conditions
!   globaldata%status=ibset(globaldata%status,GSINCON)
   ceq%status=ibset(ceq%status,EQINCON)
1000 continue
!   write(*,*)'exit set_condition, T= ',ceq%tpval(1)
   return
 end subroutine set_cond_or_exp

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine get_condition(nterm,svr,pcond)
! finds a condition record with the given state variable expression
! If nterm<0 the absolute value of nterm is condition number, svr is irrelevant
   implicit none
   integer nterm
   type(gtp_state_variable), pointer :: svr
! NOTE: pcond must have been set to ceq%lastcondition before calling this
! pcond: pointer, to a gtp_condition record for this equilibrium
   type(gtp_condition), pointer :: pcond
!\end{verbatim} %+
   type(gtp_condition), pointer :: last
   type(gtp_state_variable), pointer :: condvar
   integer j1,num
   if(.not.associated(pcond)) goto 900
!   write(*,*)'3D in get_condition: ',svr%statevarid,svr%oldstv,svr%argtyp
!   if(nterm.lt.0) write(*,*)'3D Condition number: ',-nterm
   last=>pcond
   num=0
100 continue
! search for condition abs(nterm)
      if(nterm.lt.0 .and. num+nterm.eq.0) goto 1000
      num=num+1
      if(pcond%noofterms.eq.nterm) then
         do j1=1,nterm
            condvar=>pcond%statvar(j1)
!            write(*,*)'3D get_condition: ',num,condvar%oldstv,condvar%argtyp
! dissapointment, one cannot compare two structures ... unless pointers same
!            if(condvar.ne.svr) goto 200
!            j2=1
            if(condvar%oldstv.ne.svr%oldstv) goto 200
!            j2=2
            if(condvar%argtyp.ne.svr%argtyp) goto 200
!            j2=3
            if(condvar%phase.ne.svr%phase) goto 200
!            j2=4
            if(condvar%compset.ne.svr%compset) goto 200
            if(condvar%statevarid.lt.0) goto 1000
! for fix phase the remaining have no importance
!            j2=5
            if(condvar%component.ne.svr%component) goto 200
!            j2=6
            if(condvar%constituent.ne.svr%constituent) goto 200
!            j2=7
            if(condvar%norm.ne.svr%norm) goto 200
!            j2=8
            if(condvar%unit.ne.svr%unit) goto 200
         enddo
! we have found a condition with these state variables
!         write(*,*)'3D Found condition',pcond%active
         goto 1000
      endif
200   continue
!      write(*,*)'Failed at argument: ',j2
      pcond=>pcond%next
      if(.not.associated(pcond,last)) goto 100
900 continue
!   write(*,*)'3D get_condition: No such condition'
   gx%bmperr=7779; goto 1000
1000 continue
 return
end subroutine get_condition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
      if(ocv()) write(*,*)'get_condition: ',istv,ncc,nac
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

!\begin{verbatim} %-
 subroutine extract_stvr_of_condition(pcond,nterm)
! subroutine extract_stvr_of_condition(pcond,nterm,coeffs,statevar)
! finds a condition record with the given state variable record
! returns it as a state variable record !!!
! nterm: integer, number of terms in the condition expression
! pcond: pointer, to a gtp_condition record
   implicit none
   TYPE(gtp_condition), pointer :: pcond
! ONE CANNOT HAVE ARRAYS OF POINTERS!!! STUPID
!   TYPE(gtp_state_variable), dimension(*), pointer :: statevar
!   TYPE(gtp_state_variable), dimension(*) :: statevar
   integer nterm
!   double precision coeffs(*)
!\end{verbatim}
   TYPE(gtp_condition), pointer :: current,first
!   integer, dimension(4) :: indx
   integer ncc,nac,j1,istv,iref,iunit
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

!\begin{verbatim}
 subroutine locate_condition(seqz,pcond,ceq)
! locate a condition using a sequential number
   implicit none
   integer seqz
   type(gtp_condition), pointer :: pcond
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ij
   pcond=>ceq%lastcondition
   do ij=1,seqz
      pcond=>pcond%next
      if(seqz.gt.ij .and. associated(pcond,ceq%lastcondition)) then
         write(*,*)'Locate condition called with too high index: ',seqz
         gx%bmperr=7777; goto 1000
      endif
   enddo
1000 continue
   return
 end subroutine locate_condition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine apply_condition_value(current,what,value,cmix)
! This is called when calculating an equilibrium.
! It returns a condition at each call, at first call current must be nullified?
! When all conditions done the current is nullified again
! If what=-1 then return degrees of freedoms and maybe something more
! what=0 means calculate current values of conditions
! calculate the value of a condition, used in minimizing G
   implicit none
   integer what,cmix(*)
   double precision value
!   TYPE(gtp_equilibrium_data), pointer :: ceq    
   TYPE(gtp_condition), pointer :: current
!\end{verbatim} %+
! ceq is actually redundant as current is a pointer to condition list in ceq
   integer, dimension(4) :: indices
   integer iref,iunit,jl,istv,ip
   character encoded*60,actual_arg*60
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
! cannot hanlde conditions with several terms
      write(*,*)'Found condition with several terms'
      gx%bmperr=7777; goto 900
   endif
! for debugging
   istv=current%statev
   do jl=1,4
      indices(jl)=current%indices(jl,1)
   enddo
   iref=current%iref
   iunit=current%iunit
   ip=1
   encoded=' '
   actual_arg=' '
!------------------
   if(current%statev.lt.0) then
! a fix phase cpndition has state variable equal to -iph, ics is stored in iref
      cmix(1)=4
      cmix(2)=-current%statev
      cmix(3)=current%iref
      value=current%prescribed
!      write(*,*)'3D Fix phase: ',-current%statev,current%iref,value
   elseif(current%statev.eq.1) then
! temperature
      cmix(1)=1
      value=current%prescribed
!      write(*,*)'conditon on T'
   elseif(current%statev.eq.2) then
! pressure
      cmix(1)=2
      value=current%prescribed
!      write(*,*)'conditon on P'
   elseif(current%statev.le.5) then
! potentials has statev=1..5 (T, P, MU, AC, LNAC)
      cmix(1)=3
      cmix(2)=current%statev
      cmix(3)=current%indices(1,1)
      value=current%prescribed
!      write(*,*)'condition on MU/AC/LNAC'
   elseif(current%statev.ge.10) then
! other condition must be on extensive properties (N, X, H etc)
      cmix(1)=5
!      write(*,*)'Extensive condition: ',current%statev
   else
      write(*,*)'Illegal condition',current%statev
      gx%bmperr=7777; goto 1000
   endif
   goto 900
!--------------------------------------
! Here we should return extensive condition, maybe calculate value
200 if(what.ne.0) goto 300
   cmix(1)=0
   if(current%noofterms.gt.1) then
! ignore conditions with several terms
      write(*,*)'Found condition with several terms'
      gx%bmperr=8888; goto 1000
   endif
! for debugging
   istv=current%statev
   do jl=1,4
      indices(jl)=current%indices(jl,1)
   enddo
   iref=current%iref
   iunit=current%iunit
   ip=1
   encoded=' '
   actual_arg=' '
!------------------
   if(current%statev.lt.10) goto 900
! condition must be on extensive properties (N, X, H etc)
   cmix(1)=5
   cmix(2)=current%statev
! indices are dimensioned (4,nterms)
   cmix(3)=current%indices(1,1)
   cmix(4)=current%indices(2,1)
   cmix(5)=current%indices(3,1)
   cmix(6)=current%indices(4,1)
   value=current%prescribed
   goto 900
!--------------------------------------
! this part is redundant ....
300   continue
   write(*,*)'Calling apply_condition with illegal option'
   gx%bmperr=8888; goto 1000
!-----------------------------------------------------------
! maybe something common
900 continue
!
1000 continue
   return
 end subroutine apply_condition_value

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\


!\begin{verbatim}
 subroutine amend_components(cline,last,ceq)
! enter a new set of components for equilibrium ceq
   implicit none
   integer last
   character cline*(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character symb*24
   integer lokic(maxel),ielno(10)
   double precision stoi(maxel,maxel+1),invstoi(maxel,maxel),spst(10)
   integer ic,loksp,nspel,jl,j2,ierr
   double precision spmass,qsp
! input is a list of species name, same as number of elements
   stoi=zero
   ic=0
100 continue
   call gparc('Give all components: ',cline,last,1,symb,' ',q1help)
   call find_species_record(symb,loksp)
   if(gx%bmperr.ne.0) goto 1000
! check not same species twice
   do jl=1,ic
      if(loksp.eq.lokic(jl)) then
!         write(*,*)'Same species twice'
         gx%bmperr=4162; goto 1000
      endif
   enddo
   ic=ic+1
   lokic(ic)=loksp
! get the stoichiometry and save it in row in stoi
   call get_species_data(loksp,nspel,ielno,spst,spmass,qsp)
   if(gx%bmperr.ne.0) goto 1000
   do jl=1,nspel
      stoi(ic,ielno(jl))=spst(jl)
   enddo
   if(ic.lt.noofel) goto 100
! check that stoichiometry matrix not singular, should calculate the inverse
!   do i=1,ic
!      write(*,200)(stoi(i,j),j=1,ic)
!200   format('X: ',6(1PE12.4))
!   enddo
! lukasnum routine to invert matrix
   call mdinv(maxel,maxel+1,stoi,invstoi,ic,ierr)
! check the matrix and its inverse
!   do i=1,ic
!      write(*,200)(invstoi(i,j),j=1,ic)
!   enddo
   if(ierr.eq.0) then
!      write(*,*)'Component matrix singular'
      gx%bmperr=4163; goto 1000
   endif
   if(allocated(ceq%compstoi)) then
      deallocate(ceq%compstoi)
      deallocate(ceq%invcompstoi)
   endif
   allocate(ceq%compstoi(ic,ic))
   allocate(ceq%invcompstoi(ic,ic))
!   write(*,*)(lokic(i),i=1,ic)
   do jl=1,ic
      ceq%complist(jl)%splink=lokic(jl)
! phlink=0 means no user defined reference state
      ceq%complist(jl)%phlink=0
      ceq%complist(jl)%status=0
      ceq%complist(jl)%tpref=zero
      ceq%complist(jl)%chempot=zero
      ceq%complist(jl)%mass=spmass
      do j2=1,ic
         ceq%compstoi(jl,j2)=stoi(jl,j2)
         ceq%invcompstoi(jl,j2)=invstoi(jl,j2)
      enddo
   enddo
1000 continue
   return
 end subroutine amend_components

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

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
   double precision xxx
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
! if PHNOCV set the composition is fixed
   if(btest(phlista(lokph)%status1,PHNOCV)) goto 1000
   write(*,10)
10 format('Give min or max fractions (< or negative value as max)',&
        ' or NONE for no default')
   name=' '
   ky=0
   do ll=1,phlista(lokph)%noofsubl
      if(phlista(lokph)%nooffr(ll).gt.1) then
! more than one constituent
         do iy=1,phlista(lokph)%nooffr(ll)
            ky=ky+1
            call get_phase_constituent_name(iph,ky,name,subl)
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
!
            call gparcd(quest,cline,last,1,fdef,vdef,q1help)
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
99 format('3D: ',15(f5.1))
1000 continue
   return
 end subroutine ask_default_constitution

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

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
   call gparc('Species amount as N(..) or B(...): ',&
        cline,lpos,1,species,' ',q1help)
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
         gx%bmperr=7777; goto 1000
      endif
      cval=species(k+1:)
      species=species(3:k-1)
!      write(*,*)'3D Species: ',species
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
   if(cval(1:2).eq.'= ') goto 200
   goto 300
200 continue
! the user can also give values without = or with a space before =
! but no space allowed after =
   call gparc('Amount: ',cline,lpos,1,cval,' ',q1help)
300 continue
   if(cval(1:1).eq.'=') cval(1:1)=' '
   ip=1
   call getrel(cval,ip,xval)
   if(buperr.ne.0) then
      write(*,*)'Amount must be a real number'
      goto 1000
   endif
! this return the internal code for N
!   call decode_state_variable('N ',istv,indices,iref,iunit,svr,ceq)
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
   return
 end subroutine set_input_amounts

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine get_parameter_typty(name1,lokph,typty,fractyp)
! interpret parameter identifiers like MQ&C#2 in MQ&C#2(FCC_A1,FE:C) ...
! find the property associated with this symbol
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
! no matching symbol
   gx%bmperr=7777; goto 1000
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
            write(kou,*)'Unknown element ',elnam,&
                 ' in parameter type MQ, please reenter'
            gx%bmperr=0; goto 1000
         endif
         typty=100*typty+iel
      elseif(btest(propid(typty)%status,IDCONSUFFIX)) then
! to know the constituents we must know the phase but as we do not know 
! the phase name yet but check the species exists !!!
         call find_species_by_name(elnam,isp)
         if(gx%bmperr.ne.0) then
            write(kou,*)'Unknown species ',elnam,&
                 ' in parameter type MQ, please reenter',gx%bmperr
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
         write(*,*)'Parameter specifier missing'
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
!
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

