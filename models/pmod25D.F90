!
! included in pmod25.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     8. State variable functions
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine store_putfun(name,lrot,nsymb,&
       istv,indstv,iref,iunit,idot)
! enter an expression of state variables
! name: character, name of state variable function
! lrot: pointer, to a putfun_node that is the root of the stored expression
! nsymb: integer, number of formal arguments
! istv: integer array, formal argument state variables typ
! indstv: 2D integer array, indices for the formal state variables
! iref: integer array, reference for the formal state variables
! iunit: integer array, unit of the formal state variables
   implicit none
   type(putfun_node), pointer :: lrot
   integer nsymb
   integer, dimension(*) :: istv,iref,iunit,idot
   integer, dimension(4,*) :: indstv
   character name*(*)
!\end{verbatim}
   integer jf
!    write(*,*)'store_putfun ',nsvfun
   nsvfun=nsvfun+1
   if(nsymb.gt.0) then
      allocate(svflista(nsvfun)%formal_arguments(8,nsymb))
      do jf=1,nsymb
         svflista(nsvfun)%formal_arguments(1,jf)=istv(jf)
         svflista(nsvfun)%formal_arguments(2,jf)=indstv(1,jf)
         svflista(nsvfun)%formal_arguments(3,jf)=indstv(2,jf)
         svflista(nsvfun)%formal_arguments(4,jf)=indstv(3,jf)
         svflista(nsvfun)%formal_arguments(5,jf)=indstv(4,jf)
         svflista(nsvfun)%formal_arguments(6,jf)=iref(jf)
         svflista(nsvfun)%formal_arguments(7,jf)=iunit(jf)
         svflista(nsvfun)%formal_arguments(8,jf)=idot(jf)
      enddo
   endif
   svflista(nsvfun)%name=name
   svflista(nsvfun)%linkpnode=>lrot
   svflista(nsvfun)%status=0
   svflista(nsvfun)%narg=nsymb
! this is the number of actual argument needed (like @P, @C and @S)
   svflista(nsvfun)%nactarg=0
! eqnoval indicate which equilibrium to use to get its value.
! default is 0 meaning current equilibria, can be changed by AMEND SYMBOL
   svflista(nsvfun)%eqnoval=0
1000 continue
   return
 end subroutine store_putfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine list_svfun(text,ipos,lrot,ceq)
! list a state variable function
   implicit none
   character text*(*)
   integer ipos,lrot
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! copied svflista(lrot)%formal_arguments(2..5,jt) to indices as gfortran error
   integer indstv(4),indices(4)
   character symbols(20)*32,afterdot*32
   integer js,jt,ip,istv,ii,kl
!    write(*,*)'list_svfun 1:',svflista(lrot)%narg
   if(lrot.le.0 .or. lrot.gt.nsvfun) then
      gx%bmperr=4140; goto 1000
   endif
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
!             write(*,*)'dot derivative: ',jt,&
!                  svflista(lrot)%formal_arguments(1,jt)
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
500 continue
   kl=len_trim(svflista(lrot)%name)
   text(ipos:ipos+kl+1)=svflista(lrot)%name(1:kl)//'= '
   ipos=ipos+kl+2
   call wrtfun(text,ipos,svflista(lrot)%linkpnode,symbols)
! where is pfnerr defined??
   if(pfnerr.ne.0) then
      write(kou,*)'Putfun error listing funtion ',pfnerr
      gx%bmperr=4142; goto 1000
   endif
!   text(ipos:ipos)=';'
!   ipos=ipos+1
1000 continue
   return
 end subroutine list_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine list_all_svfun(kou,ceq)
! list all state variable funtions
   implicit none
   integer kou
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character text*256
   integer ks,ipos
   write(kou,17)
17 format('List of all state variable symbols'/' No Name = expression ;')
   do ks=1,nsvfun
      ipos=1
      call list_svfun(text,ipos,ks,ceq)
      if(pfnerr.ne.0) then
         gx%bmperr=4142; pfnerr=0; goto 1000
      endif
      write(kou,76)ks,text(1:ipos-1)
76    format(i3,2x,a)
   enddo
1000 continue
   return
 end subroutine list_all_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine evaluate_all_svfun(kou,ceq)
! evaluate and list all functions
   implicit none
   integer kou
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character actual_arg(10)*24
   integer kf
   double precision val
   write(kou,75)
75 format('No  Name ',12x,'Value')
   do kf=1,nsvfun
! actual arguments needed if svflista(kf)%nactarg>0
      val=evaluate_svfun(kf,actual_arg,0,ceq)
      if(gx%bmperr.ne.0) goto 1000
      write(kou,77)kf,svflista(kf)%name,val
77    format(i3,1x,a,1x,1PE15.8)
   enddo
1000 continue
   return
 end subroutine evaluate_all_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 double precision function evaluate_svfun(lrot,actual_arg,mode,ceq)
! envaluate all funtions as they may depend on each other
! actual_arg are names of phases, components or species as @Pi, @Ci and @Si
! needed in some deferred formal parameters  (NOT IMPLEMENTED YET)
   implicit none
   integer lrot,mode
   character actual_arg(*)*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer indices1(4),indices2(4)
   double precision argval(20)
   integer jv,jt,istv,ieq,ii
   double precision value
   argval=zero
!    write(*,*)'evaluate_svfun ',lrot,svflista(lrot)%narg,svflista(lrot)%name
! locate function
   if(lrot.le.0 .or. lrot.gt.nsvfun) then
      gx%bmperr=4140; goto 1000
   endif
   if(svflista(lrot)%narg.eq.0) goto 300
! get values of arguments
   jv=0
   jt=0
100 continue
      jt=jt+1
      istv=svflista(lrot)%formal_arguments(1,jt)
      if(istv.lt.0) then
! if eqnoval nonzero it indicates from which equilibrium to get its value
         ieq=svflista(lrot)%eqnoval
         if(ieq.eq.0) then
            value=ceq%svfunres(-istv)
         else
            value=eqlista(ieq)%svfunres(-istv)
         endif
!          write(*,*)'evaluate_svfun symbol',ieq,value
      else
         do ii=1,4
            indices1(ii)=svflista(lrot)%formal_arguments(1+ii,jt)
         enddo
         if(svflista(lrot)%formal_arguments(8,jt).eq.0) then
! get state variable value
            call state_variable_val(istv,indices1,&
                 svflista(lrot)%formal_arguments(6,jt), &
                 svflista(lrot)%formal_arguments(7,jt),value,ceq)
         else
! calculate a derivative!!! need 2 state variables
            do ii=1,4
               indices2(ii)=svflista(lrot)%formal_arguments(1+ii,jt+1)
            enddo
            call state_var_value_derivative(istv,indices1, &
                 svflista(lrot)%formal_arguments(6,jt), &
                 svflista(lrot)%formal_arguments(7,jt), &
                 svflista(lrot)%formal_arguments(1,jt+1),indices2,&
                 svflista(lrot)%formal_arguments(6,jt+1), &
                 svflista(lrot)%formal_arguments(7,jt+1),value,ceq)
            jt=jt+1
         endif
      endif
      jv=jv+1
      argval(jv)=value
      if(jt.lt.svflista(lrot)%narg) goto 100
! all arguments evaluated (or no arguments needed)
300 continue
!    write(*,333)'evaluate_svfun ',svflista(lrot)%name,argval(1),argval(2)
!333 format(a,a,2(1PE15.6))
!   write(*,340)'evaluate svfun 1: ',mode,lrot
340 format(a,5i4)
   modeval: if(mode.eq.0 .and. btest(svflista(lrot)%status,SVFVAL)) then
! If mode=0 and SVFVAL set return the stored value
      value=ceq%svfunres(lrot)
!      write(*,350)'evaluate svfun 2: ',0,lrot,value
   elseif(mode.eq.0 .and. btest(svflista(lrot)%status,SVFEXT)) then
! if mode=0 and SVFEXT set use value from equilibrium eqno
      ieq=svflista(lrot)%eqnoval
      if(ceq%eqno.eq.ieq) then
         value=evalf(svflista(lrot)%linkpnode,argval)
         if(pfnerr.ne.0) then
            write(*,*)'evaluate_svfun putfunerror ',pfnerr
            gx%bmperr=4141; goto 1000
         endif
         ceq%svfunres(lrot)=value
!         write(*,350)'evaluate svfun 3: ',ieq,lrot,value
      else
         value=eqlista(ieq)%svfunres(lrot)
      endif
!      write(*,350)'evaluate svfun 4: ',ieq,lrot,value
350 format(a,2i3,1pe12.4)
   else
! if mode=1 always evaluate
      value=evalf(svflista(lrot)%linkpnode,argval)
      if(pfnerr.ne.0) then
         write(*,*)'evaluate_svfun putfunerror ',pfnerr
         gx%bmperr=4141; goto 1000
      endif
   endif modeval
! save value in current equilibrium
!   write(*,*)'25D eval_svfun: ',lrot,value,size(ceq%svfunres)
   ceq%svfunres(lrot)=value
   evaluate_svfun=value
1000 continue
   return
 end function evaluate_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>     9. Interactive things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine ask_phase_constitution(cline,last,iph,ics,lokcs,ceq)
! interactive input constitution of phase iph
   implicit none
   integer last,iph,ics,lokcs
   character cline*(*)
!\end{verbatim}
   character name1*24,quest*32
   double precision yarr(maxcons2),sites(maxsubl),qq(5),yyy,xxx,sss
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
! 1=entered, 2=fix, 3=dormant, 4=suspended, 5=hidden
   if(qph.lt.0 .and. test_phase_status(iph,ics,xxx,ceq).gt.3) goto 200
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
! ask if default
   call gparcd('Default constitution?',cline,last,1,ch1,chd,q1help)
   if(ch1.eq.'Y' .or. ch1.eq.'y') then
      call set_default_constitution(iph,ics,0,ceq)
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
      nycon: do nr=1,knl(ll)-1
         kkk=kkk+1
         loksp=phlista(lokph)%constitlist(kkk)
         line=splista(loksp)%symbol
         ip=len_trim(line)+1
         if(ll.gt.1) then
            line(ip:)='#'//char(ll+ichar('0'))
            ip=ip+2
         endif
         once=.true.
20        continue
         call gparrd(line(1:ip+2),cline,last,xxx,yarr(kkk),q1help)
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
            write(*,*)'Sum of fractions larger than unity'
            sss=sss+xxx
            yarr(kkk)=sss-1.0D-12*(knl(ll)-nr)
!            goto 20
         endif
         yarr(kkk)=max(xxx,1.0D-12)
      enddo nycon
! the last constituent is set to the rest
      kkk=kkk+1
      yarr(kkk)=max(sss,1.0D-12)
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
 subroutine enter_parameter_interactivly(cline,ip)
! enter a parameter from terminal or macro
! NOTE both for ordered and disordered fraction set !!
   implicit none
   integer ip
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
!   write(*,*)'25D: fractyp: ',kq,name1(1:kq)
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
!73       format('25D sublattice: "',a,'" position: ',i3,' in ',a,' : ',i3)
         elnam(kel:)=' '
      endif
      if(btest(propid(typty)%status,IDELSUFFIX)) then
!         write(*,*)'25D: elnam: ',kel,lk3,typty,elnam
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
!      write(*,*)'25D: conname: ',kel,lk3,typty,elnam
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
!      write(*,81)'25D: found: ',typty1,typty,lk3,k4,loksp
81    format(a,10i4)
   endif
!    write(*,*)'enter_parameter_inter 2: ',jph,lokph
! extract constituent array, remove final ) and decode
   name3=parname(lp2+1:)
   lp1=len_trim(name3)
! this removes the final )
   name3(lp1:)=' '
!
   call decode_constarr(name3,nsl,endm,nint,lint,ideg)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,83)'hej: ',name3(1:lp1),nint,(lint(2,i),i=1,nint)
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
! If parameter has no T dependendence just ask for value
   if(btest(propid(typty1)%status,IDNOTP)) then
      write(kou,*)'This parameter can only be a constant'
      call gparr('Value: ',cline,ip,xxx,zero,q1help)
      if(buperr.ne.0) then
         xxx=zero; buperr=0
      endif
! the tpfun always want a low T, expression; high T N
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
   call enter_tpfun(parname,longline,lfun)
   if(gx%bmperr.ne.0) goto 1000
!    write(*,290)'enter_par 7: ',lokph,nsl,nint,ideg,lfun,refx
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
      if(refid.eq.reflista(jl)%reference) then
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
 subroutine set_condition(cline,ip,ceq)
! decode an equilibrium condition, can be an expression with + and -
! the expression should be terminated with an = or value supplied on next line
! like "T=1000", "x(liq,s)-x(pyrrh,s)=0", "2*mu(cr)-3*mu(o)=muval"
! It can also be a "NOFIX=<phase>" or "FIX=<phase> value"
! The routine should accept conditions identified with the <number>:
! preceeding each condition in a list_condition
! It should also accept changing conditions by <number>:=new_value
   implicit none
   integer ip
   character cline*(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer nterm,kolon,iqz,krp,jp,istv,iref,iunit,jstv,jref,junit,jl,ks
   integer linkix,norem,ics,kstv,iph,nidfirst,nidlast,nidpre,qp
   character stvexp*80,stv*16,stvrest*80,textval*32,c4*4
   character svtext*60,encoded*60,defval*18
   integer indices(4),allterms(4,10),condno
   double precision coeffs(10),xxx,value
   logical inactivate
   TYPE(gtp_condition), pointer :: previous,temp,new
50 continue
! extract symbol like T, X(FCC,CR), MU(C) etc up to space, +, - or = sign
   call gparcd('State variable: ',cline,ip,ichar('='),stvexp,'T',q1help)
! can there be any error return??
   if(buperr.ne.0) goto 1000
!   write(*,*)'in set_condition: ',ip,stvexp(1:20)
! is the line empty?
   if(stvexp(1:1).eq.' ') then
      gx%bmperr=4126; goto 1000
   endif
!   write(*,*)'set_cond 0: ',stvexp(1:20)
   nterm=0
! check if there is a ":" in stvexp
!   write(*,*)'set_condition: ',stvexp(1:40),ip
   condno=0
   kolon=index(stvexp,':')
   if(kolon.gt.0) then
!      write(*,*)'25D: Found colon: ',kolon,stvexp(1:20)
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
!         write(*,*)'25D: Found condition: ',condno
         goto 155
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
!      write(*,*)'Inactivate phase condition'
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
!    write(*,*)'calling decode ',stvexp(1:20)
   call decode_state_variable(stvexp,istv,indices,iref,iunit,ceq)
   if(gx%bmperr.ne.0) goto 1000
   svtext=stvexp
   if(istv.lt.0) then
! istv < 0 means it is a parameter property symbol like TC, illegal as cond
      gx%bmperr=4127; goto 1000
   endif
! Use  kstv=(istv+1)/10+5 as check
   kstv=(istv+1)/10+5
!   write(*,*)'condition code: ',istv,kstv
!   if(kstv.eq.14 .or. kstv.eq.15 .or. kstv.ge.19) then
   if(kstv.eq.14 .or. kstv.eq.15) then
! the state variables Q and DG are illegal as conditions 
      gx%bmperr=4127; goto 1000
   endif
!   write(*,117)stvexp(1:10),istv,indices
!117 format('25D, sc: ',a,i7,5x,4i3)
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
   coeffs(nterm)=xxx
   do jl=1,4
      allterms(jl,nterm)=indices(jl)
   enddo
   stvexp=stvrest
!    write(*,*)'set_cond rest "',stvrest(1:20),'"'
   jp=1
   if(.not.eolch(stvexp,jp)) goto 100
!---------------
150 continue
! the expression (or single variable) is decoded, get value after = sign
! provide the current value as default
   encoded=' '
   call get_state_var_value(svtext,xxx,encoded,ceq)
   if(gx%bmperr.ne.0) then
! This error occurs when setting the first compositions before any calculations
!      write(*,152)gx%bmperr,svtext(1:len_trim(svtext))
!152   format('Cannot get current value of: ',a,', error: ',i5/&
!           'Setting default value to zero')
      gx%bmperr=0; xxx=zero
   endif
155 continue
   qp=1
   call wrinum(defval,qp,10,0,xxx)
   if(buperr.ne.0) then
      buperr=0; defval=' '
   endif
157 continue
   call gparcd('Value: ',cline,ip,1,textval,defval,q1help)
!   if(buperr.ne.0) then
!      gx%bmperr=4129; buperr=0; goto 1000
!   endif
! value can be a symbol or a numeric value.  Symbols not allowed yet
!    write(*,*)'set_condition value ',textval(1:20)
   inactivate=.false.
   jp=1
   c4=textval(jp:jp+3)
   call capson(c4)
   if(c4.eq.'NONE') then
! value NONE means inactivate
! Problem when setting T=NONE, ceq%tpval(1) could have any value afterwards
      inactivate=.true.
!      write(*,*)'Inactivate condition: ',condno,value,xxx
      value=xxx
   else
!      write(*,*)'Textval: ',textval
      call getrel(textval,jp,value)
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
! Some symbols must not be allowed for values .....
         gx%bmperr=4129; buperr=0; goto 1000
      else
         linkix=-1
      endif
170   continue
   endif
! check if condition already exists.  If condno>0 it must exist or error
   temp=>ceq%lastcondition
   if(condno.eq.0) then
      call get_condition(nterm,jstv,allterms,iref,iunit,temp)
   elseif(condno.lt.0) then
      if(inactivate) then
! user has given *=NONE, remove all conditions except phase status FIX
! as inactive conditions are ignored just remove first acive until error
         norem=1
210      continue
         temp=>ceq%lastcondition
         call get_condition(0,norem,allterms,iref,iunit,temp)
         if(gx%bmperr.ne.0) then
            gx%bmperr=0; goto 1000
         endif
         if(temp%statev.gt.0) then
!            write(*,*)'inactivate ',norem,temp%statev
            temp%active=1
         endif
         goto 210
! no else patth needed
      endif
   else
! user has given "5:= ..." and condition 5 must exist, otherwise error
      call get_condition(0,condno,allterms,iref,iunit,temp)
      if(gx%bmperr.ne.0) goto 1000
   endif
277 continue
!   write(*,*)'At label 277: ',gx%bmperr
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
         if(linkix.lt.0) then
            temp%prescribed=value
            temp%symlink1=-1
! special for T and P, change the local value also
!            write(*,*)'set condition ',istv,value
            if(istv.eq.1) then
               ceq%tpval(1)=value
            elseif(istv.eq.2) then
               ceq%tpval(2)=value
            endif
         else
            temp%symlink1=linkix
         endif
      endif
      goto 900
   endif newcond
   goto 500
! handle fix/nofix of a phase, a condition should be set inactive.
300 continue
   ip=ip+1
!   write(*,*)'fix phase: ',ip,cline(ip:40)
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
   temp=>ceq%lastcondition    
   call get_condition(nterm,jstv,allterms,iref,iunit,temp)
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
         temp%prescribed=value
         temp%symlink1=-1
      endif
      goto 1000
   else
! if inactivate it is an error not finding the condition
!      write(*,*)'Finding condition error ',gx%bmperr
      if(inactivate) then
         write(*,*)'We should have found a condition ',gx%bmperr
         goto 1000
      endif
      gx%bmperr=0
   endif
! get the value for the new condition
   ip=index(cline,'==')+2
   call getrel(cline,ip,value)
   if(buperr.ne.0) then
      write(*,*)'error setting fix amount ',ip,cline(1:40)
      gx%bmperr=4100; goto 1000
   endif
!   write(*,*)'Set fix phase amount: ',value
! create new condition record for this equilibrium
500 continue
   temp=>ceq%lastcondition    
!   write(*,*)'allocating condition'
   allocate(ceq%lastcondition)
   new=>ceq%lastcondition
   new%noofterms=nterm
   new%statev=jstv
   new%iunit=iunit
   new%iref=iref
   new%active=0
!    write(*,*)'allocating terms',nterm
   allocate(new%condcoeff(nterm))
   allocate(new%indices(4,nterm))
!    write(*,*)'allocations ok'
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
!      write(*,*)'set condition ',istv,jstv,value
      if(istv.eq.1) then
         ceq%tpval(1)=value
      elseif(istv.eq.2) then
         ceq%tpval(2)=value
      endif
   else
      new%symlink1=linkix
   endif
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
 end subroutine set_condition

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine apply_condition_value(current,what,type,value,cmix,ceq)
! This is called when calculating an equilibrium.
! It returns one condition at each call, at first current must be nullified
! When all conditions done the current is nullified again
! If what=-1 then return degrees of freedoms and maybe something more
! what=0 means calculate current values of conditions
! calculate the value of a condition, used in minimizing G
   implicit none
   integer what,type,cmix(*)
   double precision value
   TYPE(gtp_equilibrium_data), pointer :: ceq    
   TYPE(gtp_condition), pointer :: current
!\end{verbatim}
   integer, dimension(4) :: indices
   integer iref,iunit,degfree,jl,istv,ip
   double precision xxx
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
! ignore conditions with several terms
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
!      write(*,*)'Fix phase: ',-current%statev,current%iref
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

!\begin{verbatim}
 subroutine enter_svfun(cline,last,ceq)
! enter a state variable function
! no check that the name is unique!!!
   implicit none
   integer last
   character cline*(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer, parameter :: npfs=20
   integer ks,maxsym,ipos,jt,js,kdot,nsymb
   character name2*16,pfsym(npfs)*60,string*128,pfsymdot*60
   integer istv(npfs),indstv(4,npfs),iref(npfs),iunit(npfs),lokv(npfs)
   type(putfun_node), pointer :: lrot
!    
! maxsym is negative to allow the user to enter abs(maxs) symbols
! pfsym are the entered symbols
! lokv is only internal strage in putfun
! lrot is the root node of expression
! nsymb is the number of user entered symbols
!    write(kou,17)'enter svgun ',last,cline(1:20),nsvfun
17 format(a,i3,2x,a,i3)
   call gparc('Symbol name: ',cline,last,ichar('='),name2,' ',q1help)
   call capson(name2)
!   write(*,*)'enter_svfun: ',last,name2,':',cline(1:10)
   if(.not.proper_symbol_name(name2,0)) goto 1000
   do ks=1,nsvfun
      if(name2.eq.svflista(ks)%name) then
         gx%bmperr=4136; goto 1000
      endif
   enddo
! TO BE IMPLEMENTED: enter symbols with dummy arguments like CP(@P1)=HM(@P1).T
! where @Pi is a phase, @Ci is a component and @Si is a species
! these dummy variables must be defined in symbol name ?? why ?? maybe not
   call gparc('Expression, end with ";" :',cline,last,6,string,';',q1help)
   maxsym=-npfs
   ipos=1
   call putfun(string,ipos,maxsym,pfsym,lokv,lrot,nsymb)
   if(pfnerr.ne.0) then
      pfnerr=0; gx%bmperr=4134; goto 1000
   endif
!
! identify symbols as state variables, if derivative there is a dot
   jt=0
   do js=1,nsymb
      jt=jt+1
      lokv(jt)=0
      kdot=index(pfsym(js),'.')
      if(kdot.gt.0) then
! derivatives must be stored as two state variables
!          write(*,*)'Cannot handle: ',pfsym(js),kdot
!          goto 1000
         lokv(jt)=1
         pfsymdot=pfsym(js)(kdot+1:)
         pfsym(js)(kdot:)=' '
!          write(*,*)'dot ',pfsym(js)(1:20),pfsymdot(1:20)
         call decode_state_variable(pfsym(js),istv(jt),indstv(1,jt),&
              iref(jt),iunit(jt),ceq)
         if(gx%bmperr.ne.0) goto 1000
         jt=jt+1
         call decode_state_variable(pfsymdot,istv(jt),indstv(1,jt),&
              iref(jt),iunit(jt),ceq)
         if(gx%bmperr.ne.0) goto 1000
!          write(*,*)'enter_svfun ',istv(jt-1),istv(jt)
      else
         call decode_state_variable(pfsym(js),istv(jt),indstv(1,jt),&
              iref(jt),iunit(jt),ceq)
         if(gx%bmperr.ne.0) then
! symbol not a state variable, may be another function
!            write(*,*)'error: "',pfsym(js),'"'
            gx%bmperr=0
            do ks=1,nsvfun
               if(pfsym(js).eq.svflista(ks)%name) then
                  istv(jt)=-ks
                  goto 390
               endif
            enddo
            gx%bmperr=4135; goto 1000
         endif
      endif
390   continue
   enddo
! store name, all arguments and function, not if derivatives nt>nsymb
!    write(*,*)'enter_svfun ',jt,nsymb
   call store_putfun(name2,lrot,jt,istv,indstv,&
        iref,iunit,lokv)
1000 continue
   return
 end subroutine enter_svfun

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

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
   integer lokph,lokcs,ky,ll,iy
   real mmyfr(maxconst)
   character quest*32,name*24
   double precision xxx
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
! if PHNOCV set the composition is fixed
   if(btest(phlista(lokph)%status1,PHNOCV)) goto 1000
   write(*,10)
10 format('Give min or max fraction values (max as negative value)',&
        ' or zero for no default')
   name=' '
   ky=0
   do ll=1,phlista(lokph)%noofsubl
      do iy=1,phlista(lokph)%nooffr(ll)
         ky=ky+1
!         call get_species_name(phlista(lokph)%constitlist(ky),name)
         call get_phase_constituent_name(iph,ky,name)
         if(gx%bmperr.ne.0) then
            write(*,*)'default: ',iph,ky,iy
            goto 1000
         endif
         quest='Default for '//name(1:len_trim(name))//&
              '#'//char(ichar('0')+ll)
!         call gparrd(quest,cline,last,xxx,-1.0D-3,q1help)
         call gparrd(quest,cline,last,xxx,0.5D0,q1help)
         if(xxx.gt.one) xxx=one
         mmyfr(ky)=xxx
      enddo
   enddo
   call enter_default_constitution(iph,ics,mmyfr,ceq)
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
   TYPE(gtp_condition), pointer :: current,first,last
   character species*32,cval*16,statevar*4,condline*32
   integer ielno(10),indices(4)
   double precision addval(maxel)
   integer k,loksp,istv,jel,ip,iref,iunit
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
   call decode_state_variable('N ',istv,indices,iref,iunit,ceq)
! if B convert to N: moles of species = input_mass/mass_of_species
! moles of element = stoiciometry_of_element/total_number_of_elements
   if(statevar(1:1).eq.'B') then
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
! if condition not active then activate and zero prescibed amount
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

