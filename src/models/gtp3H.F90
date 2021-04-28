! gtp3H included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
!>     14. Additions and diffusion
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
! Additions have a unique number, given sequentially as implemented 
! These are all defined in gtp3.F90
!  integer, public, parameter :: INDENMAGNETIC=1
!  integer, public, parameter :: XIONGMAGNETIC=2
!  integer, public, parameter :: DEBYECP=3
!  integer, public, parameter :: EINSTEINCP=4
!  integer, public, parameter :: TWOSTATEMODEL1=5
!  integer, public, parameter :: ELASTICMODEL1=6
!  integer, public, parameter :: VOLMOD1=7
!  integer, public, parameter :: UNUSED_CRYSTBREAKDOWNMOD=8
!  integer, public, parameter :: SECONDEINSTEIN=9
!  integer, public, parameter :: SCHOTTKYANOMALY=10
!  integer, public, parameter :: DIFFCOEFS=11
!  integer, public, parameter :: TWOSTATEMODEL2=12
!------------------------------------
! For each addition XX there is a subroutine create_XX
! called from the add_addrecord
! and a subroutine calc_XX 
! called from the addition_selector, called from calcg_internal
! There is a common list routine
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine addition_selector
!\begin{verbatim}
 subroutine addition_selector(addrec,moded,phres,lokph,mc,ceq)
! called when finding an addition record while calculating G for a phase
! addrec is addition record
! moded is 0, 1 or 2 if no, first or 2nd order derivatives should be calculated
! phres is ?
! lokph is phase location
! mc is number of constitution fractions
! ceq is current equilibrium record
   implicit none
   type(gtp_phase_add), pointer :: addrec
   integer moded,lokph,mc
   TYPE(gtp_phase_varres), pointer :: phres
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
!   write(*,*)'3H select addition: ',addrec%type
   addition: select case(addrec%type)
   case default
      write(kou,*)'3H No such addition type ',addrec%type,lokph
      gx%bmperr=4330
! 1 Inden-Hillert magnetic
   case(indenmagnetic) ! Inden magnetic 
      addrec%propval=zero
      call calc_magnetic_inden(moded,phres,addrec,lokph,mc,ceq)
! 2 Inden-Hillert-Qing-Xiong magnetism
   case(xiongmagnetic) ! Inden-Qing-Xiong
      addrec%propval=zero
      call calc_xiongmagnetic(moded,phres,addrec,lokph,mc,ceq)
!      write(kou,*)'3H Inden-Qing-Xiong magn model not tested yet'
!      gx%bmperr=4332
! 3 Debye Cp
   case(debyecp) ! Debye Cp
      addrec%propval=zero
      call calc_debyecp(moded,phres,addrec,lokph,mc,ceq)
      write(kou,*)'3H Debye Cp model not implemented yet'
      gx%bmperr=4331
! 4 Einsten Cp
   case(einsteincp) ! Einstein Cp
      addrec%propval=zero
      call calc_einsteincp(moded,phres,addrec,lokph,mc,ceq)
!      gx%bmperr=4331
! 5  Twostate liquid
   case(twostatemodel1) ! Two state model with composition variable G2
      addrec%propval=zero
!      write(*,*)'3H selecting calc_twostate_model1: ',mc
      call calc_twostate_model1(moded,phres,addrec,lokph,mc,ceq)
! changed below not to calculate G2 as a mixing parameter
!      call calc_twostate_model2(moded,phres,addrec,lokph,mc,ceq)
! 6 Elastic model
   case(elasticmodel1) ! Elastic model !
      addrec%propval=zero
      call calc_elastica(moded,phres,addrec,lokph,mc,ceq)
      write(kou,*)' Elastic model not implemented yet'
      gx%bmperr=4399
! 7 Volume model
   case(volmod1) ! Simple volume model depending on V0, VA and VB
      addrec%propval=zero
      call calc_volmod1(moded,phres,addrec,lokph,mc,ceq)
! 8 UNUSED
!   case(crystalbreakdownmod) ! Limiting heat capacity of extrapolated solid
!      addrec%propval=zero
!      call calc_crystalbreakdownmod(moded,phres,addrec,lokph,mc,ceq)
! 9
   case(secondeinstein) ! Adding a second Einstein Cp
      addrec%propval=zero
      call calc_secondeinstein(moded,phres,addrec,lokph,mc,ceq)
! 10
   case(schottkyanomaly) ! Adding a second Schottky anomaly Cp
      addrec%propval=zero
      call calc_schottky_anomaly(moded,phres,addrec,lokph,mc,ceq)
! 11
   case(diffcoefs) ! Calculating diffusion coefficients
      addrec%propval=zero
      call calc_diffusion(moded,phres,addrec,lokph,mc,ceq)
!      gx%bmperr=4333
! 12  see also 5, twostatemodel1 !
   case(twostatemodel2) ! Two state model with composition independent G2
      addrec%propval=zero
!      write(*,*)'3H selecting calc_twostate_model1: ',mc
!      call calc_twostate_model1(moded,phres,addrec,lokph,mc,ceq)
! changed not to calculate G2 as a mixing parameter
      call calc_twostate_model2(moded,phres,addrec,lokph,mc,ceq)
   end select addition
1000 continue
   return
 end subroutine addition_selector

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\addtotable subroutine add_addrecord
!\begin{verbatim}
 subroutine add_addrecord(lokph,extra,addtyp)
! generic subroutine to add an addition typ addtyp (including Inden)
   implicit none
   integer lokph,addtyp
   character extra*(*)
!\end{verbatim}
   integer aff
   double precision xxx
   character name*24
   type(gtp_phase_add), pointer :: newadd,addrec,lastrec
   logical bcc,ibm
!
!   write(*,*)'3H creating addrecord: ',trim(extra),addtyp,lokph
! check if this addition already entered
   lastrec=>phlista(lokph)%additions
   addrec=>lastrec
   do while(associated(addrec))
      if(addrec%type.eq.addtyp) then
         write(*,*)'3H addition already entered ',trim(phlista(lokph)%name),&
              addtyp,lokph,extra              
         goto 1000
      else
         lastrec=>addrec
         addrec=>lastrec%nextadd
      endif
   enddo
! NOTE EET is not an addition, it is comparing the entropy of solid and liquid
! create addition record
!   write(*,*)'3H adding addition record',lokph,addtyp
   addition: select case(addtyp)
!-----------------------------------------
   case default
      write(kou,*)'No addtion type ',addtyp,lokph
!-----------------------------------------
   case(indenmagnetic) ! Inden magnetic
! 1
      if(extra(1:1).eq.'Y' .or. extra(1:1).eq.'y') then
! bcc model
         aff=-1
         call create_magrec_inden(newadd,aff)
      else
         aff=-3
         call create_magrec_inden(newadd,aff)
      endif
!-----------------------------------------
   case(xiongmagnetic) ! Inden-Qing-Xiong. Assume bcc if BCC part of phase name
! 2
!      bcc=.false.
!      if(index('BCC',phlista(lokph)%name).gt.0) bcc=.true.
      if(extra(1:1).eq.'Y' .or. extra(1:1).eq.'y') then
         bcc=.TRUE.
      else
         bcc=.FALSE.
      endif
      ibm=.FALSE.
      if(extra(2:2).eq.'I') ibm=.TRUE.
! lokph because we need to check if average or individual Boghr magnetons
!      call create_xiongmagnetic(newadd,lokph,bcc)
      call create_xiongmagnetic(newadd,ibm,bcc)
!-----------------------------------------
   case(debyecp) ! Debye Cp UNUSED
! 3
!      call create_debyecp(newadd)
!-----------------------------------------
   case(einsteincp) ! Einstein Cp
! 4
      call create_einsteincp(newadd)
!-----------------------------------------
   case(twostatemodel1) ! Liquid 2 state model
! 5
! NEW set bit to allow endmember parameter modification! Question asked in PMON6
!      write(*,*)'3H setting bit ph2state: ',PH2STATE
!      phlista(lokph)%status1=ibset(phlista(lokph)%status1,PH2STATE)
!      write(*,*)'3H extra "',extra,'"'
      if(extra(1:1).eq.'N') then
         phlista(lokph)%status1=ibset(phlista(lokph)%status1,PH2STATE)
         write(*,*)'3H G2 is assumed to be composition independent'
         call create_newtwostate_model1(newadd)
! return that the addition type has changed ...
         addtyp=twostatemodel2
      else
         phlista(lokph)%status1=ibclr(phlista(lokph)%status1,PH2STATE)
         call create_twostate_model1(newadd)
      endif
!-----------------------------------------
   case(elasticmodel1) ! Elastic model 1
! 6
      call create_elastic_model_a(newadd)
!-----------------------------------------
   case(volmod1) ! Volume model 1
! 7
      call create_volmod1(newadd)
!-----------------------------------------
!   case(crystalbreakdownmod) ! Crystal Breakdown model
! 8 UNUSED
!      call create_crystalbreakdownmod(newadd)
!-----------------------------------------
   case(secondeinstein) ! Second Einstein T
! 9
      call create_secondeinstein(newadd)
!-----------------------------------------
   case(schottkyanomaly) ! Schottky anomaly
! 10
      call create_schottky_anomaly(newadd)
!-----------------------------------------
   case(diffcoefs)  !  diffusion coefficients
! 11 
      call create_diffusion(newadd,lokph,extra)
!-----------------------------------------
   end select addition
!-----------------------------------------
   if(gx%bmperr.ne.0) goto 1000
! initiate status word for this addition
!   newadd%status=0
   if(associated(phlista(lokph)%additions)) then
!      write(*,*)'3H adding new addition record to phase  ',lokph,addtyp
      lastrec%nextadd=>newadd
   else
!      write(*,*)'3H adding first addition record to phase',lokph,addtyp
      phlista(lokph)%additions=>newadd
   endif
1000 return
 end subroutine add_addrecord

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine need_propertyid
!\begin{verbatim}
 subroutine need_propertyid(id,typty)
! get the index of the property needed
   implicit none
   integer typty
   character*4 id
!\end{verbatim} %+
! here the property list is searched for "id" and its index stored in addrec
   do typty=1,ndefprop
      if(propid(typty)%symbol.eq.id) then
         goto 1000
      endif
   enddo
   write(*,*)'3H Parameter id ',id,' not found'
   gx%bmperr=4335
   typty=-1
1000 continue
   return
 end subroutine need_propertyid

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine setpermolebit
!\begin{verbatim}
 subroutine setpermolebit(lokph,addtype)
! set bit in addition record that addition is per mole
! lokph is phase record
! addtype is the addtion record type
   implicit none
   integer lokph,addtype
!\end{verbatim}
   type(gtp_phase_add), pointer :: addrec
   addrec=>phlista(lokph)%additions
!   write(*,*)'3H set size bit: ',addtype
   do while(associated(addrec))
      if(addrec%type.eq.addtype) then
         write(*,*)'3H setting bit ADDPERMOL for addition type ',addtype
         addrec%status=ibset(addrec%status,ADDPERMOL)
         goto 1000
      endif
      addrec=>addrec%nextadd
   enddo
   write(*,*)'3H Cannot find addition ',addtype
1000 continue
   return
 end subroutine setpermolebit

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_magrec_inden
!\begin{verbatim}
 subroutine create_magrec_inden(addrec,aff)
! enters the magnetic model
   implicit none
   type(gtp_phase_add), pointer :: addrec
   integer aff
!\end{verbatim} %+
   integer typty,ip,nc
   character text*128
   integer, parameter :: ncc=6
   double precision coeff(ncc)
   integer koder(5,ncc)
! There is some trouble with memory leaks in expressions to fix!!!
   TYPE(tpfun_expression), target :: llow2,lhigh2
   TYPE(tpfun_expression), pointer :: llow,lhigh
!
   if(aff.eq.-1) then
! bcc, aff=-1
! Magnetic function below Curie Temperature
! problem in ct1xfn to start a function with +1 or 1
      text=' 1.0-.905299383*T**(-1)-.153008346*T**3-'//&
           '.00680037095*T**9-.00153008346*T**15 ;'
!       write(*,*)'3H emm 1: ',text(1:len_trim(text))
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
!       write(*,17)'3H emm 1B:',nc,(coeff(i),i=1,nc)
17     format(a,i3,5(1PE11.3))
      if(gx%bmperr.ne.0) goto 1000
! Trouble with memory leaks for expressions to be fixed ...
      llow=>llow2
!      call ct1mexpr(nc,coeff,koder,llow)
! Attempt to remove big memory leak
      call ct1mexpr(nc,coeff,koder,llow2)
      if(gx%bmperr.ne.0) goto 1000
! Magnetic function above Curie Temperature
      text=' -.0641731208*T**(-5)-.00203724193*T**(-15)'//&
           '-4.27820805E-04*T**(-25) ; '
!       write(*,*)'3H emm 2: ',text(1:len_trim(text))
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
!      call ct1mexpr(nc,coeff,koder,lhigh)
! Attempt to remove big memory leak
      call ct1mexpr(nc,coeff,koder,lhigh2)
      if(gx%bmperr.ne.0) goto 1000
   else
!------------
! fcc, aff=-3
! Magnetic function below Curie Temperature
      text='+1.0-.860338755*T**(-1)-.17449124*T**3-.00775516624*T**9'//&
           '-.0017449124*T**15 ; '
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
      llow=>llow2
!      call ct1mexpr(nc,coeff,koder,llow)
! Attempt to remove big memory leak
      call ct1mexpr(nc,coeff,koder,llow2)
      if(gx%bmperr.ne.0) goto 1000
! Magnetic function above Curie Temperature
      text='-.0426902268*T**(-5)-.0013552453*T**(-15)'//&
           '-2.84601512E-04*T**(-25) ; '
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
!      call ct1mexpr(nc,coeff,koder,lhigh)
 ! Attempt to remove big memory leak
      call ct1mexpr(nc,coeff,koder,lhigh2)
     if(gx%bmperr.ne.0) goto 1000
   endif
! reserve an addition record
   allocate(addrec)
! store data in record
   allocate(addrec%explink(2))
   nullify(addrec%nextadd)
   addrec%status=0
   addrec%aff=aff
   addrec%type=indenmagnetic
! attempt to remove memory leak
!   addrec%explink(1)=llow
!   addrec%explink(2)=lhigh
!   write(*,*)'3H magnetic expression links'
   addrec%explink(1)=llow2
   addrec%explink(2)=lhigh2
! addrecs declared in gtp3.F90 but I am not sure it is needed or used
   addrecs=addrecs+1
   allocate(addrec%need_property(2))
   addrec%addrecno=addrecs
   addrec%need_property=0
! here the property list is searched for TC and BM
   call need_propertyid('TC  ',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(1)=typty
   call need_propertyid('BMAG',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(2)=typty
1000 continue
   return
 end subroutine create_magrec_inden

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_magnetic_inden
!\begin{verbatim}
 subroutine calc_magnetic_inden(moded,phres,lokadd,lokph,mc,ceq)
! calculates Indens magnetic contribution
! NOTE: values for function not saved, should be done to save time.
! Gmagn = RT*f(T/Tc)*ln(beta+1)
! moded: integer, 0=only G, S, Cp; 1=G and dG/dy; 2=Gm dG/dy and d2G/dy2
! phres: pointer, to phase\_varres record
! lokadd: pointer, to addition record
! lokph: integer, phase record 
! mc: integer, number of constituents
! ceq: pointer, to gtp_equilibrium_data
   implicit none
   integer moded,lokph,mc
   TYPE(gtp_phase_varres) :: phres
   TYPE(gtp_phase_add), pointer :: lokadd
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer itc,ibm,jl,noprop,ik,k,jk,j,jxsym
   double precision logb1,invb1,iafftc,iaffbm,rgasm,rt,tao,gmagn,msize
   double precision dtaodt,dtaodp,beta,d2taodp2,d2taodtdp,tc,tv
   double precision tao2(2),ftao(6),dtao(3,mc),d2tao(mc*(mc+1)/2)
   double precision addgval(6),daddgval(3,mc),d2addgval(mc*(mc+1)/2)
   logical addpermole
! phres points to result record with gval etc for this phase
   TYPE(tpfun_expression), pointer :: exprot
! dgdt = Gmagn/T + RT*df/dtao*dtao/dT*ln(beta+1)
! dgdp = RT df/dtao*dtao/dP*ln(beta+1)
! dgdy = RT*df/dtao*dtao/dy*ln(beta+1) + RT*f/(beta+1)*dbeta/dy
! d2gdt2=2*R*df/dtao*dtao/dT*ln(beta+1) + RT*d2f/dtao2*(dtao/dT)**2*ln(beta+1)
!        +RT*df/dtao*d2tao/dT2*ln(beta+1)
! d2gdtdp= ...
! d2gdp2=
! d2gdtdy=
! d2gdpdy=
! d2gdydy=
! listprop(1) is the number of properties calculated
! listprop(2:listprop(1)) give the typty of different properties
! calculated in gval(*,i) etc
! one has to find those with typty equal for need_property in the magnetic
! record, i.e. typty=2 for TC and typty=3 for BM
! the properties needed.
!
   noprop=phres%listprop(1)-1
   itc=0; ibm=0
!   write(*,*)'3H cmi 2: ',mc,noprop,(phres%listprop(j),j=1,noprop)
! Inden magnetic need properties in need_property(1..2)
   findix: do jl=2,noprop
      if(phres%listprop(jl).eq.lokadd%need_property(1)) then
         itc=jl
      elseif(phres%listprop(jl).eq.lokadd%need_property(2)) then
         ibm=jl
      endif
   enddo findix
   if(itc.eq.0 .or. ibm.eq.0) then
! it is no error if no TC or BM but then magnetic contribution is zero
!       write(*,12)phlista(lokph)%name
12     format('3H Warning: Magnetic addition for phase ',a&
           /9x,'but no values for TC or BM, magnetic contribution zero')
      goto 1000
   endif
   tc=phres%gval(1,itc)
   beta=phres%gval(1,ibm)
!    write(*,95)'3H Magnetic values in: ',itc,ibm,tc,beta
!95 format(a,2i3,3(1PE15.6))
   if(tc.lt.zero) then
! we should take care of the case when tc and beta have different signs
! note: all derivatives of tc must be multiplied with iaff
      iafftc=one/lokadd%aff
      do ik=1,mc
         do k=1,3
            phres%dgval(k,ik,itc)=iafftc*phres%dgval(k,ik,itc)
         enddo
         do jk=ik,mc
            jxsym=kxsym(ik,jk)
            phres%d2gval(jxsym,itc)=&
                 iafftc*phres%d2gval(jxsym,itc)
!            phres%d2gval(ixsym(ik,jk),itc)=&
!                 iafftc*phres%d2gval(ixsym(ik,jk),itc)
         enddo
      enddo
      do k=1,6
         phres%gval(k,itc)=iafftc*phres%gval(k,itc)
      enddo
      tc=phres%gval(1,itc)
!      write(*,*)'3H Inden 1: ',tc,iafftc
   else
      iafftc=zero
   endif
! avoid diving with zero, tc is a temperature so 0.01 degree is small
   if(tc.lt.one) tc=1.0D-2
   if(beta.lt.zero) then
! note all derivatives of bm must be multipled by iaffbm
!      iaffbm=one/addlista(lokadd)%aff
      iaffbm=one/lokadd%aff
      do ik=1,mc
         do k=1,3
            phres%dgval(k,ik,ibm)=iaffbm*phres%dgval(k,ik,ibm)
         enddo
         do jk=ik,mc
            jxsym=kxsym(ik,jk)
            phres%d2gval(jxsym,ibm)=&
                 iaffbm*phres%d2gval(jxsym,ibm)
         enddo
      enddo
      do k=1,6
         phres%gval(k,ibm)=iaffbm*phres%gval(k,ibm)
      enddo
      beta=phres%gval(1,ibm)
!      write(*,*)'3H Inden 2: ',beta,iaffbm
   endif
!
   tv=ceq%tpval(1)
   rgasm=globaldata%rgas
   rt=rgasm*tv
   tao=tv/tc
   tao2(1)=tao
! one should save values of ftao if tao2 is the same next time ....
! but as tc depend on the constitution that is maybe not so often.
   if(tao.lt.one) then
      exprot=>lokadd%explink(1)
   else
      exprot=>lokadd%explink(2)
   endif
   call ct1efn(exprot,tao2,ftao,ceq%eq_tpres)
   logb1=log(beta+one)
   invb1=one/(beta+one)
   gmagn=rt*ftao(1)*logb1
!   if(ocv()) then
!      write(*,98)'3H m1: ',tc,beta,ftao(1),logb1,rt
!      write(*,98)'3H m2: ',rt*gmagn,rt*(gmagn+phres%gval(1,1)),iafftc
!98    format(a,5(1PE14.6))
!   endif
!
   dtaodt=one/tc
   dtaodp=-tao/tc*phres%gval(3,itc)
   addgval(1)=gmagn
   addgval(2)=gmagn/tv+rt*ftao(2)*dtaodt*logb1
   addgval(3)=rt*ftao(2)*dtaodp*logb1+rt*ftao(1)*invb1*phres%gval(3,ibm)
!      phres%gval(1,1)=phres%gval(1,1)+addgval(1)/rt
!      phres%gval(2,1)=phres%gval(2,1)+addgval(2)/rt
!      phres%gval(3,1)=phres%gval(3,1)+addgval(3)/rt
! save these in record
   do j=1,3
      lokadd%propval(j)=addgval(j)
      phres%gval(j,1)=phres%gval(j,1)+addgval(j)/rt
   enddo
!   write(*,77)lokadd%type,(lokadd%propval(j),j=1,4)
!77 format('3H Addition ',i2,': ',4(1pe12.4))
! ignore second derivatives if no derivatives wanted
   if(moded.eq.0) then
      goto 1000
   endif
! Now all derivatives
! phres%gval(*,itc) are TC and derivatives wrt T and P
! phres%dgval(*,*,itc) are derivatives of TC wrt T, P and Y
! phres%d2gval(*,itc) are derivatives of TC wrt Y1 and Y2
! phres%gval(*,ibm) are beta and dervatives etc
! TC and beta must not depend on T, only on P and Y
!    dtaodt=one/tc
!    dtaodp=-tao/tc*phres%gval(3,itc)
! d2taodt2 is zero
   d2taodtdp=-one/tc*phres%gval(3,itc)
   d2taodp2=2.0d0*tao/tc**2*phres%gval(3,itc)-tao/tc*phres%gval(6,itc)
! 1-6 means F, F.T, T.P, F.T.T, F.T.P and F.P.P
   addgval(4)=2.0d0*rgasm*ftao(2)*dtaodt*logb1+&
        rt*ftao(4)*(dtaodt)**2*logb1
   addgval(5)=rgasm*ftao(2)*dtaodp*logb1+&
        rgasm*ftao(1)*invb1*phres%gval(3,ibm)+&
        rt*ftao(4)*dtaodt*dtaodp*logb1+&
        rt*ftao(2)*d2taodtdp*logb1+&
        rt*ftao(2)*dtaodt*invb1*phres%gval(3,ibm)
   addgval(6)=rt*ftao(4)*(dtaodp)**2*logb1+&
        rt*ftao(2)*d2taodp2*logb1+rt*ftao(1)*dtaodp*invb1*phres%gval(3,ibm)+&
        rt*ftao(2)*dtaodp*invb1*phres%gval(3,ibm)-&
        rt*ftao(1)*(invb1*phres%gval(3,ibm))**2+&
        rt*ftao(1)*invb1*phres%gval(6,ibm)
! G, G.T and G.Y, G.T.Y and G.Y1.Y2 correct (no P dependence checked)
   do j=1,mc
      dtao(1,j)=-tao*phres%dgval(1,j,itc)/tc
      dtao(2,j)=-phres%dgval(1,j,itc)/tc**2
      dtao(3,j)=2.0d0*tao*phres%gval(3,itc)*phres%dgval(1,j,itc)/tc**2-&
           tao*phres%dgval(3,j,itc)/tc
      do k=j,mc
         jxsym=kxsym(j,k)
         d2tao(jxsym)=&
              2.0*tao*phres%dgval(1,j,itc)*phres%dgval(1,k,itc)/tc**2&
              -tao*phres%d2gval(jxsym,itc)/tc
      enddo
   enddo
   do j=1,mc
! first derivative wrt Y, checked for bcc in Cr-Fe-Mo, error in fcc in c-cr-fe?
      daddgval(1,j)=rt*ftao(2)*dtao(1,j)*logb1+&
           rt*ftao(1)*invb1*phres%dgval(1,j,ibm)
!      write(*,43)j,daddgval(1,j),dtao(1,j),phres%dgval(1,j,ibm)
!43    format('3H Inden 4: ',i2,6(1pe12.5))
! second derivative wrt to T and Y, checked
      daddgval(2,j)=rgasm*ftao(2)*dtao(1,j)*logb1+&
           rgasm*ftao(1)*invb1*phres%dgval(1,j,ibm)+&
           rt*ftao(4)*dtaodt*dtao(1,j)*logb1+&
           rt*ftao(2)*dtao(2,j)*logb1+&
           rt*ftao(2)*dtaodt*invb1*phres%dgval(1,j,ibm)
!       write(*,56)rgasm*ftao(2)*dtao(1,j)*logb1,&
!            rgasm*ftao(1)*invb1*phres%dgval(1,j,ibm),&
!            rt*ftao(4)*dtaodt*dtao(1,j)*logb1,&
!            rgasm*ftao(2)*dtao(2,j)*logb1,&
!            rt*ftao(2)*dtaodt*invb1*phres%dgval(1,j,ibm)
!56 format('3H calcmag : ',5(1PE13.5))
! second derivative wrt P and Y, no P dependence
      daddgval(3,j)=rt*ftao(4)*dtaodp*dtao(1,j)*logb1+&
           rt*ftao(2)*dtao(3,j)*logb1+&
           rt*ftao(2)*dtao(1,j)*invb1*phres%gval(3,ibm)-&
           rt*ftao(1)*invb1**2*phres%gval(3,ibm)*phres%dgval(1,j,ibm)+&
           rt*ftao(1)*invb1*phres%dgval(3,j,ibm)
      do k=j,mc
! second derivatives wrt Y1 and Y2, wrong
         jxsym=kxsym(j,k)
         d2addgval(jxsym)=rt*ftao(4)*dtao(1,j)*dtao(1,k)*logb1+&
              rt*ftao(2)*d2tao(jxsym)*logb1+&
              rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm)+&
              rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm)-&
              rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm)+&
              rt*ftao(1)*invb1*phres%d2gval(jxsym,ibm)
!         d2addgval(ixsym(j,k))=rt*ftao(4)*dtao(1,j)*dtao(1,k)*logb1+&
!              rt*ftao(2)*d2tao(ixsym(j,k))*logb1+&
!              rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm)+&
!              rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm)-&
!              rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm)+&
!              rt*ftao(1)*invb1*phres%d2gval(ixsym(j,k),ibm)
!          write(*,57)rt*ftao(4)*dtao(1,j)*dtao(1,k)*logb1,&
!   R            rt*ftao(2)*d2tao(ixsym(j,k))*logb1,&
!               rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm),&
!               rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm),&
!              -rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm),&
!               rt*ftao(1)*invb1*phres%d2gval(ixsym(j,k),ibm)
!57 format('3H mag2y: ',6(1PE12.4))
      enddo
   enddo
! now add all to the total G and its derivatives
! NOTE if addpermole bit set we have to multiply with derivatives of
! the size of the phase ...
   if(btest(lokadd%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize magadd 1: ',lokph,addpermole,msize
! UNFINISHED: ignoring that msize depend on fractions
   else
      addpermole=.FALSE.; msize=one
   endif
   do j=1,mc
!      write(*,99)'3H magadd 1: ',1,j,phres%dgval(1,j,1),daddgval(1,j)/rt
      do k=1,3
! first derivatives
         phres%dgval(k,j,1)=phres%dgval(k,j,1)+msize*daddgval(k,j)/rt
      enddo
99    format(a,2i3,2(1pe16.8))
      do k=j,mc
! second derivatives
!         write(*,99)'3H magadd 2: ',k,j,rt*phres%d2gval(ixsym(j,k),1),&
!              d2addgval(ixsym(j,k))
         jxsym=kxsym(j,k)
         phres%d2gval(jxsym,1)=phres%d2gval(jxsym,1)+&
              msize*d2addgval(jxsym)/rt
!         phres%d2gval(ixsym(j,k),1)=phres%d2gval(ixsym(j,k),1)+&
!              msize*d2addgval(ixsym(j,k))/rt
      enddo
   enddo
!   write(*,*)'3H cm 7: ',phres%gval(1,1),addgval(1)/rt
! note phres%gval(1..3,1) already calculated above
   do j=4,6
      lokadd%propval(j)=msize*addgval(j)
      phres%gval(j,1)=phres%gval(j,1)+msize*addgval(j)/rt
   enddo
1000 continue
   return
 end subroutine calc_magnetic_inden

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_xiongmagnetic
!\begin{verbatim}
 subroutine create_xiongmagnetic(addrec,ibm,bcc)
! adds a Xiong type magnetic record, we must separate fcc and bcc by extra
! copied from Inden magnetic model
! The difference is that it uses CTA for Curie temperature and NTA for Neel
! and individual (IBM=.TRUE.) or average Bohr magneton numbers 
! BCC is .TRUE. if it is a BCC phase
   implicit none
   logical ibm,bcc
!   integer lokph
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim} %+
   integer typty,ip,nc,jj
   character text*128
   integer, parameter :: ncc=6
   double precision coeff(ncc),dval
   integer koder(5,ncc)
!   TYPE(tpfun_expression), pointer :: llow,lhigh
   TYPE(tpfun_expression) :: llow,lhigh
!
! from W Xiong et al Calphad (2012) 11-20
!
! G = RT g(tao) ln(b* + 1)
!
! tao = T/Tc
! b*= \Pi_i (b_i + 1)**(x_i) - 1
!
! tao<0 g(tao)=0
!
! tao<1 g(tao) = 1-1/D( 0.38438376/(p*tao)+0.63570895(1/p-1)*
!
!            (tao**3/6 + tao**9/135 + tao**15/600 + tao**21/1617) )
!
! tao>1 g(tao) = tao**(-7)/D( 1/21 + tao**(-14)/630 + tao**(-28)/2975 + 
!                             tao**(-42)/8232)
!
! p=0.37 for bcc and p=0.25 for non-bcc (like fcc)
!
! for bcc:      +1-.880323235*TAO**(-1)-.152870878*TAO**3-.00679426123*TAO**9
!               -.00152870878*TAO**15-5.67238878E-04*TAO**21
!
!              -.0403514888*TAO**(-7)-.00134504963*TAO**(-21)
!              -2.84834039E-04*TAO**(-35)-1.02937472E-04*TAO**(-49)
!
! for non-bcc: 
!
!      +1-.842849633*TAO**(-1)-.174242226*TAO**3-.00774409892*TAO**9
!      -.00174242226*TAO**15-6.46538871E-04*TAO**21
!
!       -.0261039233*TAO**(-7)-8.70130777E-04*TAO**(-21)
!      -1.84262988E-04*TAO**(-35)-6.65916411E-05*TAO**(-49)
!
!   write(*,*)'3H Qing-Xiong magnetic model',bcc
!
   if(bcc) then
! Magnetic function below Curie/Neel Temperature, 
! problem in ct1xfn to start a function with +1 or 1
      text=' +1-.880323235*T**(-1)-.152870878*T**3-.00679426123*T**9'//&
           '-.00152870878*T**15-5.67238878E-04*T**21'
! CHANGE OF REFERENCE STATE OF THE ELEMENTS
!      text=' +1-.152870878*T**3-.00679426123*T**9'//&
!           '-.00152870878*T**15-5.67238878E-04*T**21'
!      write(*,*)'3H emm 1: ',trim(text)
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
!       write(*,17)'3H emm 1B:',nc,(coeff(i),i=1,nc)
17     format(a,i3,5(1PE11.3))
      if(gx%bmperr.ne.0) goto 1000
      call ct1mexpr(nc,coeff,koder,llow)
      if(gx%bmperr.ne.0) goto 1000
! Magnetic function above Curie/Neel Temperature
      text='-.0403514888*T**(-7)-.00134504963*T**(-21)'//&
           '-2.84834039E-04*T**(-35)-1.02937472E-04*T**(-49)'
! CHANGE OF REFERENCE STATE OF THE ELEMENTS
!      text=' +.880323235*T**(-1)-.0403514888*T**(-7)-.00134504963*T**(-21)'//&
!           '-2.84834039E-04*T**(-35)-1.02937472E-04*T**(-49)'
!       write(*,*)'3H emm 2: ',trim(text)
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
      call ct1mexpr(nc,coeff,koder,lhigh)
      if(gx%bmperr.ne.0) goto 1000
! this is 1/(p*D) in eq. A9 in Qing et al, p=0.37
!   dval=0.880323235D0
      dval=one/(0.49649686D0+0.37D0*(0.33461979D0-0.49649686D0))
!      write(*,*)'3H added Qing-Xiong magnetic contribution to a bcc phase'
   else
!------------
! fcc
! Magnetic function below Curie/Neel Temperature
! REFERENCE STATE AT T=0
      text=' +1-.842849633*T**(-1)-.174242226*T**3-.00774409892*T**9'//&
           '-.00174242226*T**15-6.46538871E-04*T**21'
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
      call ct1mexpr(nc,coeff,koder,llow)
      if(gx%bmperr.ne.0) goto 1000
! Magnetic function above Curie/Neel Temperature
      text=' -.0261039233*T**(-7)'//&
           '-8.70130777E-04*T**(-21)-1.84262988E-04*T**(-35)'//&
           '-6.65916411E-05*T**(-49)'
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
      call ct1mexpr(nc,coeff,koder,lhigh)
      if(gx%bmperr.ne.0) goto 1000
! this is 1/(p*D) in eq. A9 in Qing et al for FCC, p=0.25; 1/p-1 =3.0
!      dval=4.0D0/(0.33461979D0+3.0D0*0.49649686D0)
!     dval=0.842849633D0
      dval=one/(0.49649686D0+0.25D0*(0.33461979D0-0.49649686D0))
!      write(*,*)'3H added Qing-Xiong magnetic contribution to a non-bcc phase'
   endif
! reserve an addition record
!   write(*,*)'3H 1/(pD)= ',dval
   allocate(addrec)
! store data in record
   allocate(addrec%explink(2))
   allocate(addrec%constants(1))
   nullify(addrec%nextadd)
   addrec%type=xiongmagnetic
! beware of segmentation fault here !!! llow and llhigh no longer pointers
   addrec%explink(1)=llow
   addrec%explink(2)=lhigh
   addrec%constants(1)=dval
   addrecs=addrecs+1
! Set bit 1 that there are properties
   addrec%status=0
   addrec%status=ibset(addrec%status,ADDHAVEPAR)
   if(bcc) addrec%status=ibset(addrec%status,ADDBCCMAG)
!   write(*,*)'3H Qing-Xiong magnetic addition: ',addrec%status,bcc,ADDBCCMAG
   allocate(addrec%need_property(3))
   addrec%addrecno=addrecs
! here the property list is searched for CTA, NTA and IBM
   call need_propertyid('CTA ',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(1)=typty
! The individual Bohr magneton number or not set in PMON
! WHEN READ FROM A UNFORMATTED FILE, DO WE KNOW LOKPH??
!   if(btest(phlista(lokph)%status1,PHBMAV)) then
   if(.not.ibm) then
! This model use an effective Bohr magneton number b*=prod(b_i+1)**x_i -1
      call need_propertyid('BMAG ',typty)
   else
! or an individual Bohr magneton number b*=prod(b_i+1)**x_i -1
!      write(*,*)'3H using induvidual Bohr magneton numbers',&
!           btest(phlista(lokph)%status1,PHBMAV)
      call need_propertyid('IBM ',typty)
   endif
!---------------------------------------------------
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(2)=typty
! NTA is not so important, anti-magnetic contributions usually small
   call need_propertyid('NTA ',typty)
   if(gx%bmperr.ne.0) then
      gx%bmperr=0
      addrec%need_property(3)=0
   else
      addrec%need_property(3)=typty
   endif
!   write(*,*)'3H need properties: ',(addrec%need_property(jj),jj=1,3)
1000 continue
   return
 end subroutine create_xiongmagnetic

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_xiongmagnetic
!\begin{verbatim}
 subroutine calc_xiongmagnetic(moded,phres,lokadd,lokph,mc,ceq)
! calculates Indens-Qing-Xiong magnetic contribution
! 
! Gmagn = RT*f(T/Tc)*ln(beta+1)
! moded: integer, 0=only G, S, Cp; 1=G and dG/dy; 2=Gm dG/dy and d2G/dy2
! phres: pointer, to phase\_varres record
! lokadd: pointer, to addition record
! lokph: integer, phase record 
! mc: integer, number of constituents
! ceq: pointer, to gtp_equilibrium_data
   implicit none
   integer moded,lokph,mc
! phres points to result record with gval etc for this phase
   TYPE(gtp_phase_varres) :: phres
   TYPE(gtp_phase_add), pointer :: lokadd
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer itc,itn,ibm,jl,noprop,ik,k,jk,j,jxsym,ip
   double precision logb1,invb1,iafftc,iaffbm,rgasm,rt,tao,gmagn,msize
   double precision dtaodt,dtaodp,beta,d2taodp2,d2taodtdp,tc,tv,plus,fixit
   double precision tao2(2),ftao(6),dtao(3,mc),d2tao(mc*(mc+1)/2)
   double precision addgval(6),daddgval(3,mc),d2addgval(mc*(mc+1)/2)
   double precision tn,tcsave,tnsave,gmdo_inf,dgmdo_infdt,d2gmdo_infdt2
   logical addpermole
   character line*128,tps(2)*3
   TYPE(tpfun_expression), pointer :: exprot
! dgdt = Gmagn/T + RT*df/dtao*dtao/dT*ln(beta+1)
! dgdp = RT df/dtao*dtao/dP*ln(beta+1)
! dgdy = RT*df/dtao*dtao/dy*ln(beta+1) + RT*f/(beta+1)*dbeta/dy
! d2gdt2=2*R*df/dtao*dtao/dT*ln(beta+1) + RT*d2f/dtao2*(dtao/dT)**2*ln(beta+1)
!        +RT*df/dtao*d2tao/dT2*ln(beta+1)
! d2gdtdp= ...
! d2gdp2=
! d2gdtdy= suck
! d2gdpdy=
! d2gdydy=
! listprop(1) is the number of properties calculated
! listprop(2:listprop(1)) give the typty of different properties
! calculated in gval(*,i) etc
! one has to find those with typty equal for need_property in the magnetic
! record, i.e. typty for CTA/NTA and typty for BMAG/IBM
! the properties needed.
!
   noprop=phres%listprop(1)-1
   itc=0; ibm=0; itn=0
   lokadd%propval=zero
!    write(*,*)'3H cmi 2: ',noprop,(phres%listprop(i),i=1,noprop)
! Inden-Qing-Xiong magnetic need properties in need_property(1..3)
   findix: do jl=2,noprop
      if(phres%listprop(jl).eq.lokadd%need_property(1)) then
         itc=jl
      elseif(phres%listprop(jl).eq.lokadd%need_property(2)) then
! we may also use an "average" Bohr magneton number in BMAG
         ibm=jl
      elseif(phres%listprop(jl).eq.lokadd%need_property(3)) then
         itn=jl
      endif
   enddo findix
! check that the needed properties are defined
!   write(*,*)'3H found magnetic properties: ',itc,ibm,itn
   if(ibm.eq.0 .or. (itc.eq.0 .and. itn.eq.0)) then
! it is no error if no CTA, NTA or BMAG but then magnetic contribution is zero
       write(*,12)trim(phlista(lokph)%name),itc,itn,ibm
12     format('3H *** Warning: Magnetic addition for phase ',a,&
            ' not calculated as '/&
            10x,'some values for CTA, NTA or BMAG/IBM are zero',3i3)
      goto 1000
   else
      tc=-one
      tn=-one
      if(itc.gt.0) tc=phres%gval(1,itc)
      if(itn.gt.0) tn=phres%gval(1,itn)
   endif
! I am not sure I calculate correct derivatives for indivudal Bohr magnetons ...
   if(.not.btest(phlista(lokph)%status1,PHBMAV)) then
      write(*,*)'3H *** Bohr magneton number derivatives not calculated'
   endif
   beta=phres%gval(1,ibm)
!   write(*,95)'3H Magnetic values in: ',itc,itn,ibm,tc,tn,beta
95 format(a,3i3,3(1PE15.6))
   if(beta.le.zero .or. (tc.le.zero .and. tn.le.zero)) then
! no magnetic contribution
      gmagn=zero
      addgval=zero
      daddgval=zero
      d2addgval=zero
      goto 1000
   endif
! we should use the appropriate tao=t/tc or tao=t/tn
! BUT WE MAY HAVE BOTH tc>0 and tn>0 !!
! use AF model unless tc negative, both cannot be negative here (test above)
! BUT WE MAY HAVE BOTH AS POSITIVE!!
   tcsave=tc
   if(tc.le.zero) then
! no ferro but maybe antiferro. One of them must be positive here!!
! Divide by AFF=3.0?
!      tc=tn/3.0D0
!      beta=beta/3.0D0
      tc=tn
! we use this index below to extract its value
      itc=itn
!   elseif(tn.gt.zero) then
! we have both AFM and FM, use FM, i.e. tc so nothing to do
   endif
!
   tv=ceq%tpval(1)
   rgasm=globaldata%rgas
   rt=rgasm*tv
   tao=tv/tc
   tao2(1)=tao
! one should save values of ftao if tao2 is the same next time ....
! but as tc depend on the constitution that is maybe not so often.
   if(tao.lt.one) then
      exprot=>lokadd%explink(1)
! VERY CLUMSY bug for debugging
   else
      exprot=>lokadd%explink(2)
   endif
!   plus=one
   plus=zero
! calculate function and derivatives wrt T, functions already created
   call ct1efn(exprot,tao2,ftao,ceq%eq_tpres)
! copied from list_addition
!   tps(1)='tao'
!   tps(2)='err'
!   ip=1
!   line=' '
!   call ct1wfn(exprot,tps,line,ip)
!   write(*,'(a,a,a/2(1pe12.4))')'f(tao)=',trim(line),';',tao,ftao(1)
!   write(*,'(a,2(1pe12.4))')'3H tao, f(tao): ',tao,ftao(1)
!   call wrice(kou,4,8,78,line(1:ip))
! the functions entered in explink use reference state at T=infinity
! correct for using reference state at T=0
! -1.0D0+0.38438376D0*lokadd%constants(1)*T/tc
! lokadd%constants(1) = 1/(p*D) in eq. A9 in paper by Qing
! NOTE tc may depend on P, we need the dtaodp and d2taodp2
   dtaodp=-tao/tc*phres%gval(3,itc)
   d2taodtdp=-one/tc*phres%gval(3,itc)
   d2taodp2=2.0d0*tao/tc**2*phres%gval(3,itc)-tao/tc*phres%gval(6,itc)
!
   fixit=0.38438376D0*lokadd%constants(1)
! this is for BCC
!   fixit=0.880323235D0
   fixit=zero
   ftao(1)=ftao(1)+plus*(fixit*tv/tc-one)
   ftao(2)=ftao(2)+plus*fixit/tc
   ftao(3)=ftao(3)-plus*fixit*tv*dtaodp/tc**2
! this is d2ftaodT2, no change
!   ftao(4)=ftao(4)
! this is d2ftaodTdP
   ftao(5)=ftao(5)-plus*fixit*dtaodp/tc**2
   ftao(6)=ftao(6)+2.0D0*plus*fixit*tv*d2taodp2/tc**3
!   if(plus.gt.zero) then
!      write(*,'(a,e17.9,a,e12.4)')'3H f(tao) correction: -1+',&
!           fixit,'/tao;',beta
!   else
!      write(*,'(a,e14.6,a,e14.6)')'3H f(tao) correction: +1-',&
!           fixit,'/tao;',beta
!   endif      
!------------------------------------------------------   
   logb1=log(beta+one)
   invb1=one/(beta+one)
   gmagn=rt*ftao(1)*logb1
!
! Calculate Gmdo_inf/RT, which value to use for "p"?
! THERE ARE T and composition derivatives of this also!!
!   gmdo_inf=-logb1*(one-0.38438376D0*lokadd%constants(1)/tao)
!   dgmdo_infdt=-logb1/tv
!   gmdo_inf=-rt*logb1*(one-0.38438376D0*lokadd%constants(1)/tc)
!   write(*,'(a,2(1pe14.6),a/5(1pe12.4))')'3H Gmdo(inf): ',&
!        rgasm*logb1*0.38438376D0*lokadd%constants(1)*tc,rgasm*logb1,'*T',&
!        rgasm,logb1,0.38438376D0,lokadd%constants(1),tc
!   gmdo_inf=-rgasm*logb1*(tv-0.38438376D0*lokadd%constants(1)*tc)
!   dgmdo_infdt=-rgasm*logb1
!   d2gmdo_infdt2=zero
!   write(*,88)'3H gmdo_inf: ',tv,gmdo_inf,dgmdo_infdt
!88 format(a,F8.2,4(1pe12.4))
!
!    write(*,98)'3H cm 97: ',tc,beta,ftao(1),logb1,rt
!    write(*,98)'3H cm 98: ',rt*gmagn,rt*(gmagn+phres%gval(1,1)),tcx,iafftc
!98  format(a,5(1PE14.6))
!
   dtaodt=one/tc
! d2taodT2=zero
! this already calculated above
!   dtaodp=-tao/tc*phres%gval(3,itc)
!   addgval(1)=gmagn+gmdo_inf
!   addgval(2)=gmagn/tv+rt*ftao(2)*dtaodt*logb1+dgmdo_infdt
!   addgval(3)=rt*ftao(2)*dtaodp*logb1+rt*ftao(1)*invb1*phres%gval(3,ibm)+&
!        d2gmdo_infdt2
   addgval(1)=gmagn
   addgval(2)=gmagn/tv+rt*ftao(2)*dtaodt*logb1
   addgval(3)=rt*ftao(2)*dtaodp*logb1+rt*ftao(1)*invb1*phres%gval(3,ibm)
! make sure d2G/dT2 is calculated and stored so it can be listed
   addgval(4)=(2.0D0*ftao(2)+tv*ftao(4)*dtaodt)*rgasm*dtaodt*logb1
!   phres%gval(1,1)=phres%gval(1,1)+addgval(1)/rt
!   phres%gval(2,1)=phres%gval(2,1)+addgval(2)/rt
!   phres%gval(3,1)=phres%gval(3,1)+addgval(3)/rt
! save these in record
! NOTE if parallel calculation the same stored values %propval will be
! written by all threads so they must not be used!!
!   write(*,77)lokadd%type,(lokadd%propval(j),j=1,4)
!77 format('3H addition ',i2,': ',4(1pe12.4))
!   if(moded.eq.0) then
   if(moded.eq.0) then
! They are included only for listing and debugging
      if(btest(lokadd%status,ADDPERMOL)) then
         addpermole=.TRUE.; msize=phres%abnorm(1)
      else
         msize=one
      endif
      do j=1,4
         lokadd%propval(j)=msize*addgval(j)
         phres%gval(j,1)=phres%gval(j,1)+msize*addgval(j)/rt
      enddo
! ignore second derivatives if no derivatives wanted
      goto 1000
   endif
! Now all derivatives with respect to fractions ...
! phres%gval(*,itc) are TC and derivatives wrt T and P
! phres%dgval(*,*,itc) are derivatives of TC wrt T, P and Y
! phres%d2gval(*,itc) are derivatives of TC wrt Y1 and Y2
! phres%gval(*,ibm) are beta and dervatives etc
! TC and beta must not depend on T, only on P and Y
!    dtaodt=one/tc
!    dtaodp=-tao/tc*phres%gval(3,itc)
! d2taodt2 is zero
!   d2taodtdp=-one/tc*phres%gval(3,itc)
!   d2taodp2=2.0d0*tao/tc**2*phres%gval(3,itc)-tao/tc*phres%gval(6,itc)
! 1-6 means F, F.T, T.P, F.T.T, F.T.P and F.P.P
!   addgval(4)=2.0d0*rgasm*ftao(2)*dtaodt*logb1+&
!        rt*ftao(4)*(dtaodt)**2*logb1
   addgval(5)=rgasm*ftao(2)*dtaodp*logb1+&
        rgasm*ftao(1)*invb1*phres%gval(3,ibm)+&
        rt*ftao(4)*dtaodt*dtaodp*logb1+&
        rt*ftao(2)*d2taodtdp*logb1+&
        rt*ftao(2)*dtaodt*invb1*phres%gval(3,ibm)
   addgval(6)=rt*ftao(4)*(dtaodp)**2*logb1+&
        rt*ftao(2)*d2taodp2*logb1+rt*ftao(1)*dtaodp*invb1*phres%gval(3,ibm)+&
        rt*ftao(2)*dtaodp*invb1*phres%gval(3,ibm)-&
        rt*ftao(1)*(invb1*phres%gval(3,ibm))**2+&
        rt*ftao(1)*invb1*phres%gval(6,ibm)
! G, G.T and G.Y, G.T.Y and G.Y1.Y2 correct (no P dependence checked)
   do j=1,mc
      dtao(1,j)=-tao*phres%dgval(1,j,itc)/tc
      dtao(2,j)=-phres%dgval(1,j,itc)/tc**2
      dtao(3,j)=2.0d0*tao*phres%gval(3,itc)*phres%dgval(1,j,itc)/tc**2-&
           tao*phres%dgval(3,j,itc)/tc
      do k=j,mc
         jxsym=kxsym(j,k)
         d2tao(jxsym)=&
              2.0*tao*phres%dgval(1,j,itc)*phres%dgval(1,k,itc)/tc**2&
              -tao*phres%d2gval(jxsym,itc)/tc
!         d2tao(ixsym(j,k))=&
!              2.0*tao*phres%dgval(1,j,itc)*phres%dgval(1,k,itc)/tc**2&
!              -tao*phres%d2gval(ixsym(j,k),itc)/tc
      enddo
   enddo
   do j=1,mc
! first derivative wrt Y, checked for bcc in Cr-Fe-Mo, error in fcc in c-cr-fe?
      daddgval(1,j)=rt*ftao(2)*dtao(1,j)*logb1+&
           rt*ftao(1)*invb1*phres%dgval(1,j,ibm)
!      write(*,43)j,daddgval(1,j),dtao(1,j),phres%dgval(1,j,ibm)
!43    format('3H Inden 4: ',i2,6(1pe12.5))
! second derivative wrt to T and Y, checked
      daddgval(2,j)=rgasm*ftao(2)*dtao(1,j)*logb1+&
           rgasm*ftao(1)*invb1*phres%dgval(1,j,ibm)+&
           rt*ftao(4)*dtaodt*dtao(1,j)*logb1+&
           rt*ftao(2)*dtao(2,j)*logb1+&
           rt*ftao(2)*dtaodt*invb1*phres%dgval(1,j,ibm)
!       write(*,56)rgasm*ftao(2)*dtao(1,j)*logb1,&
!            rgasm*ftao(1)*invb1*phres%dgval(1,j,ibm),&
!            rt*ftao(4)*dtaodt*dtao(1,j)*logb1,&
!            rgasm*ftao(2)*dtao(2,j)*logb1,&
!            rt*ftao(2)*dtaodt*invb1*phres%dgval(1,j,ibm)
!56 format('3H calcmag : ',5(1PE13.5))
! second derivative wrt P and Y, no P dependence
      daddgval(3,j)=rt*ftao(4)*dtaodp*dtao(1,j)*logb1+&
           rt*ftao(2)*dtao(3,j)*logb1+&
           rt*ftao(2)*dtao(1,j)*invb1*phres%gval(3,ibm)-&
           rt*ftao(1)*invb1**2*phres%gval(3,ibm)*phres%dgval(1,j,ibm)+&
           rt*ftao(1)*invb1*phres%dgval(3,j,ibm)
      do k=j,mc
! second derivatives wrt Y1 and Y2, wrong ??
         jxsym=kxsym(j,k)
         d2addgval(jxsym)=rt*ftao(4)*dtao(1,j)*dtao(1,k)*logb1+&
              rt*ftao(2)*d2tao(jxsym)*logb1+&
              rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm)+&
              rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm)-&
              rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm)+&
              rt*ftao(1)*invb1*phres%d2gval(jxsym,ibm)
!         d2addgval(ixsym(j,k))=rt*ftao(4)*dtao(1,j)*dtao(1,k)*logb1+&
!              rt*ftao(2)*d2tao(ixsym(j,k))*logb1+&
!              rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm)+&
!              rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm)-&
!              rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm)+&
!              rt*ftao(1)*invb1*phres%d2gval(ixsym(j,k),ibm)
!          write(*,57)rt*ftao(4)*dtao(1,j)*dtao(1,k)*logb1,&
!               rt*ftao(2)*d2tao(ixsym(j,k))*logb1,&
!               rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm),&
!               rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm),&
!              -rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm),&
!               rt*ftao(1)*invb1*phres%d2gval(ixsym(j,k),ibm)
!57 format('3H mag2y: ',6(1PE12.4))
      enddo
   enddo
! now add all to the total G
! NOTE if addpermole bit set we have to multiply with derivatives of
! the size of the phase ...
   if(btest(lokadd%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize magadd 2: ',lokph,addpermole,msize
   else
      addpermole=.FALSE.; msize=one
   endif
   do j=1,mc
      do k=1,3
!          write(*,99)'3H magadd 1: ',k,j,rt*phres%dgval(k,j,1),daddgval(k,j)
         phres%dgval(k,j,1)=phres%dgval(k,j,1)+msize*daddgval(k,j)/rt
      enddo
!99 format(a,2i3,2(1pe16.8))
      do k=j,mc
!          write(*,99)'3H magadd 2: ',k,j,rt*phres%d2gval(ixsym(j,k),1),&
!               d2addgval(ixsym(j,k))
         jxsym=kxsym(j,k)
         phres%d2gval(jxsym,1)=phres%d2gval(jxsym,1)+&
              msize*d2addgval(jxsym)/rt
!         phres%d2gval(ixsym(j,k),1)=phres%d2gval(ixsym(j,k),1)+&
!              msize*d2addgval(ixsym(j,k))/rt
      enddo
   enddo
!    write(*,*)'3H cm 7: ',rt*phres%gval(1,1),addgval(1)
! note phres%gval(1..3,1) already calculated above, multiplied with misize??
   do j=1,6
      lokadd%propval(j)=msize*addgval(j)
      phres%gval(j,1)=phres%gval(j,1)+msize*addgval(j)/rt
   enddo
! we may have destroyed the original value of tc if we have AFM
   tc=tcsave
!   write(*,900)tc,tn,tao,beta,phres%gval(4,1),lokadd%propval(4)
900 format('3H QX magn1: ',2F9.2,2F9.3,2(1pe12.4))
! jump here if no magnetic contribution
1000 continue
!   write(*,900)tc,tn,tao,beta,phres%gval(1,1),lokadd%propval(1)
   return
 end subroutine calc_xiongmagnetic

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_volmod1
!\begin{verbatim}
 subroutine create_volmod1(addrec)
! create addition record for the simple volume model
!
! currently only V = V0 * exp(VA(T))
! V0 is property (typty) 21, VA is 22 and reserved VB (Bulk modulus) 23 
! but VB is not implemented yet
   implicit none
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim} %+
   integer typty,kk
! reserve an addition record
   allocate(addrec)
! store data in record
   nullify(addrec%nextadd)
   addrec%status=0
   addrec%type=volmod1
! addrecs declared in gtp3.F90 but I am not sure it is needed or used
   addrecs=addrecs+1
   addrec%addrecno=addrecs
   allocate(addrec%need_property(3))
! properties needed
   call need_propertyid('V0  ',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(1)=typty
   call need_propertyid('VA  ',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(2)=typty
   call need_propertyid('VB  ',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(3)=typty
!   write(*,*)'Added volume model 1'
! store zero in 6 values for propval
   addrec%propval=zero
1000 continue
!   write(*,*)'3H created volume addition',addrecs
   return
 end subroutine create_volmod1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_volmod1
!\begin{verbatim} %-
 subroutine calc_volmod1(moded,phres,lokadd,lokph,mc,ceq)
! calculate the simple volume model, CURRENTLY IGNORING COMPOSITION DEPENDENCE
!
! G = P*V0(x)*exp(VA(T,x))
! moded: integer, 0=only G, S, Cp; 1=G and dG/dy; 2=Gm dG/dy and d2G/dy2
! phres: pointer, to phase\_varres record
! lokadd: pointer, to addition record
! lokph: integer, phase record index
! mc: integer, number of constituents
! ceq: pointer, to gtp_equilibrium_data
   implicit none
   integer moded,lokph,mc
! phres points to result record with gval etc for this phase
   TYPE(gtp_phase_varres) :: phres
   TYPE(gtp_phase_add), pointer :: lokadd
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer jl,iv0,iva,ivb,noprop
   double precision v0,va,vb,vol,deltap,pvol
! propval are stored locally in addition record
   lokadd%propval=zero
! if this bit not set there are no volume parameters
   if(.not.btest(lokadd%status,ADDHAVEPAR)) goto 1000
   iv0=0; iva=0; ivb=0;
   v0=zero; va=zero; vb=zero
   noprop=phres%listprop(1)-1
   findix: do jl=2,noprop
      if(phres%listprop(jl).eq.lokadd%need_property(1)) then
         iv0=jl
         v0=phres%gval(1,iv0)
      elseif(phres%listprop(jl).eq.lokadd%need_property(2)) then
         iva=jl
         va=phres%gval(1,iva)
      elseif(phres%listprop(jl).eq.lokadd%need_property(3)) then
         ivb=jl
         vb=phres%gval(1,ivb)
      endif
   enddo findix
!   write(*,'(a,3i3)')'3H volmodel: ',iv0,iva,ivb
! if iv0 is zero there are no volume data
   if(iv0.eq.0) goto 1000
! reference pressure is 1 bar
   deltap=ceq%tpval(2)-1.0D5
   if(ivb.ne.zero) then
! if ivb not zero there are bulk modulus data ... NOT IMPLEMENTED
      write(*,*)'3H Volume model with bulk modulus not implemented'
   else
      vol=v0/ceq%rtn
      if(iva.ne.zero) then
! NOTE all values should be divided by RT
         vol=v0*exp(va)/ceq%rtn
      endif
!      write(*,*)'3H v0, va: ',v0,va
      pvol=deltap*vol
! contribtions to G and derivatives, G, G.T, G.P=V, G.T.T, G.T.P, G.P.P
! NOTE there may be other parameters which depend on P!
      phres%gval(1,1)=phres%gval(1,1)+pvol
! G.T
      phres%gval(2,1)=phres%gval(2,1)+pvol*phres%gval(2,iva)
! G.P
      phres%gval(3,1)=phres%gval(3,1)+vol
! G.T.T
      phres%gval(4,1)=phres%gval(4,1)+&
           pvol*(phres%gval(2,iva)**2+phres%gval(4,iva))
! G.T.P
      phres%gval(5,1)=phres%gval(5,1)+vol*phres%gval(2,iva)
! G.P.P
! for the moment ignore pressure and composition dependence ...
!      phres%gval(6,1)=phres%gval(6,1)
   endif
!   write(*,*)'Calculated volume ',pvol,vol,deltap
!  store some property values
   lokadd%propval(1)=pvol
   lokadd%propval(2)=pvol*phres%gval(2,iva)
   lokadd%propval(3)=vol
   lokadd%propval(4)=pvol*(phres%gval(2,iva)**2+phres%gval(4,iva))
   lokadd%propval(5)=vol*phres%gval(2,iva)
   lokadd%propval(6)=zero
1000 continue
   return
 end subroutine calc_volmod1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_elastic_model_a
!\begin{verbatim}
 subroutine create_elastic_model_a(newadd)
! addition record to calculate the elastic energy contribution
   implicit none
   type(gtp_phase_add), pointer :: newadd
!\end{verbatim} %+
   integer typty
   allocate(newadd)
   newadd%type=elasticmodel1
   allocate(newadd%need_property(5))
   newadd%status=0
! needed properties
   newadd%need_property=0
   call need_propertyid('LPX ',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(1)=typty
   call need_propertyid('EC11',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(2)=typty
   call need_propertyid('EC12',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(3)=typty
   call need_propertyid('EC44',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(4)=typty
   call need_propertyid('LPTH',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(5)=typty
! now elastica is declared as pointer, is that OK?
   allocate(newadd%elastica)
1000 continue
   return
 end subroutine create_elastic_model_a

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_elastica
!\begin{verbatim}
 subroutine calc_elastica(moded,phres,addrec,lokph,mc,ceq)
! calculates elastic contribution and adds to G and derivatives
   implicit none
   integer moded,lokph,mc
   type(gtp_phase_varres), pointer :: phres
   type(gtp_phase_add), pointer :: addrec
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer jl,ilpx,ilpth,iec11,iec12,iec44,noprop,i1,i2
   double precision sum1,sum2
! get the current lattice parameters and elastic constants
   ilpx=0; ilpth=0; iec11=0; iec12=0; iec44=0
   noprop=phres%listprop(1)-1
   findix: do jl=2,noprop
      if(phres%listprop(jl).eq.addrec%need_property(1)) then
         ilpx=jl
      elseif(phres%listprop(jl).eq.addrec%need_property(2)) then
         iec11=jl
      elseif(phres%listprop(jl).eq.addrec%need_property(3)) then
         iec12=jl
      elseif(phres%listprop(jl).eq.addrec%need_property(4)) then
         iec44=jl
      elseif(phres%listprop(jl).eq.addrec%need_property(5)) then
! this one may not be needed initially at least
         ilpth=jl
      endif
   enddo findix
   if(ilpx.eq.0 .or. iec11.eq.0 .or. iec12.eq.0 .or. iec44.eq.0) then
      write(*,11)'3H Missing elastic parameter index: ',ilpx,iec11,iec12,iec44
11    format(a,5i4)
   endif
!   write(*,11)'3H indices: ',ilpx,iec11,iec12,iec44
! take care of the special elastic record
! ignore compsition derivatives at present ...
! elastic constant matrix, Voigt notation, symetric
   addrec%elastica%cmat=zero
   addrec%elastica%cmat(1,1)=phres%gval(1,iec11)
   addrec%elastica%cmat(2,2)=phres%gval(1,iec11)
   addrec%elastica%cmat(3,3)=phres%gval(1,iec11)
   addrec%elastica%cmat(4,4)=phres%gval(1,iec44)
   addrec%elastica%cmat(5,5)=phres%gval(1,iec44)
   addrec%elastica%cmat(6,6)=phres%gval(1,iec44)
   addrec%elastica%cmat(1,2)=phres%gval(1,iec12)
   addrec%elastica%cmat(1,3)=phres%gval(1,iec12)
   addrec%elastica%cmat(2,3)=phres%gval(1,iec12)
   addrec%elastica%cmat(2,1)=phres%gval(1,iec12)
   addrec%elastica%cmat(3,1)=phres%gval(1,iec12)
   addrec%elastica%cmat(3,2)=phres%gval(1,iec12)
!   write(*,22)phres%gval(1,iec11),phres%gval(1,iec12),phres%gval(1,iec44)
22 format('Elastic constants: ',3(1pe12.4))
!   write(*,19)(addrec%elastica%cmat(1,i1),i1=1,6)
!   write(*,19)(addrec%elastica%cmat(2,i1),i1=1,6)
!   write(*,19)(addrec%elastica%cmat(3,i1),i1=1,6)
!   write(*,19)(addrec%elastica%cmat(4,i1),i1=1,6)
!   write(*,19)(addrec%elastica%cmat(5,i1),i1=1,6)
!   write(*,19)(addrec%elastica%cmat(6,i1),i1=1,6)
19 format('3H CIJ: ',6(1pe12.4))
!....................
! equilibrium lattice constant (cubic, just diagonal)
   addrec%elastica%latticepar=zero
   addrec%elastica%latticepar(1,1)=phres%gval(1,ilpx)
   addrec%elastica%latticepar(2,2)=phres%gval(1,ilpx)
   addrec%elastica%latticepar(3,3)=phres%gval(1,ilpx)
!   write(*,23)'3H Lattice parameter: ',phres%gval(1,ilpx)
!....................
! The equilibrium lattice distances are in LPX (cubic lattice)
! The current lattice parameters are in ceq%phres%curlat(3,3)
! generate epsa, Voigt notation
!   write(*,23)'3H curlat 1: ',(phres%curlat(i1,1),i1=1,3)
!   write(*,23)'3H curlat 2: ',(phres%curlat(i1,2),i1=1,3)
!   write(*,23)'3H curlat 3: ',(phres%curlat(i1,3),i1=1,3)
23 format(a,3(1pe12.4))
   addrec%elastica%epsa(1)=(phres%curlat(1,1)-addrec%elastica%latticepar(1,1))&
        /addrec%elastica%latticepar(1,1)
   addrec%elastica%epsa(2)=(phres%curlat(2,2)-addrec%elastica%latticepar(2,2))&
        /addrec%elastica%latticepar(2,2)
   addrec%elastica%epsa(3)=(phres%curlat(3,3)-addrec%elastica%latticepar(3,3))&
        /addrec%elastica%latticepar(3,3)
! as addrec%elastica%latticepar(2,3) is zero for cubic use (1,1)
   addrec%elastica%epsa(4)=&
        (2*(phres%curlat(2,3)-addrec%elastica%latticepar(2,3)))&
        /addrec%elastica%latticepar(1,1)
   addrec%elastica%epsa(5)=&
        (2*(phres%curlat(1,3)-addrec%elastica%latticepar(1,3)))&
        /addrec%elastica%latticepar(1,1)
   addrec%elastica%epsa(6)=&
        (2*(phres%curlat(1,2)-addrec%elastica%latticepar(1,2)))&
        /addrec%elastica%latticepar(1,1)
!   write(*,25)'3H ev1: ',(addrec%elastica%epsa(i1),i1=1,6)
25 format(a,6(1pe12.4))
!....................
! calculate the elastic energy ... I do not know how to use F08 matrix mult
   sum1=zero
   do i1=1,6
      sum2=zero
      do i2=1,6
         sum2=sum2+addrec%elastica%cmat(i1,i2)*addrec%elastica%epsa(i2)
      enddo
!      write(*,23)'3H sum2: ',sum2
      sum1=sum1+addrec%elastica%epsa(i1)*sum2
   enddo
   addrec%elastica%eeadd(1)=5.0D-1*sum1
   write(*,30)'3H: Elastic energy: ',addrec%elastica%eeadd(1)
30 format(a,1pe15.7)
! TYPE gtp_elastic_modela
!    double precision, dimension(3,3) :: latticepar
! epsilon in Voigt notation
!    double precision, dimension(6) :: epsa
! elastic constant matrix in Voigt notation
!    double precision, dimension(6,6) :: cmat
! calculated elastic energy addition (with derivative to T and P?)
!    double precision, dimension(6) :: eeadd
! maybe more
! end TYPE gtp_elastic_modela
   
1000 continue
   return
 end subroutine calc_elastica

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_lattice_parameters
!\begin{verbatim}
 subroutine set_lattice_parameters(iph,ics,xxx,ceq)
! temporary way to set current lattice parameters for use with elastic model a
   implicit none
   integer iph,ics
   double precision, dimension(3,3) :: xxx
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokph,lokcs
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   ceq%phase_varres(lokcs)%curlat=xxx
!   write(*,*)'3H Phase+set: ',lokph,lokcs
!   write(*,23)'3H slp 1: ',(ceq%phase_varres(lokcs)%curlat(i1,1),i1=1,3)
!   write(*,23)'3H slp 2: ',(ceq%phase_varres(lokcs)%curlat(i1,2),i1=1,3)
!   write(*,23)'3H slp 3: ',(ceq%phase_varres(lokcs)%curlat(i1,3),i1=1,3)
23 format(a,3(1pe12.4))
1000 continue
   return
 end subroutine set_lattice_parameters

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_einsteincp
!\begin{verbatim}
 subroutine create_einsteincp(newadd)
   implicit none
   type(gtp_phase_add), pointer :: newadd
!\end{verbatim} %+
   integer, parameter :: ncc=6
   integer typty
!
! G/RT = 1.5*R*THET + 3*ln( 1 - exp( -THET/T ) ) 
! No need to use TPFUN
!
! gtp_phase_add has variables:
! integer :: type,addrecno,aff
! integer, allocatable :: need_property
! type(tpfun_expression), dimension, pointer :: explink
! type(gtp_phase_add), pointer :: nextadd   
! for spme additions one may create other records but they must have
! the variables type and nextadd
!------------------------------------------
   allocate(newadd)
! Both Einstein and Debye models use THET
   newadd%type=einsteincp
   newadd%status=0
!   call need_propertyid('THET',typty)
   call need_propertyid('LNTH',typty)
   if(gx%bmperr.ne.0) goto 1000
   allocate(newadd%need_property(1))
   newadd%need_property(1)=typty
   nullify(newadd%nextadd)
1000 continue
   return
 end subroutine create_einsteincp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_einsteincp
!\begin{verbatim}
 subroutine calc_einsteincp(moded,phres,addrec,lokph,mc,ceq)
! Calculate the contibution due to Einste Cp model for low T
! moded 0, 1 or 2
! phres all results
! addrec pointer to addition record
! lokph phase record
! mc number of variable fractions
! ceq equilibrum record
!
! G = 1.5*R*THET + 3*R*T*ln( 1 - exp( -THET/T ) ) 
! This is easier to handle inside the calc routine without TPFUN
!
   implicit none
   integer moded,lokph,mc
   type(gtp_phase_varres), pointer :: phres
   type(gtp_phase_add), pointer :: addrec
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ith,noprop,extreme,j1,j2
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1,fact
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2,theta
   double precision gein,dgeindt,d2geindt2,msize,theta,test
   double precision addphm(6)
   logical addpermole
   double precision, allocatable :: dthet(:),d2thet(:),dein(:),d2ein(:)
!
   noprop=phres%listprop(1)-1
!   write(*,*)'3H thet: ',phres%listprop(2),addrec%need_property(1)
   findix: do ith=2,noprop
      if(phres%listprop(ith).eq.addrec%need_property(1)) goto 100
   enddo findix
!   write(*,*)'3H No value of THET for phase ',trim(phlista(lokph)%name)
   write(*,*)'3H No value of LNTH for phase ',trim(phlista(lokph)%name)
   gx%bmperr=4336; goto 1000
100 continue
   if(phres%gval(1,ith).le.one) then
!      write(*,69)'3H Illegal THET for phase ',trim(phlista(lokph)%name),&
      write(*,69)'3H Illegal LNTH for phase ',trim(phlista(lokph)%name),&
           phres%gval(1,ith)
69    format(a,a,1pe12.4)
      gx%bmperr=4399; goto 1000
   endif
! NOTE the parameter value is ln(thera)! take the exponential!
! ln(thet) is in gval(1,ith), derivatives in dgval(*,ith,*) and d2gval(ith,*)
! G/RT = 1.5*THET/T + 3*LN(1-exp(-THET/T)) 
! NOTE ALL VALUES CALCULATED AS FOR G/RT
! kvot=theta/T
   if(phres%gval(1,ith).gt.1.0D2) then
!      write(*,*)'3H Probably wrong value of THET, parameter should be ln(THET)'
      write(*,*)'3H Probably wrong value of LNTH, parameter should be ln(THET)'
      write(*,*)'3H error in phase: ',trim(phlista(lokph)%name)
      gx%bmperr=4399; goto 1000
   endif
! The exp( ) because the parameter value is LN(THETA)   
   theta=exp(phres%gval(1,ith))
!   kvot=exp(phres%gval(1,ith))/ceq%tpval(1)
   kvot=theta/ceq%tpval(1)
!   write(*,*)'3H LN(THET): ',trim(phlista(lokph)%name),phres%gval(1,ith),theta
!   write(*,'(a,4(1pe11.3))')'3H dTH/dyi: ',(phres%dgval(1,j1,ith),j1=1,mc)
!   write(*,'(a,4(2i2,1pe11.3))')'3H d2TH/dyidyj: ',&
!        ((j1,j2,phres%d2gval(ixsym(j1,j2),ith),j2=j1,mc),j1=1,mc)
! We must convert all derivatives to real THET ?? no ??
   allocate(dthet(mc))
   allocate(d2thet(mc*(mc+1)/2))
   do j1=1,mc
      do j2=j1,mc
         d2thet(ixsym(j1,j2))=phres%d2gval(ixsym(j1,j2),ith)
      enddo
      dthet(j1)=phres%dgval(1,j1,ith)
   enddo
!   do j1=1,mc
!      do j2=j1,mc
!         d2thet(ixsym(j1,j2))=exp(phres%d2gval(ixsym(j1,j2),ith))
!      enddo
!      dthet(j1)=exp(phres%dgval(1,j1,ith))
!   enddo
!   write(*,'(a,4(1pe11.3))')'3H A dTH/dyi: ',(dthet(j1),j1=1,mc)
!   write(*,'(a,4(2i2,1pe11.3))')'3H A d2TH/dyidyj: ',&
!        ((j1,j2,d2thet(ixsym(j1,j2)),j2=j1,mc),j1=1,mc)
! simpler .... if it is correct??  
!  dTHETA/dy_i = d/dyi(exp(LNTH))= exp(LNTH)*d/dy1(LNTH) = THETA*dLNTH/dy1 ??
!   do j1=1,mc
!      do j2=j1,mc
!         d2thet(ixsym(j1,j2))=theta*phres%d2gval(ixsym(j1,j2),ith)
!      enddo
!      dthet(j1)=theta*phres%dgval(1,j1,ith)
!   enddo
!   write(*,'(a,4(1pe11.3))')'3H B dTH/dyi: ',(dthet(j1),j1=1,mc)
!   write(*,'(a,4(2i2,1pe11.3))')'3H B d2TH/dyidyj: ',&
!        ((j1,j2,d2thet(ixsym(j1,j2)),j2=j1,mc),j1=1,mc)
!   write(*,70)'3H phres: ',ceq%tpval(1),phres%gval(1,1),phres%gval(2,1),&
!        phres%gval(3,1),phres%gval(4,1),kvot
! we should be careful with numeric overflow, for small T or large T
! no risk for overflow for exp(-kvot)
   if(kvot.gt.1.0D2) then
! T is very small, kvot very large, exp(kvot) may cause overflow, 
! exp(-kvot) is very small, ln(1-exp(-kvot)) is close to zero
! exp(kvot) may cause overflow, kvot/(exp(kvot)-1)=
! kvot*exp(-kvot)/(1-exp(-kvot)) = (1-kvot+kvot**2/2-...)/(1-kvot/2+...) = 1
      extreme=-1
!      kvotexpkvotm1=one
      expmkvot=zero
      kvotexpkvotm1=zero
      ln1mexpkvot=zero
   elseif(kvot.lt.1.0D-2) then
! T is very big, kvot is very small, exp(-kvot) approch 1, 1-exp(-kvot)=kvot
! exp(-kvot) is close unity, ln(1-exp(-kvot))=ln(1-(1-kvot+kvot**2/2+...)) =
!            ln(kvot-kvot**2/2+...)=ln(kvot)
! exp(kvot) is close to unity: exp(kvot)-1 = kvot+kvot**2/2+ ...
      extreme=1
      expmkvot=exp(-kvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
      ln1mexpkvot=log(one-expmkvot)
   else
! range of T and kvot where value varies, take care of composition derivatives
      extreme=0
      expmkvot=exp(-kvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
      ln1mexpkvot=log(one-expmkvot)
!      ln1mexpkvot=log(exp(kvot)-one)
   endif
! kvot is +THETA/T; gein is integrated cp contribution to the Gibbs energy
! gein is Einsten contribution/RT
   gein=1.5D0*kvot+3.0D0*ln1mexpkvot
!   write(*,71)'3H Cp E1:',extreme,ceq%tpval(1),gein,ln1mexpkvot,expmkvot,&
!        kvotexpkvotm1
! first derivative wrt T taking care of overflow
   dgeindt=3.0D0*(ln1mexpkvot-kvotexpkvotm1)/ceq%tpval(1)
! second derivative wrt T
! This is d2G/dT**2/(RT) = -T**2/R*(Einstein Cp/RT) (or rather Cv/RT)
   if(extreme.eq.-1) then
! take care of overflow at low T, kvotexpkvotm1=expmkvot=0 set above
      d2geindt2=zero
   else
      d2geindt2=-3.0D0*kvotexpkvotm1**2/(expmkvot*ceq%tpval(1)**2)
   endif
!   NOTE if addpermole bit set we have to multiply with derivatives of
! the size of the phase ...
   if(btest(addrec%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize lowT: ',lokph,addpermole,msize
   else
      addpermole=.FALSE.; msize=one
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize lowT: ',lokph,addpermole,msize
   endif
! BEGIN NEW CODE ---------------------------------------------------------
! wrong G^Ein/RT = gein = 1.5*THETA/T + 3*LN(exp(THETA/T) - 1) 
! G^Ein/RT = gein = 1.5*THETA/T + 3*LN(1-exp(-THETA/T)) 
! where z=ln(THETA); THETA=exp(z); z depend on composition; kvot=THETA/T
! we have dz/dy_i etc in phres%dgval(1,i,ith); (note z does not depend on T)
!         d2z/dy_idy_j in phres%d2gval(ixsym(i,j))
!
! wrong G^Ein = RT*(G^Ein/RT) = 1.5*R*T*kvot + 3*R*T*ln(exp(kvot) - 1); 
!          kvot=THETA/T
! G^Ein = 1.5*R*THETA + 3*R*T*ln(1-exp(-THETA/T)); THETA=exp(z(y_i))
!
! dG^Ein/dy_i=1.5*R*dTHETA/dy_i+3*R*exp(-THETA/T)/(1-exp(-THETA/T))*dTHETA/dy_i
!     = (1.5+3*exp(-THETA/T)/(exp(-THETA/T)-1))*R*dTHETA/dy_i
! Composition derivative of the Einstein function is
! dEin/dy_i/RT = ((1.5+3*exp(-THETA/T)/(exp(-THETA/T)-1)/T)*dTHETA/dy_i
! dTHETA/dy_i = exp(z)*dz/y_i = THETA*dz/dy_i
! dEin/dy_i/RT = ((1.5+(3/T)*exp(-THETA/T)/(1-exp(-THETA/T)))*(THETA/T)*dz/dy_i
! dTHETA/dy_i = exp(z)*dz/y_i = THETA*dz/y_i
! REMEMBER kvot=THET/T; expmkvot=exp(-kvot); gasconstant R=globaldata%rgas
!
! This is composition derivatives of THET
!
! as I want NOTE expmkvot is exp(-kvot) !! NOW IT WORKS !!! SUCK
   fact=1.5D0*(one+expmkvot)/(one-expmkvot)*kvot
! the curve below better, correct shape ...
!   write(*,77)'3H Einstein dE/dy',lokph,1,theta,phres%dgval(1,1,ith),fact,&
!        phres%dgval(1,1,1),phres%dgval(1,1,1)+fact*phres%dgval(1,1,ith)
77 format(a,2i2,5(1pe12.4))
!   allocate(dein(mc))
!   allocate(d2ein(mc*(mc+1)/2))
! VERY MESSY CODING AND I DO NOT UNDERSTAND IT/BoS 2021.04.02
   do j1=1,mc
!     write(*,77)'3H Einstein dE/dy',lokph,j1,theta,phres%dgval(1,j1,ith),fact,&
!           phres%dgval(1,j1,1),phres%dgval(1,j1,1)+fact*phres%dgval(1,j1,ith)
! we must use dthet(j1)=THETA*dgval(1,j1,ith) and not dgval(1,j1,ith) !!
! old phres%dgval(1,j1,1)=msize*(phres%dgval(1,j1,1)+fact*phres%dgval(1,j1,ith))
      phres%dgval(1,j1,1)=phres%dgval(1,j1,1)+msize*fact*dthet(j1)
!      write(*,*)'3H second derivatives missing for Einstein, SUCK'
!      dein(j1)=msize*fact*dthet(j1)
      do j2=j1,mc
! d2Ein/dy1dy2 = (1.5R*theta+3*R*(exp(kvot)-1)**(-1))*d2theta/dy1dy2 -
!                 3*R*exp(kvot)/(T*(exp(kvot)-1))**2*dtheta/dy1*dtheta/dy2
         phres%d2gval(ixsym(j1,j2),1)=phres%d2gval(ixsym(j1,j2),1)+&
              msize*(fact*d2thet(ixsym(j1,j2))-&
              3.0d0*exp(kvot)*(ceq%tpval(1)*(exp(kvot)-one))**(-2)*&
              dthet(j1)*dthet(j2))
!         d2ein(ixsym(j1,j2))=msize*(fact*d2thet(ixsym(j1,j2))-&
!              3.0d0*exp(kvot)*(ceq%tpval(1)*(exp(kvot)-one))**(-2)*&
!              dthet(j1)*dthet(j2))
      enddo
   enddo
! listing ....
!   write(*,'(a,4(1pe11.3))')'3H dein: ',(dein(j1),j1=1,mc)
!   write(*,'(a,4(2i2,1pe11.3))')'3H d2ein: ',&
!        ((j1,j2,d2ein(ixsym(j1,j2)),j2=j1,mc),j1=1,mc)
!
! END NEW CODE ----------------------------------------------------------
! debug value of G
!   write(*,77)'3H Einstein ln(theta):',lokph,0,theta,gein,test,msize
! return the values in phres%gval(*,1)
   phres%gval(1,1)=phres%gval(1,1)+msize*gein
   phres%gval(2,1)=phres%gval(2,1)+msize*dgeindt
!   phres%gval(3,1)=phres%gval(3,1)
   phres%gval(4,1)=phres%gval(4,1)+msize*d2geindt2
!   phres%gval(5,1)=phres%gval(5,1)
!   phres%gval(6,1)=phres%gval(6,1)
   addrec%propval(1)=msize*gein
   addrec%propval(2)=msize*dgeindt
   addrec%propval(4)=msize*d2geindt2
!   write(*,70)'3H Cp E3: ',ceq%tpval(1),gein,dgeindt,d2geindt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!
! NOTE Missing implementation of derivatives wrt comp.dep of THET.
! the THET parameter cannot depend on T
!   write(*,*)'3H calc_einsteincp not including composition dependence of thet'
! addphm should be G^phys, dG^phys/dT, dG^phys/dP, d2G^phys/dT^2 etc
   addphm=zero
   addphm(1)=gein
   addphm(2)=dgeindt
   addphm(4)=d2geindt2
! correct for formula unit
   call add_size_derivatives(moded,phres,addphm,lokph,mc,ceq)
!
1000 continue
   return
 end subroutine calc_einsteincp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine add_size_derivatives
!\begin{verbatim}
 subroutine add_size_derivatives(moded,phres,addphm,lokph,mc,ceq)
! Many physical models are defined per mole of atoms, as the Gibbs energy
! is calculate per mole formula unit this routine will handle the
! additional derivatives needed when M*ADD(1 mole)
! mc is number of constituent variables
! addphm(1..6) is G, dG/dt, dG/dp, d2G/dt2, d2G/dtdp, d2G/dp2 for the addition
   implicit none
   integer moded,lokph,mc
   type(gtp_phase_varres), pointer :: phres
!   type(gtp_phase_add), pointer :: addrec
   type(gtp_equilibrium_data), pointer :: ceq
   double precision addphm(6)
!\end{verbatim} %+
   integer i1,i2,j1,j2,jxsym,k,s1,s2
   double precision site1,site2
! Moles of constituents per mole formula units is:
! M = \sum_s a_s \sum_i y_si; dM/dy = a_s
! what about disordered fraction sets? ignore .... UNFINISHED??
! sites are in phres%sites
! number of constituents in siblattice s is in phlista(lokph)%nooffr(s)
   goto 1000
   write(*,*)'3H inside add_size_derivatives',lokph,addphm(1)
   s1=1
   s2=1
   site1=phres%sites(s1)
   j1=1
   do i1=1,mc
!
! G^phy_M = N*G^phy_m (already done)
! dG^phy_M/dyi = (dN/dyi)*G^phy_m  + N*(dG^phy_m/dyi) ??
! d2G/dyidyj = (dN/dyi)*(dGm/dyj)+(dN/dyj)*(dGm/dyi)+N*(d2Gm/dyidyj) ignore
! dN/dyi = a_s for sublattice s with constituent i 
! and d2N/dyidyj = 0
!
! I am not sure about this routine ....
!
! dGM/dyi = (dN/dyi)*Gm  (+N*(dGM/dyi) already done)
      do k=1,3
! this is dG/dy and d2G/dtdy and d2G/dpdy (k=1,2,3)
         phres%dgval(k,i1,1)=phres%dgval(k,i1,1)+site1*addphm(k)
      enddo
      write(*,*)'3H addition to mu',i1,site1*addphm(1)
!      write(*,*)'3H ignoring 2nd derivatives of size'
!+      j2=1
!+      site2=phres%sites(s2)
!+      do i2=i1,mc
! For the moment ignore all second derivatives !!
! no second derivatives wrt same constituents twice
!+        if(i2.gt.i1) then
!+            site2=phres%sites(s2)
!+            jxsym=kxsym(i1,i2)
! d2G/dyidyj = (dN/dyi)*(dGm/dyj)+(dN/dyj)*(dGm/dyi) (+N*(d2Gm/dyidyj) done)
!+            phres%d2gval(jxsym,1)=phres%d2gval(jxsym,1)+site1*site2*addphm(1)
!+         endif
!+         j2=j2+1
!+         if(j2.gt.phlista(lokph)%nooffr(s2)) then
!+            j2=1; s2=s2+1
!+            if(s2.le.phlista(lokph)%noofsubl) site2=phres%sites(s2)
!+         endif
!+      enddo
      j1=j1+1
      if(j1.gt.phlista(lokph)%nooffr(s1)) then
         j1=1; s1=s1+1
         if(s1.le.phlista(lokph)%noofsubl) site1=phres%sites(s1)
      endif
   enddo
1000 continue
   return
 end subroutine add_size_derivatives

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_schottky_anomaly
!\begin{verbatim}
 subroutine create_schottky_anomaly(newadd)
! Adding a Schottky anomaly to Cp
   implicit none
   type(gtp_phase_add), pointer :: newadd
!\end{verbatim} %+
   integer, parameter :: ncc=6
   integer typty
!
! G/RT =  SAM * ln( 1 + exp( SAM/T ) ) 
! No need to use TPFUN
!
! gtp_phase_add has variables:
! integer :: type,addrecno,aff
! integer, allocatable :: need_property
! type(tpfun_expression), dimension, pointer :: explink
! type(gtp_phase_add), pointer :: nextadd   
! for spme additions one may create other records but they must have
! the variables type and nextadd
!------------------------------------------
   allocate(newadd)
! Schottky anomaly uses THT2 and DCP2, same as second Einstein
   newadd%status=0
   newadd%type=schottkyanomaly
   allocate(newadd%need_property(2))
   call need_propertyid('TSCH',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(1)=typty
   call need_propertyid('CSCH',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(2)=typty
   nullify(newadd%nextadd)
1000 continue
   return
 end subroutine create_schottky_anomaly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_schottky_anomaly
!\begin{verbatim}
 subroutine calc_schottky_anomaly(moded,phres,addrec,lokph,mc,ceq)
! Calculate the contibution due to a Schottky anomaly
! moded 0, 1 or 2
! phres all results
! addrec pointer to addition record
! lokph phase record
! mc number of variable fractions
! ceq equilibrum record
!
! G = DCP2*T*ln( 1 + exp( -THT2/T ) ) 
! dG/dT = DCP2*(ln(1+exp(-THT2/T))+(THT/T)*(1+exp(+THT2/T))**(-1)
! d2G/dT2 = -DCP2*THT2**2*T**(-3)*exp(THT2/T)*(1+exp(+THT2/T))**(-2)
!
   implicit none
   integer moded,lokph,mc
   type(gtp_phase_varres), pointer :: phres
   type(gtp_phase_add), pointer :: addrec
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ith,jth,noprop,extreme,j1
   double precision kvot,expkvot,expmkvot,ln1pexpmkvot,kvotexpkvotp1,fact
   double precision gsch,dgschdt,d2gschdt2,dcp2,msize
   logical addpermole
!
   noprop=phres%listprop(1)-1
!   write(*,*)'3H lnth: ',phres%listprop(2),addrec%need_property(1)
   ith=0
   jth=0
   findix: do j1=2,noprop
      if(phres%listprop(j1).eq.addrec%need_property(1)) then
         ith=j1
      elseif(phres%listprop(j1).eq.addrec%need_property(2)) then
         jth=j1
      endif
   enddo findix
! ith is THT2 and jth is DCP2
   if(ith.eq.0 .or. jth.eq.0) then
!      write(*,*)'3H missing Schottky anomaly parameter for phase ',&
!           trim(phlista(lokph)%name)
      goto 1000
   endif
! phres%gval(1,ith) and phres(1,jth) must not depend on T
   kvot=exp(phres%gval(1,ith))/ceq%tpval(1)
   dcp2=phres%gval(1,jth)
   if(kvot.le.zero) goto 1000
! we should be careful with numeric overflow, for small T or large T
   if(kvot.gt.1.0D2) then
! T is very small, kvot very large, exp(kvot) may cause overflow, 
! exp(-kvot) is very small, ln(1-exp(-kvot)) is close to zero
! exp(kvot) may cause overflow, kvot/(exp(kvot)-1)=
! kvot*exp(-kvot)/(1-exp(-kvot)) = (1-kvot+kvot**2/2-...)/(1-kvot/2+...) = 1
      extreme=-1
      expmkvot=zero
      kvotexpkvotp1=zero
      ln1pexpmkvot=zero
   elseif(kvot.lt.1.0D-2) then
! T is very big, kvot is very small, exp(-kvot) approch 1, 1-exp(-kvot)=kvot
! exp(-kvot) is close unity, ln(1-exp(-kvot))=ln(1-(1-kvot+kvot**2/2+...)) =
!            ln(kvot-kvot**2/2+...)=ln(kvot)
! exp(kvot) is close to unity: exp(kvot)-1 = kvot+kvot**2/2+ ...
      extreme=1
      expmkvot=exp(-kvot)
      kvotexpkvotp1=kvot/(exp(kvot)+one)
      ln1pexpmkvot=log(one+expmkvot)
   else
! normal range of T and kvot
      extreme=0
      expmkvot=exp(-kvot)
      kvotexpkvotp1=kvot/(exp(kvot)+one)
      ln1pexpmkvot=log(one+expmkvot)
   endif
! 
! Note this is the G/RT value dcp2*ln(1+exp(tht2/T)
! G = DCP2*T*ln( 1 + exp( -THT2/T ) ) 
   gsch=dcp2*ln1pexpmkvot
!   write(*,71)'3H Cp E1:',extreme,ceq%tpval(1),gein,ln1mexpkvot,expmkvot,&
!        kvotexpkvotm1
! first derivative wrt T taking care of overflow
! dcp2*ln(1+exp(tht2))/T -(tht2/T**2)/exp
! dG/dT = DCP2*(ln(1+exp(-THT2/T))+(THT/T)*(1+exp(+THT2/T))**(-1)
   dgschdt=DCP2*(ln1pexpmkvot-kvotexpkvotp1)/ceq%tpval(1)
! This is d2G/dT**2/(RT) = -T**2/R*(Einstein Cp/RT) (or rather Cv/RT)
   if(extreme.eq.-1) then
! take care of overflow at low T, kvotexpkvotm1=expmkvot=0 set above
      d2gschdt2=zero
   else
! d2G/dT2 = -DCP2*THT2**2*T**(-3)*exp(THT2/T)*(1+exp(+THT2/T))**(-2)
      d2gschdt2=-DCP2*kvotexpkvotp1**2/(expmkvot*ceq%tpval(1))
   endif
! first derivative for each constituent. The parameter value is ln(theta)
! and we should divide by RT
!   fact=1.5D0*kvot+3.0D0*kvotexpkvotm1
!   do j1=1,mc
!      phres%dgval(1,j1,1)=phres%dgval(1,j1,1)+fact*phres%dgval(1,j1,ith)
!   enddo
! return the values in phres%gval(*,1)
! NOTE if addpermole bit set we have to multiply with derivatives of
! the size of the phase ...
   if(btest(addrec%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize schky: ',lokph,addpermole,msize
   else
      addpermole=.FALSE.; msize=one
   endif
!
   phres%gval(1,1)=phres%gval(1,1)+msize*gsch
   phres%gval(2,1)=phres%gval(2,1)+msize*dgschdt
!   phres%gval(3,1)=phres%gval(3,1)
   phres%gval(4,1)=phres%gval(4,1)+msize*d2gschdt2
!   phres%gval(5,1)=phres%gval(5,1)
!   phres%gval(6,1)=phres%gval(6,1)
   addrec%propval(1)=msize*gsch
   addrec%propval(2)=msize*dgschdt
   addrec%propval(4)=msize*d2gschdt2
!   write(*,70)'3H Schottky: ',ceq%tpval(1),gsch,dgschdt,d2gschdt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!
! Missing implem of derivatives wrt comp.dep of thet.  thet2 cannot depend on T
!
1000 continue
   return
 end subroutine calc_schottky_anomaly

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_secondeinstein
!\begin{verbatim}
 subroutine create_secondeinstein(newadd)
   implicit none
   type(gtp_phase_add), pointer :: newadd
!\end{verbatim} %+
   integer, parameter :: ncc=6
   integer typty
!
! G/RT = DCP2*ln( 1 - exp( -THT2/T ) ) 
! No need to use TPFUN
!
! gtp_phase_add has variables:
! integer :: type,addrecno,aff
! integer, allocatable :: need_property
! type(tpfun_expression), dimension, pointer :: explink
! type(gtp_phase_add), pointer :: nextadd   
! for spme additions one may create other records but they must have
! the variables type and nextadd
!------------------------------------------
   allocate(newadd)
   newadd%type=secondeinstein
   newadd%status=0
! The second Einstein use THT2 and DCP2
   allocate(newadd%need_property(2))
   call need_propertyid('THT2',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(1)=typty
   call need_propertyid('DCP2',typty)
   if(gx%bmperr.ne.0) goto 1000
   newadd%need_property(2)=typty
   nullify(newadd%nextadd)
   write(*,*)'3H created 2nd Einstein: ',newadd%type
1000 continue
   return
 end subroutine create_secondeinstein

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_secondeinstein
!\begin{verbatim}
 subroutine calc_secondeinstein(moded,phres,addrec,lokph,mc,ceq)
! Calculate the contibution due to Einste Cp model for low T
! moded 0, 1 or 2
! phres all results
! addrec pointer to addition record
! lokph phase record
! mc number of variable fractions
! ceq equilibrum record
!
! G = 1.5*R*THET + 3*R*T*ln( 1 - exp( -THET/T ) ) 
! This is easier to handle inside the calc routine without TPFUN
!
   implicit none
   integer moded,lokph,mc
   type(gtp_phase_varres), pointer :: phres
   type(gtp_phase_add), pointer :: addrec
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ith,jth,noprop,extreme,j1
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1,fact
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2
   double precision gein,dgeindt,d2geindt2,deltacp,msize
   logical addpermole
!
   noprop=phres%listprop(1)-1
!   write(*,*)'3H tht2: ',phres%listprop(2),addrec%need_property(1),&
!        addrec%need_property(2)
   ith=0; jth=0;
   findix: do j1=2,noprop
      if(phres%listprop(j1).eq.addrec%need_property(1)) then
         ith=j1
      elseif(phres%listprop(j1).eq.addrec%need_property(2)) then
         jth=j1
      endif
   enddo findix
   if(ith.eq.0 .or. jth.eq.0) then
      write(*,*)'3H Missing second Einstein properties for phase ',&
           trim(phlista(lokph)%name)
      gx%bmperr=4336; goto 1000
   endif
100 continue
   if(phres%gval(1,ith).le.one) then
      write(*,70)'3H Illegal LNTH for phase ',trim(phlista(lokph)%name),&
           phres%gval(1,ith)
      gx%bmperr=4399; goto 1000
   endif
! NOTE the parameter value is ln(thera)! take the exponential!
! thet is in gval(1,ith), derivatives in dgval(*,ith,*) and d2gval(ith,*)
! G/RT = phres%gval(1,jth)*R*LN(exp(THET/T) - 1) 
! NOTE ALL VALUES CALCULATED AS FOR G/RT
   kvot=exp(phres%gval(1,ith))/ceq%tpval(1)
!   write(*,70)'3H phres: ',ceq%tpval(1),phres%gval(1,1),phres%gval(2,1),&
!        phres%gval(3,1),phres%gval(4,1),kvot
! we should be careful with numeric overflow, for small T or large T
! no risk for overflow for exp(-kvot)
   if(kvot.gt.1.0D2) then
! T is very small, kvot very large, exp(kvot) may cause overflow, 
! exp(-kvot) is very small, ln(1-exp(-kvot)) is close to zero
! exp(kvot) may cause overflow, kvot/(exp(kvot)-1)=
! kvot*exp(-kvot)/(1-exp(-kvot)) = (1-kvot+kvot**2/2-...)/(1-kvot/2+...) = 1
      extreme=-1
!      kvotexpkvotm1=one
      expmkvot=zero
      kvotexpkvotm1=zero
      ln1mexpkvot=zero
   elseif(kvot.lt.1.0D-2) then
! T is very big, kvot is very small, exp(-kvot) approch 1, 1-exp(-kvot)=kvot
! exp(-kvot) is close unity, ln(1-exp(-kvot))=ln(1-(1-kvot+kvot**2/2+...)) =
!            ln(kvot-kvot**2/2+...)=ln(kvot)
! exp(kvot) is close to unity: exp(kvot)-1 = kvot+kvot**2/2+ ...
      extreme=1
      expmkvot=exp(-kvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
      ln1mexpkvot=log(one-expmkvot)
   else
! normal range of T and kvot
      extreme=0
      expmkvot=exp(-kvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
      ln1mexpkvot=log(one-expmkvot)
   endif
! 
! The Delta Cp value is given in phres%gval(1,jth)  It can be negative!
! and it can depend on P and composition !! NOT IMPLEMENTED !! BEWHERE
! In normal Einstein deltacp=3.0
   deltacp=phres%gval(1,jth)
   gein=deltacp*ln1mexpkvot
!   write(*,71)'3H Cp E1:',extreme,ceq%tpval(1),gein,ln1mexpkvot,expmkvot,&
!        kvotexpkvotm1
! first derivative wrt T taking care of overflow
   dgeindt=deltacp*(ln1mexpkvot-kvotexpkvotm1)/ceq%tpval(1)
! This is d2G/dT**2/(RT) = -T**2/R*(Einstein Cp/RT) (or rather Cv/RT)
   if(extreme.eq.-1) then
! take care of overflow at low T, kvotexpkvotm1=expmkvot=0 set above
      d2geindt2=zero
   else
      d2geindt2=-deltacp*kvotexpkvotm1**2/(expmkvot*ceq%tpval(1)**2)
   endif
!   write(*,16)'3H 2nd Einstein: ',kvot,deltacp,d2geindt2
16 format(a,6(1pe12.4))
! check if addition is per mole 
   if(btest(addrec%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize 2ndein: ',lokph,addpermole,msize
   else
      addpermole=.FALSE.; msize=one
   endif
! first derivative for each constituent. The parameter value is ln(theta)
! and we should divide by RT
   fact=deltacp*kvotexpkvotm1
   do j1=1,mc
      phres%dgval(1,j1,1)=phres%dgval(1,j1,1)+fact*phres%dgval(1,j1,ith)
   enddo
! return the values in phres%gval(*,1)
   phres%gval(1,1)=phres%gval(1,1)+msize*gein
   phres%gval(2,1)=phres%gval(2,1)+msize*dgeindt
!   phres%gval(3,1)=phres%gval(3,1)
   phres%gval(4,1)=phres%gval(4,1)+msize*d2geindt2
!   phres%gval(5,1)=phres%gval(5,1)
!   phres%gval(6,1)=phres%gval(6,1)
   addrec%propval(1)=msize*gein
   addrec%propval(2)=msize*dgeindt
   addrec%propval(4)=msize*d2geindt2
!   write(*,70)'3H Cp E3: ',ceq%tpval(1),gein,dgeindt,d2geindt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!
! Missing implem of derivatives wrt comp.dep of tht2 and dcp2.
! Neither tht2 nor dcp2 can depend on T
!
1000 continue
   return
 end subroutine calc_secondeinstein

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_twostate_model1
!\begin{verbatim}
 subroutine create_twostate_model1(addrec)
! newadd is location where pointer to new addition record should be stored
   implicit none
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim} %+
   integer typty
! this is bad programming as it cannot be deallocated but it will never be ...
! maybe pointers can be deallocated?
   allocate(addrec)
   addrec%status=0
! nullify pointer to next addition
   nullify(addrec%nextadd)
!-----------------------------
! The model consists of two contributions
! The first is the harmonic vibrations of an ideal amprthous phase
!     this requires a THETA representing the Einstein T
! The second is a term - RT*(1+exp(G2/RT))
! which represent the change from "solid like" to "liquid like"
!-----------------------------
! I am not sure what this is used for
   addrecs=addrecs+1
   addrec%addrecno=addrecs
! property needed
   allocate(addrec%need_property(2))
   call need_propertyid('G2  ',typty)
   addrec%need_property(1)=typty
!   call need_propertyid('THET ',typty)
   call need_propertyid('LNTH  ',typty)
   addrec%need_property(2)=typty
! type of addition
   addrec%type=twostatemodel1
! store zero.  Used to extract current value of this property
   addrec%propval=zero
1000 continue
!   write(*,*)'Created two state liquid record'
   return
 end subroutine create_twostate_model1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_twostate_model1
!\begin{verbatim}
 subroutine create_newtwostate_model1(addrec)
! newadd is location where pointer to new addition record should be stored
   implicit none
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim} %+
   integer typty
! this is bad programming as it cannot be deallocated but it will never be ...
! maybe pointers can be deallocated?
   allocate(addrec)
! nullify pointer to next addition
   nullify(addrec%nextadd)
   addrec%status=0
!-----------------------------
! The model consists of two contributions
! The first is the harmonic vibrations of an ideal amprthous phase
!     this requires a THETA representing the Einstein T
! The second is a term - RT*(1+exp(G2/RT))
! which represent the change from "solid like" to "liquid like"
!-----------------------------
! I am not sure what this is used for
   addrecs=addrecs+1
   addrec%addrecno=addrecs
! property needed G2 is not needed as composition independent
   allocate(addrec%need_property(1))
!   call need_propertyid('G2  ',typty)
!   addrec%need_property(1)=typty
!   call need_propertyid('THETA  ',typty)
   call need_propertyid('LNTH ',typty)
   addrec%need_property(1)=typty
! type of addition  this is 12
   addrec%type=twostatemodel2
! store zero.  Used to extract current value of this property
   addrec%propval=zero
1000 continue
!   write(*,*)'Created two state liquid record ',twostatemodel2
   return
 end subroutine create_newtwostate_model1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_twostate_model_john
!\begin{verbatim}
 subroutine calc_twostate_model_john(moded,phres,addrec,lokph,mc,ceq)
! subroutine calc_twostate_model1(moded,phres,addrec,lokph,mc,ceq)
! CURRENTLY NOT USED
! This routine works OK but I am testing a modification
! moded is 0, 1 or 2 if no, first or 2nd order derivatives should be calculated
! addrec is addition record
! phres is phase_varres record
! lokph is phase location
! mc is number of constitution fractions
! ceq is current equilibrium record
   implicit none
   integer moded,lokph,mc
   TYPE(GTP_PHASE_ADD), pointer :: addrec
   TYPE(GTP_PHASE_VARRES), pointer :: phres
   TYPE(GTP_EQUILIBRIUM_DATA), pointer :: ceq
!\end{verbatim} %+
! two state model for extrapolating liquid to low T
! DG = d(H-RT) + RT( dln(d)+(1-d)ln(1-d))
! where d is "liquid like" atoms.  H is enthalpy to form defects
! At equilibrium
!
! d = exp(-H/RT) / (1 + e(-H/RT) ) is the integrated Einstein Cp -H/R is THET
!
! G^liq - G^amorph = G^amorph - RT ln(1+exp(-DG_d/RT)
! DG_d is the enthalpy of forming 1 mole of defects in the glassy state
!
!------------------------------
! The value of Gd for the phase is calculated and added to G
   integer jj,noprop,ig2,ith,extreme,jth,kth
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2
   double precision gein,dgeindt,d2geindt2
   double precision g2ein,dg2eindt,d2g2eindt2,theta2,dcpl
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1
   double precision g2val,dg2,expg2,expmg2,rt,tv,rg,dg2dt,dgfdt,d2g2dt2
   double precision expmg2p1,msize
   logical addpermole
! This is Johns original model
   write(*,*)'3H THIS VERSION OF TWOSTATE MODEL SHOULD NOT BE USED!'
   stop
! number of properties calculatied
   noprop=phres%listprop(1)-1
! locate the THET and G2 property record 
   ig2=0
   ith=0
   jth=0
   kth=0
! check if addition is per mole 
   if(btest(addrec%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize john: ',lokph,addpermole,msize
   else
      addpermole=.FALSE.; msize=one
   endif
   findix: do jj=2,noprop
      if(phres%listprop(jj).eq.addrec%need_property(1)) then
! current values of G2 is stored in phres%gval(1,ig2)
         ig2=jj;
      elseif(phres%listprop(jj).eq.addrec%need_property(2)) then
! current value of THET are stored in phres%gval(1,ith)
         ith=jj
! SECOND EINSTEIN CP CONTRBUTION ADDED SEPARATELY
!      elseif(phres%listprop(jj).eq.14) then
! current value of LIQUID THET are stored in with index 14 VISC
!         jth=jj
!         theta2=exp(phres%gval(1,jth))
!         write(*,*)'3H found liquid theta: ',theta2
!      elseif(phres%listprop(jj).eq.27) then
! current value of LIQUID THET are stored in with index 14 VISC
!         kth=jj
!         dcpl=phres%gval(1,kth)
!         write(*,*)'3H found liquid delta-cp: ',dcpl
      endif
   enddo findix
   if(ith.eq.0) then
!      write(*,*)'3H Cannot find value for amorphous THET'
      write(*,*)'3H warning no values for amorphous LNTH'
      gein=zero; dgeindt=zero; d2geindt2=zero
      goto 300
!      gx%bmperr=4367; goto 1000
   endif
!----------------------------------
! for the moment the composition dependence is ignored
!   write(*,19)'3H 2no1: ',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
!------ this THET part copied from calc_einstein
! thet is in gval(1,ith), derivatives in dgval(*,ith,*) and d2gval(ith,*)
! G/RT = 1.5*THET/T + 3*R*LN(exp(THET/T) - 1) 
! NOTE ALL VALUES CALCULATED AS FOR G/RT
! kvot=theta/T
! NOTE the stored value is ln(theta)! !!!
   kvot=exp(phres%gval(1,ith))/ceq%tpval(1)
!   write(*,70)'3H phres: ',ceq%tpval(1),phres%gval(1,1),phres%gval(2,1),&
!        phres%gval(3,1),phres%gval(4,1),kvot
! we should be careful with numeric overflow, for small T or large T
! no risk for overflow for exp(-kvot)
!   expmkvot=exp(-kvot)
!   ln1mexpkvot=log(one-expmkvot)
   if(kvot.gt.1.0D2) then
! T is very small, kvot very large, exp(kvot) may cause overflow, 
! exp(-kvot) is very small, ln(1-exp(-kvot)) is close to zero
! exp(kvot) may cause overflow, kvot/(exp(kvot)-1)=
! kvot*exp(-kvot)/(1-exp(-kvot)) = (1-kvot+kvot**2/2-...)/(1-kvot/2+...) = 1
      extreme=-1
      expmkvot=zero
      ln1mexpkvot=zero
      kvotexpkvotm1=zero
   elseif(kvot.lt.1.0D-2) then
! T is very big, kvot is very small, exp(-kvot) approch 1, 1-exp(-kvot)=kvot
! exp(-kvot) is close unity, ln(1-exp(-kvot))=ln(1-(1-kvot+kvot**2/2+...)) =
!            ln(kvot-kvot**2/2+...)=ln(kvot)
! exp(kvot) is close to unity: exp(kvot)-1 = kvot+kvot**2/2+ ...
      extreme=1
      expmkvot=exp(-kvot)
      ln1mexpkvot=log(one-expmkvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
   else
! normal range of T and kvot
      extreme=0
      expmkvot=exp(-kvot)
      ln1mexpkvot=log(one-expmkvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
   endif
! 
   gein=1.5D0*kvot+3.0D0*ln1mexpkvot
!   write(*,71)'3H Cp E1:',extreme,ceq%tpval(1),gein,ln1mexpkvot,expmkvot,&
!        kvotexpkvotm1
! first derivative wrt T taking care of overflow
   dgeindt=3.0D0*(ln1mexpkvot-kvotexpkvotm1)/ceq%tpval(1)
! This is d2G/dT**2/(RT) = -T**2/R*(Einstein Cp/RT) (or rather Cv/RT)
   if(extreme.eq.-1) then
! take care of overflow at low T
      d2geindt2=zero
   else
      d2geindt2=-3.0D0*kvotexpkvotm1**2/(expmkvot*ceq%tpval(1)**2)
   endif
! return the values in phres%gval(*,1)
   phres%gval(1,1)=phres%gval(1,1)+msize*gein
   phres%gval(2,1)=phres%gval(2,1)+msize*dgeindt
!   phres%gval(3,1)=phres%gval(3,1)
   phres%gval(4,1)=phres%gval(4,1)+msize*d2geindt2
!
! ?????????????????????????????
!
! NO DERIVATIVES WITH RESPECT TO FRACTIONS ??????????????????
!
! ?????????????????????????????
!
!   phres%gval(5,1)=phres%gval(5,1)
!   phres%gval(6,1)=phres%gval(6,1)
!   addrec%propval(1)=gein
!   addrec%propval(2)=dgeindt
!   addrec%propval(4)=d2geindt2
!   write(*,71)'3H Cp E3: ',extreme,ceq%tpval(1),gein,dgeindt,d2geindt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!  thet cannot depend on T
! Missing implementation of derivatives wrt comp.dep of thet.
   tv=ceq%tpval(1)
!-------------------------- two state part DIVIDED BY RT
! hump was an attempt to reduce the hump due to state change entropy
! it does not seem to work ... in fact it is the same as scaling G^HT-G^LT
! ****************** This is Johns orignal model *********************
!   hump=1.0
! Jump here if no Einstein solid
300 continue
   if(ig2.eq.0) then
      write(*,*)'Cannot find value for G2 two-state parameter'
      gx%bmperr=4367; goto 1000
   endif
! NOTE g2val and derivatives in phres%gval(..) are not divided by RT !!
   g2val=phres%gval(1,ig2); dg2dt=phres%gval(2,ig2)
   dg2=zero; d2g2dt2=zero
   if(g2val.eq.zero .and. dg2dt.eq.zero) then
!      write(*,*)'3H: G2 parameter zero, ignoring bump',g2val
      goto 900
   endif
!   write(*,19)'3H +am ',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
19 format(a,6(1pe11.3))
   rt=ceq%rtn
!   tv=ceq%tpval(1)
   rg=globaldata%rgas
!   expmg2=exp(-g2val/rt)
!   if(g2val/rt.gt.2.0D2) then
!      expmg2=exp(-g2val/(rt))
!      expg2=one/expmg2
!   elseif(g2val/rt.lt.1e-30) then
!      expmg2=exp(-g2val/(rt))
!      expg2=one/expmg2
!   else
      expmg2=exp(-g2val/(rt))
      expmg2p1=expmg2+one
      expg2=one/expmg2
!   endif
!   dg2=log(one+expmg2)
   dg2=log(expmg2p1)
!   write(*,19)'3H G2: ',g2val/rt,expmg2,dg2,dg2*rt
! NOTE values added to gval(*,1) must be divided by RT
! G = G - RT*ln(1+exp(-g2/RT))
! G
!   phres%gval(1,1)=phres%gval(1,1)-dg2
! (R*ln(1+g2val) + (g2/tv-dg2/dt)/(1+exp(-g2/RT)))/RT
! G.T
!   dgfdt=(rg*dg2+(g2val/tv-dg2dt)/(expg2+one))/(rt)
   dgfdt=(rg*dg2+(g2val/tv-dg2dt)/(expg2+one))/(rt)
!   phres%gval(2,1)=phres%gval(2,1)-dgfdt
! G.P   is zero
! ****************** This is Johns orignal model *********************
!-------------------------- tentative:
! d2g2/dt2/(1+exp(g2/RT)+
!   ((g2/tv)**2+(dg2/dt)**2-2*g2/tv*dg2/dt)*exp(g2/rt)/((1+exp(g2/RT)))**2/rt
! G.T.T 
! This what my derivation gives:
!   d2g2dt2=(phres%gval(4,ig2)/(one+expg2)+&
! works after Qing checked the signs
   d2g2dt2=(phres%gval(4,ig2)/(one+expg2)-&
        ((g2val/tv)**2+(dg2dt)**2-2.0D0*(g2val/tv)*dg2dt)*expg2/&
        (rt*(one+expg2)**2))/rt
!   phres%gval(4,1)=phres%gval(4,1)+d2g2dt2
! Maybe the error is here !!  YES now it works!
!   phres%gval(4,1)=phres%gval(4,1)-d2g2dt2
!   write(*,19)'3H dg2A: ',tv,  g2val/rt, one/(one+expg2),&
!        -rg*tv**2*phres%gval(4,1),&
!        -rg*tv**2*d2g2dt2,
! The Eistein contribution is OK
!   -rg*tv**2*d2geindt2
!        -rg*tv*phres%gval(4,ig2)/(one+expg2),&
!        -rg*tv*((g2val/tv)**2+(dg2dt)**2-2.0D0*(g2val/tv)*dg2dt)*expg2/&
!        (rt*(one+expg2)**2)
!   write(*,19)'3H dg2B: ',phres%gval(4,ig2)/(one+expg2),&
!        (g2val/tv)**2+(dg2dt)**2-2.0D0*(g2val/tv)*dg2dt,&
!        rt*(one+expg2)**2
! G.T.P is zero
! G.P.P is zero
800 continue
   phres%gval(1,1)=phres%gval(1,1)-msize*dg2
   phres%gval(2,1)=phres%gval(2,1)-msize*dgfdt
   phres%gval(4,1)=phres%gval(4,1)+msize*d2g2dt2
!
!
!   write(*,19)'3H 2st:',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
! save local values divided by RT?
! THIS ROUTINE CURRENTLY NOT USED <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
900 continue
   addrec%propval=zero
   addrec%propval(1)=msize*(gein-dg2)
   addrec%propval(2)=msize*(dgeindt-dgfdt)
   addrec%propval(4)=msize*(d2geindt2-d2g2dt2)
1000 continue
   write(*,*)'3H YOU ARE USING WRONG LIQUID TWOSTATE MODEL'
   return
 end subroutine calc_twostate_model_john

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_twostate_model1
!\begin{verbatim}
 subroutine calc_twostate_model1(moded,phres,addrec,lokph,mc,ceq)
! this routine is used when G2 and LNTH are composition dependent
! subroutine calc_twostate_modelny(moded,phres,addrec,lokph,mc,ceq)
! The routine _john works OK but I am testing a modification
! moded is 0, 1 or 2 if no, first or 2nd order derivatives should be calculated
! addrec is addition record
! phres is phase_varres record
! lokph is phase location
! mc is number of constitution fractions
! ceq is current equilibrium record
   implicit none
   integer moded,lokph,mc
   TYPE(GTP_PHASE_ADD), pointer :: addrec
   TYPE(GTP_PHASE_VARRES), pointer :: phres
   TYPE(GTP_EQUILIBRIUM_DATA), pointer :: ceq
!\end{verbatim} %+
! two state model for extrapolating liquid to low T
! DG = d(H-RT) + RT( dln(d)+(1-d)ln(1-d))
! where d is "liquid like" atoms.  H is enthalpy to form defects
! At equilibrium
!
! d = exp(-H/RT) / (1 + e(-H/RT) ) is the integrated Einstein Cp -H/R is THET
!
! G^liq - G^amorph = G^amorph - RT ln(1+exp(-DG_d/RT)
! DG_d is the enthalpy of forming 1 mole of defects in the glassy state
!
!------------------------------
! The value of Gd for the phase is calculated and added to G
   integer jj,noprop,ig2,ith,extreme,jth,kth
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2
   double precision gein,dgeindt,d2geindt2
   double precision xi,hump
   double precision, parameter :: humpfact=5.0D0
   logical addpermole
!   double precision g2ein,dg2eindt,d2g2eindt2,theta2,dcpl
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1
   double precision g2val,dg2,expg2,expmg2,rt,tv,rg,dg2dt,dgfdt,d2g2dt2
   double precision expmg2p1,fact,g2sum,msize,fact2
   double precision, allocatable :: mux(:)
   integer, save :: maxwarnings=0
! number of properties calculatied
   noprop=phres%listprop(1)-1
! locate the THET and G2 property record 
   ig2=0
   ith=0
   jth=0
! check if addition is per mole 
   if(btest(addrec%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize 2-state: ',lokph,addpermole,msize
   else
      addpermole=.FALSE.; msize=one
   endif
   findix: do jj=2,noprop
      if(phres%listprop(jj).eq.addrec%need_property(1)) then
! current values of G2 is stored in phres%gval(1,ig2)
         ig2=jj;
      elseif(phres%listprop(jj).eq.addrec%need_property(2)) then
! current value of LNTH are stored in phres%gval(1,ith)
         ith=jj
!      elseif(phres%listprop(jj).eq.22) then
! current value of DCP2 are stored in phres%gval(1,ith)
!         jth=jj
      endif
   enddo findix
   if(ith.eq.0) then
!      write(*,*)'3H Cannot find value for amorphous LNTH'
      if(maxwarnings.lt.20) then
         maxwarnings=maxwarnings+1
         write(*,*)'3H twostatemodel1 no values for amorphous LNTH:',maxwarnings
      endif
      gein=zero; dgeindt=zero; d2geindt2=zero
      goto 300
!      gx%bmperr=4367; goto 1000
   endif
   if(ig2.eq.0) then
      write(*,*)'3H twostate_model1 Cannot find G2 two-state parameter'
      gx%bmperr=4367; goto 1000
   endif
!----------------------------------
! for the moment the composition dependence is ignored
!   write(*,19)'3H 2no1: ',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
!------ this THET part copied from calc_einstein
! thet is in gval(1,ith), derivatives in dgval(*,ith,*) and d2gval(ith,*)
! G/RT = 1.5*THET/T + 3*R*LN(exp(THET/T) - 1) 
! NOTE ALL VALUES CALCULATED AS FOR G/RT
! kvot=theta/T
! NOTE the stored value is ln(theta! !!!
   kvot=exp(phres%gval(1,ith))/ceq%tpval(1)
!   write(*,70)'3H phres: ',ceq%tpval(1),phres%gval(1,1),phres%gval(2,1),&
!        phres%gval(3,1),phres%gval(4,1),kvot
! we should be careful with numeric overflow, for small T or large T
! no risk for overflow for exp(-kvot)
!   expmkvot=exp(-kvot)
!   ln1mexpkvot=log(one-expmkvot)
   if(kvot.gt.1.0D2) then
! T is very small, kvot very large, exp(kvot) may cause overflow, 
! exp(-kvot) is very small, ln(1-exp(-kvot)) is close to zero
! exp(kvot) may cause overflow, kvot/(exp(kvot)-1)=
! kvot*exp(-kvot)/(1-exp(-kvot)) = (1-kvot+kvot**2/2-...)/(1-kvot/2+...) = 1
      extreme=-1
      expmkvot=zero
      ln1mexpkvot=zero
      kvotexpkvotm1=zero
   elseif(kvot.lt.1.0D-2) then
! T is very big, kvot is very small, exp(-kvot) approch 1, 1-exp(-kvot)=kvot
! exp(-kvot) is close unity, ln(1-exp(-kvot))=ln(1-(1-kvot+kvot**2/2+...)) =
!            ln(kvot-kvot**2/2+...)=ln(kvot)
! exp(kvot) is close to unity: exp(kvot)-1 = kvot+kvot**2/2+ ...
      extreme=1
      expmkvot=exp(-kvot)
      ln1mexpkvot=log(one-expmkvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
   else
! normal range of T and kvot
      extreme=0
      expmkvot=exp(-kvot)
      ln1mexpkvot=log(one-expmkvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
   endif
! 
   gein=1.5D0*kvot+3.0D0*ln1mexpkvot
!   write(*,71)'3H Cp E1:',extreme,ceq%tpval(1),gein,ln1mexpkvot,expmkvot,&
!        kvotexpkvotm1
! first derivative wrt T taking care of overflow
   dgeindt=3.0D0*(ln1mexpkvot-kvotexpkvotm1)/ceq%tpval(1)
! This is d2G/dT**2/(RT) = -T**2/R*(Einstein Cp/RT) (or rather Cv/RT)
   if(extreme.eq.-1) then
! take care of overflow at low T
      d2geindt2=zero
   else
      d2geindt2=-3.0D0*kvotexpkvotm1**2/(expmkvot*ceq%tpval(1)**2)
   endif
! composition variable Variable LNTH
   fact=1.5D0*(one+expmkvot)/(one-expmkvot)*kvot
! the curve below better, correct shape ...
   do jj=1,mc
      phres%dgval(1,jj,1)=msize*(phres%dgval(1,jj,1)+fact*phres%dgval(1,jj,ith))
   enddo
!-------------------------- jump here if no LNTH variable
! return the values in phres%gval(*,1)
300 continue
   phres%gval(1,1)=phres%gval(1,1)+msize*gein
   phres%gval(2,1)=phres%gval(2,1)+msize*dgeindt
!   phres%gval(3,1)=phres%gval(3,1)
   phres%gval(4,1)=phres%gval(4,1)+msize*d2geindt2
!   phres%gval(5,1)=phres%gval(5,1)
!   phres%gval(6,1)=phres%gval(6,1)
!   write(*,71)'3H Cp E3: ',extreme,ceq%tpval(1),gein,dgeindt,d2geindt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!  thet cannot depend on T
! include the composition dependence of the eistein contribution? DONE ABOVE
! End of Einstein part
!----------------------------------------------------------------
!  write(*,*)'3H calc_twostate_model1 not including composition dependence thet'
!-------------------------- two state part DIVIDE BY RT
! NOTE the values in phres%gval(1,ig2), phres%dgval(1,jj,ig2)
!        are not divided by T.
!
   rt=ceq%rtn
   tv=ceq%tpval(1)
   rg=globaldata%rgas
   g2val=phres%gval(1,ig2); dg2dt=phres%gval(2,ig2)
   dg2=zero; d2g2dt2=zero; expmg2=zero
!   write(*,*)'3H gval1: ',g2val
   if(g2val.eq.zero .and. dg2dt.eq.zero) then
      write(*,*)'3H: G2 parameter zero, ignoring 2-state model'
      goto 900
   endif
   d2g2dt2=phres%gval(4,ig2)
   goto 600
!------------------------------------------
600 continue
! if g2val is positive we are in the amorphous region
! if g2val is negative we are in the liquid region
! The if statements here ensure expmg2 is between 1e-60 and 1e+60
!   write(*,'(a,6(1pe12.4))')'3H g2val: ',g2val,dg2dt,-g2val/rt
   if(-g2val/rt.gt.2.0D2) then
! LIQUID REGION exp(200) >> 1, thus d2g=ln(1+exp(g2val))=g2val
! and the derivatives are those above. DIVIDED BY RT?
      dg2=g2val/rt
      dgfdt=dg2dt/rt
      d2g2dt2=d2g2dt2/rt
      goto 700
   elseif(-g2val/rt.lt.-2.0D2) then
! Low T AMORPHOUS REGION: exp(-200)=0; ln(1)=0 and everything is zero
      dg2=zero
      dg2dt=zero
      d2g2dt2=zero
      goto 900
   else
! intermediate T range, we have to calculate, exp( -200 to +200) is OK
      expmg2=exp(-g2val/rt)
      expg2=one/expmg2
      expmg2p1=expmg2+one
      dg2=log(expmg2p1)
!      write(*,'(a,4(1pe12.4))')'3H intermed: ',phres%gval(1,ig2),&
!           g2val/rt,expmg2p1,dg2
   endif
!   write(*,19)'3H gval8: ',g2val/rt,expmg2,dg2
!   write(*,19)'3H dg2: ',tv,g2val,expmg2,dg2
!   write(*,19)'3H G2: ',tv,xi,g2val/rt,expmg2,dg2,dg2*rt
! NOTE values added to phres%gval(*,1) must be divided by RT
! G = G - RT*ln(1+exp(-g2/RT))
! G
!   phres%gval(1,1)=phres%gval(1,1)-dg2
! (R*ln(1+g2val) + (g2/tv-dg2/dt)/(1+exp(-g2/RT)))/RT
! G.T
   dgfdt=(rg*dg2+(g2val/tv-dg2dt)/(expg2+one))/rt
!   dgfdt=dg2+(g2val/tv-dg2dt)/(expg2+one)
! G.P   is zero
!-------------------------- tentative:
! d2g2/dt2/(1+exp(g2/RT)+
!   ((g2/tv)**2+(dg2/dt)**2-2*g2/tv*dg2/dt)*exp(g2/rt)/((1+exp(g2/RT)))**2/rt
! G.T.T 
! Fixed sign problem
!   d2g2dt2=(phres%gval(4,ig2)/(one+expg2)-&
   d2g2dt2=(d2g2dt2/(one+expg2)-&
        ((g2val/tv)**2+(dg2dt)**2-2.0D0*(g2val/tv)*dg2dt)*expg2/&
        (rt*(one+expg2)**2))/rt
700 continue
!   write(*,705)'3H 2SL: ',g2val/rt, dg2, dgfdt, dgfdt, d2g2dt2, tv,&
!        rt, expg2, dg2dt, msize, d2g2dt2*rt
705 format(a,6(1pe12.4)/8x,6(1pe12.4))
! THIS IS THE SUBROUTINE USED FOR 2STATE LIQUID with composition dependent GD
! This should be OK/ 2020.02.27
   phres%gval(1,1)=phres%gval(1,1)-msize*dg2
   phres%gval(2,1)=phres%gval(2,1)-msize*dgfdt
   phres%gval(4,1)=phres%gval(4,1)+msize*d2g2dt2
! values of T, \xi, g, s and cp   
!
! ADDING DERIVATIVES WITH RESPECT TO FRACTIONS !!!!!!!!!!!
!
   fact=expmg2/(expmg2+one)/rt
   fact2=(fact/rt)**2/expmg2*phres%gval(2,ig2)
!   write(*,*)'3H calculating twostatemodel, wrong dg/dy?',fact2
   do jj=1,mc
      phres%dgval(1,jj,1)=phres%dgval(1,jj,1)+fact*phres%dgval(1,jj,ig2)
! d2(G2)/dydT
!      phres%dgval(1,jj,1)=phres%dgval(2,jj,1)-fact2*phres%dgval(1,jj,ig2)
! ignore other 2nd derivatives
   enddo
   goto 900
!============================================================
!----------skipping old code below
! Searching for bug when entering G2 as a comp.dependent parameter rather
! than as a part of the pure element data. Liquid with: T=1950; x(v)=.5
! 1. OC exactly same as TC when G2 is part of pure elements:
!    G=-127282 J; a(ti)=2.9601E-4; a(v)=5.1269E-4 (SER refstate)
! 2. When modeling G2 as a separate parameter we get:
!    TC: G=-127223 J; ac(ti)=2.9714E-4; ac(v)=5.1446E-4;  <<<<<<<<<<<<
!    OC: G=-127223 J; ac(ti)=3.3896E-4; ac(v)=4.5099E-4 (divided by RT)
! 3. The chemical potential wrong and give strange (wrong) phase diagram
!    OC: divide by T: G=-127223 J; ac(ti)=2.9714E-4; ac(v)=5.1446E-4 WoW
!  but the phase diagram still wrong, not same composition at BCC/LIQ minimum!!
!    BCC is identical in TC and OC, no problem with Einstein
! 4. Calculating liquid at T=1920; x(v)=.1 gives different results:
!    TC: G=-122239 J; ac(ti)=5.4656E-4; ac(v)=1.2773E-4
!    OC: G=-122239 J; ac(ti)=5.2677E-4; ac(v)=1.7801E-4
!    T=1900, x(v)=.1 (for BCC!!)
!    TC: G=-120324 J; ac(ti)=5.6252E-4; ac(v)=1.4795E-4 BCC! DGM(liq)=-4.49E-3
!    OC: G=-120324 J; ac(ti)=5.6252E-4; ac(v)=1.4795E-4 BCC! DGM(liq)=-1.02E-2
!    T=1910, x(v)=.1 (TC gives 2-phase equil, BCC+LIQ, OC just BCC)
! 5. After some more changes (introducing msize mm) which should NOT change
!    the result, it is bad again ... SUCK
!    OC: divide by T: G=-127223 J; ac(ti)=2.5623E-4; ac(v)=5.9662E-4
! does some variable have random a value??  BCC stable at this calc
!    OC: divide by RT: G=-127223 J; ac(ti)=3.3896E-4; ac(v)=4.5099E-4 (as above)
!    mulip with 0.5/T: G=-127223 J; ac(ti)=3.0040E-4; ac(v)=5.0889E-4
!    muli with 0.52/T: G=-127223 J; ac(ti)=2.9849E-4; ac(v)=5.1214E-4
!    muli with 0.53/T: G=-127223 J; ac(ti)=2.9754E-4; ac(v)=5.1377E-4
! this also gave reasonable phase diagram!!
!    mul with 0.533/T: G=-127223 J; ac(ti)=2.9726E-4; ac(v)=5.1426E-4
!    mu with 0.5335/T: G=-127223 J; ac(ti)=2.9721E-4; ac(v)=5.1434E-4
!    mu with 0.5338/T: G=-127223 J; ac(ti)=2.9818E-4; ac(v)=5.1439E-4
!    mul with 0.534/T: G=-127223 J; ac(ti)=2.9717E-4; ac(v)=5.1442E-4
!    mu with 0.5342/T: G=-127223 J; ac(ti)=2.9715E-4; ac(v)=5.1446E-4
! phase diagram correct with multiplied factor with 0.534... Wow, why?
!    correct:          G=-127223 J; ac(ti)=2.9714E-4; ac(v)=5.1446E-4
!    muli with 0.54/T: G=-127223 J; ac(ti)=2.9660E-4; ac(v)=5.1541E-4
!---------------------------------------------
! COMPLETELY NONSCIENTIFIC ... ENGINEERING 
!---------------------------------------------
! the factor /rt because phres%dgval(1,jj,ig2) is not divided by RT
!   fact=expmg2/(expmg2+one)/rt
! no RT?? No, activites zero
!   fact=expmg2/(expmg2+one)
!>>>>>>>>>>>>>> this factor 0.5336 gives correct chemical potentials
   fact=0.5342D0*msize*expmg2/(expmg2+one)/tv
!>>>>>>>>>>>>>> but I do not understad why
!   fact=msize*expmg2/(expmg2+one)/rt
! sign?? - no good; 
!   fact=10*expmg2/(expmg2+one)/rt
!!   fact=10*expmg2/(expmg2+one)/rt
!   fact=1.2E0*expmg2/(expmg2+one)/tv
!   write(*,'(a,6(1pe12.4))')'3H missing dg2/dy: ',expmg2,tv,msize,fact
! temporary debug ...; calculate contribution to chemical potential
!   allocate(mux(mc))
!   mux=zero
   g2sum=dg2
   do jj=1,mc
! check than x1*mu1+x2*mu2=g
!      g2sum=g2sum-phres%yfr(jj)*fact*phres%dgval(1,jj,ig2)
! in dgval(i,j,k) index i=1 means d/y; 2 means d2/dydT, 3 means d2/dydP
!                 index j is constituent; index k is property, k=1 is G
!      write(*,710)ig2,jj,phres%dgval(1,jj,1),phres%dgval(1,jj,ig2)/rt,&
!           fact*phres%dgval(1,jj,ig2),&
!           phres%dgval(1,jj,1)+fact*phres%dgval(1,jj,ig2),phres%yfr(jj)
!           phres%dgval(2,jj,1),phres%dgval(2,jj,ig2)/rt,&
!           phres%dgval(2,jj,1)+fact*phres%dgval(2,jj,ig2)
710   format('3H G2: ',2i2,6(1pe11.3))
      phres%dgval(1,jj,1)=phres%dgval(1,jj,1)+fact*phres%dgval(1,jj,ig2)
      phres%dgval(2,jj,1)=phres%dgval(2,jj,1)+fact*phres%dgval(2,jj,ig2)
   enddo
!   do jj=1,mc
!      mux(jj)=g2sum+fact*phres%dgval(1,jj,ig2)
!   enddo! These should be the same!
!   g2sum=zero
!   do jj=1,mc
!      g2sum=g2sum+phres%yfr(jj)*mux(jj)
!   enddo
!   write(*,*)'3H same?: ',dg2,g2sum
800 continue
!   write(*,19)'3H G9: ',tv,hump,dg2*rt,-dgfdt*rt,-d2g2dt2*rt*tv
!   phres%gval(4,1)=phres%gval(4,1)-d2g2dt2
!   write(*,19)'3H dg2A: ',tv,  g2val/rt, one/(one+expg2),&
!        -rg*tv**2*phres%gval(4,1),&
!        -rg*tv**2*d2g2dt2,
! The Eistein contribution is OK
!   -rg*tv**2*d2geindt2
!        -rg*tv*phres%gval(4,ig2)/(one+expg2),&
!        -rg*tv*((g2val/tv)**2+(dg2dt)**2-2.0D0*(g2val/tv)*dg2dt)*expg2/&
!        (rt*(one+expg2)**2)
!   write(*,19)'3H dg2B: ',phres%gval(4,ig2)/(one+expg2),&
!        (g2val/tv)**2+(dg2dt)**2-2.0D0*(g2val/tv)*dg2dt,&
!        rt*(one+expg2)**2
! G.T.P is zero
! G.P.P is zero
!   write(*,19)'3H 2st:',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
! jump here skipping old code above
!============================================================
! save Einstein and G2 values in addrec, multiply with RT?
900 continue
   addrec%propval=zero
   addrec%propval(1)=msize*(gein-dg2)
   addrec%propval(2)=msize*(dgeindt-dgfdt)
   addrec%propval(4)=msize*(d2geindt2-d2g2dt2)
!
1000 continue
   return
! this routine is used when G2 and LNTH are composition dependent
 end subroutine calc_twostate_model1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_twostate_model2
!\begin{verbatim}
 subroutine calc_twostate_model2(moded,phres,addrec,lokph,mc,ceq)
! subroutine calc_twostate_model_nomix(moded,phres,addrec,lokph,mc,ceq)
! moded is 0, 1 or 2 if no, first or 2nd order derivatives should be calculated
! addrec is addition record
!
! IN THIS VERSION G2 is treated as a composition independent parameter 
!  thus this just handles the Einsten Cp
! NOTE Einstein LNTH should be composition dependent
!
! phres is phase_varres record
! lokph is phase location
! mc is number of constitution fractions
! ceq is current equilibrium record
   implicit none
   integer moded,lokph,mc
   TYPE(GTP_PHASE_ADD), pointer :: addrec
   TYPE(GTP_PHASE_VARRES), pointer :: phres
   TYPE(GTP_EQUILIBRIUM_DATA), pointer :: ceq
!\end{verbatim} %+
! two state model for extrapolating liquid to low T
! DG = d(H-RT) + RT( dln(d)+(1-d)ln(1-d))
! where d is "liquid like" atoms.  H is enthalpy to form defects
! At equilibrium
!
! d = exp(-H/RT) / (1 + e(-H/RT) ) is the integrated Einstein Cp -H/R is THET
!
! G^liq - G^amorph = G^amorph - RT ln(1+exp(-DG_d/RT)
! DG_d is the enthalpy of forming 1 mole of defects in the glassy state
!
!------------------------------
! The value of Gd for the phase is calculated and added to G
   integer jj,noprop,ig2,ith,extreme,jth,kth
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2
   double precision gein,dgeindt,d2geindt2
   double precision xi,hump
   double precision, parameter :: humpfact=5.0D0
   logical addpermole
!   double precision g2ein,dg2eindt,d2g2eindt2,theta2,dcpl
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1
   double precision g2val,dg2,expg2,expmg2,rt,tv,rg,dg2dt,dgfdt,d2g2dt2
   double precision expmg2p1,fact,g2sum,msize
   double precision, allocatable :: mux(:)
   integer, save :: maxwarnings=0
! number of properties calculatied
   noprop=phres%listprop(1)-1
! locate the THET and G2 property record 
   ig2=0
   ith=0
   jth=0
! check if addition is per mole 
   if(btest(addrec%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize 2-state: ',lokph,addpermole,msize
   else
      addpermole=.FALSE.; msize=one
   endif
   findix: do jj=2,noprop
! start from 2 as phres%listprop(1) is always G
      if(phres%listprop(jj).eq.addrec%need_property(1)) then
! current value of THET are stored in phres%gval(1,ith)
         ith=jj
      endif
   enddo findix
   if(ith.eq.0) then
!      write(*,*)'3H Cannot find value for amorphous LNTH'
      if(maxwarnings.lt.20) then
         maxwarnings=maxwarnings+1
         write(*,*)'3H twostatemodel2 no values for amorphous LNTH:',maxwarnings
      endif
      gein=zero; dgeindt=zero; d2geindt2=zero
      goto 1000
!      goto 300
!      gx%bmperr=4367; goto 1000
   endif
!   write(*,*)'3H Using composition independent G2 values, THET=',&
!        phres%gval(1,ith)
!   if(ig2.eq.0) then
!      write(*,*)'3H Cannot find value for G2 two-state parameter'
!      gx%bmperr=4367; goto 1000
!   endif
!----------------------------------
! for the moment the composition dependence is ignored
!   write(*,19)'3H 2no1: ',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
!------ this THET part copied from calc_einstein
! thet is in gval(1,ith), derivatives in dgval(*,ith,*) and d2gval(ith,*)
! G/RT = 1.5*THET/T + 3*R*LN(exp(THET/T) - 1) 
! NOTE ALL VALUES CALCULATED AS FOR G/RT
! kvot=theta/T
! NOTE the stored value is ln(theta)! !!!
   kvot=exp(phres%gval(1,ith))/ceq%tpval(1)
!   write(*,70)'3H phres: ',ceq%tpval(1),phres%gval(1,1),phres%gval(2,1),&
!        phres%gval(3,1),phres%gval(4,1),kvot
! we should be careful with numeric overflow, for small T or large T
! no risk for overflow for exp(-kvot)
!   expmkvot=exp(-kvot)
!   ln1mexpkvot=log(one-expmkvot)
   if(kvot.gt.1.0D2) then
! T is very small, kvot very large, exp(kvot) may cause overflow, 
! exp(-kvot) is very small, ln(1-exp(-kvot)) is close to zero
! exp(kvot) may cause overflow, kvot/(exp(kvot)-1)=
! kvot*exp(-kvot)/(1-exp(-kvot)) = (1-kvot+kvot**2/2-...)/(1-kvot/2+...) = 1
      extreme=-1
      expmkvot=zero
      ln1mexpkvot=zero
      kvotexpkvotm1=zero
   elseif(kvot.lt.1.0D-2) then
! T is very big, kvot is very small, exp(-kvot) approch 1, 1-exp(-kvot)=kvot
! exp(-kvot) is close unity, ln(1-exp(-kvot))=ln(1-(1-kvot+kvot**2/2+...)) =
!            ln(kvot-kvot**2/2+...)=ln(kvot)
! exp(kvot) is close to unity: exp(kvot)-1 = kvot+kvot**2/2+ ...
      extreme=1
      expmkvot=exp(-kvot)
      ln1mexpkvot=log(one-expmkvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
   else
! normal range of T and kvot
      extreme=0
      expmkvot=exp(-kvot)
      ln1mexpkvot=log(one-expmkvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
   endif
! 
   gein=1.5D0*kvot+3.0D0*ln1mexpkvot
!   write(*,71)'3H Cp E1:',extreme,ceq%tpval(1),gein,ln1mexpkvot,expmkvot,&
!        kvotexpkvotm1
! first derivative wrt T taking care of overflow
   dgeindt=3.0D0*(ln1mexpkvot-kvotexpkvotm1)/ceq%tpval(1)
! This is d2G/dT**2/(RT) = -T**2/R*(Einstein Cp/RT) (or rather Cv/RT)
   if(extreme.eq.-1) then
! take care of overflow at low T
      d2geindt2=zero
   else
      d2geindt2=-3.0D0*kvotexpkvotm1**2/(expmkvot*ceq%tpval(1)**2)
   endif
! derivatives with respect to composition dependence of THET
   fact=1.5D0*(one+expmkvot)/(one-expmkvot)*kvot
   do jj=1,mc
      phres%dgval(1,jj,1)=msize*(phres%dgval(1,jj,1)+fact*phres%dgval(1,jj,ith))
   enddo
! Ignore second derivatives as this seems a small efect, 
!-------------------------- jump here if no THET variable
! return the values in phres%gval(*,1)
300 continue
   phres%gval(1,1)=phres%gval(1,1)+msize*gein
   phres%gval(2,1)=phres%gval(2,1)+msize*dgeindt
!   phres%gval(3,1)=phres%gval(3,1)
   phres%gval(4,1)=phres%gval(4,1)+msize*d2geindt2
!   phres%gval(5,1)=phres%gval(5,1)
!   phres%gval(6,1)=phres%gval(6,1)
!   write(*,71)'3H Cp E3: ',extreme,ceq%tpval(1),gein,dgeindt,d2geindt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!----------------------------------------------------------------
! skip the 2-state model as G2 included in the ^oG for the endmember
!   goto 1000
1000 continue
   return
 end subroutine calc_twostate_model2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_twostate_endmember
!\begin{verbatim}
 subroutine calc_twostate_model_endmember(proprec,g2values,ceq)
! This calculated G2 (GD in some papers) for a pure endmember
! No composition dependence
! Value calculated here added to ^oG for the endmember
! moded is 0, 1 or 2 if no, first or 2nd order derivatives should be calculated
! phres is phase_varres record
! lokph is phase location
! ceq is current equilibrium record
   implicit none
   TYPE(gtp_property), pointer :: proprec
   TYPE(GTP_EQUILIBRIUM_DATA), pointer :: ceq
   double precision g2values(6)
!\end{verbatim}
   TYPE(gtp_property), pointer :: propg2
   integer lokfun,typty
   double precision tv,rt,rg,dg2,dgfdt,expg2,expmg2,expmg2p1
   double precision g2val,dg2dt,d2g2dt2,vals(6)
   g2values=zero
! do not destroy the value of proprec!!
   propg2=>proprec
! At present %proptype 16 is G2 but can be changed anytime!!
! However, it will always be called G2 and need_property('G2 ',typty)
! will return is current index
   call need_propertyid('G2  ',typty)
!   write(*,*)'3H found G2 typty: ',typty
   liq2state: do while(associated(propg2))
! How to find addrec%need_property ....??
!      if(phres%listprop(jl).eq.addrec%need_property(1)) then
!         ilpx=jl
!      write(*,*)'3H property: ',propg2%proptype
      if(propg2%proptype.eq.typty) goto 77
      propg2=>propg2%nextpr
   enddo liq2state
   write(*,*)'Missing liquid twostate parameter G2'
   gx%bmperr=4399; goto 1000
! found p   
77 continue
! calculate G2 value at current T for the endmember
   lokfun=propg2%degreelink(0)
   call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,'(a,3(1pe12.4))')'3H G2 endmember: ',vals(1),vals(2),vals(4)
   g2val=vals(1)
   dg2dt=vals(2)
   d2g2dt2=vals(4)
!   
   rt=ceq%rtn
   tv=ceq%tpval(1)
   rg=globaldata%rgas
! if g2val is positive we are in the amorphous region
! if g2val is negative we are in the liquid region
! The if statements here ensure expmg2 is between 1e-60 and 1e+60
!   write(*,'(a,6(1pe12.4))')'3H g2val: ',g2val,dg2dt,-g2val/rt
   if(-g2val/rt.gt.2.0D2) then
! LIQUID REGION exp(200) >> 1, thus d2g=ln(1+exp(g2val))=g2val
! and the derivatives are those above. DIVIDED BY RT?
      dg2=g2val/rt
      dgfdt=dg2dt/rt
      d2g2dt2=d2g2dt2/rt
      goto 700
   elseif(-g2val/rt.lt.-2.0D2) then
! AMORPHOUS REGION: exp(-200)=0; ln(1)=0 and everything is zero
      dg2=zero
      dg2dt=zero
      d2g2dt2=zero
      goto 700
   else
! intermediate T range, we have to calculate, exp( -200 to +200) is OK
      expmg2=exp(-g2val/rt)
      expg2=one/expmg2
      expmg2p1=expmg2+one
      dg2=log(expmg2p1)
!      write(*,'(a,4(1pe12.4))')'3H intermed: ',g2val/rt,expmg2p1,dg2
   endif
!   write(*,19)'3H gval8: ',g2val/rt,expmg2,dg2
!   write(*,19)'3H dg2: ',tv,g2val,expmg2,dg2
!   write(*,19)'3H G2: ',tv,xi,g2val/rt,expmg2,dg2,dg2*rt
! NOTE values added to phres%gval(*,1) must be divided by RT
! G = G - RT*ln(1+exp(-g2/RT))
! G
!   phres%gval(1,1)=phres%gval(1,1)-dg2
! (R*ln(1+g2val) + (g2/tv-dg2/dt)/(1+exp(-g2/RT)))/RT
! G.T
   dgfdt=(rg*dg2+(g2val/tv-dg2dt)/(expg2+one))/rt
!   dgfdt=dg2+(g2val/tv-dg2dt)/(expg2+one)
! G.P   is zero
!-------------------------- tentative:
! d2g2/dt2/(1+exp(g2/RT)+
!   ((g2/tv)**2+(dg2/dt)**2-2*g2/tv*dg2/dt)*exp(g2/rt)/((1+exp(g2/RT)))**2/rt
! G.T.T 
! Fixed sign problem
!   d2g2dt2=(phres%gval(4,ig2)/(one+expg2)-&
   d2g2dt2=(d2g2dt2/(one+expg2)-&
        ((g2val/tv)**2+(dg2dt)**2-2.0D0*(g2val/tv)*dg2dt)*expg2/&
        (rt*(one+expg2)**2))/rt
700 continue
! return these values to be added to ^oG for the endmember
!   g2values(1)=-g2val/rt
   g2values(1)=-dg2
!   g2values(2)=-dg2dt
   g2values(2)=-dgfdt
   g2values(4)=d2g2dt2
!   write(*,'(a,3(1pe12.4))')'3H g2values: ',g2values(1),g2values(2),g2values(4)
! No P derivatives (yet)
!   write(*,705)'3H 2SL: ',g2val/rt, dg2, dgfdt, dgfdt, d2g2dt2, tv,&
!        rt, expg2, dg2dt, msize, d2g2dt2*rt
705 format(a,6(1pe12.4)/8x,6(1pe12.4))
! each endmember has its own value of G2
!   phres%gval(1,1)=phres%gval(1,1)-msize*dg2
!   phres%gval(2,1)=phres%gval(2,1)-msize*dgfdt
!   phres%gval(4,1)=phres%gval(4,1)+msize*d2g2dt2
! values of T, \xi, g, s and cp   
!
1000 continue
   return
 end subroutine calc_twostate_model_endmember

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_twostate_model_old
!\begin{verbatim}
 subroutine calc_twostate_model_old(moded,phres,addrec,lokph,mc,ceq)
! Failed attempt to decrease the hump when the g2 parameter changes sign
! moded is 0, 1 or 2 if no, first or 2nd order derivatives should be calculated
! addrec is addition record
! phres is phase_varres record
! lokph is phase location
! mc is number of constitution fractions
! ceq is current equilibrium record
   implicit none
   integer moded,lokph,mc
   TYPE(GTP_PHASE_ADD), pointer :: addrec
   TYPE(GTP_PHASE_VARRES), pointer :: phres
   TYPE(GTP_EQUILIBRIUM_DATA), pointer :: ceq
!\end{verbatim}
! two state model for extrapolating liquid to low T
! DG = d(H-RT) + RT( dln(d)+(1-d)ln(1-d))
! where d is "liquid like" atoms.  H is enthalpy to form defects
! At equilibrium
!
! d = exp(-H/RT) / (1 + e(-H/RT) ) is the integrated Einstein Cp -H/R is THET
!
! G^liq - G^amorph = G^amorph - RT ln(1+exp(-DG_d/RT)
! DG_d is the enthalpy of forming 1 mole of defects in the glassy state
!
!------------------------------
! The value of Gd for the phase is calculated and added to G
   integer jj,noprop,ig2,ith,extreme,m4
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2
   double precision gein,dgeindt,d2geindt2
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1
   double precision g2val,dg2,expg2,expmg2,rt,tv,rg,dg2dt,dgfdt,d2g2dt2
   double precision hump,fq,dfq,d2fq,addq,daddq,d2addq,dd
! number of properties calculatied
   noprop=phres%listprop(1)-1
! locate the LNTH and G2 property record 
   ig2=0
   ith=0
   findix: do jj=2,noprop
      if(phres%listprop(jj).eq.addrec%need_property(1)) then
! current values of G2 is stored in phres%gval(1,ig2)
         ig2=jj;
      elseif(phres%listprop(jj).eq.addrec%need_property(2)) then
! current value of THET are stored in phres%gval(1,ith)
         ith=jj
      endif
   enddo findix
   if(ith.eq.0) then
      write(*,*)'Cannot find value for amorphous LNTH'
      gx%bmperr=4399; goto 1000
   endif
   if(ig2.eq.0) then
      write(*,*)'Cannot find value for G2 two-state parameter'
      gx%bmperr=4399; goto 1000
   endif
!----------------------------------
! for the moment the composition dependence is ignored
!   write(*,19)'3H 2no1: ',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
!------ this THET part copied from calc_einstein
! thet is in gval(1,ith), derivatives in dgval(*,ith,*) and d2gval(ith,*)
! G/RT = 1.5*THET/T + 3*R*LN(exp(THET/T) - 1) 
! NOTE ALL VALUES CALCULATED AS FOR G/RT
! kvot=theta/T
! NOTE the stored value is ln(theta! !!!
   kvot=exp(phres%gval(1,ith))/ceq%tpval(1)
!   write(*,70)'3H phres: ',ceq%tpval(1),phres%gval(1,1),phres%gval(2,1),&
!        phres%gval(3,1),phres%gval(4,1),kvot
! we should be careful with numeric overflow, for small T or large T
! no risk for overflow for exp(-kvot)
!   expmkvot=exp(-kvot)
!   ln1mexpkvot=log(one-expmkvot)
   if(kvot.gt.1.0D2) then
! T is very small, kvot very large, exp(kvot) may cause overflow, 
! exp(-kvot) is very small, ln(1-exp(-kvot)) is close to zero
! exp(kvot) may cause overflow, kvot/(exp(kvot)-1)=
! kvot*exp(-kvot)/(1-exp(-kvot)) = (1-kvot+kvot**2/2-...)/(1-kvot/2+...) = 1
      extreme=-1
      expmkvot=zero
      ln1mexpkvot=zero
      kvotexpkvotm1=zero
   elseif(kvot.lt.1.0D-2) then
! T is very big, kvot is very small, exp(-kvot) approch 1, 1-exp(-kvot)=kvot
! exp(-kvot) is close unity, ln(1-exp(-kvot))=ln(1-(1-kvot+kvot**2/2+...)) =
!            ln(kvot-kvot**2/2+...)=ln(kvot)
! exp(kvot) is close to unity: exp(kvot)-1 = kvot+kvot**2/2+ ...
      extreme=1
      expmkvot=exp(-kvot)
      ln1mexpkvot=log(one-expmkvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
   else
! normal range of T and kvot
      extreme=0
      expmkvot=exp(-kvot)
      ln1mexpkvot=log(one-expmkvot)
      kvotexpkvotm1=kvot/(exp(kvot)-one)
   endif
! 
   gein=1.5D0*kvot+3.0D0*ln1mexpkvot
!   write(*,71)'3H Cp E1:',extreme,ceq%tpval(1),gein,ln1mexpkvot,expmkvot,&
!        kvotexpkvotm1
! first derivative wrt T taking care of overflow
   dgeindt=3.0D0*(ln1mexpkvot-kvotexpkvotm1)/ceq%tpval(1)
! This is d2G/dT**2/(RT) = -T**2/R*(Einstein Cp/RT) (or rather Cv/RT)
   if(extreme.eq.-1) then
! take care of overflow at low T
      d2geindt2=zero
   else
      d2geindt2=-3.0D0*kvotexpkvotm1**2/(expmkvot*ceq%tpval(1)**2)
   endif
! return the values in phres%gval(*,1)
   phres%gval(1,1)=phres%gval(1,1)+gein
   phres%gval(2,1)=phres%gval(2,1)+dgeindt
!   phres%gval(3,1)=phres%gval(3,1)
   phres%gval(4,1)=phres%gval(4,1)+d2geindt2
!   phres%gval(5,1)=phres%gval(5,1)
!   phres%gval(6,1)=phres%gval(6,1)
   addrec%propval(1)=gein
   addrec%propval(2)=dgeindt
   addrec%propval(4)=d2geindt2
!   write(*,71)'3H Cp E3: ',extreme,ceq%tpval(1),gein,dgeindt,d2geindt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!  thet cannot depend on T
! Missing implem of derivatives wrt comp.dep of thet.
!-------------------------- two state part DIVIDE BY RT
! NOTE g2val and derivatives not divided by RT !!
   g2val=phres%gval(1,ig2)
   dg2dt=phres%gval(2,ig2)
   dg2=zero; d2g2dt2=zero
   if(g2val.eq.zero .and. dg2dt.eq.zero) then
      write(*,*)'3H: G2 parameter zero, ignoring twostate model',g2val
      goto 900
   endif
!   write(*,19)'3H +am ',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
19 format(a,6(1pe11.3))
   rt=ceq%rtn
   tv=ceq%tpval(1)
   rg=globaldata%rgas
!   expmg2=exp(-g2val/rt)
!   if(g2val/rt.gt.2.0D2) then
!      expmg2=exp(-g2val/(rt))
!      expg2=one/expmg2
!   elseif(g2val/rt.lt.1e-30) then
!      expmg2=exp(-g2val/(rt))
!      expg2=one/expmg2
!   else
      expmg2=exp(-g2val/(rt))
      expg2=one/expmg2
!   endif
!   dg2=log(one+expmg2)
   dg2=log(one+expmg2)
!   write(*,19)'3H G2: ',g2val/rt,expmg2,dg2,dg2*rt
! NOTE values added to gval(*,1) must be divided by RT
! G = G - RT*ln(1+exp(-g2/RT))
! G
   phres%gval(1,1)=phres%gval(1,1)-dg2
! (R*ln(1+g2val) + (g2/tv-dg2/dt)/(1+exp(-g2/RT)))/RT
! G.T
!   dgfdt=(rg*dg2+(g2val/tv-dg2dt)/(expg2+one))/(rt)
   dgfdt=(rg*dg2+(g2val/tv-dg2dt)/(expg2+one))/(rt)
   phres%gval(2,1)=phres%gval(2,1)-dgfdt
! G.P   is zero
!-------------------------- tentative:
! d2g2/dt2/(1+exp(g2/RT)+
!   ((g2/tv)**2+(dg2/dt)**2-2*g2/tv*dg2/dt)*exp(g2/rt)/((1+exp(g2/RT)))**2/rt
! G.T.T 
! This what my derivation gives:
!   d2g2dt2=(phres%gval(4,ig2)/(one+expg2)+&
! Qing proposal, works after fixing the sign also below
   d2g2dt2=(phres%gval(4,ig2)/(one+expg2)-&
! This is which is the same as TC
!   d2g2dt2=(-phres%gval(4,ig2)/(one+expg2)+&
        ((g2val/tv)**2+(dg2dt)**2-2.0D0*(g2val/tv)*dg2dt)*expg2/&
        (rt*(one+expg2)**2))/rt
! Maybe the error is here !!  YES now it works!
   phres%gval(4,1)=phres%gval(4,1)+d2g2dt2
!   phres%gval(4,1)=phres%gval(4,1)-d2g2dt2
! G.T.P is zero
! G.P.P is zero
! This is my addition to the two-state model to control the size of the hump
!   goto 1000
   hump=1.0D0
   m4=2
   fq=g2val/rt
   dfq=dg2dt/rt-fq/tv
   d2fq=phres%gval(4,ith)/rt-2.0D0/(rt*tv)*(dg2dt+g2val/tv)
   dd=one+(2.0D-1*fq)**m4
   addq=hump/dd
   daddq=-m4*hump*fq**(m4-1)*dfq/dd**2
   d2addq=-m4*hump*fq**(m4-2)*((m4-1)*dfq**2+fq*d2fq)/dd**2+&
        2.0d0*m4**2*hump*fq**(2*m4-2)*dfq**2/dd**3
! ignoring T dependence
   d2addq=5.0D-5/dd
   phres%gval(1,1)=phres%gval(1,1)+addq
   phres%gval(2,1)=phres%gval(2,1)+daddq
   phres%gval(4,1)=phres%gval(4,1)+d2addq
   write(*,800)'3H added hump',tv,fq,dd,-rg*tv**2*d2addq
800 format(a,6(1pe11.3))
! save local values divided by RT?
900 continue
   addrec%propval=zero
   addrec%propval(1)=gein-dg2
   addrec%propval(2)=dgeindt-dgfdt
   addrec%propval(4)=d2geindt2-d2g2dt2
1000 continue
   return
 end subroutine calc_twostate_model_old

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_debyecp
!\begin{verbatim}
 subroutine create_debyecp(addrec)
! enters a record for the debye model
   implicit none
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim} %+
   integer typty
! reserve an addition record
   allocate(addrec)
! Set the type of addition and look for needed parameter properties
   addrec%type=debyecp
   addrec%status=0
   allocate(addrec%need_property(1))
   call need_propertyid('LNTH ',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(1)=typty
! missing things for the actual Cp function ...
!
   write(kou,*)'Not implemented yet'; gx%bmperr=4078
!
1000 continue
   return
 end subroutine create_debyecp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_debyecp
!\begin{verbatim}
 subroutine calc_debyecp(moded,phres,lokadd,lokph,mc,ceq)
! calculates Mauro Debye contribution
! NOTE: values for function not saved, should be done to save calculation time.
! moded: integer, 0=only G, S, Cp; 1=G and dG/dy; 2=Gm dG/dy and d2G/dy2
! phres: pointer, to phase\_varres record
! lokadd: pointer, to addition record
! lokph: integer, phase record 
! mc: integer, number of constituents
! ceq: pointer, to gtp_equilibrium_data
   implicit none
   integer moded,lokph,mc
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_phase_add), pointer :: lokadd
   TYPE(gtp_phase_varres) :: phres
!\end{verbatim}
   integer ith,noprop
! value of THET and derivatives have type ??
   noprop=phres%listprop(1)-1
!    write(*,*)'3H cmi 2: ',noprop,(phres%listprop(i),i=1,noprop)
! Find thet, index stored in need_property(1)
   do ith=2,noprop
      if(phres%listprop(ith).eq.lokadd%need_property(1)) goto 100
   enddo
   write(*,*)'3H No Debye temperature LNTH',lokph
   gx%bmperr=4336; goto 1000
100 continue
   write(*,*)'3H Deby low T heat capacity model not implemented'
   gx%bmperr=4078
1000 continue
   return
 end subroutine calc_debyecp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine create_diffusion
!\begin{verbatim}
 subroutine create_diffusion(addrec,lokph,text)
   implicit none
   integer lokph
   character text*(*)
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim} %+
   integer typty,jj,last,is,js,ks,loksp,loksp2,ll,nsl
   character typ*24,quest*38,spname*24
   double precision alpha
   type(gtp_diffusion_model), pointer :: diffcoef
   logical once
! initiate
   quest='Dependent constituent in sublattice X:'
! reserve an addition record
   allocate(addrec)
! nullify pointer to next addition
   nullify(addrec%nextadd)
   addrec%status=0
! Set the type of addition and look for needed parameter properties
   addrec%type=diffcoefs
! Some information is needed
   last=1
100 continue
   call gparcdx('Type of diffusion model: ',text,last,1,typ,'SIMPLE',&
        '?Amend diffusion')
   call capson(typ)
!   write(*,*)'3H typ: ',index('MAGNETIC',trim(typ)),trim(typ)
   if(index('SIMPLE',trim(typ)).eq.1) then
      write(*,*)'Simple diffusion  model selected'
      jj=2
   elseif(index('MAGNETIC',trim(typ)).eq.1) then
      write(*,*)'Magnetic diffusion model selected'
      jj=3
   else
      write(*,*)'Dilute diffusion model selected'
      jj=1
   endif
! allocate diffusion record for data
   allocate(addrec%diffcoefs)
   diffcoef=>addrec%diffcoefs
!   addrec%diffcoefs=>diffcoef
! ???????????? must we have a diffusion record for each composition set??
   diffcoef%difftypemodel=jj
   diffcoef%status=0
   nullify(diffcoef%nextcompset)
! dependent component for each sublattice
   nsl=phlista(lokph)%noofsubl
   allocate(diffcoef%depcon(nsl))
   is=1
   do ll=1,nsl
      quest(37:37)=char(ll+ichar('0'))
      once=.true.
200   continue
      loksp=phlista(lokph)%constitlist(is)
      spname=splista(loksp)%symbol
      call gparcdx(quest,text,last,1,typ,spname,'?Amend diffusion')
      call find_species_record(typ,loksp2)
      if(gx%bmperr.ne.0) then
         if(once) then
            once=.false.
            write(*,*)'No such species'
         else
            goto 1000
         endif
      endif
! we must also check this species is a constient in the sublattice!!
      if(loksp2.ne.loksp) then
         do js=is,is+phlista(lokph)%nooffr(ll)-1
            if(loksp2.eq.phlista(lokph)%constitlist(js)) goto 250
         enddo
         write(*,*)'This species is not a constituent of this sublattice'
         if(once) goto 200
         gx%bmperr=4399; goto 1000
      endif
250   continue
! is is always the first constituent in each sublattice
      diffcoef%depcon(ll)=loksp2
      is=is+phlista(lokph)%nooffr(ll)
   enddo
! for jj=3 we must ask for ALPHA and ALPHA2 (with species names)
   if(jj.eq.3) then
      allocate(diffcoef%alpha(phlista(lokph)%nooffr(2)))
      call gparrdx('Value of ALPHA: ',text,last,alpha,0.3D0,'?Amend diffusion')
      diffcoef%alpha(1)=alpha
      if(nsl.eq.2 .and. phlista(lokph)%nooffr(2).gt.1) then
         ks=2
         is=phlista(lokph)%nooffr(1)
         loop: do ll=1,phlista(lokph)%nooffr(2)
            loksp=phlista(lokph)%constitlist(is+ll)
! if constituent is Va ignore!!
            if(.not.btest(splista(loksp)%status,SPVA)) then
               spname=splista(loksp)%symbol
               quest='Value of ALPHA2&'//trim(spname)
               call gparrdx(quest,text,last,alpha,1.0D0,'?Amend diffusion')
               if(ks.le.size(diffcoef%alpha)) diffcoef%alpha(ks)=alpha
               ks=ks+1
            endif
         enddo loop
      endif
!      write(*,*)'3H alpha: ',diffcoef%alpha
   endif
!   write(*,*)'3H depcon: ',diffcoef%depcon
! This addition may use MQ, MF, MG and maybe more
   allocate(addrec%need_property(3))
   call need_propertyid('MQ  ',typty)
   addrec%need_property(1)=typty
   call need_propertyid('MF  ',typty)
   addrec%need_property(2)=typty
   call need_propertyid('MG  ',typty)
   addrec%need_property(3)=typty
   if(gx%bmperr.ne.0) goto 1000
   write(*,*)'Diffusion record created'
!   write(kou,*)'Not implemented yet'; gx%bmperr=4078
1000 continue
   return
 end subroutine create_diffusion

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine diffusion_onoff
!\begin{verbatim}
 subroutine diffusion_onoff(phasetup,bitval)
! switches the bit which calculates diffusion coefficients on/off
! if bitval 0 calculate is turned on, 1 turn off
   implicit none
   integer bitval
   type(gtp_phasetuple) :: phasetup
!\end{verbatim} %+
   integer lokph
   type(gtp_phase_add), pointer :: addrec
   lokph=phasetup%lokph
   addrec=>phlista(lokph)%additions
   loop: do while(associated(addrec))
      if(addrec%type.eq.DIFFCOEFS) then
         if(bitval.eq.0) then
            addrec%diffcoefs%status=ibclr(addrec%diffcoefs%status,0)
         else
            addrec%diffcoefs%status=ibset(addrec%diffcoefs%status,0)
         endif
         exit loop
      endif
   enddo loop
1000 continue
   return
 end subroutine diffusion_onoff

 !/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_diffusion
!\begin{verbatim}
 subroutine calc_diffusion(moded,phres,lokadd,lokph,mc,ceq)
! calculates diffusion coefficients
! NOTE: values for function not saved, should be done to save calculation time.
! moded: integer, 0=only G, S, Cp; 1=G and dG/dy; 2=Gm dG/dy and d2G/dy2
! phres: pointer, to phase\_varres record
! lokadd: pointer, to addition record
! lokph: integer, phase record 
! mc: integer, number of constituents
! ceq: pointer, to gtp_equilibrium_data
   implicit none
   integer moded,lokph,mc
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_phase_add), pointer :: lokadd
   TYPE(gtp_phase_varres) :: phres
!\end{verbatim} %+
   type(gtp_diffusion_model), pointer :: diffcoef
   diffcoef=>lokadd%diffcoefs
!   write(*,*)'Diffusion phase and model: ',trim(phlista(lokph)%name),&
!        diffcoef%difftypemodel
!   write(*,*)'Dependent const: ',diffcoef%depcon
!   write(*,*)'Alpha: ',diffcoef%alpha
!   write(*,*)'Calculation on the diffusion record not implemented'
1000 continue
   return
 end subroutine calc_diffusion

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine get_diffusion_matrix
!\begin{verbatim}
 subroutine get_diffusion_matrix(phtup,mdm,dcval,ceq)
! extracts calculated diffusion coefficients for a phase tuple
! phtup phase tuple
! dcval diffusion matrix
! ceq: pointer, to gtp_equilibrium_data
   implicit none
   integer mdm
   double precision dcval(mdm,*)
   TYPE(gtp_phasetuple) :: phtup
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   TYPE(gtp_phase_add), pointer :: lokadd
   TYPE(gtp_phase_varres) :: phres
1000 continue
   return
 end subroutine get_diffusion_matrix

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_addition
!\begin{verbatim}
 subroutine list_addition(unit,CHTD,phname,ftyp,lokadd)
! list description of an addition for a phase on unit
! used when writing databases files and phase data
   implicit none
   integer unit,ftyp
! CHTD is letter for TDB files TYPE_DEFINITION ... suck
   character CHTD*1,phname*(*)
   TYPE(gtp_phase_add), pointer :: lokadd
!\end{verbatim} %+
   integer ip
   TYPE(tpfun_expression), pointer :: exprot
   character line*256,tps(2)*3,chc*2
   double precision ff
!
   if(unit.eq.kou) then
      chc='  '
   else
      chc='$ '
   endif
!
!   if(.not.btest(lokadd%status,ADDHAVEPAR)) then
! skip additions with no parameters for this phase
! REMOVED THIS because it creates confusion
! If parameters added before addition specified then ADDHAVEPAR is not set
!      write(*,*)chc,'3H No parameters for addition: ',&
!           trim(additioname(lokadd%type))
!      goto 1000
!   else
!      write(*,*)'3H status word for this addition: ',lokadd%status
!   endif
   addition: select case(lokadd%type)
   case default
      write(unit,*)'Unknown addtion type: ',phname,lokadd%type
   case(indenmagnetic) ! Inden magnetic model
      if(ftyp.eq.2) then
! TDB file: I do not think I have saved the enthalpy factor, bcc (-1) it is 0.4
         ff=0.28D0
         if(lokadd%aff.eq.-1) ff=0.4D0
         write(unit,88)CHTD,phname(1:len_trim(phname)),lokadd%aff,ff
88       format(' TYPE_DEFINITION ',a,' GES A_P_D ',a,' MAGNETIC ',i3,F8.4,'!')
      else
         write(unit,100)lokadd%aff
100      format(2x,'+ Magnetic model by Inden, anti-ferromagnetic factor:',i3,/&
              4x,'Magnetic function below the ordering temperature TC',&
              ' with TAO=T/TC:')
         tps(1)='TAO'
         tps(2)='err'
         ip=1
         line=' '
         exprot=>lokadd%explink(1)
         call ct1wfn(exprot,tps,line,ip)
         call wrice(unit,4,8,78,line(1:ip))
         write(unit,110)
110      format(4x,'Magnetic function above the ordering temperature TC ',&
              'with TAO=T/TC:')
         ip=1
         line=' '
         exprot=>lokadd%explink(2)
         call ct1wfn(exprot,tps,line,ip)
         call wrice(unit,4,8,78,line(1:ip))
! write current values of gmagn and values
!         write(unit,120)(lokadd%propval(ip),ip=1,4)
!120      format('    Curr. contrib. G, G.T etc:',4(1pe12.4))
      endif
!---------------------------------------------
   case(debyecp) ! Debye Cp model
      write(unit,200)chc
200   format(a,'+ Debye Cp model, not implemented yet')
!---------------------------------------------
   case(xiongmagnetic) ! Inden-Qing-Xiong
      write(unit,300)chc,lokadd%status
300   format(a,'+ Inden magnetic model modified by Qing and Xiong ',Z8/&
           4x,'with separate Curie and Neel temperatures.'/&
           4x,'Magnetic function below the ordering temperature TC ',&
           ' with TAO=T/TC:')
      tps(1)='TAO'
      tps(2)='err'
      ip=1
      line=' '
      exprot=>lokadd%explink(1)
      call ct1wfn(exprot,tps,line,ip)
      call wrice(unit,4,8,78,line(1:ip))
      write(unit,110)
      ip=1
      line=' '
      exprot=>lokadd%explink(2)
      call ct1wfn(exprot,tps,line,ip)
      call wrice(unit,4,8,78,line(1:ip))
!---------------------------------------------
   case(einsteincp) ! Einstein Cp model
      write(unit,400)chc
400   format(a,'+ Einstein Cp model: 1.5R*exp(LNTH(x)) +',&
           ' 3RT*ln(1-exp(-exp(LNTH(x))/T))')
!---------------------------------------------
   case(elasticmodel1) ! Elastic model 1
      write(unit,500)
500   format(1x,'+ Elastic model 1, with P interpreted as a force in',&
           ' the X direction.')
!---------------------------------------------
   case(twostatemodel1) ! Liquid two state  model including Einstein
      write(unit,510)chc,chc
510   format(a,'+ Liquid 2 state model: G(liq)-RT*ln(1+exp(-G2(x,T)/RT))'/&
           a,'+ Einstein Cp model: 1.5R*exp(LNTH(x)) ',&
           '+ 3RT*ln(1-exp(-exp(LNTH(x))/T))')
!---------------------------------------------
   case(volmod1) ! Volume model 1
      write(unit,520)chc
520   format(a,'+ Volume model P*V0(x)*exp(VA(x,T))')
!---------------------------------------------
!   case(crystalbreakdownmod) ! Crystal breakdown model UNUSED, EET not listed
!      write(unit,530)chc
!530   format(a,'+ Crystal breakdown model used above current value of CBT')
!---------------------------------------------
   case(secondeinstein) ! Second Einstein Cp contribution
      write(unit,540)chc
540   format(a,'+ Second Einstein: DCP2(x)*RT*ln(exp(ln(THT2(x))/T)-1)')
!---------------------------------------------
   case(schottkyanomaly) ! Schottky Anomaly
      write(unit,550)chc
550   format(a,'+ Schottky anomaly DSCH(x)*RT*ln(1+exp(-ln(TSCH(x))/T)) ')
!---------------------------------------------
! THIS MODEL IS OBSOLETE
   case(twostatemodel2) ! Liquid two state  model with fix G2 and Einstein
      write(unit,511)chc,chc
511   format(a,' + wrong Liquid 2 state model: G(liq)-RT*ln(1+exp(-G2(T)/RT))'/&
           a,' + Einstein Cp model: 1.5R*exp(LNTH(x)) ',&
           '+ 3RT*ln(1-exp(exp(LNTH(x))/T))')
!---------------------------------------------
   end select addition
1000 continue
   return
 end subroutine list_addition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine list_addition_values
!\begin{verbatim} %-
 subroutine list_addition_values(unit,phres)
! lists calculated values for this addition
! Used for the command CALCULATE PHASE to inform about additions
   implicit none
   integer unit
   TYPE(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   integer lokph,j1
   TYPE(gtp_phase_add), pointer :: addrec
!
   lokph=phres%phlink
   addrec=>phlista(lokph)%additions
   do while(associated(addrec))
!      write(lut,77)addrec%type,(addrec%propval(j1),j1=1,4)
77    format('Addition type ',i2,': ',4(1pe12.4))
! ignore additions without parameters
      if(btest(addrec%status,ADDHAVEPAR)) &
           write(lut,78)additioname(addrec%type),(addrec%propval(j1),j1=1,4)
78    format('Addition/RT ',a,':',4(1pe10.2))
      addrec=>addrec%nextadd
   enddo
1000 continue
   return
 end subroutine list_addition_values

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine add_ternary_extrapol_method
!\begin{verbatim} %-
 subroutine add_ternary_extrapol_method(kou,lokph,typ,toophead,species)
! add a Toop or Kohler extrapolation method for a ternary subsystem to a phase
! interactive or from database
! if kou>0 then the command is given interactivly and questions can be asked
   implicit none
   integer kou,lokph
   character typ*1,species(3)*24
! VERY DIFFICULT TO FIX: I must use an external sequental list for all tooprec
! and allocate the "seq" pointer
! oterwise the allocated tooprecord were destroyed in an arbitrary way ..
   type(gtp_tooprec) :: toophead
!\end{verbatim}
   integer jj,kk,loksp,tint(3),toop,done,conix
   type(gtp_endmember), pointer :: endmem
   type(gtp_interaction),pointer :: intrec12,intrec13,intrec23,intrec
   type(gtp_interaction),pointer :: intrec1,intrec2
   type(gtp_tooprec), allocatable, target :: newrec
!   type(gtp_tooprec), pointer :: method,duplicate
   type(gtp_tooprec), pointer :: duplicate,next
   type(gtp_tooprec), pointer :: method
   character dummy*24
   logical checkdup
! for debugging
   integer nextmethod
   integer aheremethod12, aheremethod13,aheremethod23
! this is incremented by 1 each time a record is created in any phase
!   integer, save ::uniqid=0
!
!   write(*,8)lokph,typ,trim(species(1)),trim(species(2)),trim(species(3))
8  format('3H Add_ternary_extrapol_method ',i3,2x,a,': ',a,'-',a,'-',a)
   if(phlista(lokph)%noofsubl.gt.1) then
      write(*,*)'3H Kohler/Toop not allowed for phases with sublattices'
      gx%bmperr=4399; goto 1000
   endif
   if(associated(phlista(lokph)%disordered)) then
      write(*,*)'3H Kohler/Toop not allowed for phases with disordered set'
      gx%bmperr=4399; goto 1000
   endif
   if(.not.(typ.eq.'K' .or. typ.eq.'T')) then
      write(*,*)'3H illegal extrapolation method: ',typ
      gx%bmperr=4399; goto 1000
   endif
! find the 3 constituents in the same sublattice in the phase
! they can be in any order
   all3: do jj=1,3
      call find_species_record_noabbr(species(jj),loksp)
      if(gx%bmperr.ne.0) goto 1000
! check if it is a constituent
      isconst: do kk=1,size(phlista(lokph)%constitlist)
!         write(*,'(a,5i5)')'3H constituent?',jj,kk,&
!              phlista(lokph)%constitlist(kk),loksp
         if(loksp.eq.phlista(lokph)%constitlist(kk)) then
! I should check they are in the same sublattice ... but we have just one
            tint(jj)=kk; cycle all3
         endif
      enddo isconst
! a constituent not found
      write(*,*)'3H no such constituent: ',jj,species(jj)
      gx%bmperr=4399; goto 1000
   enddo all3
!   write(*,'(a,3i4)')'3H found constituents: ',tint
   if(tint(1).eq.tint(2) .or. tint(1).eq.tint(3) .or. tint(2).eq.tint(3)) then
      write(*,*)'3H Same element twice, not a ternary!'
      gx%bmperr=4399; goto 1000
   endif
! sort constituents in numerical order (clumsy but not time-critical)
   toop=0
   if(typ.eq.'T') toop=1
! if 1>2 change order
   if(tint(1).gt.tint(2)) then
      jj=tint(1); tint(1)=tint(2); tint(2)=jj
      dummy=species(1); species(1)=species(2); species(2)=dummy
      if(toop.eq.1) toop=2
   endif
! now 1 < 2, check if 2>3?
   if(tint(2).gt.tint(3)) then
      jj=tint(2); tint(2)=tint(3); tint(3)=jj
      dummy=species(2); species(2)=species(3); species(3)=dummy
      if(toop.eq.2) toop=3
   endif
! now 3 > (1,2) check again if 1>2
   if(tint(1).gt.tint(2)) then
      jj=tint(1); tint(1)=tint(2); tint(2)=jj
      dummy=species(1); species(1)=species(2); species(2)=dummy
      if(toop.eq.1) toop=2
   endif
! now 1 < 2 < 3 and species in same order
!   write(*,'(a,3i3)')'3H tint:',tint
! now we have to find the interaction records ... suck
! Some binary parameter in the ternary may not have an interaction record.   
! there must be a check if this ternary is not a duplicate
! disordered fraction set not allowed
! Try to give reasonable error messages
   endmem=>phlista(lokph)%ordered
   if(.not.associated(endmem)) then
      write(*,*)'3H No endmembers or interaction in this phase!'
      gx%bmperr=4399; goto 1000
   endif
! nullify pointers from binary interactions
   nullify(intrec1); nullify(intrec2)
   nullify(intrec12); nullify(intrec13); nullify(intrec23)
! look for endmember with lowest constituent index (they are ordered that way)
! REMEMBER we have a substitutional model, no sublattices
! one may add a second sublattice with same constituent in all 3 endmembers
   start: do conix=1,2
!      write(*,'(a,3i3)')'3H endmem1: ',conix,endmem%fraclinks(1,1),tint(conix)
      findem: do while(endmem%fraclinks(1,1).lt.tint(conix))
!      write(*,'(a,3i3)')'3H endmem2: ',conix,endmem%fraclinks(1,1),tint(conix)
         endmem=>endmem%nextem;
         if(associated(endmem)) then
            cycle findem
         elseif(associated(intrec1) .or. associated(intrec2)) then
            exit start
         else
            write(*,*)'3H no binary model parameters for this ternary!'
            gx%bmperr=4399; goto 1000
         endif
      enddo findem
      if(.not.associated(endmem)) exit start
!      write(*,'(a,3i3)')'3H endmem3: ',conix,endmem%fraclinks(1,1),tint(conix)
! we have found an endmember >=tint(conix)
      if(endmem%fraclinks(1,1).eq.tint(conix)) then
!         write(*,*)'3H Found endmember for ',species(conix)
         if(conix.eq.1) then
! interaction AB and AC are in same list of interactions, separate later
            write(*,*)'3H interaction list 1-2 and 1-3 found'
            intrec1=>endmem%intpointer
!            if(associated(intrec1%tooprec)) write(*,'(a,4i3)')'3H link 1 id:',&
!                 intrec1%fraclink(1),intrec1%tooprec%uniqid
         else
            write(*,*)'3H interaction list 2-3 found'
            intrec2=>endmem%intpointer
!            if(associated(intrec2%tooprec)) write(*,'(a,4i3)')'3H link 2 id:',&
!                 intrec1%fraclink(1),intrec1%tooprec%uniqid
         endif
      else
         write(*,*)'3H No endmember for constituent:',trim(species(conix))
      endif
   enddo start
   if(.not.associated(intrec1) .and. .not.associated(intrec2)) then
      write(*,*)'3H Ignoring extrapolation method as no binary interactions'
      gx%bmperr=4399; goto 1000
   endif
!   if(associated(intrec1%tooprec)) write(*,'(a,4i3)')'3H link 3 id:',&
!        intrec1%fraclink(1),intrec1%tooprec%uniqid
!   if(associated(intrec2%tooprec)) write(*,'(a,4i3)')'3H link 4 id:',&
!        intrec2%fraclink(1),intrec2%tooprec%uniqid
! We have set intrec1 to first interaction record of A and B
! Look for interaction record AB, AC and BC in these interaction lists
!--------------------------------------
! we have a link to an interaction record for tint(1) or tint(2)
100 continue
! there is no order of the interaction constituent
! one or two binary interaction parameter can be missing but not all three!!
! Remember: there is only a single set of sites
!   checkdup=.TRUE.
! skip check for duplicates
!   checkdup=.FALSE.
! We look for interaction records 1-2, 1-3 and 1-3
   write(*,*)'3H intrec: ',associated(intrec1),associated(intrec2)
   intrec=>intrec1
   done=tint(1)
!   if(associated(intrec%tooprec)) write(*,'(a,4i3)')'3H link 5 id:',&
!        intrec%fraclink(1),intrec%tooprec%uniqid
! Jump back here for next interaction list. note intrec13 set from intrec1
200 continue
   intwith1: do while(associated(intrec))
!      if(associated(intrec%tooprec)) write(*,'(a,4i3)')'3H link 6 id:',&
!           done,intrec%fraclink(1),intrec%tooprec%uniqid
! The endmember record is either tint(1) or tint(2)
      if(intrec%fraclink(1).eq.tint(2)) then
! In this case the endmember record must be tint(1)
         intrec12=>intrec
         write(*,*)'3H Found interaction 1-2: ',trim(species(1)),'-',&
              trim(species(2)),associated(intrec12%tooprec)
      elseif(intrec%fraclink(1).eq.tint(3)) then
! The endmember record can be tint(1) or tint(2)
         if(done.eq.tint(1)) then
            intrec13=>intrec
            write(*,*)'3H Found interaction 1-3: ',trim(species(1)),'-',&
                 trim(species(3)),associated(intrec13%tooprec)
         else
            intrec23=>intrec
            write(*,*)'3H Found interaction 2-3: ',trim(species(2)),'-',&
                 trim(species(3)),associated(intrec23%tooprec)
            exit intwith1
         endif
      endif
      intrec=>intrec%nextlink
   enddo intwith1
   if(done.eq.tint(1)) then
 ! we have searched the whole ternary list from intrec12 and intrec13?
     done=tint(2); intrec=>intrec2; goto 200
   endif
!--------------------------------------------------------------
! This seems to work: have a global list of tooprecords and allocate the 
! toophead%seq record to use when linking together the ternaries.
! This sequential list of tooprecord is not used for anything.
! We have not really linked the method yet, only found the binaries where we
! should find the records to set the links to the tooprecord.
! we have now links from the 3 binary interaction records (same may be null)
! Link the new record from these and copy the previos method records to
! the new method record.
   allocate(toophead%seq)
   method=>toophead%seq
   method%toop=toop
   method%const1=tint(1)
   method%const2=tint(2)
   method%const3=tint(3)
   nullify(method%next12); nullify(method%next13); nullify(method%next23)
   method%extra=0
   uniqid=uniqid+1
   method%uniqid=uniqid
! problem to insert this record in the structure ...
! seach the end of all lists of methods to add the new one
!   write(*,97)uniqid,toop,tint,(trim(species(jj)),jj=1,3),&
!        associated(intrec12),associated(intrec13),associated(intrec23)
97 format('3H created method record: ',2i3,2x,3i3," : ",a,1x,a,1x,a,2x,3l2)
!=====================================================================
! search for the last tooprec for each binary and check for duplicates
! if no tooprec records then link from intrec12%tooprec
! New method A-B-C, we have found interaction records (if null no binary param)
! intrec12: a list with methods A-B-X, A-X-B or X-A-B. Duplicate if A-B-C exist
!               next links can be A-B-X, A-X-B or X-A-B 
! intrec13: a list with methods A-C-X or A-X-C
! intrec23: a list with methods B-C-X and B-X-C
   if(.not.associated(intrec12)) then
! There is no binary interaction record, no place to store a link
!      write(*,103)trim(species(1)),trim(species(2))
103   format('3H No interaction record',a,'-',a)
   else
!---------------------------------------------------------------- A-B
! Linking from binary A-B
      duplicate=>intrec12%tooprec
      nextmethod=0
!      write(*,105)' 1-2 ',associated(duplicate),nextmethod,&
!           trim(species(1)),trim(species(2))
105   format('3H linking ',a,l2,i2,': ',a,'-',a)
      if(associated(duplicate)) then
! The constituents in each method records are always ordered alphabetically
! %const1<%const2<%const3 X<Y<Z and tint(1)<tint(2)<tint(3), A<B<C
! X can be %const1 or %const2 and Y can be %const2 or %const3
! Check also for duplicates
         next12: do while(associated(duplicate))
! save link to duplicate
            next=>duplicate
!            write(*,107)'3H The intrec12 list: ',trim(species(1)),&
!                 trim(species(2)),trim(species(3)),duplicate%uniqid,&
!                 tint(1),tint(2),&
!                 duplicate%const1,duplicate%const2,duplicate%const3,next%uniqid
107         format(a,a,'-',a,' with ',a,i3,3x,2i3,2x,3i3,5x,i3)
! In this method A<B can be any combination of %const1<%const2<%const3
! BUT THE THEY ARE ALL IN INCREASING ORDER 
! %const1 can be less than tint(1) but if so %const2 must be equal to tint(1)
! we also have to check if the new is a duplicate and which next link to follow
            if(duplicate%const1.eq.tint(1)) then
! When %const1=tint(1) then tint(2) can be %const2 or %const3
               if(duplicate%const2.eq.tint(2)) then
! Method A-B-X, Check if this is a duplicate
                  if(duplicate%const3.eq.tint(3)) then
                     write(*,*)'3H found duplicate 1'; goto 1100
                  endif
! The A-B is the first binay, use the first link, %cont3 is not tint(3)
                  nextmethod=1
                  duplicate=>duplicate%next12
               elseif(duplicate%const3.eq.tint(2)) then
! Method A-X-B with A-B as the second binary, X cannot be tint(3) as higher
                  nextmethod=2
                  duplicate=>duplicate%next13
               else
! ERROR as no binary A-B in this method 
                  write(*,*)'3H exit 1'
                  goto 1200
                  exit next12
               endif
            elseif(duplicate%const2.eq.tint(1) .and.&
                 duplicate%const3.eq.tint(2)) then
! Metod X-A-B with A-B as the third binary, X cannot be tint(3)
               nextmethod=3
               duplicate=>duplicate%next23
!               write(*,*)'3H take next23 link ',associated(duplicate)
            else
!               write(*,130)duplicate%uniqid,&
!                    duplicate%const2,tint(2),duplicate%const3,tint(3)
130            format('3H not found? ',i3,5x,2i4,2x,2i4)
! ERROR as no binary A-B in this method 
               write(*,*)'3H exit 2'
               goto 1200
               exit next12
            endif
            write(*,*)'3H loop duplicate: ',associated(next),nextmethod
         enddo next12
      endif
! we come here at the end of the list of methods including the A-B binary
! This is the first method for the binary A-B, note duplicate is null!
      if(nextmethod.eq.0) then
         intrec12%tooprec=>method
!         write(*,120)trim(species(1)),trim(species(2))
120      format('3H New method linked from interacion: ',a,'-'a)
      else
! insert the correct link, note duplicate is null !!
         if(nextmethod.eq.1) next%next12=>method
         if(nextmethod.eq.2) next%next13=>method
         if(nextmethod.eq.3) next%next23=>method
!         write(*,121)next%uniqid,nextmethod,next%const1,next%const2,next%const3
121      format('3H New method linked 1B:',i3,2x,i3,2x,3i3)
      endif
   endif
!   write(*,*)'3H linked method for 1-2'
!----------------------------------------------- A-C, tint(1)-tint(3)
! A method here X<Y<Z can have A-C as any binary
   if(.not.associated(intrec13)) then
!      write(*,103)trim(species(1)),trim(species(3))
   else
      duplicate=>intrec13%tooprec
      nextmethod=0
!      write(*,105)' 1-3 ',associated(duplicate),nextmethod,&
!           trim(species(1)),trim(species(3))
      if(associated(duplicate)) then
! The constituents in each method records are always ordered alphabetically
! %const1<%const2<%const3 X<Y<Z and tint(1)<tint(2)<tint(3), A<B<C
! X can be %const1 or %const2 and Y can be %const2 or %const3
! Check also for duplicates
         next13: do while(associated(duplicate))
            next=>duplicate
!            write(*,107)'3H The intrec13 list: ',trim(species(1)),&
!                 trim(species(2)),trim(species(3)),duplicate%uniqid
! In this method A<B can be any combination of %const1<%const2<%const3
! BUT THE THEY ARE ALL IN INCREASING ORDER 
! %const1 can be less than tint(1) but if so %const2 must be equal to tint(1)
! we also have to check if this is duplicate or which next link to follow
            if(tint(1).eq.duplicate%const1) then
! When %const1=tint(1) tint(3) must be %cont3
               if(duplicate%const2.eq.tint(3)) then
! Method A-C-X, no duplicate as tint(2)<tint(3) 
!                  if(duplicate%const3.eq.tint(3)) write(*,*)'3H duplicate 1'
!                  if(tint(2).eq.duplicate%const3) goto 1100
! The A-C is the first binary, use the first link
                  nextmethod=1
                  duplicate=>duplicate%next12
               elseif(duplicate%const3.eq.tint(3)) then
! Method A-X-C, no duplicate possible
                  nextmethod=2
                  duplicate=>duplicate%next13
               else
! Error as there is no A-C binary
                  write(*,*)'3H exit 3'  
                  goto 1200
                  exit next13
               endif
            elseif(duplicate%const2.eq.tint(1) .and.&
                 duplicate%const3.eq.tint(2)) then
               nextmethod=3
               duplicate=>duplicate%next23
            else
! No A-C in this method, it does no fit in this list!
               write(*,*)'3H exit 4'
               goto 1200
               exit next13
            endif
         enddo next13
      endif
! we come here at the end of the list of methods including the A-B binary
! This is the first method for the binary A-B
      if(nextmethod.eq.0) then
         intrec13%tooprec=>method
!         write(*,120)trim(species(1)),trim(species(3))
      else
! insert the correct link
         if(nextmethod.eq.1) next%next12=>method
         if(nextmethod.eq.2) next%next13=>method
         if(nextmethod.eq.3) next%next23=>method
!         write(*,121)next%uniqid,nextmethod,next%const1,next%const2,next%const3
      endif
   endif
!   write(*,*)'3H linked method for 1-3'
!----------------------------------------------- B-C list
   if(.not.associated(intrec23)) then
!      write(*,103)trim(species(2)),trim(species(3))
   else
      duplicate=>intrec23%tooprec
      nextmethod=0
!      write(*,105)' 2-3 ',associated(duplicate),nextmethod,&
!           trim(species(2)),trim(species(3))
      if(associated(duplicate)) then
! The constituents in each method records are always ordered alphabetically
! %const1<%const2<%const3 X<Y<Z and tint(1)<tint(2)<tint(3), A<B<C
! X can be %const1 or %const2 and Y can be %const2 or %const3
! Check also for duplicates
         next23: do while(associated(duplicate))
            next=>duplicate
!            write(*,107)'3H The intrec23 list: ',trim(species(1)),&
!                 trim(species(2)),trim(species(3)),duplicate%uniqid
! In this method A<B can be any combination of %const1<%const2<%const3
! BUT THE THEY ARE ALL IN INCREASING ORDER 
! %const1 can be less than tint(1) but if so %const2 must be equal to tint(1)
! we also have to check if this is duplicate or which next link to follow
            if(duplicate%const1.eq.tint(1)) then
               if(duplicate%const2.eq.tint(2)) then
! Method A-C-X, Check if this is a duplicate is redundant 
                  if(duplicate%const3.eq.tint(3)) write(*,*)'3H duplicate 2'
                  if(tint(3).eq.duplicate%const3) goto 1100
! The A-C is the first binay, use the first link
                  nextmethod=1
                  duplicate=>duplicate%next12
               elseif(tint(3).eq.duplicate%const3) then
! Method B-X-C Check for duplicate is redundant as tint(2)<tint(3)=%const3
!               if(tint(2).eq.duplicate%const3) goto 1100
! The A-C is the second binary, use the second link
                  nextmethod=2
                  duplicate=>duplicate%next13
               else
                  write(*,*)'3H exit 5'
                  goto 1200
                  exit next23
               endif
            elseif(tint(1).eq.duplicate%const2 .and.&
                 tint(2).eq.duplicate%const3) then
! Method X-A-C, duplicate check probably redundant
               if(duplicate%const1.eq.tint(1)) write(*,*)'3H duplicate 3'
               if(tint(1).eq.duplicate%const1) goto 1100
               nextmethod=3
               duplicate=>duplicate%next23
            else
! No A-C in this method does no fit in this list!
               write(*,*)'3H exit 6'
               goto 1200
               exit next23
            endif
         enddo next23
      endif
! we come here at the end of the list of methods including the A-B binary
! This is the first method for the binary A-B
      if(nextmethod.eq.0) then
         intrec23%tooprec=>method
!         write(*,120)trim(species(2)),trim(species(3))
      else
! insert the correct link
         if(nextmethod.eq.1) next%next12=>method
         if(nextmethod.eq.2) next%next13=>method
         if(nextmethod.eq.3) next%next23=>method
!         write(*,121)next%uniqid,nextmethod,next%const1,next%const2,next%const3
      endif
   endif
!   write(*,*)'3H linked method for 2-3'
!===================================================================
! Puuuuuuuuuuuuuuuuuhhhhhhhhhhhhhhhhhhhh
!   write(*,510)'3H Toop/Kohler record:',method%uniqid,method%toop,&
!        method%const1,method%const2,method%const3,method%extra,&
!        associated(method%next12),associated(method%next13),&
!        associated(method%next23)
510 format(a,6i4,3L2)
1000 continue
   return
! Error: Found duplicate method
1100 continue
   write(*,1110)duplicate%uniqid,duplicate%const1,duplicate%const2,&
        duplicate%const3,tint(1),tint(2),tint(3)
1110 format('3H Error found duplicate method ',i3,5x,3i3,5x,3i3)
   gx%bmperr=4399; goto 1000
! Error: Trying to enter a method with wrong set of constituents
1200 continue
   write(*,1210)trim(species(1)),trim(species(2)),trim(species(3)),tint,&
        duplicate%const1,duplicate%const2,duplicate%const3
1210 format('3H Error: ternary system with ',a,'-',a,'-',a,': ',3i3,&
          ' does not fit method: ',3i3)
 end subroutine add_ternary_extrapol_method
 
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

