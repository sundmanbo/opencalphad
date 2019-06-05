!
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
!  integer, public, parameter :: TWOSTATESMODEL1=5
!  integer, public, parameter :: ELASTICMODEL1=6
!  integer, public, parameter :: VOLMOD1=7
!  integer, public, parameter :: UNUSED_CRYSTBREAKDOWNMOD=8
!  integer, public, parameter :: SECONDEINSTEIN=9
!  integer, public, parameter :: SCHOTTKYANOMALITY=10
!  integer, public, parameter :: DIFFCOEFS=11
!------------------------------------
! For each addition XX there is a subroutine create_XX
! called from the add_addrecord
! and a subroutine calc_XX 
! called from the addition_selector, called from calcg_internal
! There is a common list routine
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

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
      write(kou,*)'No such addition type ',addrec%type,lokph
      gx%bmperr=4330
! !
   case(indenmagnetic) ! Inden magnetic 
      addrec%propval=zero
      call calc_magnetic_inden(moded,phres,addrec,lokph,mc,ceq)
! 2
   case(debyecp) ! Debye Cp
      addrec%propval=zero
      call calc_debyecp(moded,phres,addrec,lokph,mc,ceq)
      write(kou,*)' Debye Cp model not implemented yet'
      gx%bmperr=4331
! 3
   case(xiongmagnetic) ! Inden-Qing-Xiong
      addrec%propval=zero
      call calc_xiongmagnetic(moded,phres,addrec,lokph,mc,ceq)
!     write(kou,*)'Inden-Qing-Xiong magn model not implemented yet'
!      gx%bmperr=4332
! 4
   case(einsteincp) ! Einstein Cp
      addrec%propval=zero
      call calc_einsteincp(moded,phres,addrec,lokph,mc,ceq)
!      gx%bmperr=4331
! 5
   case(elasticmodel1) ! Elastic model !
      addrec%propval=zero
      call calc_elastica(moded,phres,addrec,lokph,mc,ceq)
      write(kou,*)' Elastic model not implemented yet'
      gx%bmperr=4399
! 6
   case(twostatemodel1) ! Two state model
      addrec%propval=zero
      call calc_twostate_model1(moded,phres,addrec,lokph,mc,ceq)
! 7
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
   case(schottkyanomality) ! Adding a second Schottky anomality Cp
      addrec%propval=zero
      call calc_schottky_anomality(moded,phres,addrec,lokph,mc,ceq)
! 11
   case(diffcoefs) ! Calculating diffusion coefficients
      addrec%propval=zero
      call calc_diffusion(moded,phres,addrec,lokph,mc,ceq)
!      gx%bmperr=4333
   end select addition
1000 continue
   return
 end subroutine addition_selector

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

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
   logical bcc
!
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
! lokph because we need to check if average or individual Boghr magnetons
      call create_xiongmagnetic(newadd,lokph,bcc)
!-----------------------------------------
   case(debyecp) ! Debye Cp
! 3
      call create_debyecp(newadd)
!-----------------------------------------
   case(einsteincp) ! Einstein Cp
! 4
      call create_einsteincp(newadd)
!-----------------------------------------
   case(elasticmodel1) ! Elastic model 1
! 5
      call create_elastic_model_a(newadd)
!-----------------------------------------
   case(twostatemodel1) ! Liquid 2 state model
! 6
      call create_twostate_model1(newadd)
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
   case(schottkyanomality) ! Schottky anomality
! 10
      call create_schottky_anomality(newadd)
!-----------------------------------------
   case(diffcoefs)  !  diffusion coefficients
! 11 
      call create_diffusion(newadd,lokph,extra)
!-----------------------------------------
   end select addition
!-----------------------------------------
   if(gx%bmperr.ne.0) goto 1000
! initiate status word for this addition
   newadd%status=0
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

!\begin{verbatim}
 subroutine need_propertyid(id,typty)
! get the index of the property needed
   implicit none
   integer typty
   character*4 id
!\end{verbatim}
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
   integer itc,ibm,jl,noprop,ik,k,jk,j
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
            phres%d2gval(ixsym(ik,jk),itc)=&
                 iafftc*phres%d2gval(ixsym(ik,jk),itc)
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
            phres%d2gval(ixsym(ik,jk),ibm)=&
                 iaffbm*phres%d2gval(ixsym(ik,jk),ibm)
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
         d2tao(ixsym(j,k))=&
              2.0*tao*phres%dgval(1,j,itc)*phres%dgval(1,k,itc)/tc**2&
              -tao*phres%d2gval(ixsym(j,k),itc)/tc
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
         d2addgval(ixsym(j,k))=rt*ftao(4)*dtao(1,j)*dtao(1,k)*logb1+&
              rt*ftao(2)*d2tao(ixsym(j,k))*logb1+&
              rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm)+&
              rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm)-&
              rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm)+&
              rt*ftao(1)*invb1*phres%d2gval(ixsym(j,k),ibm)
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
      write(*,'(a,i4,l2,1pe12.4)')'3H msize magadd 1: ',lokph,addpermole,msize
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
         phres%d2gval(ixsym(j,k),1)=phres%d2gval(ixsym(j,k),1)+&
              msize*d2addgval(ixsym(j,k))/rt
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

!\begin{verbatim}
 subroutine create_xiongmagnetic(addrec,lokph,bcc)
! adds a Xiong type magnetic record, we must separate fcc and bcc by extra
! copied from Inden magnetic model
! The difference is that it uses TCA for Curie temperature and TNA for Neel
! and individual Bohr magneton numbers
   implicit none
   logical bcc
   integer lokph
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim} %+
   integer typty,ip,nc
   character text*128
   integer, parameter :: ncc=6
   double precision coeff(ncc)
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
   if(bcc) then
! Magnetic function below Curie/Neel Temperature, 
! problem in ct1xfn to start a function with +1 or 1
      text=' +1-.880323235*T**(-1)-.152870878*T**3-.00679426123*T**9'//&
           '-.00152870878*T**15-5.67238878E-04*T**21'
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
!       write(*,*)'3H emm 2: ',trim(text)
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
      call ct1mexpr(nc,coeff,koder,lhigh)
      if(gx%bmperr.ne.0) goto 1000

   else
!------------
! fcc
! Magnetic function below Curie/Neel Temperature
      text='+1-.842849633*T**(-1)-.174242226*T**3-.00774409892*T**9'//&
           '-.00174242226*T**15-6.46538871E-04*T**21'
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
      call ct1mexpr(nc,coeff,koder,llow)
      if(gx%bmperr.ne.0) goto 1000
! Magnetic function above Curie/Neel Temperature
      text='-.0261039233*T**(-7)-8.70130777E-04*T**(-21)'//&
           '-1.84262988E-04*T**(-35)-6.65916411E-05*T**(-49)'
      ip=1
      nc=ncc
      call ct1xfn(text,ip,nc,coeff,koder,.FALSE.)
      if(gx%bmperr.ne.0) goto 1000
      call ct1mexpr(nc,coeff,koder,lhigh)
      if(gx%bmperr.ne.0) goto 1000
   endif
! reserve an addition record
   allocate(addrec)
! store data in record
   allocate(addrec%explink(2))
   nullify(addrec%nextadd)
   addrec%type=xiongmagnetic
! beware of segmentation fault here !!! llow and llhigh no longer pointers
   addrec%explink(1)=llow
   addrec%explink(2)=lhigh
   addrecs=addrecs+1
   allocate(addrec%need_property(3))
   addrec%addrecno=addrecs
! here the property list is searched for CTA, NTA and IBM
   call need_propertyid('CTA ',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(1)=typty
   if(btest(phlista(lokph)%status1,PHBMAV)) then
! This model can use an effective Bohr magneton number b*=prod(b_i+1)**x_i -1
      call need_propertyid('BMAG ',typty)
   else
! or an individual Bohr magneton number b*=prod(b_i+1)**x_i -1
      call need_propertyid('IBM ',typty)
   endif
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
1000 continue
   return
 end subroutine create_xiongmagnetic

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
   integer itc,itn,ibm,jl,noprop,ik,k,jk,j
   double precision logb1,invb1,iafftc,iaffbm,rgasm,rt,tao,gmagn,msize
   double precision dtaodt,dtaodp,beta,d2taodp2,d2taodtdp,tc,tv
   double precision tao2(2),ftao(6),dtao(3,mc),d2tao(mc*(mc+1)/2)
   double precision addgval(6),daddgval(3,mc),d2addgval(mc*(mc+1)/2)
   double precision tn,tcsave,tnsave
   logical addpermole
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
! record, i.e. typty=2 for TC and typty=3 for BM
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
         ibm=jl
      elseif(phres%listprop(jl).eq.lokadd%need_property(3)) then
         itn=jl
      endif
   enddo findix
! check that the needed properties are defined
   if(ibm.eq.0 .or. (itc.eq.0 .and. itn.eq.0)) then
! it is no error if no CTA, NTA or BMAG but then magnetic contribution is zero
       write(*,12)trim(phlista(lokph)%name)
12     format('3H Warning: Magnetic addition for phase ',a,&
            ' not calculated as '/10x,&
            'values for CTA, NTA or BMAG, magnetic G are zero')
      goto 1000
   else
      tc=-one
      tn=-one
      if(itc.gt.0) tc=phres%gval(1,itc)
      if(itn.gt.0) tn=phres%gval(1,itn)
   endif
! I am not sure I calculate correct derivatives for indivudal Bohr magnetons ...
   if(btest(phlista(lokph)%status1,PHBMAV)) then
      write(*,*)'3H *** Warning: individual Bohjr magneton number not checked'
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
! we should use the appropriate tao=t/tc or t/tn
! use AF model unless tc negative, both cannot be negative here
   if(tc.le.zero) then
! no ferro but antiferro.  tn>0 as both tn and tc checked against zero above
      tcsave=tc
      tc=tn
! we use this index below to extract its value
      itc=itn
! if tn negative use tc
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
   else
      exprot=>lokadd%explink(2)
   endif
   call ct1efn(exprot,tao2,ftao,ceq%eq_tpres)
   logb1=log(beta+one)
   invb1=one/(beta+one)
   gmagn=rt*ftao(1)*logb1
!    write(*,98)'3H cm 97: ',tc,beta,ftao(1),logb1,rt
!    write(*,98)'3H cm 98: ',rt*gmagn,rt*(gmagn+phres%gval(1,1)),tcx,iafftc
!98  format(a,5(1PE14.6))
!
   dtaodt=one/tc
   dtaodp=-tao/tc*phres%gval(3,itc)
   addgval(1)=gmagn
   addgval(2)=gmagn/tv+rt*ftao(2)*dtaodt*logb1
   addgval(3)=rt*ftao(2)*dtaodp*logb1+rt*ftao(1)*invb1*phres%gval(3,ibm)
!   phres%gval(1,1)=phres%gval(1,1)+addgval(1)/rt
!   phres%gval(2,1)=phres%gval(2,1)+addgval(2)/rt
!   phres%gval(3,1)=phres%gval(3,1)+addgval(3)/rt
! save these in record
! NOTE if parallel calculation the same stored values %propval will be
! written by all threads so they must not be used!!
! They are included only for listing and debugging
   do j=1,3
      lokadd%propval(j)=addgval(j)
      phres%gval(j,1)=phres%gval(j,1)+addgval(j)/rt
   enddo
!   write(*,77)lokadd%type,(lokadd%propval(j),j=1,4)
77 format('3H addition ',i2,': ',4(1pe12.4))
! ignore second derivatives if no derivatives wanted
   if(moded.eq.0) then
! make sure Cp is calculated and stored so it can be listed
      addgval(4)=2.0d0*rgasm*ftao(2)*dtaodt*logb1+&
           rt*ftao(4)*(dtaodt)**2*logb1
      lokadd%propval(4)=addgval(4)
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
         d2tao(ixsym(j,k))=&
              2.0*tao*phres%dgval(1,j,itc)*phres%dgval(1,k,itc)/tc**2&
              -tao*phres%d2gval(ixsym(j,k),itc)/tc
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
         d2addgval(ixsym(j,k))=rt*ftao(4)*dtao(1,j)*dtao(1,k)*logb1+&
              rt*ftao(2)*d2tao(ixsym(j,k))*logb1+&
              rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm)+&
              rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm)-&
              rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm)+&
              rt*ftao(1)*invb1*phres%d2gval(ixsym(j,k),ibm)
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
      write(*,'(a,i4,l2,1pe12.4)')'3H msize magadd 2: ',lokph,addpermole,msize
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
         phres%d2gval(ixsym(j,k),1)=phres%d2gval(ixsym(j,k),1)+&
              msize*d2addgval(ixsym(j,k))/rt
      enddo
   enddo
!    write(*,*)'3H cm 7: ',rt*phres%gval(1,1),addgval(1)
! note phres%gval(1..3,1) already calculated above
   do j=4,6
      lokadd%propval(j)=msize*addgval(j)
      phres%gval(j,1)=phres%gval(j,1)+msize*addgval(j)/rt
   enddo
! we may have destroyed the original value of tc if we have AFM
   tc=tcsave
! jump here if no magnetic contribution
1000 continue
   return
 end subroutine calc_xiongmagnetic

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
! store zero in 6 values for propval
   addrec%propval=zero
1000 continue
!   write(*,*)'3H created volume addition',addrecs
   return
 end subroutine create_volmod1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
! store some property values
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
   call need_propertyid('THET',typty)
   if(gx%bmperr.ne.0) goto 1000
   allocate(newadd%need_property(1))
   newadd%need_property(1)=typty
   nullify(newadd%nextadd)
1000 continue
   return
 end subroutine create_einsteincp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
   integer ith,noprop,extreme,j1
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1,fact
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2
   double precision gein,dgeindt,d2geindt2,msize
   logical addpermole
!
   noprop=phres%listprop(1)-1
!   write(*,*)'3H thet: ',phres%listprop(2),addrec%need_property(1)
   findix: do ith=2,noprop
      if(phres%listprop(ith).eq.addrec%need_property(1)) goto 100
   enddo findix
   write(*,*)'3H No value of THET for phase ',trim(phlista(lokph)%name)
   gx%bmperr=4336; goto 1000
100 continue
   if(phres%gval(1,ith).le.one) then
      write(*,69)'3H Illegal THET for phase ',trim(phlista(lokph)%name),&
           phres%gval(1,ith)
69    format(a,a,1pe12.4)
      gx%bmperr=4399; goto 1000
   endif
! NOTE the parameter value is ln(thera)! take the exponential!
! thet is in gval(1,ith), derivatives in dgval(*,ith,*) and d2gval(ith,*)
! G/RT = 1.5*THET/T + 3*R*LN(exp(THET/T) - 1) 
! NOTE ALL VALUES CALCULATED AS FOR G/RT
! kvot=theta/T
   if(phres%gval(1,ith).gt.1.0D2) then
      write(*,*)'Most likely wrong value of THET, parameter should be ln(THET)'
      gx%bmperr=4399; goto 1000
   endif
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
   gein=1.5D0*kvot+3.0D0*ln1mexpkvot
!   write(*,71)'3H Cp E1:',extreme,ceq%tpval(1),gein,ln1mexpkvot,expmkvot,&
!        kvotexpkvotm1
! first derivative wrt T taking care of overflow
   dgeindt=3.0D0*(ln1mexpkvot-kvotexpkvotm1)/ceq%tpval(1)
! This is d2G/dT**2/(RT) = -T**2/R*(Einstein Cp/RT) (or rather Cv/RT)
   if(extreme.eq.-1) then
! take care of overflow at low T, kvotexpkvotm1=expmkvot=0 set above
      d2geindt2=zero
   else
      d2geindt2=-3.0D0*kvotexpkvotm1**2/(expmkvot*ceq%tpval(1)**2)
   endif
! first derivative for each constituent. The parameter value is ln(theta)
! and we should divide by RT
   fact=1.5D0*kvot+3.0D0*kvotexpkvotm1
   do j1=1,mc
      phres%dgval(1,j1,1)=phres%dgval(1,j1,1)+fact*phres%dgval(1,j1,ith)
   enddo
! NOTE if addpermole bit set we have to multiply with derivatives of
! the size of the phase ...
   if(btest(addrec%status,ADDPERMOL)) then
      addpermole=.TRUE.; msize=phres%abnorm(1)
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize lowT: ',lokph,addpermole,msize
   else
      addpermole=.FALSE.; msize=one
!      write(*,'(a,i4,l2,1pe12.4)')'3H msize lowT: ',lokph,addpermole,msize
   endif
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
!
1000 continue
   return
 end subroutine calc_einsteincp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine create_schottky_anomality(newadd)
! Adding a Schottky anomality to Cp
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
! Schottky anomality uses THT2 and DCP2, same as second Einstein
   newadd%type=schottkyanomality
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
 end subroutine create_schottky_anomality

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_schottky_anomality(moded,phres,addrec,lokph,mc,ceq)
! Calculate the contibution due to a Schottky anomality
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
!   write(*,*)'3H thet: ',phres%listprop(2),addrec%need_property(1)
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
!      write(*,*)'3H missing Schottky anomality parameter for phase ',&
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
 end subroutine calc_schottky_anomality

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
      write(*,70)'3H Illegal THET for phase ',trim(phlista(lokph)%name),&
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
! cehck if addition is per mole 
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

!\begin{verbatim}
 subroutine create_twostate_model1(addrec)
! newadd is location where pointer to new addition record should be stored
   implicit none
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim}
   integer typty
! this is bad programming as it cannot be deallocated but it will never be ...
! maybe pointers can be deallocated?
   allocate(addrec)
! nullify pointer to next addition
   nullify(addrec%nextadd)
!-----------------------------
! The model consists of two contributions
! The first is the harmonic vibrations of an ideal amprthous phase
!     this requires a THETA representing the Einstein T
! The second is a term - RT*(1+exp(G2/RT))
! which represent the change from "solid like" to "liquid like"
!-----------------------------
! I am not sure what is is used for
   addrecs=addrecs+1
   addrec%addrecno=addrecs
! property needed
   allocate(addrec%need_property(2))
   call need_propertyid('G2  ',typty)
   addrec%need_property(1)=typty
   call need_propertyid('THETA  ',typty)
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

!\begin{verbatim}
 subroutine calc_twostate_model_john(moded,phres,addrec,lokph,mc,ceq)
! subroutine calc_twostate_model1(moded,phres,addrec,lokph,mc,ceq)
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
   integer jj,noprop,ig2,ith,extreme,jth,kth
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2
   double precision gein,dgeindt,d2geindt2
   double precision g2ein,dg2eindt,d2g2eindt2,theta2,dcpl
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1
   double precision g2val,dg2,expg2,expmg2,rt,tv,rg,dg2dt,dgfdt,d2g2dt2
! This is Johns original model
! number of properties calculatied
   noprop=phres%listprop(1)-1
! locate the THET and G2 property record 
   ig2=0
   ith=0
   jth=0
   kth=0
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
      write(*,*)'Cannot find value for amorphous THET'
      gein=zero; dgeindt=zero; d2geindt2=zero; goto 300
!      gx%bmperr=4399; goto 1000
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
!   addrec%propval(1)=gein
!   addrec%propval(2)=dgeindt
!   addrec%propval(4)=d2geindt2
!   write(*,71)'3H Cp E3: ',extreme,ceq%tpval(1),gein,dgeindt,d2geindt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!  thet cannot depend on T
! Missing implem of derivatives wrt comp.dep of thet.
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
      gx%bmperr=4399; goto 1000
   endif
! NOTE g2val and derivatives not divided by RT !!
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
      expg2=one/expmg2
!   endif
!   dg2=log(one+expmg2)
   dg2=log(one+expmg2)
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
   phres%gval(1,1)=phres%gval(1,1)-dg2
   phres%gval(2,1)=phres%gval(2,1)-dgfdt
   phres%gval(4,1)=phres%gval(4,1)+d2g2dt2
!   write(*,19)'3H 2st:',phres%gval(1,1),phres%gval(2,1),phres%gval(4,1)
! save local values divided by RT?
900 continue
   addrec%propval=zero
   addrec%propval(1)=gein-dg2
   addrec%propval(2)=dgeindt-dgfdt
   addrec%propval(4)=d2geindt2-d2g2dt2
1000 continue
   return
 end subroutine calc_twostate_model_john

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_twostate_model1(moded,phres,addrec,lokph,mc,ceq)
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
   integer jj,noprop,ig2,ith,extreme,jth,kth
!   double precision del1,del2,del3,del4,gein,dgeindt,d2geindt2
   double precision gein,dgeindt,d2geindt2
   double precision xi,hump
   double precision, parameter :: humpfact=5.0D0
!   double precision g2ein,dg2eindt,d2g2eindt2,theta2,dcpl
   double precision kvot,expkvot,expmkvot,ln1mexpkvot,kvotexpkvotm1
   double precision g2val,dg2,expg2,expmg2,rt,tv,rg,dg2dt,dgfdt,d2g2dt2
! number of properties calculatied
   noprop=phres%listprop(1)-1
! locate the THET and G2 property record 
   ig2=0
   ith=0
   jth=0
   findix: do jj=2,noprop
      if(phres%listprop(jj).eq.addrec%need_property(1)) then
! current values of G2 is stored in phres%gval(1,ig2)
         ig2=jj;
      elseif(phres%listprop(jj).eq.addrec%need_property(2)) then
! current value of THET are stored in phres%gval(1,ith)
         ith=jj
      elseif(phres%listprop(jj).eq.22) then
! current value of DCP2 are stored in phres%gval(1,ith)
         jth=jj
      endif
   enddo findix
   if(ith.eq.0) then
      write(*,*)'Cannot find value for amorphous THET'
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
!   write(*,71)'3H Cp E3: ',extreme,ceq%tpval(1),gein,dgeindt,d2geindt2
70 format(a,F7.2,5(1pe12.4))
71 format(a,i3,1x,F7.2,5(1pe12.4))
!  thet cannot depend on T
!----------------------------------------------------------------
!-------------------------- two state part DIVIDE BY RT
! NOTE g2val and derivatives not divided by RT !!
   rt=ceq%rtn
   tv=ceq%tpval(1)
   rg=globaldata%rgas
   g2val=phres%gval(1,ig2); dg2dt=phres%gval(2,ig2)
   dg2=zero; d2g2dt2=zero
   expmg2=zero
!   write(*,*)'3H gval1: ',g2val
   if(g2val.eq.zero .and. dg2dt.eq.zero) then
!      write(*,*)'3H: G2 parameter zero, ignoring bump',g2val
      goto 900
   endif
   d2g2dt2=phres%gval(4,ig2)
! hump is an attempt to reduce the hump due to state change entropy
! This is testing a modification to prevent a hump
! If G^d is positive we are in LT range and hump=0
! If G^d is negative we are in HT rannge and 
! by scaling G^d to vary between 0 and 1 when G^d is negative (xi>0.5)
   xi=zero
   hump=one
   if(jth.gt.0) then
      hump=phres%gval(1,jth)
! there should be a chack if -200 < -g2val/rt or -g2val/rt > 200
      xi=exp(-g2val/rt)/(one+exp(-g2val/rt))
      write(*,19)'3H gval6: ',tv,g2val/rt,hump,xi
      if(-g2val/rt.gt.200) then
! Fraction liquid is very large 
         xi=one
         hump=one
!      elseif(-g2val/rt.lt.-200) then
! Fraction liquid is very small and can be ignored         
!         hump=zero
      else
! This is the intermediate range when hump*xi should approach unity
         if(-g2val/rt.lt.zero) then
            hump=0.5*humpfact*hump*xi
         else
            hump=0.5*humpfact*hump*(one-xi)+(2*xi-one)
         endif
      endif
   else
! This is classical Schottky model
      hump=one
   endif
!   write(*,19)'3H gval7: ',tv,g2val/rt,hump,xi,g2val*hump/rt
   g2val=hump*g2val
   dg2dt=hump*dg2dt
   d2g2dt2=hump*d2g2dt2
19 format(a,6(1pe11.3))
! if g2val is positive we are in the amorphous region
! if g2val is negative we are in the liquid region
! The if statements here ensure expmg2 is between 1e-60 and 1e+60
   if(-g2val/rt.gt.2.0D2) then
! exp(200) >> 1, thus d2g=ln(1+exp(g2val))=g2val
! and the derivatives are those above
      dg2=g2val
      dgfdt=dg2dt
      d2g2dt2=d2g2dt2
      goto 700
   elseif(-g2val/rt.lt.-2.0D2) then
! exp(-200)=0; ln(1)=0 and everything is zero
      dg2=zero
      dg2dt=zero
      d2g2dt2=zero
      goto 800
   else
! intermediate T range, we have to calculate
      expmg2=exp(-g2val/rt)
      expg2=one/expmg2
      dg2=log(one+expmg2)
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
   phres%gval(1,1)=phres%gval(1,1)-dg2
   phres%gval(2,1)=phres%gval(2,1)-dgfdt
   phres%gval(4,1)=phres%gval(4,1)+d2g2dt2
! values of T, \xi, g, s and cp   
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
! save local values divided by RT?
900 continue
   addrec%propval=zero
   addrec%propval(1)=gein-dg2
   addrec%propval(2)=dgeindt-dgfdt
   addrec%propval(4)=d2geindt2-d2g2dt2
1000 continue
   return
 end subroutine calc_twostate_model1

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
! locate the THET and G2 property record 
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
      write(*,*)'Cannot find value for amorphous THET'
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

!\begin{verbatim}
 subroutine create_crystalbreakdownmod(addrec)
! enters a record for the crystal breakdown
   implicit none
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim} %+
   integer typty
! reserve an addition record
   allocate(addrec)
! nullify link to next   
   nullify(addrec%nextadd)
! Set the type of addition and look for needed parameter properties
!   addrec%type=crystalbreakdownmod
   write(*,*)'3H crystal breakdown not an addition'
   gx%bmperr=4399; goto 1000
   allocate(addrec%need_property(1))
   call need_propertyid('CBT ',typty)
   if(gx%bmperr.ne.0) goto 1000
   addrec%need_property(1)=typty
! store zero.  Used to extract current value of this property
   addrec%propval=zero
!
!   write(kou,*)'3H Not implemented yet'; gx%bmperr=4078
!
1000 continue
   return
 end subroutine create_crystalbreakdownmod

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_crystalbreakdown_steep(moded,phres,addrec,lokph,mc,ceq)
! calculates the metastable extrapolation above crystal breakdown T (CBT)
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
   TYPE(gtp_phase_add), pointer :: addrec
   TYPE(gtp_phase_varres) :: phres
!\end{verbatim}
   integer ith,noprop,mc2,jj,jk,saveph
   double precision cbt,rg,tval,rt,psi
! BEWHARE: NOT IN PARALLEL
   double precision, allocatable :: cbtgval(:),cbtdgval(:,:),cbtd2gval(:)
   double precision, allocatable :: gsol(:),dgsol(:,:),d2gsol(:)
   save saveph,cbtgval,cbtdgval,cbtd2gval
   double precision x1,x2,x3
! extract the current value of the crystal breakdown T
   noprop=phres%listprop(1)-1
   do ith=2,noprop
      if(phres%listprop(ith).eq.addrec%need_property(1)) goto 100
   enddo
   write(*,*)'3H No value of CBT found for phase ',trim(phlista(lokph)%name)
   gx%bmperr=4336; goto 1000
100 continue
   cbt=phres%gval(1,ith)
! if current T lower than CBT just exit   
   if(ceq%tpval(1).le.cbt) goto 1000
! if T is higher we have to calculate everything for T=CBT ... HOW?   
!   write(*,110)'3H T > CBT: ',lokph,saveph,ceq%tpval(1),cbt
110 format(a,2i3,2F10.2)
   if(lokph.ne.saveph) then
      if(allocated(cbtgval)) then
         write(*,*)'3H deallocating as ',lokph,' not same as ',saveph
         deallocate(cbtgval)
         deallocate(cbtdgval)
         deallocate(cbtd2gval)
      endif
   endif
   if(.not.allocated(cbtgval)) then
      saveph=lokph
      write(*,*)'3H Allocating cbt arrays',saveph
      allocate(cbtgval(6))
      allocate(cbtdgval(3,mc))
      mc2=mc*(mc+1)
      allocate(cbtd2gval(mc2))
! assume current T is OK to save ...
      do jj=1,6
         cbtgval(jj)=phres%gval(jj,1)
      enddo
      write(*,101)'3H cbtgsol: ',cbtgval(1),cbtgval(2),cbtgval(4)
      do jk=1,mc
         do jj=1,3
            cbtdgval(jj,mc)=phres%dgval(jj,mc,1)
         enddo
      enddo
      do jk=1,mc2
         cbtd2gval(jj)=phres%d2gval(jk,1)
      enddo
! we do not have to weight together when we stored!!
      goto 1000
   endif
! for intermediate results, they are deallocated when leaving the routine
   allocate(gsol(6))
   allocate(dgsol(3,mc))
   allocate(d2gsol(mc2))
! Now we can weight together the values at CBT and current T. psi<1.0
   tval=ceq%tpval(1)
   psi=cbt/tval
   rg=globaldata%rgas*tval
   rt=rg*tval
   x1=3.0d0*(-tval/cbt*log(tval)+tval/cbt*log(cbt)-log(tval)+&
        2.0D0*tval/cbt+log(cbt)-2.0d0)
   gsol(1)=((tval**(-2)/6.0D0+cbt/(3.0d0*tval)+0.5D0*cbt**2)*cbtgval(4)/cbt**4+&
        (tval-cbt)*cbtgval(2)+cbtgval(1))*cbt/tval+x1*cbt/tval
!        3.0d0*log(psi)+1.5d0*tval*(one/psi-psi)
!   gsol(2)=(-(cbt/tval)**(-3)*cbt*cbtgval(4)/3.0d0+cbt*cbtgval(4)/3.0d0)/rt+&
!        cbtgval(2)/rt
   x1=-(cbt/tval)**3*cbt*cbtgval(4)/3.0d0
   x2=cbt*cbtgval(4)/3.0d0
   x3=cbtgval(2)
   gsol(2)=(x1+x2+x3)*cbt/tval&
        +3.0d0*(-log(tval)/cbt+log(cbt)/cbt-one/tval+one/cbt)*(cbt/tval)
! This Cp curve is OK, but the decrease of the Cp from LT is too steep
   x1=-3.0D0*(tval-cbt)/(cbt*tval**2)*(cbt/tval)
   gsol(4)=(cbt/tval)**4*cbtgval(4)*(cbt/tval)+x1
!        -3.0D0*(one/tval-one/cbt)/tval**3
   write(*,101)'3H extra: ',tval,cbt/tval,x1,x2,x3,gsol(4)
101 format(a,6(1pe10.2))
! this is just for a pure element !!
   phres%gval(1,1)=gsol(1)
   phres%gval(2,1)=gsol(2)
   phres%gval(4,1)=gsol(4)
1000 continue
   return
 end subroutine calc_crystalbreakdown_steep

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_crystalbreakdownmod(moded,phres,addrec,lokph,mc,ceq)
! calculates the metastable extrapolation above crystal breakdown T (CBT)
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
   TYPE(gtp_phase_add), pointer :: addrec
   TYPE(gtp_phase_varres) :: phres
!\end{verbatim}
   integer ith,noprop,mc2,jj,jk,saveph
   double precision cbt,rg,tval,rt,psi
! BEWHARE: NOT IN PARALLEL
   double precision, allocatable :: cbtgval(:),cbtdgval(:,:),cbtd2gval(:)
   double precision, allocatable :: gsol(:),dgsol(:,:),d2gsol(:)
   save saveph,cbtgval,cbtdgval,cbtd2gval
   double precision x1,x2,x3,dcpb
! extract the current value of the crystal breakdown T
   noprop=phres%listprop(1)-1
   do ith=2,noprop
      if(phres%listprop(ith).eq.addrec%need_property(1)) goto 100
   enddo
   write(*,*)'3H No value of CBT found for phase ',trim(phlista(lokph)%name)
   gx%bmperr=4336; goto 1000
100 continue
   cbt=phres%gval(1,ith)
! if current T lower than CBT just exit   
   if(ceq%tpval(1).le.cbt) goto 1000
! if T is higher we have to calculate everything for T=CBT ... HOW?   
!   write(*,110)'3H T > CBT: ',lokph,saveph,ceq%tpval(1),cbt
110 format(a,2i3,2F10.2)
   if(lokph.ne.saveph) then
      if(allocated(cbtgval)) then
         write(*,*)'3H deallocating as ',lokph,' not same as ',saveph
         deallocate(cbtgval)
         deallocate(cbtdgval)
         deallocate(cbtd2gval)
      endif
   endif
   if(.not.allocated(cbtgval)) then
      saveph=lokph
      write(*,*)'3H Allocating cbt arrays',saveph
      allocate(cbtgval(6))
      allocate(cbtdgval(3,mc))
      mc2=mc*(mc+1)
      allocate(cbtd2gval(mc2))
! assume current T is OK to save ...
      do jj=1,6
         cbtgval(jj)=phres%gval(jj,1)
      enddo
!      write(*,101)'3H cbtgsol: ',cbtgval(1),cbtgval(2),cbtgval(4)
      do jk=1,mc
         do jj=1,3
            cbtdgval(jj,mc)=phres%dgval(jj,mc,1)
         enddo
      enddo
      do jk=1,mc2
         cbtd2gval(jj)=phres%d2gval(jk,1)
      enddo
! we do not have to weight together when we stored!!
      goto 1000
   endif
! for intermediate results, they are deallocated when leaving the routine
   allocate(gsol(6))
   allocate(dgsol(3,mc))
   allocate(d2gsol(mc2))
! We should maybe approach zero as the Einstein  contribution is still there??
!   dcpb=3.0d0
   dcpb=zero
! Now we can weight together the values at CBT and current T. psi<1.0
   tval=ceq%tpval(1)
   psi=cbt/tval
   rg=globaldata%rgas*tval
   rt=rg*tval
   x1=dcpb*(-tval/cbt*log(tval)+tval/cbt*log(cbt)-log(tval)+&
        2.0D0*tval/cbt+log(cbt)-2.0d0)
!  gsol(1)=((tval**(-2)/6.0D0+cbt/(3.0d0*tval)+0.5D0*cbt**2)*cbtgval(4)/cbt**4+&
!        (tval-cbt)*cbtgval(2)+cbtgval(1))*cbt/tval+x1*cbt/tval
   gsol(1)=((tval**(-1)/2.0D0+cbt*tval/2.0d0+cbt**2)*cbtgval(4)/cbt**3+&
        (tval-cbt)*cbtgval(2)+cbtgval(1))*cbt/tval+x1*cbt/tval
!-----------------------------------------
!   x1=-(cbt/tval)**3*cbt*cbtgval(4)/3.0d0
!   x2=cbt*cbtgval(4)/2.0d0
   x1=-(cbt/tval)**2*cbt*cbtgval(4)/2.0d0
   x2=cbt*cbtgval(4)/2.0d0
   x3=cbtgval(2)
   gsol(2)=(x1+x2+x3)*cbt/tval&
        +dcpb*(-log(tval)/cbt+log(cbt)/cbt-one/tval+one/cbt)*(cbt/tval)
! This Cp curve is OK, but the decrease of the Cp from LT is too steep
   x1=-dcpb*(tval-cbt)/(cbt*tval**2)*(cbt/tval)
   gsol(4)=(cbt/tval)**3*cbtgval(4)*(cbt/tval)+x1
!
!   write(*,101)'3H extra: ',tval,cbt/tval,x1,x2,x3,gsol(4)
101 format(a,6(1pe10.2))
! this is just for a pure element !!
   phres%gval(1,1)=gsol(1)
   phres%gval(2,1)=gsol(2)
   phres%gval(4,1)=gsol(4)
1000 continue
   return
 end subroutine calc_crystalbreakdownmod

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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
   allocate(addrec%need_property(1))
   call need_propertyid('THET',typty)
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
   write(*,*)'3H No Debye temperature THET',lokph
   gx%bmperr=4336; goto 1000
100 continue
   write(*,*)'3H Deby low T heat capacity model not implemented'
   gx%bmperr=4078
1000 continue
   return
 end subroutine calc_debyecp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine create_diffusion(addrec,lokph,text)
   implicit none
   integer lokph
   character text*(*)
   type(gtp_phase_add), pointer :: addrec
!\end{verbatim}
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
! Set the type of addition and look for needed parameter properties
   addrec%type=diffcoefs
! Some information is needed
   last=1
100 continue
   call gparcd('Type of diffusion model: ',text,last,1,typ,'SIMPLE',nohelp)
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
      call gparcd(quest,text,last,1,typ,spname,nohelp)
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
      call gparrd('Value of ALPHA: ',text,last,alpha,0.3D0,nohelp)
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
               call gparrd(quest,text,last,alpha,1.0D0,nohelp)
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

!\begin{verbatim}
 subroutine diffusion_onoff(phasetup,bitval)
! switches the bit which calculates diffusion coefficients on/off
! if bitval 0 calculate is turned on, 1 turn off
   implicit none
   integer bitval
   type(gtp_phasetuple) :: phasetup
!\end{verbatim}
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
!\end{verbatim}
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
   if(.not.btest(lokadd%status,ADDHAVEPAR)) then
! skip additions with no parameters for this phase
!      write(*,*)chc,'3H No parameters for addition: ',&
!           trim(additioname(lokadd%type))
      goto 1000
!   else
!      write(*,*)'3H status word for this addition: ',lokadd%status
   endif
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
100      format(2x,'+ Magnetic model by Inden, anti-ferromagntic factor:',i3,/&
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
      write(unit,300)chc
300   format(a,'+ Inden magnetic model modified by Qing and Xiong'/&
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
400   format(a,'+ Einstein Cp model: 1.5R*THET(x) +',&
           ' 3RT*ln(exp(ln(THET(x))/T)-1)')
!---------------------------------------------
   case(elasticmodel1) ! Elastic model 1
      write(unit,500)
500   format(2x,'+ Elastic model 1, with P interpreted as a force in',&
           ' the X direction.')
!---------------------------------------------
   case(twostatemodel1) ! Liquid two state  model including Einstein
      write(unit,510)chc,chc
510   format(a,' + Liquid 2 state model: G(liq)-RT*ln(1+exp(-G2(x,T)/RT))'/&
           a,' + Einstein Cp model: 1.5R*THET(x) ',&
           '+ 3RT*ln(exp(ln(THET(x))/T)-1)')
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
   case(schottkyanomality) ! Schottky Anomality
      write(unit,550)chc
550   format(a,'+ Schottky anomality DSCH(x)*RT*ln(1+exp(-ln(TSCH(x))/T)) ')
   end select addition
1000 continue
   return
 end subroutine list_addition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

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

