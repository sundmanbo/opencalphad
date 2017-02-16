!
! gtp3H included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!
!>     14. Additions
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
!------------------------------------
! For each addition XX there is a subroutine create_XX
! called from the add_addrecord
! and a subroutine calc_XX 
! called from the addition_selector, called from calcg_internal
! There is a common list routine
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!

!\begin{verbatim}
 subroutine addition_selector(addrec,moded,phres,lokph,mc,ceq)
! called when finding an addtion record while calculating G for a phase
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
   addition: select case(addrec%type)
   case default
      write(kou,*)'No such addition type ',addrec%type,lokph
      gx%bmperr=4330
   case(indenmagnetic) ! Inden magnetic
      call calc_magnetic_inden(moded,phres,addrec,lokph,mc,ceq)
   case(debyecp) ! Debye Cp
      call calc_debyecp(moded,phres,addrec,lokph,mc,ceq)
      write(kou,*)' Debye Cp model not implemented yet'
      gx%bmperr=4331
   case(xiongmagnetic) ! Inden-Xiong
      call calc_xiongmagnetic(moded,phres,addrec,lokph,mc,ceq)
!     write(kou,*)'Inden magnetic model with sep TC and TN not implemented yet'
!      gx%bmperr=4332
   case(einsteincp) ! Einstein Cp
      call calc_einsteincp(moded,phres,addrec,lokph,mc,ceq)
      write(kou,*)' Einstein Cp model not implemented yet'
      gx%bmperr=4331
   case(elasticmodel1) ! Elastic model !
      call calc_elastica(moded,phres,addrec,lokph,mc,ceq)
      write(kou,*)' Elastic model not implemented yet'
      gx%bmperr=4399
   case(twostatemodel1) ! Two state model
      write(kou,*)' Two state model not implemented yet'
      gx%bmperr=4333
   case(volmod1) ! Simple volume model depending on V0, VA and VB
      call calc_volmod1(moded,phres,addrec,lokph,mc,ceq)
!      write(kou,*)' Two state model not implemented yet'
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
   case default
      write(kou,*)'No addtion type ',addtyp,lokph
   case(indenmagnetic) ! Inden magnetic
      if(extra(1:1).eq.'Y' .or. extra(1:1).eq.'y') then
! bcc model
         aff=-1
         call create_magrec_inden(newadd,aff)
      else
         aff=-3
         call create_magrec_inden(newadd,aff)
      endif
   case(xiongmagnetic) ! Inden-Xiong.  Assume bcc if BCC part of phase name
!      bcc=.false.
!      if(index('BCC',phlista(lokph)%name).gt.0) bcc=.true.
      if(extra(1:1).eq.'Y' .or. extra(1:1).eq.'y') then
         bcc=.TRUE.
      else
         bcc=.FALSE.
      endif
      call create_xiongmagnetic(newadd,bcc)
   case(debyecp) ! Debye Cp
      call create_debyecp(newadd)
   case(einsteincp) ! Einstein Cp
      call create_einsteincp(newadd)
   case(elasticmodel1) ! Elastic model 1
      call create_elastic_model_a(newadd)
   case(twostatemodel1) ! Two state model 1
      call create_twostate_model1(newadd)
   case(volmod1) ! Volume model 1
      call create_volmod1(newadd)
   end select addition
   if(gx%bmperr.ne.0) goto 1000
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
! subroutine add_magrec_inden(lokph,addtyp,aff)
! adds a magnetic record to lokph
! lokph is phase location
! addtyp should be 1 of Inden model
! aff is antiferromagnic factor, -1 for bcc and -3 for fcc and hcp
!   implicit none
!   integer lokph,addtyp,aff
!\end{verbatim} %+
!   integer mc
!   type(gtp_phase_add), pointer :: newadd,addrec,addrec1
!   mc=phlista(lokph)%tnooffr
! check if there are other additions
!   addrec1=>phlista(lokph)%additions
!   addrec=>addrec1
!200 do while(associated(addrec))
!      if(addrec%type.eq.indenmagnetic) then
!         write(*,*)'3H Inden magnetic model already entered'
!         goto 1000
!      else
!         addrec1=addrec
!         addrec=>addrec%nextadd
!      endif
!   enddo
! add the addition record last
!   call create_magrec_inden(newadd,aff)
!   if(gx%bmperr.ne.0) goto 1000
!   addrec1%nextadd=>newadd
!1000 continue
!   return
! end subroutine add_magrec_inden

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
   double precision logb1,invb1,iafftc,iaffbm,rgasm,rt,tao,gmagn
   double precision dtaodt,dtaodp,beta,d2taodp2,d2taodtdp,tc,tv
   double precision tao2(2),ftao(6),dtao(3,mc),d2tao(mc*(mc+1)/2)
   double precision addgval(6),daddgval(3,mc),d2addgval(mc*(mc+1)/2)
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
!               rt*ftao(2)*d2tao(ixsym(j,k))*logb1,&
!               rt*ftao(2)*dtao(1,j)*invb1*phres%dgval(1,k,ibm),&
!               rt*ftao(2)*dtao(1,k)*invb1*phres%dgval(1,j,ibm),&
!              -rt*ftao(1)*invb1**2*phres%dgval(1,j,ibm)*phres%dgval(1,k,ibm),&
!               rt*ftao(1)*invb1*phres%d2gval(ixsym(j,k),ibm)
!57 format('3H mag2y: ',6(1PE12.4))
      enddo
   enddo
! now add all to the total G and its derivatives
! something wrong here, j should go from 1 to 9 in my fenix case ...
   do j=1,mc
!      write(*,99)'3H magadd 1: ',1,j,phres%dgval(1,j,1),daddgval(1,j)/rt
      do k=1,3
! first derivatives
         phres%dgval(k,j,1)=phres%dgval(k,j,1)+daddgval(k,j)/rt
      enddo
99    format(a,2i3,2(1pe16.8))
      do k=j,mc
! second derivatives
!         write(*,99)'3H magadd 2: ',k,j,rt*phres%d2gval(ixsym(j,k),1),&
!              d2addgval(ixsym(j,k))
         phres%d2gval(ixsym(j,k),1)=phres%d2gval(ixsym(j,k),1)+&
              d2addgval(ixsym(j,k))/rt
      enddo
   enddo
!   write(*,*)'3H cm 7: ',phres%gval(1,1),addgval(1)/rt
! note phres%gval(1..3,1) already calculated above
   do j=4,6
      lokadd%propval(j)=addgval(j)
      phres%gval(j,1)=phres%gval(j,1)+addgval(j)/rt
   enddo
1000 continue
   return
 end subroutine calc_magnetic_inden

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine create_xiongmagnetic(addrec,bcc)
! adds a Xiong type magnetic record, we must separate fcc and bcc by extra
! copied from Inden magnetic model
! The difference is that it uses TCA for Curie temperature and TNA for Neel
! and individual Bohr magneton numbers
   implicit none
   logical bcc
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
! This model use an effective Bohr magneton number b*=prod(b_i+1)**x_i -1
!   call need_propertyid('IBM ',typty)
   call need_propertyid('BMAG ',typty)
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
! calculates Indens-Xiong magnetic contribution
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
   double precision logb1,invb1,iafftc,iaffbm,rgasm,rt,tao,gmagn
   double precision dtaodt,dtaodp,beta,d2taodp2,d2taodtdp,tc,tv
   double precision tao2(2),ftao(6),dtao(3,mc),d2tao(mc*(mc+1)/2)
   double precision addgval(6),daddgval(3,mc),d2addgval(mc*(mc+1)/2)
   double precision tn,tcsave,tnsave
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
! Inden magnetic need properties in need_property(1..3)
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
            ' not calculated as '/10x&
            'values for CTA, NTA or BMAG, magnetic G are zero')
      goto 1000
   else
      tc=-one
      tn=-one
      if(itc.gt.0) tc=phres%gval(1,itc)
      if(itn.gt.0) tn=phres%gval(1,itn)
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
   do j=1,mc
      do k=1,3
!          write(*,99)'3H magadd 1: ',k,j,rt*phres%dgval(k,j,1),daddgval(k,j)
         phres%dgval(k,j,1)=phres%dgval(k,j,1)+daddgval(k,j)/rt
      enddo
!99 format(a,2i3,2(1pe16.8))
      do k=j,mc
!          write(*,99)'3H magadd 2: ',k,j,rt*phres%d2gval(ixsym(j,k),1),&
!               d2addgval(ixsym(j,k))
         phres%d2gval(ixsym(j,k),1)=phres%d2gval(ixsym(j,k),1)+&
              d2addgval(ixsym(j,k))/rt
      enddo
   enddo
!    write(*,*)'3H cm 7: ',rt*phres%gval(1,1),addgval(1)
! note phres%gval(1..3,1) already calculated above
   do j=4,6
      lokadd%propval(j)=addgval(j)
      phres%gval(j,1)=phres%gval(j,1)+addgval(j)/rt
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
! calculate the simple volume model
!
! G = P*V0*exp(VA(T))
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
   integer jl,iv0,iva,ivb,noprop
   double precision v0,va,vb,vol,deltap,pvol
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
! if iva is zero there are no volume data
   if(iva.eq.0) goto 1000
! reference pressure is 1 bar
   deltap=ceq%tpval(2)-1.0D5
   if(ivb.ne.zero) then
! if ivb not zero there are bulk modulus data
      write(*,*)'3H Volume model with bulk modulus not implemented'
   elseif(v0.ne.zero) then
! NOTE all values should be divided by RT
      vol=v0*exp(va)/ceq%rtn
      pvol=deltap*vol/ceq%rtn
      if(iva.gt.0) then
!         write(*,17)'3H volume: ',lokph,iv0,iva,v0,va,vol,pvol,&
!              phres%gval(2,iva),phres%gval(4,iva)
17       format(a,3i3,6(1pe11.3))
      else
!         write(*,17)'3H volume: ',lokph,iv0,iva,vol,pvol
      endif
! contribtions to G and derivatives, G, G.T, G.P=V, G.T.T, G.T.P, G.P.P
! NOTE other parameters may depend on P!
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
!      phres%gval(6,1)=phres%gval(6,1)
! for the moment ignore composition dependence ...
! store some property values
      lokadd%propval(1)=pvol
      lokadd%propval(2)=pvol*phres%gval(2,iva)
      lokadd%propval(3)=vol
      lokadd%propval(4)=pvol*(phres%gval(2,iva)**2+phres%gval(4,iva))
      lokadd%propval(5)=vol*phres%gval(2,iva)
      lokadd%propval(6)=zero
   endif
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
! G/RT = 3*ln( 1 - exp( THET/T ) ) 
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
! G = 3*R*T*ln( 1 - exp( THET/T ) ) 
! This is easier to handle inside the calc routine without TPFUN
!
   implicit none
   integer moded,lokph,mc
   type(gtp_phase_varres), pointer :: phres
   type(gtp_phase_add), pointer :: addrec
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ith,noprop
   double precision del1,del2,del3,gein,dgeindt,d2geindt2
!
   noprop=phres%listprop(1)-1
   findix: do ith=2,noprop
      if(phres%listprop(ith).eq.addrec%need_property(1)) goto 100
   enddo findix
   write(*,*)'3H No theta value. ',lokph
   gx%bmperr=4336; goto 1000
100 continue
! thet is in gval(ith,1), derivatives in dgval(*,ith,*) and d2gval(ith,*)
! G/RT = 3*ln( 1 - exp( THET/T ) ) 
! NOTE DIRIVATES CALCULATED FOR G/RT
   del1=phres%gval(ith,1)/ceq%tpval(1)
   del2=exp(del1)
   del3=1.0d0-del2
   gein=3.0D0*log(del3)
   dgeindt=3.0D0*(del1/ceq%tpval(1))*(del2/del3)
!   d2geindt2=3.0D0*(del1**2/ceq%tpval(1))*(del2/del3**2)
   d2geindt2=dgeindt*del1/del3
! Missing implem of derivatives wrt fractions of thet.  thet cannot depend on T
1000 continue
   return
 end subroutine calc_einsteincp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine create_twostate_model1(newadd)
! not implemented
   implicit none
   type(gtp_phase_add), pointer :: newadd
!\end{verbatim}
   write(kou,*)'Not implemented yet'; gx%bmperr=4078
1000 continue
   return
 end subroutine create_twostate_model1

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
   write(*,*)'3H Not implemented yet'
   gx%bmperr=4078
1000 continue
   return
 end subroutine calc_debyecp

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine list_addition(unit,CHTD,phname,ftyp,lokadd)
! list description of an addition for a phase on unit
   implicit none
   integer unit,ftyp
   character CHTD*1,phname*(*)
   TYPE(gtp_phase_add), pointer :: lokadd
!\end{verbatim} %+
   integer ip
   TYPE(tpfun_expression), pointer :: exprot
   character line*256,tps(2)*3
   double precision ff
!
   addition: select case(lokadd%type)
   case default
      write(unit,*)'Unknown addtion type: ',phname,lokadd%type
   case(indenmagnetic) ! Inden magnetic model
      if(ftyp.eq.2) then
! TDB file: a do not think I have saved the enthalpy factor, bcc (-1) it is 0.4
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
      write(unit,200)
200   format(2x,'+ Debye Cp model, not implemented yet')
!---------------------------------------------
   case(xiongmagnetic) ! Inden-Xiong
      write(unit,300)
300   format(2x,'+ Inden magnetic model modified by Xiong'/&
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
      write(unit,400)
400   format(2x,'+ Einstein Cp model:'/4x,'G = 3*R*T*LN(1-THET/T)')
!---------------------------------------------
   case(elasticmodel1) ! Elastic model 1
      write(unit,500)
500   format(2x,'+ Elastic model 1, with P interpreted as a force in',&
           ' the X direction.')
!---------------------------------------------
   case(twostatemodel1) ! Two state  model 1
      write(unit,*)'Two state model 1, not implemented yet'
!---------------------------------------------
   case(volmod1) ! Volume model 1
      if(unit.eq.kou) then
         write(unit,*)' + Volume model V=V0(x)*exp(VA(x,T))'
      else
         write(unit,*)'$  + Volume model V=V0(x)*exp(VA(x,T))'
      endif
   end select addition
1000 continue
   return
 end subroutine list_addition

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine list_addition_values(unit,phres)
! lists calculated values for this addition
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
      write(lut,77)addrec%type,(addrec%propval(j1),j1=1,4)
77    format('Addition type ',i2,': ',4(1pe12.4))
      addrec=>addrec%nextadd
   enddo
1000 continue
   return
 end subroutine list_addition_values

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

