!
! gtp3X included in gtp3.F9029
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     15. Calculate things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calcg(iph,ics,moded,lokres,ceq)
! calculates G for phase iph and composition set ics in equilibrium ceq
! checks first that phase and composition set exists
! Data taken and stored in equilibrium record ceq
! lokres is set to the phase_varres record with all fractions and results
! moded is 0, 1 or 2 depending on calculating no, first or 2nd derivarives
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer iph,ics,moded,lokres
!\end{verbatim}
   integer jcs,lokcs,lokph
   if(gx%bmperr.ne.0) then
      write(*,*)'3X Error code set when calling calcg: ',gx%bmperr
      goto 1000
   endif
   if(iph.le.0 .or. iph.gt.noofph) then
! the selected_element_reference phase with iph=0 is calculated separtely
      gx%bmperr=4050; goto 1000
   endif
   lokph=phases(iph)
   if(lokph.le.0 .or.lokph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   endif
!    write(*,*)'3X calcg 1: ',phlista(lokph)%name
! find fractions for this composition set
   if(ics.le.1) then
      jcs=1
   elseif(ics.le.phlista(lokph)%noofcs) then
      jcs=ics
   else
! no such composition set
!      write(*,*)'3X calcg 1 error 4072'
      gx%bmperr=4072; goto 1000
   endif
!   if(phlista(1)%noofcs.gt.1) then
! strange error that liquid (phase 1) has 3 composition set
!      write(*,*)'3X csbug: ',lokph,jcs,phlista(1)%noofcs
!      stop 'csbug'
!   endif
! Find fraction record this composition set
   lokcs=phlista(lokph)%linktocs(ics)
!-----
!   mcs=1
!   lokcs=phlista(lokph)%cslink
!   do while(mcs.lt.jcs)
!      mcs=mcs+1
! firsteq is the first equilibrium and a global variable in this module
!      lokcs=firsteq%phase_varres(lokcs)%next
!      if(lokcs.le.0) then
!         write(*,*)'3X calcg 2 error 4072'
!         gx%bmperr=4072; goto 1000
!      endif
!   enddo
   lokres=lokcs
!   write(*,*)'3X calcg 7: ',lokres,ceq%eqname(1:10)
! call using the local structure phase_varres
! results can be obtained through lokres
!   write(*,17)'3X calcg: ',lokph,lokres,ceq%phase_varres(lokres)%yfr(1)
!17 format(a,2i4,1pe15.6)
   call calcg_internal(lokph,moded,ceq%phase_varres(lokres),ceq)
1000 continue
   return
 end subroutine calcg

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calcg_internal(lokph,moded,cps,ceq)
! Central calculating routine calculating G and everyting else for a phase
! ceq is the equilibrium record, cps is the phase_varres record for lokph
! moded is type of calculation, 0=only G, 1 G and first derivatives
!    2=G and all second derivatives
! Can also handle the ionic liquid model now ....
   implicit none
   integer lokph,moded
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_phase_varres), target :: cps
!\end{verbatim}
! fractype defines fraction type (1=constituent fractions)
! empermut and ipermut permutation of fractions for phases with option F and B
! permrecord, maxprec and sameint to handle permutation in the interaction tree
   integer, parameter :: permstacklimit=150
   integer fractype,epermut,ipermut,typty,pmq,maxprec
   integer sameint(5)
   integer, dimension(permstacklimit) :: lastpmq,maxpmq
!   character bug*60
!   dimension sites(maxsubl),pushpop(maxpp)
   double precision, dimension(:), allocatable :: dpyq(:),d2pyq(:),d2vals(:)
   double precision, dimension(:,:), allocatable :: dvals(:,:)
   double precision vals(6)
! this array has the sum of constituents up to and including current sublattice
   integer incffr(0:maxsubl)
! in local gz: gz%intlevel level of interaction, gz%intcon and gz%intlat are
! used also in cgint when calculating interactions.
   TYPE(gtp_parcalc) :: gz
! disordered fraction set
   TYPE(gtp_fraction_set) :: fracset,dislink
   TYPE(gtp_phase_varres), pointer :: phres,phpart,phmain
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_endmember), pointer :: endmemrec
   TYPE(gtp_interaction), pointer :: intrec
   TYPE(gtp_pystack), pointer :: pystack
   TYPE(gtp_phase_add), pointer :: addrec
! for an ordered phase like FCC with a disordered contribution one must
! calculate the ordered part twice, one with original fractions and once
! with these replaced by the disordered fractions. and subdrahera.  This means
! one must have space to save fractions and results
!   double precision, dimension(:), allocatable :: savey
   double precision, dimension(maxconst) :: savey
   double precision, dimension(:,:), allocatable :: saveg
   double precision, dimension(:,:,:), allocatable :: savedg
   double precision, dimension(:,:), allocatable :: saved2g
   double precision, dimension(:,:), allocatable :: tmpd2g
! added when implicit none
   double precision rtg,pyq,ymult,add1,sum,yionva,fsites,xxx,sublf
   integer nofc2,nprop,nsl,msl,lokdiseq,ll,id,id1,id2,lm,qz,floryhuggins
   integer lokfun,itp,nz,intlat,ic,jd,jk,ic1,jpr,ipy,i1,j1
   integer i2,j2,ider,is,kk,ioff,norfc,iw,iw1,iw2,lprop,jonva,icat
   integer nsit1,nsit2
! storage for calculated Flopry Huggins volume parameters
   integer, dimension(:), allocatable :: fhlista
   double precision, dimension(:,:), allocatable :: fhv
   double precision, dimension(:,:,:), allocatable :: dfhv
   double precision, dimension(:,:), allocatable :: d2fhv
! to handle parameters with wildcard constituent and other things
   logical wildc,nevertwice,first,chkperm,ionicliq,iliqsave,iliqva,iliqneut
! debugging for partitioning and ordering
!   integer clist(4)
! calculate RT to normalize all Gibbs energies, ceq is current equilibrium
   rtg=globaldata%rgas*ceq%tpval(1)
   ceq%rtn=rtg
!-----------------------
! this is used for the Flory-Huggins model
   floryhuggins=0
   chkperm=.false.
   if(btest(phlista(lokph)%status1,PHFORD) .or. &
        btest(phlista(lokph)%status1,PHBORD)) then
      chkperm=.true.
! This is needed only once unless parameters are changed.  It numbers the
! interaction records sequentially for the permutations
      call palmtree(lokph)
      if(gx%bmperr.ne.0) goto 1000
   endif
!   if(ocv()) write(*,*)'3X in gcalc_internal: ',lokph
!-----------------------------------------------------------------
50  continue
! local work arrays for products of Y and calculated parameters are allocated
   gz%nofc=phlista(lokph)%tnooffr
   nofc2=gz%nofc*(gz%nofc+1)/2
!   write(*,17)'3X calcg, ',lokph,gz%nofc,nofc2,size(cps%d2gval),cps%nprop,&
!        cps%yfr(1)
!17 format(a,5i4,1pe15.6)
! for disordered fraction sets gz%nofc must be from disordered fraction record
! maybe these should not be allocated for moded=0 and 1
!   if(ocv()) write(*,*)'3X First allocate: ',gz%nofc,nofc2
   allocate(dpyq(gz%nofc))
   allocate(d2pyq(nofc2))
! these return values from excess parameters that may depend on constitution
   allocate(dvals(3,gz%nofc))
   allocate(d2vals(nofc2))
   nullify(pystack)
! do they have to be zeroed? YES!
   dpyq=zero
   d2pyq=zero
! dimension for number of parameter properties
   nprop=cps%nprop
! phres will point either to ordered or disordered results
! phmain will always point to record for ordered phase_varres
   phmain=>cps
   phres=>cps
! zero result arrays for all properties, maybe one should do it separately for
! each property as it is found but it may be faster to do it like this anyway
   phres%gval=zero
   if(moded.gt.0) then
      phres%dgval=zero
      if(moded.gt.1) then
         phres%d2gval=zero
      endif
   endif
! copy current values of T, P and RT from gtp_phase_varres
   gz%tpv(1)=ceq%tpval(1)
   gz%tpv(2)=ceq%tpval(2)
!   write(*,*)'3X calcg_i: ',gz%tpv
   gz%rgast=ceq%tpval(1)*globaldata%rgas
!   gz%rgast=ceq%tpval(1)*ceq%rgas
! this is used to check the number of times an ordered phase is calculated
   first=.true.
!-------------------------------------------------------------------
! calculate configurational entropy.
   nsl=phlista(lokph)%noofsubl
   ionicliq=.FALSE.
   iliqsave=.FALSE.
   if(btest(phlista(lokph)%status1,PHIONLIQ)) then
      call config_entropy_i2sl(moded,nsl,phlista(lokph)%nooffr,phres,&
           phlista(lokph)%i2slx,gz%tpv(1))
      ionicliq=.TRUE.
!      iliqsave=.FALSE.
      iliqva=.FALSE.
      jonva=0
!      write(*,*)'3X Config G 1: ',phres%gval(1,1)*rtg
!      if(phlista(lokph)%i2slx(1).gt.phlista(lokph)%tnooffr .and. &
!           phlista(lokph)%i2slx(2).gt.phlista(lokph)%tnooffr) then
!         onlyanions=.TRUE.
!      else
!         onlyanions=.FALSE.
!      endif
   elseif(btest(phlista(lokph)%status1,PHFHV)) then
! Flory-Huggins model require special treatment to calculate the molar
! volumes of the constituents.  The entropy is calculated in the end and
! a second loop through all parameters done by jumping to label 100
! check that just one sublattice and sites equal to one
      if(nsl.ne.1 .or. phres%sites(1).ne.one) then
         write(*,*)'3X Flory-Huggins model must have one lattice and site'
         gx%bmperr=4337; goto 1000
      endif
      floryhuggins=-1
   else
! NOTE: for phases with disordered fraction set this is calculated
! ONLY for the ordered original constituent fraction set
      call config_entropy(moded,nsl,phlista(lokph)%nooffr,phres,gz%tpv(1))
   endif
   if(gx%bmperr.ne.0) goto 1000
!-------------------------------------------------------------------
! start BIG LOOP for all fraction variables and parameters
! there may be several different properties in addition to G like TC, MQ& etc
! each of these are stored in separate gval(*,ipy) where ipy is an integer
! set for each property. lprop is incremented by one for each new property
! found (each phase may have different) and in listprop the original type
! of property is stored.  listprop will always be associated with phmain
100 continue
! yionva is used as indicator below if there are Va or just neutrals ...
   yionva=zero
   nevertwice=.true.
   lprop=2
   phmain%listprop(1)=1
   fractype=0
!   write(*,*)'3X calcg 99: ',lokph,cps%phtupx,cps%disfra%varreslink
!--------------------------------------------------------------------
! VERY STRANGE ERROR
! wrong results calculating with disordered fraction set disappeared
! when adding this write statements (and the one after the else statement
! and the one with the text nevertwice: further below)
!   write(*,101)'3X calcg 100 ',nsl,phlista(lokph)%nooffs,&
!        btest(phmain%status2,CSORDER),phres%gval(1,1),cps%gval(1,1)
101 format(a,2i4,1x,l,4(1pe14.4))
!--------------------------------------------------------------------
! loop for different types of fractions: site fractions, mole fractions ...
   fractyp: do while(fractype.lt.phlista(lokph)%nooffs)
105 continue
!     write(*,7)'3X label 105: ',fractype,btest(phlista(lokph)%status1,PHSUBO),&
!           btest(phmain%status2,CSORDER),btest(phlista(lokph)%status1,PHMFS),&
!           fracset%totdis,phres%gval(1,1)
7     format(a,i2,3(1x,l),i3,3(1pe12.4))
      fractype=fractype+1
! return here for calculating with disordered fractions for same fraction type
110 continue
! gz%nofc is number of fraction variables, msl is number of sublattices
! for this set of fractions!!! Ordering in FCC may have 5 sublattices with
! 4 participating in ordering and one interstitial.  The second fraction
! set may have 2 sublattices, 1 for the 4 ordering and one interstitial
      fracset=phmain%disfra
      ftype: if(fractype.eq.1) then
!---------------------------------------------- ordered (or only) fraction set
         if(btest(phlista(lokph)%status1,PHMFS)) then
! there is a disordered fractions set, we need additional fracset
            if(fracset%totdis.ne.0) then
               if(btest(phlista(lokph)%status1,PHSUBO)) then
! if phsubo set skip subtracting the ordered part as disordered, just add
                  goto 106
               endif
! the phase can totally disorder, disordered skip ordered part
! the CSORDER bit set by calc_disfrac called from set_constitution
               if(btest(phmain%status2,CSORDER)) then
! the phase is ordered, we have to calculate this part twice
!                  write(*,*)'3X Setting nevetwice false'
                  nevertwice=.false.
! independent if ordered or disordered always calculate first fraction set
               else
! the phase is disordered, skip ordered part and just calculate disordered
!                  write(*,*)'3X Skipping ordered part'
                  goto 105
               endif
            endif
         endif
106      continue
         gz%nofc=phlista(lokph)%tnooffr
         msl=nsl
         incffr(0)=0
         do qz=1,nsl
            incffr(qz)=incffr(qz-1)+phlista(lokph)%nooffr(qz)
         enddo
!         write(*,*)'3X after 106: ',fractype,nevertwice
! the results will be stored in the results arrays indicated by phres
! it was set above for the ordered fraction set. 
      else
!-------------------------------------------------
! disorderd/other fraction sets, take data from  gtp_fraction_set
!         write(*,*)'3X Fraction type: ',fractype,cps%disfra%varreslink
         msl=fracset%ndd
         gz%nofc=fracset%tnoofxfr
         incffr(0)=0
         do qz=1,msl
            incffr(qz)=incffr(qz-1)+fracset%nooffr(qz)
         enddo
! we have to deallocate and allocate local arrays, not if moded=0 or 1??
         deallocate(dpyq)
         deallocate(d2pyq)
         allocate(dpyq(gz%nofc))
         allocate(d2pyq(nofc2))
!         if(ocv()) write(*,*)'3X Allocated dpyq 2'
         dpyq=zero
         deallocate(dvals)
         deallocate(d2vals)
         allocate(dvals(3,gz%nofc))
         allocate(d2vals(nofc2))
         if(ocv()) write(*,*)'3X Allocated vals 2'
! the results will be stored in result arrays indicated by phres
! for the disordered fraction set phres must be set here and the arrays zeroed
         dislink=cps%disfra
!         write(*,*)'3X Calc internal disordred part 1A',dislink%fsites
         lokdiseq=dislink%varreslink
!         write(*,*)'3X Calc internal disordred part 1B',lokdiseq
         phres=>ceq%phase_varres(lokdiseq)
         phres%gval=zero
!         write(*,*)'3X Calc internal disordred part 1c'
         if(moded.gt.0) then
            phres%dgval=zero
            if(moded.gt.1) then
               phres%d2gval=zero
            endif
         endif
!         write(*,*)'3X Calc internal disordred part 2'
      endif ftype
!==========================================================
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! code below is an attempt to parallelize the calculation of each
! endmember tree for a single phase ...
! It is commented away as there has been som many changes
!
! there can be ordered and disordered fraction sets selected by fractype
! One endmember at a time but to speed up when having several
! CPU we give one endmamber plus its interaction tree to each tread.  
! To handle this all endmember records should be in an array
!      if(fractype.eq.1) then
! TYPE gtp_phase must be extended with these lists
!         endmemrec=>phlista(lokph)%ordered
!         oendmems: do i=1,phlista(lokph)%noemr
!            call calc_endmemtree(lokph,moded,msl,&
!                 phlista(lokph)%oendmemarr(i)%p1,phres,phmain,ceq)
!         enddo oendmems
!      else
!         endmemrec=>phlista(lokph)%disordered
!         dendmems: do i=1,phlista(lokph)%ndemr
!            call calc_endmemtree(lokph,moded,msl,&
!                 phlista(lokph)%dendmemarr(i)%p1,phres,phmain,ceq)
!         enddo dendmems
!      endif
!
! calculated for this fraction type, initiation for next in the beginning of
! loop but we may have to calculate once again with same fraction type but
! with the fractions as disordered fractions
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!==========================================================
! there can be ordered and disordered fraction sets selected by fractype
      if(fractype.eq.1) then
         endmemrec=>phlista(lokph)%ordered
      else
         endmemrec=>phlista(lokph)%disordered
      endif
!
! here we take one endmember at a time but to speed up when having several
! CPU we give one endmamber plus its interaction tree to each tread.  
! To handle this all endmember records should be in an array.  An attempt to
! implement this was made in calcg_internal2 but not updated for permutations
!
! empermut, lastpmq and maxpmq controls permutations (option F and B)
! maxpmq is set to zero for each new endmember but keep its content
! during calculation of all permutations of the same endmember and interactions
! big loop for all permutation of fractions (ordering option F and B)
! including all interaction parameters linked from this endmember
!
!      write(*,*)'3X Config G 2: ',phres%gval(1,1)*rtg
      endmemloop: do while(associated(endmemrec))
!
! The array maxpmq is used for interaction permutations.  It must be
! initialized to zero at the first endmember permutation.  It is set to
! limits for the interacton permutations for all interaction records.
         maxpmq=0
         maxprec=0
         epermut=0
         sameint=0
!         write(*,*)'3X: start endmember list'
         empermut: do while(epermut.lt.endmemrec%noofpermut)
            epermut=epermut+1
! calculate py, calculate parameter, calculate contribution to G etc
! py is product of all fractions, dpy are first derivatives and d2py second
            pyq=one
            if(moded.gt.0) then
! moded=0, only G, =1 only G and dG/dy, moded=2 all Gm dG/dy and d2G/dy2
               dpyq=zero
               if(moded.gt.1) then
                  d2pyq=zero
               endif
            endif
!-----------------------------------------------------
            pyqloop: do ll=1,msl
               id=endmemrec%fraclinks(ll,epermut)
! remove next line when all fixed
!               if(ll.lt.5) clist(ll)=id
! id negative means wildcard, independent of the fraction in this sublattice
               if(id.lt.0) then
                  gz%yfrem(ll)=one
               else
                  gz%yfrem(ll)=phres%yfr(id)
                  if(gz%yfrem(ll).lt.bmpymin) gz%yfrem(ll)=bmpymin
                  if(gz%yfrem(ll).gt.one) gz%yfrem(ll)=one
               endif
! gz%endcon is used for interaction parameters below
               gz%endcon(ll)=id
               pyq=pyq*gz%yfrem(ll)
               if(ionicliq .and. ll.eq.2) then
! For ionic liquid we must check when Va or neutral in second sublattice
! i2slx(1) is index of vacancy, i2slx(2) is first neutral
                  if(id.eq.phlista(lokph)%i2slx(1) .and. yionva.eq.zero) then
                     iliqva=.TRUE.
                     yionva=gz%yfrem(ll)
                     jonva=phlista(lokph)%i2slx(1)
! We found Va.  Save all calculated values as the follwing terms should all
! be multiplied with Q (done after finishing calculation)
!                     nprop=phmain%nprop
! we have already extracted nprop above .... 
                     allocate(saveg(6,nprop))
                     saveg=phres%gval
!                     if(ocv()) write(*,*)'3X saveg allocated 1A:',size(saveg),&
!                          gz%nofc,nofc2,nprop,moded
                     if(moded.gt.0) then
! only allocate if needed, some "out of memory" problems here calculating grid
! with just ionic liquid phase
                        allocate(savedg(3,gz%nofc,nprop))
                        allocate(saved2g(nofc2,nprop))
                        savedg=phres%dgval
                        saved2g=phres%d2gval
                     endif
!                     if(ocv()) write(*,*)'3X saveg allocated 1B: '
!                     write(*,*)'3X Config G 3A: ',phres%gval(1,1)*rtg
                     phres%gval=zero
                     phres%dgval=zero
                     phres%d2gval=zero
!                     write(*,*)'3X Config G 3B: ',phres%gval(1,1)*rtg
                     iliqsave=.TRUE.
!                     write(*,117)'3X Saved ionliq G at Va id: ',&
!                          id,yionva,saveg(1,1)
117                  format(a,i3,6(1pe12.4))
                  elseif(id.eq.phlista(lokph)%i2slx(2)) then
! we have a neutral in second sublattice
                     iliqva=.FALSE.
                     yionva=-one
                     jonva=0
                     if(.not.iliqsave) then
! We may have model without Va, for exampel (Ca+2)p(O-2,SiO4-4,SiO2)q, if so
! we must save all calculated values as the rest should be multiplied with Q
!                        nprop=phmain%nprop
! we already know nprop from above
                        allocate(saveg(6,nprop))
                        allocate(savedg(3,gz%nofc,nprop))
                        allocate(saved2g(nofc2,nprop))
!                        if(ocv()) write(*,*)'3X saveg allocated 2:',size(saveg)
                        saveg=phres%gval
                        savedg=phres%dgval
                        saved2g=phres%d2gval
                        phres%gval=zero
                        phres%dgval=zero
                        phres%d2gval=zero
                        iliqsave=.TRUE.
!                        write(*,117)'3X Saved ionliq G at neutral id: ',&
!                           id,yionva,saveg(1,1)
                     endif
                  endif
               endif
            enddo pyqloop
            if(moded.eq.0) goto 150
!---------------------------------------------------- first derivatives
            dpyqloop: do ll=1,msl
! here pyq is known, same loop as above to calculate dpyq(i)=pyq/y_i
               id=endmemrec%fraclinks(ll,epermut)
               if(id.gt.0) then
! pyq was multiplied with gz%yfrem(11) above, now divide with it
                  dpyq(id)=pyq/gz%yfrem(ll)
!                  write(*,*)'3X dpq/dy: ',ll,id,dpyq(id)
               elseif(.not.ionicliq) then
! wildcard in the sublattice and NOT ionic liquid
                  do iw=incffr(ll-1)+1,incffr(ll)
                     dpyq(iw)=pyq
                  enddo
               elseif(ll.ne.1) then
! wildcard in second subl of ionic liquid, same as for CEF
                  do iw=incffr(ll-1)+1,incffr(ll)
                     dpyq(iw)=pyq
                  enddo
!               else
! wildcard in first subl of ionic liquid then just ignore first derivatives
! with respect to constituents in first sublattice
!                  continue
               endif
            enddo dpyqloop
            if(moded.le.1) goto 150
!---------------------------------------------------- second derivatives
            d2pyqloop1: do ll=1,msl
               id1=endmemrec%fraclinks(ll,epermut)
               d2pyloop2: do lm=ll+1,msl
                  id2=endmemrec%fraclinks(lm,epermut)
                  if(id1.gt.0) then
                     if(id2.gt.0) then
                        d2pyq(ixsym(id1,id2))=dpyq(id1)/gz%yfrem(lm)
                     else
! wildcard in sublattice lm
                        do iw=incffr(lm-1)+1,incffr(lm)
                           d2pyq(ixsym(id1,iw))=dpyq(id1)
                        enddo
                     endif
                  else
! wildcard in sublattice ll
                     if(id2.gt.0) then
                        do iw=incffr(ll-1)+1,incffr(ll)
                           d2pyq(ixsym(id2,iw))=one
                        enddo
                     else
! wildcards in both sublattice ll and lm
                        do iw1=incffr(ll-1)+1,incffr(ll)
                           do iw2=incffr(lm-1)+1,incffr(lm)
                              d2pyq(ixsym(iw1,iw2))=pyq
                           enddo
                        enddo
                     endif
                  endif
               enddo d2pyloop2
            enddo d2pyqloop1
!---- jump here if moded is 0 or 1
150         continue
!            write(*,228)'3X d2pyq 0:',d2pyq
!            write(*,*)'3X Config G 4A: ',phres%gval(1,1)*rtg
!            write(*,154)'3X endmember permutation: ',epermut,(clist(i),i=1,4)
154         format(a,i5,4i4,'--------------------------------')
155         format(a,i5,10i4)
            proprec=>endmemrec%propointer
            emprop: do while(associated(proprec))
               typty=proprec%proptype
               if(typty.ne.1) then
! if property different from 1 (=G) find where to store it, use phmain link
                  do qz=2,lprop-1
                     if(phmain%listprop(qz).eq.typty) goto 170
                  enddo
! a new property, save its typty in listprop and increment lprop
! note that the property index typty is not used as index in gval etc
! as that can be very large. lprop is incremented by 1 for each property
! actually used in the model of the phase.
                  qz=lprop
                  phmain%listprop(qz)=typty
                  if(allocated(phmain%listprop)) then
!                  if(lprop.ge.nprop) then
! VERY STRANGE ERROR, nprop is suddenly zero ....
                     if(lprop.ge.size(phmain%listprop)) then
                        write(*,169)'3X Too many parameter properties ',&
                             lprop,nprop,typty,lokph,&
                             size(phmain%listprop),phlista(lokph)%name
169                     format(a,3i3,2x,2i3,2x,a)
                        gx%bmperr=4338; goto 1000
                     endif
                  else
                     write(*,*)'3X Internal error, listprop not allocated',&
                          lokph,phlista(lokph)%name
                     gx%bmperr=4339; goto 1000
                  endif
                  lprop=lprop+1
                  phmain%listprop(1)=lprop
!                  write(*,*)'3X lprop: ',lprop,typty,phmain%listprop(1)
170               continue
                  ipy=qz
               else
                  ipy=1
               endif
! calculate function and derivatives wrt T and P
! the results from eval_tpfun must also be different in different treads ...
               lokfun=proprec%degreelink(0)
               call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
!               write(*,*)'3X calcg calling eval_tpfun 2: ',gx%bmperr,vals(1)
               if(gx%bmperr.ne.0) goto 1000
               prop1: if(ipy.eq.1) then
! property 1 i.e. Gibbs energy, should be divided by RT
                  vals=vals/rtg
               endif prop1
! debug
!               write(*,173)'3X endmember: ',endmemrec%antalem,ipy,pyq,vals(1)
173            format(a,2i4,4(1pe12.4))
! multiply with py and derivatives. vals is composition independent
!               write(*,*)'3X Config G 4B: ',vals(1)*rtg
               noderz2: if(moded.gt.0) then
                  derloopz2: do id=1,gz%nofc
                     do itp=1,3
                        phres%dgval(itp,id,ipy)=phres%dgval(itp,id,ipy)+ &
                             dpyq(id)*vals(itp)
                     enddo
                     if(moded.gt.1 .and. dpyq(id).gt.zero) then
                        do jd=id+1,gz%nofc
                           phres%d2gval(ixsym(id,jd),ipy)= &
                                phres%d2gval(ixsym(id,jd),ipy)+ &
                                d2pyq(ixsym(id,jd))*vals(1)
                        enddo
                     endif
                  enddo derloopz2
               endif noderz2
               do itp=1,6
                  phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*vals(itp)
               enddo
! strange values of mobilities for ordered phases ...
!               if(ipy.ne.1) then
!               write(*,173)'3X gval:      ',phmain%listprop(ipy),ipy,&
!                    phres%gval(1,ipy),pyq,vals(1)
!               endif
               proprec=>proprec%nextpr
!               write(*,*)'3X Config G 4C: ',phres%gval(1,1)*rtg
! debug problem with mobility calculation
!               if(ipy.eq.2) then
!                  write(*,172)'3X mob: ',ipy,phmain%listprop(1),&
!                       phmain%listprop(ipy),&
!                       pyq,vals(1),pyq*vals(1),phres%gval(1,ipy)
!172               format(a,3i4,6(1pe12.4))
!               endif
            enddo emprop
!------------------------------------------------------------------
! take link to interaction records, use push and pop to save pyq etc
! pmq keeps track of the location in LASTPMQ and MAXPMQ
! for each interaction record in this binary interaction tree
            intrec=>endmemrec%intpointer
            gz%intlevel=0
            pmq=1
! pmq is initiated by palmtree above in the interaction records
!            write(*,*)'3X excess 0: ',associated(intrec),phres%gval(1,1)*rtg
            interloop: do while(associated(intrec))
!----------------------------------------------------------------
! come back here an interaction at a higher level or a poped next that must
! be pushed 
200            continue
               gz%intlevel=gz%intlevel+1
!               write(*,*)'3X excess 1: ',gz%intlevel,phres%gval(1,1)*rtg
               call push_pyval(pystack,intrec,pmq,&
                    pyq,dpyq,d2pyq,moded,gz%nofc)
! intrec%order is initiated by palmtree to set a sequential number
               pmq=intrec%order
!               write(*,155)'3X Pushed: ',pmq,gz%intlevel
!-------------------------------------------------------------------
! come back here for another permutation of same paremeter (no push needed)
220            continue
               bford: if(chkperm) then
                  setipermut: if(maxpmq(pmq).eq.0) then
! ipermut must be initiated and saved in lastpmq
                     ipermut=1; lastpmq(pmq)=ipermut
! On level 1 the number of permutation is in first location
! On level 2 it is more complicated but the first number of perm is in 2nd loc
                     maxpmq(pmq)=intrec%noofip(gz%intlevel)
                  else
! lastpmq and maxpmq already initiated (NOTE: they are used for all
! permutations of the same endmember, that is why they are stored here
! They cannot be pushed on the stack as the stack is also popped
                     ipermut=lastpmq(pmq)+1
                     plimit: if(ipermut.gt.maxpmq(pmq)) then
! maximum interaction level allowed when permutations  is 2
                        level: if(gz%intlevel.eq.1) then
! This is always simple for level 1, 
                           maxpmq(pmq)=maxpmq(pmq)+&
                                intrec%noofip(1)
!                           write(*,155)'3X new limit: ',ipermut,&
!                                maxpmq(pmq)
                           if(ipermut.le.maxpmq(pmq)) goto 230
                        elseif(gz%intlevel.gt.2) then
                           write(*,*)'3X Max level 2 interactions allowed'
                           gx%bmperr=4340; goto 1000
                        else
                           varying: if(intrec%noofip(1).eq.1) then
! If this is 1 then noofip(2) is number of permutations each time
                              maxpmq(pmq)=maxpmq(pmq)+intrec%noofip(2)
                              if(ipermut.le.maxpmq(pmq)) goto 230
                           else
! This is more complicated, different number of permutations each time
! Example: noofip=(3,2,1,0,12) means there are 3 different permutations
! first time; 2 the second time; 1 the last time none;
! 12 is the total number of permutationss (including first order)
! Example 1: end member (A:A:A:A), no permutation
! first int B in 1 with perms:  2nd int C in 2 with perms: (3,3,3,3,12)
! (AB:A:A:A)                   (AB:AC:A:A) (AB:A:AC:A) (AB:A:A:AC)
! (A:AB:A:A)                   (AC:AB:A:A) (A:AB:AC:A) (A:AB:A:AC)
! (A:A:AB:A)                   (AC:A:AB:A) (A:AC:AB:A) (A:A:AB:AC)
! (A:A:A:AB)                   (AC:A:A:AB) (A:AC:A:AB) (A:A:AC:AB)
! Example 2: end member (A:A:A:A), no permutation
! first int B in 1 with perms:  2nd int B in 2 with perms: (3,2,1,0,6)
! (AB:A:A:A)                   (AB:AB:A:A) (AB:A:AB:A) (AB:A:A:AB)
! (A:AB:A:A)                   (A:AB:AB:A) (A:AB:A:AB)
! (A:A:AB:A)                   (A:A:AB:AB)
! (A:A:A:AB)                   none
! If mod(ipermut,noofip(1)) is 0 one should start from index 2
                              nz=intrec%noofip(1)
!                              write(*,155)'3X noofip: ',ipermut,pmq,&
!                                   maxpmq(pmq),(intrec%noofip(j),j=1,nz)
                              if(maxpmq(pmq).gt.0) then
! Previous increase of limit was greater than zero, special case for noofip=2
                                 if(intrec%noofip(1).eq.2) then
                                    maxpmq(pmq)=-maxpmq(pmq)
                                 else
                                    nz=mod(ipermut-1,intrec%noofip(1))
                                    if(nz.eq.0) then
                                       maxpmq(pmq)=-maxpmq(pmq)
                                    else
                                       maxpmq(pmq)=maxpmq(pmq)+&
                                            intrec%noofip(1+nz)
                                    endif
                                 endif
                                 if(ipermut.le.maxpmq(pmq)) goto 230
                              else
! Previous increase of limit was 0, start repeating values from noofip(2..
                                 maxpmq(pmq)=intrec%noofip(2)-&
                                      maxpmq(pmq)
                                 if(ipermut.le.maxpmq(pmq)) goto 230
                              endif
!                              write(*,155)'3X noperm: ',ipermut,pmq,&
!                                   lastpmq(pmq),maxpmq(pmq)
                           endif varying
! as we have passed the limit of permutations, take higher or next interaction
!???                           if(ipermut.le.maxpmq(pmq)) goto 230
                        endif level
! We have exeeded the permutation limit, we should not go to any
! higher interaction but to a next interaction on same level (if any)
! or go down one level
                        if(associated(intrec%highlink)) then
                           if(gz%intlevel.eq.2) then
                              write(*,229)gz%intlevel
229                           format('3X Error, max 2 levels of interactions',/&
                                   ' with permutations!! ',i3)
                              gx%bmperr=4340; goto 1000
                           endif
! Take the link to higher as no more permutations here
                           goto 290
                        endif
!..............................
! No higher level, if we cannot pop we must return to endmember
                        if(gz%intlevel.eq.0) exit interloop
! we must pop lower order interaction records here to get correct permutation
                        call pop_pyval(pystack,intrec,pmq,&
                             pyq,dpyq,d2pyq,moded,gz%nofc)
                        gz%intlevel=gz%intlevel-1
                        pmq=intrec%order
!.................................
! intrec must not be associated in the popint: do-loop
                        nullify(intrec)
                        goto 295
                     endif plimit
! We have now the permutation for this interaction in ipermut
230                  continue
                  endif setipermut
! Found the permutations for option F and B, save it in lastpmq(pmq)
                  lastpmq(pmq)=ipermut
! Without permutations just set ipermut=1
               else
!                  write(*,*)'3X no permutations'
                  ipermut=1
               endif bford
!-------------------------------------------------------------------
! Code below until label 290 the same with and without permutations
! extract  sublattice, constituent and fraction of interacting constituent
! NOTE "ic" used several times below, do not manipulate it!!!
               intlat=intrec%sublattice(ipermut)
               ic=intrec%fraclink(ipermut)
               gz%intlat(gz%intlevel)=intlat
               gz%intcon(gz%intlevel)=ic
               gz%yfrint(gz%intlevel)=phres%yfr(ic)
!               write(*,*)'3X excess 2: ',ionicliq,iliqsave
               if(ionicliq .and. iliqsave) then
                  if(intlat.eq.1 .and. yionva.gt.zero) then
! iliqsave is TRUE for ionic_liquid and for excess parameters without anions
! For cation interactions multiply with yionva.  If no vacancies yionva=-1.0
                     gz%yfrint(gz%intlevel)=phres%yfr(ic)*yionva
!                     write(*,*)'3X *yionva: ',yionva,gz%yfrint(gz%intlevel)
                  endif
               endif
! calculate new PY incl derivatives. Moded to avoid unrequested derivatives
!
! IF interaction endmember is WILDCARD then the interaction is special,
! L(*,A) is y_A *(1-y_A) where 1-y_A is the sum of all fractions except A
! pyq = pyq * y_ic * (y_ix + y_iy + ... ) (all_other_in_same_sublattice))
! derivatives are calculated for all constituents in intlat
! note one can also have wildcards in other sublattices ....
               if(gz%endcon(intlat).gt.0) then
                  wildc=.FALSE.
                  ymult=gz%yfrint(gz%intlevel)
               else
                  if(iliqsave) then
! I sincerely hope wildcards are never used in 2nd subl of ionic liquids ...
        write(*,*)'3X Wildcard in second sublattice illegal for ionic liquids'
                     gx%bmperr=4341; goto 1000
                  endif
                  wildc=.TRUE.
                  write(*,*)'3X wildcard found!'
                  ymult=gz%yfrint(gz%intlevel)*(one-gz%yfrint(gz%intlevel))
               endif
!               write(*,228)'3X d2pyq 1:',d2pyq
!---------------------------------
               noder3A: if(moded.gt.0) then
! ...................................... loop for first derivatives
!                  write(*,228)'3X pyq  0:',pyq
!                  write(*,228)'3X dpyq 1:',dpyq
!                  write(*,228)'3X divers:',ymult,yionva,gz%yfrem(1)
                  iloop1: do id=1,gz%nofc
                     if(moded.gt.1) then
! ...................................... second derivatives
                        if(iliqsave .and. intlat.eq.1) then
! This IF loop is only executed when Va in second sublattice, i.e. when cation
! interactions which should also be multiplied with the power of yionva
! which is gz%intlevel+1
! jonva=phlista(lokph)%i2slx(1) is index of vacancy, i2slx(2) is first neutral
                           if(jonva.le.0) then
                              write(*,*)'3X illegal interaction'
                              gx%bmperr=4342; goto 1000
                           endif
!                           iloop2A: do jd=id,gz%nofc
! for ionic liquids pyq can contain powers of Va, for such a case jd must loop
! from id, not id+1.
!                              d2pyq(ixsym(id,jd))=&
!                                   (gz%intlevel+1)*d2pyq(ixsym(id,jd))*ymult
!                              write(*,209)'3X d2pyq 2: ',id,jd,ixsym(id,jd),&
!                                   gz%intlevel+1,d2pyq(ixsym(id,jd)),ymult
209                           format(a,4i3,6(1pe12.4))
! code above wrong
! CODE BELOW IS UNCERTAIN
                           if(id.eq.ic) then
! derivative wrt new interaction, ic, and Va
                              d2pyq(ixsym(ic,jonva))=&
                                   (gz%intlevel+1)*gz%yfrem(1)*yionva
!                              write(*,209)'3X d2pyq A: ',ic,jonva,&
!                                   ixsym(ic,jonva),gz%intlevel+1,&
!                                   d2pyq(ixsym(ic,jonva)),gz%yfrem(1),yionva
                           elseif(id.eq.jonva) then
! double derivative wrt Va
                              d2pyq(ixsym(jonva,jonva))=&
                 (gz%intlevel+1)*dpyq(jonva)*gz%yfrint(gz%intlevel)/yionva
!                              write(*,209)'3X d2pyq B: ',jonva,jonva,&
!                                   ixsym(jonva,jonva),gz%intlevel+1,&
!                                   d2pyq(ixsym(jonva,jonva)),dpyq(jonva),&
!                                   gz%yfrint(gz%intlevel),yionva
!                           elseif(id.lt.jonva) then
                           elseif(id.le.phlista(lokph)%nooffr(1)) then
! two derivatives for all other anions ...
                              d2pyq(ixsym(id,ic))=dpyq(id)*yionva
!                              write(*,209)'3X d2pyq C: ',id,ic,&
!                                   ixsym(id,ic),0,&
!                                   d2pyq(ixsym(id,ic)),dpyq(id),&
!                                   yionva
                              d2pyq(ixsym(id,jonva))=(gz%intlevel+1)*ymult
!                              write(*,209)'3X d2pyq D: ',id,jonva,&
!                                   ixsym(id,jonva),gz%intlevel+1,&
!                                   d2pyq(ixsym(id,jonva)),ymult
! no 2nd derivatves wrt neutrals
                           endif
!                           write(*,228)'3X d2pyq 7:',d2pyq
!                           enddo iloop2A
! END VERY UNCERTAIN CODE, the rest seems to work ...
                        else
! For all other models and cases it is simply ...
                           iloop2B: do jd=id+1,gz%nofc
                              d2pyq(ixsym(id,jd))=d2pyq(ixsym(id,jd))*ymult
                           enddo iloop2B
                        endif
! NOTE "ic" has been set above as the interacting constituent
                        if(iliqsave) then
                           if(intlat.eq.1 .and. yionva.gt.zero) then
                              continue
! This already done ...
! For ionic liquid model the 2nd derivatives must be multipled with yionva
!                              if(id.eq.phlista(lokph)%i2slx(1)) then
! This is the vacancy, all 2nd derivatives multiplied with a factor
!                                 d2pyq(ixsym(id,ic))=&
!                                      (gz%intlevel+1)*d2pyq(ixsym(id,ic))
!                               write(*,209)'3X d2pyq 3: ',id,ic,ixsym(id,ic),&
!                                      gz%intlevel+1,d2pyq(ixsym(id,ic))
!                              endif
                           else
                              d2pyq(ixsym(id,ic))=dpyq(id)*yionva
!                              write(*,209)'3X d2pyq 4: ',id,ic,ixsym(id,ic),&
!                                   0,d2pyq(ixsym(id,ic)),yionva
                           endif
                        else
                           d2pyq(ixsym(id,ic))=dpyq(id)
                        endif
                     endif
! ................................. this is the first derivative, must be exact
! very messy for the ionic liquid here ...
                     dpyq(id)=dpyq(id)*ymult
                     if(ionicliq .and. iliqsave) then
!                        write(*,202)'3X Extra va power: ',id,gz%intlevel,&
!                             gz%intlat(gz%intlevel),dpyq(id)
202                     format(a,3i3,4(1pe12.4))
                        if(id.eq.phlista(lokph)%i2slx(1) .and. &
                             gz%intlat(gz%intlevel).eq.1) then
! for vacancies there is an additional factor if interaction in first subl
                           dpyq(id)=(gz%intlevel+1)*dpyq(id)
!                           write(*,197)gz%intlevel,gz%intcon(gz%intlevel)
197                        format('3X: Va inter: ',5i3)
                        endif
                     endif
                  enddo iloop1
! we must check if any endmember is wildcard like L(*:A,B)
! Hopefully this works also for ionic liquid interaction between neutrals
                  do ll=1,msl
                     if(ll.ne.intlat) then
                        if(gz%endcon(ll).lt.0) then
                           do iw=incffr(ll-1)+1,incffr(ll)
                              d2pyq(ixsym(iw,ic))=pyq
                           enddo
                        endif
                     endif
                  enddo
                  wildcard: if(wildc) then
! The interacting constituent is a wildcard ... calculate the contribution
! to second derivate from all fractions in intlat, remember incffr(0)=0.
! Ionic liquids should never have wildcards as intercations ... ?
                     do iw=incffr(intlat-1)+1,incffr(intlat)
                        if(iw.ne.ic) then
                           d2pyq(ixsym(iw,ic))=dpyq(iw)
                        endif
!                        write(*,213)'3X 529: ',iw,ic,ixsym(iw,ic),&
!                             gz%intlevel,intlat,incffr(intlat)
                        dpyq(iw)=pyq*gz%yfrint(gz%intlevel)
!                        dpyq(jd)=pyq*gz%yfrint(gz%intlevel)
                     enddo
213                  format(a,10i5)
                     dpyq(ic)=pyq*(one-gz%yfrint(gz%intlevel))
                  else
! this is the normal first derivative of pyq*y(ic) with respect to y(ic)=ymult
                     dpyq(ic)=pyq
                     if(ionicliq) then
!                        write(*,214)'3X Multiply with y_va: ',&
!                             iliqsave,ic,intlat,yionva,pyq
214                     format(a,l,2i3,4(1pe12.4))
                        if(iliqsave .and. intlat.eq.1.and.yionva.gt.zero) then
! for compatibility with substitutional liquids, multiply interactions 
! of cations (in 1st subl) when vacancies in 2nd with the vacancy fraction
                           dpyq(ic)=pyq*yionva
                        endif
                     endif
                  endif wildcard
!                  write(*,228)'3X dpyq: ',(dpyq(ll),ll=1,4)
228               format(a,6(1pe12.4))
               endif noder3A
! pyq calculated identically for wildcards as ymult set differently above
! It should work for ionic liquids as ymult has been multiplied with yionva
               pyq=pyq*ymult
               proprec=>intrec%propointer
!               write(*,218)'3X pyq: ',associated(proprec),ymult,pyq
218            format(a,l,2(1pe12.4))
! list values of pyq, dpyg, d2pyg
!               write(*,228)'3X pyq:',pyq
!               write(*,228)'3X dpy:',dpyq
!               write(*,228)'3X d2py:',d2pyq
219            format(a,6(1pe12.4))
!..............................
               intprop: do while(associated(proprec))
! calculate interaction parameter, can depend on composition
                  call cgint(lokph,proprec,moded,vals,dvals,d2vals,gz,ceq)
                  if(gx%bmperr.ne.0) goto 1000
!                  write(*,228)'3X val:',vals(1),(dvals(1,id),id=1,gz%nofc)
! G parameters (ipy=1) are divided by RT inside cgint
                  typty=proprec%proptype
                  if(typty.ne.1) then
! other properties than 1 (G) must be stored in different gval(*,ipy) etc
                     do qz=2,lprop-1
                        if(phmain%listprop(qz).eq.typty) goto 250
                     enddo
! a new property, save its typty in listprop and increment lprop
                     qz=lprop
                     phmain%listprop(qz)=typty
                     lprop=lprop+1
                     phmain%listprop(1)=lprop
250                  continue
! here the value of ipy is set, 1 means G
                     ipy=qz
                  else
                     ipy=1
                  endif
! note: adding to phres%gval at the end of noder4: if(....)
                  noder4: if(moded.gt.0) then
                     iloop3: do id=1,gz%nofc
                        if(moded.gt.1) then
!                           iloop4: do jd=id+1,gz%nofc
! This loop was constructed for normal cases when pyq has each fraction once
! in ionic liquids Va can have a power so loop for all!
                           iloop4: do jd=id,gz%nofc
                              phres%d2gval(ixsym(id,jd),ipy)= &
                                   phres%d2gval(ixsym(id,jd),ipy)+ &
                                   d2pyq(ixsym(id,jd))*vals(1)
!                              write(*,251)'3X G:',id,jd,ixsym(id,jd),&
!                                   d2pyq(ixsym(id,jd)),vals(1)
251                           format(a,3i3,4(1pe12.4))
                           enddo iloop4
                        endif
                        do itp=1,3
                           phres%dgval(itp,id,ipy)=&
                                phres%dgval(itp,id,ipy)+ &
                                dpyq(id)*vals(itp)
                        enddo
                     enddo iloop3
!                     write(*,211)'3X Interactions: ',gz%iq,jonva
211                  format(a,5i3,5x,i3)
                     if(jonva.gt.0) then
!                        write(*,212)jonva,phres%dgval(1,jonva,1)*rtg
212                     format('3X with va: ',i3,6(1pe12.4))
                     endif
!...............................
! below contribution to derivatives from composition dependent parameters
! the values of gz%iq represent interacting constituents and are set in cgint
                     cdex1: if(gz%iq(5).gt.0) then
! gz%iq(5) is nonzero only for TOOP and similar models not implemented yet ...
                        gx%bmperr=4086; goto 1000
                     elseif(gz%iq(4).gt.0) then
!...............................
! composition dependent reciprocal parameter
! for ionic liquid one must consider extra vacancy fractions ...
! remember ipy is property type for this parameter, set above
!                        write(*,333)'3X comp dep reciprocal:',gz%iq,pyq,vals(1)
333                     format(a,5i4,4(1pe14.6))
                        if(moded.gt.0) then
                           do jk=1,4
                           if(moded.gt.1) then
! contribution to second derivatives with respect to 2 const previously ignored
! No second derivatives calculated in cgint for this case
                           do qz=jk,4
                              phres%d2gval(ixsym(gz%iq(jk),gz%iq(qz)),ipy)=&
                                  phres%d2gval(ixsym(gz%iq(jk),gz%iq(qz)),ipy)+&
                                  dpyq(gz%iq(jk))*dvals(1,gz%iq(qz))+&
                                  dpyq(gz%iq(qz))*dvals(1,gz%iq(jk))
                           enddo
                           endif
! first derivatives, including 2nd wrt T and P
                              do itp=1,3
! itp=1 for 1st derivative, =2 for 2nd derivative also with T, =3 also with P
                                 phres%dgval(itp,gz%iq(jk),ipy)=&
                                      phres%dgval(itp,gz%iq(jk),ipy)+&
                                      pyq*dvals(itp,gz%iq(jk))
                              enddo
                           enddo
                        endif
                     elseif(gz%iq(3).gt.0) then !cedex1
! composition dependent ternary interaction in same sublattice, Mats model
! PROBABLY ERRORS HERE as no consideration of derivatives wrt other endmember
! constituents, only to the 3 interacting
! ALSO used to indicate derivatives wrt vacancies in ionic liquid model ??? NO!
!...<<<<<<<...... indentation back 2 levels
                  if(moded.gt.1) then
                     noindent1: do jk=1,3
                        do qz=jk+1,3
! the second derivative for jk=qz calculated below as it is simpler
                           phres%d2gval(ixsym(gz%iq(jk),gz%iq(qz)),ipy)=&
                                phres%d2gval(ixsym(gz%iq(jk),gz%iq(qz)),ipy)+&
                                dpyq(gz%iq(jk))*dvals(1,gz%iq(qz))+&
                                dpyq(gz%iq(qz))*dvals(1,gz%iq(jk))
                        enddo
                     enddo noindent1
                  endif
                  do jk=1,3
                     do itp=1,3
                        phres%dgval(itp,gz%iq(jk),ipy)=&
                             phres%dgval(itp,gz%iq(jk),ipy)&
                             +pyq*dvals(itp,gz%iq(jk))
                     enddo
                     phres%d2gval(ixsym(gz%iq(jk),gz%iq(jk)),ipy)=&
                          phres%d2gval(ixsym(gz%iq(jk),gz%iq(jk)),ipy)+&
                          2.0D0*dpyq(gz%iq(jk))*dvals(1,gz%iq(jk))
                  enddo
!...>>>>>>...........indentation forward
                     elseif(gz%iq(2).gt.0) then !cedex1
! gz%iq(2) nonzero means composition dependent binary interaction parameter,
! only RK yet.
                        noder3B: if(moded.gt.1) then
! one can maybe make this loop faster by just looping throungh endmembrs
! but then one must handle wildcard endmembers ....
! and there may be other bugs here anyway ....
                           do ic1=1,gz%nofc
                              add1=dpyq(ic1)*dvals(1,gz%iq(1))+&
                                   dpyq(gz%iq(1))*dvals(1,ic1)+&
                                   pyq*d2vals(ixsym(ic1,gz%iq(1)))
                              phres%d2gval(ixsym(ic1,gz%iq(1)),ipy)=&
                                   phres%d2gval(ixsym(ic1,gz%iq(1)),ipy)+add1
                              if(ic1.ne.gz%iq(1)) then
! this IF to avoid that the second derivative gz%iq(1) and gz%iq(2) is
! calculated twice. ic1 will at some time be equal to gz%iq(1) and to gz%iq(2)
                                 add1=dpyq(ic1)*dvals(1,gz%iq(2))+&
                                      dpyq(gz%iq(2))*dvals(1,ic1)+&
                                      pyq*d2vals(ixsym(ic1,gz%iq(2)))
                                 phres%d2gval(ixsym(ic1,gz%iq(2)),ipy)=add1+&
                                      phres%d2gval(ixsym(ic1,gz%iq(2)),ipy)
                              endif
                           enddo
                        endif noder3B
                        do itp=1,3
                           phres%dgval(itp,gz%iq(1),ipy)=&
                                phres%dgval(itp,gz%iq(1),ipy)&
                                +pyq*dvals(itp,gz%iq(1))
                           phres%dgval(itp,gz%iq(2),ipy)=&
                                phres%dgval(itp,gz%iq(2),ipy)&
                                +pyq*dvals(itp,gz%iq(2))
                           if(ionicliq) then
! for ionic liquid when interactions involve cations there is a contribution
! due to the vacancy fraction multiplied with the cations yc1*yc2*yva**2
! we are dealing with binary RK interactions, gz%intlevel=1, check if 
! interaction is in first sublattice (between cations) and vacancy in second
                              if(iliqva .and. jonva.gt.0) then
                                 if(gz%intlat(1).eq.1) then
! add pyq multipled with the derivative with respect to vacancy fraction
! This should be done for d2gval also but I skip that at present ...
                                    phres%dgval(itp,jonva,ipy)=&
                                         phres%dgval(itp,jonva,ipy)+&
                                         pyq*dvals(itp,jonva)
!                                 write(*,*)'3X jonva:',jonva,pyq,dvals(1,jonva)
                                 elseif(gz%intlat(1).eq.2 .and. &
                                      gz%iq(2).gt.jonva) then
! This fixed the problem with Pd-Ru-Te in the fuel (+ fix in cgint)
!                 write(*,55)'3X (C:Va,K)',iliqva,gz%intlat(1),jonva,gz%iq(1),&
!                      gz%iq(2),gz%endcon(1),pyq,dvals(itp,gz%endcon(1))
55               format(a,l2,5i3,4(1pe12.4))
                                    icat=gz%endcon(1)
                                    phres%dgval(itp,icat,ipy)=&
                                         phres%dgval(itp,icat,ipy)+&
                                         pyq*dvals(itp,icat)
                                 endif
                              endif
                           endif
                        enddo
                     endif cdex1
! end contribution to derivates from composition dependent parameters
!......................
                  endif noder4
! finally add the contribution to G, G.T etc
                  iloop6: do itp=1,6
                     phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*vals(itp)
                  enddo iloop6
! debug problem with mobility calculation
!                  if(ipy.eq.2) then
!                     write(*,172)'3X imob:',ipy,phmain%listprop(1),&
!                          phmain%listprop(ipy),&
!                          pyq,vals(1),pyq*vals(1),phres%gval(1,ipy)
!                  endif
                  proprec=>proprec%nextpr
               enddo intprop
!               write(*,*)'3X Config G 4F: ',phres%gval(1,1)*rtg
! finished one interaction (or permutation on this level), go to higher level
! note that ipermut is saved in lastpmq(pmq).  If there are more
! permutations on this level they will be calculated later also including 
! higher order parameters.  
!------------------------------------------------------------------
! Take link to higher level records for current permutation
290            continue
               intrec=>intrec%highlink
               wrong: if(chkperm .and. associated(intrec)) then
! We must go to higher as we can have interactions with different permutations?
                  jpr=intrec%order
                  if(lastpmq(jpr).gt.0 .and.lastpmq(jpr).ge.maxpmq(jpr)) then
! if we nullify here we will take next rather than higher
!                     nullify(intrec)
!                     write(*,155)'3X Maybe skipping higer?: ',jpr,&
!                          lastpmq(jpr),maxpmq(jpr),gz%intlevel
!                     if(maxpmq(jpr).lt.0) maxpmq(jpr)=intrec%noofip(2)-&
!                          maxpmq(jpr)
                  endif
               endif wrong
! if intrec is associated then go to big "interloop: do while()" loop
295            continue
               popint: do while(.not.associated(intrec))
! No higher level, pop lower order interaction records, if no pop: endmember
                  if(gz%intlevel.eq.0) exit interloop
                  call pop_pyval(pystack,intrec,pmq,&
                       pyq,dpyq,d2pyq,moded,gz%nofc)
                  gz%intlevel=gz%intlevel-1
                  pmq=intrec%order
! check if we have more permutations for this record
                  if(chkperm) then
                     if(lastpmq(pmq).lt.maxpmq(pmq)) then
                        goto 200
                     endif
                  endif
                  intrec=>intrec%nextlink
               enddo popint
! we should loop here if we found a higher order record or 
! a lower order record with a next link
            enddo interloop
298         continue
!            write(*,*)'3X Config G 4X: ',phres%gval(1,1)*rtg
! take next permutation of the end member fractions
         enddo empermut
300      continue
! take next end member
!      write(*,155)'3X endmem: ',epermut,endmemrec%noofpermut,endmemrec%antalem
         endmemrec=>endmemrec%nextem
      enddo endmemloop
!      write(*,*)'3X Config G 5: ',phres%gval(1,1)*rtg
!------------------------------------------------------------------------
! end loop for this fraction type, initiation for next in the beginning of loop
! but we may have to calculate once again with same fraction type but
! with the fractions as disordered fractions
!      write(*,*)'3X Testing nevertwice ',nevertwice
! Jump to 400 terminates calculation for this fraction type
!      write(*,303)'3X Nevertwice: ',nevertwice,&
!           btest(phlista(lokph)%status1,phsubo),&
!           first,fractype,phres%gval(1,1)
303   format(a,3(1x,L),i3,4(1pe12.4))
!      write(*,623)'3X order/disorder: ',lprop,phres%gval(1,2),phres%gval(1,3)
      if(nevertwice) goto 400
! UNIFINISHED ??
! TEST IF WE SHOULD SUBTRACT THE ORDERED ENERGY AS DISORDERED AS IN THE
! CURRENT IMPLEMENTATION IN THERMO-CALC. BY JUMPING TO 400 WE SKIP THAT.
      if(btest(phlista(lokph)%status1,phsubo)) then
!         write(*,*)'3X phsubo bit set'
         goto 400
      endif
! PARTITION PROBLEM FOR ORDERED PHASES
!      goto 400
!------------------------------------------------
!      write(*,611)'3X ftyp1:',fractype,btest(phlista(lokph)%status1,phmfs),&
!           btest(phmain%status2,csorder),first,lokph,phres%gval(1,1)
611   format(a,i3,3(1x,L),i3,3(1pe12.4))
      disord: if(fractype.eq.1 .and. btest(phlista(lokph)%status1,phmfs) &
           .and. btest(phmain%status2,csorder)) then
! Handle additions of several fraction set ?? Additions calculated
! after both ordered and disordered fraction set calculated
!         write(*,611)'3X ftyp:',fractype,btest(phlista(lokph)%status1,phmfs),&
!              btest(phmain%status2,csorder),first,lokph,phres%gval(1,1)
         if(first) then
! calculate with next fraction type
! alternative meithod: no need to calculate with all fractions as disordered
            first=.false.
!            write(*,*)'3X: next fraction type'
!            goto 400
! we must save phres%yfr before disorder ....
!            allocate(savey(gz%nofc))
! this creates problem with pointers in disordery?? avoid by allocating savey??
!            savey=phres%yfr
            do j1=1,size(phres%yfr)
               savey(j1)=phres%yfr(j1)
            enddo
!------------ code below was removed for a while but is now reinstated
!            write(*,*)'3X cg: ',phmain%phlink,phmain%disfra%varreslink
! ??? very uncertain how to call disordery .....
!            call disordery(phmain,phmain%disfra%varreslink,ceq)
!            write(*,*)'3X At disordery: ',phmain%disfra%varreslink,&
!                 cps%disfra%varreslink
            call disordery(phmain,ceq)
! if call to disordery here no crash in disordery ...
! if call moved to after assignment of savey there is a crash (GNU fortran)
!----------
!            allocate(savey(gz%nofc))
!            savey=phres%yfr
!            nprop=phmain%nprop
! error calculating volumes for order/disorder, V0 in gval(1,2), VA in gval(1,3)
!            write(*,623)'3X V0,VA 1: ',lprop,phres%gval(1,2),phres%gval(1,3)
623         format(a,i3,6(1pe12.4))
! we already know nprop
            allocate(saveg(6,nprop))
            allocate(savedg(3,gz%nofc,nprop))
            allocate(saved2g(nofc2,nprop))
!            write(*,*)'3X saveg allocated 3: ',size(saveg)
            saveg=phres%gval
            savedg=phres%dgval
            saved2g=phres%d2gval
!            do i1=1,gz%nofc
!               write(*,602)'3X G4y: ',i1,phres%dgval(1,i1,1),savedg(1,i1,1)
!            enddo
            phres%gval=zero
            phres%dgval=zero
            phres%d2gval=zero
            goto 110
         else
! We have now calculated the 4SL model both as original and disordered
! We should now subtract the disordered from the ordered
! this is debug output
!            do i1=1,gz%nofc
!               write(*,602)'3X G4x: ',i1,phres%dgval(1,i1,1),savedg(1,i1,1)
!            enddo
602         format(a,i3,6(1pe14.6))
! Ordered part calculated with disordered fractions, subtract this
! from the first, restore fractions and deallocate
! THIS IS TRICKY
! NOTE all sublattices are identical in this case with the same number 
! of constituents
! First sum all second derivatives into tmpd2g, moded=1 means only 1st deriv
! error calculating volumes for order/disorder, V0 in gval(1,2), VA in gval(1,3)
!             write(*,623)'3X V0,VA 2: ',lprop,phres%gval(1,2),phres%gval(1,3)
            noder6A: if(moded.gt.1) then
               nz=fracset%tnoofxfr
               allocate(tmpd2g(nz*(nz+1)/2,nprop))
               tmpd2g=zero
! remove this comment to obtain old code
!               goto 666
!--------------------------------------------------------------------------
! now summation derived 2017-02-20
! d2G/dy_is dy_jt = a_s a_t (d2G/dx_is dx_js + d2G/dx_is dx_jt +
!                            d2G/dx_it dx_js + d2G/dx_it dx_jt)
! first calculate the term within ( )
               do ipy=1,lprop-1
                  do i1=1,gz%nofc
! note j1 and j2 are set to constituent index of i1,i2 in disordered sublattice
! or second sublattice if interstital ?? Maybe not
                     j1=fracset%y2x(i1)
                     do i2=i1,gz%nofc
                        j2=fracset%y2x(i2)
                        tmpd2g(ixsym(j1,j2),ipy)=tmpd2g(ixsym(j1,j2),ipy)+&
                             phres%d2gval(ixsym(i1,i2),ipy)
                     enddo
                  enddo
               enddo
! Then subtract this from saved2g
               do ipy=1,lprop-1
                  do i1=1,gz%nofc
                     j1=fracset%y2x(i1)
                     do i2=i1,gz%nofc
                        j2=fracset%y2x(i2)
!                        write(*,637)'3X sublf: ',ipy,i1,i2,nsit1,nsit2,&
!                             gz%nofc,sublf
                        phres%d2gval(ixsym(i1,i2),ipy)=&
                             saved2g(ixsym(i1,i2),ipy)-&
                             tmpd2g(ixsym(j1,j2),ipy)*&
                             fracset%dxidyj(i1)*fracset%dxidyj(i2)
                     enddo
                  enddo
               enddo
               goto 667
!----------------------- old code below not used
!               if(nsl.gt.3) goto 666
! probable BUG here with 2nd derivatives of ordered FCC calculated as disordered
! But this is necessary for 2 sublattice ordered model !! ??
! DEBUG, problem with partitioning
!               write(*,613)'3X sub: ',nz,gz%nofc,fracset%latd,fracset%y2x
!613            format(a,3i3,2x,20i3)
!               write(*,614)'3X dxi/dyj: ',fracset%dxidyj
!614            format(a,(10f7.4))
! Formula derived 2017-02-20
! d2G/dy_isdy_jt = a_s a_t (d2G/dx_isdx_is + 2 d2G/dx_isdx_jt + d2G/dx_jsdx_jt)
! THIS IS ONLY DONE FOR ORDERING ON 2 SUBLATTICES
               do ipy=1,lprop-1
                  do i1=1,gz%nofc
                     j1=fracset%y2x(i1)
                     do i2=i1,gz%nofc
                        j2=fracset%y2x(i2)
                        tmpd2g(ixsym(j1,j2),ipy)=tmpd2g(ixsym(j1,j2),ipy)+&
                             phres%d2gval(ixsym(i1,i2),ipy)
                     enddo
                  enddo
               enddo
! phres%d2gval(ixsym(i1,i2),ipy) is 2nd derivatives of the orderd part
! calculated with the same fractions in all sublattices.
!               write(*,603)'3X d2Gord(x)/dy1/dy2: ',nz,tmpd2g(1,1),&
!                    tmpd2g(2,1),tmpd2g(3,1)
603            format(a,i3,6(1pe12.4))
! tmpd2g is now d2G/dxidxj calculated with disordered fractions
! subract that from saved d2G/dyidyj saved in saved2g taking into account
! the derivatives dxi/dyj (in fracset%dxidyj)
! original
!               do ipy=1,lprop-1
!                  do i1=1,gz%nofc
!                     j1=fracset%y2x(i1)
!                     do i2=i1,gz%nofc
! subtract from saved value
!                        j2=fracset%y2x(i2)
!                        phres%d2gval(ixsym(i1,i2),ipy)=&
!                             saved2g(ixsym(i1,i2),ipy)
! just ignore contribution to 2nd derivatives from ordered calc as disordered??
!                             saved2g(ixsym(i1,i2),ipy)-&
!                             tmpd2g(ixsym(j1,j2),ipy)*&
!                             fracset%dxidyj(i1)*fracset%dxidyj(i2)
!                     enddo
!                  enddo
!               enddo
666            continue
! subract the corresponding d2G/dy1/dy2 calculated as disordered ...
! but multiplied with the product of the sites on the sublattices for
! the constituents used as derivatives
               do ipy=1,lprop-1
                  nsit1=1
                  do i1=1,gz%nofc
!                     j1=fracset%y2x(i1)
! incffr is the number of constituents up to and including nsit1
                     if(i1.gt.incffr(nsit1)) nsit1=nsit1+1
                     nsit2=nsit1
                     do i2=i1,gz%nofc
! subtract from saved value
!                        j2=fracset%y2x(i2)
                        if(i2.gt.incffr(nsit2)) nsit2=nsit2+1
                        sublf=phres%sites(nsit1)*phres%sites(nsit2)
!                        write(*,637)'3X sublf: ',ipy,i1,i2,nsit1,nsit2,&
!                             gz%nofc,sublf
637                     format(a,6i3,1pe12.4)
! evidently this attempt to improve made it worse ...
!                        sublf=one
                        phres%d2gval(ixsym(i1,i2),ipy)=&
                             saved2g(ixsym(i1,i2),ipy)-&
                             phres%d2gval(ixsym(i1,i2),ipy)
!                             sublf*phres%d2gval(ixsym(i1,i2),ipy)
! maybe not needed ??
!                             phres%d2gval(ixsym(i1,i2),ipy)*&
!                             fracset%dxidyj(i1)*fracset%dxidyj(i2)
                     enddo
                  enddo
               enddo
667            continue
               if(allocated(tmpd2g)) deallocate(tmpd2g)
            endif noder6A
!---------------------
! sum all first partial derivates to first sublattice
            noder6B: if(moded.gt.0) then
!               write(*,613)'3X dG/dx: ',fracset%ndd,fracset%nooffr
               do ipy=1,lprop-1
                  do ider=1,3
                     do is=1,fracset%nooffr(1)
                        sum=zero
                        kk=is
                        do ll=1,fracset%latd
                           sum=sum+phres%dgval(ider,kk,ipy)
! it is not really necessary to put phres%dgval it to zero, just for prudence
!                           phres%dgval(ider,kk,ipy)=zero
                           kk=kk+fracset%nooffr(1)
                        enddo
                        phres%dgval(ider,is,ipy)=sum
                     enddo
                     if(fracset%ndd.eq.2) then
! one can have 2 sets of ordered subl like (Al,Fe)(Al,Fe)...(C,Va)(C,Va)...
! BUT I doubt that works ...
                        ioff=fracset%nooffr(1)*fracset%latd
                        do is=1,fracset%nooffr(2)
                           sum=zero
                           kk=ioff+is
                           do ll=fracset%latd+1,phlista(lokph)%noofsubl
                              sum=sum+phres%dgval(ider,kk,ipy)
                              phres%dgval(ider,kk,ipy)=zero
                              kk=kk+fracset%nooffr(2)
                           enddo
                           phres%dgval(ider,ioff+is,ipy)=sum
                        enddo
                     endif
                  enddo
               enddo
!-------------------------
               if(moded.gt.0) then
                  do ipy=1,lprop-1
! loop in negative direction avoid destroy the values in phres%dgval first subl
                     do i1=gz%nofc,1,-1
! all derivatives wrt same element from all sublattices is in first sublattice
                        j1=fracset%y2x(i1)
                        do ider=1,3
! Finally subtract this contribution from saved values
!                           phres%dgval(ider,i1,ipy)=savedg(ider,i1,ipy)-&
                           xxx=savedg(ider,i1,ipy)-&
                                phres%dgval(ider,j1,ipy)*fracset%dxidyj(i1)
!                           write(*,615)'3X Gy-Gx: ',ider,i1,ipy,j1,&
!                                savedg(ider,i1,ipy),phres%dgval(ider,j1,ipy),&
!                                fracset%dxidyj(i1),xxx
!615                        format(a,4i3,4(1pe14.6))
                           phres%dgval(ider,i1,ipy)=xxx
                        enddo
                     enddo
                  enddo
               endif
            endif noder6B
! check for bug, phres%gval(1,1) must not be negative!!
!            write(*,617)'3X do=o-oasd: ',saveg(1,1),phres%gval(1,1),&
!                 saveg(1,1)-phres%gval(1,1)
617         format(a,6(1pe12.4))
            do ipy=1,lprop-1
               do ider=1,6
                  phres%gval(ider,ipy)=saveg(ider,ipy)-&
                       phres%gval(ider,ipy)
               enddo
            enddo
! error calculating volumes for order/disorder, V0 in gval(1,2), VA in gval(1,3)
!            write(*,623)'3X V0,VA 3: ',lprop,phres%gval(1,2),phres%gval(1,3)
! restore ordered fractions and deallocate save arrays why not allocate savey?
!            write(*,612)'3X yd: ',(phres%yfr(ipy),ipy=1,gz%nofc)
!            do ipy=1,gz%nofc
            phres%yfr=savey
!            enddo
!            write(*,612)'3X yo: ',(phres%yfr(ipy),ipy=1,gz%nofc)
612         format(a,6(1pe11.3)/(7x,6e11.3))
! why set to zero if I deallocate ??
!            savey=zero
!            saveg=zero
!            savedg=zero
!            saved2g=zero
!            if(ocv()) write(*,*)'3X saveg DE-allocated 1: ',size(saveg)
!            deallocate(savey)
            deallocate(saveg)
            deallocate(savedg)
            deallocate(saved2g)
         endif
! code above reinstated but has problems ....
      endif disord
! WE CAN JUMP HERE WITHOUT CALCULATING THE ORDERED PART AS DISORDERED
400 continue
   enddo fractyp
!--------------------------------------------------------------
! finished loops for all fractypes, now add together G and all
! partial derivatives for all fractypes
410 continue
! cheking for properties
!   if(ocv()) then
!      write(*,411)lprop-1,(phmain%listprop(j1),j1=2,lprop)
!      write(*,412)'Val: ',(phmain%gval(1,j1),j1=1,lprop-1)
!411   format('3X Properties: ',i3,': ',10i4)
412   format(a,(6E12.4))
!   endif
   norfc=phlista(lokph)%tnooffr
   fractionsets: if(btest(phlista(lokph)%status1,phmfs)) then
!----------------------------------------------------------------
! for disordered part of sigma we may have to multiply the disordered
! part with fsites to have correct formula unit
!      write(*,*)'3X fsites 1: ',phmain%disfra%fsites
      fsites=phmain%disfra%fsites
! add together contributions from different fractypes
! phres is last calculated part, set phpart to ordered part (phmain)
      phpart=>phmain
! loop for all second and first derivatives using chain rule
! and coefficients from fracset%dxidyj
! d2f1/dyidyj = d2f2/dxkdxl*dxk/dyi*dxl/dyj
! gz%nofc are number of disordered constituents
! norfc are number of ordered constituents
! lprop-1 is number of properties to be summed
! G(tot)    = GD(x)+(GO(y)-GO(y=x))
! G(tot).yj = dGD(x).dxi*dxdyj + (GO(y).yj - GO(y=x).yj)
! configurational entropy calculated only for GO(y)
      noder7A: if(moded.gt.0) then
         do i1=1,norfc
            j1=fracset%y2x(i1)
! second derivatives
            noder7B: if(moded.gt.1) then
               do i2=i1,norfc
! add the contributions from the disordered part
                  j2=fracset%y2x(i2)
                  do ipy=1,lprop-1
                     phpart%d2gval(ixsym(i1,i2),ipy)=&
                          phpart%d2gval(ixsym(i1,i2),ipy)+&
                          fsites*phres%d2gval(ixsym(j1,j2),ipy)*&
                          fracset%dxidyj(i1)*fracset%dxidyj(i2)
                  enddo
               enddo
            endif noder7B
! first derivatives
            do ipy=1,lprop-1
!               add1=phpart%dgval(1,i1,ipy)
               do ider=1,3
!                  phpart%dgval(ider,i1,ipy)=phpart%dgval(ider,i1,ipy)+&
                  xxx=phpart%dgval(ider,i1,ipy)+&
                    fsites*phres%dgval(ider,j1,ipy)*fracset%dxidyj(i1)
! phres have the disordred contribution
!                  write(*,413)'3X Gd+Go:',ider,i1,j1,&
!                       phpart%dgval(ider,i1,ipy),fsites,&
!                       phres%dgval(ider,j1,ipy),fracset%dxidyj(i1),xxx
!                  write(*,413)'3X Gd+Go:',ider,i1,j1,&
!                       phmain%dgval(ider,i1,ipy),fsites,&
!                       phres%dgval(ider,j1,ipy),fracset%dxidyj(i1),xxx
                  phpart%dgval(ider,i1,ipy)=xxx
               enddo
            enddo
         enddo
      endif noder7A
413   format(a,3i3,6(1pe12.4))
! Check Integral values, phpart%gval(1,1) is ordered-ordasdis, phres is disord
!      write(*,617)'3X g=do+d:   ',phpart%gval(1,1),fsites*phres%gval(1,1),&
!           phpart%gval(1,1)+fsites*phres%gval(1,1)
      do ipy=1,lprop-1
!         add1=phpart%gval(1,ipy)
         do ider=1,6
            phpart%gval(ider,ipy)=phpart%gval(ider,ipy)+&
                 fsites*phres%gval(ider,ipy)
                 
         enddo
!         if(ocv()) write(*,413)'3X G:',ipy,0,0,&
!         write(*,413)'3X G:',ipy,0,0,&
!              phpart%gval(1,ipy),add1,phres%gval(1,ipy)
      enddo
!      write(*,413)'3X 413:',ipy,0,0,&
!           phpart%gval(1,ipy),add1,phres%gval(1,1)
   endif fractionsets
! now set phres to ordered+disorded results and forget phpart
   phres=>phmain
!................................
!   write(*,*)'3X: ioliq+saved: ',ionicliq,iliqsave,phres%gval(1,1)
   ionliqsum: if(ionicliq .and. iliqsave) then
! For ionic liquid we may have to add gsave+Q*gval (with chain rule ...)
! G = saveg + Q*phres%gval with 1st and 2nd derivatives
! NOT FINISHED !!! interaction parameters above with VA must be treated
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! BEWHARE: FOR IONIC_LIQUID Thermo-Calc (version S) calculates G = Q G_M 
! if there are no end-member parameters (G_M is the Gibbs energy per
! formula unit and Q is the number of sites in second sublattice), 
! This is wrong (but all endmember parameters are never zero for a real liquid)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!      write(*,*)'3X Config G 6: ',phres%gval(1,1)*rtg
      if(moded.eq.0) goto 490
!      write(*,491)'3X ionliq: ',phlista(lokph)%i2slx,phlista(lokph)%nooffr
491   format(a,2i3,5x,2i3)
      firstd: do i1=1,norfc
         secondd: do i2=i1,norfc
            do ipy=1,lprop-1
!               write(*,497)'3X adding: ',i1,i2,ixsym(i1,i2),ipy
497            format(a,10i3)
               phres%d2gval(ixsym(i1,i2),ipy)=saved2g(ixsym(i1,i2),ipy)+&
                    phres%sites(2)*phres%d2gval(ixsym(i1,i2),ipy)
               add1=zero
! IMPORTANT note dpqdy(i1) the the charge of iq, do not confuse with dpyq ...
               if(i1.le.phlista(lokph)%nooffr(1)) then
                  add1=phres%dpqdy(i1)*phres%dgval(1,i2,ipy)
               endif
               if(i2.le.phlista(lokph)%nooffr(1)) then
                  add1=add1+phres%dpqdy(i2)*phres%dgval(1,i1,ipy)
               endif
               phres%d2gval(ixsym(i1,i2),ipy)=phres%d2gval(ixsym(i1,i2),ipy)+&
                    add1
            enddo
         enddo secondd
! hm, when debugging here phres%dgval(1,*,1)=0 so ...
         add1=savedg(1,i1,1)
         sum=phres%dgval(1,i1,1)
         if(phres%dpqdy(i1).lt.1.0D-60) phres%dpqdy(i1)=zero
         do ipy=1,lprop-1
            do ider=1,3
! this calculates the proper ionic liquid model, not Q times
               phres%dgval(ider,i1,ipy)=&
                    savedg(ider,i1,ipy)+&
                    phres%sites(2)*phres%dgval(ider,i1,ipy)
! The contribution from the derivative of Q = \sum_i nu_i y_i, dQ/dy_i = nu_i
! G = G1 + Q G2 where
! G1 = \sum_i \sum_j y_i y_j G_ij + config.entropy
! G2 = y_va\sum_i y_i G_i + Q\sum_k y_k G_k
! Above were added:               dG/dy_i = dG1/dy_i + + Q dG2/dy_i 
! For cations we must add also    dG/dy_i = dG/dy_i + nu_i G2 
               if(i1.le.phlista(lokph)%nooffr(1)) then
! nooffr(1) is the number of constituents in first sublattice
                  phres%dgval(ider,i1,ipy)=phres%dgval(ider,i1,ipy)+&
                       phres%dpqdy(i1)*phres%gval(ider,ipy)
               endif
            enddo
         enddo
!     write(*,747)'3X suming: ',i1,savedg(1,i1,1)*rtg,phres%dgval(1,i1,1)*rtg,&
!              phres%dpqdy(i1),phres%gval(1,1)
!         write(*,747)'3Xx:',i1,add1,sum,phres%dgval(1,i1,1),phres%dpqdy(i1),&
!              phres%sites(2),savedg(1,i1,1)
!747      format(a,i2,6(1pe12.4))
      enddo firstd
!      write(*,*)'3X summed: ',savedg(1,1,1)*rtg,phres%dgval(1,1,1)*rtg
! Integral values: G = saveg + Q*phres%gval with T and P derivatives
490   continue
!      write(*,492)'3X ionsum: ',saveg(1,1),phres%gval(1,1),&
!           (saveg(1,1)+phres%gval(1,1))*rtg*phres%sites(2)
492   format(a,6(1pe12.4))
!      write(*,*)'3X Config G 7A: ',phres%gval(1,1)*rtg
      do ipy=1,lprop-1
         do ider=1,6
            phres%gval(ider,ipy)=saveg(ider,ipy)+&
                 phres%sites(2)*phres%gval(ider,ipy)
         enddo
      enddo
!      write(*,*)'3X Config G 7B: ',phres%gval(1,1)*rtg,saveg(1,1)*rtg
! strange bug which changes the results for a calculation with only C1
! if the ionic liquid has been non-suspended at some previous calculation ...
      saveg=zero
!      if(ocv()) write(*,*)'3X deallocated saveg 2: ',size(saveg)
! no need to set them zero if they will be deallocated??
!      savedg=zero
!      saved2g=zero
      deallocate(saveg)
      if(moded.gt.0) then
         deallocate(savedg)
         deallocate(saved2g)
      endif
!499   continue
   endif ionliqsum
!................................
! we have now finished calculate all parameters including those 
! properties that affect the Gibbs energy indirectly like Curie T etc
! For a model like the Flory-Huggins the partial molar volume must be
! calculated before any other part of G so all must be calculated again ...
! The label here just to indicate this, there is no explict jump here
500 continue
   if(btest(phmain%status2,CSADDG)) then
! we have an addition to G, at present just a constant /RT
      if(allocated(phmain%addg)) then
         xxx=phmain%addg(1)/ceq%rtn
      else
         write(*,*)'3X not allocated addg'
         xxx=zero
      endif
!      write(*,*)'Addition to G:',xxx
! a constant addition affects G and dG/dy and d2G/dy2
      phmain%gval(1,1)=phmain%gval(1,1)+xxx
! dgval( 1/dT/dP , i , property)
      do id=1,gz%nofc
         phmain%dgval(1,id,1)=phmain%dgval(1,id,1)+xxx
         do jd=id,gz%nofc
            phmain%d2gval(ixsym(id,jd),1)=phmain%d2gval(ixsym(id,jd),1)+xxx
         enddo
      enddo
   endif
   if(floryhuggins.lt.0) then
! The Flory-Huggins entropy require that we use the volume parameters
! These have now been calculated and can be used in a second loop through
! the other parameters
! find where the molar volumes are stored, phmain%listprop(1) is no props
!      write(*,507)'3X FH: ',nofc2,phmain%listprop(1),&
!           (phmain%listprop(ipy),ipy=2,phmain%listprop(1)-1)
!507   format(a,i5,i3,20i5)
      allocate(fhlista(gz%nofc))
      fhlista=0
      ll=1
      do ipy=2,phmain%listprop(1)
         if(phmain%listprop(ipy).gt.2000) then
! NOTE each element has a Flory-Huggins volume ... 2001, 2002 etc in any order
! fhlista(i) is the index to gval(*,ipy)
            fhlista(phmain%listprop(ipy)-2000)=ipy
         endif
      enddo
! we must save the Flory-Huggins volumes as they are used in next loop
      allocate(fhv(gz%nofc,6))
      allocate(dfhv(gz%nofc,3,gz%nofc))
      allocate(d2fhv(gz%nofc,nofc2))
      dfhv=zero
      d2fhv=zero
!      fhvsum=zero
      do qz=1,gz%nofc
         ipy=fhlista(qz)
! if ipy is zero then the volume is constant
         do ll=1,6
            if(ipy.gt.0) then
               fhv(qz,ll)=phmain%gval(ll,ipy)
            endif
         enddo
         do ll=1,gz%nofc
            if(ipy.gt..0) then
               dfhv(qz,1,ll)=phmain%dgval(1,ll,ipy)
               dfhv(qz,2,ll)=phmain%dgval(2,ll,ipy)
               dfhv(qz,3,ll)=phmain%dgval(3,ll,ipy)
            endif
         enddo
         do ll=1,nofc2
            if(ipy.gt.0) then
               d2fhv(qz,ll)=phmain%d2gval(ll,ipy)
            endif
         enddo
! this is the only non-zero value for elements with no Flory-Huggins para,eter
         if(ipy.eq.0) fhv(qz,1)=one
      enddo
! the Flory-Huggins parametr for each constituent is in the "fhv" arguments
! They may be updated inside this subroutine ...
      call config_entropy_floryhuggins(moded,gz%nofc,phmain,gz%tpv(1),&
           fhv,dfhv,d2fhv)
! set floryhuggins to 1 so the entropy is not calculated again but 
! the molar volumes that have been calculated here can be used
      floryhuggins=1
! we must calculate the other parameters again using the specific molar volumes
! not implemented yet ...
!      write(*,*)'3X Flory-Huggins model only config entropy, no goto 100'
      goto 100
   endif
!................................
! calculate additions like magnetic contributions etc and add to G
   addrec=>phlista(lokph)%additions
   additions: do while(associated(addrec))
! Note for phases with a disordered fraction set, gz%nofc is equal to
! the disordered number of fractions here 
      gz%nofc=phlista(lokph)%tnooffr
! moded is 0, 1 or 2 if derivatives should be calculated, phres is pointer
! to result arrays, lokadd is the addition record, listprop is needed to
! find where TC and BM are stored, gz%nofc are number of constituents
      call addition_selector(addrec,moded,phres,lokph,gz%nofc,ceq)
      if(gx%bmperr.ne.0) goto 1000
! NOTE that the addition record is not in the dynamic data struturce
! but the values calculated are returned added to phres
! There is a temporary storage of results for listing only.
      addrec=>addrec%nextadd
   enddo additions
! there are some special properties like mobilities and similar which
! have a conmponent or constituent index like MQ&<constituent>
!   ipy=typty/100+mod(typty,100)
!   if(ipy.gt.10) then
!      write(*,*)'3X Property ',typty,ipy
1000 continue
!   ipy=phlista(lokph)%linktocs(1)
!   write(*,*)'3X exit 1: ',lokph,ipy,ceq%phase_varres(ipy)%disfra%varreslink
!   ipy=phlista(lokph)%linktocs(2)
!   if(ipy.gt.0) &
!      write(*,*)'3X exit 2: ',lokph,ipy,ceq%phase_varres(ipy)%disfra%varreslink
   if(chkperm) then
! wait for checking for errors ....
!      write(*,*)'3X Press return'
!      read(*,297)ch1
!297   format(a)
   endif
! running out of memory??
   deallocate(dpyq)
   deallocate(d2pyq)
   deallocate(dvals)
   deallocate(d2vals)
   if(allocated(fhv)) then
      deallocate(fhv)
      deallocate(dfhv)
      deallocate(d2fhv)
   endif
!    write(*,1001)gx%bmperr,(phres%gval(i,1),i=1,4)
!    write(*,1002)(phres%dgval(1,i,1),i=1,3)
!    write(*,1003)(phres%d2gval(i,1),i=1,6)
1001 format('3X calcg g: ',i5,4(1PE15.7))
1002 format('3X calcg dg:  ',3(1PE15.7))
1003 format('3X calcg d2g: ',6(1PE11.3))
   return
 end subroutine calcg_internal

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine setendmemarr(lokph,ceq)
! stores the pointers to all ordered and disordered endmemners in arrays
! UNUSED ??
   implicit none
   integer lokph
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ll,nz,noemr
   TYPE(gtp_endmember), pointer :: emrec
   TYPE(gtp_fraction_set), pointer :: disfraset
   if(allocated(phlista(lokph)%oendmemarr)) then
      deallocate(phlista(lokph)%oendmemarr)
! allways allocate place for maximum endmembers (product of constituents)
      nz=1
      do ll=1,phlista(lokph)%noofsubl
         nz=nz*phlista(lokph)%nooffr(ll)
      enddo
      allocate(phlista(lokph)%oendmemarr(nz))
      noemr=0
      emrec=>phlista(lokph)%ordered
      do while(associated(emrec))
         noemr=noemr+1
         phlista(lokph)%oendmemarr(noemr)%p1=>emrec
         emrec=>emrec%nextem
      enddo
      phlista(lokph)%noemr=noemr
   endif
! same for disordered endmembers (if any)
! Data for this is stored in phase_varres record, same index as phlista !!!
   if(allocated(phlista(lokph)%dendmemarr)) then
      deallocate(phlista(lokph)%dendmemarr)
! allways allocate place for maximum endmembers (product of constituents)
      disfraset=>ceq%phase_varres(lokph)%disfra
      nz=1
      do ll=1,disfraset%ndd
         nz=nz*disfraset%nooffr(ll)
      enddo
      allocate(phlista(lokph)%dendmemarr(nz))
      noemr=0
      emrec=>phlista(lokph)%disordered
      do while(associated(emrec))
         noemr=noemr+1
         phlista(lokph)%dendmemarr(noemr)%P1=>emrec
         emrec=>emrec%nextem
      enddo
      phlista(lokph)%ndemr=noemr
   endif
1000 continue
   return
 end subroutine setendmemarr

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine tabder(iph,ics,ceq)
! tabulate derivatives of phase iph with current constitution and T and P
   implicit none
   integer iph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character name*24
   double precision kappa,napfu,t,p,rtg,g,v,s,h,u,f,cp,alpha
   integer tnk,lokph,nsl,lokres,lokcs,ll,ll2,kk1,kk2,kk3,kk4,loksp
!
   lokph=phases(iph)
   nsl=phlista(lokph)%noofsubl
! calculate G and derivatives, lokres returns index of phase_varres
   call calcg(iph,ics,2,lokres,ceq)
   if(gx%bmperr.ne.0) then
      goto 1000
   endif
! number of moles of atoms per formula unit
   napfu=ceq%phase_varres(lokres)%abnorm(1)
   T=ceq%tpval(1)
   P=ceq%tpval(2)
   rtg=globaldata%rgas*T
   lokcs=lokres
! returned values: G, G.T=-S, G.P=V, G.T.T=-Cp/T G.T.P=V*alpha, G.P.P=-V*kappa
! all divided by RT and per mole formula unit of phase
! G=H-TS, F=U-TS, H=U+PV, S=-G.T, V=G.P
! H=G+TS=G-T*G.T, U=H-PV=(G-T*G.T)-P*G.P, CP=-T*G.T.T
! alpha= 1/V*V.T = G.T.P/V, kappa = -1/V*V.P = -G.P.P/V
   G=rtg*ceq%phase_varres(lokcs)%gval(1,1)
!    write(*,5)'3X tabder 2: ',rtg,G
   S=-rtg*ceq%phase_varres(lokcs)%gval(2,1)
   V=rtg*ceq%phase_varres(lokcs)%gval(3,1)
   H=G+T*S
   U=H-P*V
   F=U-T*S
   CP=-T*rtg*ceq%phase_varres(lokcs)%gval(4,1)
   if(V.ne.zero) then
      alpha=rtg*ceq%phase_varres(lokcs)%gval(5,1)/V
      kappa=rtg*ceq%phase_varres(lokcs)%gval(6,1)/V
   else
      alpha=zero
      kappa=zero
   endif
   write(kou,100)napfu,rtg,T,P,G,G/rtg
100 format(/'Per mole FORMULA UNIT of the phase, ',1pe12.4,' atoms/F.U., RT=',&
         1pe15.7/&
         'at T= ',0pF8.2,' K and P= ',1PE13.6,' Pa',8x,'SI units',9x,'/RT'/ &
         'Gibbs energy J/mol ',28('.'),1Pe16.8,e16.7)
   write(kou,102)F,F/rtg,H,H/rtg,U,U/rtg,S,S/rtg,V,V/rtg,&
        CP,CP/rtg,alpha,alpha/rtg,kappa,kappa/rtg
102 format('Helmholtz energy J/mol ',24('.'),1PE16.8,e16.7 &
        /'Enthalpy J/mol ',32('.'),1PE16.8,e16.7 &
        /'Internal energy J/mol ',25('.'),1PE16.8,e16.7 &
        /'Entropy J/mol/K ',31('.'),1PE16.8,e16.7 &
        /'Volume m3 ',37('.'),1PE16.8,e16.7 &
        /'Heat capacity J/mol/K ',25('.'),1PE16.8,e16.7 &
        /'Thermal expansion 1/K ',25('.'),1PE16.8,e16.7 &
        /'Bulk modulus 1/Pa ',29('.'),1PE16.8,e16.7)
   tnk=phlista(lokph)%tnooffr
   ll=1
   kk1=0
   kk2=phlista(lokph)%nooffr(ll)
   dy1loop: do while(kk1.le.tnk)
      kk1=kk1+1
      if(kk1.gt.kk2) then
!          write(*,11)'3X tabder 2: ',kk1,kk2,ll,tnk,nsl
!11 format(a,10i3)
         ll=ll+1
         if(ll.gt.nsl) exit
         kk2=kk2+phlista(lokph)%nooffr(ll)
      endif
      if(phlista(lokph)%nooffr(ll).eq.1) then
!          write(*,*)'3X tabder 1: ',kk1,kk2,ll,tnk
         ll=ll+1
         if(ll.gt.nsl) exit
         kk2=kk2+phlista(lokph)%nooffr(ll)
         cycle
      endif
      loksp=phlista(lokph)%constitlist(kk1)
      name=splista(loksp)%symbol
      write(kou,110)name(1:len_trim(name)),ll
110 format('First partial derivative with respect to ',a,&
        ' in sublattice ',i2,' of')
      write(kou,120)rtg*ceq%phase_varres(lokcs)%dgval(1,kk1,1),&
           ceq%phase_varres(lokcs)%dgval(1,kk1,1),&
           rtg*(ceq%phase_varres(lokcs)%dgval(1,kk1,1)-&
           T*ceq%phase_varres(lokcs)%dgval(2,kk1,1)),&
           ceq%phase_varres(lokcs)%dgval(1,kk1,1)-&
           T*ceq%phase_varres(lokcs)%dgval(2,kk1,1),&
           rtg*ceq%phase_varres(lokcs)%dgval(2,kk1,1),&
           ceq%phase_varres(lokcs)%dgval(2,kk1,1),&
           rtg*ceq%phase_varres(lokcs)%dgval(3,kk1,1),&
           ceq%phase_varres(lokcs)%dgval(3,kk1,1)
120    format(5x,'G ',40('.'),1PE16.8,e16.7, &
           /5x,'H ',40('.'),1PE16.8,e16.7, &
           /5x,'G.T ',38('.'),1PE16.8,e16.7, &
           /5x,'G.P ',38('.'),1PE16.8,e16.7)
      kk3=kk1
      kk4=kk2
      ll2=ll
      write(kou,150)
150 format(5x,'Second partial derivative of Gibbs energy with respect to also')
      dy2loop: do while(kk3.le.tnk)
         if(phlista(lokph)%nooffr(ll2).gt.1) then
!            write(kou,160)name(1:len_trim(name)),ll2, &
            write(kou,160)name,ll2, &
                 rtg*ceq%phase_varres(lokcs)%d2gval(ixsym(kk1,kk3),1),&
                 ceq%phase_varres(lokcs)%d2gval(ixsym(kk1,kk3),1)
160          format(10x,a,'   in ',i2,5('.'),1PE16.8,e16.7)
         endif
         kk3=kk3+1
         if(kk3.le.tnk) then
            loksp=phlista(lokph)%constitlist(kk3)
            name=splista(loksp)%symbol
         endif
         if(kk3.gt.kk4) then
            ll2=ll2+1
            if(ll2.gt.nsl) exit
            kk4=kk4+phlista(lokph)%nooffr(ll2)
         endif
      enddo dy2loop
!       write(*,*)'3X tabder 7A: ',kk1,kk2
   enddo dy1loop
900 continue
!    write(*,*)'3X tabder 7B: ',kk2
!    write(*,*)'3X tabder: ',rtg,rtg*phase_varres(lokcs)%gval(1,1)
1000 continue
   return
 end subroutine tabder

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine cgint(lokph,lokpty,moded,vals,dvals,d2vals,gz,ceq)
! calculates an excess parameter that can be composition dependent
! gz%yfrem are the site fractions in the end member record
! gz%yfrint are the site fractions in the interaction record(s)
! lokpty is the property index, lokph is the phase record
! moded=0 means only G, =1 G and dG/dy, =2 all
   implicit none
   integer moded,lokph
   TYPE(gtp_property), pointer :: lokpty
   TYPE(gtp_parcalc) :: gz
   double precision vals(6),dvals(3,gz%nofc)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! temporary data like gz%intlevel, gz%nofc etc
   double precision d2vals(gz%nofc*(gz%nofc+1)/2),valtp(6)
   double precision vv(0:2),fvv(0:2)
   integer lfun,jdeg,jint,qz,ivax,icat
   double precision rtg,dx0,dx,dx1,dx2,ct,fvs,dvax0,dvax1,dvax2,yionva
   double precision ycat0,dcat1,dcat2,dyvan1,dyvan2
   double precision, parameter :: onethird=one/3.0D0,two=2.0D0
   logical ionicliq,iliqva,iliqneut
! zeroing 5 iq, and vals, dvals and d2vals
!   write(*,*)'3X cgint 1:',gz%iq(1),gz%iq(2),gz%iq(3)
! why zero qz%iq, it has been set before calling ...
   gz%iq=0
   vals=zero
   dvals=zero
   d2vals=zero
   rtg=gz%rgast
!   write(*,*)'3X in cgint',lokph
   if(lokpty%degree.eq.0) then
!----------------------------------------------------------------------
! no composition dependence
      lfun=lokpty%degreelink(0)
      call eval_tpfun(lfun,gz%tpv,vals,ceq%eq_tpres)
      if(gx%bmperr.ne.0) goto 1000
      if(lokpty%proptype.eq.1) then
         vals=vals/rtg
      endif
      goto 1000
   endif
!----------------------------------------------------------------------
! for composition dependent param set default variables for ionic liquid
   ionicliq=.FALSE.
   iliqva=.FALSE.
   iliqneut=.FALSE.
   yionva=zero
   if(btest(phlista(lokph)%status1,PHIONLIQ)) then
! prepare for ionic liquid interactions
!      write(*,17)'3X RK: ',phlista(lokph)%i2slx(1),gz%endcon(gz%intlat(1))
17    format(a,10i4)
!      write(*,*)'3X ionicliq set true'
!      write(*,17)'3X Const in subl: ',gz%intlat(1),gz%endcon(gz%intlat(1)),&
!           gz%endcon(2),phlista(lokph)%i2slx(1),gz%intlevel
      ionicliq=.TRUE.
      if(gz%endcon(2).eq.phlista(lokph)%i2slx(1)) then
! VA endmember in the 2nd sublattice, this is the complicated case
         yionva=gz%yfrem(2)
         ivax=phlista(lokph)%i2slx(1)
!         write(*,64)'3X iliq with Va: ',ivax,yionva
64       format(a,i3,6(1pe12.4))
         if(gz%intlat(1).eq.1) then
! interaction in sublattice 1 between two cations same as substituional L_A,B
! with each cation fraction multiplied with vacancy 
! Also set TRUE for reciprocal interactions (gz%intlevel=2)
            iliqva=.TRUE.
         else
! interaction in sublattice 2 between Va and neutral (i.e. cation and neutral)
! same as substitutional L_A,B with cation fraction multiplied with vacancy
! Hm, I am not sure interactions are ordered so all interactions in first
! sublattice comes before any in second sublattice ??
            iliqneut=.TRUE.
         endif
!      else
! constituent in second sublattice is not vacancy, no particular action ??
!         write(*,17)'3X 2nd sublattice constituent not Va: ',gz%endcon(2)
      endif
   endif
   intlev: if(gz%intlevel.eq.1) then
!----------------------------------------------------------------------
! plain binary Redlich Kister. gz%endcon can be wildcard, i.e. negative
! but for the moment give error message in that case
! A binary wildcard excess parameter means y_A ( 1 - y_A) * L_A*
! most naturally gz%intcon(1) would be negative
      gz%iq(1)=gz%endcon(gz%intlat(1))
      gz%iq(2)=gz%intcon(1)
      if(gz%iq(1).lt.0 .or. gz%iq(2).lt.0) then
! composition dependent wildcard interaction not implemented
! y(1-y) ( L0 + (2y-1) L1 + (2y-1)**2 L2 + ....) ??
         gx%bmperr=4031; goto 1000
      endif
! endmember fraction minus interaction fraction
      dx0=gz%yfrem(gz%intlat(1))-gz%yfrint(1)
! ycat is one unless ionic liquid with vacancy-neutral interaction
      ycat0=one
      if(ionicliq) then
         if(iliqva) then
! interaction between cations with vacancy on second sublattice
! NOTE intraction fraction alreay multiplied with yionva before calling cgint
            dvax0=gz%yfrem(gz%intlat(1))-gz%yfrint(1)/yionva
!            dvax0=dx0
            dx0=yionva*dvax0
!            write(*,65)'3X Va on 2nd: ',gz%iq(2),gz%intlat(1),dvax0,dx0,&
!                 gz%yfrem(gz%intlat(1)),gz%yfrint(1)
65          format(a,2i3,6(1pe12.4))
         elseif(iliqneut) then
! interaction between vacancy and neutral in second sublattice
! we must know the cation
            icat=gz%endcon(1)
            ycat0=gz%yfrem(1)
! the fraction difference is between (y_cation * y_Va - y_neutral)
            dx0=gz%yfrem(1)*yionva-gz%yfrint(1)
            dvax0=ycat0
!            write(*,*)'3X dx0: ',dx0,ycat0,yionva
         endif
      endif
      vals=zero
      dx=one
      dx1=zero
      dx2=zero
      dvax1=zero
      dvax2=zero
      dyvan1=one
      dyvan2=one
!      write(*,*)'3X cgint 2:',gz%iq(1),gz%iq(2),gz%iq(3),icat
!      write(*,*)'3X c1bug: ',ionicliq,iliqva,iliqneut
! special for ionic liquid: when two cation interacts with Va in second
! sublattice the vacancy fraction is raised by power 2
      RK: do jdeg=0,lokpty%degree
         lfun=lokpty%degreelink(jdeg)
         call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
         if(gx%bmperr.ne.0) goto 1000
         if(lokpty%proptype.eq.1) then
! property type 1 is G and should be normalized by RT
            valtp=valtp/rtg
         endif
! vals and valtp are arrays with 6 elements: G, G.T, G.P, G.T.T ...
         vals=vals+dx*valtp
!         write(*,11)'3X dx: ',gz%iq(1),gz%iq(2),jdeg,vals(1),dx,valtp(1)
11       format(a,3i2,6(1pe11.3))
! no composition derivative.  if moded=0 only G, =1 G+G.Y, =2 all
         noder5: if(moded.gt.0) then
! first derivatives, qz=1: dG/dyA dG/dyB; qz=2: d2G/dTdy; qz=3: d3G/dPdy
! for iliqneut there should not be same -dx1 ... gz%iq(2) is neutral
            do qz=1,3
! For interactions between Va and neutral in ionic liguid a power of yionva
! is required for the cation derivative as we have (y_cation*yionva-y_neutral)
! In all other cases dyvan1=unity
               dvals(qz,gz%iq(1))=dvals(qz,gz%iq(1))+ycat0*dx1*valtp(qz)
               dvals(qz,gz%iq(2))=dvals(qz,gz%iq(2))-dx1*valtp(qz)
! The handling of ionic liquid parameter derivatives can be simplified ...
               if(iliqva) then
! derivative with respect to vacancy fraction for (yc1-yc2)*yva: yc1-yc2
                  dvals(qz,ivax)=dvals(qz,ivax)+dvax1*valtp(qz)
!                  if(qz.eq.1) write(*,11)'3X iliqva: ',0,0,ivax,dvax1
              elseif(iliqneut) then
! derivative with respect to cation (yc1*yva-yn):
! multiply with a power of y_Va
                  dvals(qz,icat)=dvals(qz,icat)+yionva*dx1*valtp(qz)
!                write(*,19)'3X mess:',qz,gz%iq(1),gz%iq(2),icat,yionva,ycat0,&
!                       valtp(qz),dx1
19                format(a,4i3,6(1pe12.4))
               endif
            enddo
!            write(*,11)'3X dx1:',gz%iq(1),gz%iq(2),jdeg,dvals(1,gz%iq(1)),&
!                 dvals(1,gz%iq(2)),dx1,valtp(1)
! second derivatives, d2G/dyAdyA d2G/dyAdyB d2G/dyBdyB
            if(moded.gt.1) then
               d2vals(ixsym(gz%iq(1),gz%iq(1)))=&
                    d2vals(ixsym(gz%iq(1),gz%iq(1)))+dx2*valtp(1)
               d2vals(ixsym(gz%iq(1),gz%iq(2)))=&
                    d2vals(ixsym(gz%iq(1),gz%iq(2)))-dx2*valtp(1)
               d2vals(ixsym(gz%iq(2),gz%iq(2)))=&
                    d2vals(ixsym(gz%iq(2),gz%iq(2)))+dx2*valtp(1)
!               if(iliqva) then
! UNFINISHED d2G/dyvdyv d2G/dyvdyA d2G/dyvdyB interactions two cations
!                  d2vals(ixsym(ivax,ivax))=&
!                       d2vals(ixsym(ivax,ivax))+dvax2*valtp(1)
!               elseif(iliqneut) then
! UNFINISHED also for interactions Va-neutral
!                  continue
!               endif
            endif
         endif noder5
! next power of dx
         if(iliqva) then
! interaction between two cations, dx0=y_va*(y_c1 - y_c2)
! NO CHANGE HERE WHEN FIXING ERROR FOR Va-Neutal interaction ...
            dx2=(jdeg+1)*dx1
            dvax2=(jdeg+1)*dvax1
            if(jdeg.eq.0) then
               dx1=yionva
               dvax1=dvax0
            else
               dx1=(jdeg+1)*dx1*dx0
               dvax1=(jdeg+1)*dvax1*dx0
            endif
            dx=dx*dx0
!            write(*,23)'3X iliqvb: ',jdeg,dx,dx1,dx2,dvax0,dvax1,dvax2
23          format(a,i2,6(1pe12.4))
         elseif(iliqneut) then
! interaction between Va and neutral a bit more complicated ... NOT TESTED
! NOTE 2nd derivatives ignored ...
            dx2=(jdeg+1)*dx1
            dvax2=dvax1
            if(jdeg.eq.0) then
               dx1=one
               dvax1=dvax0
            else
               dx1=(jdeg+1)*dx1*dx0
               dvax1=(jdeg+1)*dvax0*dx1
            endif
            dx=dx*dx0
         else
! normal CEF model
            dx2=(jdeg+1)*dx1
            dx1=(jdeg+1)*dx
            dx=dx*dx0
         endif
      enddo RK
   elseif(gz%intlevel.eq.2) then !intlev
!----------------------------------------------------------------------
! important to set ivax=0 here as tested below if not zero
      ivax=0
      if(ionicliq) then
!         write(*,*)'3X Comp.dep ternary ionic liquid parameter: ',iliqva
         if(iliqva) then
! the endmember constituent in second sublattice is Va, no anions!!
            if(gz%intlat(1).eq.1 .and. gz%intlat(2).eq.1) then
! we have 3 cations interacting in first sublattice and Va in second
! require treatment of extra vacancy fraction
               write(*,*)'3X 3 interacting cations in liquid not implemented'
               gx%bmperr=4343; goto 1000
            elseif(gz%intlat(1).eq.1 .and. gz%intlat(2).eq.2) then
               ivax=gz%endcon(2)
            endif
         elseif(gz%intcon(2).eq.phlista(lokph)%i2slx(1)) then
! reciprocal interaction between anion and vacancy is the second sublattice
! with 2 cations in first sublattice
!            write(*,*)'3X reciprocal with 2 cations and anion and Va'
            ivax=gz%intcon(2)
            yionva=gz%yfrint(2)
         else
            write(*,*)'3X: ionic liquid model parameter not implemented'
            gx%bmperr=4342; goto 1000
         endif
! other ternary parameters in ionic liquid OK, no extra vacancy fraction
      endif
!................................................................
! ternary composition dependent interaction
      ternary: if(gz%intlat(1).eq.gz%intlat(2)) then
! Ternary composition dependent interaction in same sublattice, Hillert form.
! The idea is that the sum of vv is always unity even in higher order systems
! whereas the sum of the constituent frations are not
! If wildcard then any of the gz%iq would be negative, not allowed
         gz%iq(1)=gz%endcon(gz%intlat(1))
         gz%iq(2)=gz%intcon(1)
         gz%iq(3)=gz%intcon(2)
         if(gz%iq(1).lt.0 .or. gz%iq(2).lt.0 .or. gz%iq(3).lt.0) then
            gx%bmperr=4031; goto 1000
         endif
         vv(0)=gz%yfrem(gz%intlat(1))
         vv(1)=gz%yfrint(1)
         vv(2)=gz%yfrint(2)
         ct=(one-vv(0)-vv(1)-vv(2))*onethird
         vv=vv+ct
         fvv(0)=two*onethird
         fvv(1)=-onethird
         fvv(2)=-onethird
         terloop: do jint=0,2
! calculate parameter
            lfun=lokpty%degreelink(jint)
            call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
            if(lokpty%proptype.eq.1) then
               valtp=valtp/rtg
            endif
! function value
            vals=vals+vv(jint)*valtp
            noder6: if(moded.gt.0) then
! first derivatives
               do qz=1,3
                  dvals(qz,gz%iq(1))=dvals(qz,gz%iq(1))+fvv(0)*valtp(qz)
                  dvals(qz,gz%iq(2))=dvals(qz,gz%iq(2))+fvv(1)*valtp(qz)
                  dvals(qz,gz%iq(3))=dvals(qz,gz%iq(3))+fvv(2)*valtp(qz)
               enddo
! there is no contribution to the second derivatives from this interaction
            endif noder6
            fvs=fvv(2)
            fvv(2)=fvv(1)
            fvv(1)=fvv(0)
            fvv(0)=fvs
         enddo terloop
      else
!.........................................................
! composition dependent reciprocal interactions here only degree 1 and 2
         if(lokpty%degree.gt.2) then
            write(*,*)'3X Composition dependent reciprocal degree max 2'
            gx%bmperr=4078; goto 1000
         else
!            write(*,32)lokph,lokpty%degree,gz%intlat(1),gz%intlat(2),&
!                 gz%iq(1),gz%iq(2),gz%iq(3),gz%iq(4)
32          format('3X Comp.dep. rec. param: ',i3,2x,i1,2x,2i2,4i5)
         endif
! Note the composition dependence is defined that 
! L = y'_Ay'_By"_Cy"_D (0L + (y"_C-y"_D)*1L + (y'_A-y'_B)*2L)
! it is a bit strange that 2nd sublattice is 1L ... but that is the definition
         gz%iq(1)=gz%endcon(1)
         gz%iq(2)=gz%intcon(1)
         gz%iq(3)=gz%endcon(2)
         gz%iq(4)=gz%intcon(2)
! degree 0 not composition dependent, vals multiplied with pyq after return
         lfun=lokpty%degreelink(0)
         if(lfun.gt.0) then
            call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
            if(gx%bmperr.ne.0) goto 1000
            if(lokpty%proptype.eq.1) then
               valtp=valtp/rtg
            endif
            vals=vals+valtp
         endif
! lokpty%degree must be 1 or 2 otherwise we would not be here
!         write(*,17)'3X composition dependent reciprocal',ivax
         lfun=lokpty%degreelink(1)
         recip1: if(lfun.gt.0) then
! degree 2 can be empty, otherwise multiplied with gz%iq(3)-gz%iq(4)
! no problem with ionic liquid except there may be values in dvals
            call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
            if(gx%bmperr.ne.0) goto 1000
!            write(*,62)'3X rp: ',valtp(1),valtp(2)
62          format(a,6(1pe12.4))
            if(lokpty%proptype.eq.1) then
               valtp=valtp/rtg
            endif
            vals=vals+(gz%yfrem(gz%intlat(2))-gz%yfrint(2))*valtp
! dvals(1,const) is the 1st derivative of the fun wrt const
! dvals(2,const) is the 2nd derivative of the fun wrt const and T
! dvals(3,const) is the 2nd derivative of the fun wrt const and P
! one dvals(*,ivax) could have been assigned a value above (for ionic liquid)
            do qz=1,3
!      write(*,63)'3X dvals: ',qz,gz%iq(3),dvals(qz,gz%iq(3)),&
!           dvals(qz,gz%iq(4)),valtp(qz)
               dvals(qz,gz%iq(3))=dvals(qz,gz%iq(3))+valtp(qz)
               dvals(qz,gz%iq(4))=dvals(qz,gz%iq(4))-valtp(qz)
!      write(*,63)'3X dvals: ',qz,gz%iq(4),dvals(qz,gz%iq(3)),dvals(qz,gz%iq(4))
            enddo
63          format(a,2i3,6(1pe12.4))
         endif recip1
! degree 2 can be empty, otherwise multiplied with y(gz%iq(1))-y(gz%iq(2))
         recip2: if(lokpty%degree.gt.1) then 
            lfun=lokpty%degreelink(2)
            if(lfun.gt.0) then
               call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
               if(gx%bmperr.ne.0) goto 1000
               if(lokpty%proptype.eq.1) then
                  valtp=valtp/rtg
               endif
               if(ivax.gt.0) then
!                  write(*,67)ivax,gz%iq(1),gz%iq(2),gz%iq(3),gz%iq(4),yionva
!67                format('3X ion liq recip: ',i3,2x,4i3,1pe12.4)
! interaction in ionic liquid with vacancy as one constituent in 2nd subl.
                  vals=vals+yionva*(gz%yfrem(gz%intlat(1))-gz%yfrint(1))*valtp
                  do qz=1,3
                     dvals(qz,gz%iq(1))=+yionva*valtp(qz)
                     dvals(qz,gz%iq(2))=-yionva*valtp(qz)
                  enddo
! we have to take into account extra derivatives wrt vacancies if vacancy
! is a constituent in second sublattice
                  do qz=1,3
                     dvals(qz,ivax)=&
                          (gz%yfrem(gz%intlat(1))-gz%yfrint(1))*valtp(qz)
                  enddo
               else
! not ionic liquid .... puuuh
                  vals=vals+(gz%yfrem(gz%intlat(1))-gz%yfrint(1))*valtp
                  do qz=1,3
                     dvals(qz,gz%iq(1))=+valtp(qz)
                     dvals(qz,gz%iq(2))=-valtp(qz)
                  enddo
               endif
            endif
         endif recip2
      endif ternary
!----------------------------------------------------------------------
   elseif(gz%intlevel.ge.3) then !intlev
! higher interaction levels have no composition dependence
      write(*,999)
999   format('Composition dependence for parameters with >2 interacting ',&
           'constituents'/'not implemented!')
      gx%bmperr=4078; goto 1000
   endif intlev
!----------------------------------------------------------------------
! finished finally ....
1000 continue
   return
 end subroutine cgint

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine config_entropy(moded,nsl,nkl,phvar,tval)
! calculates configurational entropy/R for phase lokph
   implicit none
   integer moded,nsl
   integer, dimension(nsl) :: nkl
   TYPE(gtp_phase_varres), pointer :: phvar
!\end{verbatim}
   integer ll,kk,kall,nk,jl
   double precision tval,ss,yfra,ylog
   ll=0
   kall=0
   sublatticeloop: do while (ll.lt.nsl)
      ll=ll+1
      nk=nkl(ll)
      kk=0
      ss=zero
      fractionloop: do while (kk.lt.nk)
         kk=kk+1
         kall=kall+1
         if(nk.eq.1) cycle sublatticeloop
         yfra=phvar%yfr(kall)
         if(yfra.lt.bmpymin) yfra=bmpymin
         if(yfra.gt.one) yfra=one
         ylog=log(yfra)
! gval(1:6,1) are G and derivator wrt T and P
! dgval(1,1:N,1) are derivatives of G wrt fraction 1:N
! dgval(2,1:N,1) are derivatives of G wrt fraction 1:N and T
! dgval(3,1:N,1) are derivatives of G wrt fraction 1:N and P
! d2dval(ixsym(N*(N+1)/2),1) are derivatives of G wrt fractions N and M
! this is a symmetric matrix and index givem by ixsym(M,N)
         ss=ss+yfra*ylog
         if(moded.gt.0) then
            phvar%dgval(1,kall,1)=phvar%sites(ll)*(one+ylog)
            phvar%d2gval(ixsym(kall,kall),1)=phvar%sites(ll)/yfra
         endif
      enddo fractionloop
      phvar%gval(1,1)=phvar%gval(1,1)+phvar%sites(ll)*ss
   enddo sublatticeloop
! set temperature derivative of G and dG/dy
   phvar%gval(2,1)=phvar%gval(1,1)/tval
   if(moded.gt.0) then
      do jl=1,kall
         phvar%dgval(2,jl,1)=phvar%dgval(1,jl,1)/tval
      enddo
   endif
1000 continue
   return
 end subroutine config_entropy

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine config_entropy_floryhuggins(moded,nofc,phvar,tval,fhv,dfhv,d2fhv)
! calculates configurational entropy/R for a phase with Flory-Huggins model
! moded=0 only G; =1 G and dG/dy; =2; G, dG/dy and d2G/dy2
! nofc number of constituents, phvar phase_varres record
   implicit none
   integer moded,nofc
   TYPE(gtp_phase_varres), pointer :: phvar
! fvh(1,1) is FH volume for constituent 1 etc.
   double precision tval,fhv(nofc,*),dfhv(nofc,3,*),d2fhv(nofc,*)
!\end{verbatim}
   integer kall,nofc2,k1,k2
   double precision ss,sy,st,sp,yfra1,ylog,sumq
   double precision, allocatable :: pfhv(:,:),dpfhv(:,:,:),d2pfhv(:,:)
   double precision, allocatable :: yfra(:),qfra(:),sumsy(:,:)
!
   nofc2=nofc*(nofc+1)/2
!   write(*,1)'3X Config entropy FH model: ',nofc,(fhv(kall,1),kall=1,nofc)
1  format(a,i3,5F8.2)
107 format(a,6(1pe12.4))
108 format(a,i2,6(1pe12.4))
   allocate(yfra(nofc))
   allocate(qfra(nofc))
   allocate(sumsy(nofc,nofc))
! temporary arrays, maybe all not needed??
   allocate(pfhv(nofc,6))
   allocate(dpfhv(nofc,3,nofc))
   allocate(d2pfhv(nofc,nofc2))
   dpfhv=zero
   d2pfhv=zero
!
! sum the FH volumes for current composition and use as normallizing
   sumq=zero
   sumsy=zero
   do kall=1,nofc
      yfra1=phvar%yfr(kall)
      if(yfra1.lt.bmpymin) yfra1=bmpymin
      if(yfra1.gt.one) yfra1=one
!v      sumq=sumq+fhv(kall,1)*yfra1
      sumq=sumq+yfra1
      yfra(kall)=yfra1
!      do k1=1,nofc
!         if(k1.eq.kall) then
! 1st DERIVATIVES of q_i = p_i/\sum_j p_j
! fhv(i,1) is FH volume for const i, fhv(i,2) is T deriv, fhv(i,3) is P der
! dfhv(i,1,j) derivative of FH volume for i wrt const j
! dfhv(i,2,j) 2nd derivative of FH volume for i wrt const j and T
! dfhv(i,3,j) 2nd derivative of FH volume for i wrt const j and P
! UNFINISHED ?? sumsy including T and P derivatives ??
!            dpfhv(kall,1,k1)=dfhv(kall,1,k1)*yfra(kall)+fhv(kall,1)
!            dpfhv(kall,2,k1)=dfhv(kall,2,k1)*yfra(kall)+fhv(kall,2)
!            dpfhv(kall,3,k1)=dfhv(kall,3,k1)*yfra(kall)+fhv(kall,3)
!            sumsy(1,k1)=sumsy(1,k1)+dfhv(kall,1,k1)*yfra1+fhv(kall,1)
!            write(*,106)'3X sum1: ',kall,k1,sumsy(1,k1)
!         else
!            dpfhv(kall,1,k1)=dfhv(kall,1,k1)*yfra(kall)
!            dpfhv(kall,2,k1)=dfhv(kall,2,k1)*yfra(kall)
!            dpfhv(kall,3,k1)=dfhv(kall,3,k1)*yfra(kall)
!            sumsy(1,k1)=sumsy(1,k1)+dfhv(kall,1,k1)*yfra1
!            write(*,106)'3X sum2: ',kall,k1,sumsy(1,k1)
!         endif
!      enddo
   enddo
106 format(a,2i3,1pe12.4)
! we should extract fhv for each constituent from the species record ...
!   write(*,117)'3X sumq, vi: ',sumq,(fhv(kall,1)*yfra(kall)/sumq,kall=1,nofc)
117 format(a,6(1pe12.4))
!-----------------------------------------
! fhv(i,1) is FH volume for const i, fhv(i,2) is T deriv, fhv(i,3) is P deriv
! dfhv(i,1,j) derivative of FH volume for i wrt const j
! dfhv(i,2,j) 2nd derivative of FH volume for i wrt const j and T
! dfhv(i,3,j) 2nd derivative of FH volume for i wrt const j and P
! d2fhv(i,ixsym(j,k)) 2nd derivative of FH volume for i wrt const j and k
! gval(1:6,1) are G and derivator wrt T and P
! dgval(1,1:N,1) are derivatives of G wrt fraction 1:N
! dgval(2,1:N,1) are derivatives of G wrt fraction 1:N and T
! dgval(3,1:N,1) are derivatives of G wrt fraction 1:N and P
! d2dval(ixsym(N*(N+1)/2),1) are derivatives of G wrt fractions N and M
! this is a symmetric matrix and index givem by ixsym(M,N)
! ========== IMPORTANT: We have not implemented composition dependence in the
! Gibbs energy calculations !!!  Test that is not used!
   do kall=1,nofc
      ss=dfhv(kall,1,1)
      do k1=2,nofc
         if(abs(dfhv(kall,1,k1)-ss).gt.1.0D-8) then
            write(*,77)kall,k1,dfhv(kall,1,k1)
77          format(' *** Warning, Flory-Huggins model implemented',&
                 ' for constant FHV only',2i3,1pe12.4)
         endif
      enddo
   enddo
! =========================================================================
! Calculate the confurational entropy
   ss=zero
   fractionloop: do kall=1,nofc
! We use the already calculated partial molar volumes v_i = fhv_i * y_i
! This means we cannot calculate this before calculating all parameters once!!
! so this routine must be called after a first calculation of fhv, dfhv etc
! and then we must calculate all parameters again ... as they depend of v_i
! this is ln(n_i/(\sum_j n_j)), 0<yfra(kall)<1
!v      ylog=log(fhv(kall,1)*yfra(kall)/sumq)
!      ylog=log(yfra(kall)/sumq)    .... sumq=1.0!!
      ylog=log(yfra(kall))
      if(moded.gt.0) then
! UNFINISHED derivatives wrt T, P
         loopk1: do k1=1,nofc
            loopk2: do k2=k1,nofc
! UNFINISHED all second derivatives ignored, just set as for ideal 1/y
               if(kall.eq.k1 .and. k1.eq.k2) then
                  phvar%d2gval(ixsym(k1,k2),1)=one/yfra(kall)
               endif
            enddo loopk2
! This is df_i/dz, eq. 11 in FH documentation. Note dpfhv and sumsy set
! above to include the extra term if kall=k1
!            sy=dpfhv(kall,1,k1)*(ylog+one)-fhv(kall,1)/sumq*sumsy(1,k1)
! UNFINISHED ... IGNORE THE COMPOSITION DEPENENCE OF fhv_i ...
         enddo loopk1
! 
!         phvar%dgval(1,kall,1)=fhv(kall,1)/sumq*(ylog+one)
!         phvar%dgval(1,kall,1)=(ylog+one)/sumq
! each species has now a flory-Huggins segment value in fhv(kall,1)
         phvar%dgval(1,kall,1)=(ylog+one)/fhv(kall,1)
         phvar%dgval(2,kall,1)=phvar%dgval(2,kall,1)/tval
      endif
!      ss=ss+(fhv(kall,1)/sumq)*yfra(kall)*ylog
!v      ss=ss+yfra(kall)*ylog
      ss=ss+yfra(kall)/fhv(kall,1)*ylog
!      write(*,300)'3X ss: ',ss,yfra(kall),fhv(kall,1),fhv(kall,1)*yfra/sumq,&
!           ylog
300   format(a,6(1pe12.4))
   enddo fractionloop
! each species may now have a Flory Huggins segment number ....
!   ss=ss/sumq
   ss=ss
! The integral entropy and its T and P derivatives
! UNFINISHED should include T and P derivatives of fhv ....
!   phvar%gval(1,1)=phvar%gval(1,1)+ss
   phvar%gval(1,1)=ss
! UNFINISHED add T derivates of dfhv .... and any P derivatives
   phvar%gval(2,1)=phvar%gval(1,1)/tval
1000 continue
   return
 end subroutine config_entropy_floryhuggins

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine config_entropy_i2sl(moded,nsl,nkl,phvar,i2slx,tval)
! calculates configurational entropy/R for ionic liquid model
! Always 2 sublattices, the sites depend on composition
! P = \sum_j (-v_j) y_j + Q y_Va
! Q = \sum_i v_i y_i
! where v is the charge on the ions. P and Q calculated by set_constitution
   implicit none
   integer moded,nsl,i2slx(2)
   integer, dimension(nsl) :: nkl
   TYPE(gtp_phase_varres), pointer :: phvar
!\end{verbatim}
   integer ll,kk,kall,nk,j1,j2
   double precision tval,ss,yfra,ylog,yva,spart(2)
   ll=0
   kall=0
   spart=zero
   yva=zero
   sublatticeloop: do while (ll.lt.nsl)
      ll=ll+1
      nk=nkl(ll)
      kk=0
      ss=zero
      fractionloop: do while (kk.lt.nk)
         kk=kk+1
         kall=kall+1
! no cycle as we may need values of spart and yva ...
!         if(nk.eq.1) cycle sublatticeloop
         yfra=phvar%yfr(kall)
         if(yfra.lt.bmpymin) yfra=bmpymin
         if(yfra.gt.one) yfra=one
! save current value of vacancy fraction
         if(kall.eq.i2slx(1)) yva=yfra
!         write(*,2)'3X yva: ',kall,i2slx(1),yva,yfra
!2        format(a,2i3,6(1pe12.4))
         ylog=log(yfra)
! gval(1:6,1) are G and derivator wrt T and P
! dgval(1,1:N,1) are derivatives of G wrt fraction 1:N
! dgval(2,1:N,1) are derivatives of G wrt fraction 1:N and T
! dgval(3,1:N,1) are derivatives of G wrt fraction 1:N and P
! d2dval(ixsym(N*(N+1)/2),1) are derivatives of G wrt fractions N and M
! this is a symmetric matrix and index givem by ixsym(M,N)
         ss=ss+yfra*ylog
         if(moded.gt.0) then
            phvar%dgval(1,kall,1)=phvar%sites(ll)*(one+ylog)
            phvar%d2gval(ixsym(kall,kall),1)=phvar%sites(ll)/yfra
         endif
      enddo fractionloop
      phvar%gval(1,1)=phvar%gval(1,1)+phvar%sites(ll)*ss
      if(ll.eq.1) then
         spart(1)=ss
      else
         spart(2)=ss
      endif
   enddo sublatticeloop
   if(moded.eq.0) goto 900
! convergence problem with ionic liquid, skip contribution to 2nd derivatuves
!   localmoded=moded
!   if(moded.eq.2) localmoded=1
!   write(*,*)'3X ionic config_entropy: ',i2slx,kall
! additional derivatives as sublattice sites depend on composition
! -------------------------- derivatives of config entropy
! S = P*S1 + Q*S2
! S1 = \sum_i y_i*ln(y_i)
! S2 = \sum_j y_j*ln(y_j)+y_Va*ln(y_Va)+\sum_k y_k*ln(Y_k))
! P = \sum_j (-v_j)*y_j + Q*y_Va
! Q = \sum_i v_i*y_i
! term within [...] already calculated as part of normal config.entropy
! dS/dy_i        = +v_i*S2 + v_i*y_Va*S1 + [P*(1+ln(y_i)]   ..cation   OK
! dS/dy_j        = -v_j*S1 +               [Q*(1+ln(y_j))]  ..anion    OK
! dS/dy_Va       = Q*S1 +                  [Q*(1+ln(y_Va))] ..Va       OK
! dS/dy_k        =                         [Q*(1+ln(y_k)]   ..neutral  OK
! d2S/dy_i1dy_i2 = v_i1*y_Va*(1+ln(y_i2) + v_i2*y_Va*(1+ln(y_i1) + 
!                  [P*(1/y_i1**2)]            ..last term zero unless i1=i2  OK
! d2S/dy_idy_j   = v_i*(1+ln(y_j)) + (-v_j)*(1+ln(y_i))                      OK
! d2S/dy_idy_Va  = v_i*(1+ln(y_Va)) + v_i*S1 + Q*(1+ln(y_i))                 OK
! d2S/dy_idy_k   = v_i*(1+ln(y_k))                                           OK
! d2S/dy_j1d_j2  = [only Q/y**2 if j1=j2]                                    OK
! d2S/dy_jdy_Va  = zero                                                      OK
! d2S/dy_jdy_k   = zero                                                      OK
! d2S/dy_Va2     = [only Q/y_Va**2]                                          OK
! d2S/dy_Vady_k  = zero                                                      OK
! d2S/dy_k1dy_k2 = [only Q/y_k1**2 if k1=k2]                                 OK
! ----------------------
! the coding is not optimal for speed, all the 1/y**2 term calculated above
! i2slx(1) is index of vacancy, i2slx(2) is index of first neutral
! if either (or both) are missing their index is higher than last constituent
!   write(*,102)'3X va+neutral: ',i2slx
!102 format(a,10i3)
! dpqdy is calculated in gtp3X: set_constitution ??
!   write(*,108)'3X dpqdy: ',(phvar%dpqdy(j1),j1=1,nkl(1)+nkl(2))
108 format(a,10F7.3)
   cation: do j1=1,nkl(1)
! this was an attempt to improve convergence ... it did but not enough
!      if(localmoded.eq.1) goto 109
      cation2: do j2=j1,nkl(1)
! d2S/dy_i1dy_i2 = v_i1*y_Va*(1+ln(y_i2) + v_i2*y_Va*(1+ln(y_i1) + 
!                  [P*(1/y_i1**2)]         ..last term already calculated  OK
!         write(*,103)'3X ij: ',j1,j2,ixsym(j1,j2),yva,&
!              phvar%d2gval(ixsym(j1,j2),1),&
!              phvar%dpqdy(j1),phvar%dpqdy(j2),&
!              phvar%dgval(1,j1,1),phvar%dgval(1,j2,1)
!103      format(a,3i3,6(1pe11.3))
         phvar%d2gval(ixsym(j1,j2),1)=phvar%d2gval(ixsym(j1,j2),1)+&
              (phvar%dpqdy(j1)*phvar%dgval(1,j2,1)+&
               phvar%dpqdy(j2)*phvar%dgval(1,j1,1))*yva/phvar%sites(1)
      enddo cation2
      anion2: do kk=1,nkl(2)
         j2=nkl(1)+kk
         if(j2.lt.min(i2slx(1),i2slx(2))) then
! d2S/dy_idy_j   = v_i*(1+ln(y_j)) + (-v_j)*(1+ln(y_i))    ...cation+anion OK
            phvar%d2gval(ixsym(j1,j2),1)=&
                 phvar%dpqdy(j1)*phvar%dgval(1,j2,1)/phvar%sites(2)+&
                 phvar%dpqdy(j2)*phvar%dgval(1,j1,1)/phvar%sites(1)
         elseif(j2.eq.i2slx(1)) then
! d2S/dy_idy_Va  = v_i*(1+ln(y_Va)) + v_i*S1 + Q*(1+ln(y_i))   ...cation+Va OK
            phvar%d2gval(ixsym(j1,j2),1)=&
                 phvar%dpqdy(j1)*phvar%dgval(1,j2,1)/phvar%sites(2)+&
                 phvar%dpqdy(j1)*spart(1)+&
                 phvar%sites(2)*phvar%dgval(1,j1,1)/phvar%sites(1)
         else
! d2S/dy_idy_k   = v_i*(1+ln(y_k))                        ...cation+neutral OK
!            write(*,107)'3X i,va: ',j1,j2,phvar%dpqdy(j1),phvar%dgval(1,j2,1),&
!                 phvar%sites(2)
!107         format(a,2i2,6(1pe12.4))
            phvar%d2gval(ixsym(j1,j2),1)=&
                 phvar%dpqdy(j1)*phvar%dgval(1,j2,1)/phvar%sites(2)
         endif
      enddo anion2
109   continue
! this done at the end as original dgval(1,j1,1)=P*(1+ln(y_j1))/P used above
! dS/dy_i        = +v_i*S2 + v_i*y_Va*S1 + [P*(1+ln(y_i)]   ..cation   OK
!      write(*,19)'3X c: ',j1,phvar%dgval(1,j1,1),&
!           phvar%dpqdy(j1),spart(2),phvar%dpqdy(j1),yva,spart(1)
!19    format(a,i3,6(1pe12.4))
      phvar%dgval(1,j1,1)=phvar%dgval(1,j1,1)+&
           phvar%dpqdy(j1)*spart(2)+phvar%dpqdy(j1)*yva*spart(1)
   enddo cation
! this done separately as original dgval(1,j2,1)=Q*(1+ln(y_j2))/Q used above
! kall here should be total number of constituents
   anion1: do j2=nkl(1)+1,min(i2slx(1),kall)
      if(j2.lt.min(i2slx(1),i2slx(2))) then
! dS/dy_j        = -v_j*S1 +               [Q*(1+ln(y_j))]  ..anion    OK
!         write(*,*)'3X anion1 A: ',j2
         phvar%dgval(1,j2,1)=phvar%dgval(1,j2,1)+phvar%dpqdy(j2)*spart(1)
      elseif(j2.eq.i2slx(1)) then
! dS/dy_Va       = Q*S1 +                  [Q*(1+ln(y_Va))] ..Va       OK
!         write(*,*)'3X anion1 B: ',j2
         phvar%dgval(1,j2,1)=phvar%dgval(1,j2,1)+phvar%sites(2)*spart(1)
!      else
! dS/dy_k        = nothing +               [Q*(1+ln(y_k)]   ..neutral  OK
      endif
!      write(*,*)'3X anion1 C: ',j2
   enddo anion1
! set temperature derivative of dG/dy
   do j1=1,kall
      phvar%dgval(2,j1,1)=phvar%dgval(1,j1,1)/tval
   enddo
900 continue
!  phvar%gval(1,1)=phvar%gval(1,1)+phvar%sites(ll)*ss
!   write(*,905)'3X parts: ',phvar%gval(1,1),phvar%sites,spart
!905 format(a,6(1pe12.4))
! set temperature derivative of G
   phvar%gval(2,1)=phvar%gval(1,1)/tval
1000 continue
   return
 end subroutine config_entropy_i2sl

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine push_pyval(pystack,intrec,pmq,pyq,dpyq,d2pyq,moded,iz)
! push data when entering an interaction record
   implicit none
   integer pmq,moded,iz
   double precision pyq,dpyq(iz),d2pyq(iz*(iz+1)/2)
   type(gtp_pystack), pointer :: pystack
   type(gtp_interaction), pointer :: intrec
!\end{verbatim} %+
   type(gtp_pystack), pointer :: new
!
   if(associated(pystack)) then
      allocate(new)
      new%previous=>pystack
      pystack=>new
   else
      allocate(pystack)
      nullify(pystack%previous)
   endif
! save data
   pystack%intrecsave=>intrec
   pystack%pmqsave=pmq
   pystack%pysave=pyq
   if(moded.ge.1) then
! if moded 0 there are no derivatives
      allocate(pystack%dpysave(iz))
      pystack%dpysave=dpyq
      if(moded.eq.2) then
! if moded 1 there are no second derivatives
         allocate(pystack%d2pysave(iz*(iz+1)/2))
         pystack%d2pysave=d2pyq
      endif
   endif
1000 continue
   return
 end subroutine push_pyval

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine pop_pyval(pystack,intrec,pmq,pyq,dpyq,d2pyq,moded,iz)
! pop data when entering an interaction record
   implicit none
   integer iz,pmq,moded
   double precision pyq,dpyq(iz),d2pyq(iz*(iz+1)/2)
   type(gtp_pystack), pointer :: pystack
   type(gtp_interaction), pointer :: intrec
!\end{verbatim}
   type(gtp_pystack), pointer :: old
   if(.not.associated(pystack)) then
!      write(*,*)'3X Tying to pop from an empty PY stack'
      gx%bmperr=4075; goto 1000
   endif
! restore data
   intrec=>pystack%intrecsave
   pmq=pystack%pmqsave
   pyq=pystack%pysave
   if(moded.ge.1) then
! if moded >0 there are derivatives
      dpyq=pystack%dpysave
      if(moded.eq.2) then
! if moded 2 there are second derivatives
         d2pyq=pystack%d2pysave
      endif
   endif
! release memory
   old=>pystack
   pystack=>pystack%previous
   deallocate(old)
1000 continue
   return
 end subroutine pop_pyval

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine calc_disfrac(lokph,lokcs,ceq)
! calculate and set disordered set of fractions from sitefractions
! The first derivatives are dxidyj.  There are no second derivatives
!   TYPE(gtp_fraction_set), pointer :: disrec
   implicit none
   integer lokph,lokcs
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
!   TYPE(gtp_fraction_set), pointer :: disrec
   TYPE(gtp_phase_varres), pointer :: phord
   TYPE(gtp_phase_varres), pointer :: phdis
!   logical ordered
! minimum difference in site fraction to be set as ordered
!   double precision, parameter :: yminord=1.0D-10
   integer lokdis,is
!
!   write(*,*)'3X entering calc_disfrac'
! this is the record with the ordered constitution
   phord=>ceq%phase_varres(lokcs)
!   disrec=phord%disfra
!   lokdis=disrec%varreslink
!   phdis=>disrec%phdapointer
   lokdis=ceq%phase_varres(lokcs)%disfra%varreslink
   phdis=>ceq%phase_varres(lokdis)
! this is a record within the ordered constitution record for disordered fracs
!   disrec=>phord%disfra
! to find the varres record with disordered fractions use varreslink
! this is the index to the phase_varres record with the ordered fractions ???
   lokdis=ceq%phase_varres(lokcs)%disfra%varreslink
!   write(*,*)'3X Calc disfra: ',lokph,lokcs,lokdis
!   phdis=>ceq%phase_varres(lokdis)
!   call calc_disfrac2(ceq%phase_varres(lokcs)%disfra,&
!   call calc_disfrac2(ceq%phase_varres(lokcs),ceq%phase_varres(lokdis),ceq)
   call calc_disfrac2(phord,phdis,ceq)
1000 continue
   return
 end subroutine calc_disfrac

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine calc_disfrac2(phord,phdis,ceq)
! calculate and set disordered set of fractions from sitefractions
! The first derivatives are dxidyj.  There are no second derivatives
!   TYPE(gtp_fraction_set), pointer :: disrec
   implicit none
   TYPE(gtp_phase_varres), target :: phord
   TYPE(gtp_phase_varres), pointer :: phdis
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   TYPE(gtp_fraction_set), pointer :: disrec
   logical ordered
! minimum difference in site fraction to be set as ordered
   double precision, parameter :: yminord=1.0D-10
   integer lokdis,is
!
!   write(*,*)'3X entering calc_disfrac'
!   disrec=phord%disfra
!   lokdis=disrec%varreslink
!   phdis=>disrec%phdapointer
! this is the record with the ordered constitution
!   phord=>ceq%phase_varres(lokcs)
! this is a record within the ordered constitution record for disordered fracs
   disrec=>phord%disfra
! to find the varres record with disordered fractions use varreslink
! this is the index to the phase_varres record with the ordered fractions ???
!   lokdis=disrec%varreslink
!   phdis=>ceq%phase_varres(lokdis)
!   write(*,*)'3X calc_disfrac 1A'
! check that some values are accessable
!   write(*,*)'3X calc_disfra phase index: ',phord%phlink
!   write(*,*)'3X calc_disfra disordered sublattices: ',disrec%ndd
!   write(*,*)'3X calc_disfra ordered and disordered records: ',lokcs,lokdis
!   write(*,*)'3X calc_disfra phase index via disordred record: ',phdis%phlink
!   write(*,*)'3X calc_disfrac 1B'
   phdis%yfr=zero
!   write(*,*)'3X disfrac 1: ',disrec%tnoofyfr
   do is=1,disrec%tnoofyfr
      phdis%yfr(disrec%y2x(is))=&
           phdis%yfr(disrec%y2x(is))+disrec%dxidyj(is)*phord%yfr(is)
!      write(*,77)'3X disfrac 2: ',is,disrec%y2x(is),phdis%yfr(disrec%y2x(is)),&
!           disrec%dxidyj(is),phord%yfr(is)
77    format(a,2i3,3(1pe12.4))
   enddo
!   write(*,*)'3X calc_disfrac 2'
! check if phase is really ordered, meaning that the disordered fractions
! are equal to the ordered ones
   ordered=.false.
   do is=1,disrec%tnoofyfr
      if(abs(phdis%yfr(disrec%y2x(is))-&
           phord%yfr(is)).gt.yminord) ordered=.true.
   enddo
!   write(*,*)'3X calc_disfrac 3'
   if(.not.ordered) then
! if this bit set one will not calculate the ordered part of the phase
      phord%status2=ibclr(phord%status2,csorder)
   else
! bit must be cleared as it might have been set at previous call
      phord%status2=ibset(phord%status2,csorder)
   endif
!   write(*,*)'3X calc_disfrac 4: ',phord%status2
! copy these to the phase_varres record that belongs to this fraction set
! a derivative dGD/dyj = sum_i dGD/dxi * dxidyj
! where dGD/dxi is dgval(1,y2x(j),1) and dxidyj is disrec%dxidyj(j)
! because each y constituent contributes to only one disordered x fraction
1000 continue
   return
! G(tot)    = GD(xdis)+(GO(yord)-GO(yord=xdis))
! G(tot).yj = dGD(xdis).dxi*dxdyj + GO.yj - GO.yj ... 
 end subroutine calc_disfrac2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine disordery(phvar,ceq)
! sets the ordered site fractions in FCC and other order/disordered phases
! equal to their disordered value in order to calculate and subtract this part
! phvar is pointer to phase_varres for ordered fractions
   implicit none
   TYPE(gtp_phase_varres), pointer :: phvar
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   TYPE(gtp_fraction_set), pointer :: disrec
!   TYPE(gtp_phase_varres) :: phdis
   TYPE(gtp_phase_varres), pointer :: phdis
   integer lokdcs,ii,nofc1,nofc2
   double precision xxx
!   integer lokdcs,kk,ll,is,nis,nsl
! find disordered fractions
!   lokdcs=phvar%disfra%varreslink
!   disrec=>phvar%disfra
!   write(*,*)'3X disordery: ',disrec%latd,disrec%nooffr(1),lokdcs
!   phdis=ceq%phase_varres(lokdcs)
!   write(*,*)'3X disordery: ',ceq%xconv
!   write(*,*)'3X disordery: ',phdis%yfr(1)
!   phdis=>ceq%disrec%phdapointer
! find disordered fractions
   disrec=>phvar%disfra
   lokdcs=phvar%disfra%varreslink
!   write(*,9)trim(phlista(phvar%phlink)%name),lokdcs
9  format('3X diordery: ',a,i5)
! problem that this pointer is not always ok ....???
   phdis=>ceq%phase_varres(lokdcs)
   call disordery2(lokdcs,phvar,disrec,ceq)
!   write(*,11)'3X phvary: ',(phvar%yfr(ii),ii=1,nofc1)
!   write(*,11)'3X phdisy: ',(phdis%yfr(ii),ii=1,nofc2)
!   write(*,*)'3X Done disorder'
!
1000 continue
   return
 end subroutine disordery

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim} %-
 subroutine disordery2(lokdcs,phvar,disrec,ceq)
! subroutine disordery2(phdis,phvar,disrec,ceq)
! sets the ordered site fractions in FCC and other order/disordered phases
! equal to their disordered value in order to calculate and subtract this part
! phvar is pointer to phase_varres for ordered fractions
! phdis is pointer to phase_varres for disordered fractions
   implicit none
   TYPE(gtp_phase_varres), pointer :: phvar
   TYPE(gtp_fraction_set), pointer :: disrec
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokdcs,kk,ll,is,nis,nsl
   double precision xxx
! copy fractions, loop through all ordered sublattices in phvar
! and store fraction from lokdis
!   write(*,*)'3X dis2: 1'
   kk=0
! this was never assigned!! BOS 16.11.04
   lokdcs=disrec%varreslink
!   write(*,*)'3X lokdcs: ',lokdcs,allocated(ceq%phase_varres(lokdcs)%yfr)
! here copy: 
! y(ord,1,1)=y(dis,1); y(ord,1,2)=y(dis,2); y(ord,1,3)=y(dis,3); 
! y(ord,2,1)=y(dis,1); y(ord,2,2)=y(dis,2); y(ord,2,3)=y(dis,3); 
!   write(*,*)'3X dis2: 2',disrec%latd,disrec%nooffr(1)
!   write(*,*)'3X dis2: 3',phdis%yfr(1)
!   write(*,*)'3X disordery2: ',lokdcs
   do ll=1,disrec%latd
      do is=1,disrec%nooffr(1)
         kk=kk+1
!         phvar%yfr(kk)=phdis%yfr(is)
!         xxx=phdis%yfr(is)
         xxx=ceq%phase_varres(lokdcs)%yfr(is)
         phvar%yfr(kk)=xxx
      enddo
   enddo
!   write(*,*)'3X dis2: 4',disrec%ndd
   if(disrec%ndd.eq.2) then
! one can have 2 sets of ordered subl. like (Al,Fe)(Al,Fe)...(C,Va)(C,Va)...
      nis=disrec%nooffr(1)
      nsl=size(phvar%sites)
!      write(*,*)'3X dis2: 5',nis,nsl
!      write(*,*)'3X dy: ',nis,kk,disrec%latd,nsl,disrec%nooffr(2)
      do ll=disrec%latd+1,nsl
         do is=1,disrec%nooffr(2)
            kk=kk+1
!            phvar%yfr(kk)=phdis%yfr(nis+is)
            phvar%yfr(kk)=ceq%phase_varres(lokdcs)%yfr(nis+is)
         enddo
      enddo
   endif
1000 continue
   return
 end subroutine disordery2

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine set_driving_force(iph,ics,dgm,ceq)
! set the driving force of a phase explicitly
   implicit none
   type(gtp_equilibrium_data), pointer :: ceq
   integer iph,ics
   double precision dgm
!\end{verbatim}
   integer lokph,lokcs
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   ceq%phase_varres(lokcs)%dgm=dgm
1000 continue
   return
 end subroutine set_driving_force

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine extract_massbalcond(tpval,xknown,antot,ceq)
! extract T, P,  mol fractions of all components and total number of moles
! for use when minimizing G for a closed system.  Probably redundant
   implicit none
   double precision, dimension(*) :: tpval,xknown
   double precision antot
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! This routine MUST return error 4143 or 4144 (too few or too many conditions)
! if that is the fact.  Other error codes can be returned if there are
! conditions which does not allow the grid minimizer.
   integer, dimension(4) :: indices
   double precision, dimension(maxel) :: ani,abi,xset,wset
   double precision mass,h298,s298,xxx,xsum,wsum
   double precision sumwdivm,anisum,abisum,restmass,divisor,dividend,abtot
   TYPE(gtp_condition), pointer :: current,last
   character encoded*16,actual_arg(1)*16,elsym*2,elname*16,refstat*16
   integer nox,now,nc,jl,iref,iunit,ip,idf,ie,more,numberest,istv,localerr
   logical allmassbal
!
   do ie=1,noel()
      xknown(ie)=zero
   enddo
!   write(*,*)'3X in extract_massbal'
   ani=zero; abi=zero; xset=zero; wset=zero
   antot=zero; abtot=zero
   xsum=zero; wsum=zero
   anisum=zero; abisum=zero
   nox=0; now=0
   localerr=0
!
!   write(*,*)"3X in extract massbalace 1"
   last=>ceq%lastcondition
   if(.not.associated(last)) then
      gx%bmperr=4143; goto 1000
   endif
!   write(*,*)"3X in extract massbalace 2"
   current=>last
   nc=0
   allmassbal=.TRUE.
100 continue
      current=>current%next
! ignore inactive conditions
      if(current%active.ne.0) goto 300
! if a conditions has several terms we cannot calculate x
      if(current%noofterms.gt.1) then
!         write(*,*)'3X Grid minimizer cannot be used with expressions'
         localerr=4179
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
      if(current%symlink1.gt.0) then
! the value is a symbol, the node to the expression is in
! svflista(current%symlink1)%linkpnode
! NOTE THIS IS NOT THE SAME AS meq_evaluate_svfun but OK as no derivative
         xxx=evaluate_svfun_old(current%symlink1,actual_arg,1,ceq)
      else
         xxx=current%prescribed
      endif
!      write(*,17)'3X massbal: ',encoded,istv,indices,iunit,iref,xxx
17    format(a,2x,a,2x,i3,2x,4i3,2x,2i3,1PE15.7)
! extract values of T, P, N, B, X and W
      if(current%statev.eq.1) then
! this is the temperature
         tpval(1)=xxx
         nc=nc+1
      elseif(current%statev.eq.2) then
! this is the pressure
         tpval(2)=xxx
         nc=nc+1
      elseif(current%statev.eq.110) then
! this is N=value or N(element)=value
         if(indices(2).gt.0) then
! this should mean the number of moles of a component in a phase, illegal here
!            write(*,*)'3X N with 2 indices illegal in this case'
            localerr=4179
         elseif(indices(1).gt.0) then
! N(i)=xxx
            ani(indices(1))=xxx
            anisum=anisum+xxx
         else
! N=xxx
            antot=xxx
         endif
         nc=nc+1
      elseif(current%statev.eq.111) then
         if(indices(2).gt.0) then
            localerr=4179; goto 1000
         endif
! this is X(index1)=value, CHECK UNIT if %!!!
         if(iunit.eq.100) xxx=1.0D-2*xxx
         xset(current%indices(1,1))=xxx
         xsum=xsum+xxx
         nc=nc+1
         nox=nox+1
      elseif(current%statev.eq.120) then
! this is B=value or B(i)=value
         if(indices(2).gt.0) then
! this should mean the mass of a component in a phase, illegal here
            write(*,*)'3X B with 2 indices illegal'
            localerr=4179
         elseif(indices(1).gt.0) then
! B(i)=xxx
            abi(indices(1))=xxx
            abisum=abisum+xxx
         else
! B=xxx
            abtot=xxx
         endif
         nc=nc+1
      elseif(current%statev.eq.122) then
         if(indices(2).gt.0) then
            localerr=4179
         endif
! this is W(index1)=value, CHECK UNIT if %!!!
         if(iunit.eq.100) xxx=1.0D-2*xxx
         wset(current%indices(1,1))=xxx
         wsum=wsum+xxx
         nc=nc+1
         now=now+1
      else
! this is not a massbalance condition but continue just to check how many cond
         allmassbal=.FALSE.
         nc=nc+1
      endif
! take next condition if we have not done all
300   continue
      if(ocv()) write(*,310)'3X massbal: ',current%prescribed,last%prescribed
310   format(a,6(1pe12.4))
      if(.not.associated(current,last)) goto 100
!--------------------------------------
! check if correct number of conditions found
500 continue
   idf=noofel+2-nc
   if(idf.ne.0) then
! if idf is not zero there are not enough conditions
      gx%bmperr=4144; goto 1000
   elseif(.not.allmassbal) then
! some conditions are not massbalance
      localerr=4151
   endif
!   write(*,*)'3X extract_massbal: ',localerr
! We have correct number of conditions but if localerr set we do not have
! all as massbalance conditions.  Return with that code set
   if(localerr.ne.0) then
      gx%bmperr=localerr; goto 1000
   endif
! we have extracted all conditions N, B, X, W
! check that only one value per component
   do ie=1,noel()
      if(xset(ie).gt.zero) then
         if(wset(ie).gt.zero) goto 1100
         if(ani(ie).gt.zero) goto 1100
         if(abi(ie).gt.zero) goto 1100
      elseif(wset(ie).gt.zero) then
         if(ani(ie).gt.zero) goto 1100
         if(abi(ie).gt.zero) goto 1100
      elseif(ani(ie).gt.zero) then
         if(abi(ie).gt.zero) goto 1100
      elseif(abi(ie).le.zero) then
! this can be "the rest"
!         write(*,*)'3X massbal',ie,abi(ie),antot,abtot
         if(.not.btest(globaldata%status,GSNOTELCOMP)) then
            if(antot.eq.zero .and. abtot.eq.zero) goto 1105
!         else
!            write(*,*)'3X Other components then elements 1'
         endif
      endif
   enddo
!   write(*,510)'N: ',(ani(i),i=1,noel())
!   write(*,510)'B: ',(abi(i),i=1,noel())
!   write(*,510)'x: ',(xset(i),i=1,noel())
!   write(*,510)'w: ',(wset(i),i=1,noel())
510 format(a,7F9.6)
   bigif: if(antot.gt.zero) then
! we have a value for total number of moles, N, there must not be one for B
      if(abtot.ne.zero) goto 1110
      more=0
      numberest=0
      sumwdivm=zero
! convert as much as possible to N(i).  Sum also some data needed if there
! are conditions on mass fractions
      do ie=1,noel()
         call get_element_data(ie,elsym,elname,refstat,mass,h298,s298)
         if(xset(ie).gt.zero) then
            ani(ie)=antot*xset(ie)
            anisum=anisum+ani(ie)
            abisum=abisum+mass*ani(ie)
         elseif(abi(ie).gt.zero) then
            ani(ie)=abi(ie)/mass
            anisum=anisum+ani(ie)
            abisum=abisum+mass*ani(ie)
         elseif(wset(ie).gt.zero) then
            sumwdivm=sumwdivm+wset(ie)/mass
            more=1
         elseif(ani(ie).eq.zero) then
            if(numberest.gt.0) then
               write(*,*)'3X Missing condition for two elements.'
! ??               gx%bmperr=0; goto 1000
               gx%bmperr=4151; goto 1000
            endif
            restmass=mass
            numberest=ie
         endif
      enddo
      if(numberest.eq.0) then
         write(*,*)'3X Error - condition on all elements and N??'
         gx%bmperr=0; goto 1000
      endif
      if(more.gt.0) then
! there are some mass fractions, we have to calculate B
! but first we must determine the number of moles of "the rest" element
         divisor=antot-anisum-abisum/(one-wsum)*sumwdivm
         dividend=one+restmass/(one-wsum)*sumwdivm
         ani(numberest)=divisor/dividend
         abi(numberest)=restmass*ani(numberest)
         abisum=abisum+abi(numberest)
! now calculate B
         abtot=abisum/(one-wsum)
!         write(*,520)'3X nrest: ',numberest,divisor,dividend,ani(numberest),&
!              abi(numberest),abtot
520 format(a,i3,6(1pe12.4))
! now calculate moles of elements with massfractions
         do ie=1,noel()
            if(wset(ie).gt.zero) then
               abi(ie)=abtot*wset(ie)
               call get_element_data(ie,elsym,elname,refstat,mass,h298,s298)
               ani(ie)=abi(ie)/mass
            endif
         enddo
      else
! all conditions are mole fractions, just set "the rest"
         ani(numberest)=antot-anisum
      endif
      do ie=1,noel()
         xset(ie)=ani(ie)/antot
      enddo
   elseif(abtot.gt.zero) then
! we have a value for total mass, B, not common and too complicated
!      write(*,*)'3X Cannot handle condition on total mass'
      gx%bmperr=4180
   elseif(xsum.eq.zero .and. wsum.eq.zero) then
! just N(i)= and B(i)=, no N= nor B= and no X nor W, No rest element
!      write(*,520)'3X N(i): ',0,anisum,(ani(j),j=1,noel())
      do ie=1,noel()
         if(abi(ie).gt.zero) then
            call get_element_data(ie,elsym,elname,refstat,mass,h298,s298)
            ani(ie)=abi(ie)/mass
            anisum=anisum+ani(ie)
         endif
      enddo
      antot=anisum
      do ie=1,noel()
         xset(ie)=ani(ie)/antot
         if(xset(ie).le.zero) then
            if(.not.btest(globaldata%status,GSNOTELCOMP)) then
! when other components than elements the mole fractions can be <0 or > 1
               write(*,*)'3X mass balance error: ',ie
               gx%bmperr=4181; goto 1000
!            else
!               write(*,*)'3X Other components than elements 2'
            endif
         endif
      enddo
   else
! any other combination of conditions ....
      write(*,*)'3X Cannot handle these massbalance conditions'
      gx%bmperr=4182
   endif bigif
! copy fractions to arguments
900 continue
   do ie=1,noel()
      xknown(ie)=xset(ie)
   enddo
1000 continue
   return
! errors
1100 continue
   write(*,*)'3X Two mass balance conditions for same element',ie
   gx%bmperr=4183; goto 1000
1105 continue
   write(*,*)'3X One component without condition'
   gx%bmperr=4181; goto 1000
1110 continue
   write(*,*)'3X Both N and B cannot be set'
   gx%bmperr=4184; goto 1000
!
 end subroutine extract_massbalcond

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine save_constitutions(ceq,copyofconst)
! copy the current phase amounts and constituitions to be restored
! if calculations fails during step/map
! DANGEROUS IF NEW COMPOSITION SETS CREATED
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision, allocatable, dimension(:) :: copyofconst
!\end{verbatim} %+
   integer varresx,nz,ij,syfr
! calculate dimension of copyofconst
   nz=0
! skippa varres with index 1, that is the reference phase
!   do varresx=2,csfree-1
   do varresx=2,highcs
      if(allocated(ceq%phase_varres(varresx)%yfr)) then
! NOTE size( ... ) can return reasonable value even if not allocated !!!
! BUT why is phas_varres(varresx)%yfr it not allocated ???
! evidently the composition set for varresx is created ... maybe removed??
         syfr=size(ceq%phase_varres(varresx)%yfr)
      else
         syfr=0
      endif
!      write(*,12)'3X Varres record and size: ',varresx,1+syfr,nz
12    format(a,3i5)
      nz=nz+1+syfr
   enddo
!   write(*,*)'3X In save_constitution',nz,highcs
   allocate(copyofconst(nz))
   nz=1
!   do varresx=2,csfree-1
   do varresx=2,highcs
! save 1+sfr values for each composition set
! segmentation fault in this loop for stepbug (>20 elements COST507)
! crash happends when higher composition sets are stored ...
      copyofconst(nz)=ceq%phase_varres(varresx)%amfu
      if(allocated(ceq%phase_varres(varresx)%yfr)) then
         syfr=size(ceq%phase_varres(varresx)%yfr)
      else
         syfr=0
      endif
!      write(*,16)'3X Storing varres record: ',varresx,syfr,size(copyofconst),nz
16    format(a,5i5)
! the segmentation fault seems not to be the allocation of copyofconst but
! rather that we cannot access the yfr in ceq%phase_varres(varresx)
! for the extra composition sets created by gridmin
      do ij=1,syfr
         copyofconst(nz+ij)=ceq%phase_varres(varresx)%yfr(ij)
      enddo
!      write(*,17)varresx,nz,1+syfr,(copyofconst(ij),ij=nz,nz+syfr)
17    format('3X s:',3i4,10(F6.3))
      nz=nz+1+syfr
   enddo
1000 continue
   return
 end subroutine save_constitutions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine restore_constitutions(ceq,copyofconst)
! restore the phase amounts and constitutions from copyofconst
! if calculations fails during step/map
! DANGEROUS IF NEW COMPOSITION SETS CREATED
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision copyofconst(*)
!\end{verbatim}
   integer nz,varresx,ij,syfr
   nz=1
! skippa varres with index 1, that is the reference phase
!   do varresx=2,csfree-1
   do varresx=2,highcs
      ceq%phase_varres(varresx)%amfu=copyofconst(nz)
      if(allocated(ceq%phase_varres(varresx)%yfr)) then
         syfr=size(ceq%phase_varres(varresx)%yfr)
      else
         syfr=0
      endif
      do ij=1,syfr
         ceq%phase_varres(varresx)%yfr(ij)=copyofconst(nz+ij)
      enddo
!      write(*,17)varresx,nz,syfr,ceq%phase_varres(varresx)%amfu,&
!           (ceq%phase_varres(varresx)%yfr(ij),ij=1,syfr)
17    format('3X r:',i2,2i3,6(1pe12.4))
      nz=nz+1+size(ceq%phase_varres(varresx)%yfr)
   enddo
1000 continue
   return
 end subroutine restore_constitutions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

