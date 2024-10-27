!
! gtp3X included in gtp3.F90, separate gtp3XQ for MQMQA
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     15. Section: calculate G and other things
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calcg
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
!   write(*,*)'3X in calcg',iph,ics,moded
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
!   write(*,*)'3X in calcg: ',lokcs
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
!   write(*,*)'3X calcg 7: ',lokres,trim(ceq%eqname)
! call using the local structure phase_varres
! results can be obtained through lokres
!   write(*,17)'3X calcg: ',lokph,lokres,ceq%phase_varres(lokres)%yfr(1)
17 format(a,2i4,1pe15.6)
   call calcg_internal(lokph,moded,ceq%phase_varres(lokres),ceq)
1000 continue
!   write(*,*)'3X back from calcg_internal'
   return
 end subroutine calcg

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calcg_internal
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
   integer fractype,epermut,ipermut,typty,pmq,maxprec,already
   integer sameint(5)
   integer, dimension(permstacklimit) :: lastpmq,maxpmq
!   character bug*60
!   dimension sites(maxsubl),pushpop(maxpp)
   double precision, dimension(:), allocatable :: dpyq(:),d2pyq(:),d2vals(:)
   double precision, dimension(:,:), allocatable :: dvals(:,:)
   double precision vals(6)
! this array has the sum of constituents up to and including current sublattice
   integer incffr(0:maxsubl)
! Kohler-Toop binary excess model link
   type(gtp_tooprec), pointer :: tooprec
! in local gz: gz%intlevel level of interaction, gz%intcon and gz%intlat are
! used also in cgint when calculating interactions.
   TYPE(gtp_parcalc) :: gz
! disordered fraction set
!   TYPE(gtp_fraction_set) :: fracset,dislink
   TYPE(gtp_fraction_set), pointer :: fracset,dislink
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
   integer nofc2,nprop,nsl,msl,lokdiseq,ll,id,id1,id2,lm,qz
   integer lokfun,itp,nz,intlat,ic,jd,jk,ic1,jpr,ipy,i1,j1,jj,jxsym
   integer i2,j2,ider,is,kk,ioff,norfc,iw,iw1,iw2,lprop,jonva,icat
   integer nsit1,nsit2
! cqc configurational entropy
   integer nclust
   double precision, allocatable, dimension(:,:) :: gclust
! mqmqa endmember counting and other specials
   integer mqmqj,kend
!   double precision, dimension(:,:), allocatable :: fhv
!   double precision, dimension(:,:,:), allocatable :: dfhv
!   double precision, dimension(:,:), allocatable :: d2fhv
   double precision g2val(6)
! to handle parameters with wildcard constituent and other things
   logical wildc,nevertwice,first,chkperm,ionicliq,iliqsave,iliqva,iliqneut
! mobility parameters must not have wildcard constituents
   logical liq2state,wildmob,mqmqa
! pointer to mqmqaf record with all fraction records for MQMQA
   type(gtp_mqmqa_var), pointer :: mqf
! debugging for partitioning and ordering
   integer idlist(9)
! calculate RT to normalize all Gibbs energies, ceq is current equilibrium
   rtg=globaldata%rgas*ceq%tpval(1)
   ceq%rtn=rtg
!   write(*,*)'3X in calcg_internal: ',lokph
!   if(ocv()) write(*,*)'3X in gcalc_internal: ',lokph
!-----------------------
   chkperm=.false.
   mqmqa=.false.
   already=0
   if(btest(phlista(lokph)%status1,PHMQMQA)) then
!      write(*,*)'3X phase has MQMQA model'
      mqmqa=.TRUE.
   endif
   if(btest(phlista(lokph)%status1,PHFORD) .or. &
        btest(phlista(lokph)%status1,PHBORD)) then
! PHPALM is needed for phases with permutations such as ordered FCC/BCC/HCP
      chkperm=.true.
      if(.not.btest(phlista(lokph)%status1,PHPALM)) then
!         write(*,*)'3X calling palmtree ',lokph,cps%phtupx
! This is needed only once unless parameters are changed.  It numbers the
! interaction records sequentially for the permutations
! palmtree is in gtp3Y.F90 for some unknown reason ...
         call palmtree(lokph)
         if(gx%bmperr.ne.0) goto 1000
! this must be zeroed if a new interaction parameter is added
!         phlista(lokph)%status1=ibset(phlista(lokph)%status1,PHPALM)
      endif
   endif
!-----------------------------------------------------------------
50  continue
! local work arrays for products of Y and calculated parameters are allocated
   gz%nofc=phlista(lokph)%tnooffr
   nofc2=gz%nofc*(gz%nofc+1)/2
!   write(*,*)'3X in calcginternal ',btest(phlista(lokph)%status1,PHLIQ2STATE)
!   write(*,17)'3X calcg, ',lokph,gz%nofc,nofc2,size(cps%d2gval),cps%nprop,&
!        cps%yfr(1)
!17 format(a,5i4,1pe15.6)
! for disordered fraction sets gz%nofc must be from disordered fraction record
! maybe these should not be allocated for moded=0 and 1
!   write(*,*)'3X allocate: ',gz%nofc,nofc2
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
   elseif(mqmqa) then
! MQMQA FactSage entropy model
! strange error OC dies when calling this using "c g" as first command
! OK when using c ph .... , problem with arguments
!      write(*,222)moded,lokph,phlista(lokph)%tnooffr,phres%yfr(1),gz%tpv(1)
222   format('3X call mqmqa entropy: ',3i3,2(1pe12.4))
! in gtp3_XQ
      call config_entropy_mqmqa1(phres,moded,lokph,gz%tpv(1))
! attempt to simplify the call .... 
! when we come back mqmqaf should have some arrays allocated ....
! DO NOT SET THIS POINTER BEFORE ARRAYS ARE ALLOCATED IN CONFIG_ENTROPY_MQMQA
! this pointer now obsolete, phres%mqmqaf used inside _mqmqa1
!      mqf=>phres%mqmqaf
!     write(*,'(a,1pe12.4)')'3X MQMQA cfgG:',phres%gval(1,1)*gz%rgast*phres%amfu
!      write(*,'(a,1pe12.4)')'3X MQMQA cfgG:',phres%gval(1,1)*gz%rgast
!      do mqmqj=1,mqf%npair
!         write(*,777)mqf%pair(mqmqj),(mqf%dpair(mqmqj,kk),kk=1,mqf%nconst)
!         write(*,777)mqf%pair(mqmqj),mqf%dpair(mqmqj,1),mqf%dpair(mqmqj,2),&
!              mqf%dpair(mqmqj,3)
777      format('3X pairs:',F10.6,2x,10F9.5)
!      enddo
!      write(*,*)'3X back from config_entropy_mqmqa1 '
   elseif(btest(phlista(lokph)%status1,PHQCE)) then
! this is the corrected QC, Hillert-Selleby-Sundman model
      call config_entropy_qchillert(moded,phlista(lokph)%nooffr(1),&
           phres,phlista(lokph),gz%tpv(1))
!      write(*,480)'3X dg/dt/RT: 1: ',qcmodel,phres%yfr(3),&
!           phres%gval(1,1),phres%gval(2,1)
480   format(a,i2,6(1pe12.4))
! several old versions to be deleted ...
!      call config_entropy_cqc(moded,phlista(lokph)%tnooffr(1),&
!           phres,phlista(lokph),nclust,gclust,gz%tpv(1))

   elseif(btest(phlista(lokph)%status1,PHCVMCE)) then
! the classical quasichemical or tetraherdon CVM model with LRO
      call config_entropy_qcwithlro(moded,phlista(lokph)%tnooffr,phres,&
           phlista(lokph),gz%tpv(1))
! phstate
   elseif(btest(phlista(lokph)%status1,PHTISR)) then
! the configurational model by E Kremer (Calphad 2022)
      call config_entropy_tisr(moded,phlista(lokph)%tnooffr,phres,&
           phlista(lokph),gz%tpv(1))
   elseif(btest(phlista(lokph)%status1,PHSROT)) then
! the configurational model is a modified tetrahedron quasichemical model
      call config_entropy_srot(moded,phlista(lokph)%tnooffr,phres,&
           phlista(lokph),gz%tpv(1))
! phstate
   elseif(btest(phlista(lokph)%status1,PHSSRO)) then
! CVM tetraheron SRO (no LRO) configurational entropy
      call config_entropy_ssro(moded,lokph,phres,gz%tpv(1))
! phstate
   elseif(btest(phlista(lokph)%status1,PHCVMTFL)) then
! CVM tetraheron SRO (no LRO) configurational entropy
      call config_entropy_cvmtfl(moded,lokph,phres,gz%tpv(1))
   else
!----------- the CF Bragg-Williams ideal configurational entropy per sublattice
! NOTE: for phases with disordered fraction set this is calculated
! ONLY for the original constituent fraction set with ordering sublattices
      call config_entropy(moded,nsl,phlista(lokph)%nooffr,phres,gz%tpv(1))
   endif
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3X segmentation fault 10'
!-------------------------------------------------------------------
! MQMQA separate calculation of G as well (as above) the entropy!!!
! NO no 2state model, einstein, magnetism etc.
   if(mqmqa) then
!      write(*,*)'3X Separate routine for nonconfig MQMQA'
      call calc_mqmqa(lokph,phres,ceq)
!      write(*,*)'3X back from nonconfig MQMQA'
      goto 1000
   endif
!-------------------------------------------------------------------
! start BIG LOOP for all fraction variables and parameters
! there may be several different properties in addition to G like TC, MQ& etc
! each of these are stored in separate gval(*,ipy) where ipy is an integer
! set for each property. lprop is incremented by one for each new property
! found (each phase may have different) and in listprop the original type
! of property is stored.  listprop will always be associated with phmain
!100 continue
! yionva is used as indicator below if there are Va or just neutrals ...
   yionva=zero
   nevertwice=.true.
   lprop=2
   phmain%listprop=0
   phmain%listprop(1)=lprop
!   write(*,168)'3X lprop0:',lprop,0,(phmain%listprop(jj),jj=1,10)
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
!      write(*,*)'3X segmentation fault 20',fractype
! return here for calculating with disordered fractions for same fraction type
110 continue
! gz%nofc is number of fraction variables, msl is number of sublattices
! for this set of fractions!!! Ordering in FCC may have 5 sublattices with
! 4 participating in ordering and one interstitial.  The second fraction
! set may have 2 sublattices, 1 for the 4 ordering and one interstitial
!      fracset=phmain%disfra
      fracset=>phmain%disfra
!      write(*,*)'3X segmentation fault 30',associated(fracset)
      ftype: if(fractype.eq.1) then
!---------------------------------------------- ordered (or only) fraction set
         if(btest(phlista(lokph)%status1,PHMFS)) then
! there is a disordered fractions set, we need additional fracset
            if(fracset%totdis.ne.0) then
               if(btest(phlista(lokph)%status1,PHSUBO)) then
! if phsubo set skip subtracting the ordered part as disordered, just add
                  goto 106
               endif
! the phase can totally disorder, if disordered skip ordered part
! the CSORDER bit set by calc_disfrac called from set_constitution
               if(btest(phmain%status2,CSORDER)) then
! the phase is ordered, we have to calculate this part twice
!                  write(*,*)'3X Setting nevertwice false'
                  nevertwice=.false.
! independent if ordered or disordered always calculate first fraction set
               else
! the phase is disordered, skip ordered part and just calculate disordered
! nevertwice is already set TRUE
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
         if(ocv()) write(*,*)'3X Allocated dpyq 2'
         dpyq=zero
         deallocate(dvals)
         deallocate(d2vals)
         allocate(dvals(3,gz%nofc))
         allocate(d2vals(nofc2))
         if(ocv()) write(*,*)'3X Allocated vals 2'
! the results will be stored in result arrays indicated by phres
! for the disordered fraction set phres must be set here and the arrays zeroed
!         dislink=cps%disfra
         dislink=>cps%disfra
!         write(*,*)'3X Calc internal disordred part 1A',dislink%fsites
         lokdiseq=dislink%varreslink
!         write(*,*)'3X Calc internal disordred part 1B',lokdiseq
         phres=>ceq%phase_varres(lokdiseq)
! Wow phres%gval etc not allocated !!
         if(.not.allocated(phres%gval)) then
            allocate(phres%gval(6,nprop))
         endif
         phres%gval=zero
!         write(*,*)'3X Calc internal disordred part 1c',&
!              allocated(phres%dgval),gz%nofc
         if(moded.gt.0) then
            if(.not.allocated(phres%dgval)) then
               allocate(phres%dgval(3,gz%nofc,nprop))
            endif
            phres%dgval=zero
            if(moded.gt.1) then
               nofc2=gz%nofc*(gz%nofc+1)/2
!               write(*,*)'3X segmentation fault 48: ',&
!                    allocated(phres%d2gval),nofc2
               if(.not.allocated(phres%d2gval)) then
                  allocate(phres%d2gval(nofc2,nprop))
               endif
               phres%d2gval=zero
            endif
         endif
!         write(*,*)'3X Calc internal disordred part 2'
      endif ftype
!==========================================================
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! code below is an abandonned test to parallelize the calculation of each
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
      endmemloop: do while(associated(endmemrec))
!
! The array maxpmq is used for interaction permutations.  It must be
! initialized to zero at the first endmember permutation.  It is set to
! limits for the interacton permutations for all interaction records.
         maxpmq=0
         maxprec=0
         epermut=0
         sameint=0
!--------------------------------- quick test of mqmqa reference state
         if(mqmqa) then
            stop '3X MQMQA separate routine'
         endif
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
            wildmob=.FALSE.
            pyqloop: do ll=1,msl
               id=endmemrec%fraclinks(ll,epermut)
! debugging 4SL with wildcards
               idlist(ll)=id
! remove next line when all fixed
!               if(ll.lt.5) clist(ll)=id
! id negative means wildcard, independent of the fraction in this sublattice
               if(id.lt.0) then
                  gz%yfrem(ll)=one
                  wildmob=.TRUE.
               else
                  gz%yfrem(ll)=phres%yfr(id)
                  if(gz%yfrem(ll).lt.bmpymin) gz%yfrem(ll)=bmpymin
                  if(gz%yfrem(ll).gt.one) gz%yfrem(ll)=one
               endif
! gz%endcon is used for interaction parameters below
               gz%endcon(ll)=id
               pyq=pyq*gz%yfrem(ll)
!               write(*,33)ll,epermut,id,gz%yfrem(ll),pyq
33             format('3X py: ',i3,2i5,2(1pe12.4))
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
                  elseif(id.eq.phlista(lokph)%i2slx(2) .and. jonva.eq.0) then
! we have NO vacancy but a neutral in second sublattice
                     iliqva=.FALSE.
                     yionva=-one
                     jonva=0
                     why1: if(.not.iliqsave) then
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
!                             id,yionva,saveg(1,1)
!                     else
!                        write(*,*)'3X neutral: ',jonva,yionva,iliqva
                     endif why1
                  endif
               endif
            enddo pyqloop
            if(moded.eq.0) goto 150
!---------------------------------------------------- first derivatives of py
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
!---------------------------------------------------- second derivatives of py
! searching for bug with interaction wildcards in 4SL
!            write(*,68)'3X d2P/dyi2A:',nofc2,(d2pyq(id),id=1,nofc2)
!            d2pyq is all zero here
            d2pyqloop1: do ll=1,msl
               id1=endmemrec%fraclinks(ll,epermut)
! too complicated here ...               jxsym=ixsym(ll,ll+1)
               d2pyloop2: do lm=ll+1,msl
                  id2=endmemrec%fraclinks(lm,epermut)
                  if(id1.gt.0) then
                     if(id2.gt.0) then
                        d2pyq(ixsym(id1,id2))=dpyq(id1)/gz%yfrem(lm)
                     else
! wildcard in sublattice lm, real component in ll
!                        do iw=incffr(lm)+1,incffr(lm)
!                           d2pyq(ixsym(id1,iw))=dpyq(id1)
!                        enddo
! This derivative should be zero!! /170324/BoS
                        continue
                        wildmob=.TRUE.
                     endif
                  else
! wildcard in sublattice ll, real component in lm
                     if(id2.gt.0) then
!                        do iw=incffr(ll-1)+1,incffr(ll)
!                           d2pyq(ixsym(id2,iw))=one
!                        enddo
! This should be zero!! /170324/BoS
                        continue
                     else
! wildcards in both sublattice ll and lm
!                        do iw1=incffr(ll-1)+1,incffr(ll)
!                           do iw2=incffr(lm-1)+1,incffr(lm)
!                              d2pyq(ixsym(iw1,iw2))=pyq
!                           enddo
!                        enddo
! I think this should be zero too!! /170324/BoS
                     endif
                     wildmob=.TRUE.
                  endif
               enddo d2pyloop2
            enddo d2pyqloop1
! searching for bug with interaction wildcards in 4SL
!            write(*,67)'B:' ,pyq,(idlist(iw1),iw1=1,nsl)
67          format('3X endmem',a,e12.4,9i4)
!            write(*,68)'3X d2P/dyi2B:',nofc2,(d2pyq(id),id=1,nofc2)
68          format(a,i3,5(1pe12.4)/,(16x,5e12.4))
!---- jump here if moded is 0 or 1
150         continue
!
! if debugpar nonzero add call for debug_endmemberpar(....) here
!
!-----------------------------------------------------
! d2pyq contains 2nd serivatives of endmember fractions.
!            write(*,228)'3X d2pyq 0:',d2pyq
!            write(*,*)'3X Config G 4A: ',phres%gval(1,1)*rtg
!            write(*,154)'3X endmember permutation: ',epermut,(clist(i),i=1,4)
154         format(a,i5,4i4,'--------------------------------')
155         format(a,i5,10i4)
            proprec=>endmemrec%propointer
! for liquids with twostate models first calculate the g2 parameter
            if(btest(phlista(lokph)%status1,PH2STATE)) then
               write(*,*)'3X Phase ',trim(phlista(lokph)%name),&
                    ' has PH2STATE bit set'
               call calc_twostate_model_endmember(proprec,g2val,ceq)
               if(gx%bmperr.ne.0) goto 1000
!               write(*,'(a,6(1pe12.4))')'3X g2val:',g2val
               liq2state=.true.
            else
               liq2state=.false.
               g2val=zero
            endif
            emprop: do while(associated(proprec))
               typty=proprec%proptype
               if(typty.ne.1) then
! if property different from 1 (=G) find where to store it, use phmain link
! First check if the parameter is a mobility and there are wildcrds
                  if(wildmob) then
! nowildcard(1..3) set in gtp_init in gtp3A.F90 for mobility parameters
! typty is indicator*100 + constituent index
                     do qz=1,3
                        if(typty/100.eq.nowildcard(qz)) then
                           write(*,*)&
                                '3X mobilities must not have wildcards',lokph
                           gx%bmperr=4374; goto 1000
                        endif
                     enddo
                  endif
                  do qz=2,lprop-1
                     if(phmain%listprop(qz).eq.typty) goto 170
                  enddo
! a new property, save its typty in listprop and increment lprop
! note that the property index typty is not used as index in gval etc
! as that can be very large. lprop is incremented by 1 for each property
! actually used in the model of the phase.  lprop is last free index
                  qz=lprop
                  if(qz.gt.size(phmain%listprop)) then
                     write(*,*)'Too many differnt parameter identifiers',qz
                     gx%bmperr=4338; goto 1000
                  endif
                  phmain%listprop(qz)=typty
! a bit stupid to allocate listprop, it should have fixed allocation ...
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
! listprop(1) is number of defined properties, listprop(2..) is property
!                  write(*,168)'3X lprop: ',lprop,typty,&
!                       (phmain%listprop(ipy),ipy=1,lprop)
168               format(a,2i5,': ',10i5)
! jump here is we already have found this property and know its ipy
170               continue
                  ipy=qz
               else
                  ipy=1
               endif
!================ here we calculate the endmember parameter ============
! calculate function and derivatives wrt T and P
! the results from eval_tpfun must also be different in different treads ...
               lokfun=proprec%degreelink(0)
               call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
!               write(*,167)'3X eval_tpfun: ',ipy,lokfun,pyq,vals(1),vals(1)/rtg
167            format(a,2i5,6(1pe12.4))
               if(gx%bmperr.ne.0) goto 1000
               prop1: if(ipy.eq.1) then
! property 1 i.e. Gibbs energy, should be divided by RT
                  vals=vals/rtg
                  if(liq2state) then
! if phase has liquid twostate model add g2val!!
!                     write(*,'(a,6(1pe12.4))')'3X +g2val',&
!                          vals(1),g2val(1),vals(1)+g2val(1),&
!                          vals(4),g2val(4),vals(4)+g2val(4)
                     vals=vals+g2val
                  endif
               endif prop1
!               write(*,*)'3X property type: ',typty,ipy,vals(1)
!================ now we calculated the endmember parameter ============
! take care of derivatives of fraction variables ...
!               write(*,173)'3X endmember: ',endmemrec%antalem,ipy,pyq,vals(1)
173            format(a,2i4,4(1pe12.4))
! multiply with py and derivatives. vals is composition independent
!               write(*,*)'3X Config G 4B: ',vals(1)*rtg
! segmentation fault between 64 and 65 ....
               noderz2: if(moded.gt.0) then
                  derloopz2: do id=1,gz%nofc
                     do itp=1,3
                        phres%dgval(itp,id,ipy)=phres%dgval(itp,id,ipy)+ &
                             dpyq(id)*vals(itp)
                     enddo
                     if(moded.gt.1 .and. dpyq(id).gt.zero) then
                        jxsym=kxsym(id,id+1)
                        do jd=id+1,gz%nofc
! trying to replace calls of ixsym ... OK here
                           if(ixsym(id,jd).ne.jxsym) then
                              write(*,*)'ISYM error 1',id,jd,ixsym(id,jd),jxsym
                              stop
                           endif
!                           write(*,*)'3X segfault 64C',allocated(d2pyq),&
!                                d2pyq(jxsym)
!                           write(*,*)'3X segfault 64E',allocated(phres%d2gval)
! phres%d2gval not allocated!!
!                           write(*,*)'3X segfault 64F',phres%d2gval(jxsym,ipy)
                           phres%d2gval(jxsym,ipy)= &
                                phres%d2gval(jxsym,ipy)+ &
                                d2pyq(jxsym)*vals(1)
                           jxsym=jxsym+jd
                        enddo
                     endif
                  enddo derloopz2
               endif noderz2
               do itp=1,6
                  phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*vals(itp)
               enddo
!              write(*,171)'3X phres7: ',ipy,phres%gval(1,1),phres%gval(1,ipy),&
!                    pyq,vals(1)
171            format(a,i3,6(1pe12.4))
! strange values of mobilities for ordered phases ... EINSTEIN
!               if(ipy.ne.1) then
!                  write(*,173)'3X gval:      ',phmain%listprop(ipy),ipy,&
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
! take link to first interaction records, use push and pop to save pyq etc
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
! check if there is a Kohler-Toop link (NOT YET)
!               write(*,*)'3X testing tooprec: ',associated(intrec%tooprec)
               if(associated(intrec%tooprec)) then
!                  write(*,*)'3X Toop/Kohler model: ',&
!                       associated(intrec%tooprec),chkperm,gz%intlevel
                  tooprec=>intrec%tooprec
                  if(chkperm) then
                     write(*,*)'3X Toop/Kohler and permutations is illegal'
                     gx%bmperr=4399; goto 1000
                  endif
! we need this additional information inside calc_toop
! I find it very elegant just to include a pointer to the phase_varres record!
                  tooprec%phres=>cps
               else
                  nullify(tooprec)
               endif
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
! if intlat or ic is zero or less give error message and skip
               if(intlat.le.0 .or. ic.le.0) then
                  if(already.eq.0 .or. intrec%antalint.ne.already) then
                     already=intrec%antalint
                     write(*,231)'3X error: ',phlista(lokph)%alphaindex,&
                          (idlist(iw1),iw1=1,nsl),&
                          (gz%intlat(iw1),gz%intcon(iw1),iw1=1,gz%intlevel)
                     write(*,231)'3X intp: ',intrec%antalint,gz%intlevel,&
                          ipermut,intlat,ic,pmq,maxpmq(pmq)
231                  format(a,10i5)
                  endif
                  goto 290
               endif
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
                  wildmob=.TRUE.
!                  write(*,*)'3X wildcard found!'
                  ymult=gz%yfrint(gz%intlevel)*(one-gz%yfrint(gz%intlevel))
               endif
!               write(*,228)'3X d2pyq 1:',d2pyq
!---------------------------------
!               write(*,*)'3X ionic liquid: ',iliqsave,yionva
               cationintandva: if(.not.iliqsave) then
! iliqsave is TRUE when interaction in first sublattice and Va in second
!                  write(*,228)'3X d2pyq 7:',d2pyq
                  modedx: if(moded.gt.0) then
! ...................................... loop for first derivatives
                     iloop1: do id=1,gz%nofc
                        if(moded.gt.1) then
! ...................................... second derivatives
! For all models except ionic liquids 2nd derivatives are simple ...
                           iloop2B: do jd=id+1,gz%nofc
!                              d2pyq(ixsym(id,jd))=d2pyq(ixsym(id,jd))*ymult
                              jxsym=kxsym(id,jd)
                              d2pyq(jxsym)=d2pyq(jxsym)*ymult
                           enddo iloop2B
                           d2pyq(ixsym(id,ic))=dpyq(id)
                        endif
! I FORGOT THIS LINE WHEN TRYING TO FIX IONIC LIQUID !!! TOTAL MESS !!!
                        dpyq(id)=dpyq(id)*ymult
                     enddo iloop1
                  endif modedx
               else ! here we have cation interaction with Va in second subl.
! SPECIAL FOR IONIC LIQUID
! This is needed for interactions from endmembers with Va in second sublattice
! as the model must be compatibel with a regular solution, like
! (Mo+4,Pd+2,Rh+3)p(Va)q must be identical to (Mo,Pd,Rh) and
! (Fe+2)p(Va,C)q must be identical to (Fe,C)
! This requires that each cation fracition is multiplied with fraction of Va
! Instead of just yA+yB+yVa one must have yA+yB+yVa**2
!                  write(*,228)'3X pyq  0:',pyq
!                  write(*,228)'3X dpyq 1:',dpyq
!                  write(*,228)'3X divers:',ymult,yionva,gz%yfrem(1)
                  if(jonva.le.0 .and. intlat.eq.1) then
                     write(*,*)'Illegal cation interaction with neutral'
                     gx%bmperr=4265; goto 1000
                  endif
! ...................................... loop for first derivatives
                  iliqloop1: do id=1,gz%nofc
                     seconder2: if(moded.gt.1) then
! CODE BELOW IS UNCERTAIN
! This IF loop is only executed when Va in second sublattice, i.e. when cation
! interactions which should also be multiplied with the power of yionva
! which is gz%intlevel+1
! jonva=phlista(lokph)%i2slx(1) is index of vacancy, i2slx(2) is first neutral
! index of the constituent in first sublattice is gz%endcon(1)
! index of the constituent in second sublattice is gz%endcon(2) = jonva
! index of interaction constituents are in gz%intcon(gz%intlevel+)
! pyq, dpyq and d2pyq set for the endmember
!
! NOTE: some 2nd derivatives wrong for (Fe+2)p(Va,C)q and more ...
! NOT tested (Ca+2)p(O-2,SiO4-4,SiO2)q
! ...................................... loop for second derivatives
                        iloop2X: do jd=id+1,gz%nofc
                           jxsym=kxsym(id,jd)
                           if(jd.le.phlista(lokph)%nooffr(1)) then
! both id and jd are cations, interaction must be multiplied with yionva
!                              d2pyq(ixsym(id,jd))=&
!                                   d2pyq(ixsym(id,jd))*ymult*yionva
                              d2pyq(jxsym)=d2pyq(jxsym)*ymult*yionva
!                              write(*,215)gz%intlevel,ic,id,jd,&
!                                   d2pyq(ixsym(id,jd)),ymult,&
!                                   d2pyq(ixsym(id,jd))*ymult*yionva
215                           format('3X d2pyq: ',4i3,4(1pe12.4))
                           elseif(jd.lt.jonva) then
! if jd<jonva derivative wrt anion and cation or two cations, jd must be anion
!                              d2pyq(ixsym(id,jd))=&
!                                   d2pyq(ixsym(id,jd))*ymult
                              d2pyq(jxsym)=d2pyq(jxsym)*ymult
                           elseif(jd.eq.jonva) then
! calculate also d2pyq(ixsym(jonva,jonva)) the only nonzero diagonal element
!                              d2pyq(ixsym(id,jd))=(gz%intlevel+1)/gz%intlevel*&
!                                   d2pyq(ixsym(id,jd))*ymult
                              d2pyq(jxsym)=(gz%intlevel+1)/gz%intlevel*&
                                   d2pyq(jxsym)*ymult
                           else
! second derivatives with two neutrals
!                              d2pyq(ixsym(id,jd))=&
!                                   d2pyq(ixsym(id,jd))*ymult
                              d2pyq(jxsym)=d2pyq(jxsym)*ymult
                           endif
                        end do iloop2X
!                        write(*,216)'3X dpyq: ',id,jd,dpyq
!                        write(*,217)'3X d2pyq:',d2pyq
216                     format(a,2i3,6(1pe10.2))
217                     format(a,6(1pe10.2))
                     endif seconder2
! assigning d2pyq before updating dpyq ??
9991                 continue
                     d2pyq(ixsym(ic,id))=dpyq(id)
! ........ this is the first derivative, must be exact NO CHANGE 17.12.06/BoS
! ic is the constituent index of the interaction
!                     write(*,314)'3X this dpyq1:',id,ic,jonva,ixsym(id,jonva)
314                  format(a,4i4)
                     if(dpyq(id).ne.zero) then
                        dpyq(id)=dpyq(id)*ymult
!                        write(*,216)'3X this dpyq1:',id,ic,dpyq(id),ymult
                     elseif(jonva.gt.0) then
                        if(d2pyq(ixsym(id,jonva)).ne.zero) then
! this is adding more first order derivatives ???
                           dpyq(id)=dpyq(jonva)*ymult
!                        write(*,216)'3X this dpyq2:',id,jonva,dpyq(id),ymult
                        endif
                     endif
                     if(id.eq.phlista(lokph)%i2slx(1) .and. &
                          gz%intlat(gz%intlevel).eq.1) then
! for vacancies there is an additional power  in first subl
                        dpyq(id)=(gz%intlevel+1)*dpyq(id)
! should maybe be:
!                        dpyq(id)=(gz%intlevel+1)/gzintlevel*dpyq(id)
!                        write(*,197)gz%intlevel,gz%intcon(gz%intlevel)
197                     format('3X: Va inter: ',5i3)
                     endif
!                     write(*,216)'3X d2pyq1: ',ic,id,d2pyq(ixsym(ic,id))
                  enddo iliqloop1
! This is a special 2nd derivative wrt Va twice
!                  d2pyq(ixsym(jonva,jonva))=dpyq(jonva)/yionva
                  if(jonva.gt.0) then
                     d2pyq(kxsym(jonva,jonva))=dpyq(jonva)/yionva
                  endif
!                  write(*,216)'3X all dpyq:',gz%intlevel,ic,dpyq
!                  write(*,216)'3X all d2pyq:',gz%intlevel,ic,d2pyq
! END SPECIAL FOR IONIC LIQUID
!---------------------------------------------------------------------
               endif cationintandva
! we must check if any endmember is wildcard like L(phase,*:A,B)
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
213               format(a,10i5)
                  dpyq(ic)=pyq*(one-gz%yfrint(gz%intlevel))
               else ! not a wildcard
! this is the normal first derivative of pyq*y(ic) with respect to y(ic)=ymult
                  dpyq(ic)=pyq
                  if(ionicliq) then
!                        write(*,214)'3X Multiply with y_va: ',&
!                             iliqsave,ic,intlat,yionva,pyq
214                  format(a,l2,2i3,4(1pe12.4))
                     if(iliqsave .and. intlat.eq.1.and.yionva.gt.zero) then
! for compatibility with substitutional liquids, multiply interactions 
! of cations (in 1st subl) when vacancies in 2nd with the vacancy fraction
                        dpyq(ic)=pyq*yionva
                     endif
                     endif
                  endif wildcard
!                  write(*,228)'3X dpyq: ',(dpyq(ll),ll=1,4)
228               format(a,6(1pe12.4))
! pyq calculated identically for wildcards as ymult set differently above
! It should work for ionic liquids as ymult has been multiplied with yionva
                  pyq=pyq*ymult
                  proprec=>intrec%propointer
!               write(*,218)'3X pyq: ',associated(proprec),ymult,pyq
218               format(a,l2,2(1pe12.4))
! list values of pyq, dpyg, d2pyg
!               write(*,228)'3X pyq:',pyq
!               write(*,228)'3X dpy:',dpyq
!               write(*,228)'3X d2py:',d2pyq
219               format(a,6(1pe12.4))
!..............................
! Here we finally calculate the interaction parameter .... SUCK
               intprop: do while(associated(proprec))
! calculate interaction parameter, can depend on composition
! maybe faster to zero here than inside cgint ??
                  vals=zero
                  dvals=zero
                  d2vals=zero
!                  call cgint(lokph,proprec,moded,vals,dvals,d2vals,gz,ceq)
                  call cgint(lokph,proprec,moded,&
                       vals,dvals,d2vals,gz,tooprec,ceq)
                  if(gx%bmperr.ne.0) goto 1000
!                  write(*,228)'3X val:',vals(1),(dvals(1,id),id=1,gz%nofc)
! G parameters (ipy=1) are divided by RT inside cgint
                  typty=proprec%proptype
                  if(typty.ne.1) then
! check if magnetic and wildcard ...
                     if(wildmob) then
! nowildcard(1..3) set in gtp_init in gtp3A.F90 for mobility parameters
! typty is indicator*100 + constituent index
                        do qz=1,3
                           if(typty/100.eq.nowildcard(qz)) then
                              write(*,*)&
                                   '3X mobilities must not have wildcards',lokph
                              gx%bmperr=4374; goto 1000
                           endif
                        enddo
                     endif
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
! Testing using jxsym ... OK here also
                           jxsym=kxsym(id,id)
!                           iloop4: do jd=id+1,gz%nofc
! This loop was constructed for normal cases when pyq has each fraction once
! in ionic liquids Va can have a power so loop for all!
                           iloop4: do jd=id,gz%nofc
!                              phres%d2gval(ixsym(id,jd),ipy)= &
!                                   phres%d2gval(ixsym(id,jd),ipy)+ &
!                                   d2pyq(ixsym(id,jd))*vals(1)
                              if(ixsym(id,jd).ne.jxsym) then
                                 write(*,*)'ISYM error 2',id,jd,&
                                      ixsym(id,jd),jxsym
                                 stop
                              endif
                              phres%d2gval(jxsym,ipy)= &
                                   phres%d2gval(jxsym,ipy)+ &
                                   d2pyq(jxsym)*vals(1)
                              jxsym=jxsym+jd
!                              write(*,251)'3X G:',id,jd,ixsym(id,jd),&
!                                   d2pyq(ixsym(id,jd)),vals(1)
251                           format(a,3i3,4(1pe12.4))
                           enddo iloop4
                        endif
     toop7: if(associated(tooprec)) then
! this is part of iloop3 is for all components "id", starting 30 lines above  
! Normally binay interactions depend only on the constituents gz%iq(1) 
! and gz%iq(2) but Toop/Kohler method depend also on other constituents! 
! Thus dvals may not be correctly updated but that is taken care here
!        write(*,'(a,i3,5(1pe12.4))')'3X toop7: ',&
!             id,pyq,dvals(1,id),phres%dgval(1,id,ipy)+pyq*dvals(1,id),&
!             (phres%dgval(1,id,ipy)-pyq*dvals(1,id))*gz%rgast
! For those who forget dgval1,j,1) is derivative of G wrt constituent j
! dgval(2,j,1) is derivative of G wrt constituent j and T
! dgval(3,j,1) is derivative of G wrt constituent j and P. 
! Third index is for other properties
        do itp=1,3
           phres%dgval(itp,id,ipy)=phres%dgval(itp,id,ipy)-pyq*dvals(itp,id)
        enddo
! I AM NOT SURE THIS IS CORRECT/BoS
! ignore contribution to the second derivatives phres%d2gval
! iloop3  ends just a few lines below
     endif toop7
                        do itp=1,3
                           phres%dgval(itp,id,ipy)=&
                                phres%dgval(itp,id,ipy)+dpyq(id)*vals(itp)
                        enddo
                     enddo iloop3
!                     write(*,211)'3X Interactions: ',gz%iq,jonva
211                  format(a,5i3,5x,i3)
!                     if(jonva.gt.0) then
!                        write(*,212)jonva,phres%dgval(1,jonva,1)*rtg
212                     format('3X with va: ',i3,6(1pe12.4))
!                     endif
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
! no jxsym here ... to complicated
                           do qz=jk,4
!                              phres%d2gval(ixsym(gz%iq(jk),gz%iq(qz)),ipy)=&
!                                 phres%d2gval(ixsym(gz%iq(jk),gz%iq(qz)),ipy)+&
!                                  dpyq(gz%iq(jk))*dvals(1,gz%iq(qz))+&
!                                  dpyq(gz%iq(qz))*dvals(1,gz%iq(jk))
! I do not trust optimized gfortran will eliminate 2 calls to ixsym !!!
                              jxsym=ixsym(gz%iq(jk),gz%iq(qz))
                              phres%d2gval(jxsym,ipy)=&
                                   phres%d2gval(jxsym,ipy)+&
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
!                           phres%d2gval(ixsym(gz%iq(jk),gz%iq(qz)),ipy)=&
!                                phres%d2gval(ixsym(gz%iq(jk),gz%iq(qz)),ipy)+&
!                                dpyq(gz%iq(jk))*dvals(1,gz%iq(qz))+&
!                                dpyq(gz%iq(qz))*dvals(1,gz%iq(jk))
! not trusting gfortran optimizing
                           jxsym=ixsym(gz%iq(jk),gz%iq(qz))
                           phres%d2gval(jxsym,ipy)=&
                                phres%d2gval(jxsym,ipy)+&
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
!                     phres%d2gval(ixsym(gz%iq(jk),gz%iq(jk)),ipy)=&
!                          phres%d2gval(ixsym(gz%iq(jk),gz%iq(jk)),ipy)+&
!                          2.0D0*dpyq(gz%iq(jk))*dvals(1,gz%iq(jk))
! not trusing gforntran optimizing
                     jxsym=ixsym(gz%iq(jk),gz%iq(jk))
                     phres%d2gval(jxsym,ipy)=&
                          phres%d2gval(jxsym,ipy)+&
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
!                              add1=dpyq(ic1)*dvals(1,gz%iq(1))+&
!                                   dpyq(gz%iq(1))*dvals(1,ic1)+&
!                                   pyq*d2vals(ixsym(ic1,gz%iq(1)))
!                              phres%d2gval(ixsym(ic1,gz%iq(1)),ipy)=&
!                                   phres%d2gval(ixsym(ic1,gz%iq(1)),ipy)+add1
! not trusing gfortran optimizing
                              jxsym=ixsym(ic1,gz%iq(1))
                              add1=dpyq(ic1)*dvals(1,gz%iq(1))+&
                                   dpyq(gz%iq(1))*dvals(1,ic1)+&
                                   pyq*d2vals(jxsym)
                              phres%d2gval(jxsym,ipy)=&
                                   phres%d2gval(jxsym,ipy)+add1
                              if(ic1.ne.gz%iq(1)) then
! this IF to avoid that the second derivative gz%iq(1) and gz%iq(2) is
! calculated twice. ic1 will at some time be equal to gz%iq(1) and to gz%iq(2)
!                                 add1=dpyq(ic1)*dvals(1,gz%iq(2))+&
!                                      dpyq(gz%iq(2))*dvals(1,ic1)+&
!                                      pyq*d2vals(ixsym(ic1,gz%iq(2)))
!                                 phres%d2gval(ixsym(ic1,gz%iq(2)),ipy)=add1+&
!                                      phres%d2gval(ixsym(ic1,gz%iq(2)),ipy)
! not trusting gfortran optimizing
                                 jxsym=ixsym(ic1,gz%iq(2))
                                 add1=dpyq(ic1)*dvals(1,gz%iq(2))+&
                                      dpyq(gz%iq(2))*dvals(1,ic1)+&
                                      pyq*d2vals(jxsym)
                                 phres%d2gval(jxsym,ipy)=add1+&
                                      phres%d2gval(jxsym,ipy)
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
! to many indentations--------------------------------------------
         catx: if(ionicliq) then
! for ionic liquid when interactions involve cations there is a contribution
! due to the vacancy fraction multiplied with the cations yc1*yc2*yva**2
! we are dealing with binary RK interactions, gz%intlevel=1, check if 
! interaction is in first sublattice (between cations) and vacancy in second
            if(iliqva .and. jonva.gt.0) then
               if(gz%intlat(1).eq.1) then
! add pyq multipled with the derivative with respect to vacancy fraction
! This should be done for d2gval also but I skip that at present ...
                  phres%dgval(itp,jonva,ipy)=&
                       phres%dgval(itp,jonva,ipy)+pyq*dvals(itp,jonva)
!                       write(*,*)'3X jonva:',jonva,pyq,dvals(1,jonva)
               elseif(gz%intlat(1).eq.2 .and. gz%iq(2).gt.jonva) then
! This fixed the problem with Pd-Ru-Te in the fuel (+ fix in cgint)
!                write(*,55)'3X (C:Va,K)',iliqva,gz%intlat(1),jonva,gz%iq(1),&
!                     gz%iq(2),gz%endcon(1),pyq,dvals(itp,gz%endcon(1))
55                format(a,l2,5i3,4(1pe12.4))
               icat=gz%endcon(1)
               if(icat.gt.0) then
                  phres%dgval(itp,icat,ipy)=&
                    phres%dgval(itp,icat,ipy)+pyq*dvals(itp,icat)
!               else
! wow, icat can be -99 meaning interaction between neutrals ....
! but then just skip as the assignment above is not relevant
! Error occured for ionic liquid with:
! 1    2    3   4    5    6    7    8    9    10  11
! BA+2 CE+3 CS+ GD+3 LA+3 MO+4 PD+2 PU+3 RU+4 U+4 ZR+4 :               
! I- MOO4-2 O-2 VA CEO2 CS2TE CSO2 I2 MOO3 O  PUO2 TE TEO2
! 12 13     14  15 16   17    18   19 20   21 22   23 24
! we come here with: BA+2:VA,TE; PD+2:VA,TE; RU+4:VA,TE; *:CS2TE,TE
! gz%intlat(1)=2 OK; jonva=15 OK; gz%iq(1)=17(Cs2Te); gz%iq(2)=23(Te); 
! gz%endcon(1)=-99
!                  write(*,55)'3X Instroini:',iliqva,gz%intlat(1),jonva,&
!                       gz%iq(1),gz%iq(2),gz%endcon(1),pyq
               endif
            endif
         endif
      endif catx
! increase indendation-----------------------------------------
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
! here we could maybe use cycle interloop ??/Bosse 2023.12.16
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
303   format(a,3(1x,l1),i3,4(1pe12.4))
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
         returnoradd: if(first) then
! we have calculated for the first, now calculate for second fraction type
! alternative method: no need to calculate with all fractions as disordered
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
!               allocate(tmpd2g(nz*(nz+1)/2,nprop))
!               tmpd2g=zero
!--------------------------------------------------------------------------
! simplest way of correcting 2nd deruvatives, Gord(y=x) in phres%d2gval
! phres%d2gval(i,j) = saved2g(i,j) - phres%d2gval(i,j)
               do ipy=1,lprop-1
                  do i1=1,gz%nofc
!                     jxsym=ixsym(i1,i1)
                     jxsym=kxsym(i1,i1)
! It should work with jxsym here 
                     do i2=i1,gz%nofc
!                        if(ixsym(i1,i2).ne.jxsym) then
! this ixsym test works and has run of few 1000 times, removed for speed!!
!                           write(*,*)'ISYM error 3',i1,i2,ixsym(i1,i2),jxsym
!                           stop
!                        endif
                        phres%d2gval(jxsym,ipy)=&
                             saved2g(jxsym,ipy)-&
                             phres%d2gval(jxsym,ipy)
! adding i2 to jxsym here seems correct!!
                        jxsym=jxsym+i2
!                        phres%d2gval(ixsym(i1,i2),ipy)=&
!                             saved2g(ixsym(i1,i2),ipy)-&
!                             phres%d2gval(ixsym(i1,i2),ipy)
                     enddo
                  enddo
               enddo
!               goto 667
! old code removed
!667            continue
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
         endif returnoradd
! code above reinstated but has problems ....
      endif disord
! WE CAN JUMP HERE WITHOUT CALCULATING THE ORDERED PART AS DISORDERED
400   continue
!      write(*,*)'3X calcg_internal at label 400'
   enddo fractyp
!   norfc=phlista(lokph)%tnooffr
! 4SL FCC all correct here
!   write(*,69)'3X d2G/dy2B:',norfc,(phres%d2gval(ixsym(j1,j1),1),j1=1,norfc)
69 format(a,i3,6(1pe12.4))
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
! problem using jxsym here, map13 crashed FCC 4 sublattice orering!!!
! PAY ATTENTION TO indices!! we have both i1, i2 and j1, j2
!               jxsym=ixsym(i1,i1)
               jxsym=kxsym(i1,i1)
               do i2=i1,norfc
! add the contributions from the disordered part
                  j2=fracset%y2x(i2)
!                  if(ixsym(i1,i2).ne.jxsym) then
! this ixsym test works and has run of few 1000 times, removed for speed!!
!                     write(*,*)'ISYM error 4',i1,i2,ixsym(i1,i2),jxsym
!                     stop
!                  endif
                  do ipy=1,lprop-1
                     phpart%d2gval(jxsym,ipy)=&
                          phpart%d2gval(jxsym,ipy)+&
                          fsites*phres%d2gval(ixsym(j1,j2),ipy)*&
                          fracset%dxidyj(i1)*fracset%dxidyj(i2)
!                     phpart%d2gval(ixsym(i1,i2),ipy)=&
!                          phpart%d2gval(ixsym(i1,i2),ipy)+&
!                          fsites*phres%d2gval(ixsym(j1,j2),ipy)*&
!                          fracset%dxidyj(i1)*fracset%dxidyj(i2)
                  enddo
                  jxsym=jxsym+i2
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
!         jxsym=ixsym(i1,i1)
         jxsym=kxsym(i1,i1)
         secondd: do i2=i1,norfc
            do ipy=1,lprop-1
!               write(*,497)'3X adding: ',i1,i2,ixsym(i1,i2),ipy
497            format(a,10i3)
!               if(ixsym(i1,i2).ne.jxsym) then
! this ixsym test works and has run of few 1000 times, removed for speed!!
!                  write(*,*)'ISYM error 5',i1,i2,ixsym(i1,i2),jxsym
!                  stop
!               endif
               phres%d2gval(jxsym,ipy)=saved2g(jxsym,ipy)+&
                    phres%sites(2)*phres%d2gval(jxsym,ipy)
!               phres%d2gval(ixsym(i1,i2),ipy)=saved2g(ixsym(i1,i2),ipy)+&
!                    phres%sites(2)*phres%d2gval(ixsym(i1,i2),ipy)
               add1=zero
! IMPORTANT note dpqdy(i1) the the charge of iq, do not confuse with dpyq ...
               if(i1.le.phlista(lokph)%nooffr(1)) then
                  add1=phres%dpqdy(i1)*phres%dgval(1,i2,ipy)
               endif
               if(i2.le.phlista(lokph)%nooffr(1)) then
                  add1=add1+phres%dpqdy(i2)*phres%dgval(1,i1,ipy)
               endif
               phres%d2gval(jxsym,ipy)=phres%d2gval(jxsym,ipy)+add1
!               phres%d2gval(ixsym(i1,i2),ipy)=phres%d2gval(ixsym(i1,i2),ipy)+&
!                    add1
               jxsym=jxsym+i2
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
! The label here just a label, there is no explict jump here
500 continue
!   write(*,69)'3Xa d2G/dy2C:',norfc,(phres%d2gval(ixsym(j1,j1),1),j1=1,norfc)
   if(btest(phmain%status2,CSADDG)) then
! we have an constant addition to G, at present just a constant /RT
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
! doubting gfortran optimizer ...
            jxsym=kxsym(id,jd)
!            phmain%d2gval(ixsym(id,jd),1)=phmain%d2gval(ixsym(id,jd),1)+xxx
            phmain%d2gval(jxsym,1)=phmain%d2gval(jxsym,1)+xxx
         enddo
      enddo
   endif
! uniquac model
   uniquac: if(btest(phlista(lokph)%status1,phuniquac)) then
!      write(*,'(a,6(1pe12.4))')'3X calling uniquac: ',&
!           phmain%dgval(1,1,1),phres%dgval(1,2,1)
      call uniquac_model(moded,gz%nofc,phmain,ceq)
      if(gx%bmperr.ne.0) goto 1000
   endif uniquac
!................................
! calculate additions like magnetic contributions etc and add to G
! Now also Einstein, twostate liquid, volume ...
! if liq2state is FALSE we should add that constribution
! using composition dependent G2 parameters
   addrec=>phlista(lokph)%additions
!   write(*,*)'3X check for first addrec: ',associated(addrec)
   additions: do while(associated(addrec))
! Note for phases with a disordered fraction set, gz%nofc is equal to
! the disordered number of fractions here 
      gz%nofc=phlista(lokph)%tnooffr
! moded is 0, 1 or 2 if derivatives should be calculated, phres is pointer
! to result arrays, lokadd is the addition record, listprop is needed to
! find where TC and BM are stored, gz%nofc are number of constituents
! EINSTEIN
!      write(*,*)'3X addition select: ',phres%gval(1,2),gz%nofc
!      write(*,1001)'Addto: ',gx%bmperr,(phres%gval(j1,1),j1=1,4)
      call addition_selector(addrec,moded,phres,lokph,gz%nofc,ceq)
      if(gx%bmperr.ne.0) goto 1000
! NOTE that the addition record is not in the dynamic data structure
! but the values calculated are returned added to phres which is dynamic
! There is a temporary storage of results for listing only.
      addrec=>addrec%nextadd
!      write(*,*)'3X check for next addrec: ',associated(addrec)
   enddo additions
! there are some special properties like mobilities and similar which
! have a conmponent or constituent index like MQ&<constituent>
!   ipy=typty/100+mod(typty,100)
!   if(ipy.gt.10) then
!      write(*,*)'3X Property ',typty,ipy
!   write(*,*)'3X extra 2: ',phres%gval(1,2)
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
! 4SL all correct here also!
!   write(*,69)'3Xb d2G/dy2F:',norfc,(phres%d2gval(ixsym(j1,j1),1),j1=1,norfc)
! running out of memory??
! these are locally allocated, should be deallocated automatically
! Segmentation fault if I do not write ... but ... reason somewhere else
!   write(*,*)'3X deallocate dpyq?',allocated(dpyq)
   if(allocated(dpyq)) deallocate(dpyq)
!   write(*,*)'3X deallocate d2pyq?',allocated(d2pyq)
   if(allocated(d2pyq)) deallocate(d2pyq)
!   write(*,*)'3X deallocate dvals?',allocated(dvals)
   if(allocated(dvals)) deallocate(dvals)
!   write(*,*)'3X deallocate d2vals?',allocated(d2vals)
   if(allocated(d2vals)) deallocate(d2vals)
!   write(*,*)'3X calcg_internal all deallocated'
!   if(size(phres%yfr).gt.2) then
! debug cqc:
!      write(*,480)'3X dg/dt/RT: 2: ',qcmodel,phres%yfr(3),&
!           phres%gval(1,1),phres%gval(2,1)
!   endif
!   write(*,1001)'Total: ',gx%bmperr,(phres%gval(j1,1),j1=1,4)
!    write(*,1002)(phres%dgval(1,i,1),i=1,3)
!    write(*,1003)(phres%d2gval(i,1),i=1,6)
1001 format('3X/',a,i5,4(1PE12.4))
1002 format('3X calcg dg:  ',3(1PE15.7))
1003 format('3X calcg d2g: ',6(1PE11.3))
   return
 end subroutine calcg_internal

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

!\addtotable subroutine cgint
!\begin{verbatim}
 subroutine cgint(lokph,lokpty,moded,vals,dvals,d2vals,gz,tooprec,ceq)
! calculates an excess parameter that can be composition dependent
! gz%yfrem are the site fractions in the end member record
! gz%yfrint are the site fractions in the interaction record(s)
! lokpty is the property index, lokph is the phase record
! vals, dvals, d2vals multiplied by endmember and interaction fractions outside
! moded=0 means only G, =1 G and dG/dy, =2 all
   implicit none
   integer moded,lokph
   TYPE(gtp_property), pointer :: lokpty
   TYPE(gtp_parcalc) :: gz
   double precision vals(6),dvals(3,gz%nofc)
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_tooprec), pointer :: tooprec
!\end{verbatim}
! temporary data like gz%intlevel, gz%nofc etc
   double precision d2vals(gz%nofc*(gz%nofc+1)/2),valtp(6)
   double precision vv(0:2),fvv(0:2)
   integer lfun,jdeg,jint,qz,ivax,icat
   double precision rtg,dx0,dx,dx1,dx2,ct,fvs,dvax0,dvax1,dvax2,yionva
   double precision ycat0,dcat1,dcat2,dyvan1,dyvan2
   double precision, parameter :: onethird=one/3.0D0,two=2.0D0
   logical ionicliq,iliqva,iliqneut,iliq3cat
! zeroing 5 iq, and vals, dvals and d2vals
!   write(*,*)'3X cgint 1:',gz%iq(1),gz%iq(2),gz%iq(3)
! why zero qz%iq, it has been set before calling ...
   gz%iq=0
! these 3 arrays are in the call and may already have some stored values 
!   vals=zero
!   dvals=zero
!   d2vals=zero
!-------------
   rtg=gz%rgast
! to avoid warnings from -Wmaybe-uninitiated
   icat=0
   ivax=0
   dvax0=zero
   dvax1=zero
!   write(*,*)'3X in cgint',lokph
   if(lokpty%degree.eq.0) then
!----------------------------------------------------------------------
! Easy: no composition dependence.  This applies also to Toop/Kohler parameters
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
! plain binary Redlich Kister or Toop/Kohler method
! gz%endcon can be wildcard, i.e. negative
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
      if(associated(tooprec)) then
! This is a Kohler-Toop method parameter
! only for binary interaction parameters with Kohler or Toop models
! not composition dependent we never come here as we exit 50 lines above
         call calc_toop(lokph,lokpty,moded,vals,dvals,d2vals,gz,tooprec,ceq)
! we have calculated all, skip the rest of this subroutine
         goto 1000
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
! we must know the cation (if only neutrals set to one)
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
! normal CEF model.  Note negative sign taken care of when "added"
            dx2=(jdeg+1)*dx1
            dx1=(jdeg+1)*dx
            dx=dx*dx0
         endif
      enddo RK
   elseif(gz%intlevel.eq.2) then !intlev
!----------------------------------------------------------------------
! important to set ivax=0 here as tested below if not zero
      ivax=0
      iliq3cat=.FALSE.
! it can be a ternary interaction in same sublattice or a reciprocal parameter
!      write(*,*)'3X gz%intlat: ',gz%intlat
!      write(*,*)'3X gz%intcon: ',gz%intcon
      if(ionicliq) then
!         write(*,*)'3X Comp.dep ternary ionic liquid parameter: ',iliqva
         if(gz%intlat(1).eq.2) then
! Both interacting constituents in second sublattice, this should handle these:
! TAFID problem: (5):(33,37,56) is (CA+2):(ALO2-,SIO4-4,SIO2) !!!
! TAFID problem: (12):(33,37,56) is (MG+2):(ALO2-,SIO4-4,SIO2) !!!
! TAFID problem: (9):(38,39,41)  is (Fe+2):(VA,B,C) !!
! not tested: (Fe+2):(S-2,Va,S) or similar ... but it should be OK
            continue
         elseif(iliqva) then
! the pair constituent in second sublattice is Va, no anions!!
            if(gz%intlat(1).eq.1 .and. gz%intlat(2).eq.1) then
! we have 3 cations interacting in first sublattice and Va in second
! with composition dependence .... require treatment of extra vacancy fraction
! TAFID not implemented: (Fe+2,Cr+2,Ni+1):(Va) for example ....
! ternary term: y_Va*(y_Cr*L;0 +y_Fe*L;1 +y_Ni*L;2)
!               write(*,*)'3X unimplemented comp. dep. ternary cation',&
!                    'interaction in liquid'
!               gx%bmperr=4343; goto 1000
! 
               iliq3cat=.TRUE.
               ivax=gz%intcon(2)
            elseif(gz%intlat(1).eq.1 .and. gz%intlat(2).eq.2) then
! This is a reciprocal interaction, two cations, vacancy and neutral
               ivax=gz%endcon(2)
            endif
         elseif(gz%intcon(2).eq.phlista(lokph)%i2slx(1)) then
! reciprocal interaction between two cations in first sublattice and
! an anion and vacancy is the second sublattice
!            write(*,*)'3X reciprocal with 2 cations and anion and Va'
            ivax=gz%intcon(2)
            yionva=gz%yfrint(2)
         else
! I do not know what kind of parameter this is
! THIS ERROR OCCUR ONLY IN PARALLEL
!            write(*,28)'3X: unknown I2SL parameter on level: ',&
!                 gz%intlevel,gz%endcon(1),gz%intcon(1),gz%intcon(2),&
!                 gz%endcon(2),gz%iq,iliqva
28          format(a,i2,': ',i2,',',i2,',',i2,':',i2,5x,5i3,2x,l2)
!            gx%bmperr=4342; goto 1000
            goto 1000
         endif
! other ternary parameters in ionic liquid OK, no extra vacancy fraction
      endif
!................................................................
300 continue
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
! derivatives of vv w.r.t. the 3 constituents 0, 1 and 2
         fvv(0)=two*onethird
         fvv(1)=-onethird
         fvv(2)=-onethird
         if(size(lokpty%degreelink).eq.2) then
! KRASCH if only two degrees of ternary parameter (3 MUST BE GIVEN)
! If only one it is not composition dependent!
            write(*,37)trim(phlista(lokph)%name)
37          format('3X Database error, ternary composition dependent',&
                 ' parameter in ',a/'must have 3 degrees.')
            gx%bmperr=4342; goto 1000
         endif
         terloop: do jint=0,2
! calculate parameters, there are 3 of them, jint=0, 1 and 2
            lfun=lokpty%degreelink(jint)
            call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
            if(lokpty%proptype.eq.1) then
               valtp=valtp/rtg
            endif
! function value
            if(iliq3cat) then
! this is when there are 3 cations in ionic liquid, yionva is vacancy fraction
! NOTE vals and valtp both have dimension 6!!
               vals=vals+yionva*vv(jint)*valtp
! there are also a contrbution to df/dyva!, d2f/dyvadT ... calculated below
! ivax is the index of vacancy
!               write(*,*)'3X ivax: ',ivax
               dvals(1,ivax)=dvals(1,ivax)+vv(jint)*valtp(1)
               dvals(2,ivax)=dvals(2,ivax)+vv(jint)*valtp(2)
               dvals(3,ivax)=dvals(3,ivax)+vv(jint)*valtp(3)
            else
               vals=vals+vv(jint)*valtp
            endif
            noder6: if(moded.gt.0) then
! first derivatives, qz=2 is for T and qz=3 is for P derivatives
               do qz=1,3
! for interaction with 3 cations and Va in 2nd sublattice
! valtp(1) is G; valtp(2) is dG/dT; valtp(3) is dG/dP
                  if(iliq3cat) then
! the first derivatives
                     dvals(qz,gz%iq(1))=dvals(qz,gz%iq(1))+&
                          yionva*fvv(0)*valtp(qz)
                     dvals(qz,gz%iq(2))=dvals(qz,gz%iq(2))+&
                          yionva*fvv(1)*valtp(qz)
                     dvals(qz,gz%iq(3))=dvals(qz,gz%iq(3))+&
                          yionva*fvv(2)*valtp(qz)
                  else
                     dvals(qz,gz%iq(1))=dvals(qz,gz%iq(1))+fvv(0)*valtp(qz)
                     dvals(qz,gz%iq(2))=dvals(qz,gz%iq(2))+fvv(1)*valtp(qz)
                     dvals(qz,gz%iq(3))=dvals(qz,gz%iq(3))+fvv(2)*valtp(qz)
                  endif
               enddo
            endif noder6
            if(iliq3cat) then
! with ionic liquid and 3 cations iteraction there are 2nd derivatives
! with respect to Va and the cation (but no T or P derivative)!
! gz%iq(1) is
               d2vals(ixsym(ivax,gz%iq(jint+1)))=&
                    d2vals(ixsym(ivax,gz%iq(jint+1)))+fvv(0)*valtp(1)
            endif
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

!\addtotable subroutine config_entropy
!\begin{verbatim}
 subroutine config_entropy(moded,nsl,nkl,phvar,tval)
! calculates CEF configurational entropy/R for phase lokph
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
! d2dval(ixsym(N*(M+1)/2),1) are derivatives of G wrt fractions N and M
! this is a symmetric matrix and index givem by ixsym(M,N)
         ss=ss+yfra*ylog
         if(moded.gt.0) then
            phvar%dgval(1,kall,1)=phvar%sites(ll)*(one+ylog)
!            phvar%d2gval(ixsym(kall,kall),1)=phvar%sites(ll)/yfra
! kxsym same as ixsym when first index is >= second index
            phvar%d2gval(kxsym(kall,kall),1)=phvar%sites(ll)/yfra
         endif
      enddo fractionloop
      phvar%gval(1,1)=phvar%gval(1,1)+phvar%sites(ll)*ss
   enddo sublatticeloop
! looking for error calculating 4 sublattice ordered FCC
!   write(*,69)kall,(phvar%d2gval(ixsym(jl,jl),1),jl=1,kall)
69 format('3X d2G/dy2: ',i3,6(1pe12.4))
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

!\addtotable subroutine config_entropy_ssro
!\begin{verbatim}
 subroutine config_entropy_ssro(moded,lokph,phvar,tval)
! test calculates SSRO configurational entropy/R for phase lokph
   implicit none
   integer moded,lokph
   double precision tval
   TYPE(gtp_phase_varres), pointer :: phvar
!\end{verbatim}
! CVM tetrahedron model with ONLY SRO, no sublattices
   integer ia,ib,ic,id,jj,jk,jl,jm,ne,nc,zz,pz
! this is multiplicities
   integer, allocatable :: mijkl(:)
! this is mole fractions and derivatives
   double precision, allocatable :: xf(:),dxf(:,:),d2xf(:,:)
! these are constituent fractions for entropy!!!
   double precision, allocatable :: ysro(:)
! these are site fractions, not needed when only SRO as y^s_i=x_i
!   double precision, allocatable :: yf1(:),yf2(:),yf3(:),yf4(:)
!   double precision, allocatable :: dyf1(:,:),dyf2(:,:),dyf3(:,:),yf4(:,:)
! these are pair fractions, same for all bonds!
!   double precision, allocatable :: p12(:,:),p13(:,:),p14(:,:),&
!        p23(:,:),p24(:,:),p34(:,:)
!   double precision, allocatable :: dp12(:,:,:),dp13(:,:,:),dp14(:,:,:),&
!        dp23(:,:,:),dp24(:,:,:),p34(:,:,:)
   double precision, allocatable :: pstij(:,:)
   double precision, allocatable :: dpstij(:,:,:),d2pstij(:,:,:)
   double precision, parameter :: f1=0.75D0, f2=0.5D0, f3=0.25D0
!   double precision, parameter :: f1=3.0D0, f2=2.0D0, f3=1.0D0
! auxilliary
   character dummy*80
   double precision pijx,fpf,ssumy,ssump,ssumx,ylog,yfra,ycorr,yrest
   double precision pstijtest,pstijsave
! debugging
   double precision ylog1(5),ylog2,ylog3
   double precision sylog1,sylog2,sylog3
! These factors should be 2.0, -6.0 and 5.0 according to Kikuchi
! I have no idea why I have to divide them by 10 ...
   double precision, parameter :: syfact=2.0D0, spfact=-6.0D0, sxfact=5.0D0
   logical sdebug
   save pstijsave
! use phvar%volatile to initiate ycorr at first itertation
!   sdebug=.TRUE.
   sdebug=.FALSE.
!   write(*,*)'3X sdebug: ',sdebug
   if(phvar%volatile.eq.0) then
! phvar%volatile is set to zero in matsmin: meq_calceq at first iteration
! decrease ycorr when pstijsave constant
      phvar%volatile=phvar%volatile+1
      ycorr=0.5D0
      pstijsave=0.0D0
   endif
   yrest=1.0D0-ycorr
! nc is number of constituent, ne is number of elements
   nc=phlista(lokph)%tnooffr
! using empirical rule to calcuöate ne from nc
   select case(nc)
   case default
      write(*,*)'3X SRO number of constituents not implemented',nc
   case(1)
      ne=1
      write(*,*)'3X SRO entropy zero for single element'
      goto 1000
   case(5)
! binary system
      ne=2 
   case(15) 
      ne=3
   case(35) 
      ne=4
   case(70) 
! without merging AAAB etc there would be 625 clusters instead of just 70
      ne=5
   case(126)
      ne=6
   case(210)
      ne=7
   case(330)
! without merging AAAB etc there would be 4096 clusters
      ne=8
   end select
!
!   write(*,*)'3X CVMTFS model for configurational entropy',nc,ne
!
   allocate(mijkl(nc))
   allocate(xf(ne))
   allocate(dxf(ne,nc))
   allocate(d2xf(ne,nc))
!   allocate(yf1(ne))
!   allocate(yf2(ne))
!   allocate(yf3(ne))
!   allocate(yf4(ne))
!   allocate(dyf1(ne,nc))
!   allocate(dyf2(ne,nc))
!   allocate(dyf3(ne,nc))
!   allocate(dyf4(ne,nc))
! jj incremented for each cluster
   jj=0
   mijkl=0
   xf=zero
   dxf=zero
! site fractions same as mole fractions as no LRO
!   yf1=zero
!   yf2=zero
!   yf3=zero
!   yf4=zero
!   dyf1=zero
!   dyf2=zero
!   dyf3=zero
!   dyf4=zero
! extrahera mole fractions from constituent fractions
   do ia=1,ne
      jj=jj+1
! this is AAAA or BBBB etc
      mijkl(jj)=1
      xf(ia)=xf(ia)+phvar%yfr(jj)
      dxf(ia,jj)=1
      do ib=ia+1,ne
         jj=jj+3
         mijkl(jj-2)=4
         mijkl(jj-1)=6
         mijkl(jj)=4
! jj-2 is A3B1, jj-1 is A2B2, jj is A1B3 including permutations in mijkl
         xf(ia)=xf(ia)+f1*phvar%yfr(jj-2)+f2*phvar%yfr(jj-1)+f3*phvar%yfr(jj)
         xf(ib)=xf(ib)+f3*phvar%yfr(jj-2)+f2*phvar%yfr(jj-1)+f1*phvar%yfr(jj)
         dxf(ia,jj-2)=f1; dxf(ia,jj-1)=f2; dxf(ia,jj)=f3;
         dxf(ib,jj-2)=f3; dxf(ib,jj-1)=f2; dxf(ib,jj)=f1;
         do ic=ib+1,ne
            jj=jj+3
            mijkl(jj-2)=12
            mijkl(jj-1)=12
            mijkl(jj)=12
! jj-2 is A2BC, jj-1 is AB2C, jj is ABC2
            xf(ia)=xf(ia)+f2*phvar%yfr(jj-2)+f3*phvar%yfr(jj-1)+f3*phvar%yfr(jj)
            xf(ib)=xf(ib)+f3*phvar%yfr(jj-2)+f2*phvar%yfr(jj-1)+f3*phvar%yfr(jj)
            xf(ic)=xf(ic)+f3*phvar%yfr(jj-2)+f3*phvar%yfr(jj-1)+f2*phvar%yfr(jj)
            dxf(ia,jj-2)=f2; dxf(ia,jj-1)=f3; dxf(ia,jj)=f3;
            dxf(ib,jj-2)=f3; dxf(ib,jj-1)=f2; dxf(ib,jj)=f3;
            dxf(ic,jj-2)=f3; dxf(ic,jj-1)=f3; dxf(ic,jj)=f2;
            do id=ic+1,ne
               jj=jj+1
               mijkl(jj-2)=24
! jj is ABCD
               xf(ia)=xf(ia)+f3*phvar%yfr(jj)
               xf(ib)=xf(ib)+f3*phvar%yfr(jj)
               xf(ic)=xf(ic)+f3*phvar%yfr(jj)
               xf(id)=xf(id)+f3*phvar%yfr(jj)
               dxf(ia,jj)=f3; dxf(ib,jj)=f3; dxf(ic,jj)=f3; dxf(id,jj)=f3;
            enddo
         enddo
      enddo
   enddo
! Convergence problem, when pstij constant make rest approach 1.0
   pstijtest=zero
   do ia=1,ne
      do ib=ia+1,ne
         if(xf(ia)*xf(ib).gt.pstijtest) pstijtest=xf(ia)*xf(ib)
      enddo
   enddo
! pstijtest is maximum pair fraction using mole fractions calculated
!     from cluster fractions provided by the minimizer
! if almost the same as previous decrease ycorr
   if(abs(pstijtest-pstijsave).lt.1.0D-4) then
      ycorr=max(0.5D0*ycorr,1.0D-8)
! when yrest=1 we use the fractions from the minimizer
   else
      ycorr=min(0.5d0*ycorr,0.5D0)
   endif
   yrest=1.0D0-ycorr
   if(SDEBUG) write(*,'(a,1x,l,1x,12F6.3)')'3X corr: ',sdebug,&
        pstijtest,pstijsave,yrest,ycorr
! without this stupid dummy statement the calculations does not converge
   write(dummy,'(a,1x,l,1x,12F6.3)')'3X corr: ',sdebug,&
        pstijtest,pstijsave,yrest,ycorr
! save pstij for next iteration
   pstijsave=pstijtest
! IDEA
! We have no LRO, fractions on all sublattices same and equal to molefractions
! THUS recalculate the cluster fractions from the mole fractions ....
   allocate(ysro(jj))
   jj=0
   do ia=1,ne
      jj=jj+1
      ysro(jj)=yrest*phvar%yfr(jj)+ycorr*xf(ia)**4
      do ib=ia+1,ne
         jj=jj+1
         ysro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ib)*xf(ia)**3
         jj=jj+1
         ysro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ib)**2*xf(ia)**2
         jj=jj+1
         ysro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ib)**3*xf(ia)
         do ic=ib+1,ne
            jj=jj+1
            ysro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ia)**2*xf(ib)*xf(ic)
            jj=jj+1
            ysro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ia)*xf(ib)**2*xf(ic)
            jj=jj+1
            ysro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ia)*xf(ib)*xf(ic)**2
            do id=ic+1,ne
               jj=jj+1
               ysro(jj)=yrest*phvar%yfr(jj)+&
                    ycorr*mijkl(jj)*xf(ia)*xf(ib)*xf(ic)*xf(id)
            enddo
         enddo
      enddo
   enddo
!
   if(SDEBUG) then
      write(*,'(a,12F7.4)')'3X yfr: ',(phvar%yfr(jj),jj=1,nc)
      write(*,'(a,12F6.3)')'3X ysro:',ysro
   endif
!
!   write(*,'(a,F5.2,10F6.3)')'3X xf2: ',fpf,(xf(ia),ia=1,ne)
! calculate pair fractions using the site fractions
!   allocate(p12(ne,ne))
!   allocate(p13(ne,ne))
!   allocate(p14(ne,ne))
!   allocate(p23(ne,ne))
!   allocate(p24(ne,ne))
!   allocate(p34(ne,ne))
   allocate(pstij(ne,ne))
   allocate(dpstij(ne,ne,nc))
   allocate(d2pstij(ne,ne,nc*(nc+1)/2))
   pstij=zero
   dpstij=zero
   d2pstij=zero
! calculation of pair fractions using mole fractions (same as site fractions)
! include AA, AB, AC, BB etc pairs but exclude BA, CA ??
   do ia=1,ne
      do ib=1,ne
! If we had LRO we should need yf1, yf2, yf3 and yf4
! but the site fractions are the same in all sublattces when no LRO (?)
         pstij(ia,ib)=xf(ia)*xf(ib)
         do jj=1,nc
            dpstij(ia,ib,jj)=dxf(ia,jj)*xf(ib)+dxf(ib,jj)*xf(ia)
            do jk=jj,nc
               zz=ixsym(jj,jk)
!               write(*,'(a,2i4,2x,2i4,i7)')'3X d2pstij: ',ia,ib,jj,jk,zz
!               d2pstij(ia,ib,zz)=d2pstij(ia,ib,zz)+&
!                    dxf(ia,jj)*dxf(ib,jk)+dxf(ib,jj)*dxf(ia,jk)
            enddo
         enddo
      enddo
   enddo
! All necessary fractions and derivatives calculated, now the entropy
! S = -2 \sum_ijkl y_ijkl\ln(y_ijkl) +
!     6  \sum_ij p_ij\ln(p_ij) - 5 \sum_j x_j\ln(x_j)
! As we calculate G the signs are inverse
! constituent fraction entropy
! entropy contribution from the constituent fractions ------------------
! USE ysro fractions!!! not phvar%yfr  or a mix ...
   ssumy=zero
!   sylog1=zero
   do jj=1,nc
!      yfra=phvar%yfr(jj)
      yfra=ysro(jj)
      if(yfra.lt.bmpymin) yfra=bmpymin
      if(yfra.gt.one) yfra=one
! yfra is divided by mijkl as it represent mijkl fractions
      ylog=log(yfra/mijkl(jj))
! debugging
!      if(jj.le.5) ylog1(jj)=ylog
!      sylog1=sylog1+yfra*ylog
      ssumy=ssumy+syfact*yfra*ylog
      if(moded.gt.0) then
! dgval(1,1:nc,1) are derivative of G/RT wrt fraction 1:nc
! d2gval(ixsym(jj*(jk+1)/2,1) are 2nd derivarive of G/RT wrt fraction jj and jk
! dgval and d2gval are zero before this loop
!         phvar%dgval(1,jj,1)=syfact*(mijkl(jj)+ylog)
! convergence problem test
         phvar%dgval(1,jj,1)=syfact*(one+ylog)
! ? T and y derivative
         phvar%dgval(2,jj,1)=syfact*(mijkl(jj)+ylog)/tval
! 2nd derivative, each term depend on a single y fraction
         phvar%d2gval(kxsym(jj,jj),1)=syfact*mijkl(jj)/yfra
      endif
!      write(*,'(a,3(1pe12.4))')'3X ssumy: ',yfra,ylog,ssumy
   enddo
   phvar%gval(1,1)=ssumy
   phvar%gval(2,1)=ssumy/tval
! No convergence just using ysro ... change the phvar%yfr ....
   do jj=1,nc
      phvar%yfr(jj)=ysro(jj)
   enddo
! entropy contributions from the 6 pair fractions -----------------------
   ssump=zero
! all pairs the same, no need to loop over sublattices
   fpf=spfact
!   sylog2=zero
   do ia=1,ne
!      do ib=ia+1,ne
      do ib=1,ne
! bond between atom ia and ib
         ylog=log(pstij(ia,ib))
! debugging
!         ylog2=ylog
!         sylog2=sylog2+pstij(ia,ib)*ylog
!         write(*,'(a,2i3,3(1pe12.4))')'3X ylogp: ',ia,ib,pstij(ia,ib),ylog,&
!              pstij(ia,ib)*ylog
         ssump=ssump+fpf*pstij(ia,ib)*ylog
         if(moded.gt.0) then
            do jl=1,nc
! we need derivatives of pair fractions wrt constituent fractions
               phvar%dgval(1,jl,1)=phvar%dgval(1,jl,1)+&
                    fpf*dpstij(ia,ib,jl)*(one+ylog)
               do jm=jl,nc
                  zz=ixsym(jl,jm)
                  phvar%d2gval(zz,1)=phvar%d2gval(zz,1)+&
                       fpf*(d2pstij(ia,ib,zz)*(one+ylog)+&
                       dpstij(ia,ib,jl)*dpstij(ia,ib,jm)/pstij(ia,ib))
               enddo
            enddo
         endif
!         write(*,'(a,2i3,4(1pe12.4))')'3X ssump: ',ia,ib,fpf,pstij(ia,ib),&
!              ylog,ssump
      enddo
   enddo
   phvar%gval(1,1)=phvar%gval(1,1)+ssump
   phvar%gval(2,1)=phvar%gval(2,1)+ssump/tval
! entropy contributions from the mole fractions -----------------------
   ssumx=zero
!   sylog3=zero
   do ia=1,ne
! we need derivatives of mole fractions wrt constituent fractions
      ylog=log(xf(ia))
! debugg
!      ylog3=ylog
!      sylog3=sylog3+xf(ia)*ylog
      ssumx=ssumx+sxfact*xf(ia)*ylog
      if(moded.gt.0) then
         do jl=1,nc
! we need derivatives of mole fractions wrt constituent fractions
            phvar%dgval(1,jl,1)=phvar%dgval(1,jl,1)+&
                 sxfact*dxf(ia,jl)*(ylog+one)
            do jm=jl,nc
! note d2xf/dyidyj==0
               zz=kxsym(jl,jm)
               phvar%d2gval(zz,1)=phvar%d2gval(zz,1)+&
                    sxfact*dxf(ia,jl)*dxf(ia,jm)/xf(ia)
            enddo
         enddo
      endif
   enddo
   phvar%gval(1,1)=phvar%gval(1,1)+ssumx
   phvar%gval(2,1)=phvar%gval(2,1)+ssumx/tval
! All done
   if(SDEBUG) then
      write(*,900)'3X xf: ',(xf(ia),ia=1,ne)
900   format(a,6F7.3)
      write(*,910)'3X pstij: ',((pstij(ia,ib),ib=1,ne),ia=1,ne)
910   format(a,10F7.3)
!      write(*,'(a,5(1pe12.4))')'3X sylog: ',sylog1,sylog2,sylog3
!      write(*,'(a,5(1pe12.4))')'3X ylog1: ',ylog1
!      write(*,'(a,5(1pe12.4))')'3X ylogx: ',ylog2,ylog3
      write(*,'(a,5(1pe12.4))')'3X cvmtfs: ',ssumy,ssump,ssumx,&
           phvar%gval(1,1),8.31451*phvar%gval(1,1)
!   write(*,*)'3X cvmtfs model not yet released'
   endif
1000 continue
   return
 end subroutine config_entropy_ssro

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine config_entropy_cvmtfl
!\begin{verbatim}
 subroutine config_entropy_cvmtfl(moded,lokph,phvar,tval)
! CVM tetrahedron model for FCC with LRO and SRO, all constituents mix
!
   implicit none
   integer moded,lokph
   double precision tval
   TYPE(gtp_phase_varres), pointer :: phvar
!\end{verbatim}
!========================================
! Current code 230308 if same as _SSRO, only SRO ordering
! Modified code 230311 works partially for SRO (not LRO)
!========================================
   integer ia,ib,ic,id,jj,jk,jl,jm,ne,nc,zz,pz,ja,jb
! this is multiplicities
!   integer, allocatable :: mijkl(:)
! this is mole fractions and derivatives
   double precision, allocatable :: xf(:),dxf(:,:),d2xf(:,:)
! these are constituent fractions for entropy!!!
!   double precision, allocatable :: ysro(:)
   double precision, allocatable :: ylro(:)
! these are site fractions, needed for LRO  yf(1..4,*)
   double precision, allocatable :: yf(:,:)
   double precision, allocatable :: dyf(:,:,:)
! these are pair fractions, same for all bonds!
   double precision, allocatable :: pst(:,:,:)
   double precision, allocatable :: dpst(:,:,:,:)
! use this dqst to relate derivatives pf pair directly to cluster fractions
   double precision, allocatable :: dqst(:,:,:,:)
   double precision, allocatable :: d2pst(:,:,:,:)
!   double precision, allocatable :: p12(:,:),p13(:,:),p14(:,:),&
!        p23(:,:),p24(:,:),p34(:,:)
!   double precision, allocatable :: dp12(:,:,:),dp13(:,:,:),dp14(:,:,:),&
!        dp23(:,:,:),dp24(:,:,:),dp34(:,:,:)
!   double precision, allocatable :: d2p12(:,:,:),d2p13(:,:,:),d2p14(:,:,:),&
!        d2p23(:,:,:),d2p24(:,:,:),d2p34(:,:,:)
! to be removed
!   double precision, allocatable :: pstij(:,:)
!   double precision, allocatable :: dpstij(:,:,:),d2pstij(:,:,:)
!   double precision, parameter :: f1=0.75D0, f2=0.5D0, f3=0.25D0
! auxilliary
   character dummy*80
   double precision pijx,ssumy,ssump,ssumx,ylog,yfra,pfra,ycorr,yrest,ybase
   double precision pstijtest,pstijsave
! reduce derivatives
   double precision, parameter :: dfdy=1.0D0
! debugging
   double precision ylog1(5),ylog2,ylog3,mass
   double precision sylog1,sylog2,sylog3
   integer ss,tt,mm,wb1,wb2,wa1,wa2,iphf
   character cluster*4, chel(2)*1, spname*24
! Kikuchi factors, note we calculate dG/dT, not S, thus signs inverted
! spfact -6 not used as we calculate and add 6 pair, sxfact/4 as 4 sublattices
   double precision, parameter :: syfact=2.0D0, spfact=-6.0D0, sxfact=1.250D0
   logical sdebug,s2debug
   save pstijsave
! use phvar%volatile to initiate ycorr at first itertation
!========================================
   write(*,*)
   write(*,*)'3X CVM tetrahedron FCC with LRO for testing only'
!   gx%bmperr=4399
!   goto 1000
!========================================
!   sdebug=.TRUE.
   sdebug=.FALSE.
!   s2debug=.TRUE.
   s2debug=.FALSE.
!   write(*,*)'3X sdebug: ',sdebug
   if(phvar%volatile.eq.0) then
! phvar%volatile is set to zero in matsmin: meq_calceq at first iteration
! decrease ycorr when pstijsave constant
      phvar%volatile=phvar%volatile+1
      ycorr=0.5D0
      pstijsave=0.0D0
   endif
   yrest=1.0D0-ycorr
! nc is number of constituent, ne is number of elements
   nc=phlista(lokph)%tnooffr
! using empirical rule to calcuöate ne from nc
   select case(nc)
   case default
      write(*,*)'3X SRO number of constituents not implemented',nc
   case(1)
      ne=1
      write(*,*)'3X SRO entropy zero for single element'
      goto 1000
   case(16)
! binary system 2*2*2*2=16
      ne=2 
   case(81) 
      write(*,*)'3X max 2 elements!'; gx%bmperr=4399; goto 1000
      ne=3
   case(256) 
      write(*,*)'3X max 2 elements!'; gx%bmperr=4399; goto 1000
      ne=4
!   case(625) 
!      ne=5
   end select
!
!   write(*,*)'3X CVMTFS model for configurational entropy',nc,ne
!
   allocate(xf(ne))
   allocate(dxf(ne,nc))
   allocate(d2xf(ne,nc))
   allocate(yf(4,ne))
   allocate(dyf(4,ne,nc))
   allocate(dqst(6,ne,ne,nc))
! jj incremented for each cluster
   xf=zero
   dxf=zero
! site fractions same as mole fractions as no LRO
   yf=zero
   dyf=zero
! derivative of pair fractions related to clusters
   dqst=zero
! for nice debug output ...
   chel(1)='A'
   chel(2)='B'
! phase number ... 
   iphf=1
! extrahera site fractions from cluster fractions
! order is always AAAA, AAAB_01.._04; AABB_01.._06; ABBB_01.._04, BBBB etc
   jj=0
   ialoop: do ia=1,ne
      jj=jj+1
! this is AAAA or BBBB etc
      do ss=1,4
         yf(ss,ia)=yf(ss,ia)+phvar%yfr(jj)
         dyf(ss,ia,jj)=dfdy
         cluster(ss:ss)=chel(ia)
      enddo
      xf(ia)=xf(ia)+phvar%yfr(jj)
      dxf(ia,jj)=dfdy
!      call get_constituent_name(iphf,jj,spname,mass)
!      write(*,*)'3X cluster: ',cluster,' ',trim(spname),jj
!      write(*,300)(xf(ss),ss=1,ne)
!      write(*,310)((yf(ss,zz),ss=1,4),zz=1,ne)
!      write(*,*)'3X done AAAA: ',jj
      ibloop: do ib=ia+1,ne
! handle 4 of AAAB, 6 of AABB and 4 of ABBB with elements in different order
! In first 4 there is one atom ib in sublattice 4, 3, 2 and 1
         wb1=4
         do jj=jj+1,jj+4
! These are 4 clusters ordered AAAB, AAABA, ABAA, BAAA
!            write(*,*)'3X site fractions in cluster: ',jj,phvar%yfr(jj)
            do ss=1,4
               if(ss.eq.wb1) then
                  yf(ss,ib)=yf(ss,ib)+phvar%yfr(jj)
                  dyf(ss,ib,jj)=dfdy
                  cluster(ss:ss)=chel(ib)
               else
                  yf(ss,ia)=yf(ss,ia)+phvar%yfr(jj)
                  dyf(ss,ia,jj)=dfdy
                  cluster(ss:ss)=chel(ia)
               endif
            enddo
!            call get_constituent_name(iphf,jj,spname,mass)
!            write(*,*)'3X cluster: ',cluster,' ',trim(spname),jj
            wb1=wb1-1
            xf(ia)=xf(ia)+0.75D0*phvar%yfr(jj)
            dxf(ia,jj)=0.75D0*dfdy
            xf(ib)=xf(ib)+0.25D0*phvar%yfr(jj)
            dxf(ib,jj)=0.25D0*dfdy
! derivative of pair fraction relative to cluster fractions, 3 AB bonds
            dqst(1,ia,ib,jj)=3.0d0
         enddo
! after a loop the loop variable is one higher than max limit
         jj=jj-1
!         write(*,300)(xf(ss),ss=1,ne)
!         write(*,310)((yf(ss,zz),ss=1,4),zz=1,ne)
!         write(*,*)'3X done AAAB: ',jj
! there are 3 clusters ordered AABB, ABAB, ABBA, always A in first
         wa1=2
         do jj=jj+1,jj+3
!            write(*,*)'3X site fractions in cluster: ',jj,phvar%yfr(jj)
            yf(1,ia)=yf(1,ia)+phvar%yfr(jj)
            dyf(1,ia,jj)=dfdy
            cluster(1:1)=chel(ia)
            do ss=2,4
               if(ss.eq.wa1) then
                  yf(ss,ia)=yf(ss,ia)+phvar%yfr(jj)
                  dyf(ss,ia,jj)=dfdy
                  cluster(ss:ss)=chel(ia)
               else
                  yf(ss,ib)=yf(ss,ib)+phvar%yfr(jj)
                  dyf(ss,ib,jj)=dfdy
                  cluster(ss:ss)=chel(ib)
               endif
            enddo
!            call get_constituent_name(iphf,jj,spname,mass)
!            write(*,*)'3X cluster: ',cluster,' ',trim(spname),jj
            wa1=wa1+1
            xf(ia)=xf(ia)+0.5D0*phvar%yfr(jj)
            dxf(ia,jj)=0.5d0*dfdy
            xf(ib)=xf(ib)+0.5D0*phvar%yfr(jj)
            dxf(ib,jj)=0.5D0*dfdy
! derivative of pair fraction relative to cluster fractions, 4 AB bonds
            dqst(1,ia,ib,jj)=4.0d0
         enddo
         jj=jj-1
!         write(*,300)(xf(ss),ss=1,ne)
!         write(*,310)((yf(ss,zz),ss=1,4),zz=1,ne)
!         write(*,*)'3X done first half AABB: ',jj
! these are 3 clusters ordered BAAB, BABA, BBAA
         wb1=4
         do jj=jj+1,jj+3
!            write(*,*)'3X site fractions in cluster: ',jj,phvar%yfr(jj)
            yf(1,ib)=yf(1,ib)+phvar%yfr(jj)
            dyf(1,ib,jj)=dfdy
            cluster(1:1)=chel(ib)
            do ss=2,4
               if(ss.eq.wb1) then
                  yf(ss,ib)=yf(ss,ib)+phvar%yfr(jj)
                  dyf(ss,ib,jj)=dfdy
                  cluster(ss:ss)=chel(ib)
               else
                  yf(ss,ia)=yf(ss,ia)+phvar%yfr(jj)
                  dyf(ss,ia,jj)=dfdy
                  cluster(ss:ss)=chel(ia)
               endif
            enddo
!            call get_constituent_name(iphf,jj,spname,mass)
!            write(*,*)'3X cluster: ',cluster,' ',trim(spname),jj
            wb1=wb1-1
            xf(ia)=xf(ia)+0.5D0*phvar%yfr(jj)
            dxf(ia,jj)=0.5D0*dfdy
            xf(ib)=xf(ib)+0.5D0*phvar%yfr(jj)
            dxf(ib,jj)=0.5D0*dfdy
! derivative of pair fraction relative to cluster fractions, 3 AB bonds
            dqst(1,ia,ib,jj)=4.0d0
         enddo
         jj=jj-1
!         write(*,300)(xf(ss),ss=1,ne)
!         write(*,310)((yf(ss,zz),ss=1,4),zz=1,ne)
!         write(*,*)'3X done second half AABB: ',jj
! now 4 clusters ABBB, BABB, BBAB, BBBA
! These are 4 clusters ordered ABBB, BABB, BBAB, BBBA
         wa1=1
         do jj=jj+1,jj+4
!            write(*,*)'3X site fractions in cluster: ',jj,phvar%yfr(jj)
            do ss=1,4
               if(ss.eq.wa1) then
                  yf(ss,ia)=yf(ss,ia)+phvar%yfr(jj)
                  dyf(ss,ia,jj)=dfdy
                  cluster(ss:ss)=chel(ia)
               else
                  yf(ss,ib)=yf(ss,ib)+phvar%yfr(jj)
                  dyf(ss,ib,jj)=dfdy
                  cluster(ss:ss)=chel(ib)
               endif
            enddo
!            call get_constituent_name(iphf,jj,spname,mass)
!            write(*,*)'3X cluster: ',cluster,' ',trim(spname),jj
            wa1=wa1+1
            xf(ia)=xf(ia)+0.25D0*phvar%yfr(jj)
            dxf(ia,jj)=0.25D0*dfdy
            xf(ib)=xf(ib)+0.75D0*phvar%yfr(jj)
            dxf(ib,jj)=0.75D0*dfdy
! derivative of pair fraction relative to cluster fractions, 3 AB bonds
            dqst(1,ia,ib,jj)=3.0d0
         enddo
         jj=jj-1
!         write(*,300)(xf(ss),ss=1,ne)
!         write(*,310)((yf(ss,zz),ss=1,4),zz=1,ne)
!         write(*,*)'3X done ABBB: ',jj
! constitent list
!   3 Q01                       A1A1A1A 
!   4 Q02                       A1A1A1B  
!   5 Q03                       A1A1B1A  
!   6 Q04                       A1B1A1A  
!   7 Q05                       B1A1A1A  
!   8 Q06                       A1A1B1B  
!   9 Q07                       A1B1A1B  
!  10 Q08                       A1B1B1A  
!  11 Q09                       B1A1A1B  
!  12 Q10                       B1A1B1A  
!  13 Q11                       B1B1A1A  
!  14 Q12                       A1B1B1B  
!  15 Q13                       B1A1B1B  
!  16 Q14                       B1B1A1B  
!  17 Q15                       B1B1B1A  
!  18 Q16                       B1B1B1B  
! loop below not active when ne=2 ==============================
         icloop: do ic=ib+1,ne
            write(*,*)'3X LRO not implemented for 3 elements'
            gx%bmperr=4399; goto 1000
!            mijkl(jj-2)=12
!            mijkl(jj-1)=12
!            mijkl(jj)=12
! jj-2 is A2BC, jj-1 is AB2C, jj is ABC2
!           xf(ia)=xf(ia)+f2*phvar%yfr(jj-2)+f3*phvar%yfr(jj-1)+f3*phvar%yfr(jj)
!           xf(ib)=xf(ib)+f3*phvar%yfr(jj-2)+f2*phvar%yfr(jj-1)+f3*phvar%yfr(jj)
!           xf(ic)=xf(ic)+f3*phvar%yfr(jj-2)+f3*phvar%yfr(jj-1)+f2*phvar%yfr(jj)
!           dxf(ia,jj-2)=f2; dxf(ia,jj-1)=f3; dxf(ia,jj)=f3;
!           dxf(ib,jj-2)=f3; dxf(ib,jj-1)=f2; dxf(ib,jj)=f3;
!           dxf(ic,jj-2)=f3; dxf(ic,jj-1)=f3; dxf(ic,jj)=f2;
            idloop: do id=ic+1,ne
               jj=jj+1
!               mijkl(jj-2)=24
! jj is ABCD
!               xf(ia)=xf(ia)+f3*phvar%yfr(jj)
!               xf(ib)=xf(ib)+f3*phvar%yfr(jj)
!               xf(ic)=xf(ic)+f3*phvar%yfr(jj)
!               xf(id)=xf(id)+f3*phvar%yfr(jj)
!               dxf(ia,jj)=f3; dxf(ib,jj)=f3; dxf(ic,jj)=f3; dxf(id,jj)=f3;
            enddo idloop
         enddo icloop
      enddo ibloop
   enddo ialoop
! debug   
!   if(sdebug) then
   write(*,300)(xf(ia),ia=1,ne)
   do ia=1,ne
      write(*,290)chel(ia),(dxf(ia,jj),jj=1,8)
      write(*,291)chel(ia),(dxf(ia,jj),jj=9,16)
   enddo
290 format('3X dxf "',a,'" 1-8 : ',8F6.3)
291 format('3X dxf "',a,'" 9-16: ',8F6.3)
   write(*,310)((yf(ss,ia),ss=1,4),ia=1,ne)
!   if(sdebug) then
      do ia=1,ne
         write(*,270)1,ia,(dyf(1,ia,jj),jj=1,16)
         write(*,270)2,ia,(dyf(2,ia,jj),jj=1,16)
         write(*,270)3,ia,(dyf(3,ia,jj),jj=1,16)
         write(*,270)4,ia,(dyf(4,ia,jj),jj=1,16)
      enddo
270   format('3X dyf: ',2i2,16F4.1)
280   format('3X dyf: ',2i2,8F6.3)
300   format('3X xf: ',8F6.3)
310   format('3X yf: A: ',4F6.3,'  B: ',4F6.3,' C: ',4F6.3)
!   endif
!============================================== pair fractions
   allocate(pst(6,ne,ne))
   allocate(dpst(6,ne,ne,nc))
   allocate(d2pst(6,ne,ne,nc*(nc+1)/2))
   pst=zero
   dpst=zero
   d2pst=zero
   zz=0
   ploop: do ia=1,ne
! in ternary systems pst(ss,ia,ib) can be different from pst(ib,ia) etc
      do ib=1,ne
! assume only SRO, all pair frations the same (in a binary)
!         pst(1,ia,ib)=xf(ia)*xf(ib)
!         pst(2,ia,ib)=pst(1,ia,ib)
!         pst(3,ia,ib)=pst(1,ia,ib)
!         pst(4,ia,ib)=pst(1,ia,ib)
!         pst(5,ia,ib)=pst(1,ia,ib)
!         pst(6,ia,ib)=pst(1,ia,ib)
         pst(1,ia,ib)=yf(1,ia)*yf(2,ib)
         pst(2,ia,ib)=yf(1,ia)*yf(3,ib)
         pst(3,ia,ib)=yf(1,ia)*yf(4,ib)
         pst(4,ia,ib)=yf(2,ia)*yf(3,ib)
         pst(5,ia,ib)=yf(2,ia)*yf(4,ib)
         pst(6,ia,ib)=yf(3,ia)*yf(4,ib)
! Taking the average values here improve convergence.  But not below 0.4 ...
!         p12(ia,ib)=0.5*(yf(1,ia)*yf(2,ib)+yf(2,ia)*yf(1,ib))
!         p13(ia,ib)=0.5*(yf(1,ia)*yf(3,ib)+yf(3,ia)*yf(1,ib))
!         p14(ia,ib)=0.5*(yf(1,ia)*yf(4,ib)+yf(4,ia)*yf(1,ib))
!         p23(ia,ib)=0.5*(yf(2,ia)*yf(3,ib)+yf(3,ia)*yf(2,ib))
!         p24(ia,ib)=0.5*(yf(2,ia)*yf(4,ib)+yf(4,ia)*yf(2,ib))
!         p34(ia,ib)=0.5*(yf(3,ia)*yf(4,ib)+yf(4,ia)*yf(3,ib))
         dploop: do jj=1,nc
! when ignoring LRO
!            dpst(1,ia,ib,jj)=xf(ia)*dxf(ib,jj)+dxf(ia,jj)*dxf(ib,jj)
!            dpst(2,ia,ib,jj)=dpst(1,ia,ib,jj)
!            dpst(3,ia,ib,jj)=dpst(1,ia,ib,jj)
!            dpst(4,ia,ib,jj)=dpst(1,ia,ib,jj)
!            dpst(5,ia,ib,jj)=dpst(1,ia,ib,jj)
!            dpst(6,ia,ib,jj)=dpst(1,ia,ib,jj)
! NOTE:  Make use of the fact that pstij depend direcly on cluster fractions
! Just using relations to y^s_A fauvours z_AAAA and z_BBBB, i.e. no mixing
! IDEA: Use p_AB = 3z_AAAB + 4z_ABAB + 3z_ABBB
! dp_AB/d_AAAB = 3 etc., calculated above when extracting yf
            dpst(1,ia,ib,jj)=yf(1,ia)*dyf(2,ib,jj)+dyf(1,ia,jj)*yf(2,ib)
            dpst(2,ia,ib,jj)=yf(1,ia)*dyf(3,ib,jj)+dyf(1,ia,jj)*yf(3,ib)
            dpst(3,ia,ib,jj)=yf(1,ia)*dyf(4,ib,jj)+dyf(1,ia,jj)*yf(4,ib)
            dpst(4,ia,ib,jj)=yf(2,ia)*dyf(3,ib,jj)+dyf(2,ia,jj)*yf(3,ib)
            dpst(5,ia,ib,jj)=yf(2,ia)*dyf(4,ib,jj)+dyf(2,ia,jj)*yf(4,ib)
            dpst(6,ia,ib,jj)=yf(3,ia)*dyf(4,ib,jj)+dyf(3,ia,jj)*yf(4,ib)
! using average values
!            dp12(ia,ib,jj)=0.5D0*(yf(1,ia)*dyf(2,ib,jj)+dyf(1,ia,jj)*yf(2,ib)+&
!                 dyf(1,ia,jj)*yf(2,ib)+yf(1,ia)*dyf(2,ib,jj))
!            dp13(ia,ib,jj)=0.5D0*(yf(1,ia)*dyf(3,ib,jj)+dyf(1,ia,jj)*yf(3,ib)+&
!                 dyf(1,ia,jj)*yf(3,ib)+yf(1,ia)*dyf(3,ib,jj))
!            dp14(ia,ib,jj)=0.5D0*(yf(1,ia)*dyf(4,ib,jj)+dyf(1,ia,jj)*yf(4,ib)+&
!                 dyf(1,ia,jj)*yf(4,ib)+yf(1,ia)*dyf(4,ib,jj))
!            dp23(ia,ib,jj)=0.5D0*(yf(2,ia)*dyf(3,ib,jj)+dyf(2,ia,jj)*yf(3,ib)+&
!                 dyf(2,ia,jj)*yf(3,ib)+yf(2,ia)*dyf(3,ib,jj))
!            dp24(ia,ib,jj)=0.5D0*(yf(2,ia)*dyf(4,ib,jj)+dyf(2,ia,jj)*yf(4,ib)+&
!                 dyf(2,ia,jj)*yf(4,ib)+yf(2,ia)*dyf(4,ib,jj))
!            dp34(ia,ib,jj)=0.5D0*(yf(3,ia)*dyf(4,ib,jj)+dyf(3,ia,jj)*yf(4,ib)+&
!                 dyf(3,ia,jj)*yf(4,ib)+yf(3,ia)*dyf(4,ib,jj))
! Hm 2nd derivatives ... should check ... kxsym requires ib>ia
!            zz=kxsym(ia,ib) 16*15/2=8*15=4*30=120
! let second derivatives be zeo
            cycle dploop
            do mm=jj,nc
               zz=kxsym(jj,mm)
               d2pst(1,ia,ib,zz)=dyf(1,ia,jj)*dyf(2,ib,mm)+&
                    dyf(1,ia,jj)*dyf(2,ib,mm)
            enddo
         enddo dploop
     enddo
  enddo ploop
  do ss=1,6
     write(*,320)'3X pst:',ss,((pst(ss,ia,ib),ia=1,ne),ib=1,ne)
320  format(a,i2,8F7.4)
  enddo
  do ia=1,ne
     do ib=1,ne
        write(*,*)'3X pair ',chel(ia),'-',chel(ib)
        do ss=1,6
           write(*,330)'3X dpst 1-8 : ',ss,(dpst(ss,ia,ib,jj),jj=1,8)
           write(*,330)'3X dpst 9-16: ',ss,(dpst(ss,ia,ib,jj),jj=9,16)
        enddo
     enddo
  enddo
330  format(a,i2,8F7.4)
!
!  write(*,*)'3X calculation of pair fractions done'
!  gx%bmperr=4399; goto 1000
!
! Recalculate constituent fractions from the site fractions  
! Note the ordering of the constituent important  !!!
   allocate(ylro(nc))
   ylro=zero
   jj=0
   do ia=1,ne
! this is AAAA or BBBB etc
      ybase=yf(1,ia)*yf(2,ia)*yf(3,ia)*yf(4,ia)
      jj=jj+1
      write(*,*)'3X ylro_',chel(ia)//chel(ia)//chel(ia)//chel(ia),&
           ' from yf(ss,ia): ',jj,ybase
      ylro(jj)=ybase
      do ib=ia+1,ne
         do ss=4,1,-1
            jj=jj+1
! to obtain 4 AAAB replace one yf(ss,ia) with yf(ss,ib) in sublattice ss
! NOTE it has to be in correct order: AAAB, AABA, ABAA and BAAA
            ylro(jj)=ybase*yf(ss,ib)/yf(ss,ia)
         enddo
! to obtain the 3 AABB in correct order: AABB, ABAB, ABBA
         jj=jj+1
         ylro(jj)=yf(1,ia)*yf(2,ia)*yf(3,ib)*yf(4,ib)
         jj=jj+1
         ylro(jj)=yf(1,ia)*yf(2,ib)*yf(3,ia)*yf(4,ib)
         jj=jj+1
         ylro(jj)=yf(1,ia)*yf(2,ib)*yf(3,ib)*yf(4,ia)
! to obtain the 3 more variants: BAAB, BABA, BBAA
         jj=jj+1
         ylro(jj)=yf(1,ib)*yf(2,ia)*yf(3,ia)*yf(4,ib)
         jj=jj+1
         ylro(jj)=yf(1,ib)*yf(2,ia)*yf(3,ib)*yf(4,ia)
         jj=jj+1
         ylro(jj)=yf(1,ib)*yf(2,ib)*yf(3,ia)*yf(4,ia)
! to obtain the 4 variants ABBB, similar to AAAB
         ybase=yf(1,ib)*yf(2,ib)*yf(3,ib)*yf(4,ib)
         do ss=1,4
            jj=jj+1
! to obtain 4 ABBB replace one yf(ss,ib) with yf(ss,ia) in sublattice ss
! NOTE it has to be in correct order: ABBB, BABB, BBAB and BBBA
            ylro(jj)=ybase*yf(ss,ia)/yf(ss,ib)
         enddo
! element C and D ...
! These lines not needed for binary systems ... NOT IMPLEMENTED
         do ic=ib+1,ne
            write(*,*)'3X ternary LRO not implemented'
            gx%bmperr=4399; goto 1000
            jj=jj+1
!           ylro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ia)**2*xf(ib)*xf(ic)
!            jj=jj+1
!           ylro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ia)*xf(ib)**2*xf(ic)
!            jj=jj+1
!           ylro(jj)=yrest*phvar%yfr(jj)+ycorr*mijkl(jj)*xf(ia)*xf(ib)*xf(ic)**2
            do id=ic+1,ne
               jj=jj+1
!               ylro(jj)=yrest*phvar%yfr(jj)+&
!                    ycorr*mijkl(jj)*xf(ia)*xf(ib)*xf(ic)*xf(id)
            enddo
         enddo
      enddo
   enddo
!
! handle only SRO, set all pst same
   write(*,'(a,3F10.7)')'3X >>>>>>>>>>>>> only SRO <<<<<<<<<',xf(1),xf(2)
   do ia=1,ne
      do ib=1,ne
         do ss=1,6
            pst(ss,ia,ib)=xf(ia)*xf(ib)
            do jj=1,nc
               dpst(ss,ia,ib,jj)=dxf(ia,jj)*xf(ib)+xf(ia)*dxf(ib,jj)
            enddo
         enddo
      enddo
   enddo
! SKIP CORRECTION, then converges with (almost) exactly ideal start values
! pstijtest is maximum pair fraction using mole fractions calculated
!     from cluster fractions provided by the minimizer
! if almost the same as previous decrease ycorr
   if(abs(pstijtest-pstijsave).lt.1.0D-4) then
      ycorr=max(0.5D0*ycorr,1.0D-8)
! when yrest=1 we use the fractions from the minimizer
   else
      ycorr=min(2.0d0*ycorr,0.5D0)
   endif
   ycorr=zero
   yrest=1.0D0-ycorr
   write(*,'(a,1x,l,1x,12F6.3)')'3X corr: ',sdebug,&
        pstijtest,pstijsave,yrest,ycorr
! without this stupid dummy statement the calculations does not converge SRO
   write(dummy,'(a,1x,l,1x,12F6.3)')'3X corr: ',sdebug,&
        pstijtest,pstijsave,yrest,ycorr
! save pstij for next iteration
   pstijsave=pstijtest
! We have LRO, fractions on sublattices can be different
!   do jj=1,nc
!      ylro(jj)=phvar%yfr(jj)
!   enddo
!
!   if(SDEBUG) then
   write(*,'(a,8F7.4)')'3X yfr: ',(phvar%yfr(jj),jj=1,8)
   write(*,'(a,8F7.4)')'3X ylro:',(ylro(jj),jj=1,8)
   write(*,'(a,8F7.4)')'3X yfr: ',(phvar%yfr(jj),jj=9,nc)
   write(*,'(a,8F7.4)')'3X ylro:',(ylro(jj),jj=9,16)
!   endif
!
!      write(*,*)'3X calculation of clusters from site fractions done'
!   gx%bmperr=4399; goto 1000
!
! All necessary fractions and derivatives calculated, now the entropy
! S = -2 \sum_ijkl y_ijkl\ln(y_ijkl) +
!     6  \sum_ij p_ij\ln(p_ij) - 5 \sum_j x_j\ln(x_j)
! As we calculate dG/dT the signs of Kikuchi are inverse
! entropy contribution from the constituent fractions ------------------
   ssumy=zero
!   sylog1=zero
   do jj=1,nc
!      yfra=phvar%yfr(jj)
! use cluster factions recalculated from site fractions
      yfra=ylro(jj)
      if(yfra.lt.bmpymin) yfra=bmpymin
      if(yfra.gt.one) yfra=one
! yfra is divided by mijkl as it represent mijkl fractions
!      ylog=log(yfra/mijkl(jj))
      ylog=log(yfra)
! debugging
!      if(jj.le.5) ylog1(jj)=ylog
!      sylog1=sylog1+yfra*ylog
      ssumy=ssumy+syfact*yfra*ylog
      if(moded.gt.0) then
! dgval(1,1:nc,1) are derivative of G/RT wrt fraction 1:nc
! d2gval(ixsym(jj*(jk+1)/2,1) are 2nd derivarive of G/RT wrt fraction jj and jk
! dgval and d2gval are zero before this loop
         phvar%dgval(1,jj,1)=syfact*(one+ylog)
! ? T and y derivative
!         phvar%dgval(2,jj,1)=syfact*(mijkl(jj)+ylog)/tval
         phvar%dgval(2,jj,1)=syfact*(one+ylog)/tval
! 2nd derivative, each term depend on a single y fraction
         phvar%d2gval(kxsym(jj,jj),1)=syfact/yfra
      endif
   enddo
   phvar%gval(1,1)=ssumy
   phvar%gval(2,1)=ssumy/tval
!   do jj=1,nc
! use cluster fractions recalculated from site fractions !!!
!      phvar%yfr(jj)=ylro(jj)
!   enddo
   write(*,*)'3X calculation of cluster entropy term done',ssumy
!  gx%bmperr=4399; goto 1000
! entropy contributions from the 6 pair fractions -----------------------
! each can be different pst(ss,ia,ib) can be different for each ss
   ssump=zero
   do ia=1,ne
      do ib=1,ne
         ssloop: do ss=1,6
! bond between atom ia and ib in sublattice s,t, there are 6 such paris
! we must repeat this for all 6 pairs in each cluster
! Note this term is negative in the Kikuchi summation !!!
            pfra=pst(ss,ia,ib)
            if(pfra.lt.bmpymin) pfra=bmpymin
            if(pfra.gt.one) pfra=one
            ylog=log(pfra)
            ssump=ssump-pfra*ylog
            if(moded.gt.0) then
               do jl=1,nc
! we need derivatives of pair fractions wrt constituent fractions
! Note negative sign because this is part of the 6 pln(p) term
! NOTE:  Make use of the fact that pstij depend direcly on cluster fractions
! Just using relations to y^s_A fauvours z_AAAA and z_BBBB, i.e. no mixing
! Use p_AB = 3z_AAAB + 4z_ABAB+3z_ABBB
! dp_AB/d_AAAB = 3 etc.
                  if(ia.eq.ib) then
                     phvar%dgval(1,jl,1)=phvar%dgval(1,jl,1)-&
                          dpst(ss,ia,ib,jl)*(one+ylog)
                  else
! this is derivative of AB pair related to clusters with ABBB, ABAB and ABBB
                     phvar%dgval(1,jl,1)=phvar%dgval(1,jl,1)-&
                          dqst(1,ia,ib,jl)*(one+ylog)
                  endif
                  do jm=jl,nc
                     zz=ixsym(jl,jm)
                     phvar%d2gval(zz,1)=phvar%d2gval(zz,1)-&
                          d2pst(ss,ia,ib,zz)*(one+ylog)+&
                          dpst(ss,ia,ib,jl)*dpst(ss,ia,ib,jm)/pst(ss,ia,ib)
                  enddo
               enddo
            endif
         enddo ssloop
! done all 6 pairs for ia,ib
      enddo
   enddo
! done all pairs ia,ib
! We have changed the sign when we calculate 6*p*ln(p) terms above
! The derivatives have the negative sign already
   phvar%gval(1,1)=phvar%gval(1,1)+ssump
   phvar%gval(2,1)=phvar%gval(2,1)+ssump/tval
   write(*,*)'3X calculation of 6 pair fractions entropy term done',ssump
! entropy contributions from the site fractions, factor 5/4 = 1.25
   ssumx=zero
!   sylog3=zero
   do ss=1,4
      do ia=1,ne
! we need derivatives of mole fractions wrt constituent fractions
         pfra=yf(ss,ia)
         if(pfra.lt.bmpymin) pfra=bmpymin
         if(pfra.gt.one) pfra=one
         ylog=log(pfra)
         ssumx=ssumx+sxfact*pfra*ylog
         if(moded.gt.0) then
            do jl=1,nc
! we need derivatives of site fractions wrt constituent fractions
               phvar%dgval(1,jl,1)=phvar%dgval(1,jl,1)+&
                    sxfact*dyf(ss,ia,jl)*(ylog+one)
               do jm=jl,nc
! note d2xf/dyidyj==0
                  zz=kxsym(jl,jm)
                  phvar%d2gval(zz,1)=phvar%d2gval(zz,1)+&
                       sxfact*dyf(ss,ia,jl)*dyf(ss,ia,jm)/pfra
               enddo
            enddo
         endif
      enddo
   enddo
! the sxfact +5/4 is already included above
   phvar%gval(1,1)=phvar%gval(1,1)+ssumx
   phvar%gval(2,1)=phvar%gval(2,1)+ssumx/tval
   write(*,'(a,F12.4)')'3X calculation of site fraction entropy: ',ssumx
!   write(*,*)'3X all entropy calculations made'
! All done
!   if(SDEBUG) then
      write(*,900)'3X xf: ',(xf(ia),ia=1,ne)
900   format(a,6F10.7)
      do ss=1,6
         write(*,320)'3X pst',ss,((pst(ss,ia,ib),ia=1,ne),ib=1,ne)
      enddo
910   format(a,10F7.3)
      write(*,310)((yf(ss,ia),ss=1,4),ia=1,ne)
      do ia=1,ne
         write(*,270)1,ia,(dyf(1,ia,jj),jj=1,16)
         write(*,270)2,ia,(dyf(2,ia,jj),jj=1,16)
         write(*,270)3,ia,(dyf(3,ia,jj),jj=1,16)
         write(*,270)4,ia,(dyf(4,ia,jj),jj=1,16)
      enddo
      write(*,'(a,8F7.4)')'3X yfr  1-8: ',(phvar%yfr(jj),jj=1,8)
      write(*,'(a,8F7.4)')'3X ylro 1-8: ',(ylro(jj),jj=1,8)
      write(*,'(a,8F7.4)')'3X yfr  9-16:',(phvar%yfr(jj),jj=9,nc)
      write(*,'(a,8F7.4)')'3X ylro 9-16:',(ylro(jj),jj=9,16)
      write(*,'(a,5(1pe12.4))')'3X cvmtfs: ',ssumy,ssump,ssumx,&
           phvar%gval(1,1),8.31451*phvar%gval(1,1)
!   write(*,*)'3X cvmtfs model not yet released'
!   endif
1000 continue
   return
 end subroutine config_entropy_cvmtfl

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine config_entropy_i2sl
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
   integer ll,kk,kall,nk,j1,j2,jxsym
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
!            phvar%d2gval(ixsym(kall,kall),1)=phvar%sites(ll)/yfra
            phvar%d2gval(kxsym(kall,kall),1)=phvar%sites(ll)/yfra
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
! to avoid calling ixsym ...
         jxsym=kxsym(j1,j1)
      cation2: do j2=j1,nkl(1)
! d2S/dy_i1dy_i2 = v_i1*y_Va*(1+ln(y_i2) + v_i2*y_Va*(1+ln(y_i1) + 
!                  [P*(1/y_i1**2)]         ..last term already calculated  OK
         if(ixsym(j1,j2).ne.jxsym) then
! this ixsym test works and has run of few 1000 times, removed for speed!!
            write(*,*)'3X ISYM error 5',j1,j2,ixsym(j1,j2),jxsym
            stop "3X ixsym indexing error 17"
         endif
!         phvar%d2gval(ixsym(j1,j2),1)=phvar%d2gval(ixsym(j1,j2),1)+&
!              (phvar%dpqdy(j1)*phvar%dgval(1,j2,1)+&
!               phvar%dpqdy(j2)*phvar%dgval(1,j1,1))*yva/phvar%sites(1)
         phvar%d2gval(jxsym,1)=phvar%d2gval(jxsym,1)+&
              (phvar%dpqdy(j1)*phvar%dgval(1,j2,1)+&
               phvar%dpqdy(j2)*phvar%dgval(1,j1,1))*yva/phvar%sites(1)
         jxsym=jxsym+j2
      enddo cation2
      anion2: do kk=1,nkl(2)
         j2=nkl(1)+kk
         jxsym=kxsym(j1,j2)
         if(j2.lt.min(i2slx(1),i2slx(2))) then
! d2S/dy_idy_j   = v_i*(1+ln(y_j)) + (-v_j)*(1+ln(y_i))    ...cation+anion OK
!            phvar%d2gval(ixsym(j1,j2),1)=&
!                 phvar%dpqdy(j1)*phvar%dgval(1,j2,1)/phvar%sites(2)+&
!                 phvar%dpqdy(j2)*phvar%dgval(1,j1,1)/phvar%sites(1)
            phvar%d2gval(jxsym,1)=&
                 phvar%dpqdy(j1)*phvar%dgval(1,j2,1)/phvar%sites(2)+&
                 phvar%dpqdy(j2)*phvar%dgval(1,j1,1)/phvar%sites(1)
         elseif(j2.eq.i2slx(1)) then
! d2S/dy_idy_Va  = v_i*(1+ln(y_Va)) + v_i*S1 + Q*(1+ln(y_i))   ...cation+Va OK
!            phvar%d2gval(ixsym(j1,j2),1)=&
!                 phvar%dpqdy(j1)*phvar%dgval(1,j2,1)/phvar%sites(2)+&
!                 phvar%dpqdy(j1)*spart(1)+&
!                 phvar%sites(2)*phvar%dgval(1,j1,1)/phvar%sites(1)
            phvar%d2gval(jxsym,1)=&
                 phvar%dpqdy(j1)*phvar%dgval(1,j2,1)/phvar%sites(2)+&
                 phvar%dpqdy(j1)*spart(1)+&
                 phvar%sites(2)*phvar%dgval(1,j1,1)/phvar%sites(1)
         else
! d2S/dy_idy_k   = v_i*(1+ln(y_k))                        ...cation+neutral OK
!            write(*,107)'3X i,va: ',j1,j2,phvar%dpqdy(j1),phvar%dgval(1,j2,1),&
!                 phvar%sites(2)
!107         format(a,2i2,6(1pe12.4))
!            phvar%d2gval(ixsym(j1,j2),1)=&
!                 phvar%dpqdy(j1)*phvar%dgval(1,j2,1)/phvar%sites(2)
            phvar%d2gval(jxsym,1)=&
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

!\addtotable subroutine config_entropy_qcwithlro
!\begin{verbatim}
 subroutine config_entropy_qcwithlro(moded,ncon,phvar,phrec,tval)
!
! calculates configurational entropy/R for the quasichemial liquid with LRO
!
! moded=0 only G, =1 G and dG/dy, =2 G, dG/dy and d2G/dy1/dy2
! ncon is number of constituents
! phvar is pointer to phase_varres record
! phrec is the phase record
! tval is current value of T
   implicit none
   integer moded,ncon
   TYPE(gtp_phase_varres), pointer :: phvar
   TYPE(gtp_phaserecord) :: phrec
   double precision tval
!\end{verbatim} %+
! First A=(z/2)*(\sum_i (y_ii*ln(y_ii) + \sum_(j>=i) y_ij*ln(y_ij/2))
! and calculate all x_i = y_ii + \sum_j a/(a+b)*y_ij
! Then calculate the SRO: q_ij=(y_ij/(x_i*x_j)-1)*(x_i+x_j)**2
! and B=\sum_i x_i*ln(x_i)*(1-z + \sum_(j>i) (z/2-1)*f(q_ij))
! -S = A+B
   integer icon,loksp,lokel,iel,nqij,kqij,jxsym,infirst,lat2,i,j
   double precision zhalf,yfra,ylog,cluster,sbonds,scorr,stoi1,stoi2
   double precision xp,xs,gamma,x1,x2,sumx(2),gamma2
   double precision, allocatable, dimension(:) :: qij,ycluster,&
        dgamma,d2gamma
   double precision, allocatable, dimension(:,:) :: xval
   double precision, allocatable, dimension(:,:,:) :: dxval
   integer, allocatable, dimension(:,:) :: qxij
   logical iscluster
   double precision, parameter :: half=0.5D0
!
   zhalf=half*phvar%qcbonds
   allocate(xval(noofel,2))
   allocate(dxval(noofel,ncon,2))
!   allocate(ycluster(noofel))
   xval=zero
   dxval=zero
!   write(*,*)'3X classical qc with LRO!',zhalf
!   gx%bmperr=4399; goto 1000
!
   sbonds=zero
   nqij=0
   sumx=zero
   ally: do icon=1,ncon
      yfra=phvar%yfr(icon)
      if(yfra.lt.bmpymin) yfra=bmpymin
      if(yfra.gt.one) yfra=one
! loksp is set to the index of the constituent in the species array
      loksp=phrec%constitlist(icon)
! if two elements it is an AB bond
! To identify if the cluster constituent is on the first or second sublattice
! use the alphabetical order of the species name.  If first letter<second
! then the first element is in the first sublattice:
! thus AB means first  element in first sublattice,
!      BA means second element in first sublattice
! The elements are always ordered alphabetically in splista(loksp)%ellinks
      infirst=1
      if(splista(loksp)%noofel.eq.2) then
         cluster=half
         iscluster=.TRUE.
!         write(*,*)'3X CQC classic 0: ',qcmodel,iscluster,yfra
         if(splista(loksp)%symbol(1:1).gt.splista(loksp)%symbol(2:2)) then
! this is constituent BA
            infirst=2
         endif
      elseif(splista(loksp)%noofel.eq.1) then
! same element in both sublattices
         cluster=one
         iscluster=.FALSE.
      else
         write(*,*)'3X cluster with too many elements'
         gx%bmperr=4399; goto 1000
      endif
      ylog=log(yfra)
! gval(1:6,1) are G and derivator wrt T and P
! dgval(1,1:N,1) are derivatives of G wrt fraction 1:N
! dgval(2,1:N,1) are derivatives of G wrt fraction 1:N and T
! dgval(3,1:N,1) are derivatives of G wrt fraction 1:N and P
! d2dval(ixsym(N*(N+1)/2),1) are derivatives of G wrt fractions N and M
! this is a symmetric matrix and index givem by ixsym(M,N)
      sbonds=sbonds+zhalf*yfra*ylog
      if(moded.gt.0) then
         phvar%dgval(1,icon,1)=zhalf*(one+ylog)
         phvar%d2gval(kxsym(icon,icon),1)=zhalf/(yfra)
      endif
! lokel is index in ellista of first element in alphabetical order of element
      lokel=splista(loksp)%ellinks(1)
!      write(*,17)'3X qccorr: ',trim(splista(loksp)%symbol),icon,loksp,lokel,&
!           infirst,iscluster,yfra
17    format(a,a,4i4,l3,F7.4)
      if(iscluster) then
         nqij=nqij+1
! if a bond cluster there must be two elements         
         iel=ellista(lokel)%alphaindex
         stoi1=splista(loksp)%stoichiometry(1)
         stoi2=splista(loksp)%stoichiometry(2)
         xval(iel,infirst)=xval(iel,infirst)+stoi1/(stoi1+stoi2)*yfra
         sumx(infirst)=sumx(infirst)+stoi1/(stoi1+stoi2)*yfra
         dxval(iel,icon,infirst)=stoi1/(stoi1+stoi2)
!         write(*,60)'3X qc 3A: ',infirst,iel,yfra,((xval(i,j),i=1,2),j=1,2)
60       format(a,2i4,F7.3,': ',10F7.3)
         lokel=splista(loksp)%ellinks(2)
         iel=ellista(lokel)%alphaindex
         xval(iel,3-infirst)=xval(iel,3-infirst)+stoi2/(stoi1+stoi2)*yfra
         sumx(3-infirst)=sumx(3-infirst)+stoi2/(stoi1+stoi2)*yfra
         dxval(iel,icon,3-infirst)=stoi2/(stoi1+stoi2)
!         write(*,60)'3X qc 3B: ',3-infirst,iel,yfra,((xval(i,j),i=1,2),j=1,2)
      else
! the same element in both sublattices, we already know lokel
!         lokel=splista(loksp)%ellinks(1)
         iel=ellista(lokel)%alphaindex
         xval(iel,1)=xval(iel,1)+half*yfra
         sumx(1)=sumx(1)+half*yfra
         dxval(iel,icon,1)=half
         xval(iel,2)=xval(iel,2)+half*yfra 
         sumx(2)=sumx(2)+half*yfra
         dxval(iel,icon,2)=half
!         write(*,60)'3X qc 3C: ',1,iel,yfra,((xval(i,j),i=1,2),j=1,2)
      endif
!      write(*,60)'3X sumx: ',icon,0,yfra,sumx
!      write(*,'(a,2i2,": ",8(i2,F6.2))')'3X dx:',icon,1,&
!           ((iel,dxval(iel,i,1),i=1,ncon),iel=1,noofel)
!      write(*,'(a,2i2,": ",8(i2,F6.2))')'3X dx:',icon,2,&
!           ((iel,dxval(iel,i,2),i=1,ncon),iel=1,noofel)
   enddo ally
!----------------------------------------
! Here we have all x values and derivatives
! The correction term is composition independent 1-z
!   gamma=one-2.0D0*zhalf
! factor 0.5 gives OK SRO but no LRO
   gamma=0.5D0*(one-2.0D0*zhalf)     ! OK SRO but no LRO
!   gamma=sumx(1)*(one-2.0D0*zhalf) no improvement
!   gamma=0.75D0*(one-2.0D0*zhalf) totally wrong
!   gamma=0.25D0*(one-2.0D0*zhalf) Very bad
! THIS FACTOR 2 MAKES LRO STABLE ... BUT IS IT CORRECT???
   gamma2=2.0D0*gamma
! MAYBE THERE IS SOME ERROR IN THE DERIVATIVES BELOW?
! Some elements may not be dissolved in this phase ??
!   write(*,'(a,7F7.3)')'3X x1:',sumx(1),(xval(iel,1),iel=1,noofel)
!   write(*,'(a,5F7.3)')'3X x2:',sumx(2),(xval(iel,2),iel=1,noofel)
   do lat2=1,2
      do iel=1,noofel
         xval(iel,lat2)=xval(iel,lat2)/sumx(lat2)
      enddo
   enddo
!   write(*,'(a,5F7.3)')'3X QCLRO: ',sumx,gamma,gamma2
!   write(*,'(a,7F7.3)')'3X x3:',sumx(1),(xval(iel,1),iel=1,noofel)
!   write(*,'(a,7F7.3)')'3X x4:',sumx(2),(xval(iel,2),iel=1,noofel)
!   write(*,'(a,8(i2,F7.3))')'3X x5:',((iel,dxval(iel,icon,1),&
!        icon=1,ncon),iel=1,noofel)
!   write(*,'(a,8(i2,F7.3))')'3X x6:',((iel,dxval(iel,icon,2),&
!        icon=1,ncon),iel=1,noofel)
   scorr=zero
   sub2: do lat2=1,2
      allx: do iel=1,noofel
         yfra=xval(iel,lat2)
         if(yfra.le.bmpymin) yfra=bmpymin
         if(yfra.gt.one) yfra=one
         ylog=log(yfra)
! this is the contribution to integral G, multiplied with gamma after the loop
         scorr=scorr+yfra*ylog
! WE MUST ALSO CALCULATE DERIVATIVES OF x_i wrt y USING CHAIN RULE
         if(moded.gt.0) then
            ally2: do icon=1,ncon
! dgval(1,1:N,1) are derivatives of G wrt fraction 1..N
! dgval(2,1:N,1) are derivatives of G wrt fraction 1..N and T
! dgval(3,1:N,1) are derivatives of G wrt fraction 1..N and P
! d2dval(ixsym(N*(M+1)/2),1) are derivatives of G wrt fractions N and M
               phvar%dgval(1,icon,1)=phvar%dgval(1,icon,1)+&
                    gamma2*(one+ylog)*dxval(iel,icon,lat2)
               jxsym=kxsym(icon,icon)
               do loksp=icon,ncon
                  if(ixsym(icon,loksp).ne.jxsym) then
! this ixsym test works and has run of few 1000 times, removed for speed!!
                     write(*,*)'3X KSYM error 18',ixsym(icon,loksp),jxsym
                     stop
                  endif
                  phvar%d2gval(jxsym,1)=phvar%d2gval(jxsym,1)+&
                       gamma2*dxval(iel,icon,lat2)*dxval(iel,loksp,lat2)/yfra
! this replaces call to ixsym(loksp,icon)
                  jxsym=jxsym+loksp
               enddo
            enddo ally2
         endif
      enddo allx
   enddo sub2
!
!   write(*,'(a,8(i2,F7.3))')'3X x5:',((iel,dxval(iel,ncon,1),iel=1,noofel),&
!
!- ixsym --------------- ixsym end modification
! now all is calculated gval(1,1)=G; gval(2,1)=S etc
   write(*,'(a,4(1pe12.4))')'3X scorr: ',tval,gamma2,scorr,sbonds
   phvar%gval(1,1)=sbonds+gamma*scorr
   phvar%gval(2,1)=(sbonds+gamma*scorr)/tval
!   write(*,12)'3X QCLRO: ',phvar%gval(1,1),phvar%gval(2,1),gamma,&
!        zhalf,sbonds,scorr
12 format(a,6(1pe11.3))
!
1000 continue
   return
! NO LRO .... SUCK .... but LRO by doubling gamma2, why??
 end subroutine config_entropy_qcwithlro

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine config_entropy_cqc_classicqc
!\begin{verbatim}
 subroutine config_entropy_cqc_classicqc(moded,ncon,phvar,phrec,tval)
!
! calculates configurational entropy/R for the quasichemial liquid with LRO
!
! THIS ROUTINE NOT USED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! THIS WORKS OK 2019-01-10: DO NOT CHANGE ANYTHING!! works for qcmodel=1
! routine for qcmodel=2 and 3 laret
!
! only question is the parameter, value of G = K*T*R/2;
! K=-10, T=600 means G= -10*600*R/2 = -3000*R gives same curves as in paper.   
!
! moded=0 only G, =1 G and dG/dy, =2 G, dG/dy and d2G/dy1/dy2
! ncon is number of constituents
! phvar is pointer to phase_varres record
! phrec is the phase record
! tval is current value of T
   implicit none
   integer moded,ncon
   TYPE(gtp_phase_varres), pointer :: phvar
   TYPE(gtp_phaserecord) :: phrec
   double precision tval
!\end{verbatim} %+
! First A=(z/2)*(\sum_i (y_ii*ln(y_ii) + \sum_(j>=i) y_ij*ln(y_ij/2))
! and calculate all x_i = y_ii + \sum_j a/(a+b)*y_ij
! Then calculate the SRO: q_ij=(y_ij/(x_i*x_j)-1)*(x_i+x_j)**2
! and B=\sum_i x_i*ln(x_i)*(1-z + \sum_(j>i) (z/2-1)*f(q_ij))
! -S = A+B
   integer icon,loksp,lokel,iel,nqij,kqij,jxsym
   double precision zhalf,yfra,ylog,cluster,sbonds,scorr,stoi1,stoi2
   double precision xp,xs,gamma,x1,x2
   double precision, allocatable, dimension(:) :: xval,qij,ycluster,&
        dgamma,d2gamma
   double precision, allocatable, dimension(:,:) :: dxval,dqij
   integer, allocatable, dimension(:,:) :: qxij
   logical iscluster
   double precision, parameter :: half=0.5D0
!
! qcmodel=1 is classical qc without LRO, 2 is q**2, 3 is 0.5*(1+q)*q**2
!   qcmodel=1
!
   zhalf=half*phvar%qcbonds
   write(*,*)'3X classic cqc, zhalf: ',zhalf
   allocate(xval(noofel))
   allocate(dxval(noofel,ncon))
!   allocate(ycluster(noofel))
   xval=zero
   dxval=zero
!   write(*,*)'3X classical quasichemical!',zhalf
!
   sbonds=zero
   nqij=0
   do icon=1,ncon
      yfra=phvar%yfr(icon)
      if(yfra.lt.bmpymin) yfra=bmpymin
      if(yfra.gt.one) yfra=one
! if set the constituent is a binary cluster
      if(btest(phvar%constat(icon),CONQCBOND)) then
         cluster=half
         iscluster=.TRUE.
!         write(*,*)'3X CQC classic 0: ',qcmodel,iscluster,yfra
      else
         cluster=one
         iscluster=.FALSE.
      endif
! entropy is y*ln(y) for single atoms, y*ln(y/2) for clusters
      ylog=log(cluster*yfra)
! gval(1:6,1) are G and derivator wrt T and P
! dgval(1,1:N,1) are derivatives of G wrt fraction 1:N
! dgval(2,1:N,1) are derivatives of G wrt fraction 1:N and T
! dgval(3,1:N,1) are derivatives of G wrt fraction 1:N and P
! d2dval(ixsym(N*(N+1)/2),1) are derivatives of G wrt fractions N and M
! this is a symmetric matrix and index givem by ixsym(M,N)
      sbonds=sbonds+zhalf*yfra*ylog
      if(moded.gt.0) then
         phvar%dgval(1,icon,1)=zhalf*(one+ylog)
!         phvar%d2gval(ixsym(icon,icon),1)=zhalf/(yfra)
         phvar%d2gval(kxsym(icon,icon),1)=zhalf/(yfra)
! These should be the correct derivatives but with these it does not converge!!
!         phvar%dgval(1,icon,1)=zhalf*(one/cluster+ylog)
!         phvar%d2gval(ixsym(icon,icon),1)=zhalf/(cluster*yfra)
! DO NOT CHANGE ANYTHING!!
      endif
! loksp is set to the index of the species array
      loksp=phrec%constitlist(icon)
      lokel=splista(loksp)%ellinks(1)
!      if(btest(phvar%constat(icon),CONQCBOND)) then
      if(iscluster) then
         nqij=nqij+1
!         ycluster(nqij)=yfra
! if a bond cluster there must be two elements         
! lokel is index in ellista of first element, iel is its alphabetical index
         iel=ellista(lokel)%alphaindex
         stoi1=splista(loksp)%stoichiometry(1)
         stoi2=splista(loksp)%stoichiometry(2)
         xval(iel)=xval(iel)+stoi1/(stoi1+stoi2)*yfra
         dxval(iel,icon)=stoi1/(stoi1+stoi2)
!         write(*,60)'3X qc 3A: ',iel,xval(iel),yfra
         lokel=splista(loksp)%ellinks(2)
         iel=ellista(lokel)%alphaindex
         xval(iel)=xval(iel)+stoi2/(stoi1+stoi2)*yfra
         dxval(iel,icon)=stoi2/(stoi1+stoi2)
!         write(*,60)'3X qc 3B: ',iel,xval(iel),yfra
      else
         lokel=splista(loksp)%ellinks(1)
         iel=ellista(lokel)%alphaindex
         xval(iel)=xval(iel)+yfra
         dxval(iel,icon)=one
!         write(*,60)'3X qc 3C: ',iel,xval(iel),yfra
      endif
   enddo
!----------------------------------------
! We do not need qij = y_ij/(2x_ix_j) - 1
! The correction term is composition independent 1-z
   gamma=one-2.0D0*zhalf
! Some elements may not be dissolved in this phase ??
   scorr=zero
   do iel=1,noofel
      yfra=xval(iel)
      if(yfra.le.bmpymin) yfra=bmpymin
      if(yfra.gt.one) yfra=one
      ylog=log(yfra)
      scorr=scorr+yfra*ylog
! WE MUST ALSO CALCULATE DERIVATIVES OF x_i USING CHAIN RULE
! DO NOT CHANGE ANYTHING!!
      if(moded.gt.0) then
         do icon=1,ncon
            phvar%dgval(1,icon,1)=phvar%dgval(1,icon,1)+&
                 gamma*(one+ylog)*dxval(iel,icon)
            jxsym=kxsym(icon,icon)
            do loksp=icon,ncon
!               phvar%d2gval(ixsym(icon,loksp),1)=&
!                    phvar%d2gval(ixsym(icon,loksp),1)+&
!                    gamma*dxval(iel,icon)*dxval(iel,loksp)/yfra
               if(ixsym(icon,loksp).ne.jxsym) then
! this ixsym test works and has run of few 1000 times, removed for speed!!
                  write(*,*)'3X ISYM error 18',ixsym(icon,loksp),jxsym
                  stop
               endif
               phvar%d2gval(jxsym,1)=&
                    phvar%d2gval(jxsym,1)+&
                    gamma*dxval(iel,icon)*dxval(iel,loksp)/yfra
               jxsym=jxsym+loksp
            enddo
         enddo
      endif
   enddo
!- ixsym --------------- ixsym end modification
! now all is calculated gval(1,1)=G; gval(2,1)=S etc
! DO NOT CHANGE ANYTHING!!
   phvar%gval(1,1)=sbonds+gamma*scorr
   phvar%gval(2,1)=(sbonds+gamma*scorr)/tval
   write(*,12)'3X QC1: ',qcmodel,phvar%gval(1,1),phvar%gval(2,1),gamma,&
        zhalf,sbonds,scorr
12 format(a,i2,6(1pe11.3))
!
! THIS ROUTINE NOT USED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
1000 continue
   return
 end subroutine config_entropy_cqc_classicqc

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine config_entropy_qchillert
!\begin{verbatim}
 subroutine config_entropy_qchillert(moded,ncon,phvar,phrec,tval)
!
! calculates configurational entropy/R for the corrected quasichemial liquid
! Hillert-Selleby-Sundman
! Rewritten 2019-01-12 based on cqc-classicqc which seems correct
! test1: qcmodel=1: OK for zhalf=1 and 3
! test2: qcmodel=2: OK!!
! test3: qcmodel=3: OK for SRO, problmens for miscibility gap
!
! only question is the parameter, value of G = K*T*R/2;
! K=-10, T=600 means G= -10*600*R/2 = -3000*R gives same curves as in paper.   
!
! moded=0 only G, =1 G and dG/dy, =2 G, dG/dy and d2G/dy1/dy2
! ncon is number of constituents
! phvar is pointer to phase_varres record
! phrec is the phase record
! tval is current value of T
   implicit none
   integer moded,ncon
   TYPE(gtp_phase_varres), pointer :: phvar
   TYPE(gtp_phaserecord) :: phrec
   double precision tval
!\end{verbatim}
! First A=(z/2)*(\sum_i (y_ii*ln(y_ii) + \sum_(j>=i) y_ij*ln(y_ij/2))
! and calculate all x_i = y_ii + \sum_j (a_i/(a_1+a_j))*y_ij
! and calculate all x_j = y_jj + \sum_i (a_j/(a_1+a_j))*y_ij
!         dx_i/dy_ii =1; dx_i/dy_ij = a_i/(a_i+a_j); dx_i/dy_jj =0
! Then calculate normallized sro= q_ij=(y_ij/(x_i*x_j)-1)*(x_i+x_j)**2
! gcmodel=1 : gamma = -(1-z)
! gcmodel=2 : gamma = -(1-z -(z/2-1)*sro**2)
!         gamma is multiplied with B:
! and B=\sum_i x_i*ln(x_i)
! -S = dG/dT = A+gamma*B
   integer icon,loksp,lokel,iel,nqij,kqij,jcon,kcon,maxc
   double precision zhalf,yfra,ylog,cluster,sbonds,scorr,stoi1,stoi2,temp1
   double precision xp,xs,gamma,x1,x2,sij,xnorm1,xnorm2,xprod
   double precision, allocatable, dimension(:) :: xval,qij,ycluster,&
        dgamma,d2gamma,sqz
   double precision, allocatable, dimension(:,:) :: dxval,dqij,d2qij
   integer, allocatable, dimension(:,:) :: qxij,jel
   integer, allocatable, dimension(:) :: jcluster
   logical iscluster
   double precision, parameter :: half=0.5D0
!
! qcmodel=1 is classical qc without LRO, 2 is q**2, 3 is 0.5*(1+q)*q**2
!   qcmodel=1
!   gcmodel=2
!   gcmodel=3
!
   zhalf=half*phvar%qcbonds
! ncon is number of constituents, noofel number of elements
   write(*,'(a,i2,F8.3,10i4)')'3X cqc6 start: ',qcmodel,zhalf,ncon,noofel
   allocate(xval(noofel))
   allocate(dxval(noofel,ncon))
! max antal binary cluster
   maxc=ncon*(ncon-1)/2
   allocate(ycluster(maxc))
   allocate(sqz(maxc))
   allocate(jel(maxc,2))
! this is used to indicate constituent index of a cluster
   allocate(jcluster(maxc))
   xval=zero
   dxval=zero
   sqz=zero
!   write(*,*)'3X classical quasichemical!',zhalf
!
! STEP 1: entropy for clusters
   sbonds=zero
   nqij=0
   scluster: do icon=1,ncon
      yfra=phvar%yfr(icon)
      if(yfra.lt.bmpymin) yfra=bmpymin
      if(yfra.gt.one) yfra=one
! if set the constituent is a binary cluster
      if(btest(phvar%constat(icon),CONQCBOND)) then
         cluster=half
         iscluster=.TRUE.
!         write(*,*)'3X CQC classic 0: ',qcmodel,iscluster,yfra
      else
         cluster=one
         iscluster=.FALSE.
      endif
! entropy is y*ln(y) for single atoms, y*ln(y/2) for clusters
      ylog=log(cluster*yfra)
! gval(1:6,1) are G, dG/dT, dG/dP, d2G/dT2, d2G/dTdP and d2G/dP2
! dgval(1,1:N,1) are derivatives of G wrt fraction 1:N
! dgval(2,1:N,1) are derivatives of G wrt fraction 1:N and T
! dgval(3,1:N,1) are derivatives of G wrt fraction 1:N and P
! d2dval(ixsym(N*(N+1)/2),1) are derivatives of G wrt fractions N and M
! this is a symmetric matrix and index givem by ixsym(M,N)
      sbonds=sbonds+zhalf*yfra*ylog
      if(moded.gt.0) then
! first and second derivatives for \sum_i y_ii\ln(y_ii)+y_ij\ln(y_ij/2)
         phvar%dgval(1,icon,1)=zhalf*(one+ylog)
         phvar%d2gval(ixsym(icon,icon),1)=zhalf/(yfra)
      endif
! we have to calculate the mole fractions for the correction term
! loksp is set to the index of the species array
      loksp=phrec%constitlist(icon)
      lokel=splista(loksp)%ellinks(1)
!      if(btest(phvar%constat(icon),CONQCBOND)) then
      if(iscluster) then
         nqij=nqij+1
         ycluster(nqij)=yfra
         jcluster(nqij)=icon
! if a bond cluster there must be two elements         
! lokel is index in ellista of first element, iel is its alphabetical index
         iel=ellista(lokel)%alphaindex
         stoi1=splista(loksp)%stoichiometry(1)
         stoi2=splista(loksp)%stoichiometry(2)
         xval(iel)=xval(iel)+stoi1/(stoi1+stoi2)*yfra
         dxval(iel,icon)=stoi1/(stoi1+stoi2)
! sqz=1 for equiatomic ordering
         sqz(nqij)=0.125/(stoi1*stoi2)
         jel(nqij,1)=iel
!         write(*,60)'3X qc 3A: ',iel,xval(iel),yfra
         lokel=splista(loksp)%ellinks(2)
         iel=ellista(lokel)%alphaindex
         xval(iel)=xval(iel)+stoi2/(stoi1+stoi2)*yfra
         dxval(iel,icon)=stoi2/(stoi1+stoi2)
         jel(nqij,2)=iel
!         write(*,60)'3X qc 3B: ',iel,xval(iel),yfra
      else
         lokel=splista(loksp)%ellinks(1)
         iel=ellista(lokel)%alphaindex
         xval(iel)=xval(iel)+yfra
         dxval(iel,icon)=one
!         write(*,60)'3X qc 3C: ',iel,xval(iel),yfra
      endif
   enddo scluster
! check mole fractions sum up to unity
   xs=zero
   do icon=1,noofel
      xs=xs+xval(icon)
   enddo
   if(abs(xs-one).gt.1.0D-12) then
      write(*,*)'3X cqc6: sum of molefractions not unity: ',xs
      stop
   endif
!----------------------------------------
! step 2: correction factor gamma
! NOTE sign opposite eq. 11 as this is dG/dy = -S 
   if(qcmodel.eq.1) then
! classic: composition independent gamma: 1-z
      gamma=one-2.0D0*zhalf
   else
!      write(*,*)'3X qcmodel: ',qcmodel,nqij,maxc
! we must calculate the SRO for each cluster ij (more than one cluster!)
! s_ij = 0.5 y_ij - x_ix_j/(x_i+x_j)**2
! q_ij= s_ij/x_ix_j   
! and we can have ordering at other composition than equiatoic!!
      allocate(dgamma(ncon))
      allocate(d2gamma(ncon*(ncon+1)/2))
      allocate(qij(nqij))
      allocate(dqij(nqij,ncon))
      allocate(d2qij(nqij,ncon*(ncon+1)/2))
      gamma=one-2.0D0*zhalf
!      write(*,*)'3X loop to calculate gamma'
      gammaloop: do icon=1,nqij
! jel(icon,k) is element k in cluster icon
         xnorm1=xval(jel(icon,1))/(xval(jel(icon,1))+xval(jel(icon,2)))
         xnorm2=xval(jel(icon,2))/(xval(jel(icon,1))+xval(jel(icon,2)))
         xprod=xnorm1*xnorm2
! NOTE p_AB is 0.5y_AB
!         sij=0.5D0*ycluster(icon)/xprod
         sij=sqz(icon)*ycluster(icon)/xprod
! This is the variable "q" defined by eq. 9 in the 2009 paper 
!         qij(icon)=0.5D0*ycluster(icon)/xprod-one
         qij(icon)=sij-one
         do jcon=1,ncon
! first derivatives of qij with respect to y_ij
            temp1=zero
            if(jcluster(icon).eq.jcon) temp1=sqz(icon)/xprod
            temp1=temp1-sij*(dxval(jel(icon,1),jcon)/xval(jel(icon,1))-&
                 dxval(jel(icon,2),jcon)/xval(jel(icon,2)))
            dqij(icon,jcon)=temp1
! ignore 2nd derivatives ... ???
            do kcon=jcon,ncon
               d2qij(icon,ixsym(icon,jcon))=zero
            enddo
         enddo
         if(qcmodel.eq.2) then
! this is the correction factor
            gamma=gamma+(zhalf-one)*qij(icon)**2
            do jcon=1,ncon
! THIS LINE WAS MISSING!!! works (almost)
               dgamma(jcon)=2.0D0*(zhalf-one)*qij(icon)*dqij(icon,jcon)
               do kcon=jcon,ncon
! this is approximate, no d2qij....
                  d2gamma(ixsym(jcon,kcon))=2.0D0*(zhalf-one)*&
                       dqij(icon,kcon)*dqij(icon,jcon)
               enddo
            enddo
         elseif(qcmodel.eq.3) then
! this is the correction factor for qcmodel=3
            gamma=gamma+0.5D0*(zhalf-one)*(qij(icon)+one)*qij(icon)**2
!            write(*,*)'3X qcmodel=3: ',icon,phvar%phtupx,qij(icon),gamma
            do jcon=1,ncon
               dgamma(jcon)=(zhalf-one)*qij(icon)*&
                    (1.5D0*qij(icon)+one)*dqij(icon,jcon)
               do kcon=jcon,ncon
! this is approximate, no d2qij....
                  d2gamma(ixsym(jcon,kcon))=0.5D0*(zhalf-one)*&
                       (6.0D0*qij(icon)*dqij(icon,kcon)*dqij(icon,jcon)+&
                       (3*qij(icon)+one)*d2qij(icon,ixsym(jcon,kcon)))
               enddo
            enddo
         else
            write(*,*)'3X no such qcmodel: ',qcmodel
         endif
      enddo gammaloop
   endif
!   write(*,*)'3X qcmodel',qcmodel,gamma
!----------------------------------------
! Step 3: entropy for molefractions: scorr=\sum_i x_i ln(x_i)
   scorr=zero
!   write(*,'(a,4(1pe12.4))')'3X loop scorr',sbonds,gamma,sqz
   smol: do iel=1,noofel
      yfra=xval(iel)
      if(yfra.le.bmpymin) yfra=bmpymin
      if(yfra.gt.one) yfra=one
      ylog=log(yfra)
      scorr=scorr+yfra*ylog
! WE MUST ALSO CALCULATE DERIVATIVES OF x_i, dx_i/dy_j USING CHAIN RULE
! at present ignore derivatives of gamma ....
      if(moded.gt.0) then
         do icon=1,ncon
            phvar%dgval(1,icon,1)=phvar%dgval(1,icon,1)+&
                 gamma*(one+ylog)*dxval(iel,icon)
! derivative wrt T and icon
            phvar%dgval(2,icon,1)=phvar%dgval(2,icon,1)+&
                 gamma*(one+ylog)*dxval(iel,icon)/tval
            do jcon=icon,ncon
               phvar%d2gval(ixsym(icon,jcon),1)=&
                    phvar%d2gval(ixsym(icon,jcon),1)+&
                    gamma*dxval(iel,icon)*dxval(iel,jcon)/yfra
            enddo
         enddo
      endif
   enddo smol
!   write(*,*)'3X all done, save values in phvar'
! subtract the correction which depend on qcmodel
! now all is calculated gval(1,1)=G; gval(2,1)=S etc
! Second derivates for \usm_i y_i\ln(y_i) calculated above
   phvar%gval(1,1)=sbonds+gamma*scorr
   phvar%gval(2,1)=(sbonds+gamma*scorr)/tval
   if(qcmodel.gt.1) then
! include derivatives of gamma
      do icon=1,ncon
         phvar%dgval(1,icon,1)=phvar%dgval(1,icon,1)+dgamma(icon)*scorr
         phvar%dgval(2,icon,1)=phvar%dgval(2,icon,1)+dgamma(icon)*scorr/tval
         do jcon=icon,ncon
! approximate ....
            phvar%d2gval(ixsym(icon,jcon),1)=phvar%d2gval(ixsym(icon,jcon),1)+&
                 d2gamma(ixsym(icon,jcon))*scorr
         enddo
     enddo
   endif
!               
!   write(*,12)'3X cqc6: ',qcmodel,sbonds,gamma,scorr,&
!        phvar%gval(1,1),phvar%gval(2,1)
12 format(a,i2,6(1pe11.3))
!
1000 continue
   return
 end subroutine config_entropy_qchillert !gamma, dgamma, d2gamma

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine config_entropy_cvmce
!\begin{verbatim}
 subroutine config_entropy_cvmce(moded,ncon,phvar,phrec,tval)
!
! calculates the classical QC and CVM models sith LRO
! started 2021-02-17
!
! moded=0 only G, =1 G and dG/dy, =2 G, dG/dy and d2G/dy1/dy2
! ncon is number of constituents
! phvar is pointer to phase_varres record
! phrec is the phase record
! tval is current value of T
   implicit none
   integer moded,ncon
   TYPE(gtp_phase_varres), pointer :: phvar
   TYPE(gtp_phaserecord) :: phrec
   double precision tval
!\end{verbatim}
!---------------------------------------------------------------------------1
   write(*,*)'3X classical QC model with LRO, not implemented yet'
! S = - \sum_i y_i ln(y_i) + z/2 \sum_k x_k ln(x_k)
   gx%bmperr=4399
1000 continue
   return
 end subroutine config_entropy_cvmce

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine config_entropy_tisr
!\begin{verbatim}
 subroutine config_entropy_tisr(moded,ncon,phvar,phrec,tval)
!
! calculates configurational entropy/R for the Kremer liquid SRO model
! started 2021-02-12
!
! moded=0 only G, =1 G and dG/dy, =2 G, dG/dy and d2G/dy1/dy2
! ncon is number of constituents, each cell a constituent
! phvar is pointer to phase_varres record
! phrec is the phase record
! tval is current value of T
   implicit none
   integer moded,ncon
! to obtain current fractions and store results
   TYPE(gtp_phase_varres), pointer :: phvar
! to obtain phase and constituent inforamation
   TYPE(gtp_phaserecord) :: phrec
   double precision tval
!\end{verbatim}
   integer ii,loksp
   double precision rtg
!---------------------------------------------------------------------------
! 
! rtg is set to R*T
   rtg=globaldata%rgas*tval
   write(*,10)ncon,(trim(splista(phrec%constitlist(ii))%symbol),&
        phvar%yfr(ii),ii=1,ncon)
10 format('3X config_entropy_tisr with ',i2,' constituents, fractions:'/&
        10(a,1x,F7.4,', '))
! You can enter all calculations here, nothing will be added elsewhere
   write(*,'(a,1pe16.6)')'3X TISR not implemented yet, G=',phvar%gval(1,1)
   gx%bmperr=4399
! Values returned should be:
! phvar%gval(1,1)     Gibbs energy divided by RT (G/RT below)
! phvar%gval(2,1)     derivative of G/RT wrt T
! phvar%dgval(1,ii,1) 1st derivative of G/RT wrt fraction ii
! phvar%dgval(2,ii,1) 2nd derivative of G/RT wrt T amd fraction ii
! phvar%d2gval(ixym(ii,jj),1) 2nd derivative of G/RT wrt fracton ii and jj
!    Normally sufficient to set phvar%d2gval(ixsym(ii,ii))=one/phvar%yfr(ii)
!-----------------------------------
1000 continue
   return
 end subroutine config_entropy_tisr

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine config_entropy_srot
!\begin{verbatim}
 subroutine config_entropy_srot(moded,ncon,phvar,phrec,tval)
!
! calculates configurational entropy/R for tetrahedron SRO model
! I DO NOT THINK THIS IS USED ??
!
! moded=0 only G, =1 G and dG/dy, =2 G, dG/dy and d2G/dy1/dy2
! ncon is number of constituents, each cell a constituent
! phvar is pointer to phase_varres record
! phrec is the phase record
! tval is current value of T
   implicit none
   integer moded,ncon
   TYPE(gtp_phase_varres), pointer :: phvar
   TYPE(gtp_phaserecord) :: phrec
   double precision tval
!\end{verbatim}
!---------------------------------------------------------------------------
! FOR A BINARY, with 5 SRO clusters as in tetrahedron FCC 
! 
   integer ia,ib,ja,jb,kk,mm,alpha,jxsym,loksp
   double precision rk,ggg,rtg
!   double precision s1,s2,s11,s12,s21,s22
!   double precision s111,s112,s121,s122,s211,s212,s221,s222
! model constants 
   double precision escale
   double precision, allocatable, dimension(:) :: pij,xx
   double precision rrk(0:4)
! this is a scaling factor for the entropy
   escale=0.25D0
! Coefficients to calculate the mole fractions from the clusters
! which must be ordered A, A0.75B0.25, A.5B.5, A.25B.75, B
   kk=4; rk=0.25D0
! Tetrahedron in fcc lattice, z=12; m=4.  The factors below are equal to 
! 1/permutations of the clusters. Needed to obtain ideal ordering at high T
   rrk(0)=one; rrk(1)=0.25D0; rrk(2)=one/6.0D0; rrk(3)=0.25; rrk(4)=one
   if(ncon.ne.5) then
! ncon should be 5 for a binary system
      write(*,*)'3X constituents must be 5!',ncon
      gx%bmperr=4399; goto 1000
   endif
! allocations
   allocate(xx(1:2))
   allocate(pij(0:4))
   rtg=globaldata%rgas*tval
! test without any ordering parameters, the system should be ideal ...
!-------------------------------------------------------------------
   do alpha=0,4
! these are the cluster fractions, for higher order systems must be identified
      pij(alpha)=phvar%yfr(alpha+1)
   enddo
! mole fractions, note "ia" is the number of B atoms!!
   xx=zero
   do ia=0,4
      xx(1)=xx(1)+(kk-ia)*rk*pij(ia)
      xx(2)=xx(2)+ia*rk*pij(ia)
   enddo
   if(abs(xx(1)+xx(2)-one).gt.1.0e-8) stop '3X SROT mole fraction error'
! xx(1) is fraction of A, xx(2) fraction of B.  NOT USED!!
   ggg=zero
! This is summing the SRO entropy part
   do ia=0,4
      ggg=ggg+escale*pij(ia)*log(pij(ia)*rrk(ia))
   enddo
! These are the configurational G/RT and S/R
   phvar%gval(1,1)=ggg
   phvar%gval(2,1)=ggg/tval
   do ia=0,4
! d/pij ( x1*ln(x1)+x2*log(x2)+  ...+ pij*log(pij)
! note "ia" counts the B atoms
      phvar%dgval(1,ia+1,1)=escale*(one+log(pij(ia)*rrk(ia)))
      phvar%dgval(2,ia+1,1)=phvar%dgval(1,ia+1,1)/tval
   enddo
!------------------------------------------------------
! second derivatives, symmetric, stored only upper half
! approximate with 1/pij
   jxsym=0
   jb=1
   phvar%d2gval=zero
   do ia=1,ncon
      phvar%d2gval(ixsym(ia,ia),1)=escale/(rrk(ia-1)*pij(ia-1))
   enddo
!-----------------------------------
1000 continue
   return
 end subroutine config_entropy_srot

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
 
!\addtotable subroutine push_pyval
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

!\addtotable subroutine pop_pyval
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

!\addtotable subroutine calc_disfrac
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
!   write(*,*)'3X entering calc_disfrac',lokph,lokcs
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
!   write(*,*)'3X in calc_disfrac',lokph,lokdis,&
!        associated(phord),associated(phdis)
!   write(*,*)'3X Calc disfra: ',lokph,lokcs,lokdis
!   phdis=>ceq%phase_varres(lokdis)
!   call calc_disfrac2(ceq%phase_varres(lokcs)%disfra,&
!   call calc_disfrac2(ceq%phase_varres(lokcs),ceq%phase_varres(lokdis),ceq)
   call calc_disfrac2(phord,phdis,ceq)
1000 continue
   return
 end subroutine calc_disfrac

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine calc_disfrac2
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
!   write(*,*)'3X entering calc_disfrac2'
!   disrec=phord%disfra
!   lokdis=disrec%varreslink
!   phdis=>disrec%phdapointer
! this is the record with the ordered constitution
!   phord=>ceq%phase_varres(lokcs)
! this is a record within the ordered constitution record for disordered fracs
!   write(*,*)'3X entering calc_disfrac2'
   disrec=>phord%disfra
!   write(*,*)'3X in calc_disfrac2 B',associated(disrec),associated(phdis)
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
!   write(*,*)'3X in calc_disfra2c B1: ',associated(phdis%yfr)
! Segmentation fault that phdis%yfr not always allocated !!!
!   write(*,*)'3X in calc_disfra2c B1: ',allocated(phdis%yfr)
   phdis%yfr=zero
!   write(*,*)'3X in calc_disfra2c A1: ',disrec%tnoofyfr
   do is=1,disrec%tnoofyfr
      phdis%yfr(disrec%y2x(is))=&
           phdis%yfr(disrec%y2x(is))+disrec%dxidyj(is)*phord%yfr(is)
!      write(*,77)'3X disfrac 2: ',is,disrec%y2x(is),phdis%yfr(disrec%y2x(is)),&
!           disrec%dxidyj(is),phord%yfr(is)
77    format(a,2i3,3(1pe12.4))
   enddo
!   write(*,*)'3X in calc_disfrac2 A2'
! check if phase is really ordered, meaning that the disordered fractions
! are equal to the ordered ones
   ordered=.false.
   do is=1,disrec%tnoofyfr
      if(abs(phdis%yfr(disrec%y2x(is))-&
           phord%yfr(is)).gt.yminord) ordered=.true.
   enddo
!   write(*,*)'3X calc_disfrac2 A3'
   if(.not.ordered) then
! if this bit set one will not calculate the ordered part of the phase
      phord%status2=ibclr(phord%status2,csorder)
   else
! bit must be cleared as it might have been set at previous call
      phord%status2=ibset(phord%status2,csorder)
   endif
!   write(*,*)'3X in calc_disfrac2 A4: ',phord%status2
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

!\addtotable subroutine disordery
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

!\addtotable subroutine disordery2
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
!   write(*,*)'3X dis2: 1',lokdcs
! this was never assigned!! BOS 16.11.04
   lokdcs=disrec%varreslink
!   write(*,*)'3X lokdcs: ',ceq%eqno,lokdcs,disrec%latd,&
!        allocated(ceq%phase_varres(lokdcs)%yfr)
! here copy: 
! y(ord,1,1)=y(dis,1); y(ord,1,2)=y(dis,2); y(ord,1,3)=y(dis,3); 
! y(ord,2,1)=y(dis,1); y(ord,2,2)=y(dis,2); y(ord,2,3)=y(dis,3); 
!   write(*,*)'3X dis2: 2',disrec%latd,disrec%nooffr(1)
!   write(*,*)'3X dis2: 3',phdis%yfr(1)
!   write(*,*)'3X disordery2: ',lokdcs
   kk=0
! latd is the number of sublattices to be added to first disordered sublattice
   do ll=1,disrec%latd
      do is=1,disrec%nooffr(1)
! the number of constituents in first disordered sublattice same as in ordered
         kk=kk+1
!         phvar%yfr(kk)=phdis%yfr(is)
!         xxx=phdis%yfr(is)
         xxx=ceq%phase_varres(lokdcs)%yfr(is)
! phvar is the phase_varres record of the ordered phase
         phvar%yfr(kk)=xxx
      enddo
   enddo
!   write(*,*)'3X dis2: 4',disrec%ndd
   if(disrec%ndd.eq.2) then
! one can have 2 sets of ordered subl. like (Al,Fe)(Al,Fe)...(C,Va)(C,Va)...
! BUT NEVER TESTED
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

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine uniquac_model
!\begin{verbatim}
 subroutine uniquac_model(moded,ncon,phres,ceq)
! Calculate the Gibbs energy of the UNIQUAC model (Abrams et al 1975)
! Modified 2018/Oct, Nov, Dec
! It returns UNIQUAC G and first and second derivatives of G in phres%gval etc.
! The values of q_i and r_i are be stored in species record, not identifiers
! The residual term should be stored as a UQTAU identifier 
   implicit none
   integer moded,ncon
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_phase_varres), pointer :: phres
!\end{verbatim}
   integer ia,ib,ic,id,ie,jj,nprop,nint,lokph,loksp,ii
   double precision, allocatable, dimension(:) :: theta,phi,qval,xfr,rval
   double precision, allocatable, dimension(:,:) :: tau,dtaudt
   double precision, allocatable, dimension(:) :: rho,kvottau
   double precision, allocatable, dimension(:) :: dgv,d2gv,sumdtaudt,d2gdydt
   double precision hzeta,gres,sumxq,sumxr,dgrdt,xxx,yyy,sumtt
   double precision gc,gr,term1,term2,dgr,dgc,d2gc,d2gr,sumxiqi
! we must have a property index "i" for each tau_ji 
   integer, allocatable, dimension(:) :: unqtau
!   if(moded.lt.2) then
! if moded=/=2 then some of the 2nd derivatives needed is not present ...
! moded=0 when calculating with the gridminimizer
!      write(*,*)'Skipping uniquac phase as no derivatives'
!      goto 1000
!   endif
! need theta = \sum_i q_i*x_i and Phi=\sum_i r_i*x_i
!   write(*,*)'3X in uniquac 1'
   allocate(unqtau(ncon))
   allocate(theta(ncon))
   allocate(phi(ncon))
   allocate(qval(ncon))
   allocate(rval(ncon))
   allocate(xfr(ncon))
   allocate(dgv(ncon))
   allocate(d2gdydt(ncon))
   allocate(tau(ncon,ncon))
   allocate(dtaudt(ncon,ncon))
   allocate(rho(ncon))
   allocate(kvottau(ncon))
   allocate(sumdtaudt(ncon))
! number of interactions
   nint=ncon*(ncon+1)/2
   allocate(d2gv(nint))
! we need some place to store these indices if we have no addition record ...
! UQT is a model parameter identifier with constituent index 
! tau_ji is UQT&I(LIQ,J)
   call need_propertyid('UQT ',id)
   if(gx%bmperr.ne.0) goto 1000
   nprop=phres%listprop(1)-1
!   unqq=0; unqr=0; unqt12=0; unqt21=0
   unqtau=0
   do ia=1,nprop
      if(phres%listprop(ia)/100.eq.id) then
! listprop is 2600+constituent index
! parameter syntax is UQT&UA(LIQUID,UB) for TAU_{UB,UA}
         ib=mod(phres%listprop(ia),100)
         unqtau(ib)=ia
!         write(*,*)'3X found tau: ',phres%listprop(ia),ib,unqtau(ib)
      endif
   enddo
! copy mole fractions to xfr
   xfr=phres%yfr
! extract unqq and unqr from species record
   lokph=phres%phlink
   sumxq=zero
   sumxr=zero
   do ia=1,ncon
! values of q and r for each species is stored in species record
      loksp=phlista(lokph)%constitlist(ia)
      if(btest(splista(loksp)%status,SPUQC)) then
         qval(ia)=splista(loksp)%spextra(1)
         rval(ia)=splista(loksp)%spextra(2)
      else
         qval(ia)=one
         rval(ia)=one
      endif
! calculate the sum of q and r
      sumxq=sumxq+xfr(ia)*qval(ia)
      sumxr=sumxr+xfr(ia)*rval(ia)
   enddo
! extracting residual parameters complicated ... initiate to zero
   tau=one
!
! df/dx_k = ln(phi_k) + 1 - phi_k/x_k + 
!           z/2 q_k ( ln(theta_k/phi_k) + phi_k/theta_k-1 )
!
! theta = UQQ = x_i*q_i(\sum_j q_j*x_j) 
! Phi= UQR=x_i*r_i*(\sum_i r_j*x_j)
! write(*,*)'3 X Calculate Phi, theta and some invariants for the residual term'
   dtaudt=zero
   sumdtaudt=zero
   do ia=1,ncon
      theta(ia)=xfr(ia)*qval(ia)/sumxq
      phi(ia)=xfr(ia)*rval(ia)/sumxr
!      write(*,'(a,i4,6(F7.4))')'3X theta and phi:',ia,xfr(ia),theta(ia),phi(ia)
!----------------- residual term tau_ji, may not be present!
!      write(*,'(a,2i3,3(1pe12.4))')'3X tau1: ',ia,unqtau(ia),&
!           phres%gval(1,unqtau(ia)),phres%dgval(1,3-ia,unqtau(ia))
      if(unqtau(ia).eq.0) then
! OK if zero, this means no residual parameter
         continue
      else
! there are some residual parameters, extract their values 
! MODIFIED xfr(ib)*tau_(ib,ia) stored in phres%gval(1,unqtau(ia))
! MODIFIED xfr(ia)*tau_(ia,ib) stored in phres%gval(1,unqtau(ib))
! By JING: here you need to be careful!!!
         tauloop2: do ib=1,ncon
            if(ib.eq.ia) cycle tauloop2
! This is the derivative wrt xfr(ib) of UQT&IA i.e. TAU_IB,IA
            term1=phres%dgval(1,ib,unqtau(ia))
! NOTE, default value one set above, value must be larger than zero
            if(term1.ne.zero) then
               tau(ib,ia)=term1
! The derivative wrt T is in phres%dgval(2,ib,unqtau(ia))
               dtaudt(ib,ia)=phres%dgval(2,ib,unqtau(ia))
            endif
!            write(*,'(a,3i3,4(1pe12.4))')'3X tau2: ',ia,ib,unqtau(ia),&
!                 phres%dgval(1,ib,unqtau(ia)),phres%dgval(2,ib,unqtau(ia))
         enddo tauloop2
      endif
   enddo
!-----------------------
!   write(*,10)'3X q: ',qval,phres%gval(1,unqq)
!   write(*,10)'3X UNIQUAC theta: ',xfr,theta
!   write(*,10)'3X UNIQUAC tau: ',xfr,tau
! default value of tau is unity
! OK   write(*,10)'3X tau: ',tau
10 format(a,6(1pe12.4))
! here we calculate \sum_i \theta_i \tau_{ij} stored in rho
   do ia=1,ncon
      term1=zero
      term2=zero
      do ib=1,ncon
! this is \sum_j \theta_j \tau_{ji}         
         term1=term1+theta(ib)*tau(ia,ib)
         term2=term2+theta(ib)*dtaudt(ia,ib)
! OK         write(*,'(a,2i3,6(1pe12.4))')'3X rho1: ',ia,ib,theta(ib),tau(ia,ib)
      enddo
! these values are \sum_j \theta_j\tau_ji (and the T-derivative)
      rho(ia)=term1
      sumdtaudt(ia)=term2
   enddo
! OK   write(*,10)'3X rho: ',rho
!   gx%bmperr=4399; goto 1000
   do ia=1,ncon
! need for the residual derivatives \sum_i (\theta_i \tau_ji)/\rho_i
      term1=zero
      do ib=1,ncon
! I am never sure if it should be tau(ia,ib) or tau(ib,ia) ...
         term1=term1+theta(ib)*tau(ib,ia)/rho(ib)
      enddo
! This is \sum_i (\theta_i \tau_ki)/\rho_i where index "ia" is subscript "k" 
      kvottau(ia)=term1
   enddo
!   write(*,107)'3X tau:: ',kvottau,zero,xfr(1),ceq%tpval(1)
!   gx%bmperr=4399; goto 1000
! This is z/2
   hzeta=5.0D0
! Here the UNIQUAC GIBBS ENERGY and derivatives are calculated.  
! phres%gval has the ideal configurational enntropy already
! and possibly any reference energy terms!
   gc=zero; gr=zero; dgrdt=zero
   gmloop: do ia=1,ncon
! The residual and configurational G
! ALL OK without residual term, rho_i is rho_i in abrams1-190107.pdf
      gr=gr-xfr(ia)*qval(ia)*log(rho(ia))
! NOTE gc=zero if all qval and rval equal for all components !!
      gc=gc+xfr(ia)*log(phi(ia)/xfr(ia))+&
           hzeta*xfr(ia)*qval(ia)*log(theta(ia)/phi(ia))
!      write(*,210)'3X gci: ',ia,xfr(ia),phi(ia),xfr(ia)*log(phi(ia)/xfr(ia)),&
!           qval(ia),theta(ia),hzeta*xfr(ia)*qval(ia)*log(theta(ia)/phi(ia))
210   format(a,i2,6(1pe10.2))
! The first T-derivative of residual (as in Abrams1.pdf 190105)
      dgrdt=dgrdt-xfr(ia)*qval(ia)*sumdtaudt(ia)/rho(ia)
! The second T-derivative ... NOT USED YET!
!      d2grdt2=....
   enddo gmloop
!   write(*,211)'3X GC: ',0,gr,gc,gr+gc
211   format(a,i2,6(1pe12.4))
!   gx%bmperr=4399; goto 1000
   dgv=zero
   d2gdydt=zero
   first: do ib=1,ncon
! first residual derivative with respect to component 
! kvottau was calculated above 
      dgr=qval(ib)*(one-log(rho(ib))-kvottau(ib))
! derivative with respect to T and xfr(b)
! we must sum_i theta_i/rho_i (dtau_ki/dT)
      xxx=zero; yyy=zero
      sumtt=zero
! second derivative of residual
      do ia=1,ncon
         xxx=xxx+theta(ia)*dtaudt(ib,ia)/rho(ia)
         yyy=yyy+theta(ia)*tau(ib,ia)*sumdtaudt(ia)/rho(ia)**2
      enddo
      d2gdydt(ib)=-qval(ib)*(sumdtaudt(ib)/rho(ib)+xxx-yyy)
! first derivative with respect to ib of configuration
! NOTE: 1+ln(x) already calculated, thus -log(xfr(ib))
      dgc=log(phi(ib)/xfr(ib))+one-phi(ib)/xfr(ib)+&
           hzeta*qval(ib)*(log(theta(ib)/phi(ib))-one+phi(ib)/theta(ib))
!      write(*,'(a,i3,f10.6,1pe12.4)')'3X ib xfr dgc: ',ib,xfr(ib),dgc
      dgv(ib)=dgr+dgc
!      write(*,212)'3X dgr: ',ib,qval(ib),rho(ib),kvottau(ib),dgc,dgv(ib)
212   format(a,i3,6(1pe12.4))
      second: do ic=ib,ncon
! second derivative of configuration with respect to ib and ic
! APPROXIMATE not corrected!!
         do ii=1,ncon
            sumtt=sumtt+theta(ii)*tau(ib,ii)*tau(ic,ii)/rho(ii)**2
         enddo
         d2gr=qval(ib)*qval(ic)/sumxq*(1-tau(ic,ib)/rho(ib)-&
              tau(ib,ic)/rho(ic)+sumtt)
         if(ic.eq.ib) then
            d2gc=-2.0D0*phi(ic)/xfr(ic)**2
!         else
!            d2gc=zero
!            d2gr=zero
         endif
         d2gv(ixsym(ib,ic))=d2gr+d2gc
      enddo second
! VERY APPROXIMATE SECOND DERIVATIVES 
   enddo first
!   write(*,300)'3X UQG: ',gr,gc,(dgv(ia),ia=1,ncon)
!   do ib=1,ncon
!      write(*,300)'3X D2UQG: ',(d2gv(ixsym(ia,ib)),ia=1,ncon)
!   enddo
! copy results to global arrays
! phres%gval(1,1) is Gm, %gval(2,1) is dG/dT, %gval(3,1) is dG/dP, 
!      %gval(4,1) is d2G/dT2 ...
! IMPORTANT the ideal configurational entropy is in %gval(1,1) and %gval(2,1)
! phres%dgval(1,j,1) is dG/dx_j, phres%dgval(2,j,1) is d2G/dTdx_j ... 
! phres%d2gval(ixsym(j,k),1) is d2G/dx_jdx_k stored as upper triangle
! all values divided by RT
! phres%gval(2,1)= no T dependence
! phres%gval(3,1)= no P dependence
! phres%d2gval(ixsym(j,k),1) is d2G/dx_j/dx_k
!   write(*,300)'3X G/RT: ',phres%gval(1,1),gc,phres%gval(1,1)+gc
!   write(*,'(a,i3,5(1pe12.4))')'3X UQG: ',moded,gr,gc,phres%gval(1,1),&
!        gr+gc+phres%gval(1,1)
   phres%gval(1,1)=phres%gval(1,1)+gc+gr
! add dgc/dT and T-dependence of gr.  NOTE gr is also multiplied with RT
   phres%gval(2,1)=phres%gval(2,1)+(gc+gr)/ceq%tpval(1)+dgrdt
   term1=phres%gval(2,1)
!   write(*,'(a,i2,6(1pe12.4))')'3X dG/dT:    ',phres%phtupx,term1,&
!        phres%gval(2,1),gc/ceq%tpval(1),gr/ceq%tpval(1),dgrdt
!   phres%gval(4,1)=phres%gval(2,1)+other terms T-dependent terms (Cp)
   do ia=1,ncon
      phres%dgval(1,ia,1)=phres%dgval(1,ia,1)+dgv(ia)
!      write(*,212)'3X ddy: ',ia,phres%dgval(1,ia,1),dgv(ia)
! The T-dependence of the residual term affects d2G/dydT
!      term1=phres%dgval(2,ia,1)
!      phres%dgval(2,ia,1)=phres%dgval(2,ia,1)+d2gdydt(ia)
!      write(*,'(a,i2,6(1pe12.4))')'3X d2G/dydT: ',ia,term1,&
!           phres%dgval(2,ia,1),d2gdydt
! 2nd derivatives ************ skip for the moment
!      do ib=ia,ncon
!         phres%d2gval(ixsym(ia,ib),1)=phres%d2gval(ixsym(ia,ib),1)+&
!              d2gv(ixsym(ia,ib))
!         write(*,'(a,2i3,6(1pe12.4))')'3X d2g:',ib,ic,d2gr,d2gc,&
!              phres%d2gval(ixsym(ib,ic),1)
!      enddo
   enddo
! check chemical potentials
   xxx=phres%gval(1,1)-xfr(1)*phres%dgval(1,1,1)-xfr(2)*phres%dgval(1,2,1)
!   write(*,'(a,4(1pe12.4))')'3X mu: ',xxx,xxx+phres%dgval(1,1,1),&
!        xxx+phres%dgval(1,2,1)
!   write(*,300)'3X Gm, dG/dx: ',phres%gval(1,1),(phres%dgval(1,ia,1),ia=1,ncon)
300 format(a,6(1pe12.4))
!   do ia=1,ncon
!      write(*,300)'3X d2G/dx2:   ',(phres%d2gval(ixsym(ia,ib),1),ib=1,ncon)
!   enddo
1000 continue
   return
 end subroutine uniquac_model

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine set_driving_force
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

!\addtotable subroutine extract_massbalcond
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
   integer nox,now,nc,jl,iref,iunit,ip,idf,ie,more,numberest,istv,localerr,zz
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
! BUT WE HAVE TO BE CAREFUL IF THIS MUST NOT BE EVALUATED!!
         if(btest(svflista(current%symlink1)%status,SVFVAL)) then
            xxx=ceq%svfunres(current%symlink1)
         else
            xxx=evaluate_svfun_old(current%symlink1,actual_arg,1,ceq)
         endif
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
! this is W(index1)=value, CHECK UNIT if %!!!  end x2
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
!               write(*,*)'3X Missing condition for two elements.'
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

!\addtotable subroutine save_constitutions
!\begin{verbatim}
 subroutine save_constitutions(ceq,copyofconst)
! copy the current phase amounts and constituitions to be restored
! if calculations fails during step/map
! DANGEROUS IF NEW COMPOSITION SETS CREATED
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision, allocatable, dimension(:) :: copyofconst
!\end{verbatim} %+
   integer varresx,nz,ij,syfr,allsize
! calculate dimension of copyofconst
   nz=0
! skippa varres with index 1, that is the reference phase
!   do varresx=2,csfree-1
   do varresx=2,highcs
      if(allocated(ceq%phase_varres(varresx)%yfr)) then
! NOTE size( ... ) can return reasonable value even if not allocated !!!
! BUT why is phas_varres(varresx)%yfr it not allocated ???
! evidently the composition set for varresx is created ... maybe removed??
         syfr=size(ceq%phase_varres(varresx)%yfr)+1
      else
         syfr=1
      endif
!      write(*,12)'3X Varres record and size: ',varresx,1+syfr,nz
12    format(a,3i5)
      nz=nz+2+syfr
   enddo
   allsize=nz+1
!   write(*,*)'3X In save_constitution',highcs,allsize
   allocate(copyofconst(allsize))
! modification due to problems, save allocated size in first word
   copyofconst(1)=allsize
   nz=2
!   do varresx=2,csfree-1
   do varresx=2,highcs
! save 1+syfr values for each composition set
! segmentation fault in this loop for stepbug (>20 elements COST507)
! crash happends when higher composition sets are stored ...
! SAVE also the amount of the phase, DGM and the size of yfr!!
      copyofconst(nz)=ceq%phase_varres(varresx)%amfu
      copyofconst(nz+1)=ceq%phase_varres(varresx)%dgm
! varresx is 1 higher than phase index
!      if(copyofconst(nz).gt.zero) &
!           write(*,*)'3X saving amount: ',varresx-1,nz,copyofconst(nz)
      nz=nz+1
      if(allocated(ceq%phase_varres(varresx)%yfr)) then
         syfr=size(ceq%phase_varres(varresx)%yfr)
      else
!         write(*,*)'3X no fractions for: ',varresx-1,nz+1
         syfr=0
      endif
!      write(*,16)'3X Storing varres record: ',varresx,syfr,size(copyofconst),nz
16    format(a,5i5)
! the segmentation fault seems not to be the allocation of copyofconst but
! rather that we cannot access the yfr in ceq%phase_varres(varresx)
! for the extra composition sets created by gridmin
      copyofconst(nz+1)=syfr
      nz=nz+1
      do ij=1,syfr
         copyofconst(nz+ij)=ceq%phase_varres(varresx)%yfr(ij)
      enddo
      nz=nz+1+syfr
   enddo
!   write(*,*)'3x saved size in word 1: ',highcs,allsize,nz-1
1000 continue
   return
 end subroutine save_constitutions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine restore_constitutions
!\begin{verbatim} %-
 subroutine restore_constitutions(ceq,copyofconst)
! restore the phase amounts and constitutions from copyofconst
! if calculations fails during step/map
! DANGEROUS IF NEW COMPOSITION SETS CREATED
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision copyofconst(*)
!\end{verbatim}
   integer nz,varresx,ij,syfr,savedsyfr,sizeofcopy
! size of copyofconst in first word
   sizeofcopy=int(copyofconst(1))
!   write(*,*)'3X restoring amounts: ',highcs,sizeofcopy
! skippa varres with index 1, that is the reference phase
!   do varresx=2,csfree-1
   nz=2
   do varresx=2,highcs
! note varresx is index of phase_varres, always 1 bigger than phase index
!      if(copyofconst(nz).gt.zero) &
!           write(*,*)'3X restore amount: ',varresx-1,nz,copyofconst(nz)
      ceq%phase_varres(varresx)%amfu=copyofconst(nz)
      ceq%phase_varres(varresx)%dgm=copyofconst(nz+1)
      if(allocated(ceq%phase_varres(varresx)%yfr)) then
         syfr=size(ceq%phase_varres(varresx)%yfr)
      else
         syfr=0
      endif
! fraction records may have been allocated!! use saved syfr
      nz=nz+2
      savedsyfr=int(copyofconst(nz))
      if(savedsyfr.eq.0 .or. savedsyfr.ne.syfr) then
         ceq%phase_varres(varresx)%dgm=-one
!         write(*,12)'Restore saved size for phase: ',varresx-1,nz-2,syfr,&
!              int(copyofconst(nz)),copyofconst(nz-2),&
!              ceq%phase_varres(varresx)%dgm
12       format(a,4i5,2(1pe12.4))
         syfr=savedsyfr
      endif
      do ij=1,syfr
         ceq%phase_varres(varresx)%yfr(ij)=copyofconst(nz+ij)
      enddo
!      write(*,17)varresx-1,nz,syfr,ceq%phase_varres(varresx)%amfu,&
!           (ceq%phase_varres(varresx)%yfr(ij),ij=1,syfr)
17    format('3X r:',i2,2i3,6(1pe12.4))
      nz=nz+1+syfr
      if(nz-1.gt.sizeofcopy) write(*,*)'3X problem restore:',varresx,nz
   enddo
1000 continue
!   gx%bmperr=4399
   return
 end subroutine restore_constitutions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine save_phase_constitutions
!\begin{verbatim}
 subroutine save_phase_constitutions(rw,ceq,copyofconst)
! copy the current phase amounts and constituitions to be restored
! trying to fix problems with saving invariants 
! compared to reoutines above here abnorm is also saved ...
! rw=0 if save, 1 if restore
! NOTE different ceq may be used for save and restore!
   implicit none
   integer rw
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision, allocatable, dimension(:) :: copyofconst
!\end{verbatim} %+
   integer varresx,nz,ij,syfr,allsize,savedsyfr,sizeofcopy
   if(rw.eq.0) then
! calculate dimension of copyofconst
      nz=0
! Calculate space needed
! All phases saved idependent of status
! skippa varres with index 1, that is the reference phase
      do varresx=2,highcs
         if(allocated(ceq%phase_varres(varresx)%yfr)) then
            syfr=size(ceq%phase_varres(varresx)%yfr)+1
         else
            syfr=1
         endif
! ionic liquid model require more data saved (see set_constitution)
         ij=ceq%phase_varres(varresx)%phlink
         if(btest(phlista(ij)%status1,PHIONLIQ)) then
            write(*,*)'3X cannot save ionic liquid constitutions'
            gx%bmperr=4399; goto 1000
         endif
! we should save 5 reals in addition to the fractions
         nz=nz+5+syfr
      enddo
      allsize=nz+2
      allocate(copyofconst(allsize))
! modification due to problems, save allocated size in first word
      copyofconst(1)=allsize
      copyofconst(2)=highcs
      nz=3
      do varresx=2,highcs
! save 1+syfr values for each composition set
! SAVE also the amount of the phase, DGM and the size of yfr!!
         copyofconst(nz)=ceq%phase_varres(varresx)%amfu
         copyofconst(nz+1)=ceq%phase_varres(varresx)%abnorm(1)
         copyofconst(nz+2)=ceq%phase_varres(varresx)%abnorm(2)
         copyofconst(nz+3)=ceq%phase_varres(varresx)%abnorm(3)
         copyofconst(nz+4)=ceq%phase_varres(varresx)%dgm
         if(allocated(ceq%phase_varres(varresx)%yfr)) then
            syfr=size(ceq%phase_varres(varresx)%yfr)
         else
            syfr=0
         endif
         copyofconst(nz+5)=syfr
         nz=nz+5
         do ij=1,syfr
            copyofconst(nz+ij)=ceq%phase_varres(varresx)%yfr(ij)
         enddo
         nz=nz+syfr+1
      enddo
      write(*,*)'3X saved constitution: ',allsize,nz
   else
! restore saved amounts and fractions      
      if(.not.allocated(copyofconst)) then
         write(*,*)'3X no constitutions saved!'
         gx%bmperr=4399; goto 1000
      endif
      sizeofcopy=int(copyofconst(1))
      if(copyofconst(2).ne.highcs) then
         write(*,*)'3X number of phase tuples not the same'
         gx%bmperr=4399; goto 1000
      endif
      nz=3
      do varresx=2,highcs
! note varresx is index of phase_varres, always 1 bigger than phase index
         ceq%phase_varres(varresx)%amfu=copyofconst(nz)
         ceq%phase_varres(varresx)%abnorm(1)=copyofconst(nz+1)
         ceq%phase_varres(varresx)%abnorm(2)=copyofconst(nz+2)
         ceq%phase_varres(varresx)%abnorm(3)=copyofconst(nz+3)
         ceq%phase_varres(varresx)%dgm=copyofconst(nz+4)
         if(allocated(ceq%phase_varres(varresx)%yfr)) then
            syfr=size(ceq%phase_varres(varresx)%yfr)
         else
            syfr=0
         endif
! fraction records may have been allocated!! use saved syfr
         nz=nz+5
         savedsyfr=int(copyofconst(nz))
         if(savedsyfr.eq.0 .or. savedsyfr.ne.syfr) then
            write(*,*)'3X phase with zero saved fractions'
            ceq%phase_varres(varresx)%dgm=-one
            syfr=savedsyfr
         endif
         do ij=1,syfr
            ceq%phase_varres(varresx)%yfr(ij)=copyofconst(nz+ij)
         enddo
         nz=nz+1+syfr
      enddo
      if(nz-1.gt.sizeofcopy) write(*,*)'3X problem restore:',sizeofcopy,nz
   endif
1000 continue
   return
 end subroutine save_phase_constitutions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine calc_eec_gibbsenergy
!\begin{verbatim}
 subroutine calc_eec_gibbsenergy(phres,ceq)
! calculate an ideal Gibbs energy with just configurational entropy
! phres is pointer to phase_varres record for the phase
! for a solid phase with higher entropy than the liquid
! G = RT \sum_s a_s \sum_i y_si \ln(y_si)
! dG/dy_si = RT a_s (1+ln(y_si))
! d2G/dy_si^2 = RT a_s/y_si            all other 2nd derivatives zero   
   implicit none
   TYPE(gtp_phase_varres), pointer :: phres
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,sl,i1,i2,kk
   double precision sconf,tval,as
! this is the index of the phase in phlista with phase structure
   lokph=phres%phlink
! zero all second derivatives of G, the diagonal added below
   kk=phlista(lokph)%tnooffr*(phlista(lokph)%tnooffr+1)/2
   do i1=1,kk
      phres%d2gval(1,i1)=zero
   enddo
   kk=0
   sconf=zero
   tval=ceq%tpval(1)
   do sl=1,phlista(lokph)%noofsubl
      as=phres%sites(sl)
      do i1=1,phlista(lokph)%nooffr(sl)
         kk=kk+1
         sconf=sconf+as*phres%yfr(kk)
         phres%dgval(1,kk,1)=as*(one+log(phres%yfr(kk)))
         phres%dgval(2,kk,1)=phres%dgval(2,kk,1)/tval
         phres%d2gval(kxsym(kk,kk),1)=as/phres%yfr(kk)
      enddo
   enddo
! return values divided by RT
   phres%gval(1,1)=sconf
1000 continue
 end subroutine calc_eec_gibbsenergy

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine setendmemarr
!\begin{verbatim}
 subroutine setendmemarr(lokph,ceq)
! stores the pointers to all ordered and disordered endmemners in arrays
! intended to allow parallel calculation of parameters
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

!\addtotable subroutine tabder
!\begin{verbatim}
 subroutine tabder(iph,ics,times,ceq)
! tabulate derivatives of phase iph with current constitution and T and P
   implicit none
   integer iph,ics,times
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   character name*24
   double precision kappa,napfu,t,p,rtg,g,v,s,h,u,f,cp,alpha,cpu1,cpu2
   integer tnk,lokph,nsl,lokres,lokcs,ll,ll2,kk1,kk2,kk3,kk4,loksp
!
! For time measurements
!   lokph=len(name)
   lokph=phases(iph)
   nsl=phlista(lokph)%noofsubl
! calculate G and derivatives, lokres returns index of phase_varres
   call cpu_time(cpu1)
   do loksp=1,times
      call calcg(iph,ics,2,lokres,ceq)
      if(gx%bmperr.ne.0) then
         goto 1000
      endif
   enddo
   call cpu_time(cpu2)
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
         'Gibbs energy J/FU  ',28('.'),1Pe16.8,e16.7)
   write(kou,102)F,F/rtg,H,H/rtg,U,U/rtg,S,S/rtg,V,V/rtg,&
        CP,CP/rtg,alpha,alpha/rtg,kappa,kappa/rtg
102 format('Helmholtz energy J/FU  ',24('.'),1PE16.8,e16.7 &
        /'Enthalpy J/FU  ',32('.'),1PE16.8,e16.7 &
        /'Internal energy J/FU  ',25('.'),1PE16.8,e16.7 &
        /'Entropy J/FU/K  ',31('.'),1PE16.8,e16.7 &
        /'Volume m3/FU ',34('.'),1PE16.8,e16.7 &
        /'Heat capacity J/FU/K  ',25('.'),1PE16.8,e16.7 &
        /'Thermal expansion 1/K ',25('.'),1PE16.8,e16.7 &
        /'Bulk modulus m2/N ',29('.'),1PE16.8,e16.7)
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
   if(times.gt.1) then
      write(*,11)times,cpu2-cpu1,1.0D6*(cpu2-cpu1)/dble(times)
11    format('CPU times for ',i6,' calculations: ',1pe15.7,' s, ',1pe15.7,' ms')
   endif
!    write(*,*)'3X tabder 7B: ',kk2
!    write(*,*)'3X tabder: ',rtg,rtg*phase_varres(lokcs)%gval(1,1)
1000 continue
   return
 end subroutine tabder

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!

