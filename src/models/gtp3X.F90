!
! gtp3X included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     15. Calculate things
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
!   write(*,*)'3X calcg: ',lokcs
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
!   if(ocv()) write(*,*)'3X in gcalc_internal: ',lokph
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
   elseif(mqmqa) then
! MQMQA FactSage entropy model
!      write(*,*)'3X calling MQMQA liquid entropy model'
      call config_entropy_mqmqa(moded,phlista(lokph)%tnooffr,phres,&
           phlista(lokph),gz%tpv(1))
! when we come back mqmqaf should have some arrays allocated ....
! DO NOT SET THIS POINTER BEFORE ARRAYS ARE ALLOCATED IN CONFIG_ENTROPY_MQMQA
      mqf=>phres%mqmqaf
!     write(*,'(a,1pe12.4)')'3X MQMQA cfgG:',phres%gval(1,1)*gz%rgast*phres%amfu
!      write(*,'(a,1pe12.4)')'3X MQMQA cfgG:',phres%gval(1,1)*gz%rgast
!      do mqmqj=1,mqf%npair
!         write(*,777)mqf%pair(mqmqj),(mqf%dpair(mqmqj,kk),kk=1,mqf%nconst)
!         write(*,777)mqf%pair(mqmqj),mqf%dpair(mqmqj,1),mqf%dpair(mqmqj,2),&
!              mqf%dpair(mqmqj,3)
777      format('3X pairs:',F10.6,2x,10F9.5)
!      enddo
   elseif(btest(phlista(lokph)%status1,PHSSRO)) then
! Simple short range order entropy model
      write(*,*)'3X calling SSRO liquid model'
      call config_entropy_ssro(moded,lokph,phres,gz%tpv(1))
   else
!----------- the CF Bragg-Williams ideal configurational entropy per sublattice
! NOTE: for phases with disordered fraction set this is calculated
! ONLY for the original constituent fraction set with ordering sublattices
      call config_entropy(moded,nsl,phlista(lokph)%nooffr,phres,gz%tpv(1))
   endif
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3X segmentation fault 10'
!-------------------------------------------------------------------
! MQMQA separate calculation of G as well as entropy!!!
! NO no 2state model, einstein, magnetism etc.
   if(mqmqa) then
!      write(*,*)'3X Separate routine for MQMQA'
      call calc_mqmqa(lokph,phres,mqf,ceq)
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
                     write(*,*)'3X Toop/Kohler and permutations illegal'
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
! this iloop3 is for all components "id".  Normally binay interactions depend
! only on the constituents gz%iq(1) and gz%iq(2) but Toop/Kohler method depend
! also on other constituents!  That is taken care of inside this IF
        write(*,'(a,i3,5(1pe12.4))')'3X toop dvals',&
             id,pyq,dvals(1,id),phres%dgval(1,id,ipy)+pyq*dvals(1,id),&
             (phres%dgval(1,id,ipy)-pyq*dvals(1,id))*gz%rgast
        do itp=1,3
           phres%dgval(itp,id,ipy)=phres%dgval(itp,id,ipy)-pyq*dvals(itp,id)
        enddo
! ignore contribution to the second derivatives phres%d2gval
     endif toop7
                        do itp=1,3
                           phres%dgval(itp,id,ipy)=&
                                phres%dgval(itp,id,ipy)+ &
                                dpyq(id)*vals(itp)
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

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

 !\addtotable subroutine calc_mqmqa
!\begin{verbatim}
 subroutine calc_mqmqa(lokph,phres,mqf,ceq)
! Called from calcg_internal, calculates G for the mqmqa phase
! separate subroutine for the entropy which calculates all data in phres%mqf
   integer lokph
   type(gtp_phase_varres), pointer :: phres
   type(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_mqmqa_var), pointer :: mqf
!\end{verbatim}
! Most variables here are the same as in calcg_internal ...
   integer, parameter :: f1=50
   integer mqmqj,kend,s1,s2,s3,id,nofc2,ipy,lokfun,typty,itp,zp
   double precision vals(6),pyq,rtg,aff
   double precision, dimension(:), allocatable :: dpyq(:),d2pyq(:),d2vals(:)
   double precision, dimension(:,:), allocatable :: dvals(:,:),affarr(:)
! for saving FNN reference energies
   double precision refg(f1,f1)
   double precision dummy1,dummy2
   TYPE(gtp_parcalc) :: gz
   TYPE(gtp_property), pointer :: proprec
   TYPE(gtp_endmember), pointer :: endmemrec
   TYPE(gtp_interaction), pointer :: intrec
   TYPE(gtp_pystack), pointer :: pystack
   TYPE(gtp_phase_add), pointer :: addrec
!------------------------------------- 
! allocate arrays
!   write(*,*)'3X in calc_mqmqa'
   gz%nofc=phlista(lokph)%tnooffr
   nofc2=gz%nofc*(gz%nofc+1)/2
   allocate(dpyq(gz%nofc))
   allocate(d2pyq(nofc2))
   allocate(dvals(3,gz%nofc))
   allocate(d2vals(nofc2))
   allocate(affarr(mqf%npair))
   affarr=zero
   nullify(pystack)
   rtg=globaldata%rgas*ceq%tpval(1)
   mqf=>phres%mqmqaf
!   refg=zero
   dummy2=zero
! list %pp
! %pp( quad , FNN index )
!   do mqmqj=1,mqmqa_data%nconst
!      write(*,17)'3X %pp: ',mqmqj,(mqmqa_data%pp(s1,mqmqj),s1=1,4)
!   enddo
17 format(a,i3,4(1pe12.4))
!--------------------------------------
! first loop over ALL endmembers
   mqmqj=0
   endmemrec=>phlista(lokph)%ordered
! This should be number of atoms for scaling G
!   dummy1=phres%abnorm(1)/rtg       this was OK before ...
   dummy1=one/rtg
! %amfu * %abnorm(1) is number of moles in the liquid
! in the test case we have 6 atoms in the liquid phase
!   dummy1=6.0D0/rtg
!   dummy1=one/(phres%abnorm(1)*rtg)
!   write(*,'(a,3(1pe14.6))')'3X mqmqa scaling: ',dummy1,&
!        phres%amfu,phres%abnorm(1)
! This first loop: all endmember parameters
! this can give SRO contribution and excess from SNN parameters
! or it makes it possible to calculate the G for the FNN parameters
   endmemloop1: do while(associated(endmemrec))
      mqmqj=mqmqj+1
      if(mqmqj.gt.mqmqa_data%nconst) exit endmemloop1
      kend=mqmqa_data%contyp(5,mqmqj)
      if(kend.le.0) then
! This is an SNN parameter we calculate and add SNN energy and interactions ...
!         write(*,*)'3X SNN endmember record found',mqmqj
         proprec=>endmemrec%propointer
         mqsnn: do while(associated(proprec))
! This loop is not really necessay, in mqmqa the only property is G at present
            typty=proprec%proptype
            if(typty.ne.1) stop '3X illegal typty in mqmqa model'
            ipy=1
            lokfun=proprec%degreelink(0)
            call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
            if(gx%bmperr.ne.0) goto 1000
!            write(*,'(a,6(1Pe12.4))')'3X vals1:',vals
            if(ipy.eq.1) then
               vals=vals*dummy1
! This is an SNN ordering parameter, reference state addel in second loop
            endif
            pyq=phres%yfr(mqmqj)
! Should I use any factor??
!         aff=mqmqa_data%pp(1,mqmqj)
            aff=one
! NOTE the reference state contribution to this SNN added in next loop
! for all quads!!
            do itp=1,3
               phres%dgval(itp,mqmqj,ipy)=phres%dgval(itp,mqmqj,ipy)+vals(itp)
            enddo
! Initially ignore 2nd derivatives, d2G/dy2=1/y set by entropy calculation
! ipy is property, ipy=1 means G, ipy=2 means Curie T etc.
! %gval(1,1) is total G, %gval(2,1) is total dG/dT  etc.
            do itp=1,6
               phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*vals(itp)
            enddo
!            write(*,210)'3X SRO G, dG/dqi: ',mqmqj,mqmqj,pyq,aff,&
!                 phres%gval(1,1),(phres%dgval(1,s1,1),s1=1,gz%nofc)
            proprec=>proprec%nextpr
         enddo mqsnn
!600      continue
!         write(*,*)'3X any excess parameters will be handled in 3rd loop'
         endmemrec=>endmemrec%nextem
         cycle endmemloop1
      endif
! This is an FNN parameter, we calculate and save the value for later use
      proprec=>endmemrec%propointer
      aff=one/mqmqa_data%pp(1,mqmqj)
      mq1: do while(associated(proprec))
         typty=proprec%proptype
         if(typty.ne.1) stop 'illegal typty in mqmqa model'
         ipy=1
         lokfun=proprec%degreelink(0)
         call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
         if(gx%bmperr.ne.0) goto 1000
!         write(*,'(a,i3,F7.4,3(1Pe10.2))')'3X refg:',mqmqj,aff,vals(1),vals(2)
! we should divide this by the aff of this pair and we will multiply this
! FNN same aff but SNN fractions linking to this pair uase another aff
         if(ipy.eq.1) then
            vals=vals*dummy1*aff
! save values of reference state for use with SNN parameters ??
! kend is FNN (pair) index 
            do s1=1,6
               refg(kend,s1)=vals(s1)
            enddo
         endif
! next property record (should not be any ...)
         proprec=>proprec%nextpr
         if(associated(proprec)) then
            write(*,*)'3X Warning: ignoring second mqmqa property recotd!'
         endif
!         write(*,200)'3X FNN G, dG/dqi: ',phres%gval(1,1),&
!              (phres%dgval(1,s1,1),s1=1,gz%nofc)
200      format(a,1pe12.4,2x,6(1pe12.4))
      enddo mq1
      endmemrec=>endmemrec%nextem
   enddo endmemloop1
!--------------------------------------------------- end first endmember loop
! second loop over all quads, ignore endmember records
! but add reference state parameters to SNN energies
   ipy=1
   qloop: do mqmqj=1,gz%nofc
! this is quad fraction, multiply with all FNN reference energies
      pyq=phres%yfr(mqmqj)
      zp=mqmqa_data%contyp(5,mqmqj)
      pair: if(zp.gt.0) then
! this is an FNN  pair, reference energy in refg(zp,1..6), only one y derivative
! %pp(1..4,mqmqj) is stoichiometric factor for the pair
         aff=mqmqa_data%pp(1,mqmqj)
         do itp=1,3
            phres%dgval(itp,mqmqj,ipy)=phres%dgval(itp,mqmqj,ipy)+&
                 aff*refg(zp,itp)
         enddo
! Initially ignore 2nd derivatives, d2G/dy2=1/y set by entropy calculation
         do itp=1,6
            phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*aff*refg(zp,itp)
         enddo
!         write(*,205)'3X FNN: qix, FNN, aff, pyq, fun, DG: ',mqmqj,zp,aff,&
!              pyq,refg(zp,1),pyq*aff*refg(zp,1)
205      format(a,2i3,F8.5,2x,3(1pe12.4))
!         write(*,210)'3X FNN G, dG/dqi: ',mqmqj,mqmqj,pyq,aff,phres%gval(1,1),&
!              (phres%dgval(1,s1,1),s1=1,gz%nofc)
210      format(a,2i3,2F8.5,1pe12.4,2x,6(1pe10.2))
      else
! this is an SNN with two or more pairs
! For each SNN pair add the contribution to the reference state
!         write(*,'(a,i3,4F10.7)')'3X pp: ',mqmqj,&
!              (mqmqa_data%pp(s1,mqmqj),s1=1,4)
         snnloop: do s1=6,9
! zp is index to an FNN record, there can be 2 or 4 FNN records
            zp=mqmqa_data%contyp(s1,mqmqj)
            if(zp.eq.0) exit snnloop
! %pp(1..4,mqmqj) is stoichiometric factor for the pair
            aff=mqmqa_data%pp(s1-5,mqmqj)
!          write(*,212)zp,mqmqj,ipy,phres%dgval(1,mqmqj,ipy),aff,aff*refg(zp,1)
212         format('3X adding reference state to SNN parameter:',3i3,3(1pe12.4))
            do itp=1,3
               phres%dgval(itp,mqmqj,ipy)=phres%dgval(itp,mqmqj,ipy)+&
                    aff*refg(zp,itp)
            enddo
! Initially ignore 2nd derivatives, d2G/dy2=1/y set by entropy calculation
!            write(*,213)zp,mqmqj,ipy,phres%gval(1,ipy),pyq,aff,aff*refg(zp,1)
213         format('3X adding reference state to SNN parameter:',3i3,4(1pe12.4))
            do itp=1,6
               phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*aff*refg(zp,itp)
            enddo
!            write(*,205)'3X SNN: qix, FNN, aff, pyq, fun, DG: ',mqmqj,zp,aff,&
!                 pyq,refg(zp,1),pyq*aff*refg(zp,1)
!                 (phres%dgval(1,s2,1),s2=1,gz%nofc)
         enddo snnloop
      endif pair
   enddo qloop
!   goto 800
!---------------------------------------------------------------------
!--------------------- code below needed for excess parameters
! NOTE many of them may be missing
   mqmqj=0
   endmemrec=>phlista(lokph)%ordered
   endmemloop2: do while(associated(endmemrec))
      mqmqj=mqmqj+1
      kend=mqmqa_data%contyp(5,mqmqj)
!      write(*,311)associated(endmemrec),mqmqj,mqmqa_data%nconst,kend
311   format('3X second loop for endmember records: ',l2,3i5)
      if(mqmqj.gt.mqmqa_data%nconst) exit endmemloop2
      intrec=>endmemrec%intpointer
! interaction parameters are NOT linked from SNN endmembers!!
!      write(*,*)'3X Check for interaction parameters 1',associated(endmemrec),&
!           associated(intrec),mqmqj,kend
      if(.not.associated(intrec)) cycle endmemloop2
! this is an endmember parameter with possible excess parameters
      write(*,'(a,2i3)')'3X endmember with excess parameter:',mqmqj
! just excess parameters, we must calculate product of fractions
      dpyq=zero
! the fraction variable is found from the endmember record
      id=endmemrec%fraclinks(1,1)
      pyq=phres%yfr(id)
      if(pyq.lt.bmpymin) pyq=bmpymin
      if(pyq.gt.one) pyq=one
! the endmember parameter already calculated
      dpyq(id)=one
! BRANCH for intrec%highlink and intrec%nexlink
!      write(*,'(a,i2,F10.6,6(1pe12.4))')'3X SNN df/dy: ',id,pyq,&
!           (dpyq(itp),itp=1,gz%nofc)
      proprec=>intrec%propointer
      typty=proprec%proptype
      if(typty.ne.1) stop 'illegal typty in mqmqa model'
      ipy=1
      lokfun=proprec%degreelink(0)
      call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
      if(gx%bmperr.ne.0) goto 1000
      write(*,'(a,i4,6(1Pe12.4))')'3X excess vals1:',lokfun,vals
      if(ipy.eq.1) then
         vals=vals/rtg
      endif
! ---------------------------------
! unfinished below
! ---------------------------------
      cycle endmemloop2
! multiply with fractions and derivatives of fractions                     
! gz%nofc is total number of independent fractions
         do s1=1,gz%nofc
            do itp=1,3
               phres%dgval(itp,s1,ipy)=phres%dgval(itp,s1,ipy)+&
                    dpyq(s1)*vals(itp)
            enddo
         enddo
! Initially ignore 2nd derivatives, d2G/dy2=1/y set by entropy calculation
! finally add to integral properties
         do itp=1,6
            phres%gval(itp,ipy)=phres%gval(itp,ipy)+pyq*vals(itp)
         enddo
         proprec=>proprec%nextpr
!         write(*,200)'3X SNN G, dG/dq: ',phres%gval(1,1),&
!              (phres%dgval(1,s1,1),s1=1,gz%nofc)
!      enddo mq2
! check for excess parameters ....
!      intrec=>endmemrec%intpointer
!      write(*,*)'3X Check for interaction parameters 2',associated(endmemrec),&
!           associated(intrec),mqmqj
!      do while(associated(intrec))
! calculate the excess parameter!!!
!         write(*,*)'3X there is an interaction parameter!!'
!      enddo
! next endmember ....
      endmemrec=>endmemrec%nextem
   enddo endmemloop2
!----------------------------------------------------- end SNN loop
800 continue
!   write(*,990)'3X exit calc_mqmqa G:',phres%gval(1,1),&
!        (phres%dgval(1,s1,1),s1=1,gz%nofc)
990 format(a,5(1pe14.6))
1000 return
 end subroutine calc_mqmqa

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

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine cgint
!\begin{verbatim}
 subroutine cgint(lokph,lokpty,moded,vals,dvals,d2vals,gz,tooprec,ceq)
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
   TYPE(gtp_tooprec), pointer :: toopx
! zeroing 5 iq, and vals, dvals and d2vals
!   write(*,*)'3X cgint 1:',gz%iq(1),gz%iq(2),gz%iq(3)
! why zero qz%iq, it has been set before calling ...
   gz%iq=0
!   vals=zero
!   dvals=zero
!   d2vals=zero
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
!      if(associated(tooprec)) &
!           write(*,*)'3X Toop/Kohler binary parameter constant'
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
! copy tooprec as we must not change tooprec inside calc_toop
! Toop/Kohler require recalculating the binary compostion used to
! describe the composition dependence of the parameter.  If it is 
! not composition dependent we never come here as we exit 50 lines above
         toopx=>tooprec
! copy tooprec to toopx as toopx will be changed inside calc_toop
         call calc_toop(lokph,lokpty,moded,vals,dvals,d2vals,gz,toopx,ceq)
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

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine calc_toop
! called from subroutine cgint(lokph,lokpty,moded,vals,dvals,d2vals,gz,ceq)
!\begin{verbatim}
 subroutine calc_toop(lokph,lokpty,moded,vals,dvals,d2vals,gz,toopx,ceq)
! Handle a binary interaction that is in a Toop or Kohler model
! toop is the link to the kohler-Toop record
! We come here only if the parameter is composition dependent using RK series
   implicit none
   integer moded,lokph
   TYPE(gtp_property), pointer :: lokpty
   TYPE(gtp_parcalc) :: gz
   double precision vals(6),dvals(3,gz%nofc)
   double precision d2vals(gz%nofc*(gz%nofc+1)/2)
   TYPE(gtp_equilibrium_data), pointer :: ceq
   TYPE(gtp_tooprec), pointer :: toopx
!\end{verbatim}
! we need to save this pointer from toopx
   TYPE(gtp_phase_varres), pointer :: phres
! fraction values to be used in RK series
   double precision x12,x21,sigma,dxrk
   integer, allocatable, dimension(:) :: dsigma, dx12, dx21
! ternary fraction index
   integer jj(3),j1,j2,j3,link,count,toopconst,limit,jdeg,lfun,nyfr
! for the RK calculation with Toop/Kohler fractions!
   double precision valtp(6)
   double precision dx,dx0,dx1,dx2,dxi,dxj,fff,rtg
! The first part here is to modify the fractions to be used in the RK series
! the gz record has information which elements involved
! gz%iq(1) and gz%iq(2) are index of the binary constituents
! We must also handle first and second derivatives wrt all fractions.
! we come here from a binary interaction recird but we may have to follow
! links to several other toopx records with other third elements.
! Use the phres passed on via toopx%phres if there are more toopx records
   phres=>toopx%phres
!
   write(*,10)gz%iq(1),gz%iq(2),lokpty%degree
10 format('3X in calc_toop with binary; degree:',2i3,'; ',i2,2x,20('<'))
! note vals, dvals and d2vals hare zero here
   rtg=gz%rgast
   if(lokpty%degree.eq.0) then
! quick exit if no composition dependence
      lfun=lokpty%degreelink(0)
      call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
      if(gx%bmperr.ne.0) goto 1000
      if(lokpty%proptype.eq.1) then
         valtp=valtp/rtg
      endif
! this is multiplied with y_i y_j (and derivatives) at the return
      vals=vals+valtp
      goto 1000
   endif
! We have to calculate the reduced fractions
   nyfr=size(phres%yfr)
   allocate(dsigma(nyfr))
   allocate(dx12(nyfr))
   allocate(dx21(nyfr))
! max number of binaries ...
   limit=nyfr*(nyfr-1)
! initially dsigma is set to 1, i.e. derivatives with respect to all
! constituents.  Those subtracted will be set to zero
   dsigma=one; dx12=0; dx21=0
   sigma=one
   x12=phres%yfr(gz%iq(1))
   x21=phres%yfr(gz%iq(2))
   dx12(gz%iq(1))=1
   dx21(gz%iq(2))=1
! This is the RK Muggianu fraction difference
   dxrk=x12-x21
   count=0
!-----------------------------------------------------------------
! See gtp documentation, Appendix A for the algorithm used here
!-----------------------------------------------------------------
   method: do while(associated(toopx))
! we may have several Toop/Kohler ternary methods for this binary
      count=count+1
      if(count.gt.limit) then
! something wrong in the data structure, eternal loop!
         write(*,*)'3X data structure error 1 in calc_toop',count
         gx%bmperr=4399; goto 1000
      endif
! in the toopx record there are indices of the fractions needed   
!     integer toop,const1,const2,const3,extra,uniqid
! Note const1 < const2 < const3; toop is 1, 2 or 3
      jj(1)=toopx%const1
      jj(2)=toopx%const2
      jj(3)=toopx%const3
      toopconst=0
      if(toopx%toop.gt.0) toopconst=jj(toopx%toop)
      write(*,30)count,toopconst,jj,toopx%uniqid
30    format('3X Ternary method record:',i2,', T/K: ',i1,5x,3i3,', ID: ',i2)
! we have to figure out which constituent is neither gx%iq(1) or iq(2)
! and we have find the idex for thr next toopx record (if any).
! The toopx record is linked from all 3 binaries and there can be links from
! this record that must be followed.  The constituents are in increasing order 
! Could this be simplified using select case ??
      ett: if(jj(1).eq.gz%iq(1)) then
         tva: if(jj(2).eq.gz%iq(2)) then
!---------------------------------------------------------------
! we are dealing with the binary: 1-2, any link to next is next12
!            link=1
            toopx=>toopx%next12
            if(toopconst.eq.0) then
! the binary 1-2 is a Kohler extrapolation ftom 3
! There is no derivative in dsigma with respect to this constituent
               sigma=sigma-phres%yfr(jj(3))
               dsigma(jj(3))=0
            elseif(toopconst.eq.jj(1)) then
! constituent 1 is a Toop element in 1-2-3
               x21=x21+phres%yfr(jj(3))
               dx21(jj(3))=1
            elseif(toopconst.eq.jj(2)) then
! constituent 2 is a Toop element in 1-2-3
               x12=x12+phres%yfr(jj(3))
               dx12(jj(3))=1
            endif
! if toopconst.eq.jj(3) do nothing
         elseif(jj(3).eq.gz%iq(2)) then
!---------------------------------------------------------------
! we are dealing with the binary: 1-3, any link to next is next13
!            link=2
            toopx=>toopx%next13
            if(toopconst.eq.0) then
! the binary 1-3 is a Kohler extrapolation ftom 2
! There is no derivative in dsigma with respect to this constituent
               sigma=sigma-phres%yfr(jj(2))
               dsigma(jj(2))=0
            elseif(toopconst.eq.jj(1)) then
! constituent 1 is a Toop element in 1-2-3
               x21=x21+phres%yfr(jj(2))
               dx21(jj(2))=1
            elseif(toopconst.eq.jj(3)) then
! constituent 2 is a Toop element in 1-2-3
               x12=x12+phres%yfr(jj(2))
               dx12(jj(2))=1
            endif
! if toopconst.eq.jj(2) do nothing
         else
! something is wrong in the data structure
            write(*,*)'3X data structure error 2 in calc_toop'
            gx%bmperr=4399; exit method
         endif tva
      elseif(jj(2).eq.gz%iq(1)) then
!---------------------------------------------------------------
! we are dealing with the binary: 2-3, any link is next23
!         link=3
         toopx=>toopx%next23
         if(toopconst.eq.0) then
! the binary 2-3 is a Kohler extrapolation ftom 1
! There is no derivative in dsigma with respect to this constituent
            sigma=sigma-phres%yfr(jj(1))
            dsigma(jj(1))=0
         elseif(toopconst.eq.jj(2)) then
! constituent 2 is a Toop element in 1-2-3
            x21=x21+phres%yfr(jj(1))
            dx21(jj(1))=1
         elseif(toopconst.eq.jj(3)) then
! constituent 3 is a Toop element in 1-2-3
            x12=x12+phres%yfr(jj(1))
            dx12(jj(1))=1
         endif
! if toopconst.eq.jj(1) do nothing
      else
! something wrong in the data structure
         write(*,*)'3X data structure error 3 in calc_toop'
         gx%bmperr=4399; exit method
      endif ett
!-----------------------------------------------------
! if toopx associated here extract information from that ternary method
! if toopx not associated exit here
   enddo method
! when we come here when we calculated x12, x21 and sigma
   write(*,50)'3X done:   ',x12,x21,sigma,(x12-x21)/sigma,dxrk
50 format(a,9F8.4)
51 format(a,9i8)
   write(*,50)'3X yfr:    ',phres%yfr
   write(*,51)'3X dsigma: ',dsigma
   write(*,51)'3X dx12:   ',dx12
   write(*,51)'3X dx21:   ',dx21
! Now we use (x12-x21)/sigma in the RK series ...   
! derivatives of fractions in dx12, dx21 and dsigma
   dx0=(x12-x21)/sigma
   rtg=globaldata%rgas*ceq%tpval(1)
   dx=one
   dx1=zero
   dx2=zero
   RK: do jdeg=0,lokpty%degree
      lfun=lokpty%degreelink(jdeg)
      call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
      if(gx%bmperr.ne.0) goto 1000
      if(lokpty%proptype.eq.1) then
         valtp=valtp/rtg
      endif
! the integral property values
      vals=vals+dx*valtp
      write(*,71)'3X toop: ',jdeg,gz%iq(1),gz%iq(2),nyfr,vals(1),dx,valtp(1)
71    format(a,4i2,6(1pe12.4))
      dercal1:if(moded.gt.0) then
         if(jdeg.gt.0) then
            ktloop1: do j1=1,nyfr
! derivatives with respect to all constituents!!
! dxij, dxji and dsigma non-zero for the relevant constituents!
! This must be exact!!!
               fff=(dx12(j1)-dx21(j1)-dx0*dsigma(j1)/sigma)
               write(*,8)'3X fff1:',j1,dx12(j1),dx12(j1),dsigma(j1),dx0,dx1,dx2
8              format(a,4i3,3(1pe12.4))
               dvals(1,j1)=dvals(1,j1)+fff*dx1*valtp(1)
               dvals(2,j1)=dvals(2,j1)+fff*dx1*valtp(2)
               dvals(3,j1)=dvals(3,j1)+fff*dx1*valtp(3)
               dercal2: if(moded.gt.1) then
                  ktloop2: do j2=j1,nyfr
! 2nd derivatives wrt j1 and j2 using dx12(), dx21() and dsigma()
! This need only to be approximate ...
                     fff=dsigma(j1)*(-dx12(j2)+dx21(j2)+dx0*dsigma(j2))/sigma**2
                     write(*,8)'3X fff2:',j1,j2,0,0,fff
                     d2vals(ixsym(j1,j2))=d2vals(ixsym(j1,j2))+fff*dx2*valtp(1)
                  enddo ktloop2
               endif dercal2
            enddo ktloop1
         endif
      endif dercal1
      dx2=(jdeg+1)*dx1
      dx1=(jdeg+1)*dx0
      dx=dx*dx0
      write(*,'(a,i2,4(1pe12.4))')'3X dx etc:',jdeg,dx0,dx,dx1,dx2
   enddo RK
!      
1000 continue
   write(*,1001)vals,(dvals(1,j1),j1=1,nyfr),d2vals
1001 format('3X KT: ',6(1pe11.3)/7x,6(1pe11.3)/7x,6(1pe11.3))
! Normally only derivatives wrt gz%iq(1) and gz%gq(2) but with Kohler/Toop
! there can be derivatives wrt any constituent in the same sublattice
   return
 end subroutine calc_toop

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
! d2dval(ixsym(N*(N+1)/2),1) are derivatives of G wrt fractions N and M
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
!   integer, dimension(nsl) :: nkl
   TYPE(gtp_phase_varres), pointer :: phvar
!\end{verbatim}
! SRO species AB with stoichiometry AaBb FOR BINARY TEST ONLY
! -S/R = yAln(yA)+yBln(yB)+yABln(yAB/(yA^abyB^b))
!      =(yA-ayAB)ln(yA) + (yB-byAB)ln(yB) + yABln(yAB)
! some simple way to associate yAB with yA and yB needed
! at present assume yAB is first species
   integer ll,kk,kall,nk,jl,loksp,nsl
   double precision tval,ss,yfra,ylog,ypair,ratios(2)
   ll=0
   kall=0
   nsl=phlista(lokph)%noofsubl
   if(nsl.ne.1) stop 'Illegal call to SSRO, more than one set of sites'
   sublatticeloop: do while (ll.lt.nsl)
! Hm ... we need an array with the relation between the species ...
      ll=ll+1
!      nk=nkl(ll)
      nk=phlista(lokph)%nooffr(ll)
      if(nk.ne.3) stop 'Illegal call to SSRO, not 3 constituents'
      kk=0
      ss=zero
      fractionloop: do while (kk.lt.nk)
         kk=kk+1
         kall=kall+1
! cycle if a single constituent in this sublattice
         if(nk.eq.1) cycle sublatticeloop
         yfra=phvar%yfr(kall)
         if(yfra.lt.bmpymin) yfra=bmpymin
         if(yfra.gt.one) yfra=one
         ylog=log(yfra)
! if the constituent is the first species (SRO pair) additional action needed
         if(kk.eq.1) then
! the first constituent should be the pair, 
            loksp=phlista(lokph)%constitlist(kk)
            if(splista(loksp)%noofel.eq.2) then
! conversion of rations to powers
! ratios   GNUPLOT for same y; manual for same G as ideal
! A1/2B1/2 0.8:0.8             5:5  does not converge when ordering ... SUCK
! A1/3B2/3 0.6:1.0
! A1/4B1/4 0.5:1.1               
! these are for A1/2B1/2
               ratios(1)=5*splista(loksp)%stoichiometry(1)
               ratios(2)=5*splista(loksp)%stoichiometry(2)
            else
               stop 'Illegal call to SSRO, first constituent not a pair'
            endif
            ypair=yfra
         else
! second and third constituent are elements, subtract contribution from pair
            write(*,'(a,i2,F10.6)')'3X ssro 2: ',kk,-ratios(kk-1)*ypair*ylog
            ss=ss-ratios(kk-1)*ypair*ylog
         endif
! gval(1:6,1) are G and derivator wrt T and P
! dgval(1,1:N,1) are derivatives of G wrt fraction 1:N
! dgval(2,1:N,1) are derivatives of G wrt fraction 1:N and T
! dgval(3,1:N,1) are derivatives of G wrt fraction 1:N and P
! d2dval(ixsym(N*(N+1)/2),1) are derivatives of G wrt fractions N and M
! this is a symmetric matrix and index givem by ixsym(M,N)
         ss=ss+yfra*ylog
         if(moded.gt.0) then
! there are contributions from all 3 constituents
            phvar%dgval(1,kall,1)=phvar%dgval(1,kall,1)+&
                 phvar%sites(ll)*(one+ylog)
! 2nd derivative, kxsym same as ixsym when first index is >= second index
            phvar%d2gval(kxsym(kall,kall),1)=phvar%sites(ll)/yfra
         endif
         if(kk.gt.1) then
! extra terms, assume %sites(ll) is one
            phvar%dgval(1,1,1)=phvar%dgval(1,1,1)-ratios(kk-1)*ylog
            phvar%dgval(1,kk,1)=phvar%dgval(1,kk,1)-ratios(kk-1)*ypair/yfra
! assume extra terms in second derivatives not needed
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
 end subroutine config_entropy_ssro

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

!\addtotable subroutine config_entropy_tisr provided by Ed Kremer
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
      write(*,*)'3X constituents not 5!',ncon
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

!\addtotable subroutine config_entropy_mqmqa
!\begin{verbatim}
 subroutine config_entropy_mqmqa(moded,ncon,phvar,phrec,tval)
!
! calculates configurational entropy/R for the MQMQA liquid
! started 2021-01-11
! Finally understood the model 2021-10-11
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
!--------------------------------------------------------------------------
! Problem with the formula unit of phase.  The MQMQA equations are defined
! for the number of moles of the atomes in the phase.
! in OC I work with formula units (FU) of a phase "%sites(1)=1" for MQMQA
! and the number of moles of atoms/FU, "%abnorm(1)" varies with the composition.
! Aditionaly the number of formula units "%amfu" may vary by the minimizer.
! The product "%amfu*%abnorm(1)" should be equal to the MQMQA amount of phase
! To achieve this all values should be multipled by %abnorm(1) ???
! collaboratuon with Vaishnvi Tiwari
!---------------------------------------------------------------------------
! separate documentation, i,j in first subl, k,l in second subl
! p_ijkl is cluster fraction; x_i site fraction; v_ik pair fraction
! w_i coordination equivalent site fraction;
! \sum_i x_i ln(x_i) + \sum_i x_k ln(x_k)+            subattice fractions
! \sum_i\sum_k v_ik ln(v_ik/w_i w_k)+                 pair fractions
! \sum_i\sum_k p_iikk ln(p_iikk/((v^4_ik/(w^2_i w^2_k)))+           
! \sum_i\sum_j\sum_k p_ijkk ln(p_ijkk/(2(v^2_ik v^2_jk)/(w_i w_j w^2_k)))+  
! \sum_i\sum_k\sum_l p_iikl ln(p_iikl/(2(v^2_ik v^2_il)/(w^2_i w_k w_l)))+
! \sum_i\sum_j\sum_k\sum_l p_ijkl ln(
!                         p_ijkl/(4(v_ik v_il v_jk v_jl)/(w_i w_j w_k w_l)))
!---------------------------------------------------------------------------
   type(gtp_species) :: sprec
! fq is max number of quads
! fz max number of constituents on a sublattice
! f1 for other allocated arrays
   integer, parameter :: fq=50, fz=20, f1=50
! number of pairs and sublattice fractions
   integer noofpair,ncons1,ncons2
! not needed ....
!   integer loksp,nspel,ielno(10),nextra,ncation
! these are used to as index of species on sublatte 1 (ee,ff) and 2 (gg,hh)
   integer ee,ff,gg,hh
! loop variables
   integer s1,s2,s3,s4,em,c1
! pointer to mqmqaf record with all fraction records
   type(gtp_mqmqa_var), pointer :: mqf
! site fractions and amounts
!   double precision yy1(fz),yy2(fz),nn1(fz),nn2(fz)
   double precision yy1(fz),yy2(fz)
! fractions in sublattices
   double precision sum1,sum2,sum3,sum4,half
! contyp(1-4,i) specify sublattice +/- of element and if alone or mixing 2/1
! contyp(5,i) is is the pair index for a quadrupols that is a pair
! contyp(6-7,i) for a pair are species index
! contyp(8-9,i) for a pair are ZERO
! contyp(6-9,i) for other quadrupols are pair indices (2 or 4 indicies)
! contyp(10,i)  should be i ... just as a check
! contyp(11,i)  for a pair is constituent index in sublattice 1
! contyp(12,i)  for a pair is constituent index in sublattice 1
! contyp(11-12,i) for other quadrupoles are zero
   integer em1,em2,em3,em4,mpj
! %pinq(pair) is index in %contyp for a pair
! cridx(pair_index) is the index of corresponding quad %contyp(5,q) is pair
!   integer cridx(f1) REDUNDANT
! Index to the 2-4 sublattice fractions associated with a quad
!  integer fyqix(2,fq),fyqix2(2,fq)
! pair and coord.equiv fractions for pairs in a  quad
   double precision pair(fq),ceqf1(fq),ceqf2(fq)
! test correct way to calculate pair fraction
   double precision cpair(fq),dcpair(fq,fq),cpairsum,dcpairsum(fq),dp(fq,fq)
   double precision spair(fq)
! various factors
   double precision sm1,term,fffy,fff1,fff2
   double precision ffem,ffceq1,ffceq2,once1,once2
! indicate which species that are involved in a quadrupole
!   integer involved(noofsp),stoix1,stoix2
! species in sublatice 1 and 2
   integer nspin(2),eesub,ggsub
! first and second derivatives wrt constituents ...
!   double precision dma1 is coefficent of site fraction in subl 1 for quad 
!   double precision dms1 is sum of coefficents in subl 1 for a quad i
! sum each part separately
   double precision tsub,dvvv(fq,fq),lsub(fq),tend
   double precision ssub,dssub(fq),send,dsend(fq),squad,dsquad(fq)
   double precision d2ssub(fq*(fq+1)/2)
! first index is sublattice constituent, second is quad index
   double precision b1iA(fz,fq),b2iX(fz,fq),b1iAB(fq),b2iXY(fq),sum1AB,sum2XY
! this should give stoichiometry of (species,quad) on the two sublattices
   double precision dmy1(fz,fq),dmy2(fz,fq)
! second derivatives d2xx of site fractions ...
   double precision d2yy1(fz,fq*(fq+1)/2),d2yy2(fz,fq*(fq+1)/2)
   double precision dpair(f1,fq),dceqf(f1,fq),yfrac,dummy1,dq1,dq2,dq3
   double precision dyy1(fz,fq),dyy2(fz,fq),dceqf1(f1,fq),dceqf2(f1,fq)
   double precision dsm1(fq),d2sm1(fq*(fq-1)/2),dterm(fq),ojoj,alone1,alone2
   character conname*24,endname*24,spname1*24,spname2*24,connames(fq)*24
   double precision endkvot(fq),dendkvot(fq,fq),d2endkvot(fq,fq*(fq+1)/2)
   double precision mulceq(f1),dmulceq(f1,fq*(fq+1)/2),divisor
! this is a scaling with total amount of atoms
   double precision invnorm,fqq,pairceq
! quad entropy rewritten ...
   integer pair1,pair2,pair3,pair4,e2,f2,g2,h2
! save here the indices of constituents in sublattices of pairs
! needed for the charge equivalent fractions, ceqf1 and ceqf2
! MAYBE NOT NEEDED when %contyp(11..12,quad) is are constituent indices?
   integer eij(2,fq),nomix,all2,q1,s7
! modfied AB/XY loop requires, pq is pair indices, subcon is sublattices indices
! pqq is pair index in %contyp ...
! fq is index in corresponding %constoi
   integer pq(4),pqq(4),sq1(2),sq2(2),fq1(2),fq2(2)
! to avoid adding quadrupols twice
   logical done,ddebug
! The mqmqa_data provides information about the constituents.
! In set_constititution abnorm is set to the real number of atoms in the phase
! Here we set this to unity and divides all values by abnorm
! Do we need to save the original value for model parameters??
   ddebug=.FALSE.
   invnorm=phvar%abnorm(1)
!   invnorm=one
!   phvar%abnorm(1)=one
!   phvar%abnorm(1)=one
! We should probably update abnorm(2) and (3) also ...
!   phvar%abnorm(2)=invnorm*phvar%abnorm(2)
!   phvar%abnorm(3)=invnorm*phvar%abnorm(3)
!
!   write(*,'(a,i3,1pe12.3)')'3X in MQMQA, version 5: ',ncon,one/invnorm
!
   if(.not.allocated(mqmqa_data%contyp)) then
      write(*,*)'3X MQMQA missing constituent information'
      gx%bmperr=4399; goto 1000
   endif
   if(ncon.ne.mqmqa_data%nconst) then
      write(*,*)'3Xncon, %nconst: ',ncon,mqmqa_data%nconst
      stop '3X some problems with mqmqa ...'
   endif
11    format(a,4(F5.2,2x))
!   write(*,*)'3X error return as unfinished'
!   gx%bmperr=4399; goto 1000
!   if(.not.allocated(phvar%mqmqaf%yy1)) then
! THIS MOVED BELOW BUT SHOULD EVENTUALLY BE HERE
! allocating fraction arrayes for use in entropy an excess calculations
!      write(*,*)'3X allocating phase_varres%mqmqaf arrays'
!      allocate(phvar%mqmqaf%yy1(20))
!      allocate(phvar%mqmqaf%yy2(20))
!... add more ...
!   endif
! to avoid typing too much (should mqmqaf be a target? compiler not complains)
! problem allocating arrays to this pointer !!!
   mqf=>phvar%mqmqaf
   do s1=1,ncon
! wow, using phase constituent order to find quad name !! Keep it at present
      conname=splista(phrec%constitlist(mqmqa_data%contyp(10,s1)))%symbol
      connames(s1)=conname
      if(ddebug) write(*,3)s1,(mqmqa_data%contyp(s2,s1),s2=1,14),&
           (mqmqa_data%constoi(s2,s1),s2=1,4),phvar%yfr(s1),trim(conname)
3     format('3X mq:',i2,1x,4i2,1x,i3,1x,4i2,1x,i2,4i3,4F5.1,F5.2,1x,a)
   enddo
   if(ddebug) then
      do s1=1,ncon
         write(*,4)s1,(mqmqa_data%pp(s2,s1),s2=1,4),trim(connames(s1))
4        format('3X pp:',i2,4(F8.5),2x,a)
      enddo
   endif
!   write(*,'(a,20i3)')'3X pinq: ',(mqmqa_data%pinq(s1),s1=1,mqmqa_data%npair)
!   write(*,6)phvar%yfr
6  format('3X y: ',9F7.4)
! maybe use mqf variables?  Need allocation
! local variables can be replaced by those stored in phvar
! local fraction variables and derivatives
   yy1=zero; yy2=zero; pair=zero; ceqf1=zero; ceqf2=zero;
   b1iA=zero; b2iX=zero; dpair=zero; dceqf1=zero; dceqf2=zero
   b1iAB=zero; b2iXY=zero; dmy1=zero; dmy2=zero
   cpair=zero; dcpair=zero
! write(*,431)'3X d2S/Rx:',((phvar%d2gval(ixsym(s2,s3),1),s3=s2,ncon),s2=1,ncon)
! any species used below is indicated by a 1 or 2 depending on sublattice
!   do s1=1,ncon
!      if(mqmqa_data%contyp(5,s1).ne.s1) then
!         write(*,*)'3X *** Warning %contyp index 10 not correct'
!      endif
!   enddo
! these count the sum of element and pair stoichiometries for a quad
!   fyp=zero
!----------------------------------------------------
! the array species in each sublattice will have missing values
!----------------------------------------------------
! we must calculate a number of auxilliary fraction variables from the
! site fractions using mqmqa_data%contyp
!   do s1=1,ncon
!      write(*,14)'3X %%contyp: ',s1,(mqmqa_data%contyp(s7,s1),s7=1,14),&
!           trim(connames(s1))
14    format(a,i3,1x,4i2,1x,i3,1x,4i3,1x,i3,1x,4i3,1x,a)
!   enddo
   mpj=mqmqa_data%npair
!   write(*,15)'3X pinq: ',mpj,(mqmqa_data%pinq(s1),s1=1,mpj)
15 format(a,i3,2x,20i4)
   nspin(1)=mqmqa_data%ncon1
   nspin(2)=mqmqa_data%ncon2
!   noofpair=mqmqa_data%npair
!   write(*,'(a,2i3,2x,i3)')'3X subl const and pairs: ',nspin,mpj
   noofpair=0
! BIG LOOP OVER ALL QUADS, calculating fracions of pairs, sublattices etc
   sumfrac: do s1=1,ncon
      conname=connames(s1)
      if(mqmqa_data%contyp(10,s1).ne.s1) then
         write(*,*)'3X mqmqa_data%contyp(10,s1) not redundant'
      endif
      s3=mqmqa_data%contyp(5,s1)
      typ: if(s3.gt.0) then
! AN PAIR quadrupol AA:XX, increment the pair counter
! the index of the quadrupole fraction is in ALSO in contyp(10,s1) 
!         yfrac=phvar%yfr(mqmqa_data%contyp(10,s1))
         yfrac=phvar%yfr(s1)
! Pair fractions has to be normallized later, here multiply with %pp
         noofpair=noofpair+1
         pair(s3)=pair(s3)+yfrac
         dpair(s3,s1)=one
         cpair(s3)=cpair(s3)+yfrac*mqmqa_data%pp(1,s1)
! dcpair( pairindex, quadindex )
         dcpair(s3,s1)=mqmqa_data%pp(1,s1)
! Calculating pairs
!         write(*,'(a,2i3,7F8.4)')'3X pair1: ',s1,s3,one,yfrac,pair(s3),&
!              mqmqa_data%pp(1,s1),cpair(s3),dcpair(s3,s1),dpair(s3,s1)
! all second derivatives of pair is zero
! the index of constituent in first sublattice is in %contyp(11,s1)
! ee is species index, eesub is index of species in sublattice
         ee=mqmqa_data%contyp(6,s1)
         eesub=mqmqa_data%contyp(11,s1)
         eij(1,noofpair)=eesub
! the index of constituent in second sublattice is in %contyp(12,s1)
         gg=mqmqa_data%contyp(7,s1)
         ggsub=mqmqa_data%contyp(12,s1)
         eij(2,noofpair)=ggsub
! ee and gg are pair indices, eesub and ggsub sublatice const. indices
!         write(*,50)'3X decode1: ',s1,ee,gg,eesub,ggsub
50       format(a,i3,5x,2i4,5x,2i4)
         spname1=splista(ee)%symbol
         spname2=splista(gg)%symbol
! remember which species that are used by marking them (only needed for pairs)
! this is the stoichiometric factors of the species in the pair
         fff1=2.0d0/mqmqa_data%constoi(1,s1)
         fff2=2.0d0/mqmqa_data%constoi(2,s1)
!         else
!            write(*,*)'3X contyp error 1: ',mqmqa_data%contyp(1,s1)
!            gx%bmperr=4399; goto 1000
!         endif
! SAVE the location in sublattice array of species eesub in quad s1
! eesub and ggsub are sublattice indices
! >>>>>>>>>>>>>>>>>>>>>>> .............. EQUATION B15 part 1
         
         yy1(eesub)=yy1(eesub)+fff1*yfrac
         b1iA(eesub,s1)=fff1
!         write(*,12)'3X yy1 add1:',s1,eesub,1,yy1(eesub),fff1,yfrac
12       format(a,3i3,5F10.6)
! there is a single contribution from this quad to the site fractions
         b1iAB(s1)=fff1
         yy2(-ggsub)=yy2(-ggsub)+fff2*yfrac
         b2iX(-ggsub,s1)=fff2
         b2iXY(s1)=fff2
!         write(*,12)'3X yy2 add1:',s1,ggsub,2,yy2(-ggsub),fff2,yfrac
! equivalent sublattice fraction for the sublattice constituents
! >>>>>>>>>>>>>>>>>>>>>>>> .............. EQUATION B19 part 1
         ceqf1(eesub)=ceqf1(eesub)+yfrac
         ceqf2(-ggsub)=ceqf2(-ggsub)+yfrac
         dceqf1(eesub,s1)=one
         dceqf2(-ggsub,s1)=one
! Calculating ceq
!         write(*,333)'3X ceqf1e:',s1,0,0,1,eesub,ceqf1(eesub),&
!              yfrac,one,1,trim(spname1),trim(connames(s1))
!         write(*,333)'3X ceqf2e:',s1,0,0,2,ggsub,ceqf2(-ggsub),&
!              yfrac,one,2,trim(spname2),trim(connames(s1))
333      format(a,1x,5i3,3F10.6,' ceq',i1,'(',a,')  ',a)
! end of pair summations
      else
!--------------------------------------------------------------
! this is a quadrupol AB:XY consisting of 2 or 4 pairs typ A:X and B:Y
! the pair indices in %contyp are indicated in contyp(6..9,s1)
! IT IS A BIT INVOLVED AND CAN (certainly) BE SIMPLIFIED ....
         ffem=0.5D0
         fffy=one
         yfrac=phvar%yfr(s1)
         if(mqmqa_data%contyp(9,s1).gt.0) then
! contyp(9,s1) nonzero for quadrupoles with 4 pairs A:X, A:Y, B:X, B:Y
! set ffem=0.25 if 4 pairs
            ffem=0.25D0
! set fffy=0.5 to avoid adding same fraction twice
            fffy=0.5D0
         endif
! these refer to constituent species, ff, gg in first; gg hh in second
!         ee=0; ff=0; gg=0; hh=0
! s2 loops over the species involved in the quadrupol, it can be 3 or 4
! in %contyp(1..4,s1) is indicated if same species twice (2) or not (1)
! in %constoi(1..4,s1) is the coordination number
! s2 loops positions 1..4 in contyp and constoi
! these are used to find correct stoichimetry index
! which constoi to use? AA:XY should have (1,2) and (1,3) for AA:XX and AA:YY
! which constoi to use? AB:XX should have (1,3) and (2,3) for AA:XX and BB:XX
! which constoi to use? AB:XY should have (1,3), (1,4), (2,3) and (2,4) for ...
! position 6, 7, 8, 9 are indices to pairs, s2 incremented at loop end
! in the pairs %contyp(11,pairindex) and %contyp(12,pairindex) and subl index
         once1=one; once2=one; alone1=2.0d0; alone2=2.0d0
         ffceq1=0.5D0; ffceq2=0.5D0
         pq=0; sq1=0; sq2=0
! pq are the pair indices, 2 or 4
! but below we use pq as indices to mqmqa_data ... we need pinq(pq(j))
         pq(1)=mqmqa_data%contyp(6,s1)
         pq(2)=mqmqa_data%contyp(7,s1)
         pqq(1)=mqmqa_data%pinq(pq(2))
         pqq(2)=mqmqa_data%pinq(pq(2))
! here we saved A and X assuming mixing in first sublattice
! we must also save the stoichiometric factors of the sublattice species
         sq1(1)=mqmqa_data%contyp(11,pqq(1))
! fq1 this is index to %constoi for this sublattice constituent
         fq1(1)=1
         sq2(1)=mqmqa_data%contyp(12,pqq(1))
         fq2(1)=3
         if((mqmqa_data%contyp(1,s1).eq.2)) then
! quadruplet AA:XY, pairs AA:XX and AA:YY
! Same constituents in first sublattice, indices in %contyp(11, %contyp(6,s1))
!                                               and %contyp(12, %contyp(7,s1))
! mixing in second sublattice, same constituent twice in first
            sq1(2)=sq1(1)
            fq1(2)=fq1(1)
! replace stoichiometric factor
            fq2(1)=2
            sq2(2)=mqmqa_data%contyp(12,pqq(2))
            fq2(2)=3
            alone2=one
            nomix=1
!            write(*,'(a,2i3,2x,2i3)')'3X mixing in 2: ',sq1,sq2
         elseif(abs(mqmqa_data%contyp(3,s1)).eq.2) then
! quadrupole AB:XX, first pair AA:XX, second BB:XX
! Same constituents in second sublattice, indices in %contyp(11, %contyp(6,s1))
!                                                and %contyp(12, %contyp(7,s1))
! mixing in first sublattice, same constituent twice in second
            sq2(2)=sq2(1)
            fq2(2)=fq2(1)
! add second sublattice constituent twice
            sq1(2)=mqmqa_data%contyp(11,pqq(2))
            fq1(2)=2
            alone1=one
            nomix=2
!            write(*,'(a,2i3,2x,2i3)')'3X mixing in 1: ',sq1,sq2
         else
! quadupole AB:XY, 4 pairs used, AA:XX; AA:YY; BB:XX BB:YY
! 4 pairs, we have to add 2 more
            pq(3)=mqmqa_data%contyp(8,s1)
            pq(4)=mqmqa_data%contyp(9,s1)
            pqq(3)=mqmqa_data%pinq(pq(3))
            pqq(4)=mqmqa_data%pinq(pq(4))
            sq1(2)=mqmqa_data%contyp(11,pqq(3))
! 
            fq1(2)=2
            fq2(2)=4
! I am not sure how the pairs are arranged 
! but testing 3 pairs the sublattice constituent must be different
            if(sq1(2).eq.sq1(1)) sq1(2)=mqmqa_data%contyp(11,pqq(2))
            sq2(2)=mqmqa_data%contyp(12,pqq(2))
            if(sq2(2).eq.sq2(1)) sq2(2)=mqmqa_data%contyp(12,pqq(2))
            alone1=one; alone2=one
            nomix=4
            write(*,*)'3X reciprocal cluster',mqmqa_data%contyp(2,s1)
         endif
!         write(*,313)'3X pq mm: ',s1,pq,pqq,sq1,sq2,fq1,fq2
313      format(a,i3,2x,4i2,2x,4i2,4x,2i2,2x,2i2,4x,2i2,2x,2i2)
! contribution from all pairs included in this quadruple, nonzero pq
         pqloop: do s2=1,4
!            write(*,'(a,2i3)')'3x pqloop: ',s2,pq(s2)
            if(pq(s2).eq.0) exit pqloop
            pair(pq(s2))=pair(pq(s2))+ffem*yfrac
            dpair(pq(s2),s1)=ffem
! EMERGENCY, how to know which %pp to use for each pair???
! modified in gtp3B to ensure that pairs are correlated with %constoi ??
! s2 is assumed to be %pp index, pq(s2) constittuent index ...
            cpair(pq(s2))=cpair(pq(s2))+yfrac*mqmqa_data%pp(s2,s1)
! dcpair( pairindex, quadindex )
            dcpair(pq(s2),s1)=mqmqa_data%pp(s2,s1)
!            write(*,'(a,3i3,2F10.7)')'3X dcpair2: ',pq(s2),s1,s2,&
!                 yfrac,dcpair(pq(s2),s1)
! Calculating pairs in SNN
!            write(*,'(a,3i3,6F10.6)')'3X pair2: ',s1,s2,pq(s2),ffem,yfrac,&
!                 pair(pq(s2)),mqmqa_data%pp(s2,s1),cpair(pq(s2)),&
!                 dcpair(pq(s2),s1)
         enddo pqloop
!         write(*,'(a,i3,2x,2i3,2x,2i3)')'3X sqi: ',s1,sq1,sq2
         s7=0
         subloop: do s2=1,2
! For the site fractions and equivalent fraction ceqfi we have to
! extract all constituent species of the quadrupol s1 using the pair s3
! divided by with the coordination factor in s2 for quadrupol s1
! the species in first sublattice of the pair
            if(sq1(s2).le.0) then
               write(*,*)'3X no constituent in first sublattice!!!',s1,s2,sq1
               stop
            else
! We have to use the correct sublattice index and coordination factor !!
! eesub should be in mqmqa_data%contyp(10+s2,s1) ??  What is fq(s2)?
               eesub=mqmqa_data%contyp(10+s2,s1)
!               write(*,'(a,3i3,F8.3)')'3X sublattice index: ',&
!                    eesub,sq1(s2),fq1(s2),mqmqa_data%constoi(s2,s1)
!               eesub=sq1(s2)
! SAVE the sublattice location of species eesub for quad s1
               fff1=fffy*alone1/mqmqa_data%constoi(fq1(s2),s1)
               yy1(eesub)=yy1(eesub)+fff1*yfrac
               b1iA(eesub,s1)=fff1
!               write(*,13)'3X yy1 add2:',s1,s2,eesub,yy1(eesub),fff1*yfrac,&
!                    fff1,yfrac,fffy,alone1,mqmqa_data%constoi(fq1(s2),s1)
13             format(a,3i3,3F10.6,5(F6.3))
! there can be more than one contribution to site fraction from this quad
! nomix=1 if single in 1
               if(nomix.ne.1) then
                  b1iAB(s1)=b1iAB(s1)+fff1
               else
                  b1iAB(s1)=fff1
               endif
               ceqf1(eesub)=ceqf1(eesub)+fffy*ffceq1*yfrac
               dceqf1(eesub,s1)=fffy*ffceq1
            endif
!---------- second sublattice
            if(sq2(s2).gt.0) then
! constituent index is negative in second sublattice!!
               write(*,*)'3X no constituent in second sublattice!!!',s1,s2,sq2
               write(*,14)'3X %contyp: ',s1,(mqmqa_data%contyp(s7,s1),s7=1,14)
               gx%bmperr=4399; goto 1000
            else
! NOW the species in second sublattice of the pair NOTE negative
               ggsub=sq2(s2)
! SAVE the sublattice location of species eesub  and ggsub for quad s1
! fq1(s2) specify stoichiometry index of const. in 1st sublattice in AB/XY
! the species indexing in %contyp(11..14) is the same as for %constoi(1..4)
! fq2(s2) specify stoichiometry index of const. in 2nd sublattice in AB/XY
               fff2=fffy*alone2/mqmqa_data%constoi(fq2(s2),s1)
               yy2(-ggsub)=yy2(-ggsub)+fff2*yfrac
!               write(*,13)'3X yy2 add2:',s1,ggsub,s2,yy2(-ggsub),fff2,yfrac,&
!                    fffy,alone2,mqmqa_data%constoi(fq2(s2),s1)
52             format(a,i3,5F8.5)
               b2iX(-ggsub,s1)=fff2
! nomix=2 if single in sublattice 2
               if(nomix.ne.2) then
                  b2iXY(s1)=b2iXY(s1)+fff2
               else
                  b2iXY(s1)=fff2
               endif
331         format('3Xq n(',a2,'): ',3i3,2i3,4F7.4,2x,a)
! equivalent site fraction, each mixing element will be counted twice
! for quadrupole with 4 pairs fffy=0.25; otherwice 0.5
! >>>>>>>>>>>>>>>>>>> ................ EQUATION B17 part 2
               ceqf2(-ggsub)=ceqf2(-ggsub)+fffy*ffceq2*yfrac
               dceqf2(-ggsub,s1)=fffy*ffceq2
!               write(*,333)'3X ceqf1q:',s1,0,s3,1,eesub,ceqf1(eesub),&
!                    yfrac,fffy*ffceq1,1,trim(spname1),trim(connames(s1))
!               write(*,333)'3X ceqf2q:',s1,0,s3,2,ggsub,ceqf2(-ggsub),&
!                    yfrac,fffy*ffceq2,2,trim(spname2),trim(connames(s1))
! increment s2 for next pair in quadrupole s1
               endif
            enddo subloop
         endif typ
! problem with pair fractions ...
!         do s3=1,mpj
!            write(*,'(a,2i3,5F10.7)')'3X loop:',s1,s3,(dpair(s3,s2),s2=1,ncon)
!            write(*,'(a,2i3,5F10.7)')'3X loop:',s1,s3,(dcpair(s3,s2),s2=1,ncon)
!         enddo
      enddo sumfrac
!------------------------------ end BIG LOOP over all quads
!      do s3=1,nspin(1)
!         write(*,342)'3X b1iA(m,n):',s3,s1,(b1iA(s3,s4),s4=1,ncon)
!      enddo
!      write(*,341)'3X b1iAB(n)    :',s1,(b1iAB(s4),s4=1,ncon)
!      do s3=1,nspin(2)
!         write(*,342)'3X b2iX(m,n):',s3,s1,(b2iX(s3,s4),s4=1,ncon)
!      enddo
!      write(*,341)'3X b2iXY(n)    :',s1,(b2iXY(s4),s4=1,ncon)
!      write(*,340)'3X yy1: ',(yy1(s4),s4=1,3)
!      write(*,340)'3X yy2: ',(yy2(s4),s4=1,3)
340   format(a,7F10.7)
342   format(a,2i2,7F7.4)
720   format(a,i3,4(4I3,2x))

! debug listings:
!   write(*,*)'3X summed all amounts, next normallize'
!   write(*,720)'3X contyp:  ',0,((mqmqa_data%contyp(s2,s1),s2=11,14),s1=1,ncon)
!      write(*,200)'3X p_AB/XY:',(phvar%yfr(s1),s1=1,ncon)
!      write(*,200)'3X n1     :',(yy1(s1),s1=1,nspin(1))
!      write(*,200)'3X n2     :',(yy2(s1),s1=1,nspin(2))
!      write(*,200)'3X pairs  :',(pair(s1),s1=1,noofpair)
!      write(*,200)'3X cpairs :',(cpair(s1),s1=1,noofpair)
!      do s1=1,noofpair
!         write(*,200)'3X dcpairs:',(dcpair(s1,s2),s2=1,ncon)
!      enddo
!      write(*,200)'3X ceqf1  :',(ceqf1(s1),s1=1,nspin(1))
!      write(*,200)'3X ceqf2  :',(ceqf2(s1),s1=1,nspin(2))
!   stop
!   do s3=1,nspin(1)
!      write(*,342)'3X b1iA(m,n):',s3,s1,(b1iA(s3,s4),s4=1,ncon)
!   enddo
!   write(*,341)'3X b1iAB(n)    :',s1,(b1iAB(s4),s4=1,ncon)
!   do s3=1,nspin(2)
!      write(*,342)'3X b2iX(m,n):',s3,s1,(b2iX(s3,s4),s4=1,ncon)
!   enddo
!   write(*,341)'3X b2iXY(n)    :',s1,(b2iXY(s4),s4=1,ncon)
341 format(a,i2,7F10.7)
!-------------- we have extracted all comp.variables and their deriv wrt quads
! Now sum amounts and normallize
!
! NOTE in b1iA and b1iA the first index is subl.const, second is quad 
!    sometimes I mix them up ...
!
!   write(*,*)'Sublattice fractions and detivatives:
! first sublattice
   sum1AB=zero
!   write(*,*)'3X nspin: ',nspin
   do s1=1,nspin(1)
      sum1AB=sum1AB+yy1(s1)
!      write(*,88)'3X subl: ',s1,yy1(s1),(b1iA(s1,s2),s2=1,ncon)
   enddo
88 format(a,i2,F7.3,2x,9(F8.4))
!   write(*,'(a,F7.3,a)')'3X sum1AB: ',sum1AB
   do s1=1,nspin(1)
      yy1(s1)=yy1(s1)/sum1AB
   enddo
! second sublattice
   sum2XY=zero
   do s1=1,nspin(2)
      sum2XY=sum2XY+yy2(s1)
!      write(*,88)'3X sub2: ',s1,yy2(s1),(b2iX(s1,s2),s2=1,ncon)
   enddo
   do s1=1,nspin(2)
      yy2(s1)=yy2(s1)/sum2XY
   enddo
!   write(*,*)'3X nspin2: ',nspin
! derivatives of sublattice fractions wrt quads
   all2=ncon*(ncon+1)/2
   d2yy1=zero
   dummy1=one/sum1AB**2
   yder1: do s1=1,nspin(1)
      do s2=1,ncon
! b1iAB may contain contributions from two constituents in same quad
         dyy1(s1,s2)=(b1iA(s1,s2)-yy1(s1)*b1iAB(s2))/sum1AB
!         cycle yder1
! this gives phase matrix singuler
         do s3=1,ncon
            d2yy1(s1,ixsym(s2,s3))=&
                 (-b1iA(s1,s2)*b1iAB(s3)-b1iA(s1,s3)*b1iAB(s2)+&
                 2.0D0*yy1(s1)*b1iAB(s2)*b1iAB(s3))*dummy1
!            write(*,19)'3X dyy: ',s1,s2,s3,b1iA(s1,s2),b1iAB(s3),&
!                 b1iA(s1,s3),b1iAB(s2),2.0D0*yy1(s1),d2yy1(s1,ixsym(s2,s3))
19          format(a,3i2,6(1pe10.2))
         enddo
! try ... gives also phase matrix singular ...
!         d2yy1(s1,s1)=one/yy1(s1)
      enddo
   enddo yder1
!   do s1=1,nspin(1)
!      write(*,87)'3X d2yyj: ',1,s1,(d2yy1(s1,s2),s2=1,all2)
!   enddo
87 format(a,2i3,6(1pe10.2))
   d2yy2=zero
   dummy1=one/sum1AB**2
   yder2: do s1=1,nspin(2)
! the line below works when there are no SRO quads (species)
!      dyy2(s1,s1)=one; cycle yder2
! below needed when yy2 calculated from quads
      do s2=1,ncon
! b2iXY may contain contributions from two constituents in same quads
         dyy2(s1,s2)=(b2iX(s1,s2)-yy2(s1)*b2iXY(s2))/sum2XY
         cycle yder2
         do s3=1,ncon
            if(nspin(2).eq.1) then
! single sublattice fractions should not have any second derivaties ??
               d2yy2(s1,ixsym(s2,s3))=zero
            else
! appoximate ...
               d2yy2(s1,ixsym(s2,s3))=&
                    (-b2iX(s1,s2)*b2iXY(s3)-b2iX(s1,s3)*b2iXY(s2)+&
                    2.0D0*yy2(s1)*b2iXY(s2)*b2iXY(s3))*dummy1
            endif
         enddo
      enddo
   enddo yder2
!   do s1=1,nspin(2)
!      write(*,87)'3X d2yyj: ',2,s1,(d2yy2(s1,s2),s2=1,all2)
!   enddo
! ------------------------------------------
! calculate sublattice sites related to formula units
!   dummy1=invnorm/(sum1AB+sum2XY)
!   sum1AB=sum1AB*dummy1
!   sum2XY=sum2XY*dummy1
!   sum1AB=invnorm*sum1AB
!   sum2XY=invnorm*sum2XY
! We have to sum and normalize cpair
   cpairsum=zero
   dcpairsum=zero
   dp=zero
   do s1=1,noofpair
      spair(s1)=cpair(s1)
      cpairsum=cpairsum+cpair(s1)
      do s2=1,ncon
         dcpairsum(s2)=dcpairsum(s2)+dcpair(s1,s2)
         dp(s1,s2)=dcpair(s1,s2)
      enddo
   enddo
!   write(*,'(a,F10.6,2x,10(F8.4))')'3X cpsum:',cpairsum,&
!        (dcpairsum(s2),s2=1,ncon)
   do s1=1,noofpair
      cpair(s1)=cpair(s1)/cpairsum
      do s2=1,ncon
         dcpair(s1,s2)=(cpairsum*dp(s1,s2)-spair(s1)*dcpairsum(s2))/cpairsum**2
      enddo
! replacing pair here creates problems .... do it later
!      pair(s1)=cpair(s1)
! Calculate derivatives of pairs wrt quads, NEEDED FOR REFERENCE STATE
   enddo
!   do s1=1,noofpair
!      write(*,119)'3X cpair: ',s1,cpair(s1),(dcpair(s1,s2),s2=1,ncon)
!   enddo
119 format(a,i2,F10.7,2x,8F10.6)
!
!   check pairs are unity ... this pair fraction is wrong anyway ...
!   write(*,*)'3X pair fractions and derivatives:'
   dummy1=zero
! loop over all pairs
   do s1=1,noofpair
! Check sum is unity
      dummy1=dummy1+pair(s1)
!      write(*,120)s1,pair(s1),(dpair(s1,s2),s2=1,ncon)
   enddo
120 format('3X pairs:',i3,F7.4,1x,10F6.3)
   if(abs(dummy1-one).gt.1.0D-12) then
      write(*,*)'3X pair fractions does not add up to unity',dummy1
      write(*,'(a,10F7.4)')'3X pf: ',(pair(s1),s1=1,noofpair)
      gx%bmperr=4399; goto 1000
   endif
!
! NOW list the Charge Equivalent Fractions, related to sublattices
!   write(*,*)'3X Charge Equivalent fractions and derivatives:'
   dummy1=zero
   do s1=1,nspin(1)
! Check sum is unity
      dummy1=dummy1+ceqf1(s1)
!      write(*,81)'3X ceqf:',1,ceqf1(s1),(dceqf1(s1,s2),s2=1,ncon)
   enddo
   if(abs(dummy1-one).gt.1.0D-12) then
      write(*,*)'3X Sum of charge equivalent fractions not 1 on subl 1',dummy1
      gx%bmperr=4399; goto 1000
   endif
   dummy1=zero
   do s1=1,nspin(2)
! Check sum is unity
      dummy1=dummy1+ceqf2(s1)
!      write(*,81)'3X ceqf:',2,ceqf2(s1),(dceqf2(s1,s2),s2=1,ncon)
   enddo
   if(abs(dummy1-one).gt.1.0D-12) then
      write(*,*)'3X Sum of charge equivalent fractions not 1 on subl 2',dummy1
      gx%bmperr=4399; goto 1000
   endif
81 format(a,i2,F7.4,1x,(10F7.4))
!   write(*,*)'3X all normallized fractions calculated'
!   write(*,*)'3X error return as unfinished'
!   gx%bmperr=4399
!   goto 1000
!---------------------------------------------------------------------------
! 2021.08.24 derivatives of site fractions wrt quadrupoles??
!---------------------------------------------------------------------------
!   write(*,*)'3X quitting as not finished below'
!   gx%bmperr=4399
!   goto 1000
! fraction listings
!   write(*,200)'3X p_AB/XY:',(phvar%yfr(s1),s1=1,ncon)
!   write(*,200)'3X sites/FU  :',sum1AB,sum2XY
!   write(*,200)'3X y1     :',(yy1(s1),s1=1,nspin(1))
!   write(*,200)'3X y2     :',(yy2(s1),s1=1,nspin(2))
!   do s1=1,nspin(1)
!      write(*,202)'3X dy1/dpi:',s1,(dyy1(s1,s2),s2=1,ncon)
!   enddo
!   do s1=1,nspin(2)
!      write(*,202)'3X dy2/dpi:',s1,(dyy2(s1,s2),s2=1,ncon)
!   enddo
! same as above
!   write(*,200)'3X x_A/B  :',(pair(s1),s1=1,noofpair)
!   write(*,200)'3X ceqf1  :',(ceqf1(s1),s1=1,nspin(1))
!   write(*,200)'3X ceqf2  :',(ceqf2(s1),s1=1,nspin(2))
200 format(a,(10F7.4))
202 format(a,i2,(10F7.4))
!   write(*,*)'3X now the entropy: >>>>>>>>>>>>>'
!--------------------------------------------------------------------------
! COPY ALL FRACTIONS VARIABLES AND DERIVATIVES TO MQMQAF for use in parameters
! allocate all arrays
   if(.not.allocated(phvar%mqmqaf%yy1)) then
! allocate first time only!!
      mqf%nquad=ncon; mqf%npair=noofpair; mqf%ns1=nspin(1); mqf%ns2=nspin(2)
      write(*,*)'3X allocating mqf arrays'
      allocate(mqf%yy1(nspin(1)))
      allocate(mqf%yy2(nspin(2)))
      allocate(mqf%dyy1(nspin(1),ncon))
      allocate(mqf%dyy2(nspin(1),ncon))
      allocate(mqf%d2yy1(nspin(1),ncon*(ncon+1)/2))
      allocate(mqf%d2yy2(nspin(1),ncon*(ncon+1)/2))
      allocate(mqf%ceqf1(nspin(1)))
      allocate(mqf%ceqf2(nspin(1)))
      allocate(mqf%dceqf1(nspin(1),ncon))
      allocate(mqf%dceqf2(nspin(2),ncon))
      allocate(mqf%pair(noofpair))
      allocate(mqf%dpair(noofpair,ncon))
!   else
!      write(*,*)'3X copying data to mqf arrays'
   endif
!
   do s1=1,nspin(1)
      mqf%yy1(s1)=yy1(s1)
      mqf%ceqf1(s1)=ceqf1(s1)
      do s2=1,ncon
         mqf%dyy1(s1,s2)=dyy1(s1,s2)
         mqf%dceqf1(s1,s2)=dceqf1(s1,s2)
         do s3=s2,ncon
            mqf%d2yy1(s1,ixsym(s3,s2))=d2yy1(s1,ixsym(s2,s2))
         enddo
      enddo
   enddo
   do s1=1,nspin(2)
      mqf%yy2(s1)=yy2(s1)
      mqf%ceqf2(s1)=ceqf2(s1)
      do s2=1,ncon
         mqf%dyy2(s1,s2)=dyy2(s1,s2)
         mqf%dceqf2(s1,s2)=dceqf2(s1,s2)
         do s3=s2,ncon
            mqf%d2yy2(s1,ixsym(s3,s2))=d2yy2(s1,ixsym(s2,s2))
         enddo
      enddo
   enddo
!   do s1=1,noofpair
! this will later be replaced by cpair!! for entropy the old pair works better
!      mqf%pair(s1)=pair(s1)
!      do s2=1,ncon
!         mqf%dpair(s1,s2)=dcpair(s1,s2)
! try using dpair ....
!         mqf%dpair=dpair(s1,s2)
!      enddo
!   enddo
!   write(*,777)'3X mqf sub1 1 copied:',(mqf%yy1(s1),s1=1,nspin(1))
!   write(*,777)'3X mqf sub1 2 copied:',(mqf%yy2(s1),s1=1,nspin(2))
!   write(*,777)'3X mqf pair copied:',(mqf%pair(s1),s1=1,noofpair)
!   do s1=1,noofpair
!      write(*,777)'3X mqf dpair:',mqf%pair(s1),(mqf%dpair(s1,s2),s2=1,ncon)
!   enddo
777 format(a,F10.7,2x,5(F10.6),(/5x,6F10.6))
!---------------------------------------------------------------------------
! ENTROPY CALCULATION
!---------------------------------------------------------------------------
! separate documentation, i,j in first subl, k,l in second subl
! p_ijkl is cluster fraction; x_i site fraction; v_ik pair fraction
! w_i coordination equivalent site fraction;
! \sum_i y'_i ln(y'_i) + \sum_j y"_j ln(y"_j)+        subattice fractions
!
! \sum_i\sum_k v_ik ln(v_ik/(w_i w_k))+                 pair fractions
!
! \sum_i\sum_k p_iikk ln(p_iikk/((v^4_ik/(w^2_i w^2_k)))+           
! \sum_i\sum_j\sum_k p_ijkk ln(p_ijkk/(2(v^2_ik v^2_jk)/(w_i w_j w^2_k)))+  
! \sum_i\sum_k\sum_l p_iikl ln(p_iikl/(2(v^2_ik v^2_il)/(w^2_i w_k w_l)))+
! \sum_i\sum_j\sum_k\sum_l p_ijkl ln(
!                         p_ijkl/(4(v_ik v_il v_jk v_jl)/(w_i w_j w_k w_l)))
!---------------------------------------------------------------------------
! Discovered 21/10/20 with help by Mac Poschmann:
! The entropy is distributed on the quads, dS/dquad is the sum of
! the entropy contribution from sublattices, pairs and the quads
! is related to each separate quad!  Use the dyy1(*,quadindex) etc
!-----------------------------------------------------------------
! Here we calculte for one formula unit (FU) of the phase
! at the end we multiply with current number of atomes/FU
!-----------------------------------------------------------------
   ssub=zero; dssub=zero
   dvvv=zero
!   write(*,'(a,6(1pe12.4))')'3X quads: ',(phvar%yfr(q1),q1=1,ncon)
! NEW CODE, loop over all quads
   qsub: do q1=1,ncon
! Entropy from sublattices
      tsub=zero
! replace dsub with dvvv
      s7=0
      quady: do s1=1,4
! Entropy contribution from sublattice constituents for the quad
         s7=s7+1
         s2=mqmqa_data%contyp(10+s1,q1)
         fqq=one
         if(s2.gt.0) then
! Specie in first sublattice >0, if a single species fqq=2
            if(mqmqa_data%contyp(1,q1).eq.2) fqq=2.0d0
            tsub=tsub+fqq*log(yy1(s2))/mqmqa_data%constoi(s7,q1)
!            write(*,700)'3X ssub1: ',q1,s1,s2,s7,tsub,&
!                 fqq*log(yy1(s2))/mqmqa_data%constoi(s7,q1),fqq,yy1(s2),&
!                 mqmqa_data%constoi(s7,q1)
700         format(a,4i3,2(1pe12.4),4(0PF10.6))
! the derivative of fqq*log(yy1(s2))/mqmqa_data%constoi wrt all quads!
            do s3=1,ncon
               dvvv(s3,q1)=dvvv(s3,q1)+&
                    fqq*dyy1(s2,s3)/(yy1(s2)*mqmqa_data%constoi(s7,q1))
!               write(*,706)'3X dvvv1: ',q1,s2,s3,dvvv(s3,q1)
706            format(a,3i3,4(1pe12.4))
            enddo
         elseif(s2.lt.0) then
! if a single species in second sublattice fqq=2
            if(mqmqa_data%contyp(s1,q1).eq.2) fqq=2.0d0
            tsub=tsub+log(yy2(-s2))/mqmqa_data%constoi(s7,q1)
! the derivative of fqq*log(yy2(s2))/mqmqa_data%constoi wrt all quads!
            do s3=1,ncon
               dvvv(s3,q1)=dvvv(s3,q1)+&
                    fqq*dyy2(-s2,s3)/(yy2(-s2)*mqmqa_data%constoi(s7,q1))
!               write(*,706)'3X dvvv2: ',q1,s2,s3,dvvv(s3,q1)
            enddo
         else
! no more sublattice constituents
            exit quady
         endif
! exit if this is a pair
         if(mqmqa_data%contyp(5,q1).gt.0) exit quady
      enddo quady
! first derivatives, dSsub/dquad
      lsub(q1)=tsub
      ssub=ssub+phvar%yfr(q1)*tsub
!      write(*,702)'3X ssub2: ',q1,ssub,phvar%yfr(q1),tsub
702   format(a,i3,5(1pe12.4))
   enddo qsub
! correct first derivatives with respect to quads using dvvv
!   do q1=1,ncon
!      write(*,701)'3X dvvv: ',(dvvv(s1,q1),s1=1,ncon)
!   enddo
   do q1=1,ncon
      dssub(q1)=lsub(q1)
! add on all derivatives wrt q1 from other entropy terms 
      do s1=1,ncon
         dssub(q1)=dssub(q1)+phvar%yfr(s1)*dvvv(s1,q1)
      enddo
   enddo
!   write(*,701)'3X dssub: ',(dssub(q1),q1=1,ncon)
!   write(*,701)'3X SSUB:',ssub,ssub*phvar%amfu,phvar%amfu,phvar%abnorm(1),&
!        phvar%amfu*phvar%abnorm(1)
701 format(a,5(1pe12.4))
!   stop
600 format(a,1pe12.4,2x,6(1pe10.2))
!===============
! skip the pair and quad contributions
!   write(*,*)'3X Done sublattice entropy, skipping rest',squad
!   goto 900
!
!-------------------------------------------------------
! pair entropy
   send=zero; dsend=zero
   quadcef: do q1=1,ncon
      tend=zero
! loop of all pairs of this quad
      s1=5
!      allpairs: do while(.TRUE. .and. s1.lt.10)
      allpairs: do while(.TRUE. .and. s1.lt.9)
! mqmqa_data%contyp(5,q1) is nonzero if the quad is a pair
         s2=mqmqa_data%contyp(s1,q1)
         if(s1.eq.5 .and. s2.ne.0) then
! the quad q1 is a pair with index s2, only one calculation with s2=q1
            fqq=4.0D0
            s1=10
         else
            s1=s1+1
            s2=mqmqa_data%contyp(s1,q1)
! s2 is now the index a pair in this SNN quad is in %contyp(6..9,q1) 
! exit here ifthere is no pair
            if(s2.eq.0) exit allpairs
! fqq depends on q1
            fqq=1.0D0
            if(mqmqa_data%contyp(1,q1).eq.2) then
               fqq=2.0D0
            elseif(mqmqa_data%contyp(3,q1).eq.-2) then
               fqq=2.0D0
            endif
         endif
! Here s2 is a pair of the quadrupole q1.  The pair fraction is pair(s2)
! which should be divided by ceqf1(1,s2)*ceqf2(2,s2)
! The logarithm should be multiplied by qfnnsnn for the pair.  no more??
! Entropy: quadfrac*\sum_s2 fqq*ln( pair(s2)/v_s2k/(w_i w_k))/%qfnnsnn(s2)
! MAYBE save values of "pair/(ceqf1*ceqf2)" and derivaties for later use??
! REMEMBER ceqf1 is equivalent sublattice fraction ... what is eij(1,s2)??
! eij(1..2,s2) are species in first and second sublattice of the pair
! BUT they are now in %contype(11,s2) and %contyp(12,s2) ???
! KEEP eij as it is used as link from pair to sublattice constituents
!         write(*,'(a,i3,2x,2i3,2x,2i3)')'3X keep eij?: ',s2,eij(1,s2),&
!              eij(2,s2),mqmqa_data%contyp(11,s2),mqmqa_data%contyp(12,s2)
         ee=eij(1,s2); gg=-eij(2,s2)
         dq1=ceqf1(ee)*ceqf2(gg)
         mulceq(s2)=dq1
         endkvot(s2)=pair(s2)/dq1
         fqq=fqq/mqmqa_data%qfnnsnn(s2)
! >>>>>>>>>>>>>>>>  ............. EQUATION B21 2nd line first half
! This is the entropy contribution from a pair of this quad
! %qfnnsnn is read from database
! %dfnnsnn can be different for different pairs, composition dependence???
! But it should be a sum? or is that taken care of by the sum over p_AB/XY ??
         tend=tend+fqq*log(endkvot(s2))
!         write(*,421)'3X pairs: ',q1,s1,s2,tend,endkvot(s2),&
!              fqq/mqmqa_data%qfnnsnn(s2),fqq,mqmqa_data%qfnnsnn(s2)
421      format(a,3i3,5(1pe11.3))
! first derivatives, note multiplied by p_AB/XY ....
         do s3=1,ncon
            if(s3.eq.q1) dsend(s3)=dsend(s3)+fqq*log(endkvot(s2))
            dsend(s3)=dsend(s3)+fqq/endkvot(s2)*(&
                 dpair(s2,q1)/(mulceq(s2))**2-&
                 2.0d0*pair(s2)/mulceq(s2)**4*(&
                 ceqf1(ee)*dceqf2(gg,q1)+dceqf1(ee,q1)*ceqf2(gg)))
! skip 2nd derivatives ...
         enddo
      enddo allpairs
! Finally we must multiply the tend with the quad fraction
      send=send+phvar%yfr(q1)*tend
! derivatives of send wrt quad
   enddo quadcef
!
!   write(*,600)'3X SEND: ',send,(dsend(s1),s1=1,ncon)
!========================================================================
! skip quad entropies
!   write(*,*)'3X skipping quad entropies'
!   goto 900
!========================== begin loop for all quads
!   write(*,*)'3X quadropole entropies:'
!   do s1=1,noofpair
!      write(*,440)'3X dpair/dq: ',q1,(dpair(s1,s2),s2=1,ncon)
!   enddo
440 format(a,i2,6(1pe10.2),(/20x,6e10.2))
   squad=zero; dsquad=zero
! replaced s1 by q1
   quadloop: do q1=1,ncon
      if(q1.ne.mqmqa_data%contyp(10,q1)) then
! TEST: the value in contyp(10,q1) should be q1 ...
         write(*,*)'3X problems in %contyp with quad indexing 7'
         gx%bmperr=4399; goto 1000
      endif
      lsub=zero
! New code for the general case
!                  p_i
! p_i * log( ------------------------------------)
!               xi_A/X*xi_B/X*xi_B/X*xi_B/Y
!               ---------------------------
!                  w_A * w_B * w_X * w_Y
!
      s1=mqmqa_data%contyp(5,q1)
      if(s1.gt.0) then
! this is a pair
         pair1=s1
         pair2=pair1
         pair3=pair1
         pair4=pair1
         ee=eij(1,pair1)
         ff=ee
         gg=-eij(2,pair1)
         hh=gg
! before adding this write statement hh was sometines not same as gg
! as it should be SUCK
!         write(*,*)'3X gg hh: ',gg,hh,ceqf2(gg),ceqf2(hh)
         fqq=one
!         write(*,'(a,10i3)')'3X quad1: ',q1,pair1,pair2,pair3,pair4,ee,ff,gg,hh
      elseif(mqmqa_data%contyp(9,q1).eq.0) then
! here either ee=ff or gg=hh
         pair1=mqmqa_data%contyp(6,q1)
         pair2=pair1
         ee=eij(1,pair1)
         gg=-eij(2,pair1)
         pair3=mqmqa_data%contyp(7,q1)
         pair4=pair3
         ff=eij(1,pair3)
         hh=-eij(2,pair3)
         fqq=2.0d0
!         write(*,'(a,10i3)')'3X quad2: ',q1,pair1,pair2,pair3,pair4,ee,ff,gg,hh
      else
! all ee, ff, gg, hh should be different, not certain if they are
         pair1=mqmqa_data%contyp(6,q1)
         ee=eij(1,pair1)
         gg=-eij(2,pair1)
         pair2=mqmqa_data%contyp(7,q1)
         ff=eij(1,pair2)
         hh=-eij(2,pair2)
         pair3=mqmqa_data%contyp(8,q1)
         if(ee.eq.ff) ff=eij(1,pair3)
         if(gg.eq.hh) hh=-eij(2,pair3)
         pair4=mqmqa_data%contyp(9,q1)
         fqq=4.0D0
!         write(*,'(a,10i3)')'3X quad4: ',q1,pair1,pair2,pair3,pair4,ee,ff,gg,hh
      endif
!
!      write(*,'(a,8F8.4)')'3X quadx: ',pair(pair1),ceqf1(ee),&
!           pair(pair2),ceqf1(ff),pair(pair3),ceqf2(gg),pair(pair4),ceqf2(hh)
      pairceq=fqq*pair(pair1)/ceqf1(ee)*pair(pair2)/ceqf1(ff)*&
           pair(pair3)/ceqf2(gg)*pair(pair4)/ceqf2(hh)
!      write(*,'(a,9i3,1pe12.4)')'3X quadx: ',q1,pair1,pair2,pair3,pair4,&
!           ee,ff,gg,hh,pairceq
!
      squad=squad+phvar%yfr(q1)*log(phvar%yfr(q1)/pairceq)
!      write(*,440)'3X squad: ',q1,squad,phvar%yfr(q1),pairceq
!
! New code for the general case
!                  p_i
! p_i * log( ------------------------------------)
!               xi_A/X*xi_B/X*xi_B/X*xi_B/Y
!               ---------------------------
!                  w_A * w_B * w_X * w_Y
!
! loop for derivatives
      do s1=1,ncon
         if(s1.eq.q1) lsub(s1)=log(phvar%yfr(q1)/pairceq)+one
         if(s1.eq.q1) dsquad(s1)=dsquad(s1)+log(phvar%yfr(q1)/pairceq)+one
! derivative for just q1 is OK
         lsub(s1)=lsub(s1)-phvar%yfr(q1)*&
              (dpair(pair1,s1)/pair(pair1)+dpair(pair2,s1)/pair(pair2)+&
              dpair(pair3,s1)/pair(pair3)+dpair(pair4,s1)/pair(pair4)-&
              dceqf1(ee,s1)/ceqf1(ee)-dceqf1(ff,s1)/ceqf1(ff)-&
              dceqf2(gg,s1)/ceqf2(gg)-dceqf2(hh,s1)/ceqf2(hh))
! Skipping this means I ignore effect of variable fracrion on pair and ceqf
!         dsquad(s1)=dsquad(s1)-phvar%yfr(q1)*&
!              (dpair(pair1,s1)/pair(pair1)+dpair(pair2,s1)/pair(pair2)+&
!              dpair(pair3,s1)/pair(pair3)+dpair(pair4,s1)/pair(pair4)-&
!              dceqf1(ee,s1)/ceqf1(ee)-dceqf1(ff,s1)/ceqf1(ff)-&
!              dceqf2(gg,s1)/ceqf2(gg)-dceqf2(hh,s1)/ceqf2(hh))
! skip 2nd derivatives
!         write(*,440)'3X lsub: ',s1,(lsub(s2),s2=1,ncon)
      enddo
!      write(*,440)'3X SQUAD: ',q1,squad,(dsquad(s1),s1=1,ncon)
   enddo quadloop
!
!   write(*,600)'3X SQUAD: ',squad,(dsquad(s1),s1=1,ncon)
! first derivatives are wrong ....
!   dsquad=zero
!   write(*,*)'3X quad derivatives zero'
!   goto 900
!
!***********************************************************************
900 continue
! we have multiplied with amounts above, (?) set invnorm=one
!   write(*,*)'3X second derivatives are approximate.  Atoms/FU: ',invnorm
! Values should be per formula unit!
   invnorm=one
! store results in appropriate places, values divided by RT
! This is G/RT
   phvar%gval(1,1)=phvar%gval(1,1)+invnorm*(ssub+send+squad)
! derivative of G wrt T, i.e. -S/R
   phvar%gval(2,1)=phvar%gval(2,1)+invnorm*(ssub+send+squad)/tval
   if(moded.gt.0) then
! This is if first derivatives are requested (must be exact)
      do s1=1,ncon
         phvar%dgval(1,s1,1)=phvar%dgval(1,s1,1)+&
              invnorm*(dssub(s1)+dsend(s1)+dsquad(s1))
         phvar%dgval(2,s1,1)=phvar%dgval(2,s1,1)+&
              invnorm*(dssub(s1)+dsend(s1)+dsquad(s1))/tval
         if(moded.gt.1) then
! this is if second derivatives are requested
!            do s2=s1,ncon
!               phvar%d2gval(ixsym(s1,s2),1)=phvar%d2gval(ixsym(s1,s2),1)+&
!                    invnorm*d2sm1(ixsym(s1,s2))
!            enddo
! We just set 1/quad
            dummy1=phvar%yfr(s1)
            if(dummy1.lt.1.0D-12) dummy1=1.0D-12
            phvar%d2gval(ixsym(s1,s1),1)=one/dummy1
         endif
      enddo
!      write(*,431)'3X dS/Rq  :',(phvar%dgval(1,s1,1),s1=1,ncon)
!      write(*,431)'3X d2S/Rq2:',(phvar%d2gval(s1,1),s1=1,all2)
431   format(a,6(1pe12.4),(/6x,6e12.4))
   endif
!   write(*,'(a,3(1pe14.6))')'3X MQMQA:',phvar%gval(1,1),&
!        phvar%gval(1,1)*8.31451,phvar%gval(1,1)*8.31451*phvar%amfu
! replace pair by cpair to handle endmembers
! Creates problems calculating the entropy in this routine ... SUCK
   do s1=1,mqf%npair
      mqf%pair(s1)=cpair(s1)
      do s2=1,ncon
         mqf%dpair(s1,s2)=dcpair(s1,s2)
! converge problems, maybe use dp?
!         mqf%dpair(s1,s2)=dp(s1,s2)
      enddo
   enddo
! TEST temporary fix
!   do s1=1,mqf%npair
!      write(*,'(a,F9.6,2x,10F10.6)')'3X cpair: ',mqf%pair(s1),&
!           (mqf%dpair(s1,s2),s2=1,mqf%nquad)
!   enddo
1000 continue
   return
! calculates configurational entropy/R for the FactSage MQMQA model
 end subroutine config_entropy_mqmqa

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
   phres%gval(2,1)=phres%gval(2,1)+(gc+gr)/ceq%tpval(1)+dgrdt
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

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!
