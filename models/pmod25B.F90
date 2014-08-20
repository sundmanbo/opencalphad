!
! included in pmod25.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>     6. Calculate things
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
   if(iph.le.0 .or. iph.gt.noofph) then
! the selected_element_reference phase with iph=0 is calculated separtely
      gx%bmperr=4050; goto 1000
   endif
   lokph=phases(iph)
   if(lokph.le.0 .or.lokph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   endif
!    write(*,*)'calcg 1: ',phlista(lokph)%name
! find fractions for this composition set
   if(ics.le.1) then
      jcs=1
   elseif(ics.le.phlista(lokph)%noofcs) then
      jcs=ics
   else
! no such composition set
!      write(*,*)'calcg 1 error 4072'
      gx%bmperr=4072; goto 1000
   endif
!   if(phlista(1)%noofcs.gt.1) then
! strange error that liquid (phase 1) has 3 composition set
!      write(*,*)'csbug: ',lokph,jcs,phlista(1)%noofcs
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
!         write(*,*)'calcg 2 error 4072'
!         gx%bmperr=4072; goto 1000
!      endif
!   enddo
   lokres=lokcs
!   write(*,*)'calcg 7: ',lokres,ceq%eqname(1:10)
! call using the local structure phase_varres
! results can be obtained through lokres
!   write(*,17)'calcg: ',lokph,lokres,ceq%phase_varres(lokres)%yfr(1)
17 format(a,2i4,1pe15.6)
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
   character bug*60,ch1*1,pdone*4
!   dimension sites(maxsubl),pushpop(maxpp)
   double precision, dimension(:), allocatable :: dpyq(:),d2pyq(:),d2vals(:)
   double precision, dimension(:,:), allocatable :: dvals(:,:)
   double precision vals(6),dummy(6),sum3(3)
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
   double precision, dimension(:), allocatable :: savey
   double precision, dimension(:,:), allocatable :: saveg
   double precision, dimension(:,:,:), allocatable :: savedg
   double precision, dimension(:,:), allocatable :: saved2g
   double precision, dimension(:,:), allocatable :: tmpd2g
! added when implicit none
   double precision rt,pyq,ymult,add1,sum
   integer nofc2,nprop,nsl,msl,iprop,lokdiseq,ll,id,id1,id2,lm,jl
   integer lokfun,itp,nz,intlat,ic,jd,noprop,jk,ic1,jpr,ipy,i1,j1
   integer i2,j2,ider,is,kk,ioff,norfc,iw,iw1,iw2,lprop
! to handle parameters with wildcard constituent
   logical wildc,nevertwice,first,mcs,chkperm
! debugging
   integer clist(4)
! calculate RT to normalize all Gibbs energies, ceq is current equilibrium
   rt=globaldata%rgas*ceq%tpval(1)
   ceq%rtn=rt
!-----------------------
   chkperm=.false.
   if(btest(phlista(lokph)%status1,PHFORD) .or. &
        btest(phlista(lokph)%status1,PHBORD)) then
      chkperm=.true.
! This is needed only once unless parameters are changed.  It numbers the
! interaction records sequentially for the permutations
      call palmtree(lokph)
      if(gx%bmperr.ne.0) goto 1000
   endif
!-----------------------------------------------------------------
50  continue
! local work arrays for products of Y and calculated parameters are allocated
   gz%nofc=phlista(lokph)%tnooffr
   nofc2=gz%nofc*(gz%nofc+1)/2
!   write(*,17)'calcg, ',lokph,gz%nofc,cps%yfr(1)
17 format(a,2i4,1pe15.6)
! for disordered fraction sets gz%nofc must be from disordered fraction record
! maybe these should not be allocated for moded=0 and 1
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
!   write(*,*)'calcg_i: ',gz%tpv
   gz%rgast=ceq%tpval(1)*globaldata%rgas
!   gz%rgast=ceq%tpval(1)*ceq%rgas
! this is used to check the number of times an ordered phase is calculated
   first=.true.
!-------------------------------------------------------------------
! calculate configurational entropy.
! NOTE: This is done for ordered original (site) fraction set only
   nsl=phlista(lokph)%noofsubl
   call config_entropy(moded,nsl,phlista(lokph)%nooffr,phres,gz%tpv(1))
   if(gx%bmperr.ne.0) goto 1000
!-------------------------------------------------------------------
! start BIG LOOP for all fraction variables and parameters
! there may be several different properties in addition to G like TC, MQ& etc
! each of these are stored in separate gval(*,ipy) where ipy is an integer
! set for each property. lprop is incremented by one for each new property
! found (each phase may have different) and in listprop the original type
! of property is stored.  listprop will always be associated with phmain
100 continue
   nevertwice=.true.
   lprop=2
   phmain%listprop(1)=1
   fractype=0
!   write(*,101)'calcg 100 ',nsl,phres%gval(1,1),cps%gval(1,1)
101 format(a,i4,4(1pe14.4))
!--------------------------------------------------------------------
! loop for different types of fractions: site fractions, mole fractions ...
   fractyp: do while(fractype.lt.phlista(lokph)%nooffs)
105 continue
      fractype=fractype+1
! return here for cxalculating with disordered fractions for same fraction type
110 continue
! gz%nofc is number of fraction variables, msl is number of sublattices
! for this set of fractions!!! Ordering in FCC may have 5 sublattices with
! 4 participating in ordering and one interstitial.  The second fraction
! set may have 2 sublattices, 1 for the 4 ordering and one interstitial
      fracset=phmain%disfra
      ftype: if(fractype.eq.1) then
!---------------------------------------------- ordered fraction set
         if(btest(phlista(lokph)%status1,phmfs)) then
! there is a disordered fractions set, we need fracset later
            if(fracset%totdis.ne.0) then
! the phase can totally disorder, if disordered skip ordered part
               if(btest(phmain%status2,csorder)) then
! the phase is ordered, we have to calculate this part twice
                  nevertwice=.false.
              else
! the phase is disordered, skip ordered part and just calculate disordered
                  goto 105
               endif
            endif
         endif
         gz%nofc=phlista(lokph)%tnooffr
         msl=nsl
         incffr(0)=0
         do jl=1,nsl
            incffr(jl)=incffr(jl-1)+phlista(lokph)%nooffr(jl)
         enddo
! the results will be stored in the results arrays indicated by phres
! it was set above for the ordered fraction set. 
      else
!-------------------------------------------------
! disorderd/other fraction sets, take data from  gtp_fraction_set
         msl=fracset%ndd
         gz%nofc=fracset%tnoofxfr
         incffr(0)=0
         do jl=1,msl
            incffr(jl)=incffr(jl-1)+fracset%nooffr(jl)
         enddo
! we have to deallocate and allocate local arrays, not if moded=0 or 1??
         deallocate(dpyq)
         deallocate(d2pyq)
         allocate(dpyq(gz%nofc))
         allocate(d2pyq(nofc2))
         dpyq=zero
         deallocate(dvals)
         deallocate(d2vals)
         allocate(dvals(3,gz%nofc))
         allocate(d2vals(nofc2))
! the results will be stored in result arrays indicated by phres
! for the disordered fraction set phres must be set here and the arrays zeroed
         dislink=cps%disfra
!         write(*,*)'Calc internal disordred part 1A'
         lokdiseq=dislink%varreslink
!         write(*,*)'Calc internal disordred part 1B',lokdiseq
         phres=>ceq%phase_varres(lokdiseq)
         phres%gval=zero
!         write(*,*)'Calc internal disordred part 1c'
         if(moded.gt.0) then
            phres%dgval=zero
            if(moded.gt.1) then
               phres%d2gval=zero
            endif
         endif
!         write(*,*)'Calc internal disordred part 2'
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
      endmemloop: do while(associated(endmemrec))
! The arrsy maxpmq is used for interaction permutations.  It must be
! initialized to zero at the first endmember permutation.  It is set to
! limits for the interacton permutations for all interaction records.
         maxpmq=0
         maxprec=0
         epermut=0
         sameint=0
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
               if(ll.lt.5) clist(ll)=id
!               id=emlista(lokem)%fraclinks(ll,epermut)
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
            enddo pyqloop
!----------------------------------------------------
            if(moded.eq.0) goto 150
            dpyqloop: do ll=1,msl
! here pyq is known, same loop as above to calculate dpyq(i)=pyq/y_i
               id=endmemrec%fraclinks(ll,epermut)
               if(id.gt.0) then
                  dpyq(id)=pyq/gz%yfrem(ll)
               else
                  do iw=incffr(ll-1)+1,incffr(ll)
                     dpyq(iw)=pyq
                  enddo
               endif
            enddo dpyqloop
!----------------------------------------------------
            if(moded.le.1) goto 150
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
!---- jump here if moder is 0 or 1
150         continue
!            write(*,154)'endmember permutation: ',epermut,(clist(i),i=1,4)
154         format(a,i5,4i4,'--------------------------------')
155         format(a,i5,10i4)
            proprec=>endmemrec%propointer
            emprop: do while(associated(proprec))
               typty=proprec%proptype
               if(typty.ne.1) then
! if property different from 1 (=G) find where to store it, use phmain link
                  do jl=2,lprop-1
                     if(phmain%listprop(jl).eq.typty) goto 170
                  enddo
! a new property, save its typty in listprop and increment lprop
! note that the property index typty is not used as index in gval etc
! as that can be very large. lprop is incremented by 1 for each property
! actually used in the model of the phase.
                  jl=lprop
                  phmain%listprop(jl)=typty
                  if(lprop.ge.nprop) then
                     write(*,*)'Too many parameter properties ',&
                          lprop,nprop,typty
                     gx%bmperr=7777; goto 1000
                  endif
                  lprop=lprop+1
                  phmain%listprop(1)=lprop
170               continue
                  ipy=jl
               else
                  ipy=1
               endif
! calculate function and derivatives wrt T and P
! the results from eval_tpfun must also be different in different treads ...
               lokfun=proprec%degreelink(0)
               call eval_tpfun(lokfun,ceq%tpval,vals,ceq%eq_tpres)
!               write(*,*)'calcg calling eval_tpfun 2: ',gx%bmperr,vals(1)
               if(gx%bmperr.ne.0) goto 1000
               prop1: if(ipy.eq.1) then
! property 1 i.e. Gibbs energy, should be divided by RT
                  vals=vals/rt
               endif prop1
! debug
!   write(*,173)'endmember: ',endmemrec%antalem,pyq,vals(1)
!173 format(a,i3,4(1pe12.4))
! multiply with py and derivatives. vals is composition independent
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
               proprec=>proprec%nextpr
            enddo emprop
!------------------------------------------------------------------
! take link to interaction records, use push and pop to save pyq etc
! pmq keeps track of the location in LASTPMQ and MAXPMQ
! for each interaction record in this binary interaction tree
            intrec=>endmemrec%intpointer
            gz%intlevel=0
            pmq=1
! pmq is initiated by palmtree above in the interaction records
            interloop: do while(associated(intrec))
!----------------------------------------------------------------
! come back here an interaction at a higher level or a poped next that must
! be pushed 
200            continue
               gz%intlevel=gz%intlevel+1
!---------------------------------------------
! come back here for another interaction on same level ... not used??
210            continue
               call push_pyval(pystack,intrec,pmq,&
                    pyq,dpyq,d2pyq,moded,gz%nofc)
! intrec%order is initiated by palmtree to set a sequential number
               pmq=intrec%order
!               write(*,155)'Pushed: ',pmq,gz%intlevel
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
!                           write(*,155)'new limit: ',ipermut,&
!                                maxpmq(pmq)
                           if(ipermut.le.maxpmq(pmq)) goto 230
                        elseif(gz%intlevel.gt.2) then
                           write(*,*)'Max level 2 interactions allowed'
                           gx%bmperr=7777; goto 1000
                        else
                           varying: if(intrec%noofip(1).eq.1) then
! If this is 1 then noofip(2) is number of permutations each time
                              maxpmq(pmq)=maxpmq(pmq)+intrec%noofip(2)
                              if(ipermut.le.maxpmq(pmq)) goto 230
                           else
! This is more complicated, different number of permutations each time
! Example: noofip=(3,2,1,0,12) means there are 3 different permutations
! first time 2, second time 1, last time none; 12 is total number of perms
! If mod(ipermut,noofip(1)) is 0 one should start from index 2
                              nz=intrec%noofip(1)
!                              write(*,155)'noofip: ',ipermut,pmq,maxpmq(pmq),&
!                                   (intrec%noofip(j),j=1,nz)
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
!                              write(*,155)'noperm: ',ipermut,pmq,&
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
229                           format('Error, max 2 levels of interactions',/&
                                   ' with permutations!! ',i3)
                              gx%bmperr=7777; goto 1000
                           endif
! Take the link to higher as no more permutations here
                           goto 290
!                           goto 200
                        endif
!-------------------------------
! No higher level, if we cannot pop we must return to endmember
                        if(gz%intlevel.eq.0) exit interloop
! we must pop lower order interaction records here to get correct permutation
                        call pop_pyval(pystack,intrec,pmq,&
                             pyq,dpyq,d2pyq,moded,gz%nofc)
                        gz%intlevel=gz%intlevel-1
                        pmq=intrec%order
!-----------------------------------------------------------
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
                  ipermut=1
               endif bford
! Code below until label 290 the same with and without permutations
! extract  sublattice, constituent and fraction of interacting constituent
               intlat=intrec%sublattice(ipermut)
               ic=intrec%fraclink(ipermut)
               gz%intlat(gz%intlevel)=intlat
               gz%intcon(gz%intlevel)=ic
               gz%yfrint(gz%intlevel)=phres%yfr(ic)
! calculate new PY incl derivatives. Moded to avoid unrequested derivatives
! IF interaction endmember is WILDCARD then the interaction is special,
! L(*,A) is y_A *(1-y_A) where 1-y_A is the sum of all fractions except A
! pyq = pyq * y_ic * (y_ix + y_iy + ... ) (all_other_in_same_sublattice))
! derivatives are calculated for all constituents in intlat
! note one can also have wildcards in other sublattices
               if(gz%endcon(intlat).gt.0) then
                  wildc=.FALSE.
                  ymult=gz%yfrint(gz%intlevel)
               else
                  wildc=.TRUE.
                  ymult=gz%yfrint(gz%intlevel)*(one-gz%yfrint(gz%intlevel))
               endif
               noder3: if(moded.gt.0) then
                  iloop1: do id=1,gz%nofc
                     if(moded.gt.1) then
                        iloop2: do jd=id+1,gz%nofc
                           d2pyq(ixsym(id,jd))=d2pyq(ixsym(id,jd))*ymult
                        enddo iloop2
                        d2pyq(ixsym(id,ic))=dpyq(id)
                     endif
                     dpyq(id)=dpyq(id)*ymult
                  enddo iloop1
! we must check if any endmember is wildcard like L(A,B:*)
                  do ll=1,msl
                     if(ll.ne.intlat) then
                        if(gz%endcon(ll).lt.0) then
                           do iw=incffr(ll-1)+1,incffr(ll)
                              d2pyq(ixsym(iw,ic))=pyq
                           enddo
                        endif
                     endif
                  enddo
                  if(wildc) then
! wildcard contribution to second derivate from all fractions in intlat
                     do iw=incffr(intlat-1)+1,incffr(intlat)
                        if(iw.ne.ic) then
                           d2pyq(ixsym(iw,ic))=dpyq(iw)
                        endif
!                        write(*,213)'529: ',iw,ic,ixsym(iw,ic),&
!                             gz%intlevel,intlat,incffr(intlat)
                        dpyq(iw)=pyq*gz%yfrint(gz%intlevel)
!                        dpyq(jd)=pyq*gz%yfrint(gz%intlevel)
                     enddo
213                  format(a,10i5)
                     dpyq(ic)=pyq*(one-gz%yfrint(gz%intlevel))
                  else
                     dpyq(ic)=pyq
                  endif
               endif noder3
! pyq calculated identically for wildcards as ymult set differently above
               pyq=pyq*ymult
               proprec=>intrec%propointer
! this was just for debugging
!               noprop=0
!               if(associated(proprec)) then
!                  noprop=1
!               endif
!..............................
               intprop: do while(associated(proprec))
! calculate interaction parameter, can depend on composition
                  call cgint(proprec,moded,vals,dvals,d2vals,gz,ceq)
                  if(gx%bmperr.ne.0) goto 1000
! G parameters are divided by RT inside cgint
                  typty=proprec%proptype
                  if(typty.ne.1) then
! other properties than 1 (G) must be stored in different gval(*,ipy) etc
                     do jl=2,lprop-1
                        if(phmain%listprop(jl).eq.typty) goto 250
                     enddo
! a new property, save its typty in listprop and increment lprop
                     jl=lprop
                     phmain%listprop(jl)=typty
                     lprop=lprop+1
                     phmain%listprop(1)=lprop
250                  continue
                     ipy=jl
                  else
                     ipy=1
                  endif
                  noder4: if(moded.gt.0) then
                     iloop3: do id=1,gz%nofc
                        if(moded.gt.1) then
                           iloop4: do jd=id+1,gz%nofc
                              phres%d2gval(ixsym(id,jd),ipy)= &
                                   phres%d2gval(ixsym(id,jd),ipy)+ &
                                   d2pyq(ixsym(id,jd))*vals(1)
                           enddo iloop4
                        endif
                        do itp=1,3
                           phres%dgval(itp,id,ipy)=&
                                phres%dgval(itp,id,ipy)+ &
                                dpyq(id)*vals(itp)
                        enddo
                     enddo iloop3
!...............................
! below contribution to derivatives from composition dependent parameters
! the values of gz%iq represent interacting constituents and are set in cgint
                     cdex1: if(gz%iq(5).gt.0) then
! gz%iq(5) is nonzero only for TOOP and similar models not implemented yet ...
                        gx%bmperr=4086; goto 1000
                     elseif(gz%iq(4).gt.0) then
!...............................
! composition dependent reciprocal parameter, not implemented yet ...
                        gx%bmperr=4086; goto 1000
                     elseif(gz%iq(3).gt.0) then !cedex1
!...............................
! composition dependent ternary interaction in same sublattice, Mats model
! PROBABLY ERRORS HERE as no consideration of derivatives wrt other endmember
! constituents, only to the 3 interacting
                        if(moded.gt.1) then
                        noindent1: do jk=1,3
                        do jl=jk+1,3
! the second derivative for jk=jl calculated below as it is simpler
                           phres%d2gval(ixsym(gz%iq(jk),gz%iq(jl)),ipy)=&
                                phres%d2gval(ixsym(gz%iq(jk),gz%iq(jl)),ipy)+&
                                dpyq(gz%iq(jk))*dvals(1,gz%iq(jl))+&
                                dpyq(gz%iq(jl))*dvals(1,gz%iq(jk))
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
                     elseif(gz%iq(2).gt.0) then !cedex1
!............................
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
                                 phres%d2gval(ixsym(ic1,gz%iq(2)),ipy)=&
                                 phres%d2gval(ixsym(ic1,gz%iq(2)),ipy)+add1
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
                        enddo
                     endif cdex1
! end contribution to derivates from composition dependent parameters
!......................
                  endif noder4
! finally add the contribution to G, G.T etc
                  iloop6: do itp=1,6
                     phres%gval(itp,ipy)=phres%gval(itp,ipy)+&
                          pyq*vals(itp)
                  enddo iloop6
                  proprec=>proprec%nextpr
! check for strange bug in CrFeMo
!                     write(*,144)gz%iq(1),gz%iq(2),pyq,vals(1),&
!                          (dpyq(i),i=1,3),(dvals(1,i),i=1,3)
!144 format('Interaction: ',2i3,2(1pe12.4)/3(1pe12.4),2x,3(1pe12.4))
               enddo intprop
!               if(chkperm) then
! debugg output
!                  if(noprop.eq.0) then
!                     xyz=0.0D0; pdone='no  '
!                  else
!                     xyz=vals(1); pdone='yes '
!                  endif
!                  write(*,211)pdone,xyz,gz%intlevel,intrec%antalint,&
!                       pmq,maxpmq(pmq),ipermut,intlat,ic
!211               format('int permut: ',a,1PE12.2,10i4)
!               endif
! finished one interaction (or permutation on this level), go to higher level
! note that ipermut is saved in lastpmq(pmq).  If there are more
! permutations on this level they will be calculated later also including 
! higher order parameters.  
! Take link to higher level records for current permutation
290            continue
               intrec=>intrec%highlink
               wrong: if(chkperm .and. associated(intrec)) then
! We must go to higher as we can have interactions with different permutations?
                  jpr=intrec%order
                  if(lastpmq(jpr).gt.0 .and.lastpmq(jpr).ge.maxpmq(jpr)) then
! if we nullify here we will take next rather than higher
!                     nullify(intrec)
!                     write(*,155)'Maybe skipping higer?: ',jpr,&
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
! take next permutation of the end member fractions
         enddo empermut
300      continue
! take next end member
!         write(*,155)'endmem: ',epermut,endmemrec%noofpermut,endmemrec%antalem
         endmemrec=>endmemrec%nextem
      enddo endmemloop
!------------------------------------------------------------------------
! end loop for this fraction type, initiation for next in the beginning of loop
! but we may have to calculate once again with same fraction type but
! with the fractions as disordered fractions
      if(nevertwice) then
         goto 400
      endif
      disord: if(fractype.eq.1 .and. btest(phlista(lokph)%status1,phmfs) &
           .and. btest(phmain%status2,csorder)) then
! Handle additions of several fraction set ?? Additions calculated
! after both ordered and disordered fraction set calculated
         if(first) then
! prepare to calculate again with all fractions as disordered
            first=.false.
            allocate(savey(gz%nofc))
            savey=phres%yfr
!            write(*,*)'cg: ',phmain%phlink,phmain%disfra%varreslink
! ??? very uncertain how to call disordery .....
!            call disordery(phmain,phmain%disfra%varreslink,ceq)
            call disordery(phmain,ceq)
            nprop=phmain%nprop
            allocate(saveg(6,nprop))
            allocate(savedg(3,gz%nofc,nprop))
            allocate(saved2g(nofc2,nprop))
            saveg=phres%gval
            savedg=phres%dgval
            saved2g=phres%d2gval
            phres%gval=zero
            phres%dgval=zero
            phres%d2gval=zero
            goto 110
         else
! Ordered part calculated with disordered fractions, subtract this
! from the first, restore fractions and deallocate
! THIS IS TRICKY
! NOTE all sublattices are identical in this case with the same number 
! of constituents
! First sum all second derivatives into tmpd2g
            noder6A: if(moded.gt.1) then
               nz=fracset%tnoofxfr
               allocate(tmpd2g(nz*(nz+1)/2,nprop))
               tmpd2g=zero
               do ipy=1,lprop
                  do i1=1,gz%nofc
                     j1=fracset%y2x(i1)
                     do i2=i1,gz%nofc
                        j2=fracset%y2x(i2)
                        tmpd2g(ixsym(j1,j2),ipy)=tmpd2g(ixsym(j1,j2),ipy)+&
                             phres%d2gval(ixsym(i1,i2),ipy)
                     enddo
                  enddo
               enddo
               do ipy=1,lprop
                  do i1=1,gz%nofc
                     j1=fracset%y2x(i1)
                     do i2=i1,gz%nofc
! subtract from saved value
                        j2=fracset%y2x(i2)
                        phres%d2gval(ixsym(i1,i2),ipy)=&
                             saved2g(ixsym(i1,i2),ipy)-&
                             tmpd2g(ixsym(j1,j2),ipy)*&
                             fracset%dxidyj(i1)*fracset%dxidyj(i2)
                     enddo
                  enddo
               enddo
               deallocate(tmpd2g)
            endif noder6A
!---------------------
! sum all partial derivates to first sublattice,
            noder6B: if(moded.gt.0) then
               do ipy=1,lprop
                  do ider=1,3
                     do is=1,fracset%nooffr(1)
                        sum=zero
                        kk=is
                        do ll=1,fracset%latd
                           sum=sum+phres%dgval(ider,kk,ipy)
! it is not really necessary to put phres%dgval it to zero, just for cleanness
                           phres%dgval(ider,kk,ipy)=zero
                           kk=kk+fracset%nooffr(1)
                        enddo
                        phres%dgval(ider,is,ipy)=sum
                     enddo
                     if(fracset%ndd.eq.2) then
! one can have 2 sets of ordered subl like (Al,Fe)(Al,Fe)...(C,Va)(C,Va)...
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
               do ipy=1,lprop
! loop in negative direction avoid destroy the values in phres%dgval first subl
                  do i1=gz%nofc,1,-1
! all derivatives wrt same element from all sublattices is in first sublattice
                     do ider=1,3
! Finally subtract this contribution from saved values
                        j1=fracset%y2x(i1)
                        phres%dgval(ider,i1,ipy)=savedg(ider,i1,ipy)-&
                             phres%dgval(ider,j1,ipy)*fracset%dxidyj(i1)
                     enddo
                  enddo
               enddo
            endif noder6B
            do ipy=1,lprop
               do ider=1,6
                  phres%gval(ider,ipy)=saveg(ider,ipy)-&
                       phres%gval(ider,ipy)
               enddo
            enddo
! restore ordered fractions and deallocate save arrays
            phres%yfr=savey
            deallocate(savey)
            deallocate(saveg)
            deallocate(savedg)
            deallocate(saved2g)
         endif
      endif disord
400 continue
   enddo fractyp
!--------------------------------------------------------------
! finished loops for all fractypes, now add together G and delta-G(ord-disord)
410 continue
! cheking for properties
!   write(*,411)lprop,(phmain%listprop(j),j=1,lprop)
!   write(*,412)(phmain%gval(1,j),j=1,lprop-1)
!411 format('Properties: ',i3,': ',10i4)
!412 format('Val: ',(6E12.4))
   fractionsets: if(btest(phlista(lokph)%status1,phmfs)) then
!----------------------------------------------------------------
! add together contributions from different fractypes
! phres is last calculated part, set phpart to ordered part (phmain)
      phpart=>phmain
      norfc=phlista(lokph)%tnooffr
! loop for all second and frist derivatives using chain rule
! and coefficients from fracset%dxidyj
! d2f1/dyidyj = d2f2/dxkdxl*dxk/dyi*dxl/dyj
! gz%nofc are number of disordered constituents
! norfc are number of ordered constituents
! lprop is number of properties to be summed
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
                  do ipy=1,lprop
                     phpart%d2gval(ixsym(i1,i2),ipy)=&
                          phpart%d2gval(ixsym(i1,i2),ipy)+&
                          phres%d2gval(ixsym(j1,j2),ipy)*&
                          fracset%dxidyj(i1)*fracset%dxidyj(i2)
                  enddo
               enddo
            endif noder7B
! first derivatives
            do ipy=1,lprop
               do ider=1,3
                  phpart%dgval(ider,i1,ipy)=phpart%dgval(ider,i1,ipy)+&
                    phres%dgval(ider,j1,ipy)*fracset%dxidyj(i1)
               enddo
            enddo
         enddo
      endif noder7A
! Integral values
      do ipy=1,lprop
         do ider=1,6
            phpart%gval(ider,ipy)=phres%gval(ider,ipy)+phpart%gval(ider,ipy)
         enddo
      enddo
   endif fractionsets
! now set phres to ordered results and forget phpart
   phres=>phmain
! calculate additions like magetic contributions etc and add to G
   addrec=>phlista(lokph)%additions
   additions: do while(associated(addrec))
!      if(addlista(lokadd)%type.eq.1) then
      call addition_selector(addrec,moded,phres,lokph,gz%nofc,ceq)
      if(gx%bmperr.ne.0) goto 1000
!z      if(addrec%type.eq.1) then
! moded is 0, 1 or 2 if derivatives should be calculated, phres is pointer
! to result arrays, lokadd is the addition record, listprop is needed to
! find where TC and BM are stored, gz%nofc are number of constituents
!z         call calc_magnetic_inden(moded,phres,addrec,lokph,gz%nofc,ceq)
!z         if(gx%bmperr.ne.0) goto 1000
!z      elseif(addrec%type.eq.2) then
! moded is 0, 1 or 2 if derivatives should be calculated, phres is pointer
! to result arrays, lokadd is the addition record, listprop is needed to
! find where TC and BM are stored, gz%nofc are number of constituents
!z         call calc_debye_mauro(moded,phres,addrec,lokph,gz%nofc,ceq)
!z         if(gx%bmperr.ne.0) goto 1000
!z      else
!         write(*,*)'Unknown addition type ',addlista(lokadd)%type
!z         gx%bmperr=4145; goto 1000
!z      endif
      addrec=>addrec%nextadd
   enddo additions
! there are some special properties like mobilities and similar which
! have a conmponent or constituent index like MQ&<constituent>
!   ipy=typty/100+mod(typty,100)
!   if(ipy.gt.10) then
!      write(*,*)'Property ',typty,ipy
1000 continue
   if(chkperm) then
! wait for checking for errors ....
!      write(*,*)'Press return'
!      read(*,297)ch1
297   format(a)
   endif
!    write(*,1001)gx%bmperr,(phres%gval(i,1),i=1,4)
!    write(*,1002)(phres%dgval(1,i,1),i=1,3)
!    write(*,1003)(phres%d2gval(i,1),i=1,6)
1001 format('calcg g: ',i5,4(1PE15.7))
1002 format('calcg dg:  ',3(1PE15.7))
1003 format('calcg d2g: ',6(1PE11.3))
   return
 end subroutine calcg_internal

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine setendmemarr(lokph,ceq)
! stores the pointers to all ordered and disordered endmemners in arrays
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
   character name*24,blaj*70
   double precision kappa,napfu,t,p,rt,g,v,s,h,u,f,cp,alpha
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
   rt=globaldata%rgas*T
   lokcs=lokres
! returned values: G, G.T=-S, G.P=V, G.T.T=-Cp/T G.T.P=V*alpha, G.P.P=-V*kappa
! all divided by RT and per mole formula unit of phase
! G=H-TS, F=U-TS, H=U+PV, S=-G.T, V=G.P
! H=G+TS=G-T*G.T, U=H-PV=(G-T*G.T)-P*G.P, CP=-T*G.T.T
! alpha= 1/V*V.T = G.T.P/V, kappa = -1/V*V.P = -G.P.P/V
   G=rt*ceq%phase_varres(lokcs)%gval(1,1)
!    write(*,5)'tabder 2: ',rt,G
   S=-rt*ceq%phase_varres(lokcs)%gval(2,1)
   V=rt*ceq%phase_varres(lokcs)%gval(3,1)
   H=G+T*S
   U=H-P*V
   F=U-T*S
   CP=-T*rt*ceq%phase_varres(lokcs)%gval(4,1)
   if(V.ne.zero) then
      alpha=rt*ceq%phase_varres(lokcs)%gval(5,1)/V
      kappa=rt*ceq%phase_varres(lokcs)%gval(6,1)/V
   else
      alpha=zero
      kappa=zero
   endif
   write(kou,100)napfu,T,P,G
100 format(/'Per mole FORMULA UNIT of the phase, ',1pe12.4,' atoms'/&
         'at T= ',F8.2,' K and P= ',1PE13.6,' Pa',/ &
         'Gibbs energy J/mol ',28('.'),1Pe16.8)
   write(kou,102)F,H,U,S,V,CP,alpha,kappa
102 format('Helmholtz energy J/mol ',24('.'),1PE16.8 &
        /'Enthalpy J/mol ',32('.'),1PE16.8 &
        /'Internal energy J/mol ',25('.'),1PE16.8 &
        /'Entropy J/mol/K ',31('.'),1PE16.8 &
        /'Volume m3 ',37('.'),1PE16.8 &
        /'Heat capacity J/mol/K ',25('.'),1PE16.8 &
        /'Thermal expansion 1/K ',25('.'),1PE16.8 &
        /'Bulk modulus 1/Pa ',29('.'),1PE16.8)
   tnk=phlista(lokph)%tnooffr
   ll=1
   kk1=0
   kk2=phlista(lokph)%nooffr(ll)
   dy1loop: do while(kk1.le.tnk)
      kk1=kk1+1
      if(kk1.gt.kk2) then
!          write(*,11)'tabder 2: ',kk1,kk2,ll,tnk,nsl
!11 format(a,10i3)
         ll=ll+1
         if(ll.gt.nsl) exit
         kk2=kk2+phlista(lokph)%nooffr(ll)
      endif
      if(phlista(lokph)%nooffr(ll).eq.1) then
!          write(*,*)'tabder 1: ',kk1,kk2,ll,tnk
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
      write(kou,120)rt*ceq%phase_varres(lokcs)%dgval(1,kk1,1),&
           rt*(ceq%phase_varres(lokcs)%dgval(1,kk1,1)-&
           T*ceq%phase_varres(lokcs)%dgval(2,kk1,1)),&
           rt*ceq%phase_varres(lokcs)%dgval(2,kk1,1),&
           rt*ceq%phase_varres(lokcs)%dgval(3,kk1,1)
120    format(5x,'G ',40('.'),1PE16.8, &
           /5x,'H ',40('.'),1PE16.8, &
           /5x,'G.T ',38('.'),1PE16.8, &
           /5x,'G.P ',38('.'),1PE16.8)
      kk3=kk1
      kk4=kk2
      ll2=ll
      write(kou,150)
150 format(5x,'Second partial derivative of Gibbs energy with respect to also')
      dy2loop: do while(kk3.le.tnk)
         if(phlista(lokph)%nooffr(ll2).gt.1) then
            write(kou,160)name(1:len_trim(name)),ll2, &
                 rt*ceq%phase_varres(lokcs)%d2gval(ixsym(kk1,kk3),1)
160          format(10x,a,' in sublattice ',i2,18('.'),1PE16.8)
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
!       write(*,*)'tabder 7A: ',kk1,kk2
   enddo dy1loop
900 continue
!    write(*,*)'tabder 7B: ',kk2
!    write(*,*)'tabder: ',rt,rt*phase_varres(lokcs)%gval(1,1)
1000 continue
   return
 end subroutine tabder

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine cgint(lokpty,moded,vals,dvals,d2vals,gz,ceq)
! calculates an excess parameter that can be composition dependent
! gz%yfrem are the site fractions in the end member record
! gz%yfrint are the site fractions in the interaction record(s)
! lokpty is the property index
! moded=0 means only G, =1 G and dG/dy, =2 all
   implicit none
   integer moded
   TYPE(gtp_property), pointer :: lokpty
   TYPE(gtp_parcalc) :: gz
   double precision vals(6),dvals(3,gz%nofc)
   TYPE(gtp_equilibrium_data) :: ceq
!\end{verbatim}
! temporary data like gz%intlevel, gz%nofc etc
   double precision d2vals(gz%nofc*(gz%nofc+1)/2),valtp(6)
   double precision vv(0:2),fvv(0:2)
   integer lfun,jdeg,jint,jl
   double precision rt,dx0,dx,dx1,dx2,ct,fvs
   double precision, parameter :: onethird=one/3.0D0,two=2.0D0
! zeroing 5 iq, and various vals, dvals and d2vals
   gz%iq=0
   vals=0
   dvals=0
   d2vals=0
   rt=gz%rgast
   if(lokpty%degree.eq.0) then
!----------------------------------------------------------------------
! no composition dependence
      lfun=lokpty%degreelink(0)
      call eval_tpfun(lfun,gz%tpv,vals,ceq%eq_tpres)
      if(gx%bmperr.ne.0) goto 1000
      if(lokpty%proptype.eq.1) then
         vals=vals/rt
      endif
      goto 1000
   endif
   intlev: if(gz%intlevel.eq.1) then
!----------------------------------------------------------------------
! plain binary Redlich Kister. gz%endcon can be wildcard, i.e. negative
! but for the moment give error message in that case
! A binary wildcard excess means y_A ( 1 - y_A) * L_A*
! most naturally gz%intcon(1) would be negative
      gz%iq(1)=gz%endcon(gz%intlat(1))
      gz%iq(2)=gz%intcon(1)
      if(gz%iq(1).lt.0 .or.gz%iq(2).lt.0) then
! composition dependent wildcard interaction not implemented
! y(1-y) ( L0 + (2y-1) L1 + (2y-1)**2 L2 + ....)
         gx%bmperr=4031; goto 1000
      endif
      dx0=gz%yfrem(gz%intlat(1))-gz%yfrint(1)
      vals=zero
      dx=one
      dx1=zero
      dx2=zero
      RK: do jdeg=0,lokpty%degree
         lfun=lokpty%degreelink(jdeg)
         call eval_tpfun(lfun,gz%tpv,valtp,ceq%eq_tpres)
         if(gx%bmperr.ne.0) goto 1000
         if(lokpty%proptype.eq.1) then
! property type 1 is G and should be normalized by RT
            valtp=valtp/rt
         endif
! no composition derivative.  if moded=0 only G, =1 G+G.Y, =2 all
         vals=vals+dx*valtp
         noder5: if(moded.gt.0) then
! first derivatives
            do jl=1,3
               dvals(jl,gz%iq(1))=dvals(jl,gz%iq(1))+dx1*valtp(jl)
               dvals(jl,gz%iq(2))=dvals(jl,gz%iq(2))-dx1*valtp(jl)
            enddo
! second derivatives
            if(moded.gt.1) then
               d2vals(ixsym(gz%iq(1),gz%iq(1)))=&
                    d2vals(ixsym(gz%iq(1),gz%iq(1)))+dx2*valtp(1)
               d2vals(ixsym(gz%iq(1),gz%iq(2)))=&
                    d2vals(ixsym(gz%iq(1),gz%iq(2)))-dx2*valtp(1)
               d2vals(ixsym(gz%iq(2),gz%iq(2)))=&
                    d2vals(ixsym(gz%iq(2),gz%iq(2)))+dx2*valtp(1)
            endif
         endif noder5
! next power of dx
         dx2=(jdeg+1)*dx1
         dx1=(jdeg+1)*dx
         dx=dx*dx0
      enddo RK
   elseif(gz%intlevel.eq.2) then !intlev
!----------------------------------------------------------------------
! ternary interaction
      ternary: if(gz%intlat(1).eq.gz%intlat(2)) then
! Ternary composition dependent interaction in same sublattice, Hillert form.
! The idea is that the sum of vv is always unity even in higher order systems.
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
               valtp=valtp/rt
            endif
! function value
            vals=vals+vv(jint)*valtp
            noder6: if(moded.gt.0) then
! first derivatives
               do jl=1,3
                  dvals(jl,gz%iq(1))=dvals(jl,gz%iq(1))+fvv(0)*valtp(jl)
                  dvals(jl,gz%iq(2))=dvals(jl,gz%iq(2))+fvv(1)*valtp(jl)
                  dvals(jl,gz%iq(3))=dvals(jl,gz%iq(3))+fvv(2)*valtp(jl)
               enddo
! there is no contribution to the second derivatives from this interaction
            endif noder6
            fvs=fvv(2)
            fvv(2)=fvv(1)
            fvv(1)=fvv(0)
            fvv(0)=fvs
         enddo terloop
      else
!----------------------------------------------------------------------
! composition dependent reciprocal interactions here
!         write(*,32)gz%intlevel
!32       format('Property record at intlevel ',i3/&
!              'Composition dependent reciprocal interaction not implemented')
         gx%bmperr=4078; goto 1000
      endif ternary
   else !intlev
!----------------------------------------------------------------------
! much more to be implemented here .... ionic liquids, quasichem ....
      write(*,*)'cgint 99 not implemented: ',gz%intlevel,lokpty%degree
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
! uses information about sublattices etc in phlista(lokph)
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
!      write(*,*)'Tying to pop from an empty PY stack'
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
!\end{verbatim}
   TYPE(gtp_fraction_set), pointer :: disrec
   TYPE(gtp_phase_varres), pointer :: phord
   TYPE(gtp_phase_varres), pointer :: phdis
   logical ordered
! minimum difference in site fraction to be set as ordered
   double precision, parameter :: yminord=1.0D-10
   integer lokdis,is,kk
!
!   write(*,*)'entering calc_disfrac'
!   disrec=phord%disfra
!   lokdis=disrec%varreslink
!   phdis=>disrec%phdapointer
! this is the record with the ordered constitution
   phord=>ceq%phase_varres(lokcs)
! this is a record within the ordered constitution record for disordered fracs
   disrec=>phord%disfra
! to find the varres record with disordered fractions use varreslink
! this is the index to the phase_varres record with the ordered fractions ???
   lokdis=disrec%varreslink
   phdis=>ceq%phase_varres(lokdis)
!   write(*,*)'calc_disfrac 1A'
! check that some values are accessable
!   write(*,*)'calc_disfra phase index: ',phord%phlink
!   write(*,*)'calc_disfra disordered sublattices: ',disrec%ndd
!   write(*,*)'calc_disfra ordered and disordered records: ',lokcs,lokdis
!   write(*,*)'calc_disfra phase index via disordred record: ',phdis%phlink
!   write(*,*)'calc_disfrac 1B'
   phdis%yfr=zero
!   write(*,*)'disfrac 1: ',disrec%tnoofyfr
   do is=1,disrec%tnoofyfr
      phdis%yfr(disrec%y2x(is))=&
           phdis%yfr(disrec%y2x(is))+disrec%dxidyj(is)*phord%yfr(is)
!      write(*,77)'disfrac 2: ',is,disrec%y2x(is),phdis%yfr(disrec%y2x(is)),&
!           disrec%dxidyj(is),phord%yfr(is)
77    format(a,2i3,3(1pe12.4))
   enddo
!   write(*,*)'calc_disfrac 2'
! check if phase is really ordered, meaning that the disordered fractions
! are equal to the ordered ones
   ordered=.false.
   do is=1,disrec%tnoofyfr
      if(abs(phdis%yfr(disrec%y2x(is))-&
           phord%yfr(is)).gt.yminord) ordered=.true.
   enddo
!   write(*,*)'calc_disfrac 3'
   if(.not.ordered) then
! if this bit set one will not calculate the ordered part of the phase
      phord%status2=ibclr(phord%status2,csorder)
!      write(*,*)'calc_disfrac: disordered, clear ordered bit',lokph
   else
! bit must be cleared as it might have been set at previous call
      phord%status2=ibset(phord%status2,csorder)
!      write(*,*)'calc_disfrac: ordered, set ordered bit',lokph
   endif
!   write(*,*)'calc_disfrac 4'
! copy these to the phase_varres record that belongs to this fraction set
! a derivative dGD/dyj = sum_i dGD/dxi * dxidyj
! where dGD/dxi is dgval(1,y2x(j),1) and dxidyj is disrec%dxidyj(j)
! because each y constituent contributes to only one disordered x fraction
1000 continue
   return
! G(tot)    = GD(xdis)+(GO(yord)-GO(yord=xdis))
! G(tot).yj = dGD(xdis).dxi*dxdyj + GO.yj - GO.yj ... 
 end subroutine calc_disfrac

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine disordery(phvar,ceq)
! sets the ordered site fractions in FCC and other order/disordered phases
! equal to their disordered value in order to calculate and subtract this part
! phvar is index to phase_varres for ordered fractions
   implicit none
   TYPE(gtp_phase_varres), pointer :: phvar
   TYPE(gtp_equilibrium_data) :: ceq
!\end{verbatim}
   TYPE(gtp_fraction_set), pointer :: disrec
   TYPE(gtp_phase_varres) :: phdis
   integer lokdcs,kk,ll,is,nis,nsl
! find disordered fractions
   lokdcs=phvar%disfra%varreslink
   disrec=>phvar%disfra
!   write(*,*)'disordery: ',disrec%latd,disrec%nooffr(1),lokdcs
   phdis=ceq%phase_varres(lokdcs)
!   write(*,*)'disordery: ',ceq%xconv
!   write(*,*)'disordery: ',phdis%yfr(1)
!   phdis=>ceq%disrec%phdapointer
! copy fractions, loop through all ordered sublattices in phvar
! and store fraction from lokdis
   kk=0
! here copy: 
! y(ord,1,1)=y(dis,1); y(ord,1,2)=y(dis,2); y(ord,1,3)=y(dis,3); 
! y(ord,2,1)=y(dis,1); y(ord,2,2)=y(dis,2); y(ord,2,3)=y(dis,3); 
   do ll=1,disrec%latd
      do is=1,disrec%nooffr(1)
         kk=kk+1
         phvar%yfr(kk)=phdis%yfr(is)
      enddo
   enddo
   if(disrec%ndd.eq.2) then
! one can have 2 sets of ordered subl. like (Al,Fe)(Al,Fe)...(C,Va)(C,Va)...
      nis=disrec%nooffr(1)
      nsl=size(phvar%sites)
!      write(*,*)'dy: ',nis,kk,disrec%latd,nsl,disrec%nooffr(2)
      do ll=disrec%latd+1,nsl
         do is=1,disrec%nooffr(2)
            kk=kk+1
            phvar%yfr(kk)=phdis%yfr(nis+is)
         enddo
      enddo
   endif
1000 continue
   return
 end subroutine disordery

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine scale_phase_amounts(antot,ceq)
! multiply all phase amounts in moles by antot.  Probably redundant
! NOTE: MASS updated with same factor, is that correct?
! is this routine used? needed?
   implicit none
   double precision antot
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer iph,lokph,ics,lokcs
   do iph=1,noofph
      lokph=phases(iph)
      do ics=1,phlista(lokph)%noofcs
!         lokcs=phlista(iph)%cslink
         lokcs=phlista(iph)%linktocs(ics)
         if(ceq%phase_varres(lokcs)%amfu.gt.zero) then
            ceq%phase_varres(lokcs)%amfu=&
                 ceq%phase_varres(lokcs)%amfu*antot
            ceq%phase_varres(lokcs)%netcharge=&
                 ceq%phase_varres(lokcs)%netcharge*antot
         endif
      enddo
   enddo
1000 continue
   return
 end subroutine scale_phase_amounts

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

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
   integer, dimension(4) :: indices
   double precision, dimension(maxel) :: ani,abi,xset,wset
   double precision btot,bni(maxel),mass,h298,s298,xxx,xsum,wsum
   double precision sumwdivm,anisum,abisum,restmass,divisor,dividend,abtot
   TYPE(gtp_condition), pointer :: current,last
   character encoded*16,actual_arg(1)*16,elsym*2,elname*16,refstat*16
   integer nox,now,nc,jl,iref,iunit,ip,idf,ie,more,numberest,istv
   logical allmassbal
!
   ani=zero; abi=zero; xset=zero; wset=zero
   antot=zero; abtot=zero
   xsum=zero; wsum=zero
   anisum=zero; abisum=zero
   nox=0; now=0
!
   last=>ceq%lastcondition
   if(.not.associated(last)) then
      gx%bmperr=4143; goto 1000
   endif
   current=>last
   nc=0
   allmassbal=.TRUE.
100 continue
      current=>current%next
! ignore inactive conditions
      if(current%active.ne.0) goto 300
! ignore conditions with several terms
      if(current%noofterms.gt.1) goto 300
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
         xxx=evaluate_svfun(current%symlink1,actual_arg,1,ceq)
      else
         xxx=current%prescribed
      endif
!      write(*,17)'massbal: ',encoded,istv,indices,iunit,iref,xxx
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
!            write(*,*)'N with 2 indices illegal in this case'
            gx%bmperr=4179; goto 1000
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
            write(*,*)'B with 2 indices illegal'
            gx%bmperr=4179; goto 1000
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
      gx%bmperr=4151; goto 1000
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
         if(antot.eq.zero .and. abtot.eq.zero) goto 1105
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
               write(*,*)'Missing condition for two elements.'
               gx%bmperr=0; goto 1000
            endif
            restmass=mass
            numberest=ie
         endif
      enddo
      if(numberest.eq.0) then
         write(*,*)'Error - condition on all elements and N??'
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
!         write(*,520)'nrest: ',numberest,divisor,dividend,ani(numberest),&
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
!      write(*,*)'Cannot handle condition on total mass'
      gx%bmperr=4180
   elseif(xsum.eq.zero .and. wsum.eq.zero) then
! just N(i)= and B(i)=, no N= nor B= and no X nor W, No rest element
!      write(*,520)'N(i): ',0,anisum,(ani(j),j=1,noel())
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
            write(*,*)'mass balance error: ',ie
            gx%bmperr=4181; goto 1000
         endif
      enddo
   else
! any other combination of conditions ....
      write(*,*)'Cannot handle these massbalance conditions'
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
   write(*,*)'Two mass balance conditions for same element',ie
   gx%bmperr=4183; goto 1000
1105 continue
   write(*,*)'One component without condition'
   gx%bmperr=4181; goto 1000
1110 continue
   write(*,*)'Both N and B cannot be set'
   gx%bmperr=4184; goto 1000
!
 end subroutine extract_massbalcond

