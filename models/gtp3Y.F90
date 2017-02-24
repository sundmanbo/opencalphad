!
! gtp3Y included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>      16. Grid minimizer
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine global_gridmin(what,tp,xknown,nvsph,iphl,icsl,aphl,&
      nyphl,cmu,ceq)
!
! Starting rewriting 2017-02-01
!
! finds a set of phases that is a global start point for an equilibrium 
! calculation at T and P values in tp and known mole fraction in xknown
! It is intentional that this routine is independent of current conditions
! It returns: nvsph stable phases, list of phases in iphl, amounts in aphl, 
! nyphl(i) is redundant, cmu are element chemical potentials of solution
! WHAT determine what to do with the results, 0=just return solution,
! 1=enter stable set and constitution of all phases in gtp datastructure
! and create composition sets if necessary (and allowed)
! what=-1 will check if any gridpoint is below current calculated equilibrium
   implicit none
! nyphl(j) is the start position of the constitiuent fractions of phase j in
   integer, dimension(*) :: iphl,nyphl,icsl
   integer what,nvsph

   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision, dimension(2) :: tp
! cmu(1..nrel) is the chemical potentials of the solution
   double precision, dimension(*) :: xknown,aphl,cmu
!\end{verbatim}
! yarr is used for fractions in the call to generate_grid
   double precision, dimension(maxconst) :: yarr
   integer, parameter :: maxgrid=400000,maxy=2000,maxph=500
   integer :: starttid,endoftime
   real finish2
! removed yphl as argument as it not needed outside global_gridmin
! dimensioning can be problematic if many phases with many constituents as it
! contains all constituent fractions of all gridpoints (before merge)
!   double precision, dimension(10* amounts of maxconst) :: yphl
   double precision, dimension(10*maxconst) :: yphl
   double precision amount,sum,gmax
   integer ibias,ics,ics2,icsno,icsx,ie,iph,iv,j1,j2,jp,kkz,kp,kph,jbias
   integer lokcs,lokph,mode,ng,nocsets,noofgridpoints,nr,nrel,nrph,ny,nyz
   integer preveq
! nphl(iph) is last gridpoint that belongs to  phase iph, nphl(0)=0
! xarr(nrel,i) is the composition of gridpoint i
! garr(i) is the Gibbs energy of gridpoint i
! jgrid(j) is a gridpoint in the solution
! phfrac(j) is the amount of the phase of that gridpoint
! ngrid deleted
   integer, dimension(0:maxph) :: nphl
   integer, dimension(maxel) :: jgrid
   real garr(maxgrid),starting,finished
   real, dimension (:,:), allocatable :: xarr
   real, dimension (maxel,maxel) :: xsol
   double precision, dimension(maxel) :: phfrac,phsave,xdum
   double precision qq(5),savetp(2),xbase,totam
   integer, dimension(maxph) :: iphx
   character name1*24
!   integer idum(*)
! debug
   logical trace,toomany
! sort phases depending on number of gridpoints
   integer, dimension(:), allocatable :: gridpoints,phord,starttup
! pph is set to number of phases participating, some may be suspended
   integer pph,zph,nystph,order(maxel),tbase,qbase,wbase,jj
!
!   write(*,*)'3Y in global_gridmin'
   if(btest(globaldata%status,GSNOGLOB)) then
      write(*,*)'3Y Grid minimization not allowed'
      gx%bmperr=4173; goto 1000
   endif
   call cpu_time(starting)
   call system_clock(count=starttid)
   nphl=0
! Trace turn on output of grid on a file ocgrid.dat
!   trace=.true.
   toomany=.false.
   trace=.FALSE.
   if(trace) write(*,*)'3Y Trace set TRUE'
   savetp=ceq%tpval
   ceq%tpval=tp
   nrph=noph()
   if(nrph.gt.maxph) then
! too many phases
      write(*,*)'3Y Too many phases for gridmin'
      gx%bmperr=4344; goto 1000
   endif
   nrel=noel()
   sum=zero
! problem that extract_massbalcond did not object to condition x(fcc,a)=
   do ie=1,nrel
      if(xknown(ie).le.zero .or. xknown(ie).ge.one) then
         if(.not.btest(globaldata%status,GSNOTELCOMP)) then
            write(*,*)'3Y Gridmin cannot handle these composition conditions'
            gx%bmperr=4174; goto 1000
         else
! we have other components than elements, fractions can be negative and >1
            write(kou,10)
10          format('3Y Trying to use gridmininmizer whith other components',&
                 ' than the elements'/'   Can give warnings and error messages')
            gx%bmperr=4174; goto 1000
         endif
      endif
      sum=sum+xknown(ie)
   enddo
   if(ocv()) write(*,12)'3Y gridmin: ',sum,(xknown(ie),ie=1,nrel)
12 format(a,1pe12.4,10(f8.4))
   if(abs(sum-one).gt.1.0D-8) then
      write(*,*)'3Y Sum of fractions not unity when calling global_gridmin'
      gx%bmperr=4174; goto 1000
   endif
   kp=1
   pph=0
   allocate(gridpoints(nrph))
   allocate(phord(nrph))
   ggloop: do iph=1,nrph
! include all phases with any composition set entered (but only once!)
      do ics=1,noofcs(iph)
! new: -3 suspended, -2 dormant, -1,0,1 entered, 2 fixed
! ignore phases whith no composition set entered
! If a phase+compset FIX one should never be here as conditions wrong
         if(test_phase_status(iph,ics,amount,ceq).gt.PHDORM) then
            pph=pph+1
            iphx(pph)=iph
            cycle ggloop
         endif
      enddo
   enddo ggloop
! we will generate a grid for pph phases, the phase index for phase 1..pph
! is in iphx(1..pph)
! always allocate a grid for maxgrid points
   allocate(xarr(nrel,maxgrid))
   gmax=zero
! just to be sure
   nphl(0)=0
!
   phloop: do zph=1,pph
! for phase iphx(zph) the gridpoints will be stored from position nphl(zph-1)+1
! ng should be set to number of remaining points, ny and yphl is not used
      iv=nphl(zph-1)+1
      ng=maxgrid-iv
! this call will calculate gridpoints in phase zph, that may take time ...
! ng is set to remaining dimension of garr, on return the number of generated
!    gridpoints, returned as xarr composition of these and
! ny and yarr not used here
      call generic_grid_generator(0,iphx(zph),ng,nrel,xarr(1,iv),garr(iv),&
           ny,yarr,gmax,ceq)
      if(gx%bmperr.ne.0) then
         write(*,*)'3Y grid error ',zph,gx%bmperr
         exit phloop
      endif
!      write(*,*)'3Y done generic_grid: ',zph,iphx(zph),ng
! nphl(zph) is last gridpoint in phase zph
      nphl(zph)=nphl(zph-1)+ng
   enddo phloop
   if(gx%bmperr.ne.0) goto 1000
! We should add the current set of stable phases in the grid if we have made
! a successful calculation
!   if(.not.btest(ceq%status,EQNOEQCAL)) then
! add the current set of stable phases and their constitution as gridpoints
!      write(*,*)'3Y Not yet adding current stable phases as gridpoints',&
!           what
      preveq=0
!   endif
! we may be generating a list with all gridpoints ...
   if(trace) then
      write(*,*)'3Y Closing gridgen.dat'
      close(33)
   endif
   call system_clock(count=endoftime)
   call cpu_time(finished)
! kp set to total number of grispoints in all phases
   kp=nphl(pph)
   noofgridpoints=kp
! If WHAT is -1 then just compare all gridpoints with plane defined by
! the chemical potentials cmu to see if any is below.
! If so insert the gridpoint furtherst below the plane and set WHAT 10*iph+ics
!   write(*,*)'3Y global_gridmin what: ',what
   if(what.eq.-1) then
      call gridmin_check(nystph,kp,nrel,xarr,garr,xknown,nphl,pph,&
           cmu,yarr,iphx,ceq)
      goto 1000
   endif
!-----------------------------------------------
!    write(*,109)ngrid(pph),finished-starting,endoftime-starttid
109 format('3Y Calculated ',i6,' gridpoints in ',1pe12.4,' seconds, ',&
         i7,' clockcycles')
! find the minimum of nrel gridpoints among the kp-1 gridpoint
! for current overall composition, xknown
!    write(*,*)'3Y globm 4: ',kp,garr(kp),xarr(1,kp)
!   phfrac=zero
! start with all elements having chemical potential equal to gmax
   cmu(1)=gmax
   if(ocv()) write(*,*)'3Y Finding the gridpoints for the minimum: ',kp
!   write(*,*)'3Y Finding the gridpoints for the minimum: ',kp
   call find_gridmin(kp,nrel,xarr,garr,xknown,jgrid,phfrac,cmu,trace)
   if(gx%bmperr.ne.0) goto 1000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   write(*,*)'3Y gridpoints: ',kp
!   gx%bmperr=4399; goto 1000
! The solution with nrel gridpoints are in jgrid, the amount of each in phfrac
! We later want the phases in ascending order and as the gridpoints are
! in ascending order of the phases we sort the gridpoints (and amounts)
! There must be one gridpoint per component (element)
!   write(*,62)(jgrid(jp),jp=1,nrel)
62 format('3Y Gridp: ',10i6)
   call sortin(jgrid,nrel,order)
   do nyz=1,nrel
      phsave(nyz)=phfrac(order(nyz))
   enddo
   phfrac=phsave
! get the phase and constitution for each
   nyz=1
!  Extracting constititution of the gridpoints in the solution
   if(trace) then
      write(31,745)
745   format(/'Solution: ')
   endif
   do jp=1,nrel
      iphl(jp)=0
   enddo
   solloop: do jp=1,nrel
! jgrid(jp) is a grid point in the solution, find which phase it is
! and its constituent fractions
      mode=jgrid(jp)
713   continue
      do zph=1,pph
!          write(*,*)'3Y mode and ibias 1: ',mode,ibias
! nphl(zph) is the last gridpoint of phase zph, nphl(0) is 0
         if(mode.le.nphl(zph)) then
            mode=mode-nphl(zph-1)
            ibias=nphl(zph-1)
            goto 315
         else
         endif
      enddo
! nphl(pph) is number of generated gridpoints
      if(mode-nphl(pph).le.preveq) then
! gridpoint outside generated gridpoints, should be from previous solution
         write(*,*)'3Y previous stable phase included in solution',mode,preveq
      endif
! gridpoint outside range should never occur
      write(*,*)'3Y gridpoint outside range ',jgrid(jp),nphl(pph)
! It means element je=jgrid(jp)-nphl(pph) has no chemical potential
! and possibly no composition.  Find the gripoint with max of with this 
! component and add a small amont if it to avoid that an element has 
! no phase in which is can dissolve ...
      qbase=jgrid(jp)-nphl(pph)
      xbase=zero
      wbase=0
      do tbase=1,nphl(pph)
         if(xarr(qbase,tbase).gt.xbase) wbase=tbase
      enddo
      if(wbase.eq.0) then
! we have failed to find a gridpoint with this element
         gx%bmperr=4147; goto 1000
      else
         write(*,*)'3Y using point: ',wbase
         phfrac(jp)=1.0D-4
         mode=wbase
         goto 713
      endif
315   continue
      jbias=ibias
! this call is to obtain the constitution of a phase in the solution
! mode gives in grid point index in phase iphx(zph), ibias irrelevant (?)
! NOTE jbias is changed by subroutine ??
! ny is number of constituent fractions, yarr have the constituent fractions
!      write(*,317)'3Y point: ',mode,jp,iphl(jp),(iphl(nr),nr=1,jp)
      call generic_grid_generator(mode,iphx(zph),jbias,nrel,xarr,garr,&
           ny,yarr,gmax,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,317)'3Y after0: ',mode,jp,nyz,ibias,jbias,iphl(jp),&
!           (iphl(nr),nr=1,jp)
      iphl(jp)=iphx(zph)
      aphl(jp)=phfrac(jp)
      nyphl(jp)=ny
! copy the constitution of all gridpoints to yphl, needed for possible merge
      do ie=1,ny
         yphl(nyz+ie-1)=yarr(ie)
      enddo
      nyz=nyz+ny
!      write(*,317)'3Y after1: ',mode,jp,nyz,ibias,jbias,iphl(jp),&
!           (iphl(nr),nr=1,jp)
317   format(a,i6,4i4,i3,20i3)
! finally copy the mole fractions to xsol, needed for possible merging
      do ie=1,nrel
         xsol(ie,jp)=xarr(ie,mode+ibias)
      enddo
      if(trace) then
         write(31,750)jp,jgrid(jp),iphl(jp),aphl(jp),(xsol(ie,jp),ie=1,nrel)
         write(31,760)(yphl(ie),ie=nyz-ny,nyz-1)
750      format('Point: ',i2,', gridpoint: ',i5,' phase ',i3,&
              ' amount: ',1pe12.4,', Mole fractions:'/9(0pF8.5))
760      format('Constitution:'/9(0pF8.5))
      endif
   enddo solloop
! we have now start values from the gridminimizer
   if(trace) then
      write(*,*)'3Y Closing ocgrid.dat file'
      close(31)
   endif
! there must be as many gridpoints (phases) as there are elements
   nvsph=nrel
   nr=nrel
!   if(.not.btest(globaldata%status,GSNOMERGE)) then
! For the moment we will only merge grid points in the gas phase
   call merge_gridpoints(nr,iphl,aphl,nyphl,yphl,trace,nrel,xsol,cmu,ceq)
   if(gx%bmperr.ne.0) goto 1000
!   endif
!-------------------------------------------
! number of gridpoints, nr, may have changed
!   write(*,*)'3Y After merge_gripoints: ',nr,nvsph
   nvsph=nr
! if what=-1 or 0 do nothing more, just exit
   if(what.le.0) goto 1000
!------------------------------------------------------------
! prepare for storing result in the ceq data structure
! zero all phase amounts and driving forces
   do iph=1,nrph
      lokph=phases(iph)
!      lokcs=phlista(lokph)%cslink
      do ics=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(ics)
         ceq%phase_varres(lokcs)%dgm=zero
         ceq%phase_varres(lokcs)%amfu=zero
         ceq%phase_varres(lokcs)%netcharge=zero
         if(ceq%phase_varres(lokcs)%phstate.eq.phentstab) then
! reset status of "entered and stable" to just "entered" 
            ceq%phase_varres(lokcs)%phstate=phentered
         endif
      enddo
   enddo
! store chemical potentials multiplied with RT if what not -1
   ceq%rtn=globaldata%rgas*ceq%tpval(1)
   do ie=1,nrel
!      write(*,*)'3Y grid chemical potential: ',ie,cmu(ie)*ceq%rtn
! do not care about reference state for chempot(2)
      ceq%complist(ie)%chempot(1)=cmu(ie)*ceq%rtn
      ceq%complist(ie)%chempot(2)=cmu(ie)*ceq%rtn
   enddo
! set driving force 0 for stable phases
   do ie=1,nvsph
      call set_driving_force(iphl(ie),1,zero,ceq)
      if(gx%bmperr.ne.0) goto 1000
   enddo
! store the most favourable constitution of the metastable phase
!  write(*,29)'3Y set constitution of metastable phases',pph,(iphx(ie),ie=1,pph)
!   write(*,29)'3Y gps: ',(nphl(ie),ie=0,pph)
29 format(a,(20i5))
   call set_metastable_constitutions2(pph,nrel,nphl,iphx,xarr,garr,&
        nvsph,iphl,cmu,ceq)
   if(gx%bmperr.ne.0) goto 1000
   write(*,*)'3Y Constitution of metastable phases set'
! maybe more composition sets needed
   do ie=1,nvsph
      icsl(ie)=0
   enddo
   nocsets=0
!   write(*,*)'3Y before loop1: ',nvsph,ceq%eqname
! loop for all gripoints to store them in composition sets
   loop1: do j1=1,nvsph
      if(icsl(j1).eq.0) then
! if non-zero a composition set has already been assigned in loop2
         icsl(j1)=1
         icsx=1
         loop2: do j2=j1+1,nvsph
            nextgridp: if(iphl(j1).eq.iphl(j2)) then
! one more composition set needed, does it exist?
               icsx=icsx+1
               ics2=icsx
!               write(*,31)'3Y compset needed for phase ',j1,j2,iphl(j1),ics2
31             format(a,10i4)
               call get_phase_compset(iphl(j1),ics2,lokph,lokcs)
               newset: if(gx%bmperr.ne.0) then
! there is no such composition set, is automatic creation allowed?
! NOTE: there is a EQNOACS bit also???
                  if(btest(globaldata%status,GSNOACS) .or. &
                       btest(ceq%status,EQNOACS)) then
                     write(*,*)'3Y Not allowed to create composition sets'
                     gx%bmperr=4177; goto 1000
                  endif
                  gx%bmperr=0
! >>>>>>>>>>>>>>>>>>><
! BEWARE >>> not only must this be done in all threads at the same time
! one must also avoid that it is done when some thread is working on a set
! of phase+composition sets transformed to EQCALC arrays.  If so the
! indices to lokcs etc will be incorrect ... ???
! I think OMP has "secure" points where the treads can be stopped to wait
! <<<<<<<<<<<<<<<<<<<<<<
                  kph=iphl(j1)
!                  write(*,*)'3Y extra composition set for phase: ',kph,j1,j2
! It must be done in all equilibrium records, no equilibrium record needed!!!
! one must be careful with the status word when creating comp.sets
                  call enter_composition_set(kph,'    ','AUTO',icsno)
                  if(gx%bmperr.ne.0) then
!                   write(*,*)'3Y Error entering composition set ',j1,gx%bmperr
                     if(gx%bmperr.eq.4092) then
! skip entering this set, it may work anyway ...
                        if(.not.toomany) then
                           write(kou,298)iphl(j1)
298                        format('Cannot enter enough composition sets',&
                                ' for phase',i4,' but gridmin struggles on')
                           toomany=.true.
                        endif
                        gx%bmperr=0
! to avoid later trouble we should mark there is no compset for this gridp!!
                        iphl(j2)=-kph
                        icsl(j2)=-1
                        cycle loop2
                     else
                        goto 1000
                     endif
                  endif
                  call get_phase_compset(kph,icsno,lokph,lokcs)
                  if(gx%bmperr.ne.0) goto 1000
                  ceq%phase_varres(lokcs)%status2=&
                       ibset(ceq%phase_varres(lokcs)%status2,CSAUTO)
                  nocsets=nocsets+1
                  icsl(j2)=icsno
               else
! here we should check which composition set that should have which 
! constitution for example one fcc is metallic and another is cubic carbide
                  call get_phase_name(iphl(j1),ics2,name1)
                  icsl(j2)=ics2
! check if the composition set is fix (2), dormant (2) or suspended (3)
                  kkz=test_phase_status(iphl(j1),ics2,amount,ceq)
! old kkz=2 means fix
                  if(kkz.eq.PHFIXED) then
                     write(*,*)'3Y Global minimization with fix phase!'
                     gx%bmperr=4346; goto 1000
                  elseif(kkz.lt.PHENTUNST) then
                     write(*,*)'3Y Changing status for phase ',name1
                  endif
! this means status entered PHSTATE
                  ceq%phase_varres(lokcs)%phstate=0
               endif newset
            endif nextgridp
         enddo loop2
      endif
   enddo loop1
   if(nocsets.gt.0) then
      if(.not.btest(globaldata%status,GSSILENT)) then
         write(*,*)'3Y Composition set(s) created: ',nocsets
      endif
   endif
! Above one should consider if some user created compsets are dedicated to
! certain cases (MC carbides or L1_2 ordered).  These should have
! a default constitution and CSDEFCON set)
! finally store stable phase amounts and constitutions into ceq%phase_varres
   j1=1
   allocate(starttup(nvsph))
   call extract_massbalcond(ceq%tpval,xdum,totam,ceq)
   if(gx%bmperr.ne.0) goto 1000
   sum=zero
   do iph=1,nvsph
      sum=sum+aphl(iph)
   enddo
   gmax=zero
! If there is a gas nvsph may be less than number of elements
   ceqstore: do iph=1,nvsph
      if(iphl(iph).lt.0) then
! this gripoint has no composition set because too many gridpoints in same phase
         starttup(iph)=0
         continue
      else
         call set_constitution(iphl(iph),icsl(iph),yphl(j1),qq,ceq)
         if(gx%bmperr.ne.0) goto 1000
         call get_phase_compset(iphl(iph),icsl(iph),lokph,lokcs)
         if(gx%bmperr.ne.0) goto 1000
! This is a bit quiestionable but seems to work
         amount=aphl(iph)/ceq%phase_varres(lokcs)%abnorm(1)
         gmax=gmax+amount
         aphl(iph)=amount
1789     format(a,2i4,5(1pe12.4))
         ceq%phase_varres(lokcs)%amfu=aphl(iph)
         ceq%phase_varres(lokcs)%phstate=PHENTSTAB
         starttup(iph)=ceq%phase_varres(lokcs)%phtupx
         j1=j1+nyphl(iph)
      endif
   enddo ceqstore
!-----------------------------------------------------------------------
! debug listing of tuples at gridpoints
!   write(*,1411)(starttup(iph),iph=1,nvsph)
!1411 format('3Y tupl:',18i4)
!-----------------------------------------------------------------------
! iv is total number of constituent fractions
   iv=j1
! For safty, if any iphl is negative shift all values down for all gridpoints
! I do not think yphl is used any more but ...
   j1=1
   iph=1
870 continue
      if(iphl(iph).lt.0) then
880      format(a,(20i4))
         do kph=iph,nvsph-1
            iphl(kph)=iphl(kph+1)
            icsl(kph)=icsl(kph+1)
            aphl(kph)=aphl(kph+1)
            kkz=nyphl(kph)
            nyphl(kph)=nyphl(kph+1)
         enddo
         iphl(nvsph)=-9
!         write(*,880)'3Y After:  ',(iphl(j2),j2=1,nvsph)
         if(kph.lt.nvsph) then
            do j2=j1,iv
               yphl(j2)=yphl(j2+kkz)
            enddo
            iv=iv-kkz
         endif
! we shifted all down, fewer gridpoints and do not increment iph below
         iph=iph-1
         nvsph=nvsph-1
      endif
      j1=j1+nyphl(iph)
      iph=iph+1
      if(iph.le.nvsph) goto 870
!---------------------------------------
!   write(*,*)'3Y gridpoints: ',nvsph,iph
1000 continue
!   write(*,*)'3Y at 1000: ',phlista(1)%noofcs
! restore tpval in ceq
   ceq%tpval=savetp
   call cpu_time(finish2)
   if(allocated(xarr)) deallocate(xarr)
   if(gx%bmperr.ne.0) then
      ceq%status=ibset(ceq%status,EQFAIL)
   elseif(what.eq.-1) then
      if(nystph.gt.0) what=nystph
   elseif(.not.btest(globaldata%status,GSSILENT)) then
      write(*,1010)noofgridpoints,finish2-starting,&
           endoftime-starttid,ceq%tpval(1)
1010  format('Gridmin: ',i7,' points ',1pe10.2,' s and ',&
           i7,' clockcycles, T=',0pF8.2)
! set the global bit that this is not a full equilibrium
      ceq%status=ibset(ceq%status,EQGRIDCAL)
   endif
! deallocate 
   if(allocated(gridpoints)) then
      deallocate(gridpoints)
      deallocate(phord)
   endif
   if(ocv()) write(*,*)'3Y leaving global_gridmin'
!   write(*,*)'3Y leaving global_gridmin'
   return
 end subroutine global_gridmin

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine generate_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! Different action depending of the value of mode, 
! for mode<0:  (will no longer be used ... )
!    return the number of gridpoints that will be generated for phase iph in ngg
! for mode=0:
!    on entry ngg is dimension of garr
!    on exit ngg is number of generated gridpoints ...
!    return garr(i) gibbs energy and xarr(1,i) the compositions of gridpoint i
! for mode>0:
!    return site fractions of gridpoint mode in yarr, number of fractions in ny
!    iph is phase number, ngg is number of gridpoints, nrel number of elements,
! if mode=0:
!    return xarr mole fractions of gridpoints, garr Gibbs energy of gridpoints,
!    ngg is dimension of garr, gmax maximum G (start value for chem.pot)
! if mode>0:
!   "mode" is a gridpoint of this phase in solution, return number of 
!   constituent fractions in ny and fractions in yarr for this gridpoint
! The current constitution is restored at the end of the subroutine
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*),gmax
!\end{verbatim} %+
!
!   integer idum(*)
   integer lokph,errsave
   double precision, parameter :: yzero=1.0D-12
   integer abrakadabra,i,ibas,ibin,iend,is,iter,je,jend,kend,ll,ls,nend,nsl
! used to save and restore constituent fractions
   double precision ydum(maxconst)
  integer nkl(maxsubl),knr(maxconst),inkl(0:maxsubl),nofy
   double precision, dimension(:), allocatable :: yfra
   double precision sites(maxsubl),qq(5)
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
   integer, dimension(:,:), allocatable :: endm
   integer maxdim
!--------------------------------
! grid is generated by combining end endmembers
! Number of endmemers is N
! For endmember E=1..N set fraction of enmember 
!    0.99*Y_E + 0.01*Y_all                            N of these
!    0.89*Y_E + 0.10*Y_F,F=/=E + 0.01*Y_all           N*(N-1)
!    0.74*Y_E + 0.25*Y_F,F=/=E + 0.01*Y_all           N*(N-1)+N*(N-1)*(N-2)
!             + 0.15*Y_F + 0.1*Y_G,G=/=(E,F) + 0.01*Y_all  (3 or more endmemb)
!    0.61*Y_E + 0.38*Y_F,F=/=E + 0.01*Y_all
!             + 0.25*Y_F + 0.13*Y_G,G=/=(E,F) + 0.01*Y_all (3 or more endmemb)
! added:
!    0.45*Y_E + 0.35*Y_F + 0.19*Y_G + 0.01*Y_all      N*(N-1)*(N-2) (4 or more)
!----- N=2: total 2+2+2+2=8
!----- N>2: total N*(1+(N-1)*(3+2*(N-2))); N=3:33, N=20:
! with 2 endmembers: 2*(1+3)=2*4=8
! (1.00,0.00)
! (0.89,0.11) (0.74,0.26) (0.61,39)
! (0.00,1.00) ...
! with 3 endmembers: 3*11=33 gridpoints
! (1.00,0.00,0.00) 
! (0.89,0.11,0.00)(0.89,0.00,0.11)
! (0.74,0.26,0.00)(0.74,0.00.0.26)(0.74,0.15,0.11)(0.74,0.11,0.15)
! (0.61,0.38,0.00)(0.61,0.00,0.38)(0.61,0.25,0.14)(0.61,0.14,0.25)
! (0.00,1.00,0.00)
! (0.11,0.89,0.00)(0.00,0.89.0.11)
! with 4 endmembers: 
! (0.9925,0.0025,0.0025.0.0025)
! (0.8925,0.1025,0.0025,0.0025) (-,0.0025,0.1025,-) (-,.0025,.0025,.1025) ...
!---------
! for n>50 only endmember: 51:51, N:N
! for n=31-50 only one binary combination: 
! for n=26-30 only two binary combinations: 
! for n=2 and n=15-25 three binary cobinations: 
! for n=11-14 three binary and one ternary combination
! for n<=10 use full grid: 2 binar and 2 ternar combinarions
! NOTE for ybas=0.45 never add same endmember !!!
   double precision, dimension(5), parameter:: ybas=&
        [1.0D0,0.89D0,0.74D0,0.61D0,0.45D0]
   double precision, dimension(4), parameter :: ybin=&
        [0.11D0,0.26D0,0.39D0,0.15D0]
   double precision, dimension(3), parameter :: yter=[0.0D0,0.11D0,0.13D0]
! added: not here ... just for the dense grid
!   double precision, dimension(2), parameter :: yqrt=[0.35D0,0.19D0]
! for output of gridpoints
   integer jbas,sumngg,loksp,wrongngg
   logical trace,isendmem
   save sumngg,wrongngg
!
   write(*,*)'Illegal call to generate_grid: ',mode
   stop
!   if(mode.gt.0) write(*,*)'3Y entering generate_grid: ',mode,iph,ngg
!---------------------------------------------------------
! save current constitution in ydum
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1020
!---------------------------------------------------------
   if(test_phase_status_bit(iph,PHEXCB)) then
! This phase has charged endmembers, generate neutral gridpoints (also dense)
      call generate_charged_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
      goto 1000
   elseif(test_phase_status_bit(iph,PHIONLIQ)) then
! This is the ionic liquid, requires a special grid, also used for dense
!      if(mode.gt.0) write(*,*)'3Y gridpoint in the liquid',mode
      call generate_ionliq_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
!      if(mode.le.0) write(*,*)'3Y exit ionliq ',mode,ngg,ny
      if(mode.eq.-1) then
         wrongngg=ngg
      elseif(mode.eq.0) then
! ionliq -1 makes a bad estimate of the number of gridpoints generated 
! give a warning if it may be too wrong ...
         if(ngg-wrongngg.gt.1000) then
            write(*,*)'3Y warning: ionic liquid gridpoints: ',ngg,wrongngg
         endif
      endif
      goto 1000
   elseif(test_phase_status_bit(iph,PHFORD)) then
! this phase has 4 sublattice fcc/hcp tetrahedral ordering,
! this reduces the number of gridpoints UNFINISHED: NOT IMPLEMENTED YET
      call generate_fccord_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! do not jump to 1000 until the fccord routine implemented correctly
!      goto 1000
   elseif((btest(globaldata%status,GSXGRID) .or. & 
            test_phase_status_bit(iph,PHXGRID)) .and. &
        .not.test_phase_status_bit(iph,PHGAS)) then
! Generate extra gridpoints for all phases or a special phase but never for gas
      call generate_dense_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
      goto 1000
   endif
!
! mode=0 means generate grid (-1 means estimate size of grid for allocation), 
!     >0 means find constituent fractions for gridpoint in solution)
!
   if(mode.eq.0) then
!      write(*,*)'3Y Generating grid for phase: ',iph
! trace TRUE means generate outpt for each gridpoint
!      trace=.TRUE.
      trace=.FALSE.
      if(iph.eq.1 .and. trace) then
! unit 33 is opened before calling this routine
         sumngg=0
         write(33,43)
43       format('The constituent fractions, y, enclosed within parentheses',&
              'for each sublattice'/'Mole fractions after x:, Gibbs energies',&
              ' after G:'/)
      endif
      if(trace) then
         call get_phase_record(iph,lokph)
         write(33,44)iph,phlista(lokph)%name
44       format('Endmembers (EM) and gridpoints (GP) for phase: ',i3,1x,a)
      endif
   else
      trace=.FALSE.
   endif
!---------------------------------------------------------
! calculate the number of endmembers and index of first constituent in subl ll
   nend=1
   inkl(0)=0
   do ll=1,nsl
      nend=nend*nkl(ll)
      inkl(ll)=inkl(ll-1)+nkl(ll)
   enddo
!   iliqneut=0
! ionic liquids with neutrals ....
!   if(test_phase_status_bit(iph,PHIONLIQ)) then
!      loksp=0
!      do ny=nkl(1)+1,inkl(2)
!         loksp=knr(ny)
!         write(*,63)'3Y species: ',ny,knr(ny),loksp,&
!              splista(loksp)%charge,splista(loksp)%symbol
63       format(a,3i4,F10.5,2x,a)
!         if(.not.btest(splista(loksp)%status,SPVA) .and. &
!              abs(splista(loksp)%charge).eq.zero) then
! we have a neutral (vacancies has no mass), add an endmember for that
!            iliqneut=iliqneut+1
!            write(*,*)'3Y check for neutral: ',ny,iliqneut
!         endif
!      enddo
!   endif
!
   ny=inkl(nsl)
!   write(*,1010)'3Y Saved   ',iph,(ydum(i),i=1,ny)
!
! mode<0 means calculate size of arrays to allocate
!
   negmode: if(mode.lt.0) then
!---------------------------------------------------------
! just determine the number of gridpoints for this phase for global minimimum
! ideal gases should just have the endmembers ....
      ngg=nend
      lokph=phases(iph)
      if(nend.eq.1 .or. nend.gt.50 .or. &
           btest(phlista(lokph)%status1,PHID)) then
! >50 or 1 endmember or ideal phase: only endmembers
         ngg=nend
      elseif(nend.gt.30) then
! 31-50: only one binary combination
         ngg=nend*nend
      elseif(nend.gt.25) then
! 26-30: two binary combinations
         ngg=nend*(1+2*(nend-1))
      elseif(nend.eq.2 .or. nend.ge.15) then
! 2 or 15-25: three binary combinarions
         ngg=nend*(1+3*(nend-1))
      elseif(nend.gt.10) then
! 11-14: three binary and one ternary combinarion 
!         ngg=nend*(1+(nend-1)*(3+nend-2)) ! (ternary combination skipped)
!         ny=ngg
         ngg=nend*(1+3*(nend-1))
! added yqrt
!         je=nend*(nend-1)*(nend-2)
!         write(*,*)'3Y endmemers, ngg and je: ',nend,ngg,je
      else
! 3-10: three binary and two ternary combinarions (all)
! and the added (quaternary) combinatin
!         ngg=nend*(1+(nend-1)*(3+2*(nend-2))) ! (ternary combinations skipped)
!         ny=ngg
         ngg=nend*(1+3*(nend-1))
! added ygrt
         je=nend*(nend-1)*(nend-2)
!         write(*,*)'3Y endmemers, ngg and je: ',nend,ngg,je
      endif
!      write(*,*)'3Y endmembers and gridpoints: ',nend,ngg
!      read(*,11)ch1
11    format(a)
!      ngg=ngg+iliqneut
!      ngg=ngg
      if(ocv()) write(*,*)'3Y Generate grid: ',nend,ngg
      ny=nend
      goto 1001
   endif negmode
!------------------------------------------------------------
! for mode=0:
!    set gridpoint sitefractions and calculate G
! for mode>0:
!   return sitefractions (for mode=gridpoint number (part of the solution))
!   BUT: The only way to find the site fraction of a gripoint is to generate
!   all gridpoints up the one specified by the value of mode (no G calculation)
!   write(*,*)'3Y ggy: ',mode,iph,nsl,nend,inkl(nsl)
!   if(mode.gt.0) then
!      write(*,*)'3Y looking for allocate error: ',nsl,nend,inkl(nsl)
!   endif
   allocate(endm(nsl,nend))
   allocate(yfra(inkl(nsl)))
   nofy=inkl(nsl)
! generate endmembers, endm(ll,ie) is set to consituent index in sublattice ll
   je=1
   do ll=1,nsl
      endm(ll,je)=inkl(ll-1)+1
   enddo
100 continue
   je=je+1
   if(je.gt.nend) goto 120
   do ls=1,nsl
      endm(ls,je)=endm(ls,je-1)
   enddo
   ll=0
110 ll=ll+1
   if(endm(ll,je).lt.inkl(ll)) then
      endm(ll,je)=endm(ll,je)+1
   elseif(ll.lt.nsl) then
      endm(ll,je)=inkl(ll-1)+1
      goto 110
   else
      gx%bmperr=4148; goto 1000
   endif
   goto 100
120 continue
!   if(trace) then
!      do i=1,nend
!         write(33,125)i,(endm(ls,i),ls=1,nsl)
!125      format('endmem: ',i4,2x,10i3)
!      enddo
!   endif
150 continue
!---------------------------------------
! now generate all combinations of endmembers
!   write(*,*)'3Y endmembers and gridpoints: ',nend,ngg
!   read(*,11)ch1
   ngg=0
   lokph=phases(iph)
   endmem: do iend=1,nend
      yfra=yzero
      do ls=1,nsl
         yfra(endm(ls,iend))=ybas(1)
      enddo
!      write(*,180)'3Y yfra: ',1,iend,nkl(1),endm(1,iend),yfra(endm(1,iend))
180   format(a,4i3,6(1pe16.7))
      isendmem=.TRUE.
! initiate the loop veriables below for endmembers and fractions
      ibas=2
      ibin=1
      iter=1
      jend=0
      kend=0
200   continue
      ngg=ngg+1
      if(mode.gt.0) then
! if ngg=mode we have found the gridpoint! store y and x and quit
         if(ngg.eq.mode) goto 500
      else
! calculate G and composition and save
!         write(*,201)ibas,ngg,(yfra(is),is=1,inkl(nsl))
201      format('3Y ggz: ',i2,i4,5(F10.6))
         if(ocv()) write(*,*)'3Y Calculating gridpoint: ',ngg
         call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
         if(gx%bmperr.ne.0) goto 1000
         if(garr(ngg).gt.gmax) gmax=garr(ngg)
!         if(ngg.eq.15) then
!            write(*,520)'3Y cgx: ',(xarr(is,ngg),is=1,nrel)
!         endif
         if(trace) then
            if(isendmem) then
               write(33,153,advance='no')sumngg+ngg
153            format('EM:',i4,' y: ')
            else
               write(33,154,advance='no')sumngg+ngg
154            format('GP:',i4,' y: ')
            endif
            jbas=0
            do ls=1,nsl
               write(33,155,advance='no')(yfra(jbas+is),is=1,nkl(ls)-1)
155            format('(',10(F4.2,','))
               write(33,156,advance='no')yfra(jbas+nkl(ls))
156            format(F4.2,')')
               jbas=jbas+nkl(ls)
            enddo
            write(33,157,advance='no')(xarr(is,ngg),is=1,nrel)
157         format(' x:',8(f8.5))
            write(33,158)garr(ngg)
158         format(' G:',1pe12.4)
         endif
         isendmem=.FALSE.
      endif
! depending on nend value or ideal generate combinations
      if(nend.eq.1 .or. nend.gt.50 .or. &
           btest(phlista(lokph)%status1,PHID)) cycle
      yfra=yzero
      combend: if(nend.gt.30) then
! if nend=31..50, one binary combination, 961-2500
!    0.89*Y_E + 0.11*Y_F,F=/=E
         jend=jend+1
         if(jend.eq.iend) jend=jend+1
         if(jend.gt.nend) cycle
         do ls=1,nsl
            yfra(endm(ls,iend))=ybas(ibas)
            yfra(endm(ls,jend))=yfra(endm(ls,jend))+ybin(ibin)
         enddo
!      write(*,180)'3Y yfra: ',ibas,iend,nkl(1),endm(1,jend),yfra(endm(1,jend))
         goto 200
      elseif(nend.gt.25) then
! nend=26..30 two binary combinations, 1326-1770
!    0.89*Y_E + 0.11*Y_F,F=/=E
!    0.74*Y_E + 0.26*Y_F,F=/=E
         jend=jend+1
         if(jend.eq.iend) jend=jend+1
         if(jend.gt.nend) then
            if(ibas.eq.3) cycle
            jend=1
            ibas=3; ibin=2
         endif
         do ls=1,nsl
            yfra(endm(ls,iend))=ybas(ibas)
            yfra(endm(ls,jend))=yfra(endm(ls,jend))+ybin(ibin)
         enddo
!      write(*,180)'3Y yfra: ',ibas,jend,nkl(1),endm(1,jend),yfra(endm(1,jend))
         goto 200
      elseif(nend.eq.2 .or. nend.ge.15) then
! nend=2 or nend=15..25, three binary combinations, ??-1825
!    0.89*Y_E + 0.11*Y_F,F=/=E
!    0.74*Y_E + 0.25*Y_F,F=/=E
!    0.61*Y_E + 0.39*Y_F,F=/=E
         jend=jend+1
         if(jend.eq.iend) jend=jend+1
         if(jend.gt.nend) then
            if(ibas.eq.4) cycle
            ibas=ibas+1; ibin=ibin+1
            jend=1
            if(jend.eq.iend) jend=jend+1
         endif
         do ls=1,nsl
            yfra(endm(ls,iend))=ybas(ibas)
            yfra(endm(ls,jend))=yfra(endm(ls,jend))+ybin(ibin)
         enddo
!      write(*,180)'3Y yfra: ',ibas,jend,nkl(1),endm(1,jend),yfra(endm(1,jend))
         goto 200
      elseif(nend.gt.10) then
! complicated here, iterating in both binary and ternary combinations ....
! nend=11..14, 3 binary and one ternary combination, 1331-2744
!    0.89*Y_E + 0.11*Y_F,F=/=E
!    0.74*Y_E + 0.26*Y_F,F=/=E
!             + 0.15*Y_F + 0.11*Y_G,G=/=(E,F)
!    0.61*Y_E + 0.39*Y_F,F=/=E
         if(iter.eq.2) then
! we are interating in the ternary endmember
            stop 'no ternary for 10<nend<15'
!            write(*,*)'3Y Ternary combinations for 10<nend<15'
            kend=kend+1
            if(kend.eq.iend) kend=kend+1
            if(kend.eq.jend) kend=kend+1
            if(kend.gt.nend) then
               kend=1
               jend=jend+1
               if(jend.eq.iend) jend=jend+1
               if(jend.gt.nend) then
! all ternary combinations done .... ???
                  jend=1
                  ibas=4; ibin=3
                  iter=1
               endif
            endif
         else
! we are iterating in the binary endmembers
            jend=jend+1
            if(jend.eq.iend) jend=jend+1
            if(jend.gt.nend) then
               if(ibas.eq.4) cycle
!               if(ibas.eq.2) then
!                  if(iter.eq.1) then
!                     iter=2
!                     ibin=3
!                  endif
!               endif
               ibas=ibas+1; ibin=ibin+1
               jend=1
            endif
         endif
         do ls=1,nsl
            yfra(endm(ls,iend))=ybas(ibas)
            yfra(endm(ls,jend))=yfra(endm(ls,jend))+ybin(ibin)
            if(kend.gt.0) then
               yfra(endm(ls,kend))=yfra(endm(ls,kend))+yter(iter)
            endif
         enddo
!      write(*,180)'3Y yfra: ',ibas,jend,nkl(1),endm(1,jend),yfra(endm(1,jend))
         goto 200
      else
! nend=3..10, 3 binary and 2 ternary combinations, 33-1720
!    0.89*Y_E + 0.11*Y_F,F=/=E
!    0.74*Y_E + 0.26*Y_F,F=/=E
!             + 0.15*Y_F + 0.11*Y_G,G=/=(E,F)
!    0.61*Y_E + 0.39*Y_F,F=/=E
!             + 0.25*Y_F + 0.14*Y_G,G=/=(E,F)
         if(iter.eq.2 .or. iter.eq.3) then
! we are interating in the ternary endmember
            stop 'no ternary for nend<10'
!            write(*,*)'3Y Ternary combinations for 2<nend<10'
            kend=kend+1
            if(kend.eq.iend) kend=kend+1
            if(kend.eq.jend) kend=kend+1
            if(kend.gt.nend) then
               kend=1
               jend=jend+1
               if(jend.eq.iend) jend=jend+1
               if(jend.gt.nend) then
! all second ternary combinations done .... then finished !!!
                  if(iter.eq.3) cycle
                  jend=1
                  ibas=4; ibin=3
                  iter=1
                  abrakadabra=1
               endif
            endif
         else
! we are iterating in the binary endmembers
460         continue
            jend=jend+1
            if(jend.eq.iend) jend=jend+1
            if(jend.gt.nend) then
!               if(ibas.eq.2) then
!                  abrakadabra=0
!                  iter=2; ibin=2
!                  kend=1
!               else
!                  abrakadabra=1
!                  iter=3; ibin=3; ibas=3
!                  kend=1
!               endif
               if(ibas.eq.4) cycle
               ibas=ibas+1; ibin=ibin+1
               jend=0
               goto 460
            endif
         endif
         do ls=1,nsl
! attempt to generate better start values for fcc-protototype ordering
            yfra(endm(ls,iend))=ybas(ibas)
            yfra(endm(ls,jend))=yfra(endm(ls,jend))+ybin(ibin)
            if(kend.gt.0) then
               yfra(endm(ls,kend))=yfra(endm(ls,kend))+yter(iter)
            endif
         enddo
!      write(*,180)'3Y yfra: ',ibas,kend,nkl(1),endm(1,kend),yfra(endm(1,kend))
         goto 200
      endif combend
   enddo endmem
! finished all calculations
   ny=nend
   if(trace) sumngg=sumngg+ngg
   goto 1000
! all binary and ternary combination of endmember for the grid above
!
! This should be modeified to take into account option F and B as those
! sublattices are identical ... and charged constituents ... have fun ...
!
!----------------------------------------
! here we return the constitution for gridpoint "mode" in the solution
! We must also return the mole fractions ... NO??
500 continue
!   write(*,510)'3Y ggy: ',mode,ngdim,ny,(yfra(i),i=1,ny)
510 format(a,3i5,10(F6.3))
! return values of yfra in yarr.  Note xarr calculated above also returned
    do i=1,ny
       yarr(i)=yfra(i)
    enddo
!    write(*,520)'3Y ggx: ',mode,(xarr(i,mode),i=1,nrel)
520 format(a,i4,10(f8.5))
!----------------------------------------------------------
1000 continue
   if(allocated(endm)) then
      deallocate(endm)
      deallocate(yfra)
   endif
1001 continue
!------------------------------------------------------------
! IMPORTANT !! restore original constituent fractions
!   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   errsave=gx%bmperr
   gx%bmperr=0
!   write(*,1010)'3Y Restore ',iph,(ydum(i),i=1,ny)
1010 format(a,'const for ',i3,10(f6.3))
   call set_constitution(iph,1,ydum,qq,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'3Y Error restoring constitution for phase: ',iph,gx%bmperr
   endif
   gx%bmperr=errsave
! jump here if error saving constitution!!
1020 continue
   if(gx%bmperr.ne.0) write(*,*)'3Y gengrid error: ',gx%bmperr
   return
 end subroutine generate_grid

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine generic_grid_generator(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! This generates grid for any phase
! mode=0 generate grid for phase iph and mole fraction and G for all points
!    ngg on entry is max number of gridpoints, on exit number of gridpoints
!    nrel is number of elements
!    xarr(1..nrel,gp) is composition of gripoint gp, garr(iv) its G
!    ny,yarr,gmax not used
! mode>0 return constitution for gridpoint number mode in yarr
!    iph is returned as phase index for gripoint mode
!    xarr(1..,nrel) the composition at the gripoint, garr not nused
!    ny is number of constituent fractions
!    yarr are the constituent fractions
!    gmax not used ??
!
   implicit none
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*),gmax
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
! local loop variables etc
   integer ii,ij,ik,il,im,in,is,ie,incl(0:maxsubl),maxng,ng,ncon
   integer nend,nendj,nendk,nendl,nendm
   integer ijs,iks,ils,ims
! these are for call of get_phase_data
   integer nsl,nkl(maxsubl),knr(maxconst)
   double precision ydum(maxconst),sites(maxsubl),qq(5)
! this is for generating endmembers
   integer, dimension(:,:), allocatable :: endm
! these are for generating the constituion of gridpoint
   double precision, dimension(:,:), allocatable :: yendm
   double precision, dimension(:), allocatable :: yfra
! these are factors to generate gridpoints
   double precision, dimension(5), parameter :: &
!        yf=[0.33D0,0.28D0,0.18D0,0.08D0,0.03D0]
!        yf=[0.33D0,0.28D0,0.18D0,0.08D0,0.03D0] to test with map3
        yf=[0.33D0,0.28D0,0.18D0,0.14D0,0.11D0]
!        yf=[0.33D0,0.28D0,0.18D0,0.14D0,0.11D0] OK for fuel
!        yf=[0.11D0,0.13D0,0.18D0,0.23D0,0.35D0]
   logical gas,dense,verydense,gles
!---------------
! handle special phases like ionic crystals, ionic liquids and ordere/disorder
!   write(*,*)'3Y in generic_grid_generator'
   gas=.FALSE.
   if(test_phase_status_bit(iph,PHEXCB)) then
! crystalline phase with charged endmembers
      call generate_charged_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
!      write(*,*)'3Y charged grid: ',ngg
      goto 1000
   elseif(test_phase_status_bit(iph,PHIONLIQ)) then
! This is the ionic liquid, requires a special grid, also used for dense
      call generate_ionliq_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
      goto 1000
!   elseif(test_phase_status_bit(iph,PHFORD)) then
! this phase has 4 sublattice fcc/hcp tetrahedral ordering,
! this reduces the number of gridpoints UNFINISHED: NOT IMPLEMENTED YET
!      call generate_fccord_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! do not jump to 1000 until the fccord routine implemented correctly
!      goto 1000
   elseif(test_phase_status_bit(iph,PHGAS)) then
      gas=.TRUE.
   elseif((btest(globaldata%status,GSXGRID) .or. & 
            test_phase_status_bit(iph,PHXGRID)) .and. &
        .not.test_phase_status_bit(iph,PHGAS)) then
! Generate extra gridpoints for all phases or a special phase but never for gas
!      call generate_dense_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
      dense=.TRUE.
   else
      dense=.FALSE.
   endif
!      goto 1000
!----------------
! get phase model
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! max number of gridpoints allowed, ngg returned as number of gridpoints...
   maxng=ngg
   ngg=0
! incl(ii) set to number of constituents up to including sublattice ii
   incl(0)=0
   nend=1
   do ii=1,nsl
      nend=nend*nkl(ii)
      incl(ii)=incl(ii-1)+nkl(ii)
   enddo
   ncon=incl(nsl)
! nend is number of endmembers, endm(1..nsl,ii) are constituent index of ii
! yendm(1..nsl,ii) has the constituent fractions for endmember ii
! yfra is used to generate a constitutuon from a combination of endmembers
   allocate(endm(nsl,nend))
   allocate(yendm(ncon,nend))
   allocate(yfra(ncon))
   yendm=1.0D-12
! set endm(1..nsl,1) to first constituent index for each sublattice
   do ij=1,nsl
      endm(ij,1)=incl(ij-1)+1
      yendm(endm(ij,1),1)=one
   enddo
! loop to increment the constituents to generate all endmembers
   newend: do ii=2,nend
      ij=1
! copy previous endmember
      do ij=1,nsl
         endm(ij,ii)=endm(ij,ii-1)
      enddo
! increment one constituent starting from first sublattice
      do ij=1,nsl
         if(endm(ij,ii).lt.incl(ij)) then
            endm(ij,ii)=endm(ij,ii)+1
            do ik=1,nsl
               yendm(endm(ik,ii),ii)=one
            enddo
            cycle newend
         else
            endm(ij,ii)=endm(ij,1)
         endif
      enddo
   enddo newend
! output seems OK here
!   write(*,20)'3Y gend1: ',iph,nsl,(nkl(ii),ii=1,nsl)
20 format(a,2i3,2x,10i3)
!   write(*,20)'3Y gend2: ',ncon,nend,(incl(ii),ii=0,nsl)
!   do ii=1,nsl
!      write(*,21)(endm(ii,ij),ij=1,nend)
21    format(26i3)
!   enddo
!   do ii=1,nend
!      write(*,22)ii,(yendm(ij,ii),ij=1,incl(nsl))
22    format(i3,20F4.1)
!   enddo
! now generate an grid depending on nend mixing up to 5 different endmembers.
! up to for 4 endmembers 4*4*4*4*4=1024
! 5 to 7 endmembers      7*7*7*7 =2401
! 7 to 13 endmembers     13*13*13=2197
! max 50 endmembers      50*50 = 2500
! for N>50 endmembers          = N
! FOR DENSE about 10 times more
! up to 7 endmembers 7*7*7*7*7 = 16807
! 8 to 12 endmembers 12*12*12*12 = 20736
! 13 to 15 endmembers 15*15*15  =33750
! max 150 endmembers 150*150 = 22500
! for N>150                  = N
!--------------------------------------
   ng=0
   if(gas) then
      do ii=1,nend
         do is=1,ncon
            yfra(is)=yendm(is,ii)
         enddo
         ng=ng+1
         if(mode.eq.0) then
            if(ng.gt.maxng) then
               write(*,*)'3Y Too many gridpoints 7',ng,maxng,iph
               gx%bmperr=4399; goto 1000
            endif
            call calc_gridpoint(iph,yfra,nrel,xarr(1,ng),garr(ng),ceq)
            if(gx%bmperr.ne.0) then
               write(*,*)'3Y error calculating gridpoint: ',&
                    iph,gx%bmperr
               goto 1000
            endif
         elseif(mode.eq.ng) then
! when mode>0 we just want to know the constituent fractions
            goto 900
         endif
      enddo
      goto 800
   endif
!-----------------------
!   write(*,*)'3Y grid1: ',nend,dense
   gles=.not.dense
   iiloop: do ii=1,nend
      ijs=1
      nendj=nend
      if(nend.gt.150 .or. (gles .and. nend.gt.50)) then
         nendj=ii
         ijs=ii
      endif
      ijloop: do ij=ijs,nendj
         iks=1
         nendk=nend
         if(nend.gt.15 .or. (gles .and. nend.gt.13)) then
            nendk=ij
            iks=ij
         endif
         ikloop: do ik=1,nendk
            ils=1
            nendl=nend
            if(nend.gt.12 .or. (gles .and. nend.gt.7)) then
               nendl=ik
               ils=ik
            endif
            illoop: do il=ils,nendl
               ims=1
               nendm=nend
               if(nend.gt.7 .or. (gles .and. nend.gt.4)) then
! with 4 endmembers 1024 gridpoints
                  nendm=il
                  ims=il
               endif
               imloop: do im=ims,nendm
! sum up the weighted fractions from the different endmembers
                  do is=1,ncon
                     yfra(is)=yf(1)*yendm(is,ii)+yf(2)*yendm(is,ij)+&
                          yf(3)*yendm(is,ik)+yf(4)*yendm(is,il)+&
                          yf(5)*yendm(is,im)
                  enddo
                  ng=ng+1
!                  write(*,23)'3Y imloop1: ',ng,ii,ij,ik,il,im,0.0D0,yfra
23                format(a,i5,5i3,': ',1pe12.4,0p10F5.2)
                  if(mode.eq.0) then
! strange bug in map3, maxng was zero sometimes ...
                     if(ng.gt.maxng) then
                        if(maxng.lt.100) then
                           write(*,*)'3Y max gripoints wrong 6: ',maxng,iph,mode
                        else
                           write(*,*)'3Y Too many gridpoints 6',ng,maxng,iph
                           gx%bmperr=4399; goto 1000
                        endif
                     endif
                     call calc_gridpoint(iph,yfra,nrel,xarr(1,ng),garr(ng),ceq)
                     if(gx%bmperr.ne.0) then
                        write(*,*)'3Y error calculating gridpoint: ',&
                             iph,gx%bmperr
                        goto 1000
                     endif
!                     write(*,23)'3Y imloop1: ',ng,ii,ij,ik,il,im,garr(ng),yfra
                  elseif(mode.eq.ng) then
! when mode>0 we just want to know the constituent fractions
                     goto 900
                  endif
               enddo imloop
            enddo illoop
         enddo ikloop
      enddo ijloop
   enddo iiloop
800 continue
   if(mode.gt.0) then
      write(*,*)'3Y could not find gridpoint ',mode,' in phase, ',iph
      gx%bmperr=4399
   else
      ngg=ng
   endif
   goto 1000
!--------------------------------------------
! we found the gridpoint we were looking for
900 continue
   ny=ncon
   do ii=1,ny
      yarr(ii)=yfra(ii)
   enddo
1000 continue
!   write(*,*)'3Y finished generic mode ',mode,ngg
   return
 end subroutine generic_grid_generator

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine generate_dense_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! generates more gridpoints than default generate_grid
! Different action depending of the value of mode, 
! for mode<0:  
!    return the number of gridpoints that will be generated for phase iph in ngg
! for mode=0:
!    return garr(i) gibbs energy and xarr(1,i) the compositions of gridpoint i
! for mode>0:
!    return site fractions of gridpoint mode in yarr, number of fractions in ny
!    iph is phase number, ngg is number of gridpoints, nrel number of elements,
! if mode=0:
!    return xarr mole fractions of gridpoints, garr Gibbs energy of gridpoints,
!    ngg is dimension of garr
! if mode>0:
!   "mode" is a gridpoint of this phase in solution, return number of 
!   constituent fractions in ny and fractions in yarr for this gridpoint
! The current constitution is restored at the end of the subroutine
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*),gmax
!\end{verbatim} %+
!
   integer lokph,errsave
   double precision, parameter :: yzero=1.0D-12
   integer abrakadabra,i,ibas,ibin,iend,is,iter,je,jend,kend,ll,ls,nend
   double precision ydum(maxconst)
   integer ngdim,nsl
   integer nkl(maxsubl),knr(maxconst),inkl(0:maxsubl),nofy
   double precision, dimension(:), allocatable :: yfra
   double precision sites(maxsubl),qq(5)
   real, allocatable :: xbrr(:)
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
   integer, dimension(:,:), allocatable :: endm
!--------------------------------
! grid is generated by combining end endmembers
! Number of endmemers is N
! For endmember E=1..N set fraction of enmember 
!    0.99*Y_E + 0.01*Y_all                             N*N of these
!    0.95*Y_E + 0.03*Y_F + 0.01*Y_all, F.ne.E          N*N*(N-1)
!    0.91*Y_E + 0.07*Y_F + 0.02*Y_all, F.ne.E          N*N*(N-1)
!    0.80*Y_E + 0.15*Y_F + 0.05*Y_all, F.ne.E          N*N*(N-1)
!    0.68*Y_E + 0.25*Y_F + 0.07*Y_all, F.ne.E          N*N*(N-1)
! or 0.68*Y_E + 0.16*Y_F + 0.16*Y_all, F.ne.E          N*N*(N-1)
!    0.54*Y_E + 0.36*Y_F + 0.10*Y_all, F.ne.E          N*N*(N-1)
!    0.42*Y_E + 0.35*Y_F + 0.23*Y_G, F.ne.E.ne.G       N*(N-1)*(N-2)
! or 0.42*Y_E + 0.40*Y_F + 0.18*Y_G, F.ne.E.ne.G       N*(N-1)*(N-2)
! or 0.48*Y_E + 0.40*Y_F + 0.12*Y_G, F.ne.E.ne.G       N*(N-1)*(N-2)
! with 2 endmembers: 24 gridpoints
! (1.00,0.00) (0.99,0.01) (0.01,0.99) (0.00,1.00)
! (0.96,0.04) (0.04,0.96) *2
! (0.93,0.07) (0.91,0.09) *2
! (0.85,0.15) (0.80,0.20) *2
! (0.75,0.25) (0.68,0.32) *2
! (0.64,0.36) (0.57,0.43) *2
! with 3 endmembers: 9+5*9*2+6=15+90=105
! (1.00,0,0) (0.99,0.01,0) (0.99,0,0.01)    *3
! (0.97,0.03,0) (0.96,0.04,0) (0.96,0.03,0.01) *3   <50
! (0.92,0.08,0) (0.90,0.10,0) (0.90,0.08,0.02)      <25
! (0.75,0.25,0) (0.68,0.32.0) (0.68,0.25,0.07)      <20
! (0.85,0.15,0) (0.80,0.20,0) (0.80,0,15,0.05)      <15
! (0.64,0.36,0) (0.54,0.46,0) (0.54,0.36,0.10)      <11
! (0.42,0.35,0.23) (0.42,0.23,0.35) (0.35,0.42,0.23) ... 6
!----------
! M=   N*N  + 5*N*N*(N-1) + N*(N-1)*(N-2)
! N=10  100 + 5*100*9     +   10*9*8      = 4600+720 = 5320
! N=20  400 + 400*19      +   0           >8000
! ...
   integer, parameter :: breaks6=50
   integer, parameter, dimension(5) :: breaks=[9,12,15,18,21]
!   integer, parameter, dimension(5) :: breaks=[9,12,15,18,55]
   double precision, dimension(-1:6), parameter:: ybas=&
        [1.00D0,0.99D0,0.96D0,0.91D0,0.68D0,0.80D0,0.54D0,0.44D0]
   double precision, dimension(6), parameter :: ybin=&
                      [0.03D0,0.07D0,0.16D0,0.15D0,0.36D0,0.44D0]
!                      [0.03D0,0.07D0,0.25D0,0.15D0,0.36D0,0.35D0]
   double precision, dimension(6), parameter :: yter=&
                      [0.01D0,0.02D0,0.16D0,0.05D0,0.10D0,0.12D0]
!                      [0.01D0,0.02D0,0.07D0,0.05D0,0.10D0,0.23D0]
! for output of gridpoints
   integer jbas,sumngg,loksp,l0,l1,ncon,jj,anion,isp
   logical trace,isendmem
   double precision ysum
   save sumngg
   character ch1*1
!
!   write(*,17)mode,iph,ngg
17 format('3Y entering generate_dense_grid: ',i2,i3,i10)
   if(mode.eq.0) then
!      write(*,*)'3Y Generating grid for phase: ',iph
! trace TRUE means generate outpt for each gridpoint
!      trace=.TRUE.
      trace=.FALSE.
      if(iph.eq.1 .and. trace) then
! unit 33 is opened before calling this routine
         sumngg=0
         write(33,43)
43       format('The constituent fractions, y, enclosed within parentheses',&
              'for each sublattice'/'Mole fractions after x:, Gibbs energies',&
              ' after G:'/)
      endif
      if(trace) then
         call get_phase_record(iph,nend)
!         write(33,44)iph,phlista(nend)%name
44       format('Endmembers (EM) and gridpoints (GP) for phase: ',i3,1x,a)
      endif
   else
      trace=.FALSE.
   endif
!   write(*,*)'3Y Getting phase data',mode
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! calculate the number of endmembers and index of first constituent in subl ll
   nend=1
   inkl(0)=0
   lokph=phases(iph)
   do ll=1,nsl
      if(btest(phlista(lokph)%status1,PHIONLIQ) .and. ll.eq.2) then
! multiply with charged anions and Va only, add neutrals
         do jj=1,nkl(2)
! knr(i) is species(i) location but I use constitlist as I have access to it
            isp=phlista(lokph)%constitlist(nkl(1)+jj)
            if(btest(splista(isp)%status,SPION) .or. &
                 btest(splista(isp)%status,SPVA)) then
               anion=anion+1
               cycle
            endif
         enddo
         nend=nend*anion+phlista(lokph)%nooffr(2)-anion
      else
! this is the "normal" number of endmembers
         nend=nend*nkl(ll)
      endif
      inkl(ll)=inkl(ll-1)+nkl(ll)
   enddo
!   if(btest(phlista(lokph)%status1,PHIONLIQ)) then
!      write(*,*)'3Y ionic liq: ',anion,nend
!   endif
   ny=inkl(nsl)
   ncon=inkl(nsl)
   negmode: if(mode.lt.0) then
!---------------------------------------------------------
! just determine the number of gridpoints for this phase for global minimimum
! ideal gases should just have the endmembers ....
! Hm, gases with ions??
      ngdim=ngg
      ngg=nend
      lokph=phases(iph)
!      write(*,*)'3Y nend 1: ',mode,ngg
      if(nend.eq.1 .or. nend.gt.100 .or. &
           btest(phlista(lokph)%status1,PHID)) then
! >100 or 1 endmember or ideal phase: only endmembers
         ngg=nend
      else
         ngg=nend
! The limits for various combinations will be adjusted when testing ...
! Max about 20000 gridpoints per phase
!         if(nend.ge.50) then
!         if(nend.ge.100) then
         ngg=ngg+nend*(nend-1)
!         write(*,*)'3Y dense -1A: ',iph,nend,ngg,breaks(5)
! ATTENTION
! The calculation below is not correct, it overestimates a bit the number of 
! gridpoints actually generated but it should not matter so much ... I hope
! When matching a gridpoint in the solution the code to generate the
! gridpoint is used, the code below is just an estimate for allocation
!         if(nend.le.50) then
! Try 60 to handle 53 endmembers in liquid noc2500.TDB from TAF-ID
!         if(nend.le.60) then
         if(nend.le.breaks6) then
            if(nend.gt.breaks(5)) then
               ngg=ngg+nend*(nend-1)+nend*nend*(nend-1)
!               write(*,*)'3Y dense -1B: ',iph,nend,ngg,breaks(4)
            elseif(nend.gt.breaks(5)) then
               ngg=ngg+nend*(nend-1)+2*nend*nend*(nend-1)
            elseif(nend.gt.breaks(4)) then
               ngg=ngg+nend*(nend-1)+3*nend*nend*(nend-1)
            elseif(nend.gt.breaks(3)) then
               ngg=ngg+nend*(nend-1)+4*nend*nend*(nend-1)
            elseif(nend.gt.breaks(2)) then
               ngg=ngg+nend*(nend-1)+5*nend*nend*(nend-1)
!            elseif(nend.gt.breaks(1)) then
            else
               ngg=ngg+nend*(nend-1)+6*nend*nend*(nend-1)
            endif
         endif
!         write(*,*)'3Y dense -1X: ',iph,nend,ngg
      endif
!      write(*,*)'3Y endmembers and gridpoints: ',nend,ngg
!      read(*,11)ch1
11    format(a)
      ny=nend
      goto 1001
   endif negmode
!------------------------------------------------------------
! for mode=0:
!    set gridpoint sitefractions and calculate G
! for mode>0:
!   return sitefractions (for mode=gridpoint number (part of the solution))
!   BUT: The only way to find the site fraction of a gripoint is to generate
!   all gridpoints up the one specified by the value of mode (no G calculation)
!   write(*,*)'3Y ggy: ',mode,iph,nsl,nend,inkl(nsl)
!
   allocate(yfra(inkl(nsl)))
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
! endm(1,1) is constituent in sublattice 1 of first endmember
! endm(2,1) is constituent in sublattice 2 of first endmember
! endm(nsl,2) is constituent in sublattice nsl of second endmember
! endm(1..nsl,nend) are constituents in all sublattices of last endmember
   allocate(endm(nsl,nend))
! inkl(nsl) is the number of fraction variables in the phase
!   allocate(yfra(inkl(nsl)))
   allocate(xbrr(noofel))
!   nofy=inkl(nsl)
! generate endmembers, endm(ll,ie) is set to consituent index in sublattice ll
   je=1
   do ll=1,nsl
      endm(ll,je)=inkl(ll-1)+1
   enddo
100 continue
   je=je+1
! if je>nend we are finished ...
   if(je.gt.nend) goto 120
   do ls=1,nsl
      endm(ls,je)=endm(ls,je-1)
   enddo
   ll=0
110 ll=ll+1
   if(endm(ll,je).lt.inkl(ll)) then
      endm(ll,je)=endm(ll,je)+1
   elseif(ll.lt.nsl) then
      endm(ll,je)=inkl(ll-1)+1
      goto 110
   else
      gx%bmperr=4148; goto 1000
   endif
   goto 100
!---------------------------------------
! We have now generated endm(1..nsl,j)
120 continue
!   write(*,202)'3Y special 1: ',nsl,nend,inkl(nsl),endm(1,2),endm(2,2),&
!        endm(1,nend),endm(1,3),endm(2,3),endm(1,4),endm(2,4)
! now generate all unary, binary and ternary combinations of endmember fractions
! Note the sum of constituent fractions in all sublattices must be unity
! By combining endmember fractions weighted according to ybas, ybin and yter
! we can ensure that
   ngg=0
   l0=0
   l1=0
   lokph=phases(iph)
   endmem1: do iend=1,nend
! we start with a new endmember iend, ybas(1) is 1.00
      l0=ngg
      yfra=yzero
      do ls=1,nsl
         yfra(endm(ls,iend))=ybas(-1)
      enddo
      ngg=ngg+1
      if(mode.eq.0) then
! this is fór a single endmember
!         write(*,201)'3Y end: ',ngg,(yfra(is),is=1,inkl(nsl))
         call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
         if(gx%bmperr.ne.0) goto 1000
         if(garr(ngg).gt.gmax) gmax=garr(ngg)
201      format(a,i5,20(F5.2))
      elseif(ngg.eq.mode) then
         goto 500
      endif
! binary combinations 0.99*y1 + 0.01*y2
      endmem2a: do jend=1,nend
         if(jend.eq.iend) cycle endmem2a
         yfra=zero
!         write(*,202)'3Y special 3: ',endm(1,2)
         do ls=1,nsl
! to generate better start values for fcc-protototype ordering
!            write(*,202)'3Y ls iend endm: ',ls,iend,jend,endm(ls,iend)
202         format(a,10i6)
            yfra(endm(ls,iend))=ybas(0)
            yfra(endm(ls,jend))=yfra(endm(ls,jend))+yter(1)
         enddo
         ngg=ngg+1
         if(mode.eq.0) then
! this is for 0.99*y1 + 0.01*y2
! STRANGE error that destroyed endm after the call to calc_gridpoint!!
! the error was due to wrong size allocated to xarr which is strange as it
! is done elsewhere but the error disapperared when I allocated a larger
! xarr although the allocated one did not seem too small. 
            call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
!            call calc_gridpoint2(iph,yfra,endm,nrel,xarr(1,ngg),garr(ngg),ceq)
!            call calc_gridpoint2(iph,yfra,endm,nrel,xbrr,garr(ngg),ceq)
            if(gx%bmperr.ne.0) goto 1000
            if(garr(ngg).gt.gmax) gmax=garr(ngg)
!            write(*,201)'3Y bin: ',ngg,(yfra(is),is=1,inkl(nsl))
         elseif(ngg.eq.mode) then
            goto 500
         endif
      enddo endmem2a
! ternary combinations
! ibas   -1    0      1      2      3      4      5      6  
! ybas 1.00D0,0.99D0,0.96D0,0.91D0,0.68D0,0.80D0,0.54D0,0.42D0
! ybin               0.03D0,0.07D0,0.25D0,0.15D0,0.36D0,0.35D0
! yter               0.01D0,0.02D0,0.07D0,0.05D0,0.10D0,0.23D0
!      if(nend.gt.50) cycle endmem1
!      if(nend.gt.60) cycle endmem1
      if(nend.gt.breaks6) cycle endmem1
      ibasloop: do ibas=1,6
         if(nend.ge.breaks(5) .and. ibas.eq.2) cycle endmem1
         if(nend.ge.breaks(4) .and. ibas.eq.3) cycle endmem1
         if(nend.ge.breaks(3) .and. ibas.eq.4) cycle endmem1
         if(nend.ge.breaks(2) .and. ibas.eq.5) cycle endmem1
         if(nend.ge.breaks(1) .and. ibas.eq.6) cycle endmem1
         endmem2b: do jend=1,nend
            if(jend.eq.iend) cycle endmem2b
            endmem3: do kend=1,nend
               yfra=zero
               do ls=1,nsl
                  yfra(endm(ls,iend))=ybas(ibas)
                  yfra(endm(ls,jend))=yfra(endm(ls,jend))+ybin(ibas)
                  yfra(endm(ls,kend))=yfra(endm(ls,kend))+yter(ibas)
               enddo
               ysum=zero
               do ls=1,inkl(nsl)
                  ysum=ysum+yfra(ls)
               enddo
               ngg=ngg+1
               if(mode.eq.0) then
! this is for 0.96*y1 + 0.03*y2+0.01*y3
                  call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  if(garr(ngg).gt.gmax) gmax=garr(ngg)
!                  write(*,201)'3Y ter: ',ngg,(yfra(is),is=1,inkl(nsl)),ysum
               elseif(ngg.eq.mode) then
                  goto 500
               endif
            enddo endmem3
         enddo endmem2b
         l1=ngg
      enddo ibasloop
   enddo endmem1
!   write(*,*)'3Y Calculated points: ',ngg
   goto 1000
!----------------------------------------
! here we return the constitution for gridpoint "mode" in the solution
! We must also return the mole fractions ... NO??
500 continue
!    write(*,505)'3Y ggg: ',mode,iph,nsl,inkl(nsl),ny
505 format(a,i7,i4,2x,i2,i5,i10)
!   write(*,510)'3Y ggy: ',mode,ngdim,ny,(yfra(i),i=1,ny)
510 format(a,2i5,i3,10(F6.3))
   do i=1,ny
      yarr(i)=yfra(i)
   enddo
!   write(*,520)'3Y ggx: ',(xarr(is,ngg+ngdim),is=1,nrel)
520 format(a,10(f8.5))
1000 continue
! these will be deallocated by default when exit this subroutine ...
   if(allocated(endm)) then
      deallocate(endm)
      deallocate(yfra)
   endif
1001 continue
! restore original constituent fractions
!   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   errsave=gx%bmperr
   gx%bmperr=0
!   write(*,1010)'3Y Restore ',iph,(ydum(i),i=1,ny)
1010 format(a,'const for ',i3,10(f6.3))
   call set_constitution(iph,1,ydum,qq,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'3Y Error restoring constitution for phase: ',iph,gx%bmperr
   endif
   gx%bmperr=errsave
   if(gx%bmperr.ne.0) write(*,*)'3Y gengrid error: ',gx%bmperr
   return
 end subroutine generate_dense_grid

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine generate_dense_gridX(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! generates more gridpoints than default generate_grid
! Different action depending of the value of mode, 
! for mode<0:  
!    return the number of gridpoints that will be generated for phase iph in ngg
! for mode=0:
!    return garr(i) gibbs energy and xarr(1,i) the compositions of gridpoint i
! for mode>0:
!    return site fractions of gridpoint mode in yarr, number of fractions in ny
!    iph is phase number, ngg is number of gridpoints, nrel number of elements,
! if mode=0:
!    return xarr mole fractions of gridpoints, garr Gibbs energy of gridpoints,
!    ngg is dimension of garr
! if mode>0:
!   "mode" is a gridpoint of this phase in solution, return number of 
!   constituent fractions in ny and fractions in yarr for this gridpoint
! The current constitution is restored at the end of the subroutine
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*),gmax
!\end{verbatim} %+
!
   integer lokph,errsave
   double precision, parameter :: yzero=1.0D-12
   integer abrakadabra,i,ibas,ibin,iend,is,iter,je,jend,kend,ll,ls,nend
   double precision ydum(maxconst)
   integer ngdim,nsl
   integer nkl(maxsubl),knr(maxconst),inkl(0:maxsubl),nofy
   double precision, dimension(:), allocatable :: yfra
   double precision sites(maxsubl),qq(5)
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
   integer, dimension(:,:), allocatable :: endm
!--------------------------------
! grid is generated by combining end endmembers
! Number of endmemers is N
! For endmember E=1..N set fraction of enmember 
!    0.99*Y_E + 0.01*Y_all                             N*N of these
!    0.95*Y_E + 0.03*Y_F + 0.01*Y_all, F.ne.E          N*N*(N-1)
!    0.91*Y_E + 0.07*Y_F + 0.02*Y_all, F.ne.E          N*N*(N-1)
!    0.80*Y_E + 0.15*Y_F + 0.05*Y_all, F.ne.E          N*N*(N-1)
!    0.68*Y_E + 0.25*Y_F + 0.07*Y_all, F.ne.E          N*N*(N-1)
! or 0.68*Y_E + 0.16*Y_F + 0.16*Y_all, F.ne.E          N*N*(N-1)
!    0.54*Y_E + 0.36*Y_F + 0.10*Y_all, F.ne.E          N*N*(N-1)
!    0.42*Y_E + 0.35*Y_F + 0.23*Y_G, F.ne.E.ne.G       N*(N-1)*(N-2)
! or 0.42*Y_E + 0.40*Y_F + 0.18*Y_G, F.ne.E.ne.G       N*(N-1)*(N-2)
! or 0.48*Y_E + 0.40*Y_F + 0.12*Y_G, F.ne.E.ne.G       N*(N-1)*(N-2)
! with 2 endmembers: 24 gridpoints
! (1.00,0.00) (0.99,0.01) (0.01,0.99) (0.00,1.00)
! (0.96,0.04) (0.04,0.96) *2
! (0.93,0.07) (0.91,0.09) *2
! (0.85,0.15) (0.80,0.20) *2
! (0.75,0.25) (0.68,0.32) *2
! (0.64,0.36) (0.57,0.43) *2
! with 3 endmembers: 9+5*9*2+6=15+90=105
! (1.00,0,0) (0.99,0.01,0) (0.99,0,0.01)    *3
! (0.97,0.03,0) (0.96,0.04,0) (0.96,0.03,0.01) *3   <50
! (0.92,0.08,0) (0.90,0.10,0) (0.90,0.08,0.02)      <25
! (0.75,0.25,0) (0.68,0.32.0) (0.68,0.25,0.07)      <20
! (0.85,0.15,0) (0.80,0.20,0) (0.80,0,15,0.05)      <15
! (0.64,0.36,0) (0.54,0.46,0) (0.54,0.36,0.10)      <11
! (0.42,0.35,0.23) (0.42,0.23,0.35) (0.35,0.42,0.23) ... 6
!----------
! M=   N*N  + 5*N*N*(N-1) + N*(N-1)*(N-2)
! N=10  100 + 5*100*9     +   10*9*8      = 4600+720 = 5320
! N=20  400 + 400*19      +   0           >8000
! ...
   integer, parameter :: breaks6=50
   integer, parameter, dimension(5) :: breaks=[9,12,15,18,21]
   double precision, dimension(-1:6), parameter:: ybas=&
        [1.00D0,0.99D0,0.54D0,0.91D0,0.68D0,0.80D0,0.96D0,0.44D0]
   double precision, dimension(6), parameter :: ybin=&
                      [0.03D0,0.07D0,0.16D0,0.15D0,0.36D0,0.44D0]
!                      [0.03D0,0.07D0,0.25D0,0.15D0,0.36D0,0.35D0]
   double precision, dimension(6), parameter :: yter=&
                      [0.01D0,0.02D0,0.16D0,0.05D0,0.10D0,0.12D0]
!                      [0.01D0,0.02D0,0.07D0,0.05D0,0.10D0,0.23D0]
! for output of gridpoints
   integer jbas,sumngg,loksp,l0,l1,ncon,jj,anion,isp
   logical trace,isendmem
   double precision ysum
   save sumngg
   character ch1*1
!
!   write(*,17)mode,iph,ngg
17 format('3Y entering generate_dense_grid: ',i2,i3,i10)
   if(mode.eq.0) then
!      write(*,*)'3Y Generating grid for phase: ',iph
! trace TRUE means generate outpt for each gridpoint
!      trace=.TRUE.
      trace=.FALSE.
      if(iph.eq.1 .and. trace) then
! unit 33 is opened before calling this routine
         sumngg=0
         write(33,43)
43       format('The constituent fractions, y, enclosed within parentheses',&
              'for each sublattice'/'Mole fractions after x:, Gibbs energies',&
              ' after G:'/)
      endif
      if(trace) then
         call get_phase_record(iph,nend)
!         write(33,44)iph,phlista(nend)%name
44       format('Endmembers (EM) and gridpoints (GP) for phase: ',i3,1x,a)
      endif
   else
      trace=.FALSE.
   endif
!   write(*,*)'3Y Getting phase data',mode
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! calculate the number of endmembers and index of first constituent in subl ll
   if(btest(phlista(lokph)%status1,PHIONLIQ)) then
      write(*,*)'We should not be in dense_grid!'
      gx%bmperr=4399; goto 1000
   endif
   nend=1
   inkl(0)=0
   lokph=phases(iph)
   do ll=1,nsl
!      if(btest(phlista(lokph)%status1,PHIONLIQ) .and. ll.eq.2) then
! multiply with charged anions and Va only, add neutrals
!         do jj=1,nkl(2)
! knr(i) is species(i) location but I use constitlist as I have access to it
!            isp=phlista(lokph)%constitlist(nkl(1)+jj)
!            if(btest(splista(isp)%status,SPION) .or. &
!                 btest(splista(isp)%status,SPVA)) then
!               anion=anion+1
!               cycle
!            endif
!         enddo
!         nend=nend*anion+phlista(lokph)%nooffr(2)-anion
!      else
! this is the "normal" number of endmembers
      nend=nend*nkl(ll)
      inkl(ll)=inkl(ll-1)+nkl(ll)
   enddo
!   if(btest(phlista(lokph)%status1,PHIONLIQ)) then
!      write(*,*)'3Y ionic liq: ',anion,nend
!   endif
   ny=inkl(nsl)
   ncon=inkl(nsl)
   negmode: if(mode.lt.0) then
!---------------------------------------------------------
! just determine the number of gridpoints for this phase for global minimimum
! ideal gases should just have the endmembers ....
! Hm, gases with ions??
      ngdim=ngg
      ngg=nend
      lokph=phases(iph)
!      write(*,*)'3Y nend 1: ',mode,ngg
      if(nend.eq.1 .or. nend.gt.100 .or. &
           btest(phlista(lokph)%status1,PHID)) then
! >100 or 1 endmember or ideal phase: only endmembers
         ngg=nend
      else
         ngg=nend
! The limits for various combinations will be adjusted when testing ...
! Max about 20000 gridpoints per phase
!         if(nend.ge.50) then
!         if(nend.ge.100) then
         ngg=ngg+nend*(nend-1)
!         write(*,*)'3Y dense -1A: ',iph,nend,ngg,breaks(5)
! ATTENTION
! The calculation below is not correct, it overestimates a bit the number of 
! gridpoints actually generated but it should not matter so much ... I hope
! When matching a gridpoint in the solution the code to generate the
! gridpoint is used, the code below is just an estimate for allocation
!         if(nend.le.50) then
! Try 60 to handle 53 endmembers in liquid noc2500.TDB from TAF-ID
!         if(nend.le.60) then
         if(nend.le.breaks6) then
            if(nend.gt.breaks(5)) then
               ngg=ngg+nend*(nend-1)+nend*nend*(nend-1)
!               write(*,*)'3Y dense -1B: ',iph,nend,ngg,breaks(4)
            elseif(nend.gt.breaks(5)) then
               ngg=ngg+nend*(nend-1)+2*nend*nend*(nend-1)
            elseif(nend.gt.breaks(4)) then
               ngg=ngg+nend*(nend-1)+3*nend*nend*(nend-1)
            elseif(nend.gt.breaks(3)) then
               ngg=ngg+nend*(nend-1)+4*nend*nend*(nend-1)
            elseif(nend.gt.breaks(2)) then
               ngg=ngg+nend*(nend-1)+5*nend*nend*(nend-1)
!            elseif(nend.gt.breaks(1)) then
            else
               ngg=ngg+nend*(nend-1)+6*nend*nend*(nend-1)
            endif
         endif
!         write(*,*)'3Y dense -1X: ',iph,nend,ngg
      endif
!      write(*,*)'3Y endmembers and gridpoints: ',nend,ngg
!      read(*,11)ch1
11    format(a)
      ny=nend
      goto 1001
   endif negmode
!------------------------------------------------------------
! for mode=0:
!    set gridpoint sitefractions and calculate G
! for mode>0:
!   return sitefractions (for mode=gridpoint number (part of the solution))
!   BUT: The only way to find the site fraction of a gripoint is to generate
!   all gridpoints up the one specified by the value of mode (no G calculation)
!   write(*,*)'3Y ggy: ',mode,iph,nsl,nend,inkl(nsl)
!
   allocate(yfra(inkl(nsl)))
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
! endm(1,1) is constituent in sublattice 1 of first endmember
! endm(2,1) is constituent in sublattice 2 of first endmember
! endm(nsl,2) is constituent in sublattice nsl of second endmember
! endm(1..nsl,nend) are constituents in all sublattices of last endmember
   allocate(endm(nsl,nend))
! inkl(nsl) is the number of fraction variables in the phase
!   allocate(yfra(inkl(nsl)))
!   nofy=inkl(nsl)
! generate endmembers, endm(ll,ie) is set to consituent index in sublattice ll
   je=1
   do ll=1,nsl
      endm(ll,je)=inkl(ll-1)+1
   enddo
100 continue
   je=je+1
! if je>nend we are finished ...
   if(je.gt.nend) goto 120
   do ls=1,nsl
      endm(ls,je)=endm(ls,je-1)
   enddo
   ll=0
110 ll=ll+1
   if(endm(ll,je).lt.inkl(ll)) then
      endm(ll,je)=endm(ll,je)+1
   elseif(ll.lt.nsl) then
      endm(ll,je)=inkl(ll-1)+1
      goto 110
   else
      gx%bmperr=4148; goto 1000
   endif
   goto 100
!---------------------------------------
! We have now generated endm(1..nsl,j)
120 continue
!   write(*,202)'3Y special 1: ',nsl,nend,inkl(nsl),endm(1,2),endm(2,2),&
!        endm(1,nend),endm(1,3),endm(2,3),endm(1,4),endm(2,4)
! now generate all unary, binary and ternary combinations of endmember fractions
! Note the sum of constituent fractions in all sublattices must be unity
! By combining endmember fractions weighted according to ybas, ybin and yter
! we can ensure that
   ngg=0
   l0=0
   l1=0
   lokph=phases(iph)
   endmem1: do iend=1,nend
! we start with a new endmember iend, ybas(1) is 1.00
      l0=ngg
      yfra=yzero
      do ls=1,nsl
         yfra(endm(ls,iend))=ybas(-1)
      enddo
      ngg=ngg+1
      if(mode.eq.0) then
! this is fór a single endmember
!         write(*,201)'3Y end: ',ngg,(yfra(is),is=1,inkl(nsl))
         call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
         if(gx%bmperr.ne.0) goto 1000
         if(garr(ngg).gt.gmax) gmax=garr(ngg)
201      format(a,i5,20(F5.2))
      elseif(ngg.eq.mode) then
         goto 500
      endif
! binary combinations 0.99*y1 + 0.01*y2
      endmem2a: do jend=1,nend
         if(jend.eq.iend) cycle endmem2a
         yfra=zero
!         write(*,202)'3Y special 3: ',endm(1,2)
         do ls=1,nsl
! to generate better start values for fcc-protototype ordering
!            write(*,202)'3Y ls iend endm: ',ls,iend,jend,endm(ls,iend)
202         format(a,10i6)
            yfra(endm(ls,iend))=ybas(0)
            yfra(endm(ls,jend))=yfra(endm(ls,jend))+yter(1)
         enddo
         ngg=ngg+1
         if(mode.eq.0) then
! this is for 0.99*y1 + 0.01*y2
! STRANGE error that destroyed endm after the call to calc_gridpoint!!
! the error was due to wrong size allocated to xarr which is strange as it
! is done elsewhere but the error disapperared when I allocated a larger
! xarr although the allocated one did not seem too small. 
            call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
            if(gx%bmperr.ne.0) goto 1000
            if(garr(ngg).gt.gmax) gmax=garr(ngg)
!            write(*,201)'3Y bin: ',ngg,(yfra(is),is=1,inkl(nsl))
         elseif(ngg.eq.mode) then
            goto 500
         endif
      enddo endmem2a
! ternary combinations
! ibas   -1    0      1      2      3      4      5      6  
! ybas 1.00D0,0.99D0,0.96D0,0.91D0,0.68D0,0.80D0,0.54D0,0.42D0
! ybin               0.03D0,0.07D0,0.25D0,0.15D0,0.36D0,0.35D0
! yter               0.01D0,0.02D0,0.07D0,0.05D0,0.10D0,0.23D0
!      if(nend.gt.50) cycle endmem1
!      if(nend.gt.60) cycle endmem1
      if(nend.gt.breaks6) cycle endmem1
      ibasloop: do ibas=1,6
         if(nend.ge.breaks(5) .and. ibas.eq.2) cycle endmem1
         if(nend.ge.breaks(4) .and. ibas.eq.3) cycle endmem1
         if(nend.ge.breaks(3) .and. ibas.eq.4) cycle endmem1
         if(nend.ge.breaks(2) .and. ibas.eq.5) cycle endmem1
         if(nend.ge.breaks(1) .and. ibas.eq.6) cycle endmem1
         endmem2b: do jend=1,nend
            if(jend.eq.iend) cycle endmem2b
            endmem3: do kend=1,nend
               yfra=zero
               do ls=1,nsl
                  yfra(endm(ls,iend))=ybas(ibas)
                  yfra(endm(ls,jend))=yfra(endm(ls,jend))+ybin(ibas)
                  yfra(endm(ls,kend))=yfra(endm(ls,kend))+yter(ibas)
               enddo
               ysum=zero
               do ls=1,inkl(nsl)
                  ysum=ysum+yfra(ls)
               enddo
               ngg=ngg+1
               if(mode.eq.0) then
! this is for 0.96*y1 + 0.03*y2+0.01*y3
                  call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  if(garr(ngg).gt.gmax) gmax=garr(ngg)
!                  write(*,201)'3Y ter: ',ngg,(yfra(is),is=1,inkl(nsl)),ysum
               elseif(ngg.eq.mode) then
                  goto 500
               endif
            enddo endmem3
         enddo endmem2b
         l1=ngg
      enddo ibasloop
   enddo endmem1
!   write(*,*)'3Y Calculated points: ',ngg
   goto 1000
!----------------------------------------
! here we return the constitution for gridpoint "mode" in the solution
! We must also return the mole fractions ... NO??
500 continue
!    write(*,505)'3Y ggg: ',mode,iph,nsl,inkl(nsl),ny
505 format(a,i7,i4,2x,i2,i5,i10)
!   write(*,510)'3Y ggy: ',mode,ngdim,ny,(yfra(i),i=1,ny)
510 format(a,2i5,i3,10(F6.3))
   do i=1,ny
      yarr(i)=yfra(i)
   enddo
!   write(*,520)'3Y ggx: ',(xarr(is,ngg+ngdim),is=1,nrel)
520 format(a,10(f8.5))
1000 continue
! these will be deallocated by default when exit this subroutine ...
   if(allocated(endm)) then
      deallocate(endm)
      deallocate(yfra)
   endif
1001 continue
! restore original constituent fractions
!   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   errsave=gx%bmperr
   gx%bmperr=0
!   write(*,1010)'3Y Restore ',iph,(ydum(i),i=1,ny)
1010 format(a,'const for ',i3,10(f6.3))
   call set_constitution(iph,1,ydum,qq,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'3Y Error restoring constitution for phase: ',iph,gx%bmperr
   endif
   gx%bmperr=errsave
   if(gx%bmperr.ne.0) write(*,*)'3Y gengrid error: ',gx%bmperr
   return
 end subroutine generate_dense_gridX

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine generate_ionliq_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! generates gridpoints for ionic liquid (also dense)
! Different action depending of the value of mode, 
! for mode<0:  
!    return the number of gridpoints that will be generated for phase iph in ngg
! for mode=0:
!    return garr(i) gibbs energy and xarr(1,i) the compositions of gridpoint i
! for mode>0:
!    return site fractions of gridpoint mode in yarr, number of fractions in ny
!    iph is phase number, ngg is number of gridpoints, nrel number of elements,
! if mode=0:
!    return xarr mole fractions of gridpoints, garr Gibbs energy of gridpoints,
!    ngg is dimension of garr
! if mode>0:
!   "mode" is a gridpoint of this phase in solution, return number of 
!   constituent fractions in ny and fractions in yarr for this gridpoint
! The current constitution is restored at the end of the subroutine
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*),gmax
!\end{verbatim} %+
!
   integer lokph,errsave
   double precision, parameter :: yzero=1.0D-12
   integer je,iend,jend,kend,lend,ll,nend
   double precision ydum(maxconst)
   integer ngdim,nsl
   integer nkl(2),knr(maxconst),inkl(0:2)
   double precision, dimension(:), allocatable :: yfra
   double precision, dimension(:,:), allocatable :: yendm
   double precision sites(2),qq(5)
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
   integer, dimension(:,:), allocatable :: endm
!--------------------------------
! grid is generated by combining end endmembers
! Number of endmemers is N
! First level generate for 4 endmembers including same N**4 pemutations
! 4    1+4  1+3  1+3+4 1+2  1+2+4 1+2+3 1+2+3+4
! 0.52 0.54 0.63 0.65  0.87 0.81  0.98  1.00
! 0.35 0.34 0.34 0.34  0.19 0.19  0.02  -
! 0.11 0.19 0.02 -     0.02 -     -     -
! 0.02 -    -    -     -    -     -     -
! for 2: 1.00 0.98 0.87 0.81 0.65 0.63 0.54 0.52 ...          =16
! for 3: 1.00 0.98/2 0.98/3 ...                               =81
! for N=4..7: N*N*N*N                                         =256-2401
! 1 when 8..15: skip 0.02 except first and last (not dense)              
! 2 when 16-25: 0.02 and 0.11 same except first and last (not dense) N*N*N+..
! 3 when 26-60: 0.02, 0.11 and 0.35 same except first and last N*N+..
! 4 >60 only endmembers
!----------
! IDE: 
! 1) binary liquid of all endmembers N*N (53*53=2809) (incl pure endmembers)
! 2) ternary liquid of all cations with sanme anion or Va (4*11*10*9=2880)
! 3) neutrals?
   integer, parameter, dimension(4) :: breaks=[8,15,25,60]
   double precision, dimension(1:4), parameter:: yf=&
        [0.52D0,0.35D0,0.11D0,0.02D0]
! These are fractions of mixed cations for same anion, not all variants
! Used only when there are many endmembers >15
   double precision, dimension(1:3), parameter:: yfc=&
        [0.42D0,0.33D0,0.25D0]
! for output of gridpoints
   integer l1,ncon,jj,cation,anion,isp,iva,catloop,neutral1
   logical trace,dense
   character ch1*1
!
   if(mode.eq.0) then
!      write(*,*)'3Y Calculating number of gridpoints for ionic liquid'
! trace TRUE means generate outpt for each gridpoint
!      trace=.TRUE.
      trace=.FALSE.
      if(trace) then
! unit 33 is opened before calling this routine
         write(33,43)
43       format('The constituent fractions, y, enclosed within parentheses',&
              'for each sublattice'/'Mole fractions after x:, Gibbs energies',&
              ' after G:'/)
      endif
      if(trace) then
         call get_phase_record(iph,lokph)
!         write(33,44)iph,phlista(lokph)%name
44       format('Endmembers (EM) and gridpoints (GP) for phase: ',i3,1x,a)
      endif
   else
      trace=.FALSE.
   endif
   if(btest(globaldata%status,GSXGRID) .or. & 
            test_phase_status_bit(iph,PHXGRID)) then
!      write(*,*)'Dense grid set'
      dense=.TRUE.
   else
      dense=.FALSE.
   endif
!!   write(*,*)'3Y Getting phase data',mode
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! calculate the number of endmembers and index of first constituent in subl ll
   nend=1
   inkl(0)=0
   inkl(1)=nkl(1)
   cation=nkl(1)
   inkl(2)=nkl(1)+nkl(2)
   lokph=phases(iph)
   if(.not.btest(phlista(lokph)%status1,PHIONLIQ)) then
      write(*,*)'3Y internal error, this phase has not ionic liquid model!'
      gx%bmperr=4399; goto 1000
   endif
! multiply with charged anions and Va only, add neutrals, nsl=2
   anion=0
   do jj=1,nkl(2)
! knr(i) is species(i) location but I use constitlist as I have access to it
      isp=phlista(lokph)%constitlist(nkl(1)+jj)
      if(btest(splista(isp)%status,SPION) .or. &
           btest(splista(isp)%status,SPVA)) then
         anion=anion+1
         cycle
      endif
   enddo
! error when compiling with  -O2
!   nend=inkl(1)*anion+phlista(lokph)%nooffr(2)-anion
   nend=nkl(1)*anion+nkl(2)-anion
   ny=inkl(nsl)
   ncon=inkl(nsl)
!   write(*,45)'3Y liquid endmembers: ',mode,nkl(1),nkl(2),anion,nend
45 format(a,5i5)
   negmode: if(mode.lt.0) then
!---------------------------------------------------------
! just estimate the number of gridpoints for the ionic liquid phase
! pairs of cation+anion, cation+Va, neutrals
      ngdim=ngg
      ngg=nend
      lokph=phases(iph)
!      write(*,*)'3Y nend 1: ',mode,ngg,breaks
      if(nend.eq.1 .or. nend.gt.breaks(4)) then
! Normally about 2000 gridpoints per phase, for dense 10 times more ...
! 1 or >60 endmembers: only endmembers
! if more dense cation grid do not divide cation loop by 2
         ngg=nend+(cation-2)*cation*(cation+1)*anion/2
      elseif(.not.dense .and. nend.gt.breaks(3)) then
! 26..60: between 676-3600
! if more dense cation grid do not divide cation loop by 2
         ngg=nend*nend+(cation-2)*cation*(cation+1)*anion/2
!         write(*,*)'3Y catloop 17: ',ngg,nend,cation,anion,&
!              (cation-2)*cation*(cation+1)*anion/2
      elseif(.not.dense .and. nend.gt.breaks(2)) then
! 16..25: ??
         ngg=nend*nend*2+(cation-2)*cation*(cation+1)*anion/2
!         write(*,*)'3Y catoop 18: ',ngg,nend,cation,anion
      elseif(nend.gt.breaks(1)) then
! 8..15: ??
         ngg=nend*nend*nend
      else
! 2..7, all combinations
         ngg=nend*nend*nend*nend
      endif
!      write(*,*)'3Y endmembers and gridpoints: ',nend,ngg
!      read(*,11)ch1
11    format(a)
      ny=nend
      goto 1001
   endif negmode
!------------------------------------------------------------
! for mode=0:
!    set gridpoint sitefractions and calculate G
! for mode>0:
!   return sitefractions (for mode=gridpoint number (part of the solution))
!   BUT: The only way to find the site fraction of a gripoint is to generate
!   all gridpoints up the one specified by the value of mode (no G calculation)
!   write(*,*)'3Y ggy: ',mode,iph,nsl,nend,inkl(nsl)
!
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
! endm(1,1) is constituent in sublattice 1 of first endmember
! endm(2,1) is constituent in sublattice 2 of first endmember
! endm(nsl,2) is constituent in sublattice nsl of second endmember
! endm(1..nsl,nend) are constituents in all sublattices of last endmember
!   if(mode.gt.0) write(*,*)'3Y allocate endm: ',nsl,nend
   allocate(endm(nsl,nend))
! inkl(nsl) is the number of fraction variables in the phase
!   allocate(yfra(inkl(nsl)))
!   nofy=inkl(nsl)
! generate endmembers, endm(ll,ie) is set to consituent index in sublattice ll
! For neutrals in sublattice 2, sublattice 1 has -99 as constituent
   je=1
   do ll=1,nsl
      endm(ll,je)=inkl(ll-1)+1
   enddo
! we can have an ionic liquid without any cation
   if(knr(1).eq.-99) endm(1,1)=-99
   genend: do while(je.lt.nend)
      je=je+1
! next endmember is first equal to previous
      do ll=1,nsl
         endm(ll,je)=endm(ll,je-1)
      enddo
! increment the constituent in the first sublattice
      if(endm(1,je).lt.inkl(1)) then
         endm(1,je)=endm(1,je)+1
      else
         endm(1,je)=1
         isp=endm(2,je)+1
         if(splista(knr(isp))%charge.eq.zero .and. &
              .not.btest(splista(knr(isp))%status,SPVA)) then
! The next constituent in second sublattice is not Va or a neutral
            exit genend
         else
! the next constituent is an anion or Va
            endm(2,je)=endm(2,je)+1
         endif
      endif
   enddo genend
! we must generate endmembers for neutrals, wildcard -99 in first sublattice
   do iend=je,nend
      endm(1,iend)=-99
      endm(2,iend)=isp
      isp=isp+1
   enddo
! debug check
!   write(*,*)'3Y mode, endmembers, gridpoints: ',mode,nend,ngg
!   write(*,111)(endm(1,je),endm(2,je),je=1,nend)
111 format('3Y   : ',10(i4,i3)/11(i4,i3))
!   gx%bmperr=4399; goto 1000
!---------------------------------------
! We have now generated endm(1..nsl,j)
!120 continue
!   write(*,202)'3Y special 1: ',nsl,nend,inkl(nsl),endm(1,2),endm(2,2),&
!        endm(1,nend),endm(1,3),endm(2,3),endm(1,4),endm(2,4)
! we must allocate and set endmember fractions both for mode 0 and 1
!   if(mode.gt.0) write(*,*)'3Y allocate yendm: ',inkl(2),nend
   allocate(yendm(inkl(2),nend))
   yendm=zero
!   write(*,*)'3Y endmember fractions:',mode,je
   do je=1,nend
      if(endm(1,je).lt.0) then
! for neutrals set all fractions in first sublatice to unity/number of constit
! The G value calculated for the endmember is same as if all were zero
! but when mixing endmembers with cation fractions it may differ ....
         do l1=1,inkl(1)
            yendm(l1,je)=one/inkl(1)
         enddo
      else
         yendm(endm(1,je),je)=one
      endif
      yendm(endm(2,je),je)=one
!      write(*,213)je,(yendm(l1,je),l1=1,inkl(2))
213   format('3Y#',i2,14F5.2/(15F5.2))
   enddo
!   if(mode.gt.0) write(*,*)'3Y allocate yfra: ',nsl,inkl(nsl)
   allocate(yfra(inkl(nsl)))
!---------------------------------------------
! now generate combinations of endmember fractions
! Note the sum of constituent fractions in all sublattices should be unity
   ngg=0
   lokph=phases(iph)
   endmem1: do iend=1,nend
      endmem2: do jend=1,nend
         if(nend.gt.breaks(4) .and. jend.ne.iend) cycle endmem2
         endmem3: do kend=1,nend
            if(nend.gt.breaks(3) .and. kend.ne.jend) cycle endmem3
            endmem4: do lend=1,nend
               if(nend.gt.breaks(2) .and. lend.ne.kend) cycle endmem4
               do jj=1,ncon
                  yfra(jj)=yf(1)*yendm(jj,iend)+yf(2)*yendm(jj,jend)+&
                       yf(3)*yendm(jj,kend)+yf(4)*yendm(jj,lend)
               enddo
               ngg=ngg+1
               if(mode.eq.0) then
                  call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  if(garr(ngg).gt.gmax) gmax=garr(ngg)
                  if(mod(ngg,10000).eq.0) write(*,*)'Calculated ',ngg,&
                       ' gridpoints, more to go'
!                  write(*,211)'3Y ny:',ngg,garr(ngg),(xarr(jj,ngg),jj=1,nrel)
!                  write(*,212)'3Yy: ',(yfra(jj),jj=1,ncon)
211               format(a,i7,1pe12.4,0pf7.4,6f7.4,(3x,10f7.4))
212               format(a,15F5.2,(16F5.2))
               elseif(ngg.eq.mode) then
! when mode>0 we are searching for the constitution of a grid point
! and we must know the yfra here!!
                  goto 500
               endif
            enddo endmem4
         enddo endmem3
      enddo endmem2
   enddo endmem1
!   write(*,*)'3Y Calculated points 1: ',ngg,nend,breaks(2)
!   goto 1000
!   write(*,*)'3Y special cation loop: '
   if(nend.le.breaks(2)) goto 1000
! combinations 3 different cations with same anion
   anion2: do lend=0,anion-1
      catloop=lend*cation+1
!      write(*,*)'3Y catloop: ',lend+1,cation,catloop,ngg
! calculating "c g" followed by "c n" gives better result than "c e", why??
! The reason was that I had forgotten to scale phase amounts with total moles
      endmem1b: do iend=catloop,catloop+cation-3
         endmem2b: do jend=iend+1,catloop+cation-2
            endmem3b: do kend=jend+1,catloop+cation-1
! these loops generate a  more dense grid but give same results ...
! NOTE to use these also increase gridpoints generated with mode=-1
!      endmem1b: do iend=catloop,catloop+cation-1
!         endmem2b: do jend=catloop,catloop+cation-1
!            if(jend.eq.iend) cycle endmem2b
!            endmem3b: do kend=catloop,catloop+cation-1
!               if(kend.eq.jend .or. kend.eq.iend) cycle endmem3b
               if(.not.dense .and. cation.gt.breaks(3)) then
                  write(*,*)'skipping ternary cationloop'
                  cycle endmem3b
!               elseif(mode.eq.0) then
!                  write(*,480)'3Y cations: ',lend+1,iend,jend,kend,ngg
!480               format(a,10i5)
               endif
! mixing of cations with the same anion
               do jj=1,ncon
                  yfra(jj)=yfc(1)*yendm(jj,iend)+yfc(2)*yendm(jj,jend)+&
                       yfc(3)*yendm(jj,kend)
               enddo
               ngg=ngg+1
               if(mode.eq.0) then
                  call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  if(garr(ngg).gt.gmax) gmax=garr(ngg)
!                  if(mod(ngg,10000).eq.0) write(*,*)'Calculated ',ngg,&
!                       ' gridpoints, more to go'
!                  write(*,211)'3Y ny:',ngg,garr(ngg),(xarr(jj,ngg),jj=1,nrel)
!                  write(*,212)'3Yy: ',(yfra(jj),jj=1,ncon)
               elseif(ngg.eq.mode) then
! when mode>0 we are searching for the constitution of a grid point
! and we must know the yfra here!!
                  goto 500
               endif
            enddo endmem3b
         enddo endmem2b
      enddo endmem1b
   enddo anion2
!--------
! skip next loop for the moment
   goto 1000
! combination of 2 different cations with same anion and a neutral
! there are cation*anion endmembers (incl Va as anion), neutrals follow
   neutral1=cation*anion+1
   iva=ngg
   write(*,*)'3Y ionliqgrid3: ',neutral1,ngg
   recip: do lend=0,anion-1
      catloop=lend*cation+1
!      write(*,*)'3Y catloop: ',lend+1,cation,catloop,ngg
! calculating "c g" followed by "c n" gives better result than "c e", why??
! The reason was that I had forgotten to scale phase amounts with total moles
      endmem1c: do iend=catloop,catloop+cation-2
         endmem2c: do jend=iend,catloop+cation-1
! now a neutral, 
            endmem3c: do kend=neutral1,nend
               if(.not.dense .and. cation.gt.breaks(3)) then
                  write(*,*)'skipping ternary cationloop'
                  cycle endmem3c
               elseif(mode.eq.0) then
                  write(*,480)'3Y cations: ',lend+1,iend,jend,kend,ngg
480               format(a,10i5)
               endif
! mixing of cations with the same anion and a neutral
               do jj=1,ncon
                  yfra(jj)=yfc(1)*yendm(jj,iend)+yfc(2)*yendm(jj,jend)+&
                       yfc(3)*yendm(jj,kend)
               enddo
               ngg=ngg+1
               if(mode.eq.0) then
                  call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  if(garr(ngg).gt.gmax) gmax=garr(ngg)
!                  if(mod(ngg,10000).eq.0) write(*,*)'Calculated ',ngg,&
!                       ' gridpoints, more to go'
!                  write(*,211)'3Y ny:',ngg,garr(ngg),(xarr(jj,ngg),jj=1,nrel)
!                  write(*,212)'3Yy: ',(yfra(jj),jj=1,ncon)
               elseif(ngg.eq.mode) then
! when mode>0 we are searching for the constitution of a grid point
! and we must know the yfra here!!
                  goto 500
               endif
            enddo endmem3c
         enddo endmem2c
      enddo endmem1c
   enddo recip
   write(*,*)'3Y ionliqgrid7: ',neutral1,ngg,iva
!
!   write(*,*)'3Y Calculated points 2: ',ngg
! generate combinations of ternary anions if not done above
   goto 1000
!----------------------------------------
! jump here to return the constitution for gridpoint "mode" in the solution
500 continue
   do jj=1,ny
      yarr(jj)=yfra(jj)
   enddo
1000 continue
! these should be deallocated by default when exit this subroutine ...
   if(allocated(endm)) then
      deallocate(endm)
      deallocate(yfra)
      deallocate(yendm)
   endif
1001 continue
! restore original constituent fractions, also if error in this routine
!   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   errsave=gx%bmperr
   gx%bmperr=0
   call set_constitution(iph,1,ydum,qq,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'3Y Error restoring constitution for phase: ',iph,gx%bmperr
   endif
   gx%bmperr=errsave
   if(gx%bmperr.ne.0) write(*,*)'3Y ionliq_grid error: ',gx%bmperr
   return
 end subroutine generate_ionliq_grid

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine generate_fccord_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! This generates grid for a phase with 4 sublattice fcc/hcp ordering
! mode<0 just number of gridpoints in ngg, needed for allocations
! mode=0 calculate mole fraction and G for all gridpoints
! mode>0 return constitution for gridpoint mode in yarr
   implicit none
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*),gmax
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   logical, save :: once=.TRUE.
   integer i,j,k,l,m,n,e
   integer, dimension(:,:), allocatable :: endm
! NOTHING IMPLEMENTED YET
   if(once) write(*,17)
17 format('3Y Special grid for FCC/HCP 4SL ordering not implemented yet')
   once=.FALSE.
! number of constituents in the first 4 sublattices
   n=5
! first calculate the middle value in a Pascal triangle for row n+3
   j=1
   k=1
   l=1
! n+3 can be odd or even
   if(mod(n+3,2).eq.0) then
      m=(n+3)/2
   else
      m=(n+2)/2
   endif
   do i=1,n+3
      j=j*i
      if(i.le.m) then
         k=k*i
      else
         l=l*(n+4-i)
      endif
   enddo
! This is the number of unique permutations
   e=j/(k*l)
   allocate(endm(4,e))
   e=0
   do i=1,n
      do j=i,n
         do k=j,n
            do l=k,n
               e=e+1
               endm(1,e)=i
               endm(2,e)=j
               endm(3,e)=k
               endm(4,e)=l
            enddo
         enddo
      enddo
   enddo
!  write(*,*)'Endmembers ',e
!  do i=1,e
!     write(*,30)(endm(j,i),j=1,4)
!  enddo
30 format(4i3)
1000 continue
   return
 end subroutine generate_fccord_grid

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine generate_charged_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! This generates grid for a phase with charged constituents
! mode<0 just number of gridpoints in ngg, needed for allocations
! mode=0 calculate mole fraction and G for all gridpoints
! mode>0 return constitution for gridpoint mode in yarr
   implicit none
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*),gmax
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer nkl(maxsubl),knr(maxconst),inkl(0:maxsubl)
!   double precision, dimension(:), allocatable :: yfra
   double precision sites(maxsubl),ydum(maxconst),qq(5)
   integer nend,ll,nsl,i1,i2,i3,loksp,mm,lokph,lokcs,np,nm,nn,ncc,iz,loopf
   integer, dimension(:,:), allocatable :: neutral
!   integer, dimension(:), allocatable :: savengg
!   integer ielno(10)
!   double precision stoi(10),smass,qsp
   double precision charge,ratio1,ratio2
   double precision, dimension(:), allocatable :: y1,y2,y3,y4,y5
   real xdum(nrel),gdum
   integer, parameter :: ncf5=5,ncf3=3,alloneut=40000
   integer ncf,maxngg
! These are used to combine endmembers
   double precision, dimension(7), parameter :: nfact=&
        [0.01D0,0.1D0,0.33D0,0.51D0,0.67D0,0.9D0,0.99D0]
   double precision, dimension(ncf5), parameter :: cfact5=&
        [0.05D0,0.3D0,0.5D0,0.7D0,0.95D0]
   double precision, dimension(ncf3), parameter :: cfact3=&
        [0.1D0,0.5D0,0.9D0]
   logical single
! all endmembers will have a record of this type
   type gtp_charged_endmem
! one species number for each sublattice
      integer, dimension(:), allocatable :: constit
      double precision charge
   end type gtp_charged_endmem
   type(gtp_charged_endmem), dimension(:), allocatable :: endmem
! this should be saved or passed as argument
!   save savengg
! we will select 5 or 3 gripoints below
!   write(*,*)'3Y charged grid',iph
!   ncf=ncf5
!   if(.not.allocated(savengg)) then
!      allocate(savengg(noofph))
!      savengg=0
!   endif
   maxngg=ngg
   ngg=0
! get the phase data
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
!
! I will handle this in a very clumsy way by generate all endmembers
! with their charge and then try to combine them to get neutral gridpoints.
   nend=1
   inkl(0)=0
   do ll=1,nsl
      nend=nend*nkl(ll)
! inkl(ll) is the number of constituents up to and including sublattice ll
      inkl(ll)=inkl(ll-1)+nkl(ll)
   enddo
!   write(*,*)'3Y Charged grid for phase ',iph,mode,nend
   if(nend.eq.1) then
! a single endmember, just check it is neutral
      ngg=1
      charge=zero
      do ll=1,nsl
         loksp=knr(ll)
         charge=charge+sites(ll)*splista(loksp)%charge
      enddo
      if(charge.eq.zero) then
         np=ngg
         if(ngg.gt.maxngg) then
            write(*,*)'3Y too many gripoints 2',ngg,maxngg,iph
            gx%bmperr=4399; goto 1000
         endif
         if(mode.eq.0) then
! if mode=0 calculate G for this endmember
!         write(*,*)'3Y a single neutral endmember for ',iph,mode
            call calc_gridpoint(iph,ydum,nrel,xarr(1,ngg),garr(ngg),ceq)
            if(gx%bmperr.ne.0) goto 1000
            if(garr(ngg).gt.gmax) gmax=garr(ngg)
         endif
! finally remove the request for external charge balance !!!
!         write(*,*)'3Y No external charge balance for phase:',iph,lokcs,mode
         call get_phase_compset(iph,1,lokph,lokcs)
         phlista(lokph)%status1=ibclr(phlista(lokph)%status1,PHEXCB)
         goto 1000
      endif
      call get_phase_compset(iph,1,lokph,lokcs)
      write(*,*)'3Y Phase suspended as net charge: ',phlista(lokph)%name
! suspend all composition sets
      do mm=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(mm)
         ceq%phase_varres(lokcs)%phstate=PHSUS
      enddo
      goto 1000
   endif
   np=0
   nm=0
   nn=0
!   write(*,10)'3Y nend: ',iph,nend,0.0D0,(nkl(ll),ll=1,nsl)
10 format(a,2i4,5x,1pe12.4,10i3)
! allocate a record for each endmembers
   allocate(endmem(nend))
   allocate(endmem(1)%constit(nsl))
   charge=zero
   do ll=1,nsl
      endmem(1)%constit(ll)=inkl(ll-1)+1
      loksp=knr(endmem(1)%constit(ll))
!      write(*,*)'3Y species location: ',loksp
!      call get_species_data(loksp,mm,ielno,stoi,smass,qsp)
!      if(gx%bmperr.ne.0) goto 1000
      charge=charge+sites(ll)*splista(loksp)%charge
   enddo
   endmem(1)%charge=charge
!   write(*,15)'3Y end1: ',mode,iph,nsl,charge,1,endmem(1)%constit
15 format(a,3i3,1pe12.4,i4,2x,8i3)
   if(charge.gt.zero) then
      np=np+1
   elseif(charge.lt.zero) then
      nm=nm+1
   else
      nn=nn+1
   endif
!   write(*,10)'3Y endmem: ',1,0,charge,endmem(1)%constit
   emloop: do i2=2,nend
      allocate(endmem(i2)%constit(nsl))
      endmem(i2)%constit=endmem(i2-1)%constit
      sloop: do ll=1,nsl
         if(endmem(i2)%constit(ll).lt.inkl(ll)) then
            exit sloop
         elseif(ll.lt.nsl) then
            endmem(i2)%constit(ll)=endmem(1)%constit(ll)
         else
            exit emloop
         endif
      enddo sloop
      endmem(i2)%constit(ll)=endmem(i2)%constit(ll)+1
      charge=zero
      do mm=1,nsl
         loksp=knr(endmem(i2)%constit(mm))
         charge=charge+sites(mm)*splista(loksp)%charge
      enddo
      endmem(i2)%charge=charge
!      write(*,15)'3Y endx: ',mode,iph,nsl,charge,i2,endmem(i2)%constit
      if(charge.gt.zero) then
         np=np+1
      elseif(charge.lt.zero) then
         nm=nm+1
      else
         nn=nn+1
      endif
   enddo emloop
   mm=nn*nn+np*nm*(nn+np+nm)
! select the number of gridpoints here
   if(mm.gt.2000) then
      ncf=ncf3
   else
      ncf=ncf5
   endif
!      write(*,10)'3Y endmem: ',i2,endmem(i2)%charge,&
!           (splista(knr(endmem(i2)%constit(ll)))%alphaindex,ll=1,nsl)
!   enddo
! calculate the number of gridpoints, consider single endmembers, 
! binary and ternary combinations in a triple loop
   np=0
   nn=0
!   if(mode.ge.0) then
! we have saved the number of gridpoints from the mode=-1 call here
!      np=savengg(iph)
!   write(*,*)'3Y allocate neutral: ',mode,alloneut
! guess a safe value ...
   allocate(neutral(alloneut,0:3))
   neutral=0
   np=0
   loop1: do i1=1,nend
      charge1A: if(endmem(i1)%charge.eq.zero) then
! first endmember neutral, one gridpoint
         np=np+1
         if(mode.ge.0) then
! for generating Y and G we save which endmembers to combine in neutral(*,0)
            neutral(np,0)=0
            neutral(np,1)=i1
         endif
!         write(*,298)'3Y generating 1 gp:  ',np,1,mode,0,i1,0,0
298      format(a,i5,i2,i5,i2,2x,3i3)
      endif charge1A
      loop2: do i2=i1+1,nend
         charge1: if(endmem(i1)%charge.eq.zero) then
! first endmember neutral, that gridpoint already created
            charge2A: if(endmem(i2)%charge.eq.zero) then
!-----------------------------------------------------------------------
! second endmember neutral, generate 7 points between them, (the point at pure
! i2 that will be generated later): 0.01; 0.1; 0.34; 0.51; 0.67; 0.9; 0.99
               if(mode.ge.0) then
                  do ll=1,7
                     neutral(np+ll,0)=1
                     neutral(np+ll,1)=i1
                     neutral(np+ll,2)=i2
                  enddo
               endif
!               write(*,298)'3Y generating 7 gps: ',np+1,3,mode,1,i1,i2,0
               np=np+7
            else
!-----------------------------------------------------------------------
! second endmember has charge, a third endmember needed with opposite charge
               loop3A: do i3=i2+1,nend
                  if(endmem(i2)%charge*endmem(i3)%charge.lt.zero) then
! second and third endmembers have opposite charge, we have ncf gridpoints
! I1_n(I2_(1/c2)I3_(1/c3)_(1-n) where c2 is charge of i2 and c3 charge of i3
                     if(mode.ge.0) then
                        do ll=1,ncf
                           neutral(np+ll,0)=2
                           neutral(np+ll,1)=i1
                           neutral(np+ll,2)=i2
                           neutral(np+ll,3)=i3
                        enddo
                     endif
!                     write(*,298)'3Y generating 3 gps: ',np+1,3,mode,2,i1,i2,i3
                     np=np+ncf
                  endif
               enddo loop3A
            endif charge2A
!=======================================================================
! first endmember has a charge
         elseif(endmem(i2)%charge.eq.zero) then
! second endmember is neutral, we need a third with opposite charge to first
            loop3B: do i3=i2+1,nend
               if(endmem(i1)%charge*endmem(i3)%charge.lt.zero) then
! first and third endmembers have opposite charge, we have ncf gridpoints
! (I1_(1/c1)I3_(1/c3))_n(I2)_(1-n) where c1 is charge of i1 and c3 charge of i3
! where n is 0.1; 0.5; 0.9
                  if(mode.ge.0) then
                     do ll=1,ncf
                        neutral(np+ll,0)=3
                        neutral(np+ll,1)=i1
                        neutral(np+ll,2)=i2
                        neutral(np+ll,3)=i3
                     enddo
                  endif
!                  write(*,298)'3Y generating 3 gps: ',np+1,3,mode,3,i1,i2,i3
                  np=np+ncf
               endif
            enddo loop3B
!-----------------------------------------------------------------------
! first and second endmembers have charge with opposite sign
         elseif(endmem(i1)%charge*endmem(i2)%charge.lt.zero) then 
! we have one gridpoint I1_(1/c1)I2_(1/c2)
            np=np+1
            if(mode.ge.0) then
               neutral(np,0)=4
               neutral(np,1)=i1
               neutral(np,2)=i2
            endif
!            write(*,298)'3Y generating 1 gp:  ',np,1,mode,4,i1,i2,0
!-----------------------------------------------------------------------
            loop3C: do i3=i2+1,nend
               charge3A: if(endmem(i3)%charge.eq.zero) then
! third is neutral, we have ncf more gripoints
! at (I1_(1/c1)I2_(1/c2))_n(I3)_(1-n)
                  if(mode.ge.0) then
                     do ll=1,ncf
                        neutral(np+ll,0)=5
                        neutral(np+ll,1)=i1
                        neutral(np+ll,2)=i2
                        neutral(np+ll,3)=i3
                     enddo
                  endif
!                  write(*,298)'3Y generating 3 gps: ',np+1,3,mode,5,i1,i2,i3
                  np=np+ncf
               elseif(endmem(i1)%charge*endmem(i3)%charge.lt.zero) then
!-------------------------------------------------------------
! all 3 endmembers are charged, those of i2 and i3 have same sign, ncf gridp
! (I1_(1/c1)I2_(1/c2))_n(I1_(1/c1)I3_(1/c3))_(1-n)
                  if(mode.ge.0) then
                     do ll=1,ncf
                        neutral(np+ll,0)=6
                        neutral(np+ll,1)=i1
                        neutral(np+ll,2)=i2
                        neutral(np+ll,3)=i3
                     enddo
                  endif
!                  write(*,298)'3Y generating 3 gps: ',np+1,3,mode,6,i1,i2,i3
                  np=np+ncf
               else
!-------------------------------------------------------------
! all 3 endmembers are charged, those of i1 and i3 have same sign, ncf gridp
! (I1_(1/c1)I2_(1/c2))_n(I2_(1/c2)I3_(1/c3))_(1-n)
                  if(mode.ge.0) then
                     do ll=1,ncf
                        neutral(np+ll,0)=7
                        neutral(np+ll,1)=i1
                        neutral(np+ll,2)=i2
                        neutral(np+ll,3)=i3
                     enddo
                  endif
!                  write(*,298)'3Y generating 3 gps: ',np+1,3,mode,7,i1,i2,i3
                  np=np+ncf
               endif charge3A
            enddo loop3C
!-----------------------------------------------------------------------
! first and second endmembers have charge with same sign
         else
! we need a third endmember with opposite charge
            loop3D: do i3=i2+1,nend
               if(endmem(i1)%charge*endmem(i3)%charge.lt.zero) then
! all 3 endmembers are charged, those of i1 and i2 have same sign, ncf gridp
! (I1_(1/c1)I3_(1/c2))_n(I2_(1/c1)I3_(1/c3))_(1-n)
                  if(mode.ge.0) then
                     do ll=1,ncf
                        neutral(np+ll,0)=8
                        neutral(np+ll,1)=i1
                        neutral(np+ll,2)=i2
                        neutral(np+ll,3)=i3
                     enddo
                  endif
!                  write(*,298)'3Y generating 3 gps: ',np+1,3,mode,8,i1,i2,i3
                  np=np+ncf
               endif
            enddo loop3D
         endif charge1
      enddo loop2
   enddo loop1
!=======================================================================
!   if(mode.eq.0) then
!      write(*,*)'3Y ionic crystal: ',iph,np
!   endif
!   if(mode.lt.0) then
! we have just calculated the number of gridpoints, save and exit
!      write(*,*)'3Y neutral gridpoints: ',np
!      ngg=np
!      savengg(iph)=ngg
!   else
! Generate the composition of the gridpoints from 1-3 endmembers and
! if mode=0 calculate the composition and Gibbs energy for the gridpoints
! if mode>0 return the constitution of gridpoint mode.
! How do I know mode is mode gridpoint in this phase??
!      write(*,29)'3Y we are here?',iph,mode,np,nsl,inkl(nsl)
!29    format(a,10i5)
   ncc=inkl(nsl)
   allocate(y1(ncc))
   allocate(y2(ncc))
   allocate(y3(ncc))
   allocate(y4(ncc))
!   write(*,*)'3Y allocate neutral: ',mode,alloneut,np
! loopf keeps track if several gridpoints belong together
!   if(.not.(allocated(neutral))) then
!      write(*,*)'3Y ionic liquid has not allocated neutral array'
!      allocate(neutral(alloneut,0:3))
!      neutral=0
!   endif
   loopf=0
   ygen: do nm=1,np
! neutral(nm,0) is endmember combination (0 to 8), ,1..3) is endmember index
      nn=neutral(nm,0)
      i1=neutral(nm,1)
      i2=neutral(nm,2)
      i3=neutral(nm,3)
      if(loopf.eq.0) then
! when loopf=0 we have a new set of endmembers, zero yi
         y1=zero
         y2=zero
         y3=zero
      endif
! we must generate all gridpoints to have corrrect loopf
!         if(mode.gt.0) then
!            if(mode.ne.nm) exit
!            cycle
!         endif
! now we must generate correct constituent fractions and calculate G (mode=0)
!         write(*,*)'3Y select case',iph,nn
      select case(nn)
      case default
         write(*,*)'3Y case error in generate_charged_grid!!'
!----------------------- first endmember is neutral, 1 gridpoint
! single neutral endmember
      case(0)
         do ll=1,nsl
            y1(endmem(i1)%constit(ll))=one
         enddo
         y4=y1
!            write(*,300)'3Y gp0 ',nm,nn,loopf,i1,i2,i3,zero,y4
300      format(a,i5,2i2,3i3,1pe10.2,7(0pf6.3),13(f6.3))
!----------------------- first and second endmembers are neutral, 7 gridpoints
! combine with factors: 0.01; 0.10; 0.33; 0.51; 0.67; 0.9; 0.01
      case(1)
         if(loopf.eq.0) then
            do ll=1,nsl
               y1(endmem(i1)%constit(ll))=one
               y2(endmem(i2)%constit(ll))=one
            enddo
         endif
         loopf=loopf+1
         do iz=1,ncc
            y4(iz)=nfact(loopf)*y1(iz)+nfact(8-loopf)*y2(iz)
         enddo
         if(loopf.ge.7) loopf=0
!            write(*,300)'3Y gp1 ',nm,nn,loopf,i1,i2,i3,zero,y4
!----------------------- first endmember is neutral, 2 and 3 charged, 3 gridp
! ratio 2/3 depend on charge, ratio 1/(2+3)
      case(2)
         if(loopf.eq.0) then
            do ll=1,nsl
               y1(endmem(i1)%constit(ll))=one
               y2(endmem(i2)%constit(ll))=one
               y3(endmem(i3)%constit(ll))=one
            enddo
            ratio1=abs(endmem(i3)%charge)/&
                 (abs(endmem(i2)%charge)+abs(endmem(i3)%charge))
            ratio2=abs(endmem(i2)%charge)/&
                 (abs(endmem(i2)%charge)+abs(endmem(i3)%charge))
            do iz=1,ncc
               y2(iz)=ratio1*y2(iz)+ratio2*y3(iz)
            enddo
            charge=ratio1*endmem(i2)%charge+ratio2*endmem(i3)%charge
         endif
         loopf=loopf+1
         do iz=1,ncc
            y4(iz)=cfact5(loopf)*y1(iz)+cfact5(ncf+1-loopf)*y2(iz)
         enddo
         if(loopf.ge.ncf) loopf=0
!            write(*,300)'3Y gp2 ',nm,nn,loopf,i1,i2,i3,charge,y4
!----------------------- first charged, second neutral, third charged, 3 gridp
! ratio 1/3 depend on charge, ratio 2/(1+3): 0.1; 0.5; 0.9
      case(3)
         if(loopf.eq.0) then
            do ll=1,nsl
               y1(endmem(i1)%constit(ll))=one
               y2(endmem(i2)%constit(ll))=one
               y3(endmem(i3)%constit(ll))=one
            enddo
! neutral combination of 1 and 3
            ratio1=abs(endmem(i3)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i1)%charge))
            ratio2=abs(endmem(i1)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i1)%charge))
            do iz=1,ncc
               y1(iz)=ratio1*y1(iz)+ratio2*y3(iz)
            enddo
            charge=ratio1*endmem(i1)%charge+ratio2*endmem(i3)%charge
         endif
         loopf=loopf+1
         do iz=1,ncc
            y4(iz)=cfact5(loopf)*y1(iz)+cfact5(ncf+1-loopf)*y2(iz)
         enddo
         if(loopf.ge.ncf) loopf=0
!            write(*,300)'3Y gp3 ',nm,nn,loopf,i1,i2,i3,charge,y4
!----------------------- first charged, second opposite, 1 gridp
! ratio 1/2 depend on charge
      case(4)
         do ll=1,nsl
            y1(endmem(i1)%constit(ll))=one
            y2(endmem(i2)%constit(ll))=one
         enddo
! neutral combination of 1 and 2
         ratio1=abs(endmem(i2)%charge)/&
              (abs(endmem(i1)%charge)+abs(endmem(i2)%charge))
         ratio2=abs(endmem(i1)%charge)/&
              (abs(endmem(i1)%charge)+abs(endmem(i2)%charge))
         do iz=1,ncc
            y4(iz)=ratio1*y1(iz)+ratio2*y2(iz)
         enddo
         charge=ratio1*endmem(i1)%charge+ratio2*endmem(i2)%charge
!            write(*,300)'3Y gp4 ',nm,nn,loopf,i1,i2,i3,charge,y4
!----------------------- first charged, second opposite, third neutral, 3 gridp
! ratio 1/2 depend on charge, ratio 3(1+2): 0.1; 0.5; 0.9
      case(5)
         if(loopf.eq.0) then
            do ll=1,nsl
               y1(endmem(i1)%constit(ll))=one
               y2(endmem(i2)%constit(ll))=one
               y3(endmem(i3)%constit(ll))=one
            enddo
! neutral combination of 1 and 2
            ratio1=abs(endmem(i2)%charge)/&
                 (abs(endmem(i1)%charge)+abs(endmem(i2)%charge))
            ratio2=abs(endmem(i1)%charge)/&
                 (abs(endmem(i1)%charge)+abs(endmem(i2)%charge))
            do iz=1,ncc
               y1(iz)=ratio1*y1(iz)+ratio2*y2(iz)
            enddo
            charge=ratio1*endmem(i1)%charge+ratio2*endmem(i2)%charge
         endif
         loopf=loopf+1
         do iz=1,ncc
            y4(iz)=cfact5(loopf)*y1(iz)+cfact5(ncf+1-loopf)*y3(iz)
         enddo
         if(loopf.ge.ncf) loopf=0
!            write(*,300)'3Y gp5 ',nm,nn,loopf,i1,i2,i3,charge,y4
!----------------------- all charged, 2 and 3 same sign, 3 gridp
! ratio depend on charge
      case(6)
         if(loopf.eq.0) then
            do ll=1,nsl
               y1(endmem(i1)%constit(ll))=one
               y2(endmem(i2)%constit(ll))=one
               y3(endmem(i3)%constit(ll))=one
            enddo
! neutral combination of 1 and 3
            ratio1=abs(endmem(i3)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i1)%charge))
            ratio2=abs(endmem(i1)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i1)%charge))
            do iz=1,ncc
               y3(iz)=ratio1*y1(iz)+ratio2*y3(iz)
            enddo
            charge=ratio1*endmem(i1)%charge+ratio2*endmem(i3)%charge
!               write(*,410)'3Y gp charge 1+3: ',nm,i1,i2,i3,&
!                    endmem(i1)%charge,endmem(i2)%charge,endmem(i3)%charge,&
!                    ratio1,ratio2,charge
! neutral combination of 1 and 2
            ratio1=abs(endmem(i2)%charge)/&
                 (abs(endmem(i1)%charge)+abs(endmem(i2)%charge))
            ratio2=abs(endmem(i1)%charge)/&
                 (abs(endmem(i1)%charge)+abs(endmem(i2)%charge))
            do iz=1,ncc
               y1(iz)=ratio1*y1(iz)+ratio2*y2(iz)
            enddo
            charge=ratio1*endmem(i1)%charge+ratio2*endmem(i2)%charge
         endif
         loopf=loopf+1
         do iz=1,ncc
            y4(iz)=cfact5(loopf)*y1(iz)+cfact5(ncf+1-loopf)*y3(iz)
         enddo
         if(loopf.ge.ncf) loopf=0
!            write(*,300)'3Y gp6 ',nm,nn,loopf,i1,i2,i3,charge,y4
!----------------------- all charged, 1 and 3 same sign, 3 gridp
! ratio depend on charge
      case(7)
         if(loopf.eq.0) then
            do ll=1,nsl
               y1(endmem(i1)%constit(ll))=one
               y2(endmem(i2)%constit(ll))=one
               y3(endmem(i3)%constit(ll))=one
            enddo
! neutral combination of 1 and 2
            ratio1=abs(endmem(i2)%charge)/&
                 (abs(endmem(i2)%charge)+abs(endmem(i1)%charge))
            ratio2=abs(endmem(i1)%charge)/&
                 (abs(endmem(i2)%charge)+abs(endmem(i1)%charge))
            do iz=1,ncc
               y1(iz)=ratio1*y1(iz)+ratio2*y2(iz)
            enddo
            charge=ratio1*endmem(i1)%charge+ratio2*endmem(i2)%charge
!               write(*,410)'3Y gp charge 1+2: ',nm,i1,i2,i3,&
!                    endmem(i1)%charge,endmem(i2)%charge,endmem(i3)%charge,&
!                    ratio1,ratio2,charge
! neutral combination of 2 and 3
            ratio1=abs(endmem(i3)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i2)%charge))
            ratio2=abs(endmem(i2)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i2)%charge))
            do iz=1,ncc
               y2(iz)=ratio1*y2(iz)+ratio2*y3(iz)
            enddo
            charge=ratio1*endmem(i2)%charge+ratio2*endmem(i3)%charge
         endif
         loopf=loopf+1
         do iz=1,ncc
            y4(iz)=cfact5(loopf)*y1(iz)+cfact5(ncf+1-loopf)*y2(iz)
         enddo
         if(loopf.ge.ncf) loopf=0
!            write(*,300)'3Y gp7 ',nm,nn,loopf,i1,i2,i3,charge,y4
!----------------------- all charged, 1 and 2 same sign, 3 gridp
! ratio depend on charge
      case(8)
         if(loopf.eq.0) then
            do ll=1,nsl
               y1(endmem(i1)%constit(ll))=one
               y2(endmem(i2)%constit(ll))=one
               y3(endmem(i3)%constit(ll))=one
            enddo
! neutral combination of 1 and 3
            ratio1=abs(endmem(i3)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i1)%charge))
            ratio2=abs(endmem(i1)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i1)%charge))
            do iz=1,ncc
               y1(iz)=ratio1*y1(iz)+ratio2*y3(iz)
            enddo
            charge=ratio1*endmem(i1)%charge+ratio2*endmem(i3)%charge
!               write(*,410)'3Y gp charge 1+3: ',nm,i1,i2,i3,&
!                    endmem(i1)%charge,endmem(i2)%charge,endmem(i3)%charge,&
!                    ratio1,ratio2,charge
410         format(a,i4,3i3,6(1pe10.2))
! neutral combination of 2 and 3
            ratio1=abs(endmem(i3)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i2)%charge))
            ratio2=abs(endmem(i2)%charge)/&
                 (abs(endmem(i3)%charge)+abs(endmem(i2)%charge))
            do iz=1,ncc
               y2(iz)=ratio1*y2(iz)+ratio2*y3(iz)
            enddo
            charge=ratio1*endmem(i2)%charge+ratio2*endmem(i3)%charge
         endif
         loopf=loopf+1
         do iz=1,ncc
            y4(iz)=cfact5(loopf)*y1(iz)+cfact5(ncf+1-loopf)*y2(iz)
         enddo
         if(loopf.ge.ncf) loopf=0
!            write(*,300)'3Y gp8 ',nm,nn,loopf,i1,i2,i3,charge,y4
!----------------------- 
      end select
!         if(iph.ge.72) write(*,*)'3Y end select',iph,mode
!===============================================================
! Here we have the neutral constituent fraction in y4
! if mode>0 we have found the requested constitution
!      if(mode.lt.0) then
!         write(*,*)'3Y We should never be here ...'
!         goto 1000
      if(mode.gt.0) then
         if(mode.eq.nm) then
            ny=ncc
            do ll=1,ny
               yarr(ll)=y4(ll)
            enddo
!               if(iph.ge.72) write(*,507)'3Y Solution gp: ',mode,iph,y4
507         format(a,i5,i4,10F7.4)
            goto 1000
         endif
! continue searching for correct gridpoint of the solution
      else
         ngg=ngg+1
         if(ngg.gt.maxngg) then
            write(*,*)'3Y too many gripoints 3',ngg,maxngg,iph
            gx%bmperr=4399; goto 1000
         endif
         call calc_gridpoint(iph,y4,nrel,xarr(1,ngg),garr(ngg),ceq)
         if(gx%bmperr.ne.0) goto 1000
! created a bug here, used ngg instead of nm .... suck
!            if(garr(ngg).gt.gmax) gmax=garr(ngg)
         if(garr(ngg).gt.gmax) gmax=garr(ngg)
!            write(*,512)nm,qq(2),gdum,(xarr(ll,nm),ll=1,nrel)
512      format('3Y gridpoint: ',i5,2(1pe12.4),7(0pF5.2),14(F5.2))
!            if(iph.ge.72) then
!               write(*,*)'3Y calling done'
!               write(*,515)(xarr(ll,nm),ll=1,nrel)
!515            format('3Y yx: ',10F6.3)
!            endif
      endif
   enddo ygen
!
1000 continue
! deallocate creates problems ...
!   if(allocated(savengg)) then
!      deallocate(savengg)
!      deallocate(endmem)
!   endif
!   if(allocated(neutral)) then
!      deallocate(neutral)
!      deallocate(y1)
!      deallocate(y2)
!      deallocate(y3)
!      deallocate(y4)
!   endif
! restore original constitution
!   write(*,*)'3Y Gridpoints for: ',iph,mode,np
   call set_constitution(iph,1,ydum,qq,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'Error restoring constitution for: ',iph,gx%bmperr
   endif
   return
 end subroutine generate_charged_grid

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine calc_gridpoint(iph,yfra,nrel,xarr,gval,ceq)
! called by global minimization routine
! Not adopted to charged crystalline phases as gridpoints have net charge
! but charged gripoints have high energy, better to look for neutral ones ...
   implicit none
   real xarr(*),gval
   integer iph,nrel
   double precision yfra(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! ny just needed for debugging ...
   integer i,lokres
   double precision qq(5),xmol(nrel),ytemp(maxconst),sum
! set constitution and calculate G per mole atoms and composition
!
! BEWARE must be tested for parallel processing
!
   call set_constitution(iph,1,yfra,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
   call calcg(iph,1,0,lokres,ceq)
   if(gx%bmperr.ne.0) goto 1000
   call calc_phase_mol(iph,xmol,ceq)
   if(gx%bmperr.ne.0) goto 1000
!    write(*,15)'3Y gd2: ',iph,lokres,qq(1),&
!         ceq%phase_varres(lokres)%gval(1,1),(xmol(i),i=1,nrel)
!15  format(a,2i4,2e12.4,2x,5(F9.5))
!   sum=zero
!   do i=1,nrel
!      sum=sum+xmol(i)
!   enddo
   do i=1,nrel
      xarr(i)=real(xmol(i))
!     xarr(i)=real(xmol(i)/sum)
   enddo
!   if(qq(1).lt.2.0D-1) then
! number of real atoms less than 20%, a gridpoint with just vacancies ....
   if(qq(1).lt.5.0D-1) then
! number of real atoms less than 50%, a gridpoint with mainly vacancies ....
!      gval=1.0E5
!      gval=1.0E1
      gval=max(1.0E2,real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)))
   elseif(abs(qq(2)).gt.1.0D-14) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! the gridpoint has net charge, qq(2), make gval more positive. 
! Note gval(1,1) is divided by RT so around -5<0
! A better method is needed by combining charged gripoints!!!!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+20*qq(2)**2)
!      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+5*qq(2)**2)
      write(*,*)'3Y Problem with net charge ',iph,qq(2)
      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+qq(2)**2)
      if(ocv()) write(*,66)'3Y charged gp: ',&
           ceq%phase_varres(lokres)%gval(1,1)/qq(1),qq(1),abs(qq(2))
66    format(a,6(1pe12.4))
   else
      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1))
   endif
!    read(*,20)ch1
20  format(a)
1000 continue
! check for parallel
!    jip=omp_get_thread_num()
!    write(*,1010)jip,gval,gx%bmperr
1010 format('3Y Thread ',i3,', gval: ',1pe15.6,', error: ',i6)
   return
 end subroutine calc_gridpoint

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine calc_gridpoint2(iph,yfra,endm,nrel,xarr,gval,ceq)
! called by global minimization routine
! Not adopted to charged crystalline phases as gridpoints have net charge
! but charged gripoints have high energy, better to look for neutral ones ...
! FAILED ATTEMPT TO FIND BUG THAT DESTROYED LOCAL ARRAY IN GENERATE_DENSE_GRID
   implicit none
   real xarr(*),gval
   integer iph,nrel,endm(2,*)
   double precision yfra(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! ny just needed for debugging ...
   integer i,lokres,ibas
   double precision qq(5),xmol(nrel)
   double precision, allocatable :: ydum(:)
! set constitution and calculate G per mole atoms and composition
!
! BEWARE must be tested for parallel processing
!
   call set_constitution(iph,1,yfra,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
   call calcg(iph,1,0,lokres,ceq)
   if(gx%bmperr.ne.0) goto 1000
   call calc_phase_mol(iph,xmol,ceq)
   if(gx%bmperr.ne.0) goto 1000
!    write(*,15)'3Y gd2: ',iph,lokres,qq(1),&
!         ceq%phase_varres(lokres)%gval(1,1),(xmol(i),i=1,nrel)
!15  format(a,2i4,2e12.4,2x,5(F9.5))
   do i=1,nrel
      xarr(i)=real(xmol(i))
   enddo
!   if(qq(1).lt.2.0D-1) then
! number of real atoms less than 20%, a gridpoint with just vacancies ....
   if(qq(1).lt.5.0D-1) then
! number of real atoms less than 50%, a gridpoint with mainly vacancies ....
!      gval=1.0E5
!      gval=1.0E1
      gval=1.0E2
   elseif(abs(qq(2)).gt.1.0D-14) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! the gridpoint has net charge, qq(2), make gval more positive. 
! Note gval(1,1) is divided by RT so around -5<0
! A better method is needed by combining charged gripoints!!!!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+20*qq(2)**2)
!      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+5*qq(2)**2)
      write(*,*)'3Y Problem with net charge ',iph,qq(2)
      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+qq(2)**2)
      if(ocv()) write(*,66)'3Y charged gp: ',&
           ceq%phase_varres(lokres)%gval(1,1)/qq(1),qq(1),abs(qq(2))
66    format(a,6(1pe12.4))
   else
      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1))
   endif
!    read(*,20)ch1
20  format(a)
1000 continue
! check for parallel
!    jip=omp_get_thread_num()
!    write(*,1010)jip,gval,gx%bmperr
1010 format('3Y Thread ',i3,', gval: ',1pe15.6,', error: ',i6)
   return
 end subroutine calc_gridpoint2

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine calcg_endmember(iphx,endmember,gval,ceq)
! calculates G for one mole of real atoms for a single end member
! used for reference states. Restores current composition (but not G or deriv)
! endmember contains indices in the constituent array, not species index
! one for each sublattice
   implicit none
   integer iphx
   double precision gval
   integer endmember(maxsubl)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer iph,ierr,kk0,ll,lokres,nsl,lokph
   integer nkl(maxsubl),knr(maxconst)
   double precision savey(maxconst),sites(maxsubl),yfra(maxconst)
   double precision qq(5),saveg(6)
! when called by matsmin negative iph should be interpreted as index to
! phlista, convert to phase index ... suck
   if(iphx.lt.0) then
      iph=phlista(-iphx)%alphaindex
   else
      iph=iphx
   endif
!   write(*,*)'3Y calcg_endmember: ',iphx,' ',trim(phlista(abs(iphx))%name),&
!        iph,' ',trim(phlista(iph)%name)
!
   call get_phase_data(iph,1,nsl,nkl,knr,savey,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1100
! set constitution to be just the endmember
! It is difficult to make this simpler as one can have magnetic contributions
! to G, this it is not sufficient just to calculate the G function, one must
! calculate TC etc.
!   write(*,16)'3Y call endmember: ',iph,nsl,(endmember(ll),ll=1,nsl)
   yfra=zero
   kk0=0
   do ll=1,nsl
      if(endmember(ll).gt.kk0 .and. endmember(ll).le.kk0+nkl(ll)) then
         yfra(endmember(ll))=one
      else
         write(*,16)'3Y endmember outside range 1: ',iph,ll,endmember(ll),&
              kk0,kk0+nkl(ll)
16       format(a,10i5)
         gx%bmperr=4160; goto 1100
      endif
      kk0=kk0+nkl(ll)
   enddo
!   write(*,17)'3Y set: ',kk0,(yfra(ll),ll=1,kk0)
17 format(a,i3,5(1pe12.4))
   call set_constitution(iph,1,yfra,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! this was necessary using this routine when reference states are
! defined for components
! The calcg below returns lokres but we need it to save G values first!!!
   call get_phase_compset(iph,1,lokph,lokres)
   if(gx%bmperr.ne.0) goto 1000
!   write(*,*)'3Y saving gval: ',lokres,iph
   do ll=1,6
      saveg(ll)=ceq%phase_varres(lokres)%gval(ll,1)
   enddo
! just calculate Gm no derivatives!
   call calcg(iph,1,0,lokres,ceq)
   if(gx%bmperr.ne.0) goto 1000
!   if(qq(1).ge.1.0D-3) then
   if(qq(1).ge.1.0D-2) then
! avoid calculating endmembers with too many vacancies. gval is divided by RT
      gval=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
!      write(*,*)'3Y gval: ',gval,qq(1)
   else
!      write(*,*)'3Y End member has no atoms'
      gx%bmperr=4161; goto 1000
   endif
1000 continue
! restore constitution and gval even if there has been an error flag!!
   ierr=gx%bmperr
   if(gx%bmperr.ne.0) gx%bmperr=0
!   write(*,17)'3Y res: ',kk0,(savey(i),i=1,kk0)
   do ll=1,6
      ceq%phase_varres(lokres)%gval(ll,1)=saveg(ll)
   enddo
   call set_constitution(iph,1,savey,qq,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'3Y Error resetting constitution: ',ierr,gx%bmperr
   endif
! return first error if any
   if(ierr.ne.0) gx%bmperr=ierr
1100 continue
   return
 end subroutine calcg_endmember

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine calcg_endmember6(iph,endmember,gval,ceq)
! calculates G and all derivatevs wrt T and P for one mole of real atoms
! for a single end member, used for reference states. 
! Restores current composition and G (but not deriv)
! endmember contains indices in the constituent array, not species index
! one for each sublattice
   implicit none
   integer iph
   double precision gval(6)
   integer endmember(maxsubl)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ierr,kk0,ll,lokres,lokph,nsl
   integer nkl(maxsubl),knr(maxconst),ics
   double precision savey(maxconst),sites(maxsubl),qq(5),yfra(maxconst)
   double precision saveg(6),savedabnorm(3)
!
   call get_phase_data(iph,1,nsl,nkl,knr,savey,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1100
! set constitution to be just the endmember
! It is difficult to make this simpler as one can have magnetic contributions
! to G, this it is not sufficient just to calculate the G function, one must
! calculate TC etc.
   yfra=zero
   kk0=0
! NOTE abnorm(1) not restored by setting constitution, why?
   ics=1
   call get_phase_compset(iph,ics,lokph,lokres)
   if(gx%bmperr.ne.0) goto 1000
   savedabnorm=ceq%phase_varres(lokres)%abnorm
!   write(*,432)'3Y em6a: ',ceq%phase_varres(lokres)%gval(3,1),&
!        ceq%phase_varres(lokres)%abnorm(1),ceq%phase_varres(lokres)%amfu
432 format(a,6(1pe12.4))
!   write(*,11)'3Y refstate: ',iph,nsl,nkl(1),endmember(1)
11 format(a,10i5)
   do ll=1,nsl
      if(endmember(ll).gt.kk0 .and. endmember(ll).le.kk0+nkl(ll)) then
         yfra(endmember(ll))=one
      else
         write(*,16)'3Y endmember outside range 2',ll,endmember(ll),&
              kk0,nkl(ll)
16       format(a,10i5)
         gx%bmperr=4160; goto 1100
      endif
      kk0=kk0+nkl(ll)
   enddo
!   write(*,17)'set: ',kk0,(yfra(i),i=1,kk0)
17 format(a,i3,5(1pe12.4))
   call set_constitution(iph,1,yfra,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! we do not know lokres here !! YES we do now
!   ics=1
!   call get_phase_compset(iph,ics,lokph,lokres)
!   if(gx%bmperr.ne.0) goto 1000
   do ll=1,6
! Why dividing with qq(1)???
!      saveg(ll)=ceq%phase_varres(lokres)%gval(ll,1)/qq(1)
      saveg(ll)=ceq%phase_varres(lokres)%gval(ll,1)
   enddo
! third argument to calcg is 2 to calculate all derivatives
   call calcg(iph,1,2,lokres,ceq)
   if(gx%bmperr.ne.0) goto 1000
   if(qq(1).ge.1.0D-2) then
! avoid calculating endmembers with too many vacancies. gval is divided by RT
! gval(1..6,1) are G, G.T, G.P, G.T.T, G.T.P and G.P.P
      do ll=1,6
         gval(ll)=ceq%phase_varres(lokres)%gval(ll,1)/qq(1)
      enddo
   else
!      write(*,*)'End member has no atoms'
      gx%bmperr=4161; goto 1000
   endif
! we do not restore values of other properties like TC BMAGN etc
   do ll=1,6
      ceq%phase_varres(lokres)%gval(ll,1)=saveg(ll)
   enddo
1000 continue
   ierr=gx%bmperr
   if(gx%bmperr.ne.0) gx%bmperr=0
! restore constitution
!   write(*,17)'res: ',kk0,(savey(i),i=1,kk0)
   call set_constitution(iph,1,savey,qq,ceq)
   ceq%phase_varres(lokres)%abnorm=savedabnorm
!   write(*,432)'3Y em6b: ',ceq%phase_varres(lokres)%gval(3,1),&
!        ceq%phase_varres(lokres)%abnorm(1),ceq%phase_varres(lokres)%amfu
   if(gx%bmperr.ne.0) then
      if(ierr.ne.0) then
         write(*,*)'Double errors in calcg_endmember: ',ierr,gx%bmperr
      endif
   endif
! return first error if any
   if(ierr.ne.0) gx%bmperr=ierr
1100 continue
   return
 end subroutine calcg_endmember6

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!-\begin{verbatim}
 subroutine calcg_endmember2(lokph,endmember,tpref,gval,ceq)
! calculates G for one mole of real atoms for a single end member
! used for reference states. Restores current composition (but not G or deriv)
! endmember contains indices in the constituent array, not species index
! THIS ONE NOT USED
   implicit none
   integer lokph
   double precision gval,tpref(2)
   integer endmember(maxsubl)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!-\end{verbatim}
   integer ierr,kk0,ll,lokres,nsl,iph
   integer nkl(maxsubl),knr(maxconst)
   double precision savey(maxconst),sites(maxsubl),qq(5),yfra(maxconst),tps(2)
!
   iph=phlista(lokph)%alphaindex
   call get_phase_data(iph,1,nsl,nkl,knr,savey,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1100
! set constitution to be just the endmember
   yfra=zero
   kk0=0
   do ll=1,nsl
      if(endmember(ll).gt.kk0 .and. endmember(ll).le.nkl(ll)) then
         yfra(endmember(ll))=one
      else
         write(*,11)'3Y endmember outside range 3',ll,endmember(ll),nkl(ll)
11       format(a,10i4)
         gx%bmperr=4160; goto 1100
      endif
      kk0=kk0+nkl(ll)
   enddo
!   write(*,17)'3Y set: ',kk0,(yfra(i),i=1,kk0)
!17 format(a,i3,5(1pe12.4))
   call set_constitution(iph,1,yfra,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! if tpref(1) is negative use current T, else use T and P as in tpref
   tps=ceq%tpval
   if(tpref(1).gt.zero) then
      ceq%tpval=tpref
   endif
   call calcg(iph,1,0,lokres,ceq)
   if(gx%bmperr.ne.0) then
      gval=zero
   else
      gval=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
   endif
!--------------------------------
1000 continue
   ceq%tpval=tps
   ierr=0
   if(gx%bmperr.ne.0) then
      ierr=gx%bmperr; gx%bmperr=0
   endif
! restore constitution and T and P
!   write(*,17)'3Y res: ',kk0,(savey(i),i=1,kk0)
   call set_constitution(iph,1,savey,qq,ceq)
   if(gx%bmperr.ne.0) then
      if(ierr.ne.0) then
         write(*,*)'3Y Double errors: ',ierr,gx%bmperr
      endif
! return first error
      gx%bmperr=ierr
   endif
1100 continue
   return
 end subroutine calcg_endmember2

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine find_gridmin(kp,nrel,xarr,garr,xknown,jgrid,phfrac,cmu,trace)
! there are kp gridpoints, nrel is number of components
! composition of each gridpoint in xarr, G in garr
! xknown is the known overall composition
! return the gridpoints of the solution in jgrid, the phase fraction in phfrac
! cmu are the final chemical potentials
   implicit none
   integer, parameter :: jerr=50
   integer kp,nrel
   integer, dimension(*) :: jgrid
   real xarr(nrel,*),garr(*)
   double precision xknown(*),phfrac(*),cmu(nrel)
   logical trace
!\end{verbatim}
!------------------------------------------------------------------------
! How to include conditions on chemical potentials (activities) ??
! Assume mix of conditions N(A)=value and MU(B)=value
! 1. find gridpoints with phase alpha for pure A with highest MU(A),
!    set NP(alpha)=N(A) 
! 2. Tangent plane is \sum_A MU(A)+\sum_B MU(B)
! 2. seach gridpoint with composition N(C) most below tangent plane 
! 3. setup matrix with rows \sum_alpha N(alpha,A)*NP(alpha) = N(A) 
!     replacing the gridpoints one by onw with the new
!    to find a set of gridpoints with positive NP(alpha)
! 4. replace the gridpoint that gives positive NP(alpha) with the new
! 5. repeat from 3 intil no gridpoint lower
! The set of gridpoints should now fullfill the massbalance for A and MU(B)
!---------------------------------------------------------------------------
! inverting just a C1_MO2 phase I got phfmain around 1.0e-13 as smallest
   double precision, parameter :: phfmin=1.0D-15
   real xmat(nrel,nrel),xmatsave(nrel,nrel),xmaxx(nrel)
! used to solve the linear system of equations
   double precision qmat(nrel,nrel+1),qmatsave(nrel,nrel+1)
   double precision zmat(nrel,nrel+1),cmusave(nrel)
   integer notuse(kp),i,ie,iel,ierr,iesave,inerr,inuse,ip,je,jj,jp,jsave,nj
   integer nrel1,nyp,griter,nopure(nrel),gpfail
   double precision phfsave(nrel)
   integer, dimension(jerr) :: removed
   real gmin(nrel),dg,dgmin,gplan,gy,gvvp
! gridpoints that has less difference with the plane than this limit is ignored
   real, parameter :: dgminlim=1.0D-6
   logical checkremoved,linglderr,grindingon,failadd
   character ch1*1
! if trace then open file to write grid
   linglderr=.FALSE.
   grindingon=.TRUE.
   failadd=.TRUE.
   if(trace) then
      write(*,*)'3Y Opening ocgrid.dat to write grid solution'
      open(31,file='ocgrid.dat ',access='sequential')
      write(31,700)nrel,kp,(xknown(inuse),inuse=1,nrel)
700   format('Output from OC gridmin'/' Elements: ',i2,', gridpoints: 'i5,&
           ', composition: '/6(F7.4))
      write(31,*)' Gridpoints in use: '
      do inuse=1,kp
         write(31,710)inuse,(xarr(inerr,inuse),inerr=1,nrel),garr(inuse)
710      format(i6,6(1pe12.4))
      enddo
   endif
! initiallize local arrays
   inuse=kp
   inerr=0
   removed=0
   notuse=0
   cmu=zero
   xmat=zero
   qmat=zero
   do je=1,nrel
      xmat(je,je)=9.9D-1
      jgrid(je)=0
   enddo
   nrel1=nrel+1
   checkremoved=.false.
!   write(*,11)'3Y fm8: ',(xknown(i),i=1,nrel)
!11 format(a,7(F8.4))
! Find the lowest Gibbs energy as close as possible to each pure element
! or with max content
   nopure=0
! skip code below until label 88
   goto 88
!----------------------------------------------------------------------
   do ip=1,kp
!      write(*,118)'3Y pure: ',ip,(xarr(je,ip),je=1,nrel)
118   format(a,i5,10F6.3)
      do je=1,nrel
         if(xarr(je,ip).ge.xmat(je,je)) then
            if(jgrid(je).gt.0) then
               if(garr(ip).gt.gmin(je)) goto 120
!               if(gmin(je).lt.garr(ip) .and. &
!                    xarr(je,ip).eq.xmat(je,je)) then
!                  goto 120
!               endif
            endif
            xmat(je,je)=xarr(je,ip)
            jgrid(je)=ip
            gmin(je)=garr(ip)
!            write(*,*)'3Y pure: ',je,ip,gmin(je)
!         elseif(jgrid(je).eq.0 .and. xarr(je,ip).gt.xmaxx(je)) then
! failed attempt to handle cases with no gridpoint for a pure element
!            xmaxx(je)=xarr(je,ip)
!            nopure(je)=ip
!            gmin=garr(ip)
!            write(*,*)'3Y nopure: ',je,ip,xarr(je,ip)
         endif
120      continue
      enddo
   enddo
! check that we have nrel gridpoints for the pure elements
   do je=1,nrel
      if(jgrid(je).eq.0) then
! no gridpoint assigned to this element!! error (note C in pure fcc has no gp)
!         gx%bmperr=4149; goto 1000
         write(*,122)'3Y Warning, no gridpoint for pure element ',je
!              nopure(je),xmaxx(je)
122      format(a,2i5,2F7.4)
         if(nopure(je).eq.0) then
!            write(*,122)'3Y No solubility in any phase for element ',je
            gx%bmperr=4149; goto 1000
         elseif(xarr(je,nopure(je)).gt.xknown(je)) then
! accept gripoint with highest content of element je outside known composition
            do ie=1,nrel
               xmat(ie,je)=xarr(ie,nopure(je))
            enddo
            gmin(je)=garr(ip)
            phfrac(je)=xknown(je)
         else
           write(*,122)'3Y Composition outside phase compositions for element',&
                 je,nopure(je),xmaxx(je),xknown(je)
            gx%bmperr=4149; goto 1000
         endif
      else
         ip=jgrid(je)
         do ie=1,nrel
            xmat(ie,je)=xarr(ie,ip)
         enddo
         gmin(je)=garr(ip)
         phfrac(je)=xknown(je)
      endif
   enddo
! skip code above-----------------------------------------
88 continue
! set start matrix with chemical potential equal to cmu(1) (max G all gridpoint)
! for all components
   do iel=1,nrel
      phfrac(iel)=xknown(iel)
      xmat(iel,iel)=one
      cmu(iel)=1.0D8
! jgrid value here is dummy ...
      jgrid(iel)=kp+iel
! we must also set gridpoint enegies!!! maybe 1.0D20 better ....
      gmin(iel)=1.0D8
   enddo
! check inial chemical potentials and gripoint energies...
!   write(*,63)'3Ymu: ',(cmu(ie),ie=1,nrel)
!63 format(a,8(1pe10.2))
!   write(*,63)'3Ygm: ',(gmin(ie),ie=1,nrel)
! Add nrel "gridpoints" for the pure elements
!   kp=kp+nrel
! output of start matrix
!   do ip=1,nrel
!      write(*,121)ip,phfrac(ip),cmu(ip),(xmat(je,ip),je=1,nrel)
!   enddo
!   write(*,119)'3Y start: ',0,0,kp,zero,(jgrid(ip),ip=1,nrel)
119 format(a,3i6,1pe12.4/(12i6))
121 format('3Y: ',i2,2(1pe10.2),12(0pf5.2))
123 format('3Y: ',i2,1pe12.4,10(0pf6.3))
! looking for tbase calculation error
!   if(trace) write(*,770)(jgrid(je),je=1,nrel)
!770 format('Initial set of gridpoints: '(/15i5))
   do je=1,nrel
      if(one-xmat(je,je).lt.1.0d-12) then
         cmu(je)=dble(gmin(je))
      else
! we should have a composition for an almost pure element
         gx%bmperr=4150; goto 1000
      endif
   enddo
! copy this into qmat (double precision)
   do ie=1,nrel
      do je=1,nrel
         qmat(je,ie)=dble(xmat(je,ie))
      enddo
   enddo
   qmatsave=qmat
! debug output
!   do je=1,nrel
!      write(*,177)'3Y fm4: ',jgrid(je),phfrac(je),(xmat(ie,je),ie=1,nrel)
!   enddo
177 format(a,i5,1pe11.3,2x,5(1pe11.3))
   gvvp=zero
   do ie=1,nrel
      gvvp=gvvp+xknown(ie)*cmu(ie)
   enddo
   if(trace) then
      write(31,715)nrel
715   format(/'3Y Initial matrix:',i3)
      do je=1,nrel
         write(31,720)'3Y1:',xknown(je),xknown(je),(xmat(ie,je),ie=1,nrel)
      enddo
720   format(a,2F7.4,1x,8f7.3)
      write(31,730)gvvp,(cmu(je),je=1,nrel)
730   format('3Y Gibbs energy: ',1pe14.6/'Chemical potentials: '/6(1pe12.4))
   endif
   griter=0
   gpfail=0
!   write(*,175)'3Y ini: ',gvvp,(cmu(ie),ie=1,nrel)
175 format(a,(1e12.4),2x,6(1pe12.4))
!   write(*,*)'3Y gvvp: ',gvvp
! check we have the correct global composition
!    call chocx('fgm1 ',nrel,jgrid,phfrac,xmat)
!    if(gx%bmperr.ne.0) goto 1000
!    write(*,173)gvvp,(jgrid(i),i=1,nrel)
173 format('3Y fms: ',1pe12.4,10i5)
!   read(*,174)ch1
!174 format(a)
!----------------------------------------------------------
! All setup for starting the search
! search the gridpoint most below the current hyperplane, cmu are 
! the chemical potentials of each pure element for the current lowest plane.
! set notuse nonzero for all points above so they can be skipped next time
! TBASE problem, notuse suspended as a point may fall below later ... ???
200 continue
   griter=griter+1
   dgmin=zero
   nyp=0
!   write(*,*)'3Y Gridpoints in use: ',inuse
!   ovall=zero
!   do i=1,nrel
!      ovall=ovall+xknown(i)*cmu(i)
!   enddo
!   write(*,203)'3Y ff:',inuse,ovall,(cmu(je),je=1,nrel)
203 format(a,i4,1pe12.4,6(1pe11.3))
   pointloop: do jp=1,kp
      included: if(notuse(jp).eq.0) then
         gplan=zero
! first index in xarr is component, second is gridpoint
         do iel=1,nrel
            gplan=gplan+xarr(iel,jp)*cmu(iel)
         enddo
         dg=garr(jp)-gplan
!         write(*,209)'3Y fmz: ',dg,garr(jp),gplan
209      format(a,3(1pe12.4))
         if(dg.gt.zero) then
!            inuse=inuse-1
! we cannot be sure that a point that has a positive value now will always be
! above the surface of the chemical potentials!!!
!            notuse(jp)=1
         else
! if this is the most negative dg we should include it in the solution
            if(dg.lt.dgmin) then
               dgmin=dg; nyp=jp
!               write(*,*)'3Y Lower G: ',griter,nyp,kp,dgmin
            endif
! debugging LC_CsI (61) and SC_CsI (94)
!            if(jp.eq.61 .or. jp.eq.94) &
!                 write(*,44)'3Y extra: ',jp,dg,dgmin,garr(jp),gplan
44          format(a,i5,5(1pe12.4))
         endif
!      else
!         write(*,*)'3X Excluded: ',griter,jp
      endif included
   enddo pointloop
!-----------------------------------------------------------
! OUTPUT AFTER EACH SEARCH
! if lower gridpoint nyp>0
!   write(*,43)griter,nyp,kp,dgmin,(jgrid(ie),ie=1,nrel)
43 format('3Y Finished loop ',i6,' for all gridpoints: ',2i6,1pe12.4/12i6)
! TBASE bug------------------------
!   jp=94
!   do iel=1,nrel
!      gplan=gplan+xarr(iel,jp)*cmu(iel)
!   enddo
!   dg=garr(jp)-gplan
!   write(*,7677)jp,gplan,garr(jp),dg,(xarr(iel,jp),iel=1,nrel)
!7677 format('3Y Gridpoint: ',i5,3(1pe12.4)/(10f7.4))
! TBASE bug------------------------end
! if nyp=0 we have found the lowest tangent plane including the composition
   if(nyp.eq.0 .or. abs(dgmin).lt.dgminlim) then
     if(trace) write(31,*)'Found the solution after iterations: ',griter,dgmin
!      write(31,*)'Found the solution after iterations: ',griter,dgmin
      goto 900
   else
      if(trace) write(*,*)'3Y new gridpoint: ',griter,nyp,dgmin
   endif
!   inuse=inuse-1
   notuse(nyp)=1
!   write(*,211)'3Y ny:',nyp,dgmin,(xarr(ie,nyp),ie=1,nrel)
   if(trace) write(*,212)'3Y Found gridpoint ',nyp,inuse,dgmin,garr(nyp)
211 format(a,i7,1pe12.4,0pf7.4,6f7.4,(3x,10f7.4))
212 format(a,2i8,6(1pe11.3))
!-------------------------------------------------------------------------
   qmat=qmatsave
   do i=1,nrel
      phfsave(i)=phfrac(i)
   enddo
   ie=0
! loop to try to replace an old gridpoint by nyp.  Try to replace all.
300 continue
   ie=ie+1
   if(ie.gt.nrel) then
! tried to change all columns but no solution, error
!      write(*,301)'3Y Failed gp: ',nyp,gpfail,(xarr(i,nyp),i=1,nrel)
301   format(a,i7,i3,1pe10.2,2x,8(0pF5.2))
      gpfail=gpfail+1
      if(griter.gt.10*nrel .and. gpfail.gt.8*nrel) then
! this must be wrong!!  Maybe someone can understand it ...
         if(grindingon) then
            write(*,*)'3Y Grid minimizer problem but grinding on '
            grindingon=.false.
         endif
!         gx%bmperr=4346; goto 1000
      endif
! listing restored solution ......
!      xtx=zero
!      do jjq=1,nrel
!         write(*,177)'3Y flp: ',jgrid(jjq),phfrac(jjq),(xmat(ie,jjq),ie=1,nrel)
!         do jjz=1,nrel
!            xtx(jjz)=xtx(jjz)+phfrac(jjq)*xmat(jjz,jjq)
!         enddo
!      enddo
!      gvv=zero
!      do jjq=1,nrel
!         gvv=gvv+xtx(jjq)*cmu(jjq)
!      enddo
!      write(*,175)'3Y cur: ',gvv,(cmu(ie),ie=1,nrel)
!
! >>>> problem with gas phase test case cho1 with x(c)=.2 x(o)=x(H)=.4
! The gridpoints returned not good, probably due to too many gridpoints ...
!
      if(trace) write(*,*)'3Y Failed when trying to add gridpoint ',nyp
      if(failadd) then
         write(*,*)'3Y Failed trying to use some gridpoints '         
         failadd=.false.
      endif
      if(checkremoved) goto 950
! just ignore this gridpoint and continue, it has been added to notuse
! and will be checked again later as "removed"
      inerr=inerr+1
      if(inerr.gt.jerr) then
         inerr=1
      endif
      removed(inerr)=nyp
      goto 200
   endif
! replace one column in qmat by new composition
   do je=1,nrel
      qmat(je,ie)=dble(xarr(je,nyp))
   enddo
! right hand side are the known composition
   do je=1,nrel
      qmat(je,nrel1)=xknown(je)
   enddo
! solver, note qmat is destroyed inside lingld, nrel is dimension
! qmat matrix with left hand side as additional column i.e. QMAT(1..ND1,ND2)
! phfrac(ND1) is result array, nz number of unknown, ierr nonzero=error
!    do ik=1,nrel1
!       write(*,317)'3Y fm6A: ',(qmat(je,ik),je=1,nrel)
!    enddo
!   do je=1,nrel
!      write(*,55)(qmat(je,iel),iel=1,nrel+1)
!   enddo
!55 format('3Yq:',7(1pe11.3))
   call lingld(nrel,nrel1,qmat,phfrac,nrel,ierr)
   if(ierr.ne.0) then
! error may occur and is not fatal, just try to replace next column
!      write(*,*)'3Y failed replace: ',dgmin
      if(.not.linglderr) then
         if(ocv()) write(*,*)'3Y gridmin warning(s) using lingld: ',ierr,nyp
         linglderr=.TRUE.
      endif
      qmat=qmatsave
      do i=1,nrel
         phfrac(i)=phfsave(i)
      enddo
      goto 300
   endif
!   write(*,299)'3Y q: ',ierr,(phfrac(iel),iel=1,nrel)
!299 format(a,i5,7(1pe10.2))
!   read(*,302)ch1
302 format(a)
!   write(*,*)'3Y fm6B: ',ie,ierr
!   write(*,317)'3Y fm6C: ',(phfrac(i),i=1,nrel)
317 format(a,6(1pe12.4))
!-----------------------
! if solution has only positive values accept this, ierr nonzero if singular
   do je=1,nrel
      if(phfrac(je).le.phfmin .or. phfrac(je).gt.one) then
! maybe problems if known composition have almost zero of some components?
! restore qmat
!          write(*,*)'3Y fm6D: ',je
         qmat=qmatsave
         do i=1,nrel
            phfrac(i)=phfsave(i)
         enddo
         goto 300
      endif
   enddo
!   write(*,*)'3Y Replaced column: ',ie,nyp
! we have found that column ie should be replaced
!--------------------------------------------------
! update xmat, qmatsave and gmin
! as we may fail to find the solution for the chemical potentials later
! keep a copy that can be restored
   iesave=ie
   jsave=jgrid(iesave)
! mark that the replaced gridpoint should be checked again ....
!   write(*,*)'3Y Putting gridpoint back: ',jgrid(ie)
! DO NOT SAVE THE PURE ELEMENT POINTS ... >k
!   write(*,*)'3Y for notuse: ',ie,jgrid(ie),size(notuse)
   if(jgrid(ie).le.size(notuse)) then
      notuse(jgrid(ie))=0
   endif
   jgrid(ie)=nyp
   xmatsave=xmat
   do je=1,nrel
      xmat(je,ie)=xarr(je,nyp)
      qmatsave(je,ie)=dble(xarr(je,nyp))
   enddo
   gmin(ie)=garr(nyp)
!   do ik=1,nrel
!      write(*,317)'3Y fm6F: ',(xmat(je,ik),je=1,nrel)
!   enddo
!   write(*,317)'3Y fm6G: ',(gmin(je),je=1,nrel)
! to solve for the chemical potentials we have ro replace the rows by
! columns, there is a TRANSPOSE command for symmetrical matrices
   do ie=1,nrel
      do je=1,nrel
         zmat(ie,je)=qmatsave(je,ie)
      enddo
   enddo
! we have changed the solution, calculate new chemical potentials
   do je=1,nrel
      zmat(je,nrel1)=gmin(je)
   enddo
!    do ik=1,nrel1
!       write(*,317)'3Y fm8A: ',(zmat(je,ik),je=1,nrel)
!    enddo
   cmusave=cmu
   call lingld(nrel,nrel1,zmat,cmu,nrel,ierr)
   if(ierr.ne.0) then
! this should also be handelled by ignoring the new gridpoint but
! here we must restore the xmat, qmatsave and cmu.
      write(*,*)'3Y Failed to calculate chemical potentials',ierr
!      if(trace) write(*,*)'3Y Error from LINGLD for chem.pot.: ',ierr,nyp
      if(checkremoved) goto 950
      inerr=inerr+1
      if(inerr.gt.jerr) then
         inerr=1
      endif
      removed(inerr)=nyp
      jgrid(iesave)=jsave
      cmu=cmusave
      xmat=xmatsave
      do ie=1,nrel
         do je=1,nrel
            qmatsave(ie,je)=dble(xmat(ie,je))
         enddo
      enddo
! we may have successfully added a removed gridpoint
      if(checkremoved) then
         goto 950
      endif
      goto 200
   endif
! check new chemical potentials ...
!   write(*,63)'3Yny: ',(cmu(ie),ie=1,nrel)
! calculate total G
!   gvv=zero
!   do ie=1,nrel
!      do je=1,nrel
! first index is component, second is species
!         gvv=gvv+xmat(je,ie)*cmu(je)
!      enddo
!   enddo
!   if(trace) write(*,*)'3Y New total G: ',gvv,gvvp
! check if gvv is lower than previous
!   if(gvv.gt.gvvp) then
!      write(*,*)'3Y *** Gibbs energy increased, restore!'
!   endif
!   gvvp=gvv
!----------------------------------------------------------
! debug output as we have changed one gridpoint
!   xtx=zero
!   do jjq=1,nrel
!      write(*,177)'3Y gpf: ',jgrid(jjq),phfrac(jjq),(xmat(ie,jjq),ie=1,nrel)
!      do jjz=1,nrel
!         xtx(jjz)=xtx(jjz)+phfrac(jjq)*xmat(jjz,jjq)
!      enddo
!   enddo
!   gvv=zero
!   do jjq=1,nrel
!      gvv=gvv+xtx(jjq)*cmu(jjq)
!   enddo
!   write(*,175)'3Y ny4: ',gvv,(cmu(ie),ie=1,nrel)
!   write(*,317)'3Y new cmu: ',(cmu(je),je=1,nrel)
!   read(*,321)ch1
!321 format(a)
   gy=zero
   do ie=1,nrel
      gy=gy+xknown(ie)*cmu(ie)
   enddo
!   write(*,199)griter,gvvp,gy
199 format('3Y Gibbs energy changed: ',i5,2(1pe15.6))
   gvvp=gy
!
   if(trace) then
      write(31,740)griter,nyp
740   format(/'Iteration ',i6,' found gridpoint: ',i6,', new matrix:')
      do je=1,nrel
         write(*,720)'3Yz:',phfrac(je),xknown(je),(xmat(je,ie),ie=1,nrel)
      enddo
      write(31,730)gvvp,(cmu(je),je=1,nrel)
   endif
   if(checkremoved) then
      write(*,198)nyp
198   format('3Y Added previously removed gridpoint ',i6)
      goto 950
   endif
!----------------------------------------------
! here we go back to loop through all gridpoints again
!   write(*,*)'3Y New search: ',griter
   goto 200
!==============================================
900 continue
   if(gpfail.gt.0) then
! NOTE is a H-O gas with N(H)=2 N(O)=1 it fails to find H2O because then
! there is just one gridpoint stable!
      write(*,906)gpfail
906   format('3Y Gridmin could not make use of ',i7,' gridpoint(s)')
   endif
!   write(*,*)'3Y Gridmin has found a solution'
!   write(*,316)'3Y fm9A: ',(jgrid(i),i=1,nrel)
!   do ik=1,nrel
!      write(*,317)'3Y fm9B: ',(xmat(je,ik),je=1,nrel)
!   enddo
!   write(*,317)'3Y fm9C: ',(garr(je),je=1,nrel)
!   write(*,317)'3Y fm9D: ',(cmu(je),je=1,nrel)
!   write(*,317)'3Y fm9E: ',(phfrac(je),je=1,nrel)
316 format(a,10i5)
   nj=0
!    do j=1,jerr
!       if(removed(j).gt.0) then
!          write(*,*)'3Y Failed testing gridpoint ',removed(j)
!          nj=nj+1
!       endif
!    enddo
950 continue
   nj=0
   checkremoved=.true.
!   write(*,*)'3Y Checking removed gridpoints',inerr
!   xtx=zero
!   do jjq=1,nrel
!      write(*,177)'3Y flp: ',jgrid(jjq),phfrac(jjq),(xmat(ie,jjq),ie=1,nrel)
!      do jjz=1,nrel
!         xtx(jjz)=xtx(jjz)+phfrac(jjq)*xmat(jjz,jjq)
!      enddo
!   enddo
!   gvv=zero
!   do jjq=1,nrel
!      gvv=gvv+xtx(jjq)*cmu(jjq)
!   enddo
!   write(*,175)'3Y cur: ',gvv,(cmu(ie),ie=1,nrel)
!----------------
   testloop: do jj=1,inerr
      jp=removed(jj)
!      write(*,*)'3Y Checking removed gridpoint: ',jj,jp
      if(jp.gt.0) then
         gplan=zero
         do iel=1,nrel
            gplan=gplan+xarr(iel,jp)*cmu(iel)
         enddo
         dg=garr(jp)-gplan
         if(dg.lt.zero) then
!            if(trace) write(*,985)jp,dg,garr(jp),gplan
!            write(*,982)jp,dg,garr(jp),gplan
982         format('3Y Removed gridpoint ',i5,' is below surface ',3(1pe12.4))
! try to include it ....
            ie=0
            removed(jj)=-jp
            nyp=jp
            goto 300
         else
!            write(*,983)jp,dg
983         format('3Y Removed gridpoint ',i5,' above surface ',1pe12.4)
            removed(jj)=-jp
         endif
      endif
   enddo testloop
   if(inerr.gt.0 .and. nj.eq.0) then
!      if(trace) write(*,986)inerr
986   format('3Y None of the ',i3,' removed gridpoints below final surface')
   endif
   if(trace) write(*,771)(jgrid(je),je=1,nrel)
771 format('3Y Final set of gridpoints: '(/15i5))
!   xtx=0
!   do iii=1,nrel
!      write(*,987)jgrid(iii),phfrac(iii),(xarr(i,jgrid(iii)),i=1,nrel)
!987   format('3Y GP: ',i5,F7.4,2x,6F9.6)
!      do j=1,nrel
!         xtx(j)=xtx(j)+phfrac(iii)*xarr(j,jgrid(iii))
!      enddo
!   enddo
!   write(*,988)(xtx(i),i=1,nrel)
!988 format('3Y MF: ',6F9.6)
!
!    call chocx('fgme ',nrel,jgrid,phfrac,xmat)
1000 continue
!   write(*,*)'3Y exit find_gridmin'
   return
 end subroutine find_gridmin

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine merge_gridpoints(nv,iphl,aphl,nyphl,yphl,trace,nrel,xsol,cmu,ceq)
!
! BEWARE not adopted for parallel processing
!
! if the same phase has several gridpoints check if they are really separate
! (miscibility gaps) or if they can be murged.  Compare them two by two
! nv is the number of phases, iphl(i) is the index of phase i, aphl(i) is the
! amount of phase i, nyphl is the number of site fractions for phase i, 
! and yphl is the site fractions packed together
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer nv,nrel
   integer, dimension(*) :: iphl,nyphl
   double precision, dimension(*) :: aphl,yphl,cmu
   logical trace
   real xsol(maxel,*)
!\end{verbatim}
   integer i,ip,iph,jp,jump,kk,klast,kp,lokres,nm,jj,mj,lokph,j
   integer notuse(nv),incy(nv)
   double precision ycheck(maxconst),qq(5),xerr(maxel),xfromy(maxel)
   double precision summu,sumam
   logical igen
   real xmix(maxel)
   double precision a1,a2,gdf,gval1,gval2,gval3,gval4,gval5,gmindif
!
! gmindif is the value to accept to merge two gridpoints
! It should be a variable that can be set by the user for finetuning
!   write(*,*)'3Y Now we try to merge gridpoints in the same phase'
!   write(*,7)'3Y Merge_gridpoints is dissabled for the moment',nv
!7  format(a,i3)
! NOTE, always merge gripoints in ideal phases like gas
   goto 1100
! UNFINISHED
!   write(*,*)'3Y Entering merge_gridpoints',nv,ceq%gmindif
!---------------------
   gmindif=ceq%gmindif
   notuse=0
   nm=0
   write(*,67)'3Y yl: ',(aphl(i),nyphl(i),i=1,nv)
67 format(a,20(F7.3,i4))
   incy(1)=1
   do i=2,nv
      incy(i)=incy(i-1)+nyphl(i-1)
   enddo
! start points of fractions for all gridpoints
   write(*,68)'3Y ys: ',(incy(i),i=1,nv)
68 format(a,20i5)
   summu=zero
   xerr=zero
! constitution of solution gridpoints
   do jp=1,nv
      write(*,69)'3Y y:',(yphl(incy(jp)+i-1),i=1,nyphl(jp))
   enddo
69 format(a,(12F6.2))
! this calculate the overall composition from gridpoints
   do jp=1,nv
      summu=summu+aphl(jp)
      do i=1,nrel
         xerr(i)=xerr(i)+aphl(jp)*xsol(i,jp)
      enddo
   enddo
   write(*,73)'3Y in1: ',summu,(xerr(i),i=1,nrel)
73 format(a,F5.2,2x,9(f7.4))
!----------------------------------------------
100 continue
   igen=.false.
   firstgp: do jp=1,nv-1
      write(*,*)'3Y notuse: ',jp,notuse(jp)
      if(notuse(jp).ne.0) cycle firstgp
      secondgp: do kp=jp+1,nv
         write(*,*)'3Y notuse: ',kp,notuse(kp)
         if(notuse(kp).ne.0) cycle secondgp
         sameph: if(iphl(jp).eq.iphl(kp)) then
            iph=iphl(jp)
            lokph=phases(iph)
            if(btest(phlista(lokph)%status1,PHID)) then
! always merge gridpoints in ideal gas
               goto 200
            endif
! do not merge gridpoints in other phases
            if(btest(globaldata%status,GSNOMERGE)) cycle secondgp
! calculate G at 0 and 1 and  0.25, 0.5, 0.75 mix of gridpoints
! if any of these abouve the line between any two others do not merge
! as merged gridpoints are below the initial common tahnegt plane we cannot
! use that as reference
            call set_constitution(iph,1,yphl(incy(jp)),qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call calcg(iph,1,0,lokres,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gval1=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
! second point
            call set_constitution(iph,1,yphl(incy(kp)),qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call calcg(iph,1,0,lokres,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gval5=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
! take middle point
            a1=5.0D-01
            a2=5.0D-01
            do i=0,nyphl(jp)-1
               ycheck(i+1)=a1*yphl(incy(jp)+i)+a2*yphl(incy(kp)+i)
            enddo
            call set_constitution(iph,1,ycheck,qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call calcg(iph,1,0,lokres,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gval3=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
! Check if this is above the mean of gval1 and gval5, if so quit merge
            gdf=gval3-a1*gval1-a2*gval5
! merge require that difference is less than gmindif or phase ideal
            if(gdf.gt.gmindif) then
! middle is higher, no merge 1-3-5
               write(*,830)'3Y not merged: ',jp,kp,gdf,iphl(jp),gmindif
               cycle secondgp
            endif
! calculate G at 0.25
            a1=7.5D-01
            a2=2.5D-01
            do i=0,nyphl(jp)-1
               ycheck(i+1)=a1*yphl(incy(jp)+i)+a2*yphl(incy(kp)+i)
            enddo
            call set_constitution(iph,1,ycheck,qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call calcg(iph,1,0,lokres,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gval2=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
! Check if this is above the mean of gval1 and gval5, if so quit merge
            gdf=gval2-a1*gval1-a2*gval5
            if(gdf.gt.gmindif) then
! middle is higher, no merge 1-2-5
               write(*,830)'3Y not merged: ',jp,kp,gdf,iphl(jp),gmindif
               cycle secondgp
            else
! also compare between gval1 and gval3
               gdf=gval2-a1*gval1-a2*gval3
               if(gdf.gt.gmindif) then
! gval2 is s higher, no merge 1-2-3
                  write(*,830)'3Y not merged: ',jp,kp,gdf,iphl(jp),gmindif
                  cycle secondgp
               endif
            endif
! finally calculate at 0.75
! calculate G at 0.25
            a1=2.5D-01
            a2=7.5D-01
            do i=0,nyphl(jp)-1
               ycheck(i+1)=a1*yphl(incy(jp)+i)+a2*yphl(incy(kp)+i)
            enddo
            call set_constitution(iph,1,ycheck,qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call calcg(iph,1,0,lokres,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gval4=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
! Check if this is above the mean of gval1 and gval5, if so quit merge
            gdf=gval4-a1*gval1-a2*gval5
            if(gdf.gt.gmindif) then
! gval4 is higer mo merge 1-4-5
               write(*,830)'3Y not merged ',jp,kp,gdf,iphl(jp),gmindif
               cycle secondgp
            else
               gdf=gval4-a1*gval1-a2*gval5
               if(gdf.gt.gmindif) then
! gval4 is higer mo merge 1-4-5
                  write(*,830)'3Y not merged ',jp,kp,gdf,iphl(jp),gmindif
                  cycle secondgp
               endif
               gdf=gval4-a1*gval1-a2*gval5
               if(gdf.gt.gmindif) then
! gval4 is higer mo merge 1-4-5
                  write(*,830)'3Y not merged ',jp,kp,gdf,iphl(jp),gmindif
                  cycle secondgp
               endif
            endif
! compared 1-3-5, 1-2-5, 1-2-3, 1-4-5, 3-4-5, 2-3-5               
! in no case the middle point was above, that means merge
!--------------------------------------------- here we merge !!
200         continue
! gridpoint in ideal phase or point in between has lower G, merge
            write(*,830)'3Y merging:  ',jp,kp,gdf,iphl(jp),&
                 aphl(jp),aphl(kp)
830         format('3Y GridP ',a,2i3,1pe12.4,' in phase ',i3,2(e12.4))
! If merging use correct phase amounts
            a1=aphl(jp)/(aphl(jp)+aphl(kp))
            a2=aphl(kp)/(aphl(jp)+aphl(kp))
            write(*,162)'3Y p1:',a2,(yphl(incy(jp)+j),j=0,nyphl(jp)-1)
            write(*,162)'3Y p2:',a2,(yphl(incy(kp)+j),j=0,nyphl(kp)-1)
162         format(a,1pe12.4,12(0pF5.2))
! The gridpoint jp has new amount, composition and constitution
! SURPRISE: adding together constituent fractions does not reproduce
! the correct molefractions if the constituents are molecules .... ????
            aphl(jp)=aphl(jp)+aphl(kp)
            do i=0,nyphl(jp)-1
               yphl(incy(jp)+i)=a1*yphl(incy(jp)+i)+a2*yphl(incy(kp)+i)
            enddo
            call set_constitution(iph,1,yphl(incy(jp)),qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
! extract correct mole fractions
            call calc_phase_mol(iph,xfromy,ceq)
            write(*,162)'3Y ym:',0.0D0,(yphl(incy(jp)+i),i=0,nyphl(jp)-1)
            write(*,162)'3Y xy:',0.0D0,(xfromy(i),i=1,nrel)
! calculate mole fractions from xsol to compare!!
            do i=1,nrel
               xerr(i)=a1*xsol(jp,i)+a2*xsol(kp,i)
            enddo
            write(*,162)'3Y xj+xk:',0.0D0,(xerr(i),i=1,nrel)
            do i=1,nrel
               xsol(i,jp)=xerr(i)
            enddo
            igen=.true.
            nm=nm+1
! Mark the gripoint that has disappeared
            iphl(kp)=-iphl(kp)
            notuse(kp)=1
! check overall composition of solution ...
            summu=zero
            xerr=zero
            do i=1,nrel
               if(iphl(i).lt.0) cycle
               summu=summu+aphl(i)
!               write(*,*)'3Y point: ',i,aphl(i)
               do jj=1,nrel
                  xerr(jj)=xerr(jj)+aphl(i)*xsol(jj,i)
               enddo
            enddo
            write(*,73)'3Y nu: ',summu,(xerr(jj),jj=1,nrel)
! the chemical potentials has changed but how?  Approximate the change by
! making gmindif more negative for each merge (does not affect ideal phases)
            gmindif=2.0D0*gmindif
! after merging always restart loop
            goto 100
         endif sameph
      enddo secondgp
   enddo firstgp
! if two gridpoints merged compare all grispoints again
   if(igen) goto 100
!----------------------------------------
! shift fractions for the removed phases
450 continue
!   write(*,*)'3Y at label 450: ',nm
   klast=0
   do jp=1,nv
      klast=klast+nyphl(jp)
   enddo
!
! uncomment listing here if error moving fractions
!    write(*,502)nv,(iphl(i),i=1,nv)
!    write(*,502)0,(incy(i),i=1,nv)
!    write(*,502)klast,(nyphl(i),i=1,nv)
502 format('3Y check1: ',i3,2x,20i4)
!    kk=0
!    do j=1,nv
!       write(*,510)j,(yphl(i),i=kk+1,kk+nyphl(j))
!       kk=kk+nyphl(j)
!    enddo
!
   kk=0
   jp=1
   do while(jp.lt.nv)
      if(iphl(jp).lt.0) then
! shift all fractions down.  klast should be updated each shift but ...
         jump=nyphl(jp)
!          write(*,503)jp,kk,klast,jump
503      format('3Y check3: ',5i5)
!          write(*,555)'3Y nyy1: ',(yphl(ip),ip=kk+1,kk+jump)
555      format(a,6(1pe12.4))
         do ip=kk+1,klast-jump
            yphl(ip)=yphl(ip+jump)
         enddo
!          write(*,555)'3Y nyy2: ',(yphl(ip),ip=kk+1,kk+jump)
         do kp=jp,nv-1
            iphl(kp)=iphl(kp+1)
            aphl(kp)=aphl(kp+1)
            nyphl(kp)=nyphl(kp+1)
         enddo
         nv=nv-1
      else
         kk=kk+nyphl(jp)
         jp=jp+1
      endif
500   continue
   enddo
   if(iphl(nv).lt.0) nv=nv-1
!   write(*,*)'3Y final number of gridpoints: ',nv
! list overall composition with merged gridpoints
   summu=zero
   xerr=zero
! this calculate the overall composition from gridpoints
!   write(*,87)'3Y aphl: ',nv,(aphl(jp),jp=1,nv)
87 format(a,i2,7(1pe10.2))
   do jp=1,nv
      summu=summu+aphl(jp)
      do i=1,nrel
         xerr(i)=xerr(i)+aphl(jp)*xsol(i,jp)
      enddo
   enddo
!   write(*,73)'3Y in2: ',summu,(xerr(i),i=1,nrel)
! uncomment here if problems shifting fractions
!    write(*,502)nv,(iphl(i),i=1,nv)
!    write(*,502)0,(incy(i),i=1,nv)
!    write(*,502)klast,(nyphl(i),i=1,nv)
!    kk=0
!    do j=1,nv
!       write(*,510)j,(yphl(i),i=kk+1,kk+nyphl(j))
!       kk=kk+nyphl(j)
!    enddo
! if there are two or more gripoints in the same phase we have a 
! miscibility gap and may have to create miscibility gaps.
!
! >>>> unfinished
!
510 format(i3,':',6(1pe12.4))
1000 continue
   if(ocv()) write(*,*)'3Y At return from merge_gridpoints: ',nv
   return
!------------------------------------------
! temporary fix to avoid creating several composition sets in ideal gas
1100 continue
!   write(*,1102)'3Y merge ideal: ',nv,(iphl(jp),jp=1,nv)
1102 format(a,i2,20i3)
   nm=0
   notuse=0
   incy(1)=1
   do i=2,nv
      incy(i)=incy(i-1)+nyphl(i-1)
   enddo
1110 continue
   igen=.FALSE.
   do jp=1,nv-1
      do kp=jp+1,nv
         if(notuse(kp).ne.0) cycle
         if(iphl(jp).eq.iphl(kp)) then
            iph=iphl(jp)
            lokph=phases(iph)
            if(btest(phlista(lokph)%status1,PHID)) then
! add together gridpoints in ideal phases (gas)
!               write(*,*)'3Y merging gridpoints in ideal phase'
               sumam=aphl(jp)+aphl(kp)
               a1=aphl(jp)/sumam
               a2=aphl(kp)/sumam
               aphl(jp)=aphl(jp)+aphl(kp)
! sum the constituent fractions 
               do i=0,nyphl(jp)-1
                  yphl(incy(jp)+i)=a1*yphl(incy(jp)+i)+a2*yphl(incy(kp)+i)
               enddo
! sum also mole fractions!!
               do i=1,nrel
                  xsol(i,jp)=a1*xsol(i,jp)+a2*xsol(i,kp)
               enddo
               notuse(kp)=1
               igen=.TRUE.
               nm=nm+1
               iphl(kp)=-iphl(kp)
            endif
         endif
      enddo
   enddo
   if(igen) goto 1110
   if(nm.eq.0) goto 1000
   goto 450
!
 end subroutine merge_gridpoints

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine set_metastable_constitutions(ngg,nrel,kphl,ngrid,xarr,garr,&
      nr,iphl,cmu,ceq)
! this subroutine goes through all the metastable phases
! after a global minimization and sets the constituion to the most
! favourable one.  Later care should be taken that composition set 2 
! and higher are not set identical or equal to the stable
! nrel  number of components
! ngg   number of gridpoints
! kphl  array with first points calculated for phase(i) in garr
! ngrid array with last  points calculated for phase(i) in garr
! garr  array with Gibbs energy/RT for each gridpoint
! xarr  matix with composition in all gridpoints
! nr    is the number of stable phases in the solution
! iphl  array with the phase numbers of the stable phases (not ordered)
! cmu   are the chemical potentials/RT of the solution
! ceq   equilibrium record
! called by global_gridmin
   implicit none
   integer ngg,nrel,nr
   integer, dimension(*) :: kphl,ngrid,iphl
   double precision, dimension(*) :: cmu
   real garr(*),xarr(nrel,*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ig1,ign,ip,iph,ics,jph,lokcs,lokph,mode,ny,ie,ig,kp,i
   double precision yarr(maxconst),qq(5),xxx,dgmin,gmax
   real dg,gplan
!   integer idum(1000)
!   write(*,*)'3Y Entering set_metastable'
! loop through the gridpoints for all unstable phases and insert the
! constitution that is closest to be stable
!   write(*,7)'3Y set_meta: ',kp,nrel,nr,(iphl(i),i=1,nr)
7  format(a,i9,2i4,2x,10i3)
!
! WOW kphl and ngrid must be incremented by 1 ???? why BUG??
!
!   do ig=1,noofph
!      write(*,*)'3Y phase grid: ',ig,kphl(ig+1),ngrid(ig+1)
!   enddo
   phloop: do iph=1,noofph
      do jph=1,nr
         if(iph.eq.iphl(jph)) then
!            write(*,*)'3Y Phase is stable',iph
            cycle phloop
         endif
      enddo
      call get_phase_record(iph,lokph)
      if(gx%bmperr.ne.0) goto 1000
! check if all composition sets are suspended
      do ics=1,phlista(lokph)%noofcs
         call get_phase_compset(iph,ics,lokph,lokcs)
         if(gx%bmperr.ne.0) goto 1000
! new -4 hidden, -3 susp, -2 dorm, -1,0,1 entered, 2 fixed
         if(test_phase_status(iph,ics,xxx,ceq).ge.PHENTUNST) goto 60
      enddo
      cycle phloop
! this phase is not suspended and not stable, find gridpoints
60    continue
! these are the grid points belonging to phase iph
! NOTE add 1 to phase index!!
      ig1=kphl(iph+1)
      ign=ngrid(iph+1)
!      if(ocv()) write(*,69)'3Y Searching gridpoints for: ',iph,ics,ig1,ign
69    format(a,2(i3,1x),2x,3(i6,1x))
!      write(*,68)'3Y Chempot: ',(cmu(ie),ie=1,nrel)
68    format(a,6(1pe12.4))
! if ig1=0 there are no gridpoints for this phase, it is suspended or dormant
      if(ig1.le.0) cycle phloop
! if ign=ig1 there is a single gridpoint, just stet dgm
      dgmin=-1.0d12
      ip=0
! search for gripoint closeset to stable plane
      do ig=ig1,ign
         gplan=zero
         do ie=1,nrel
            gplan=gplan+xarr(ie,ig)*cmu(ie)
         enddo
         dg=gplan-garr(ig)
!         write(*,74)'3Y dgx: ',ig,gplan,garr(ig),dg,dgmin,(xarr(i,ig),i=1,nrel)
74       format(a,i5,2(1pe10.2),2x,6(0pf5.2))
         if(ign.eq.ig1) then
!            write(*,*)'3Y Single gridpoint, dgm: ',dg
            dgmin=dg
            call set_driving_force(iph,1,dgmin,ceq)
            cycle phloop
         endif
         if(abs(dg).lt.abs(dgmin)) then
            ip=ig
            dgmin=dg
!            write(*,77)'3Y low: ',ig,gplan,garr(ig),dg,dgmin
77          format(a,i5,4(1pe12.4))
         endif
      enddo
      if(ocv()) write(*,79)'3Y Least unstable gridpoint: ',iph,ics,ig1,ign,dgmin
79    format(a,4(i6,1x),1pe12.4)
!      if(ip.eq.0 .or. dgmin.gt.zero) then
!         write(*,*)'3Y This gridpoint stable: ',ip,dgmin
!         write(*,*)'3Y data: ',ip,dgmin
!      endif
! retrieve constitution for this gridpoint and insert it in phase
! must provide mode and iph. The subroutine returns ny and yarr
! mode is the gridpoint in the phase, subtract ig1-1
      mode=ip-ig1+1
! 
!      if(ocv()) write(*,78)'3Y calling gengrid: ',iph,ig1,ip,ign,mode,dgmin
78    format(a,5i7,1pe12.4)
! find the constitution of this gridpoint
!      call generate_grid(mode,iph,ign,nrel,xarr,garr,ny,yarr,idum,ceq)
!      call generate_grid(mode,iph,ign,nrel,xarr,garr,ny,yarr,gmax,ceq)
      call generic_grid_generator(mode,iph,ign,nrel,xarr,garr,ny,yarr,gmax,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,451)(yarr(i),i=1,ny)
451   format('3Y fractions: ',6(F10.6))
      call set_constitution(iph,1,yarr,qq,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,452)iph,ics,(yarr(i),i=1,ny)
452   format('3Y my: ',2i2,6(1pe10.3))
! set driving force also ...
      call set_driving_force(iph,1,dgmin,ceq)
500   continue
   enddo phloop
1000 continue
!   write(*,*)'3Y Finished set_metastable'
   return
 end subroutine set_metastable_constitutions
 
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine set_metastable_constitutions2(pph,nrel,nphl,iphx,xarr,garr,&
      nr,iphl,cmu,ceq)
! this subroutine goes through all the metastable phases
! after a global minimization and sets the constituion to the most
! favourable one.  Later care should be taken that exiting higher composition
! sets are not set equal to the stable
! pph   number of phases for which a grid has been calculated
! nrel  number of components
! nphl(p) is last gridpoint for phase(p), nphl(0)=0, p=1,pph
! iphx(p) phase number of phase(p) (skipping dormant and suspended phases)
! xarr(1..nrel,i)  composition of gridpoint i
! garr(i)  Gibbs energy/RT for gridpoint i
! nr    is the number of stable phases in the solution
! iphl(s) the phase number of the stable phases s (not ordered)
! cmu   are the chemical potentials/RT of the solution
! ceq   equilibrium record
! called by global_gridmin
   implicit none
   integer pph,nrel,nr
   integer, dimension(0:*) :: nphl
   integer, dimension(*) :: iphl,iphx
   double precision, dimension(*) :: cmu
   real garr(*),xarr(nrel,*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ig1,ign,ip,iph,ics,jph,lokcs,lokph,mode,ny,ie,ig,kp,i,zph
   double precision yarr(maxconst),qq(5),xxx,dgmin,gmax
   real dg,gplan
!   write(*,*)'3Y In set_metastable constitution'
!   goto 1000
! The phases that have gridpoints calculated are in iphx(1..pph)
   phloop: do zph=1,pph
      iph=iphx(zph)
      do jph=1,nr
         if(iph.eq.iphl(jph)) then
! this phase is stable, skip
            cycle phloop
         endif
      enddo
! this phase is metastable, find its gridpoint closesed to the tangent plane
! the grid points belonging to phase iph is between nphl(zph-1) and nphl(zph)
! NOTE nphl(0)=0
      ig1=nphl(zph-1)+1
      ign=nphl(zph)
! if ign=ig1 there is a single gridpoint, otherwise seach for minimim
      dgmin=-1.0d12
      ip=0
! search for gripoint closeset to stable plane defined by cmu
      do ig=ig1,ign
         gplan=zero
         do ie=1,nrel
            gplan=gplan+xarr(ie,ig)*cmu(ie)
         enddo
         dg=gplan-garr(ig)
         if(abs(dg).lt.abs(dgmin)) then
            ip=ig
            dgmin=dg
         endif
      enddo
!      write(*,79)'3Y Least unstable gridpoint: ',zph,iph,ig1,ign,ip,dgmin
79    format(a,5(i6,1x),1pe12.4)
      if(ign.gt.ig1) then
! retrieve constitution for this gridpoint and insert it in phase
! must provide mode and iph. The subroutine returns ny and yarr
! mode is the gridpoint in the phase
         mode=ip-nphl(zph)
! find the constitution of this gridpoint
!      call generate_grid(mode,iph,ign,nrel,xarr,garr,ny,yarr,gmax,ceq)
!         write(*,*)'3Y Get constitution of metastable phase ',iph,mode
         if(mode.gt.0) then
            call generic_grid_generator(mode,iph,ign,nrel,xarr,garr,&
                 ny,yarr,gmax,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call set_constitution(iph,1,yarr,qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
         endif
      endif
! set driving force also for phases with fix composition
      call set_driving_force(iph,1,dgmin,ceq)
500   continue
   enddo phloop
1000 continue
   return
 end subroutine set_metastable_constitutions2
 
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 logical function global_equil_check1(mode,ceq)
! subroutine global_equil_check(ceq,newceq)
!
! THIS MUST BE REWRITTEN
!
! This subroutine checks there are any gridpoints below the calculated solution
! if not it is taken as a correct global equilibrium
! This avoids creating any new composition sets but may fail in some cases
! to detect that the equilibrium is not global.
! mode=1 means try to recalculate equilibrium if not global (not implemented)
   implicit none
   integer mode
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   TYPE(gtp_equilibrium_data), target :: cceq
   TYPE(gtp_equilibrium_data), pointer :: pceq
   logical global,newcs,notglobwarning1,notglobwarning2,wrongfrac
   real, allocatable :: xarr(:,:),garr(:)
   real sumx
   double precision, dimension(maxconst) :: yarr
   double precision totmol,totmass,amount,gmax
   integer, allocatable :: kphl(:),iphx(:)
   integer gmode,iph,ngg,nrel,ny,ifri,firstpoint,sumng,nrph,ii,jj,nz,lokcs
   integer, parameter :: maxgrid=400000
!
!   write(*,*)'3Y In global_equil_check1',mode
   global=.TRUE.
   notglobwarning1=.TRUE.
   notglobwarning2=.TRUE.
! COPY the whole equilibrium record to avoid destroying anything!!
! otherwise I had strange problems with amounts of phases ??
   cceq=ceq
   pceq=>cceq
   nrph=noofph
   allocate(kphl(0:nrph+1))
   allocate(iphx(nrph+1))
!
   sumng=0
   ifri=0
   firstpoint=1
   iphx=0
   kphl=0
   ifri=0
   wrongfrac=.true.
   loop1: do iph=1,nrph
! this loop just to calculate how many points will be generated for allocation
      if(test_phase_status(iph,1,amount,pceq).le.PHDORM) cycle loop1
      ifri=ifri+1
      iphx(ifri)=iph
!      kphl(ifri)=firstpoint
! not used: nrel, xarr, garr, yarr, gmax
!      call generate_grid(-1,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,pceq)
!      if(gx%bmperr.ne.0) goto 1000
!      kphl(ifri)=ngg
!      sumng=sumng+ngg
   enddo loop1
!   write(*,*)'3Y Generated gridpoints for check:  ',sumng
!
   nrel=noofel
! allocate arrays, added 1 to avoid a segmenentation fault ....
   allocate(xarr(nrel,maxgrid))
   allocate(garr(maxgrid))
! calculate the composition and G for the gridpoints
   ii=1
   loop2: do ifri=1,nrph
      iph=iphx(ifri)
      if(iph.gt.0) then
!         write(*,*)'3y GG: ',ifri,iph,ii
! generate_grid restore original composition before ending
!         call generate_grid(0,iph,ngg,nrel,xarr(1,ii),garr(ii),&
!              ny,yarr,gmax,pceq)
         ngg=maxgrid-ii
!         write(*,*)'3Y calling generic_grid 2: ',iph,ngg
         call generic_grid_generator(0,iph,ngg,nrel,xarr(1,ii),garr(ii),&
              ny,yarr,gmax,pceq)
         if(gx%bmperr.ne.0) goto 1000
         kphl(ifri)=kphl(ifri-1)+ngg
         ii=kphl(ifri)+1
         sumng=kphl(ifri)
      endif
   enddo loop2
!   write(*,*)'3Y Calculated ',sumng,' gridpoints for check.'
! check if any gridpoint is below the G surface defined by cmuval
   iph=0
   ny=0
   loop4: do ifri=1,sumng
! keep track of the phase the gridpoint belongs to
      if(ifri.gt.ny) then
         iph=iph+1
         ny=ny+kphl(iph)
!         write(*,*)'3Y Incrementing iph: ',iph,kphl(iph),ifri,ny
      endif
      gmax=zero
      sumx=0.0E0
      do ngg=1,nrel
         gmax=gmax+dble(xarr(ngg,ifri))*pceq%cmuval(ngg)
         sumx=sumx+xarr(ngg,ifri)
      enddo
!      if(ifri.eq.sumng) write(*,*)'3Y OK ',ifri,iph
      if(abs(sumx-1.0E0).gt.1.0E-4) then
! skip this warning message ...
!         if(wrongfrac) then
!            write(*,69)'3Y Some gridpoint fractions wrong',&
!                 iph,ifri,sumx,garr(ifri)
!69          format(a,i3,i7,2(1pe14.6))
!            wrongfrac=.false.
!         endif
         cycle loop4
      endif
!      write(*,75)'3Y check: ',ifri,iph,garr(ifri),gmax,garr(ifri)-gmax
75    format(a,i6,i4,5(1pe12.4))
      if(dble(garr(ifri))-gmax.lt.-1.0D-7*abs(gmax)) then
!      if(dble(garr(ifri))-gmax.lt.-1.0D-7) then
! if the phase is stoichiometric and stable this is no error
! find the phase record using the phase tuple
         lokcs=phasetuple(iph)%lokvares
         if(ceq%phase_varres(lokcs)%phstate.ge.PHENTSTAB) then
! if number of constituent fractions equal to sublattice the composition is fix
            nz=size(ceq%phase_varres(lokcs)%sites)-&
                 size(ceq%phase_varres(lokcs)%yfr)
            if(nz.eq.0) then
!               write(*,*)'3Y No problem, phase stable with fix composition'
               cycle loop4
            else
! This gridpoint is for a solution phase that is stable and has a grid point
! below current equilibrium plane.  Could be the cubic carbide/austenit case 
               global=.FALSE.
               if(notglobwarning2) then
                  write(kou,87)iph,(xarr(ngg,ifri),ngg=1,nrel)
                  gx%bmperr=4352
                  notglobwarning2=.FALSE.
               endif
            endif
         else
! This gridpoint is for a phase that is not stable but has a grid point
! below the current equilibrium plane 
            global=.FALSE.
            if(notglobwarning1) then
               write(kou,87)ceq%phase_varres(lokcs)%phtupx,&
                    (xarr(ngg,ifri),ngg=1,nrel)
87             format(/' *** Equilibrium not global, phasetuple ',&
                    i3,' is stable with mole fractions:'/(2x,9F8.5))
               gx%bmperr=4352
               notglobwarning1=.FALSE.
            endif
         endif
      endif
!      if(ifri.eq.sumng) write(*,*)'OK ',ifri
   enddo loop4
! no gridpoint below current G surface
!   write(*,*)'3Y Still OK?'
   goto 1000
! Found gridpoint below gmax, if mode=/=1 just return error message
500 continue
!   write(*,*)'3Y Sorry I have not yet implemented automatic recalculation!'
   if(mode.eq.1) then
! Here you should implement automatic recalculation including new phase ...
      continue
!      gx%bmperr=0
!   else
!      write(*,*)'3Y Please include this phase and recalculate equilibrium'
   endif
!
1000 continue
!   write(*,*)'3Y Deallocating, check due to segmentation fault ...'
   if(allocated(xarr)) then
      deallocate(xarr)
      deallocate(garr)
      deallocate(kphl)
      deallocate(iphx)
   endif
1010 continue
   global_equil_check1=global
   return
 end function global_equil_check1

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine gridmin_check(nystph,kp,nrel,xarr,garr,xknown,ngrid,pph,&
      cmu,yphl,iphx,ceq)
! This subroutine checks if a calculated solution is correct by  
! checking if there are any gridpoints below the surface defined
! by the chemical potentials cmu
! The gridpoints have already been calculated
! nystph return 0 or the phase tuple index for a new stable phase
! kp are the number of gridpoints, nrel is number of components
! composition of each gridpoint in xarr, G in garr
! xknown is the known overall composition
! ngrid is last calculated gridpoint point for a phase jj
! pph is number of phases for which there is a grid
! iphx is phase numbers
! cmu are the final chemical potentials
! yphl is just needed as a dummy
! ceq is current equilibrium record
   implicit none
   integer kp,nrel,jp,ie,mode,pph,nystph
   double precision, parameter :: phfmin=1.0D-8
   real xarr(nrel,*),garr(*)
   double precision xknown(*)
   integer, dimension(*) :: ngrid,iphx
   double precision cmu(*),gsurf,gstable,gd,yphl(*),qq(5),rtn,gdmin
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokph,lokcs,zph,ibias,ics,iph,ny
! setting a value of addph forces that gridpoint to be added, used for test
!   integer :: addph=8 
!   integer :: addph=100
   integer :: addph=0
   double precision gmax
!   integer idum(1000)
   save addph
!
!   write(*,*)'3Y Entering gridmin_check',addph
   gstable=zero
   nystph=0
   mode=0
   rtn=globaldata%rgas*ceq%tpval(1)
   gdmin=-1.0D5
   do jp=1,kp
      gsurf=zero
      do ie=1,nrel
         gsurf=gsurf+xarr(ie,jp)*cmu(ie)
      enddo
      gsurf=gsurf/rtn
! If garr(jp) more negative than gsurf (gd>0) this gridpoint is stable
! mixing real and double precision is not a numerical problem here
      gd=gsurf-garr(jp)
      if(gd.gt.gdmin) gdmin=gd
!      write(*,17)'3Y grid comparison: ',gd,garr(jp),gsurf
!17    format(a,3(1pe12.4))
      if(gd.gt.gstable) then
! this gridpoint should be set as stable and equilibriuum recalculated
         gstable=gd
         mode=jp
      endif
   enddo
! if mode nonzero there is a gridpoint below the calculated surface
! just for test using the FeOUZr case set mode=25, that should be in liquid
! just for test using the FeOUZr case set mode=7, that should be O2 in gas
   mode=addph
   addph=0
   if(mode.gt.0) then
! we have to find which phase it is, this strange loop should find that
      ibias=0
      do zph=1,pph
! ngrid(zph) is the first gridpoints of phase zph
         if(mode.le.ngrid(zph)) then
            mode=mode-ibias
            goto 115
         else
            ibias=ngrid(zph)
         endif
      enddo
115   continue
!      write(*,*)'3Y mode, ibias and phase: ',mode,ibias,iphx(zph)
! This call return constitution of new gridpoint
!      call generate_grid(mode,iphx(zph),ibias,nrel,xarr,garr,ny,yphl,gmax,ceq)
! ibias, xarr, garr not uses
         call generic_grid_generator(mode,iphx(zph),ibias,nrel,&
              xarr,garr,ny,yphl,gmax,ceq)
      if(gx%bmperr.ne.0) goto 1000
      iph=iphx(zph)
      write(*,*)'3Y Gridmin check found new stable phase: ',iph
!-------------------
! new stable phase is iph, constituent fractions in yphl
! check if compset 1 of phase is already stable, if so maybe create compset
      lokph=phases(iph)
      lokcs=phlista(lokph)%linktocs(1)
      ics=1
200   continue
      if(ceq%phase_varres(lokcs)%dgm.lt.zero) then
! this composition set not stable, set it as stable with fractions yphl
! finetunig needed here ....
         ceq%phase_varres(lokcs)%dgm=zero
! strange, new calculation failed with small amount but worked with a lot ...
!         ceq%phase_varres(lokcs)%amount(1)=one
         ceq%phase_varres(lokcs)%amfu=one
         write(*,222)iph,ics,ny,(yphl(ie),ie=1,ny)
222      format('3Y added: ',3i3,10(f6.3))
         call set_constitution(iph,ics,yphl,qq,ceq)
      else
         ics=ics+1
         if(ics.gt.phlista(lokph)%noofcs) then
! create new composition set if allowed
            if(btest(globaldata%status,GSNOACS)) then
               gx%bmperr=4177; goto 1000
            endif
            write(*,*)'3Y Creating new composition set for ',iph
            call enter_composition_set(iph,'    ','AUTO',ics)
            if(gx%bmperr.ne.0) goto 1000
! link to new compositiin set stored here
! set a negative zero driving force
            lokcs=phlista(lokph)%linktocs(ics)
            ceq%phase_varres(lokcs)%dgm=-one
         else
            lokcs=phlista(lokph)%linktocs(ics)
         endif
! jump back to label 200 to test if this composition set is free
         goto 200
      endif
      nystph=10*iph+ics
   else
! no new phase found, just to see some values
      write(*,*)'3Y DG at least unstable gridpoint: ',gdmin
   endif
1000 continue
   return
 end subroutine gridmin_check

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine separate_constitutions(ceq)
! Go through all entered phases and if there are two composition sets
! that have similar constitutions then separate them
! Used during mapping of for example Fe-Cr to detect the miscibility gap
    implicit none
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer phtup,nextset,lokcs1,lokcs2,ic,ll,lokph,ss
    double precision, allocatable :: yarr(:)
    allph: do phtup=1,nooftup()
       nextset=phasetuple(phtup)%nextcs
       if(nextset.eq.0) cycle allph
       lokcs1=phasetuple(phtup)%lokvares
       lokcs2=phasetuple(nextset)%lokvares
       do ic=1,size(ceq%phase_varres(lokcs1)%yfr)
          if(abs(ceq%phase_varres(lokcs1)%yfr(ic)-&
               ceq%phase_varres(lokcs2)%yfr(ic)).gt.1.0D-2) then
!             write(*,66)'3Y two compsets not same:',lokcs1,lokcs2,ceq%tpval(1)
             cycle allph
          endif
       enddo
! These two composition sets have identical compositions, skip if both stable
!       write(*,66)'3Y two compsets same:',lokcs1,lokcs2,ceq%tpval(1)
66     format(a,2i3,f10.2)
       if(ceq%phase_varres(lokcs1)%phstate.ge.PHENTSTAB) then
          if(ceq%phase_varres(lokcs2)%phstate.ge.PHENTSTAB) cycle allph
! set the constitution of lokcs2 to one-the stable
! or maybe to its default??
          lokph=phasetuple(phtup)%lokph
          ic=0
          phsubl1: do ll=1,phlista(lokph)%noofsubl
             if(phlista(lokph)%nooffr(ll).eq.1) then
                ic=ic+1; cycle phsubl1
             endif
             do ss=1,phlista(lokph)%nooffr(ll)
                ic=ic+1
!                if(ceq%phase_varres(lokcs1)%yfr(ic).gt.5.0D-1) then
!                   ceq%phase_varres(lokcs2)%yfr(ic)=1.0D-1
!                else
!                   ceq%phase_varres(lokcs2)%yfr(ic)=9.0D-1
!                endif
                ceq%phase_varres(lokcs2)%yfr(ic)=&
                     one-ceq%phase_varres(lokcs1)%yfr(ic)
             enddo
          enddo phsubl1
!          write(*,77)'3Y: nyy:',ic,(ceq%phase_varres(lokcs2)%yfr(ss),ss=1,ic)
77        format(a,i3,8F6.3)
       elseif(ceq%phase_varres(lokcs2)%phstate.ge.PHENTSTAB) then
! set the constitution of lokcs2 to one-the stable
! or maybe to its default??
          lokph=phasetuple(phtup)%lokph
          ic=0
          phsubl2: do ll=1,phlista(lokph)%noofsubl
             if(phlista(lokph)%nooffr(ll).eq.1) then
                ic=ic+1; cycle phsubl2
             endif
             do ss=1,phlista(lokph)%nooffr(ll)
                ic=ic+1
!                if(ceq%phase_varres(lokcs2)%yfr(ic).gt.5.0D-1) then
!                   ceq%phase_varres(lokcs1)%yfr(ic)=1.0D-1
!                else
!                   ceq%phase_varres(lokcs1)%yfr(ic)=9.0D-1
!                endif
                ceq%phase_varres(lokcs1)%yfr(ic)=&
                     one-ceq%phase_varres(lokcs2)%yfr(ic)
             enddo
          enddo phsubl2
!          write(*,77)'3Y: nyy:',ic,(ceq%phase_varres(lokcs1)%yfr(ss),ss=1,ic)
!       else
!          write(*,*)'Both compsets unstable'
       endif
    enddo allph
1000 continue
  end subroutine separate_constitutions

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>      17. Miscellaneous
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 integer function phvarlok(lokph)
! return index of the first phase_varres record for phase with location lokph
! needed for external routines as phlista is private
   implicit none
   integer lokph
!\end{verbatim}
   phvarlok=phlista(lokph)%linktocs(1)
   return
 end function phvarlok

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine palmtree(lokph)
! Initiates a numbering of all interaction trees of an endmember of a phase
   implicit none
   integer lokph
!\end{verbatim}
   integer seq,level
   type(gtp_endmember), pointer :: endm
   type(gtp_interaction), pointer :: intrec
   type stack
      type(gtp_interaction), pointer :: p1
   end type stack
   type(stack), dimension(5) :: int_stack
   logical both
   both=.false.
   endm=>phlista(lokph)%ordered
70 continue
   emloop:do while(associated(endm))
      intrec=>endm%intpointer
      seq=0
      level=0
100   continue
      do while(associated(intrec))
         level=level+1
         if(level.gt.5) then
            write(*,*)'3Y Interaction more than 5 levels deep!'
            gx%bmperr=4347; goto 1000
         endif
         int_stack(level)%p1=>intrec
         seq=seq+1
         intrec%order=seq
         intrec=>intrec%highlink
      enddo
      if(level.gt.0) then
         intrec=>int_stack(level)%p1
         level=level-1
         intrec=>intrec%nextlink
         goto 100
      endif
      endm=>endm%nextem
   enddo emloop
   if(.not.both .and. associated(phlista(lokph)%disordered)) then
      endm=>phlista(lokph)%disordered
      both=.true.
      goto 70
   endif
1000 continue
   return
 end subroutine palmtree

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 logical function allowenter(mode)
! Check if certain commands are allowed
! mode=1 means entering an element or species
! mode=2 means entering a phase
! mode=3 means entering an equilibrium
! returns TRUE if command can be executed
   implicit none
   integer mode
!\end{verbatim}
!   write(*,*)'3Y In allowenter: ',mode
   logical yesorno
   yesorno=.FALSE.
   if(mode.le.0 .or. mode.gt.3) goto 1000
   if(mode.eq.1) then
! enter element or species not allowed after entering first phase
      if(noofph.gt.0) goto 1000
      yesorno=.TRUE.
   elseif(mode.eq.2) then
! enter phases of a disordred fraction set not allowed
! if there are no elements or after entering a second equilibrium
!      write(*,*)'3Y allowenter ',mode,noofel,eqfree,noofph
      if(noofel.eq.0) goto 1000
      if(eqfree.gt.2) goto 1000
      yesorno=.TRUE.
   elseif(mode.eq.3) then
! there must be at lease one phase before entering a second equilibrium
! Note this is tested also for entering the default equilibrium
!      write(*,*)'3Y mode 3: ',eqfree,noofph
      if(eqfree.ge.2 .and. noofph.eq.0) goto 1000
      yesorno=.TRUE.
   endif
1000 continue
   allowenter=yesorno
!   write(*,*)'3Y: allowenter:',yesorno,mode
   return
 end function allowenter

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 logical function proper_symbol_name(name,typ)
! checks that name is a proper name for a symbol
! A proper name must start with a letter A-Z
! for typ=0 it must contain only letters, digits and underscore
! for typ=1 it may contain also +, - maybe ?
! It must not be equal to a state variable
   implicit none
   integer typ
   character name*(*)
!\end{verbatim}
   character name2*64,ch1*1,chx*1
   integer jl
   logical korrekt
!   write(*,*)'3Y entering proper_symbol_name: ',name,typ
   korrekt=.FALSE.
   if(typ.lt.0 .or. typ.gt.0) then
      gx%bmperr=4139; goto 1000
   endif
   name2=name
   call capson(name2)
   if(.not.ucletter(name2(1:1))) then
! the first character of a symbol must always be a letter A-Z
!      write(*,*)'3Y Wrong first letter of symbol: ',name2(1:1),':',name2(1:5)
      gx%bmperr=4137; goto 1000
   endif
   jl=1
!   write(*,*)'3Y check name: ',name2
100 continue
   jl=jl+1
   ch1=name2(jl:jl)
! always finish when fining a space 
   if(ch1.eq.' ') then
! any symbol with at least 3 characters OK
      if(jl.le.2) then
! A single letter must be a state variable
         if(name2(1:1).eq.'A' .or. name2(1:1).eq.'B' .or. &
              name2(1:1).eq.'G' .or. name2(1:1).eq.'H' .or. &
              name2(1:1).eq.'M' .or. name2(1:1).eq.'N' .or.&
              name2(1:1).eq.'P' .or. name2(1:1).eq.'Q' .or. &
              name2(1:1).eq.'S' .or. &
              name2(1:1).eq.'T' .or. name2(1:1).eq.'U' .or.&
              name2(1:1).eq.'V' .or. name2(1:1).eq.'W' .or. &
              name2(1:1).eq.'X' .or. name2(1:1).eq.'Y') then
            if(jl.eq.2) then
               gx%bmperr=4137; goto 1000
            elseif(name2(2:2).eq.'F' .or. name2(2:2).eq.'M' .or. &
                 name2(2:2).eq.'P' .or. name2(2:2).eq.'U' .or. &
                 name2(2:2).eq.'V' .or. name2(2:2).eq.'W') then
! A two letter name must not have certain second letter
               gx%bmperr=4137; goto 1000
            endif
         endif
      endif
      korrekt=.TRUE.
      name(jl:)=' '
      goto 1000
   endif
   if(typ.eq.0) then
      if(ch1.ge.'0' .and. ch1.le.'9') goto 100
      if(ch1.ge.'A' .and. ch1.le.'Z') goto 100
      if(ch1.eq.'_') goto 100
      gx%bmperr=4138
!   else
! unknown type of symbol
!      gx%bmperr=4139
   endif
!
1000 continue
!
   if(.not.korrekt) write(*,*)'3Y Illegal name: ',name2,jl
   proper_symbol_name=korrekt
   return
 end function proper_symbol_name

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!-\begin{verbatim}
 subroutine list_free_lists(lut)
! for debugging the free lists and routines using them
   implicit none
   integer lut
!-\end{verbatim}
   integer lok,last
   write(lut,1007)noofel,noofsp,noofph,noofem,noofint,noofprop,&
        notpf(),highcs,eqfree-1,nsvfun,reffree-1,addrecs
1007 format('Records for elements, species, phases:           ',3i5/&
            'end members, interactions, properties:           ',3i5/&
            'TP-funs, composition sets, equilibria:           ',3i5/&
            'state variable functions, references, additions: ',3i5)
!----------------------------
! first free is csfree, free list is only in equilibrium firsteq
600 continue
   write(lut,610)csfree,highcs
610 format('Phase_varres first free/highcs: ',2i5)
! NOTE csfree can be higher than highcs ... after deletion pointers can go back
! UNFINISHED??
!   lok=csfree
! list free list for composition sets
!   write(*,*)'3Y csfree and highcs: ',csfree,highcs
!611 continue
!   last=lok
!   lok=firsteq%phase_varres(last)%nextfree
!   write(*,*)'3Y lok: ',last,lok
!   if(lok.gt.0) goto 611
!
   lok=csfree
620 continue
   if(lok+5.lt.highcs) then
      last=lok
      lok=firsteq%phase_varres(last)%nextfree
      write(*,*)'3Y free varres record at: ',lok,last
      if(lok.le.0) then
         write(lut,*)'Error in phase_varres free list',last,lok
         goto 1000
      else
         goto 620
      endif
   endif
! no more
630 continue
1000 continue
   return
 end subroutine list_free_lists

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine set_phase_amounts(jph,ics,val,ceq)
! set the amount formula units of a phase. Called from user i/f
! iph can be -1 meaning all phases, all composition sets
   implicit none
   integer jph,ics
   double precision val
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer iph,lokph,lokcs
   double precision amount
   if(jph.lt.0) then
      iph=1; ics=1
   else
      iph=jph
   endif
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
100 continue
   if(test_phase_status(iph,ics,amount,ceq).gt.3) goto 700
!   ceq%phase_varres(lokcs)%amount(1)=val
   ceq%phase_varres(lokcs)%amfu=val
700 continue
   if(jph.lt.0) then
      ics=ics+1
710   continue
      call get_phase_compset(iph,ics,lokph,lokcs)
      if(gx%bmperr.ne.0) then
         gx%bmperr=0;
         iph=iph+1
         if(iph.gt.noofph) goto 1000
         ics=1; goto 710
      endif
      goto 100
   endif
1000 continue
   return
 end subroutine set_phase_amounts

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine set_default_constitution(iph,ics,ceq)
! the current constitution of (iph#ics) is set to its default constitution
!  (if any), otherwise a random value.  The amount of the phase not changed
   implicit none
   integer iph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokph,lokcs,ll,jj,kk,kk0
   type(gtp_phase_varres), pointer :: cset
   double precision, allocatable :: yarr(:)
   double precision sum, qq(5),var
!
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   cset=>ceq%phase_varres(lokcs)
! we must use set_constitution at the end to update various internal variables
   allocate(yarr(phlista(lokph)%tnooffr))
   if(btest(cset%status2,CSDEFCON)) then
! there is a preset default constitution
!      write(*,12)'3Y mmfy: ',(cset%mmyfr(kk),kk=1,phlista(lokph)%tnooffr)
      kk=0
      subl1: do ll=1,phlista(lokph)%noofsubl
         kk0=kk
         sum=zero
         if(phlista(lokph)%nooffr(ll).gt.1) then
            do jj=1,phlista(lokph)%nooffr(ll)
! negative mmy(kk) means < , a maximum, set a small value
               kk=kk+1
               if(cset%mmyfr(kk).lt.0.0E0) then
                  yarr(kk)=0.01D0
               else
                  yarr(kk)=one
               endif
               sum=sum+yarr(kk)
            enddo
            kk=kk0
! the sum of fractions should be unity, hm done in set_constitution also ...
            do jj=1,phlista(lokph)%nooffr(ll)
               kk=kk+1
               yarr(kk)=yarr(kk)/sum
            enddo
         else
! a single constituent, just increment kk and leave fraction as unity
            kk=kk+1
            yarr(kk)=one
         endif
      enddo subl1
!      write(*,12)'3Y defy: ',(yarr(kk),kk=1,phlista(lokph)%tnooffr)
12    format(a,10F6.3)
   else
! there is no default constitution, set equal amount of all fractions
! with some randomness
!      write(*,*)'3Y No default constituition for: ',iph,ics
      kk=0
      subl2: do ll=1,phlista(lokph)%noofsubl
         if(phlista(lokph)%nooffr(ll).gt.1) then
! set equal amount of all fractions with some variation
            sum=one/real(phlista(lokph)%nooffr(ll))
            var=0.1D0*sum
            do jj=1,phlista(lokph)%nooffr(ll)
               kk=kk+1
               yarr(kk)=sum+var
               var=-0.9D0*var
            enddo
         else
! a single constituent, just increment kk and ensure the fraction is unity
            kk=kk+1
            yarr(kk)=one
         endif
      enddo subl2
   endif
!   write(*,411)yarr
411 format('3Y set_def_const: ',8F7.4,(10f7.4))
! in this routine the fractions in each sublattice is normallized to be unity
   call set_constitution(iph,ics,yarr,qq,ceq)
   deallocate(yarr)
1000 continue
   return
 end subroutine set_default_constitution

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine todo_before(mode,ceq)
! this could be called before an equilibrium calculation
! It should remove any phase amounts and clears CSSTABLE
! DUMMY
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode
!\end{verbatim}
   integer iph,ics,lokph,lokcs
!
!   write(*,*)'3Y Todo_before ... not implemented'
   goto 1000
!
   phloop: do iph=1,noph()
      lokph=phases(iph)
! skip hidden phases
      if(btest(phlista(lokph)%status1,PHHID)) cycle
300      csloop: do ics=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(ics)
!         ceq%phase_varres(lokcs)%amount(1)=zero
         ceq%phase_varres(lokcs)%amfu=zero
!         ceq%phase_varres(lokcs)%status2=&
!              ibclr(ceq%phase_varres(lokcs)%status2,CSSTABLE)
      enddo csloop
   enddo phloop
!
1000 continue
   return
 end subroutine todo_before

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine todo_after_found_equilibrium(mode,ceq)
! this is called after an equilibrium calculation by calceq2 and calceq3
! It marks stable phase (set CSSTABLE and remove any CSAUTO)
! It removes redundant unstable composition sets created automatically
! (CSAUTO set).  It will also shift stable composition sets to loweest 
! possible (it will take into account if there are default constituent 
! fractions, CSDEFCON set).
! mode determine some of the actions, at present only >0 or <0 matters
!
! >>>>>>>>>>> THIS IS DANGEROUS IN PARALLEL PROCESSING
! It should work in step and map as a composition set that once been stable
! will never be removed except if one does global minimization during the
! step and map. The function global_equil_check works on a copy of the
! ceq record and creates only a grid, it does not create any composition sets.
! Automatically entered metallic-FCC and MC-carbides may shift composition sets.
! Such shifts can be avoided by manual entering composition sets with
! default constitutions, but comparing a stable constitution with a
! default is not trivial ...
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode
!\end{verbatim}
   integer iph,ics,lokph,lokics,jcs,lokjcs,lastset,lokkcs,kzz,jtup,qq
   integer jstat2,fit,phs,haha1,haha2,disfravares
   double precision val,xj1,xj2
   logical notok,noremove,globalok,once
   character jpre*4,jsuf*4
   real, dimension(:), allocatable :: tmmyfr
!
!   write(*,*)'3Y in todo_after',mode
!----------------------------------------------------------------
   if(btest(globaldata%status,GSNOAFTEREQ)) goto 1000
   nostart: if(mode.lt.0 .or. btest(globaldata%status,GSTGRID)) then
! if mode<0 the conditions did not allow gridmin before use it after
! Problems with this calculation so global_equil_check is disabled inside ...
      write(*,3)
3     format('Testing if any gridpoint is below the calculated equilibrium')
      if(btest(globaldata%status,GSNOTELCOMP)) then
         write(*,*)'3Y Cannot test global equilibrium when these components'
         goto 1000
      endif
      qq=1
      globalok=global_equil_check1(qq,ceq)
!      write(*,*)'3Y Back from GEC1'
      if(globalok) then
! if TRUE equilibrium OK or it could not be tested
         if(gx%bmperr.ne.0) then
            write(*,*)'3Y Testing equilibrium with gridminimizer failed'
            goto 1000
         endif
!         write(*,*)'3Y Grid minimizer test of equilibrium OK'
      else
! if FALSE the test showed this is not a global equilibrium
         write(*,*)'3Y Grid minimizer found equilibrium wrong',gx%bmperr
         goto 1000
      endif
   endif nostart
!--------------------------------------------------------------------
! Shift all stable composition down to lower comp.sets
200 continue
!   write(*,*)'3Y Shifting composition sets'
   phloop1: do iph=1,noph()
      lokph=phases(iph)
      if(btest(phlista(lokph)%status1,PHHID)) cycle
      csloop1: do ics=2,phlista(lokph)%noofcs
         lokics=phlista(lokph)%linktocs(ics)
!         write(*,*)'3Y shift down: ',ics,lokics,&
!              ceq%phase_varres(lokics)%phstate,&
!              btest(ceq%phase_varres(lokics)%status2,CSAUTO),&
!              btest(ceq%phase_varres(lokics)%status2,CSTEMPAR)
         if(ceq%phase_varres(lokics)%phstate.eq.PHENTSTAB .and. &
              btest(ceq%phase_varres(lokics)%status2,CSTEMPAR)) then
!              btest(ceq%phase_varres(lokics)%status2,CSAUTO)) then
            fit=100
! This comp.set is stable, check if a lower compset is unstable
            csloop2: do jcs=1,ics-1
               lokjcs=phlista(lokph)%linktocs(jcs)
               if(ceq%phase_varres(lokjcs)%phstate.le.PHENTERED) then
! do not bother if composition of lokics fits defaults in lokjcs
!                  if(.not.checkdefcon(lokics,lokjcs,fit,ceq)) cycle csloop2
!                  write(*,*)'3Y Moving stable comp.set ',ics,' down to ',jcs
                  goto 500
               elseif(jcs.eq.ics-1) then
                  if(fit.gt.2) then
! No lower unstable comp.set, or no one which almost fit default const,
! lokics must remain stable, remove CSAUTO bit
! Do not remove the suffix _AUTO
!                     write(*,*)'3Y Keeping AUTO comp.set ',ics,lokics
                     ceq%phase_varres(lokics)%status2=&
                          ibclr(ceq%phase_varres(lokics)%status2,CSAUTO)
                     exit csloop2
                  endif
               else
                  cycle csloop2
               endif
! Accept a default consitution which almost fits the default
!               write(*,*)'3Y Accept fit to default: ',fit,lokics,lokjcs
500            continue
!               write(*,*)'3Y Move stable to lower unstable compsets'
! move STABLE lokics to UNSTABLE lokjcs
!               switch_compsets(... ???
! save some jcs values of amount, dgm, status, pre&suffix and tuple index
               xj1=ceq%phase_varres(lokjcs)%amfu
               xj2=ceq%phase_varres(lokjcs)%dgm
               jtup=ceq%phase_varres(lokjcs)%phtupx
               jstat2=ceq%phase_varres(lokjcs)%status2
               jpre=ceq%phase_varres(lokjcs)%prefix
               jsuf=ceq%phase_varres(lokjcs)%suffix
               phs=ceq%phase_varres(lokjcs)%phstate
!               write(*,489)lokics,lokjcs
               if(ocv()) write(*,489)ceq%phase_varres(lokics)%phtupx,jtup
489            format('3Y move results from tuplet ',i4,' to ',i4)
!                  write(*,501)lokics,ceq%phase_varres(lokics)%mmyfr
!                  write(*,501)lokjcs,ceq%phase_varres(lokjcs)%mmyfr
501            format('3Y 501: ',i5,10F5.1)
! copy main content of the phase_varres(lokics) record to phase_varres(lokjcs)
! BEWARE mmyfr must be kept!
! BEWARE disordered fraction set!!!!
                  disfravares=ceq%phase_varres(lokjcs)%disfra%varreslink
                  if(allocated(ceq%phase_varres(lokjcs)%mmyfr)) then
                     allocate(tmmyfr(size(ceq%phase_varres(lokjcs)%mmyfr)))
                     tmmyfr=ceq%phase_varres(lokjcs)%mmyfr
                     ceq%phase_varres(lokjcs)=ceq%phase_varres(lokics)
                     ceq%phase_varres(lokics)%mmyfr=tmmyfr
                     deallocate(tmmyfr)
                  endif
                  ceq%phase_varres(lokjcs)=ceq%phase_varres(lokics)
! Some content in jcs must be set or restorted separately
                  ceq%phase_varres(lokjcs)%phtupx=jtup
                  ceq%phase_varres(lokjcs)%status2=jstat2
                  ceq%phase_varres(lokjcs)%prefix=jpre
                  ceq%phase_varres(lokjcs)%suffix=jsuf
                  ceq%phase_varres(lokjcs)%phstate=PHENTSTAB
!                  ceq%phase_varres(lokjcs)%status2=&
!                       ibset(ceq%phase_varres(lokjcs)%status2,CSSTABLE)
!                  write(*,501)lokics,ceq%phase_varres(lokics)%mmyfr
!                  write(*,501)lokjcs,ceq%phase_varres(lokjcs)%mmyfr
! maybe CSAUTO bit set, always remove it!
!                  write(*,*)'3Y Ensure CSAUTO cleared in ',jcs
                  ceq%phase_varres(lokjcs)%status2=&
                       ibclr(ceq%phase_varres(lokjcs)%status2,CSAUTO)
! Some content in ics must be set separately from saved values of jcs
                  ceq%phase_varres(lokics)%amfu=xj1
                  ceq%phase_varres(lokics)%dgm=xj2
                  ceq%phase_varres(lokics)%phstate=phs
! clear the stable bit and set AUTO of ics ??
!                  ceq%phase_varres(lokics)%status2=&
!                       ibclr(ceq%phase_varres(lokics)%status2,CSSTABLE)
!                  if(btest(ceq%phase_varres(lokics)%status2,CSAUTO)) &
!                       write(*,*)'3Y AUTO bit already set in ',ics
                  ceq%phase_varres(lokics)%status2=&
                       ibset(ceq%phase_varres(lokics)%status2,CSAUTO)
! move the link to the disordered fraction set
                  ceq%phase_varres(lokjcs)%disfra%varreslink=&
                       ceq%phase_varres(lokics)%disfra%varreslink
                  ceq%phase_varres(lokics)%disfra%varreslink=disfravares
                  exit csloop2
            enddo csloop2
         endif
      enddo csloop1
   enddo phloop1
!   haha2=phlista(lokph)%linktocs(1)
!   write(*,*)'3Y mitt 1:',lokph,haha2,ceq%phase_varres(haha2)%disfra%varreslink
!   haha2=phlista(lokph)%linktocs(2)
!   if(haha2.gt.0) &
!   write(*,*)'3Y mitt 2:',lokph,haha2,ceq%phase_varres(haha2)%disfra%varreslink
! Here we may try to ensure that the stable comp.sets fits the
! default constitutions of their current set
!   write(*,*)'3Y Try to shift stable comp.sets. to match default const.'
   call shiftcompsets(ceq)
!
! upto now is safe ... now remove CSAUTO comp.sets if allowed
! check if allowed to remove
   if(btest(globaldata%status,GSNOREMCS)) goto 1000
!
! Now try to remove unstable composition sets with CSTEMPAR bit set
!   write(*,*)'3Y loop to remove comp sets and auto bits'
   once=.TRUE.
   phloop: do iph=1,noph()
      noremove=.FALSE.
      lokph=phases(iph)
      if(btest(phlista(lokph)%status1,PHHID)) cycle
! loop backwards for compsets to remove unstable with CSAUTO set
      lastset=phlista(lokph)%noofcs
      csloopdown: do ics=lastset,2,-1
         lokics=phlista(lokph)%linktocs(ics)
!         write(*,*)'3Y Checking comp.set ',ics
!         auto: if(btest(ceq%phase_varres(lokics)%status2,CSAUTO)) then
         auto: if(btest(ceq%phase_varres(lokics)%status2,CSTEMPAR)) then
            if(ceq%phase_varres(lokics)%phstate.le.PHENTERED) then
! comp.set was created automatically but is not stable, it can be removed
               if(noeq().eq.1) then
! we have just one equilibrium, OK to remove even in parallel ...
                  if(once) then
                     if(ocv()) write(*,801)lokics
801                  format('3Y Removing unstable phase tuple(s)',i5)
                     once=.FALSE.
                  endif
!                  write(*,802)'3Y removing unstable phase tuple/compset ',&
!                       ceq%phase_varres(lokics)%phtupx,lokics
802               format(a,3i5)
! remove the higherst composition set
                  call remove_composition_set(iph,.FALSE.)
                  if(gx%bmperr.ne.0) then
                     write(*,*)'3Y failed to remove tuplet:',&
                          ceq%phase_varres(lokics)%phtupx
! reset the error code but exit the attempt to clean up
                     gx%bmperr=0; goto 1000
                  endif
!                  write(*,*)'3Y Phase tuple removed for phase: ',iph
!$               elseif(omp_get_num_threads().gt.1) then
! we are running with several threads, just suspend the compset for the
! equilibrium in this thread
!$                  call suspend_composition_set(iph,.TRUE.,ceq)
               else
! when more than one equilibria in sequential eexecution suspend the compset
! in all equilibria where it is not stable
                  call suspend_composition_set(iph,.FALSE.,ceq)
               endif
            else
! the comp.set is stable, clear the CSAUTO and CSTEMPAR bits
!               write(*,*)'3Y this comp.set. should never be removed'
               ceq%phase_varres(lokics)%status2=&
                    ibclr(ceq%phase_varres(lokics)%status2,CSAUTO)
               ceq%phase_varres(lokics)%status2=&
                    ibclr(ceq%phase_varres(lokics)%status2,CSTEMPAR)
            endif
         endif auto
      enddo csloopdown
   enddo phloop
!
1000 continue
!   write(*,*)'3Y Leaving todo_after!'
!   lokph=1
!   jcs=phlista(lokph)%linktocs(1)
!   write(*,*)'3Y after 1: ',lokph,jcs,ceq%phase_varres(jcs)%disfra%varreslink
!   jcs=phlista(lokph)%linktocs(2)
!   write(*,*)'3Y Leaving todo_after'
!   if(jcs.gt.0) &
!        write(*,*)'after 2: ',lokph,jcs,ceq%phase_varres(jcs)%disfra%varreslink
   return
 end subroutine todo_after_found_equilibrium

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine checkdefcon(lokics,lokjcs,fit,ceq)
! check if composition of lokics fits default constitution in lokjcs
! return TRUE if lokics fits default in lokjcs
! NOTE lokics and lokjcs can be the same!!
   integer lokics,lokjcs,fit
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer kk
   real xdef
!   write(*,*)'3Y in checkdefcon: ',lokics,lokjcs
   if(btest(ceq%phase_varres(lokjcs)%status2,CSDEFCON)) then
!      write(*,9)(ceq%phase_varres(lokjcs)%mmyfr(fit),&
!           fit=1,size(ceq%phase_varres(lokjcs)%yfr))
9     format('3Y default: ',10F6.2)
      fit=1
      do kk=1,size(ceq%phase_varres(lokjcs)%mmyfr)
!      do kk=1,size(ceq%phase_varres(lokjcs)%yfr)
         xdef=ceq%phase_varres(lokjcs)%mmyfr(kk)
         if(xdef.eq.0) then
! no default for this constitution
            fit=fit+1
         elseif(xdef.lt.0.0) then
! A fraction with a maximum set (mmyfr<0) must be below mmyfr(kk)
            if(ceq%phase_varres(lokics)%yfr(kk).lt.abs(xdef)) fit=fit+1
!            write(*,11)ceq%phase_varres(lokics)%yfr(kk),' < ',xdef,kk,fit
11          format('3Y If ',F10.6,a,F10.6,' increment ',2i3)
         else
! A fraction with a minimum set (mmyfr>0) should be above mmyfr(kk)
            if(ceq%phase_varres(lokics)%yfr(kk).gt.xdef) fit=fit+1
!            write(*,11)ceq%phase_varres(lokics)%yfr(kk),' > ',xdef,kk,fit
         endif
      enddo
!      write(*,*)'3Y checkdefcon fit: ',fit,kk
   else
! no default constitution, perfect fit!!
      fit=size(ceq%phase_varres(lokjcs)%yfr)
   endif
1000 continue
   return
 end subroutine checkdefcon

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine shiftcompsets(ceq)
! check phase with several composition sets if they should be shifted
! to fit the default constitution better
! IGNORE UNSTABLE COMP.SETS
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,iph,ics,lokics,jcs,lokjcs,bestfit(9,9),jj,kk,kki,kkj
   integer moveto(9)
   character ch1*1
!   write(*,*)'3Y in shiftcompsets'
   phloop: do iph=1,noofph
      lokph=phases(iph)
      manycs: if(phlista(lokph)%noofcs.gt.1) then
! seach all compset which default constitution that fits best a stable one
         bestfit=0
         csloop1: do ics=1,phlista(lokph)%noofcs
            lokics=phlista(lokph)%linktocs(ics)
! ignore UNSTABLE compsets with CSAUTO set ??
            if(ceq%phase_varres(lokics)%phstate.le.PHENTERED) cycle csloop1
            call checkdefcon(lokics,lokics,kk,ceq)
!            write(*,*)'3Y fit 1: ',kk,phlista(lokph)%tnooffr
            bestfit(ics,ics)=kk
            if(kk.eq.phlista(lokph)%tnooffr) cycle csloop1
! if no default or not perfect fit compare with other compsets
!            write(*,*)'3Y compare with next compset'
            csloop2: do jcs=1,phlista(lokph)%noofcs
               if(jcs.eq.ics) cycle csloop2
               lokjcs=phlista(lokph)%linktocs(jcs)
               if(ceq%phase_varres(lokjcs)%phstate.le.PHENTERED) cycle csloop2
               call checkdefcon(lokics,lokjcs,kk,ceq)
!               write(*,*)'3Y fit 2: ',kk,phlista(lokph)%tnooffr
               bestfit(jcs,ics)=kk
            enddo csloop2
         enddo csloop1
!         do ics=1,phlista(lokph)%noofcs
!            write(*,17)(bestfit(jcs,ics),jcs=1,phlista(lokph)%noofcs)
!         enddo
17       format('3Y bestfit: ',9i5)
! when we are here whe can use bestfit to shift constitutions
         moveto=0
         shiftfrom: do ics=1,phlista(lokph)%noofcs
            kk=bestfit(ics,ics)
            if(kk.eq.phlista(lokph)%tnooffr) cycle shiftfrom
            shiftto: do jcs=2,phlista(lokph)%noofcs
               if(bestfit(jcs,ics).gt.kk) then
                  kk=bestfit(jcs,ics)
                  write(*,*)'3Y shifting: ',ics,jcs
                  call switch_compsets2(lokph,ics,jcs,ceq)
               endif
            enddo shiftto
         enddo shiftfrom
! just check do nothing for the moment ....
! if moveto(ics) is zero do not move.  otherwise moveto moveto(ics)
! but if moveto(moveto(ics)) is zero look for a moveto() that is negative ...
!         call switch_compsets2(lokph,ics,jcs,ceq)
      endif manycs
   enddo phloop
! 
1000 continue
   return
 end subroutine shiftcompsets

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine switch_compsets(iph,ics1,ics2,ceq)
! copy constitution and results from ic2 to ic1 and vice versa
   integer iph,ics1,ics2
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,lokcs1,lokcs2
! check indices are correct
   call get_phase_compset(iph,ics1,lokph,lokcs1)
   call get_phase_compset(iph,ics2,lokph,lokcs2)
   if(gx%bmperr.ne.0) goto 1000
   call switch_compsets2(lokph,ics1,ics2,ceq)
1000 continue
   return
 end subroutine switch_compsets

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine switch_compsets2(lokph,ics1,ics2,ceq)
! copy constitution and results from ic2 to ic1 and vice versa
   integer lokph,ics1,ics2
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer iph,lokcs1,lokcs2,ncon,idum,ncc
   double precision, dimension(:), allocatable :: val
   double precision, dimension(:,:), allocatable :: gval,d2gval
   double precision, dimension(:,:,:), allocatable :: dgval
   double precision qq(5),xdum
!
!   write(*,*)'3Y In switch_compsets ',lokph,ics1,ics2
   lokcs1=phlista(lokph)%linktocs(ics1)
   lokcs2=phlista(lokph)%linktocs(ics2)
! save current constitution of lokcs1 in val
   ncon=size(ceq%phase_varres(lokcs1)%yfr)
   allocate(val(ncon))
   val=ceq%phase_varres(lokcs1)%yfr
! set the constitution in lokcs1 equal to that in lokcs2.  This call 
! also updates a number of other variables in the record
   iph=phlista(lokph)%alphaindex
   call set_constitution(iph,ics1,ceq%phase_varres(lokcs2)%yfr,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
   call set_constitution(iph,ics2,val,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! copy some variables: phstate, amfu and dgm
   idum=ceq%phase_varres(lokcs1)%phstate
   ceq%phase_varres(lokcs1)%phstate=ceq%phase_varres(lokcs2)%phstate
   ceq%phase_varres(lokcs2)%phstate=idum
   xdum=ceq%phase_varres(lokcs1)%amfu
   ceq%phase_varres(lokcs1)%amfu=ceq%phase_varres(lokcs2)%amfu
   ceq%phase_varres(lokcs2)%amfu=xdum
   xdum=ceq%phase_varres(lokcs1)%dgm
   ceq%phase_varres(lokcs1)%dgm=ceq%phase_varres(lokcs2)%dgm
   ceq%phase_varres(lokcs2)%dgm=xdum
! listprop will be the same
! Now copy result arrays
   ncon=ceq%phase_varres(lokcs1)%nprop
   allocate(gval(6,ncon))
   gval=ceq%phase_varres(lokcs1)%gval
   ceq%phase_varres(lokcs1)%gval=ceq%phase_varres(lokcs2)%gval
   ceq%phase_varres(lokcs2)%gval=gval
! ceq%phase_varres(lokph1)%ncc is not the dimension of dgval, why??
   ncc=size(ceq%phase_varres(lokcs1)%yfr)
   allocate(dgval(3,ncc,ncon))
!   write(*,77)'3Y copycomp: ',ncc,ncon,&
!        size(dgval),size(ceq%phase_varres(lokcs1)%dgval),&
!        ceq%phase_varres(lokcs2)%ncc,
!77 format(a,10i5)
   dgval=ceq%phase_varres(lokcs1)%dgval
   ceq%phase_varres(lokcs1)%dgval=ceq%phase_varres(lokcs2)%dgval
   ceq%phase_varres(lokcs2)%dgval=dgval
   allocate(d2gval(ncc*(ncc+1)/2,ncon))
   d2gval=ceq%phase_varres(lokcs1)%d2gval
   ceq%phase_varres(lokcs1)%d2gval=ceq%phase_varres(lokcs2)%d2gval
   ceq%phase_varres(lokcs2)%d2gval=d2gval
! addg!!
!   if(btest(ceq%phase_varres(lokcs1)%status2
   if(allocated(ceq%phase_varres(lokcs1)%addg)) then
      val(1)=ceq%phase_varres(lokcs1)%addg(1)
      ceq%phase_varres(lokcs1)%addg(1)=ceq%phase_varres(lokcs2)%addg(1)
      ceq%phase_varres(lokcs2)%addg(1)=val(1)
   endif
! curlat, cinvy, cxmol, cdxmol?
1000 continue
! deallocate
   deallocate(val)
   deallocate(gval)
   deallocate(dgval)
   deallocate(d2gval)
   return
 end subroutine switch_compsets2

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine shiftcompsets2(lokph,ceq)
! check if the composition sets of phase lokph
! should be shifted to fit the default constitution better
   integer lokph
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer kk,lokics,lokjcs,fit
! A fraction with a maximum set (mmyfr>0) must be below that value
! A fraction with a minimum set (mmyfr<0) should be above that value
   write(*,*)'3Y Not implemented testing defaults'
   fit=0
   do kk=1,size(ceq%phase_varres(lokics)%yfr)
      write(*,*)'3Y defconst: ',kk,ceq%phase_varres(lokics)%mmyfr(kk)
      if(ceq%phase_varres(lokics)%mmyfr(kk).gt.0.0D0) then
! A fraction with a maximum set (mmyfr>0) must be below that value
         if(ceq%phase_varres(lokjcs)%yfr(kk).gt.&
              ceq%phase_varres(lokics)%mmyfr(kk)) fit=fit+5
! A fraction with a minimum set (mmyfr<0) should be above that value
      elseif(ceq%phase_varres(lokics)%mmyfr(kk).lt.0.0D0) then
         if(ceq%phase_varres(lokjcs)%yfr(kk).lt.&
              abs(ceq%phase_varres(lokics)%mmyfr(kk))) fit=fit+1
      endif
   enddo
1000 continue
   return
 end subroutine shiftcompsets2

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine select_composition_set(iph,ics,yarr,ceq)
! if phase iph wants to become stable and there are several composition sets
! this subroutine selects the one with default composition set that fits best.
! For example if an FCC phase that could be an austenite (low carbon content)
! or a cubic carbo-nitride (high carbon or nitrogen content, low vacancy)
! Less easy to handle ordered phases like B2 or L1_2 as ordering can be
! in any sublatittice ... option B and F needed
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision, dimension(*) :: yarr
   integer iph,ics
!\end{verbatim}
   double precision, parameter :: yl=0.1D0,yh=0.5D0
   integer best,lokph,maxnh,ncc,jcs,lokcs,nh,jl
   lokph=phases(iph)
   best=1
   maxnh=0
   ncc=phlista(lokph)%tnooffr
   do jcs=1,phlista(lokph)%noofcs 
! loop through all composition sets
      lokcs=phlista(lokph)%linktocs(jcs)
! compare yarr with ceq%phase_varres(lokcs)%mmyfr
      nh=0
      do jl=1,ncc
         if(ceq%phase_varres(lokcs)%mmyfr(jl).lt.zero) then
            if(yarr(jl).lt.yl) nh=nh+1
         elseif(ceq%phase_varres(lokcs)%mmyfr(jl).gt.zero) then
            if(yarr(jl).gt.yh) nh=nh+1
         endif
      enddo
      if(nh.gt.maxnh) then
         maxnh=nh
         best=jcs
      endif
   enddo
! if only one compset return 1
   ics=best
!
1000 continue
   return
 end subroutine select_composition_set

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine verify_phase_varres_array(ieq,verbose)
! This subroutine checks that the phase varres array is consistent
! in equilibrium ieq.  For ieq=1 it also checks the free list
   implicit none
   integer ieq,verbose
!\end{verbatim}
   integer free,lokcs,lokph   
   type(gtp_phase_varres), pointer :: vares
   type(gtp_equilibrium_data), pointer :: ceq
   ceq=>eqlista(ieq)
   if(ieq.eq.1) then
! check free list inside phase_varres records
      if(csfree.lt.1 .or. csfree.gt.size(ceq%phase_varres)) then
         write(*,*)'3Y ERROR: csfree value outside limits: ',csfree
         goto 1000
      endif
      lokcs=csfree
50    continue
      if(lokcs.lt.1 .or. lokcs.gt.size(ceq%phase_varres)) then
         write(*,*)'3Y ERROR: varres free list index outside limits: ',lokcs
         goto 1000
      endif
      lokcs=ceq%phase_varres(lokcs)%nextfree
      if(lokcs.lt.size(ceq%phase_varres)) goto 50
!-------
      write(*,*)'3Y varres free list seems OK.'
   endif
!-------
! check each used varres record that it has a correct phase pointer etc.
! UNFINISHED
1000 continue
      return
 end subroutine verify_phase_varres_array

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!-\begin{verbatim}
 subroutine sort_multidim_array(nrel,ng,xarr)
! sorts values in xarr in decending order for each column (element)
! xarr(1..nrel,jp) is the composition of gridpoint jp
! UNFINISHED
   implicit none
   integer nrel,ng
   real xarr(nrel,*)
!-\end{verbatim}
   integer i1,j1,k1,l1,m1,bounds(nrel,10),nb,b1,c1
   real, dimension(:,:), allocatable :: xord
   real xx,xsame
! allocate ordered array and zero. All values in xarr>0
   allocate(xord(nrel,ng))
   xord=zero
! sort first column, increment is nrel, very brute force 
   loop1: do j1=1,ng,nrel
      xx=xarr(1,j1)
      do k1=1,j1-1
         if(xx.gt.xord(1,k1)) then
! store xarr here but first shift all values in xord after k1 up, loop bacwards
            do l1=j1,k1,-1
               do m1=1,nrel
                  xord(m1,l1+1)=xord(m1,l1)
               enddo
            enddo
            do m1=1,nrel
               xord(m1,k1)=xarr(m1,j1)
            enddo
            cycle loop1
         endif
         do m1=1,nrel
            xord(m1,j1)=xarr(m1,j1)
         enddo
      enddo
   enddo loop1
! detect the bounds
   xsame=xord(1,1)
   bounds(1,1)=1
   j1=1
   do i1=1,ng
      if(xord(1,i1).lt.xsame) then
         xsame=xord(1,i1)
         j1=j1+1
         bounds(j1,1)=j1
      endif
   enddo
   nb=j1
   write(*,11)'3Y bounds column 1: ',nb,(bounds(k1,1),k1=1,nb)
11 format(a,i3,10i5)
! now the first column sorted and all values are in xord, sort columns 2 etc
! separately for each lower column bounds
! keep the ordering in the lower columns. No need to sort Last column
   column: do c1=1,nrel
      boundloop: do b1=2,nb
         loop2: do i1=bounds(b1-1,c1),bounds(b1,c1)
            loop3: do j1=2,ng,nrel
               if(xord(i1-1,j1).lt.xsame) then
                  xsame=xord(i1-1,1)
                  xx=xord(i1,j1)
                  do k1=1,ng
                     if(xx.gt.xord(i1,k1)) then
! store xx here, shift up (all values in xord)
                        do l1=j1,k1,-1
                           do m1=1,nrel
                              xord(m1,l1+1)=xord(m1,l1)
                           enddo
                        enddo
                     endif
                     xord(i1,k1)=xx
                     cycle loop3
                  enddo
               endif
            enddo loop3
         enddo loop2
      enddo boundloop
   enddo column
1000 continue
   deallocate(xord)
   return
! UNFINISHED
 end subroutine sort_multidim_array
!   
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

