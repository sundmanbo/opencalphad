!
! gtp3Y included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>      16. Grid minimizer
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine global_gridmin
!\begin{verbatim}
 subroutine global_gridmin(what,tp,xknown,nvsph,iphl,icsl,aphl,&
      nyphl,cmu,ceq)
!
! Starting rewriting 2017-02-01
!
! finds a set of phasey cons that is a global start point for an equilibrium 
! calculation at T and P values in tp and known mole fraction in xknown
! It is intentional that this routine is independent of current conditions
! It returns: nvsph stable phases, list of phases in iphl, amounts in aphl, 
! nyphl(i) is redundant, cmu are element chemical potentials of solution
! WHAT determine what to do with the results, 0=just return solution,
! 1=enter stable set and constitution of all phases in gtp datastructure
! and create composition sets if necessary (and allowed)
! what=-1 will check if any gridpoint is below current calculated equilibrium
! ?? removed what -1 170428/BoS
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
   integer pph,zph,nystph,order(maxel),tbase,qbase,wbase,jj,zz,errall,eecliq
!
!   write(*,*)'3Y in global_gridmin'
!   nystph=0
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
!   trace=.TRUE.
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
!   write(*,*)'3Y allocating gridpoints 1',nrph
   allocate(gridpoints(nrph),stat=errall)
   allocate(phord(nrph),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 1: ',errall
      gx%bmperr=4370; goto 1100
   endif
   eecliq=0
   if(globaldata%sysreal(1).gt.one) then
! we must initiate EEC data for the liquid
      sliqmin=zero; sliqmax=zero; gliqeec=zero
   endif
!   write(*,*)'3Y loop for all phases',nrph,globaldata%sysreal(1)
   ggloop: do iph=1,nrph
! include all phases with any composition set entered (but only once!)
      do ics=1,noofcs(iph)
! new: -3 suspended, -2 dormant, -1,0,1 entered, 2 fixed
! ignore phases whith no entered composition sets
! If a phase+compset FIX one should never be here as conditions wrong
         if(test_phase_status(iph,ics,amount,ceq).gt.PHDORM) then
            pph=pph+1
            iphx(pph)=iph
            if(eecliq.eq.0 .and. globaldata%sysreal(1).gt.one) then
!               write(*,*)'3Y looking for liquid: ',pph,iph
               if(btest(phlista(iph)%status1,PHLIQ)) then
                  eecliq=phlista(iph)%alphaindex
               endif
            endif
            cycle ggloop
         endif
      enddo
   enddo ggloop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! if lutbug>0 open a file for grid graphics
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!   lutbug=37
! There is some other gridoutput associated with trace
   if(lutbug.gt.0) then
      write(*,*)'3Y Opening gridmap.dat'
      open(lutbug,file='gridmap.dat',access='sequential',status='unknown')
   endif
! we will generate a grid for pph phases, the phase index for phase 1..pph
! is in iphx(1..pph)
! always allocate a grid for maxgrid points
!   write(*,*)'3Y allocating gridpoints 2',nrel,maxgrid
   allocate(xarr(nrel,maxgrid),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 2: ',errall
      gx%bmperr=4370; goto 1000
   endif
   gmax=zero
!   write(*,11)'3Y gp1:',pph,(iphx(iph),iph=1,pph)
! just to be sure
   nphl(0)=0
!
!------------------------------------------------------
! For EEC we must always calculate the liquid phase first
! it means its grid will be calculated twice because I do not want to change
   eec: if(eecliq.gt.0) then
! if eecliq=1 the liquid is first and this is not needed
! ng should be set to number of remaining points, ny and yphl is not used
      iv=1
      ng=maxgrid
! possible output on gridmap.dat
!            if(lutbug.gt.0) then
!               lokph=phases(iphx(zph))
!               write(lutbug,16)trim(phlista(lokph)%name),zph,lokph,ng
!16             format(/'Phase: ',a,3i7)
!            endif
!      write(*,*)'3Y calculate EEC data for liquid',eecliq,iphx(eecliq),iv,ng
      if(btest(globaldata%status,GSOGRID)) then
! The possibility to use the old grid tested
         call generate_grid(0,eecliq,ng,nrel,&
              xarr(1,iv),garr(iv),ny,yarr,gmax,ceq)
      else
         call generic_grid_generator(0,eecliq,ng,nrel,&
              xarr(1,iv),garr(iv),ny,yarr,gmax,ceq)
      endif
      if(gx%bmperr.ne.0) then
         write(*,*)'3Y grid error ',eecliq,gx%bmperr
         exit eec
      endif
!      write(*,*)'3Y sliqmax: ',sliqmax
   elseif(globaldata%sysreal(1).gt.one) then
      write(*,*)'3Y EEC will not work because there is no liquid'
   endif eec
!
!----------------------------------------------------------
   phloop: do zph=1,pph
! for phase iphx(zph) the gridpoints will be stored from position nphl(zph-1)+1
! ng should be set to number of remaining points, ny and yphl is not used
      iv=nphl(zph-1)+1
      ng=maxgrid-iv
!      write(*,*)'3Y generating grid for phase ',zph
! possible output on gridmap.dat
      if(lutbug.gt.0) then
         lokph=phases(iphx(zph))
         write(lutbug,16)trim(phlista(lokph)%name),zph,lokph,ng
16       format(/'Phase: ',a,3i7)
      endif
!      write(*,*)'3Y gridgen ',zph,iv,ng
! this call will calculate gridpoints in phase zph, that may take time ...
! ng is set to remaining dimension of garr, on return the number of generated
!    gridpoints, returned as xarr composition of these and
! ny and yarr not used here
!>>>>>>> important: changes here must be made also in global_equil_check
!      write(*,*)'3Y grid for phase: ',zph,phlista(phases(iphx(zph)))%name
      if(btest(globaldata%status,GSOGRID)) then
! The possibility to use the old grid tested
         call generate_grid(0,iphx(zph),ng,nrel,xarr(1,iv),garr(iv),&
              ny,yarr,gmax,ceq)
      else
         call generic_grid_generator(0,iphx(zph),ng,nrel,xarr(1,iv),garr(iv),&
              ny,yarr,gmax,ceq)
      endif
!>>>>>>>> impportant end!
!      write(*,*)'3Y grid done'
      if(gx%bmperr.ne.0) then
         write(*,*)'3Y grid error ',iphx(zph),gx%bmperr
         exit phloop
      endif
! nphl(zph) is last gridpoint in phase zph
      nphl(zph)=nphl(zph-1)+ng
   enddo phloop
!----------------------------------------------------------
! if lutbug>0 close it
   if(lutbug.gt.0) then
      close(lutbug)
      write(*,*)'3Y closed gridmap.dat with gridpoints'
   endif
!++   write(*,11)'3Y gp2:',(nphl(iph),iph=1,nrph)
11 format(a,10i7/(7x,10i7))
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
      write(*,*)'3Y Calling global_grimin with -1 no longer supported'
      stop
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
!   write(*,*)'3Y total gridpoints: ',kp
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
      if(btest(globaldata%status,GSOGRID)) then
! The possibility to use the old grid tested
         call generate_grid(mode,iphx(zph),ng,nrel,xarr(1,iv),garr(iv),&
              ny,yarr,gmax,ceq)
      else
         call generic_grid_generator(mode,iphx(zph),jbias,nrel,xarr,garr,&
              ny,yarr,gmax,ceq)
      endif
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
!   write(*,*)'3Y merge in global?',btest(globaldata%status,GSNOMERGE)
   if(.not.btest(globaldata%status,GSNOMERGE)) then
! For the moment we will only merge grid points in the gas phase
      call merge_gridpoints(nr,iphl,aphl,nyphl,yphl,trace,nrel,xsol,cmu,ceq)
      if(gx%bmperr.ne.0) goto 1000
   endif
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
!   write(*,*)'3Y Constitution of metastable phases set'
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
! this means status entered and unknown state. PHSTATE
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
!   write(*,*)'3Y allocating startup 3',nvsph
   allocate(starttup(nvsph),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 3: ',errall
      gx%bmperr=4370; goto 1000
   endif
   call extract_massbalcond(ceq%tpval,xdum,totam,ceq)
   if(gx%bmperr.ne.0) goto 1000
   sum=zero
   do iph=1,nvsph
      sum=sum+aphl(iph)
   enddo
   gmax=zero
!   write(*,*)'3Y no segmentation fault 1'
! If there is a gas nvsph may be less than number of elements
   ceqstore: do iph=1,nvsph
!      write(*,*)'3Y no segmentation fault 2',iph
      if(iphl(iph).lt.0) then
! this gripoint has no composition set because too many gridpoints in same phase
         starttup(iph)=0
         continue
!         write(*,*)'3Y no segmentation fault 3',iph
      else
!         write(*,*)'3Y segmentation fault 5',iph,iphl(iph),icsl(iph),j1
         call set_constitution(iphl(iph),icsl(iph),yphl(j1),qq,ceq)
         if(gx%bmperr.ne.0) goto 1000
!         write(*,*)'3Y no segmentation fault 6',lokph,lokcs
         call get_phase_compset(iphl(iph),icsl(iph),lokph,lokcs)
         if(gx%bmperr.ne.0) goto 1000
!         write(*,*)'3Y no segmentation fault 7'
! This is a bit quiestionable but seems to work
         amount=aphl(iph)/ceq%phase_varres(lokcs)%abnorm(1)
         gmax=gmax+amount
         aphl(iph)=amount
1789     format(a,2i4,5(1pe12.4))
!         write(*,*)'3Y no segmentation fault 8',lokcs,iph
         ceq%phase_varres(lokcs)%amfu=aphl(iph)
         ceq%phase_varres(lokcs)%phstate=PHENTSTAB
         starttup(iph)=ceq%phase_varres(lokcs)%phtupx
         j1=j1+nyphl(iph)
      endif
!      write(*,*)'3Y no segmentation fault 9',iph,nvsph
   enddo ceqstore
!-----------------------------------------------------------------------
! debug listing of tuples at gridpoints
!   write(*,*)'3Y no segmentation fault 10'
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
!   write(*,*)'3Y no segmentation fault 20'
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
1100 continue
   if(ocv()) write(*,*)'3Y leaving global_gridmin'
   write(*,*)'3Y leaving global_gridmin'
   return
 end subroutine global_gridmin

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine generate_grid
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
   integer jbas,sumngg,loksp,wrongngg,errall
   logical trace,isendmem
   save sumngg,wrongngg
!
!   write(*,*)'Illegal call to generate_grid: ',mode
!   stop
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
   elseif(test_phase_status_bit(iph,PHFORD) .or. &
        test_phase_status_bit(iph,PHBORD)) then
! this phase has 4 sublattice fcc/hcp tetrahedral ordering,
! this reduces the number of gridpoints UNFINISHED: NOT IMPLEMENTED YET
!      write(*,*)'3Y calling ordered grid 1'
      call generate_fccord_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! do not jump to 1000 until the fccord routine implemented correctly
!      write(*,*)'3Y back from fccord_grid 1, jump to 1000',ngg
! This routine return gx%bmperr=-1 if if cannot handle the gridgenerating
      if(gx%bmperr.eq.-1) then
         gx%bmperr=0
      else
         goto 1000
      endif
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
!   write(*,*)'3Y allocating endmem: ',nsl,nend,inkl(nsl)
   allocate(endm(nsl,nend),stat=errall)
   allocate(yfra(inkl(nsl)),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 4: ',errall
      gx%bmperr=4370; goto 1001
   endif
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
! jump here from generate_fccord_grid  ... not any more ...
170 continue
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
! initiate the loop variables below for endmembers and fractions
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
         if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
              write(*,*)'3Y calculated ',ngg,' gridpoints for ',&
              trim(phlista(lokph)%name)
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
! we are iterating in the ternary endmember
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

!\addtotable subroutine generic_grid_generator
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
   integer ii,ij,ik,il,im,in,is,ie,iz,incl(0:maxsubl),maxng,ng,ncon
   integer nend,nendj,nendk,nendl,nendm
   integer ijs,iks,ils,ims,lokph,errall
! these are for call of get_phase_data
   integer nsl,nkl(maxsubl),knr(maxconst)
   double precision ydum(maxconst),sites(maxsubl),qq(5)
! this is for generating endmembers
   integer, dimension(:,:), allocatable :: endm
! these are for generating the constituion of gridpoint
   double precision, dimension(:,:), allocatable :: yendm
   double precision, dimension(:), allocatable :: yfra
   integer :: warning1const=0
   save warning1const
! ----------------------------------------------------------------
! these are the factors to generate gridpoints from endmember fractions
   double precision, dimension(5), parameter :: &
!        yf=[0.33D0,0.28D0,0.18D0,0.08D0,0.03D0]
!        yf=[0.33D0,0.28D0,0.18D0,0.08D0,0.03D0] to test with map3
! ok        yf=[0.33D0,0.28D0,0.18D0,0.14D0,0.11D0]
! better but fails Fe-C at 1100 K and w(c)=0.03
!        yf=[0.11D0,0.33D0,0.14D0,0.28D0,0.18D0]
! include a small factor
!        yf=[0.11D0,0.37D0,0.04D0,0.30D0,0.18D0]
! Try to avoid several identical compositions
        yf=[0.07D0,0.28D0,0.16D0,0.45D0,0.04D0]
!        yf=[0.33D0,0.28D0,0.18D0,0.14D0,0.11D0] OK for fuel
!        yf=[0.11D0,0.13D0,0.18D0,0.23D0,0.35D0]
!------------------------------------------------------------------
   logical gas,dense,verydense,gles,trace
! bugfix by Clement Instroini 18.02.14
   if(iph.lt.1 .or. iph.gt.noofph) then
      gx%bmperr=4050; goto 1000
   else
      lokph=phases(iph)
   endif
! handle special phases like ionic crystals, ionic liquids and order/disorder
!   write(*,*)'3Y in generic_grid_generator',iph
   gas=.FALSE.
! to have some output
!   if(mode.gt.0 .and. ny.eq.-100) then
!   if(mode.gt.0) then
!      write(*,*)'3Y turn on trace',mode,iph,ngg,ny
!      trace=.TRUE.
!   else
!      write(*,*)'3Y searching for y: ',mode,iph,ngg,ny
!   endif
   if(test_phase_status_bit(iph,PHEXCB)) then
! crystalline phase with charged endmembers
!      write(*,*)'3Y charged grid ngg: ',ngg
      call generate_charged_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
      goto 1000
   elseif(test_phase_status_bit(iph,PHIONLIQ)) then
! This is the ionic liquid, requires a special grid, also used for dense
      call generate_ionliq_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
      goto 1000
   elseif(test_phase_status_bit(iph,PHGAS)) then
      gas=.TRUE.
   elseif(test_phase_status_bit(iph,PHFORD) .or. &
        test_phase_status_bit(iph,PHBORD)) then
!      write(*,*)'3Y calling ordered grid 2'
      call generate_fccord_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
!      write(*,*)'3Y back from fccord_grid 2, jump to 1000',ngg
!      goto 200
      if(gx%bmperr.eq.-1) then
! if gx%bmperr is -1 means problems in fccord_grid, use default grindgenerator
         gx%bmperr=0
      else
         goto 1000
      endif
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
!   write(*,*)'3Y allocating endmem 2',nsl,nend,ncon,nend
   allocate(endm(nsl,nend),stat=errall)
   allocate(yendm(ncon,nend),stat=errall)
   allocate(yfra(ncon))
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 5: ',errall
      gx%bmperr=4370; goto 1000
   endif
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
! jump here from generate_fccord_grid ... not any longer ...
200 continue
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
            if(ng.gt.0 .and. mod(ng,30000).eq.0) then
               lokph=phases(iph)
               write(*,*)'3Y calculated ',ng,' gridpoints for ',&
                    trim(phlista(lokph)%name)
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
   gles=.not.dense
!   write(*,*)'3Y grid1: ',nend,mode,dense,gles
   iiloop: do ii=1,nend
      ijs=1
      nendj=nend
      if(nend.gt.150 .or. (gles .and. nend.gt.50)) then
         nendj=ii
         ijs=ii
      endif
!      write(*,*)'3Y ii:',ii,ijs,nendj
      ijloop: do ij=ijs,nendj
         iks=1
         nendk=nend
         if(nend.gt.15 .or. (gles .and. nend.gt.13)) then
            nendk=ij
            iks=ij
         endif
         ikloop: do ik=iks,nendk
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
!                  write(*,323)'3Y imloop1: ',ng,ii,ij,ik,il,im,0.0D0,yfra
323                format(a,i5,5i3,': ',1pe12.4,0p10F5.2)
                  if(mode.eq.0) then
! strange bug in map3, maxng was zero sometimes ...
                     if(ng.gt.maxng) then
                        if(maxng.lt.100) then
                           if(warning1const.ne.iph) then
                              write(*,*)'3Y max gripoints wrong 6: ',&
                                   maxng,iph,mode
                              warning1const=iph
                           endif
                        else
                           write(*,*)'3Y Too many gridpoints 6',ng,maxng,iph
                           gx%bmperr=4399; goto 1000
                        endif
                     endif
                     if(ng.gt.0 .and. mod(ng,30000).eq.0) &
                          write(*,*)'3Y calculated ',ng,' gridpoints for ',&
                          trim(phlista(lokph)%name)
                     call calc_gridpoint(iph,yfra,nrel,xarr(1,ng),garr(ng),ceq)
! generate a GNUPLOT graph for dense grid
                     if(mode.eq.0 .and. lutbug.gt.0) then
!                        write(*,710)ng,nrel,garr(ng),(xarr(iz,ng),iz=1,nrel)
                        write(lutbug,710)ng,nrel,garr(ng),&
                             (xarr(iz,ng),iz=1,nrel)
710                     format(i5,i3,1pe10.2,10(0pF6.3))
                     endif
                     if(gx%bmperr.ne.0) then
                        write(*,*)'3Y error calculating gridpoint: ',&
                             iph,gx%bmperr
                        goto 1000
                     endif
!                     write(*,323)'3Y imloop2: ',ng,ii,ij,ik,il,im,garr(ng),yfra
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
      write(*,*)'3Y Could not find gridpoint ',mode,' in phase ',iph,ng
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

!\addtotable subroutine generate_dense_grid
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
   integer jbas,sumngg,loksp,l0,l1,ncon,jj,anion,isp,errall
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
!   write(*,*)'3Y allocating yfra mm',inkl(nsl),nsl,nend
   allocate(yfra(inkl(nsl)),stat=errall)
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
! endm(1,1) is constituent in sublattice 1 of first endmember
! endm(2,1) is constituent in sublattice 2 of first endmember
! endm(nsl,2) is constituent in sublattice nsl of second endmember
! endm(1..nsl,nend) are constituents in all sublattices of last endmember
   allocate(endm(nsl,nend),stat=errall)
! inkl(nsl) is the number of fraction variables in the phase
!   allocate(yfra(inkl(nsl)))
   allocate(xbrr(noofel),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 6: ',errall
      gx%bmperr=4370; goto 1000
   endif
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
! this is for a single endmember
!         write(*,201)'3Y end: ',ngg,(yfra(is),is=1,inkl(nsl))
         if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
              write(*,*)'3Y calculated ',ngg,' gridpoints 4',&
              trim(phlista(lokph)%name)
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
            if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
                 write(*,*)'3Y calculated ',ngg,' gridpoints for ',&
                 trim(phlista(lokph)%name)
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
                  if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
                       write(*,*)'3Y calculated ',ngg,' gridpoints for ',&
                       trim(phlista(lokph)%name)
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

!\addtotable subroutine generate_ionliq_grid
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
! try to have denser cation grid when Va as cation
   integer anionva,constva,hasva
   integer, allocatable, dimension(:) :: endwithva
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
! Used when mixing cations with same anion
   double precision, dimension(1:3), parameter:: yfx=&
        [0.42D0,0.14D0,0.08D0]
! for output of gridpoints
   integer l1,ncon,jj,cation,anion,isp,iva,catloop,neutral1,errall
   integer looplim1,looplim2
   logical trace,dense
   character ch1*1,dummy*128
!
   if(mode.eq.0) then
      trace=.FALSE.
!      trace=.TRUE.
!     write(*,*)'3Y Calculating the number of gridpoints for ionic liquid'
! trace TRUE means generate outpt for each gridpoint
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
!   write(*,*)'3Y Getting phase data',mode
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! calculate the number of endmembers and index of first constituent in subl ll
   nend=1
   inkl(0)=0
   inkl(1)=nkl(1)
   cation=nkl(1)
! Why is inkl(2) set like this?  I have changed ncon below
   inkl(2)=nkl(1)+nkl(2)
   lokph=phases(iph)
   if(.not.btest(phlista(lokph)%status1,PHIONLIQ)) then
      write(*,*)'3Y internal error, this phase has not ionic liquid model!'
      gx%bmperr=4399; goto 1000
   endif
! multiply with charged anions and Va only, add neutrals, nsl=2
   anion=0
   anionva=0
   do jj=1,nkl(2)
! knr(i) is species(i) location but I use constitlist as I have access to it
      isp=phlista(lokph)%constitlist(nkl(1)+jj)
      if(btest(splista(isp)%status,SPVA)) anionva=jj
      if(btest(splista(isp)%status,SPION) .or. &
           btest(splista(isp)%status,SPVA)) then
         anion=anion+1
         cycle
      endif
   enddo
! If no Va in anion sublattice anionva=0
   if(anionva.gt.0) then
! if anaionva>0 allocate array for endmembers with Va
      allocate(endwithva(nkl(1)))
! save constituent index for Va
      constva=nkl(1)+anionva
      hasva=1
!      write(*,*)'3Y anionva: ',mode,anionva,nkl(1),constva
   endif
! error when compiling with  -O2
!   nend=inkl(1)*anion+phlista(lokph)%nooffr(2)-anion
   nend=nkl(1)*anion+nkl(2)-anion
   ny=inkl(nsl)
!   ncon=inkl(nsl)
! BoS corrected 2019/04/13: (U+4,Zr+4)(O-2,Va,O) has 5 constituents not 6
   ncon=inkl(1)+nkl(2)
!   write(*,45)'3Y liquid endmembers: ',mode,nkl(1),nkl(2),anion,nend
45 format(a,5i5)
   negmode: if(mode.lt.0) then
!---------------------------------------------------------
! this is never executed as mode<0 no longer used
! just estimate the number of gridpoints for the ionic liquid phase
! pairs of cation+anion, cation+Va, neutrals
      ngdim=ngg
      ngg=nend
      lokph=phases(iph)
      write(*,*)'3Y nend 1: ',mode,ngg,breaks
      if(nend.eq.1 .or. nend.gt.breaks(4)) then
!
! NOTE mode<0 is NO LONGER USED, this code not used <<<<<<<<<<<<!!!
!
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
!      read(*,11)ch1
11    format(a)
      ny=nend
      goto 1001
   endif negmode
! the negmode if statement above no longer used ^^^^^^^^^^^^^^^^^^^^
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
!   write(*,*)'3Y allocating endmembers 5:',nsl,nend,inkl(nsl)
   allocate(endm(nsl,nend),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 7: ',errall
      gx%bmperr=4370; goto 1000
   endif
! inkl(nsl) is the number of fraction variables in the phase
!   allocate(yfra(inkl(nsl)))
!   nofy=inkl(nsl)
! generate endmembers, endm(ll,ie) is set to consituent index in sublattice ll
! For neutrals in sublattice 2, sublattice 1 has -99 as constituent
   je=1
   do ll=1,nsl
      endm(ll,je)=inkl(ll-1)+1
   enddo
! we may not have any anions, just Va!
   isp=endm(2,je)
   if(btest(splista(knr(isp))%status,SPVA)) then
! save index for endmembers with va as anion
!      write(*,*)'3Y ionic liquid endmember with Va 1: ',je,hasva,nkl(1)
      endwithva(hasva)=je
      hasva=hasva+1
   endif
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
         if(hasva.gt.1) then
!            write(*,*)'3Y ionic liquid endmember with Va 3: ',je,hasva,nkl(1)
            endwithva(hasva)=je
            hasva=hasva+1
         endif
      else
         endm(1,je)=1
         isp=endm(2,je)+1
         if(btest(splista(knr(isp))%status,SPVA)) then
! save index for endmembers with va as anion
!            write(*,*)'3Y ionic liquid endmember with Va 2: ',je,hasva,nkl(1)
            endwithva(hasva)=je
            hasva=hasva+1
         endif
         if(splista(knr(isp))%charge.eq.zero .and. &
              .not.btest(splista(knr(isp))%status,SPVA)) then
! The next constituent in second sublattice is not Va or a neutral
            exit genend
         else
! the next constituent is an anion or Va
            endm(2,je)=endm(2,je)+1
         endif
      endif
!      write(*,171)'3Y endmember 1: ',je,endm(1,je),endm(2,je)
   enddo genend
171 format(a,i3,'  (',i2,':',i2,')')
! we must generate endmembers for neutrals, wildcard -99 in first sublattice
   do iend=je,nend
      endm(1,iend)=-99
      endm(2,iend)=isp
      isp=isp+1
!      write(*,171)'3Y endmember 2: ',je,endm(1,je),endm(2,je)
   enddo
! debug check
   if(mode.eq.0) then
!      write(*,*)'3Y NEW: mode, endmembers, gridpoints: ',mode,nend,ngg
!      write(*,111)(endm(1,je),endm(2,je),je=1,nend)
111   format('3Y list: ',10(i4,i3)/11(i4,i3))
   endif
!   gx%bmperr=4399; goto 1000
!---------------------------------------
! We have now generated endm(1..nsl,j)
!120 continue
!   write(*,202)'3Y special 1: ',nsl,nend,inkl(nsl),endm(1,2),endm(2,2),&
!        endm(1,nend),endm(1,3),endm(2,3),endm(1,4),endm(2,4)
! we must allocate and set endmember fractions both for mode 0 and >0
!   if(mode.gt.0) write(*,*)'3Y allocate yendm: ',inkl(2),nend
!   write(*,*)'3Y allocating endmembers 6:',inkl(2),nend
   allocate(yendm(inkl(2),nend),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 8: ',errall
      gx%bmperr=4370; goto 1000
   endif
   yendm=zero
! SOMETIMES IT CRASHED IF THIS LINE IS REMOVED
   write(dummy,*)'3Y endmember fractions:',mode,je
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
!   allocate(yfra(inkl(nsl)))
! this is a small allocation, max 1000 double
   allocate(yfra(ncon))
!---------------------------------------------
! now generate combinations of endmember fractions
! Note the sum of constituent fractions in all sublattices should be unity
   ngg=0
   looplim1=breaks(4)
   looplim2=breaks(3)
   if(dense) then
!      write(*,*)'3Y dense ionic liquid grid',nend
      looplim1=nend+1
      looplim2=breaks(4)
   endif
   lokph=phases(iph)
   endmem1: do iend=1,nend
      endmem2: do jend=1,nend
         if(nend.gt.looplim1 .and. jend.ne.iend) cycle endmem2
         endmem3: do kend=1,nend
            if(nend.gt.looplim2 .and. kend.ne.jend) cycle endmem3
            endmem4: do lend=1,nend
               if(nend.gt.breaks(2) .and. lend.ne.kend) cycle endmem4
               do jj=1,ncon
                  yfra(jj)=yf(1)*yendm(jj,iend)+yf(2)*yendm(jj,jend)+&
                       yf(3)*yendm(jj,kend)+yf(4)*yendm(jj,lend)
               enddo
               ngg=ngg+1
               if(mode.eq.0) then
                  if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
                       write(*,*)'3Y calculated ',ngg,' gridpoints for ',&
                       trim(phlista(lokph)%name)              
                  call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
                  if(gx%bmperr.ne.0) goto 1000
                  if(lutbug.gt.0) then
! debug output for NEW version of ionic liquid, the grid was quite strange
                     write(lutbug,710)'A: ',ngg,nrel,ncon,garr(ngg),&
                          (xarr(jj,ngg),jj=1,nrel),(yfra(jj),jj=1,ncon)
710                  format(a,i5,2i3,1pe10.2,10(0pF6.3))
                  endif
                  if(garr(ngg).gt.gmax) gmax=garr(ngg)
!                  if(mod(ngg,10000).eq.0) write(*,*)'Calculated ',ngg,&
!                       ' gridpoints, more to go'
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
!   write(*,*)'3Y Calculated points 1: ',ngg,nend,breaks(2),dense
!   goto 1000
!   write(*,*)'3Y special cation loop: '
!   if(nend.le.breaks(2)) goto 1000
!   if(.not.dense .and. nend.le.breaks(2)) goto 1000
   if(.not.dense .and. cation.gt.breaks(3)) goto 1000
! combinations 3 different cations with same anion
   anion2: do lend=0,anion-1
      catloop=lend*cation+1
! REMEMBER: endmembers with same cation are ordered sequentially!!! 
!      write(*,*)'3Y catloop: ',lend+1,cation,catloop,ngg
! calculating "c g" followed by "c n" gives better result than "c e", why??
! The reason was that I had forgotten to scale phase amounts with total moles
      endmem1b: do iend=catloop,catloop+cation-3
         endmem2b: do jend=iend+1,catloop+cation-2
            endmem3b: do kend=jend+1,catloop+cation-1
! these loops generate a  more dense grid but give same results in my test
!      endmem1b: do iend=catloop,catloop+cation-1
!         endmem2b: do jend=catloop,catloop+cation-1
!            if(jend.eq.iend) cycle endmem2b
!            endmem3b: do kend=catloop,catloop+cation-1
!               if(kend.eq.jend .or. kend.eq.iend) cycle endmem3b
               if(.not.dense .and. cation.gt.breaks(3)) then
!                  write(*,*)'3Y skipping ternary cationloop'
                  cycle endmem3b
               endif
! mixing of 3 cations with the same anion
               do jj=1,ncon
                  yfra(jj)=yfx(1)*yendm(jj,iend)+yfx(2)*yendm(jj,jend)+&
                       yfx(3)*yendm(jj,kend)
               enddo
               ngg=ngg+1
               if(mode.eq.0) then
                  if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
                       write(*,*)'3Y calculated ',ngg,' gridpoints for ',&
                       trim(phlista(lokph)%name)              
                  call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
                  if(lutbug.gt.0) then
! debug output for NEW version of ionic liquid, the grid was quite strange
                     write(lutbug,710)'B: ',ngg,nrel,ncon,garr(ngg),&
                          (xarr(jj,ngg),jj=1,nrel),(yfra(jj),jj=1,ncon)
                  endif
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
!   write(*,*)'3Y Calculated points 2: ',ngg,nend,breaks(2)
! skip next loop for the moment ....... still
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
                  if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
                       write(*,*)'3Y calculated ',ngg,' gridpoints for ',&
                       trim(phlista(lokph)%name)
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
   if(ny.ne.ncon) write(*,*)'3Y ny and ncon: ',ny,ncon
!   do jj=1,ny
   do jj=1,ncon
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

!\addtotable subroutine generate_fccord_grid
!\begin{verbatim} %-
 subroutine generate_fccord_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,ceq)
! This generates grid for a phase with 4 sublattice fcc/bcc/hcp ordering
! NO LONGER USED: mode<0 just number of gridpoints in ngg, for allocations
! mode=0 calculate mole fraction and G for all gridpoints
! mode>0 return constitution for gridpoint mode in yarr
   implicit none
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*),gmax
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   logical, save :: once=.TRUE.
   integer nsl,maxng,mend,nend,kend,ncon,i1,i2,i3,i4,i5,ij,ik,iz,ls,gridlimit
   integer nkl(1000),knr(1000),incl(0:9),loksp,lokph,jj
   integer, allocatable, dimension(:,:) :: endm
   double precision, allocatable, dimension(:,:) :: yendm
   double precision, allocatable, dimension(:) :: yfra
   double precision, allocatable, dimension(:) :: ysave
   double precision ydum(1000),ysame(1000),sites(9),qq(5)
! this generates >3000 gridpoints for a binary (A,B)(A,B)(A,B)(A,B)
!   integer, parameter, dimension(3,4) :: &
!        limits=reshape([150,50,20, 100,30,15, 20,10,7, 12,7,4],shape(limits))
   integer, parameter, dimension(3,4) :: &
        limits=reshape([150,50,20, 100,30,15, 20,10,7, 12,7,4],shape(limits))
! from generic_grid_generator
!   logical dense,gles,defgrid
   double precision, dimension(5), parameter :: &
        yf=[0.07D0,0.28D0,0.16D0,0.45D0,0.04D0]
   integer ii,ijs,iks,il,ils,im,ims,is,nendj,nendk,nendl,nendm,ng,errall
   character phname*32
! NOTHING IMPLEMENTED YET oh yes it is ...
!   write(*,*)'3Y in generate_fccord_grid ',ngg
   if(mode.lt.0) then
      write(*,*)'3Y mode <0 not allowed'
      gx%bmperr=4399; goto 1000
   endif
! check that F or B bit set
   if(.not.(test_phase_status_bit(iph,PHFORD) .or. &
        test_phase_status_bit(iph,PHBORD))) then
      write(*,*)'3Y calling ordered grid without F or B bit'
      gx%bmperr=4399; goto 1000
   endif
!
! get phase model
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1010
! max number of gridpoints allowed, ngg returned as number of gridpoints ???
! maxng is zero ...
   maxng=ngg
   ngg=0
! incl(ii) set to number of constituents up to including sublattice ii
   incl(0)=0
   nend=1
   do ij=1,nsl
      nend=nend*nkl(ij)
      incl(ij)=incl(ij-1)+nkl(ij)
   enddo
   ncon=incl(nsl)
! if nend<15 there is a single constituent on the ordered sublattices
   if(nend.lt.16) then
      gx%bmperr=-1
      goto 1010
   endif
! nend is number of endmembers, endm(1..nsl,ii) are constituent index of ii
! yendm(1..nsl,ii) has the constituent fractions for endmember ii
! yfra is used to generate a constitutuon from a combination of endmembers
!   write(*,*)'3Y allocating endmembers 8:',nsl,nend,ncon
   if(nsl*nend.gt.100000) then
      call get_phase_name(iph,1,phname)
      write(*,*)'3Y Limiting gridpoints in ',trim(phname),nend,30000
! I am not sure if nend is checked in the loops below, may cause segmentation
! fault
      nend=30000
   endif
   allocate(endm(nsl,nend),stat=errall)
   allocate(yendm(ncon,nend),stat=errall)
   allocate(yfra(ncon),stat=errall)
   allocate(ysave(ncon),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 9: ',errall
      gx%bmperr=4370; goto 1000
   endif
   yendm=1.0D-12
   do ij=1,ncon
      ysave(ij)=ydum(ij)
   enddo
! set endm(1..nsl,1) to first constituent in each sublattice
   do ij=1,nsl
      endm(ij,1)=incl(ij)
!      endm(ij,1)=incl(ij-1)+1
!      yendm(endm(ij,1),1)=one
   enddo
! loop to increment the constituents to generate all endmembers
! We should avoid all permutations according to BCC
! A:A:A:A 
! A:A:A:B   ingnore all permutations of B in sublattices
! A:A:A:C - A:A:A:X
! A:A:B:B - A:A:B:X
! A:A:C:C - A:A:X:X
! A:B:A:B - A:B:A:X <<<<<<< special for BCC
! A:C:A:C - A:X:A:X
! A:B:B:C - A:B:B:X
! A:B:C:C - A:B:X:X
! A:C:C:C - A:X:X:X
! B:B:B:B 
! B:B:B:C - B:B:B:X 
! etc
! X:X:X:X
! endm(1..nsl,1) set to first in all sublattices
! ordering always on the first 4 sublattices
! endm(1..nsl,jj) contains constituent indices in sublattice 1..nsl
! skip vacancies in the ordered sublattices ...
   kend=1
   sl1: do i1=1,nkl(1)
      endm(1,kend)=incl(0)+i1
      if(btest(splista(knr(endm(1,kend)))%status,SPVA)) then
!         write(*,*)'Skipping vacancies in ordered sublattices: ',&
!              trim(splista(knr(endm(1,kend)))%symbol),kend,knr(endm(1,kend))
         cycle sl1
      endif
      sl2: do i2=i1,nkl(2)
         endm(2,kend)=incl(1)+i2
         if(btest(splista(knr(endm(2,kend)))%status,SPVA)) cycle sl2
         sl3: do i3=i2,nkl(3) 
            endm(3,kend)=incl(2)+i3
            if(btest(splista(knr(endm(3,kend)))%status,SPVA)) cycle sl3
            sl4: do i4=i3,nkl(4)
               endm(4,kend)=incl(3)+i4
               if(btest(splista(knr(endm(4,kend)))%status,SPVA)) cycle sl4
               extra: if(nsl.gt.4) then
!                  write(*,16)'3Y endm 1: ',0,kend,(endm(ik,kend),ik=1,nsl)
                  rest: do ls=5,nsl
! Hm, problems to loop over constituents in sublattices nsl>4
                     sl5: do i5=1,nkl(ls)
!                        write(*,16)'3Y endm 2: ',i5,ls,kend,&
!                             (endm(ik,kend),ik=1,nsl),incl(ls)
                        if(endm(ls,kend).ge.incl(ls)) then
! reset the constiuent in sublattice ls to the first in this sublattice
                           endm(ls,kend)=incl(ls-1)+1
                        else
                           endm(ls,kend)=endm(ls,1)+1
                        endif
                        do ik=1,nsl
                           yendm(endm(ik,kend),kend)=one
                        enddo
!                        write(*,16)'3Y endm 3: ',i5,ls,kend,&
!                             (endm(ik,kend),ik=1,nsl),0
16                      format(a,3i5,4i4,2i7)
                        kend=kend+1
! this can be ver very many so maybe nend set to lower value above
                        if(kend.eq.nend) then
                           write(*,*)'3Y limiting grid in ',trim(phname),kend
                           goto 1000
                        endif
                        do iz=1,nsl
                           endm(iz,kend)=endm(iz,kend-1)
                        enddo
                     enddo sl5
                  enddo rest
               else
!                  write(*,16)'3Y yendm 3: ',kend,(endm(ik,kend),ik=1,nsl)
                  do ik=1,nsl
                     yendm(endm(ik,kend),kend)=one
                  enddo
!                  write(*,16)'3Y endm 4: ',0,ls,kend,&
!                       (endm(ik,kend),ik=1,nsl),0
                  kend=kend+1
! this can be ver very many so maybe nend set to lower value above
                  if(kend.eq.nend) then
                     write(*,*)'3Y limiting grid in ',trim(phname),kend
                     goto 1000
                  endif
                  do iz=1,nsl
                     endm(iz,kend)=endm(iz,kend-1)
                  enddo
               endif extra
            enddo sl4
         enddo sl3
      enddo sl2
   enddo sl1
   if(mode.eq.0 .and. test_phase_status_bit(iph,PHBORD)) then
! for BCC ordered phase add endmember with same constituents in first and third
! sublattices and loop in the others like A:B-X:A:B-X and B:C-X:B:C-X
!     write(*,*)'3Y Grid minimizer has no gridpoints for B32 ordering',kend-1
!      stop 'too many gridpoints'
   endif
! kend has been incremented one too much
   nend=kend-1
!   write(*,*)'3Y ordered endmemb: ',nend
!   if(mode.eq.0) then
! output adapted to 5 sublattices (interstitial)
!      if(nsl.eq.5) then
!         write(*,17)'3Y orded:',nend,((endm(ls,mend),ls=1,nsl),mend=1,nend)
17       format(a,i3,4(i4,4i3)/,(12x,i4,4i3,i4,4i3,i4,4i3,i4,4i3))
!      do i2=1,nend
!         write(*,18)i2,(endm(ls,i2),ls=1,nsl),(yendm(i1,i2),i1=1,nsl)
18       format('3Y yendm: ',i5,2x,4i3,2x,4F6.3)
!      enddo
!      elseif(nsl.eq.4) then
! output adapted to 4 sublattices
!         write(*,19)'3Y ordend: ',((endm(ls,mend),ls=1,nsl),mend=1,nend)
19       format(a,4(i4,3i3)/,(11x,i4,3i3,i4,3i3,i4,3i3,i4,3i3))
!      endif
!   endif
!
! copied from generic_grid_generator
!
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
!   dense=.FALSE.
!   gles=.FALSE.
!   defgrid=.TRUE.
   if(btest(globaldata%status,GSOGRID)) then
      gridlimit=3
   elseif(btest(globaldata%status,GSXGRID) .or. & 
        test_phase_status_bit(iph,PHXGRID)) then
      gridlimit=1
   else
      gridlimit=2
   endif
!   write(*,*)'3Y in generate_ordered_grid ',iph,nend,gridlimit,&
!        btest(globaldata%status,GSOGRID)
   ng=0
!-----------------------
   iiloop: do ii=1,nend
      ijs=1
      nendj=nend
      if(nend.ge.limits(gridlimit,1)) then
!      if(nend.eq.150 .or. (gles .and. nend.gt.40)) then
         nendj=ii
         ijs=ii
      endif
!      write(*,*)'3Y ii:',ii,ijs,nendj
      ijloop: do ij=ijs,nendj
         iks=1
         nendk=nend
         if(nend.ge.limits(gridlimit,2)) then
!         if(nend.gt.15 .or. (gles .and. nend.gt.13)) then
            nendk=ij
            iks=ij
         endif
         ikloop: do ik=iks,nendk
            ils=1
            nendl=nend
            if(nend.ge.limits(gridlimit,3)) then
!            if((nend.gt.12 .or. (gles .and. nend.gt.7)) then
               nendl=ik
               ils=ik
            endif
            illoop: do il=ils,nendl
!            illoop: do il=1,nendl
               ims=1
               nendm=nend
               if(nend.ge.limits(gridlimit,4)) then
!               if(nend.gt.7 .or. (gles .and. nend.gt.4)) then
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
                  if(mode.eq.0) then
!                     write(*,322)'3Y imloop: ',ng,ii,ij,ik,il,im
322                  format(a,i8,5i4)
!                     write(*,323)'3Y imloop1: ',ng,ii,ij,ik,il,im,0.0D0,yfra
323                  format(a,i5,5i3,': ',1pe12.4,0p20F5.2)
! strange bug in map3, maxng was zero sometimes ...
                     if(ng.gt.maxng) then
                        if(maxng.lt.100) then
                           write(*,*)'3Y max gripoints wrong 6: ',maxng,iph,mode
                        else
                           write(*,*)'3Y Too many gridpoints 6',ng,maxng,iph
                           gx%bmperr=4399; goto 1000
                        endif
                     endif
                     if(ng.gt.0 .and. mod(ng,30000).eq.0) then
                          lokph=phases(iph)
                          write(*,*)'3Y calculated ',ng,' gridpoints for ',&
                          trim(phlista(lokph)%name)
                       endif
                     call calc_gridpoint(iph,yfra,nrel,xarr(1,ng),garr(ng),ceq)
                     if(lutbug.gt.0) then
                        write(lutbug,710)'F: ',ngg,nrel,ncon,garr(ngg),&
                             (xarr(jj,ngg),jj=1,nrel),(yfra(jj),jj=1,ncon)
710                  format(a,i5,2i3,1pe10.2,10(0pF6.3))
                     endif
                     if(gx%bmperr.ne.0) then
                        write(*,*)'3Y error calculating gridpoint: ',&
                             iph,gx%bmperr
                        goto 1000
                     endif
!                     write(*,323)'3Y imloop2: ',ng,ii,ij,ik,il,im,garr(ng),yfra
                  elseif(mode.eq.ng) then
! when mode>0 we just want to know the constituent fractions
                     goto 900
                  endif
               enddo imloop
            enddo illoop
         enddo ikloop
      enddo ijloop
   enddo iiloop
!   do ii=1,ng
!      write(*,700)ii,garr(ii),(xarr(ij,ii),ij=1,nrel)
700   format('3Y gp: ',i5,1pe12.4,9(0pF6.3))
!   enddo
800 continue
   if(mode.gt.0) then
      write(*,*)'3Y could not find gridpoint ',mode,' in phase ',iph,ng
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
!
!
1000 continue
   if(mode.eq.0) then
! restore the composition
      call set_constitution(iph,1,ysave,qq,ceq)
   endif
! nothing done, just exit
1010 continue
   return
! dense gles
 end subroutine generate_fccord_grid

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine generate_charged_grid
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
   integer nend,ll,nsl,i1,i2,i3,loksp,mm,lokph,lokcs,np,nm,nn,ncc,iz,loopf,jj
   integer, dimension(:,:), allocatable :: neutral
   integer, dimension(10) :: gtype
!   integer, dimension(:), allocatable :: savengg
!   integer ielno(10)
!   double precision stoi(10),smass,qsp
   double precision charge,ratio1,ratio2
   double precision, dimension(:), allocatable :: y1,y2,y3,y4,y5
   real xdum(nrel),gdum
   integer, parameter :: ncf5=5,ncf3=3,alloneut=300000
   integer ncf,maxngg,ncon,maxgp1,errall
   integer, parameter :: maxgp2=10000,maxgp3=20000
! These are used to combine endmembers
   double precision, dimension(7), parameter :: nfact=&
        [0.01D0,0.1D0,0.33D0,0.51D0,0.67D0,0.9D0,0.99D0]
   double precision, dimension(ncf5), parameter :: cfact5=&
        [0.05D0,0.3D0,0.5D0,0.7D0,0.95D0]
   double precision, dimension(ncf3), parameter :: cfact3=&
        [0.1D0,0.5D0,0.9D0]
   logical single,endout,dense,skipped
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
!   endout=.true.
   endout=.FALSE.
!   skipped=.TRUE.
   skipped=.FALSE.
   if(endout) write(*,*)'3Y charged grid phase:',iph
!   ncf=ncf5
!   if(.not.allocated(savengg)) then
!      allocate(savengg(noofph))
!      savengg=0
!   endif
   maxngg=ngg
   ngg=0
   gtype=0
! get the phase data
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! Clement Introini bugfix
   lokph=phases(iph)
!
! I will handle this by first generate all endmembers with their charge
! and then try to combine them to get neutral gridpoints.
   nend=1
   inkl(0)=0
   do ll=1,nsl
      nend=nend*nkl(ll)
! inkl(ll) is the number of constituents up to and including sublattice ll
      inkl(ll)=inkl(ll-1)+nkl(ll)
   enddo
! ncon is the total number of constituents
   ncon=inkl(nsl)
!   write(*,*)'3Y Charged grid for phase ',iph,mode,nend,ncon
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
            if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
                 write(*,*)'3Y calculated ',ngg,' gridpoints for ',&
                 trim(phlista(lokph)%name)
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
!      write(*,*)'3Y Phase suspended as net charge: ',phlista(lokph)%name
! suspend all composition sets
      do mm=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(mm)
         ceq%phase_varres(lokcs)%phstate=PHSUS
      enddo
      goto 1000
   endif
!
   np=0
   nm=0
   nn=0
! Problem with CU2ZN1SN1S4 maybe because of sublattce with just VA ??
!   write(*,10)'3Y nend: ',iph,nend,0.0D0,(nkl(ll),ll=1,nsl)
10 format(a,2i4,5x,1pe12.4,10i3)
! allocate a record for each endmembers
   allocate(endmem(nend),stat=errall)
   allocate(endmem(1)%constit(nsl),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 10: ',errall
      gx%bmperr=4370; goto 1000
   endif
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
! a small allocation
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
! not using phase grid bit: test_phase_status_bit(iph,PHXGRID)
! default grid density with this big system
   if(mm.gt.30000) then
! very many components, minimum number of loops (no low density option)
      ncf=1; maxgp1=maxgp2
      if(btest(globaldata%status,GSXGRID) .or. &
           test_phase_status_bit(iph,PHXGRID)) then
! higher density requested
         ncf=ncf3; maxgp1=maxgp3
      elseif(btest(globaldata%status,GSYGRID)) then
! maximum density requested (may cause grid overflow ...)
         ncf=ncf5; maxgp1=maxgp3
      endif
   elseif(mm.gt.10000) then
! default grid density with a medium size system
      ncf=ncf3; maxgp1=maxgp2
      if(btest(globaldata%status,GSOGRID)) then
! lower density requested
         ncf=1
      elseif(btest(globaldata%status,GSXGRID) .or. &
           test_phase_status_bit(iph,PHXGRID)) then
! higher density requested
         ncf=ncf5
      elseif(btest(globaldata%status,GSYGRID))then
! maximum density requested
         ncf=ncf5; maxgp1=maxgp3
      endif
   else
! default grid for a small system
      ncf=ncf5; maxgp1=maxgp2
      if(btest(globaldata%status,GSOGRID)) then
! lower density requested
         ncf=ncf3
      elseif(btest(globaldata%status,GSXGRID) .or. &
           btest(globaldata%status,GSYGRID) .or. &
           test_phase_status_bit(iph,PHXGRID)) then
! higher density requested, this is maximum
         ncf=ncf5; maxgp1=maxgp3
      endif
   endif
! the statements below replaced by if statements above
!   if(.not.dense .and. mm.gt.5000) then
! maxgp1 used to skip some combinations
!      ncf=1; maxgp1=maxgp2
!   elseif(dense) then
! testing ...
!      ncf=ncf5; maxgp1=maxgp4
!   elseif(dense .or. mm.gt.2000) then
!      ncf=ncf3; maxgp1=maxgp3
!   else
! maximum dense grid
!      ncf=ncf5; maxgp1=maxgp3
!   endif
! debug output
!   if(mode.eq.0) then
!      write(*,'(a,10i7)')'3Y Generating charged grid: ',mode,iph,&
!           nend,mm,ncf,maxgp1
! noc2500 with just C1_MO2 ((12 * 2 * 2)=48 endmem) and GAS gives
!               level  C1_MO2  nend  mm       ncf  maxgpl   total gridpoints
! low density       0   22     48    21012    1    10000       11998
! default           1                         3    10000       33284
! high              2                         5    10000       42651
! maximum           3                         5    20000       57446
!
!   endif
! now calculate the number of gridpoints, consider single endmembers, 
! binary and ternary combinations in a triple loop
   np=0
   nn=0
!   if(mode.ge.0) then
! we have saved the number of gridpoints from the mode=-1 call here
!      np=savengg(iph)
!   write(*,*)'3Y allocate neutral: ',mode,alloneut
! guess a safe value ...
   allocate(neutral(alloneut,0:3),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 11: ',errall
      gx%bmperr=4370; goto 1000
   endif
   neutral=0
   np=0
   if(endout) write(*,*)'3Y starting generating grid in ionic solid phase',nend
   loop1: do i1=1,nend
      charge1A: if(endmem(i1)%charge.eq.zero) then
! first endmember neutral, one gridpoint
         np=np+1
         if(mode.ge.0) then
! for generating Y and G we save which endmembers to combine in neutral(*,0)
            neutral(np,0)=0
            neutral(np,1)=i1
         endif
         if(endout) write(*,298)'3Y generated 1 gp:  ',np,mode,1,0,i1,0,0
298      format(a,2i7,2i2,2x,3i4)
      endif charge1A
      gtype(1)=gtype(1)+1
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
               np=np+7
               if(endout) write(*,298)'3Y generated 7 gps: ',np,mode,3,1,i1,i2,0
               gtype(2)=gtype(2)+7
            else
!-----------------------------------------------------------------------
! second endmember has charge, a third endmember needed with opposite charge
               loop3A: do i3=i2+1,nend
                  if(endmem(i2)%charge*endmem(i3)%charge.lt.zero) then
! second and third endmembers have opposite charge, we have ncf gridpoints
! I1_n(I2_(1/c2)I3_(1/c3)_(1-n) where c2 is charge of i2 and c3 charge of i3
                     if(gtype(3).gt.maxgp1) then
                        if(skipped) &
                             write(*,*)'Skipping gridpoints type 3',gtype(3)
                        exit charge2A
                     endif
                     if(mode.ge.0) then
                        do ll=1,ncf
                           neutral(np+ll,0)=2
                           neutral(np+ll,1)=i1
                           neutral(np+ll,2)=i2
                           neutral(np+ll,3)=i3
                        enddo
                     endif
                     np=np+ncf
                     if(endout) write(*,298)'3Y generated 3A gps: ',&
                          np,mode,3,2,i1,i2,i3
                     gtype(3)=gtype(3)+ncf
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
                  if(gtype(4).gt.maxgp1) then
                     if(endout) write(*,*)'Skipping gridpoints type 4',gtype(4)
                     exit loop3B
                  endif
                  if(mode.ge.0) then
                     do ll=1,ncf
                        neutral(np+ll,0)=3
                        neutral(np+ll,1)=i1
                        neutral(np+ll,2)=i2
                        neutral(np+ll,3)=i3
                     enddo
                  endif
                  np=np+ncf
                  if(endout) write(*,298)'3Y generated 3B gps: ',&
                       np,mode,3,3,i1,i2,i3
                  gtype(4)=gtype(4)+ncf
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
            if(endout) write(*,298)'3Y generated 1 gp:  ',np,mode,1,4,i1,i2,0
            gtype(5)=gtype(5)+1
!-----------------------------------------------------------------------
            loop3C: do i3=i2+1,nend
               charge3A: if(endmem(i3)%charge.eq.zero) then
! third is neutral, we have ncf more gripoints
! at (I1_(1/c1)I2_(1/c2))_n(I3)_(1-n)
                  if(gtype(6).gt.maxgp1) then
                     if(skipped) write(*,*)'Skipping gridpoints type 6',gtype(6)
                     exit charge3A
                  endif
                  if(mode.ge.0) then
                     do ll=1,ncf
                        neutral(np+ll,0)=5
                        neutral(np+ll,1)=i1
                        neutral(np+ll,2)=i2
                        neutral(np+ll,3)=i3
                     enddo
                 endif
                 np=np+ncf
                 if(endout) write(*,298)'3Y generated 3C gps: ',&
                      np,mode,3,5,i1,i2,i3
                 gtype(6)=gtype(6)+ncf
               elseif(endmem(i1)%charge*endmem(i3)%charge.lt.zero) then
!-------------------------------------------------------------
! all 3 endmembers are charged, those of i2 and i3 have same sign, ncf gridp
! (I1_(1/c1)I2_(1/c2))_n(I1_(1/c1)I3_(1/c3))_(1-n)
                  if(gtype(7).gt.maxgp1) then
                     if(skipped) write(*,*)'Skipping gridpoints type 7',gtype(7)
                     exit charge3A
                  endif
                  if(mode.ge.0) then
                     do ll=1,ncf
                        neutral(np+ll,0)=6
                        neutral(np+ll,1)=i1
                        neutral(np+ll,2)=i2
                        neutral(np+ll,3)=i3
                     enddo
                  endif
                  np=np+ncf
                  if(endout) write(*,298)'3Y generated 3D gps: ',&
                       np,mode,3,6,i1,i2,i3
                  gtype(7)=gtype(7)+ncf
               elseif(gtype(8).lt.maxgp1) then
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
                  np=np+ncf
                  if(endout) write(*,298)'3Y generated 3E gps: ',&
                       np,mode,3,7,i1,i2,i3
                  gtype(8)=gtype(8)+ncf
               else
                  if(skipped) write(*,*)'Skipping gridpoints type 8,',gtype(8)
               endif charge3A
            enddo loop3C
!-----------------------------------------------------------------------
! first and second endmembers have charge with same sign
         elseif(gtype(9).lt.maxgp1) then
! we need a third endmember with opposite charge unless too many endmembers
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
                  np=np+ncf
                  if(endout) write(*,298)'3Y generated 3F gps: ',&
                       np,mode,3,8,i1,i2,i3
                  gtype(9)=gtype(9)+ncf
               endif
            enddo loop3D
            if(endout) write(*,777)'3Y loop3',gtype
         else
            if(skipped) write(*,*)'3Y skipping gridpoints type 9',gtype(9)
         endif charge1
         if(endout) write(*,777)'3Y loop2 ',gtype
      enddo loop2
      if(endout) write(*,777)'3Y loop1 ',gtype
777   format(a,10i6)
   enddo loop1
!=======================================================================
   if(endout) write(*,*)'3Y finished all loops for ionic phase: ',ngg
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
   allocate(y1(ncc),stat=errall)
   allocate(y2(ncc),stat=errall)
   allocate(y3(ncc),stat=errall)
   allocate(y4(ncc),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 12: ',errall
      gx%bmperr=4370; goto 1000
   endif
!   write(*,*)'3Y Charged grid: ',mode,nend,ncon,ncc
!   write(*,*)'3Y allocated neutral: ',mode,alloneut,np
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
! ncc and ncon should be identical, ny is returned as number of constituents
            ny=ncc
            do ll=1,ny
               yarr(ll)=y4(ll)
            enddo
!            write(*,507)'3Y Solution gp: ',mode,ny,y4
507         format(a,2i5,10F7.4)
            goto 1000
         endif
! continue searching for correct gridpoint of the solution
      else
         ngg=ngg+1
         if(ngg.gt.maxngg) then
            write(*,*)'3Y too many gripoints 3',ngg,maxngg,iph
            gx%bmperr=4399; goto 1000
         endif
         if(ngg.gt.0 .and. mod(ngg,30000).eq.0) &
              write(*,*)'3Y calculated ',ngg,' gridpoints for ',&
              trim(phlista(lokph)%name)
         call calc_gridpoint(iph,y4,nrel,xarr(1,ngg),garr(ngg),ceq)
         if(lutbug.gt.0) then
! ny, ncon, ncc ??
            write(lutbug,710)'I: ',ngg,nrel,ncon,garr(ngg),&
                 (xarr(jj,ngg),jj=1,nrel),(y4(jj),jj=1,ncon)
710         format(a,i5,2i3,1pe10.2,20(0pF6.3))
         endif
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

!\addtotable subroutine calc_gridpoint
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
   integer i,lokres,lokph
   double precision qq(5),xmol(nrel),ytemp(maxconst),gg,ss
   TYPE(gtp_phase_varres), pointer :: varres
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
!-------------------------------------------------
   if(globaldata%sysreal(1).gt.one) then
! if eec check if this is the liquid and if so select the maximum S ??
      varres=>ceq%phase_varres(lokres)
      lokph=varres%phlink
! value per mole component gres(2,1) is the -entropy NOTE ss is postive!!
      ss=-varres%gval(2,1)/qq(1)
      gg=varres%gval(1,1)
! I do not understand why iph and lokph is not the same!!
!      write(*,*)'3Y liqtest: ',iph,lokph
      eecheck: if(btest(phlista(lokph)%status1,PHLIQ)) then
!         write(*,*)'In calc_gridpoint liquid: ',sliqmax
!         if(sliqmax.gt.zero) exit eecheck
!         neecgrid=neecgrid+1
! note varres%gval(2,1) is -S !! Should we test max or min of liquid S ??
!         write(*,*)'Determining sliqmax: ',ss,sliqmax
         if(ss.gt.sliqmax) then
!            write(*,'(a,i3,5(1pe10.2))')'Liquid gridpoint: ',iph,&
!                 ss,sliqmin,gliqeec,varres%gval(1,1)/qq(1)
            sliqmax=ss; gliqeec=varres%gval(1,1)/qq(1)
         endif
! note sliqmax is - dg/dt and always positive.  It is the max entropy
! for the liquid at any gridpoint.  If the solid has higher entropy
! it should not be allowed to be stable.
      elseif(sliqmax.gt.zero) then
! this is a solid and we have a value for sliqmax for EEC to work better
!         write(*,*)'In calc_gridpoint solid: ',varres%gval(2,1)/qq(1),sliqmax
         if(ss.gt.sliqmax) then
! the solid s=-dG/dt is larger than sliqmax make the G more positive
! Note values are divided by RT.  Multiply with the number of atoms?
            varres%gval(1,1)=(gg+1.5D1)
!            varres%gval(1,1)=(gliqeec+one)
!            write(*,'(a,2i3,5(1pe10.2))')'3Y Solid corr: ',iph,lokres,&
!                 ss,sliqmax,gg,varres%gval(1,1),qq(1)
         endif
      else
         write(*,*)'3Y Gridminimizer has no Sliqmax for solid',iph
      endif eecheck
! list EEC values
!      write(*,11)'3Y EEC set: ',iph,xmol(1),ss,sliqmax,gg,varres%gval(1,1)
11    format(a,i3,F8.4,5(1pe12.4))
   endif
!--------------------------------------
   do i=1,nrel
      xarr(i)=real(xmol(i))
   enddo
! handle special problems
   if(qq(1).lt.5.0D-1) then
! number of real atoms less than 50%, a gridpoint with mainly vacancies ....
!      write(*,12)'3Y real atoms less than 0.5',lokres,qq(1),&
!           ceq%phase_varres(lokres)%gval(1,1)/qq(1)
12    format(a,i5,3(1pe12.4))
      gval=1.0E3
!      gval=max(1.0E2,real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)))
   elseif(abs(qq(2)).gt.1.0D-14) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! the gridpoint has net charge, qq(2), make gval more positive. 
! Note gval(1,1) is divided by RT so around -5<0
! There is special grid generator combining charged gripoints!!!!
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
!   write(*,12)'All gridpoints: ',lokres,qq(1),gval
!    read(*,20)ch1
20  format(a)
1000 continue
   return
 end subroutine calc_gridpoint

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine calcg_endmember
!\begin{verbatim}
 subroutine calcg_endmember(iphx,endmember,gval,ceq)
! calculates G for one mole of real atoms for a single end member
! used for reference states. Restores current composition (but not G or deriv)
! endmember contains indices in the constituent array, not species index
! one for each sublattice
! HERE G is divided by the number of atoms in the endmember
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
! to G, thus it is not sufficient just to calculate the G function, one must
! calculate TC etc.
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

!\addtotable subroutine calcg_endmemberx
!\begin{verbatim}
 subroutine calcg_endmemberx(iphx,endmember,gval,ceq)
! calculates G for single end member with current number of atoms
! used for reference states. Restores current composition (but not G or deriv)
! endmember contains indices in the constituent array, not species index
! one for each sublattice
! THIS ONE does not divide with the number of atoms
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
! to G, thus it is not sufficient just to calculate the G function, one must
! calculate TC etc.
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
! DO NOT DIVIDE WITH QQ
   gval=ceq%phase_varres(lokres)%gval(1,1)
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
 end subroutine calcg_endmemberx

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine calcg_endmember6
!\begin{verbatim} %-
 subroutine calcg_endmember6(iph,endmember,gval,ceq)
! calculates G AND ALL DERIVATEVS wrt T and P for one mole of real atoms
! for a single end member, used for reference states. 
! Restores current composition and G (but not deriv)
! endmember contains indices in the constituent array, not species index
! one for each sublattice
! THIS ONE returns 6 values: G, dG/dT; dG/dP; d2G/dT2; d2G/dTdP; d2G/dP2
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
! this is probably redundant ...
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

!\addtotable subroutine find_gridmin
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
   integer nrel1,nyp,griter,nopure(nrel),gpfail,evigloop
   double precision phfsave(nrel)
   integer, dimension(jerr) :: removed
   real gmin(nrel),dg,dgmin,gplan,gy,gvvp
! gridpoints that has less difference with the plane than this limit is ignored
   real, parameter :: dgminlim=1.0E-6
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
!770 format('Initial set of gridpoints: ',(/15i5))
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
   evigloop=0
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
!   if(trace) write(*,212)'3Y Found gridpoint ',nyp,inuse,dgmin,garr(nyp)
! evigloop happends when two gridpoints are exchanged
! uncomment the line below indicated gridmin can be improved ...
!   write(*,212)'3Y Found gridpoint ',nyp,inuse,dgmin,garr(nyp)
211 format(a,i7,1pe12.4,0pf7.4,6f7.4,(3x,10f7.4))
212 format(a,2i8,6(1pe11.3))
!-------------------------------------------------------------------------
! A case found when this seach never enden
   evigloop=evigloop+1
   if(evigloop.gt.500) then
      write(*,*)'3Y Gridmin gives up finding minimal set of gridpoints',evigloop
      goto 900 
   endif
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
771 format('3Y Final set of gridpoints: ',(/15i5))
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

!\addtotable subroutine merge_gridpoints
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
   integer i,ip,iph,jp,jump,kk,klast,kp,lokres,nm,jj,mj,lokph,j,npm
   integer notuse(nv),incy(nv)
   double precision ycheck(maxconst),qq(5),xerr(maxel),xfromy(maxel)
   double precision summu,sumam
   logical igen
   real xmix(maxel)
   double precision a1,a2,gdf,gval1,gval2,gval3,gval4,gval5,gmindif
   character phname*24
!
! gmindif is the value to accept to merge two gridpoints
! It should be a variable that can be set by the user for finetuning
!   write(*,*)'3Y Now we try to merge gridpoints in the same phase'
!   write(*,7)'3Y Merge_gridpoints is dissabled for the moment',nv
!7  format(a,i3)
! NOTE, always merge gripoints in ideal phases like gas
! UNFINISHED
   gmindif=ceq%gmindif
! used for testing 190603/BoS
!   gmindif=-1.0D-2
!   write(*,'(a,i3,1pe12.4)')'3Y Entering merge_gridpoints',nv,ceq%gmindif
!   goto 1100
!---------------------
   notuse=0
   nm=0
   npm=0
!   write(*,67)'3Y yl: ',(aphl(i),nyphl(i),i=1,nv)
67 format(a,20(F7.3,i4))
   incy(1)=1
   do i=2,nv
      incy(i)=incy(i-1)+nyphl(i-1)
   enddo
! start points of fractions for all gridpoints
!   write(*,68)'3Y ys: ',(incy(i),i=1,nv)
68 format(a,20i5)
   summu=zero
   xerr=zero
! constitution of solution gridpoints
!   do jp=1,nv
!      write(*,69)'3Y y:',(yphl(incy(jp)+i-1),i=1,nyphl(jp))
!   enddo
69 format(a,(12F6.2))
! this calculate the overall composition from gridpoints
   do jp=1,nv
      summu=summu+aphl(jp)
      do i=1,nrel
         xerr(i)=xerr(i)+aphl(jp)*xsol(i,jp)
      enddo
   enddo
!   write(*,73)'3Y in1: ',summu,(xerr(i),i=1,nrel)
73 format(a,F5.2,2x,9(f7.4))
!----------------------------------------------
100 continue
   igen=.false.
   firstgp: do jp=1,nv-1
!      write(*,*)'3Y notuse 1: ',jp,notuse(jp)
      if(notuse(jp).ne.0) cycle firstgp
      secondgp: do kp=jp+1,nv
!         write(*,*)'3Y notuse 2: ',kp,notuse(kp)
         if(notuse(kp).ne.0) cycle secondgp
         sameph: if(iphl(jp).eq.iphl(kp)) then
            gdf=zero
            iph=iphl(jp)
            lokph=phases(iph)
            if(btest(phlista(lokph)%status1,PHID)) then
! always merge gridpoints in ideal phases for example gas
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
!               write(*,830)'3Y not merged 9: ',jp,kp,gdf,iphl(jp),gmindif
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
!               write(*,830)'3Y not merged 1: ',jp,kp,gdf,iphl(jp),gmindif
               cycle secondgp
            else
! also compare between gval1 and gval3
               gdf=gval2-a1*gval1-a2*gval3
               if(gdf.gt.gmindif) then
! gval2 is s higher, no merge 1-2-3
!                  write(*,830)'3Y not merged 2: ',jp,kp,gdf,iphl(jp),gmindif
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
!               write(*,830)'3Y not merged 3',jp,kp,gdf,iphl(jp),gmindif
               cycle secondgp
            else
               gdf=gval4-a1*gval1-a2*gval5
               if(gdf.gt.gmindif) then
! gval4 is higer mo merge 1-4-5
!                  write(*,830)'3Y not merged 4',jp,kp,gdf,iphl(jp),gmindif
                  cycle secondgp
               endif
               gdf=gval4-a1*gval1-a2*gval5
               if(gdf.gt.gmindif) then
! gval4 is higer mo merge 1-4-5
!                  write(*,830)'3Y not merged 5',jp,kp,gdf,iphl(jp),gmindif
                  cycle secondgp
               endif
            endif
! compared 1-3-5, 1-2-5, 1-2-3, 1-4-5, 3-4-5, 2-3-5               
! in no case the middle point was above, that means merge
!--------------------------------------------- here we merge !!
200         continue
! gridpoint in ideal phase or point in between has lower G, merge
            call get_phase_name(iphl(jp),1,phname)
            if(gx%bmperr.ne.0) then
               phname='UNKNOWN'; gx%bmperr=0
            endif
!            write(*,830)'3Y merging:',jp,kp,gdf,aphl(jp),aphl(kp),trim(phname)
830         format(a,2i4,3(1pe12.4),' in ',a)
! If merging use correct phase amounts
            npm=npm+1
            a1=aphl(jp)/(aphl(jp)+aphl(kp))
            a2=aphl(kp)/(aphl(jp)+aphl(kp))
!            write(*,162)'3Y p1:',a2,(yphl(incy(jp)+j),j=0,nyphl(jp)-1)
!            write(*,162)'3Y p2:',a2,(yphl(incy(kp)+j),j=0,nyphl(kp)-1)
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
!            write(*,162)'3Y ym:',0.0D0,(yphl(incy(jp)+i),i=0,nyphl(jp)-1)
!            write(*,162)'3Y xy:',0.0D0,(xfromy(i),i=1,nrel)
! calculate mole fractions from xsol to compare!!
            do i=1,nrel
               xerr(i)=a1*xsol(jp,i)+a2*xsol(kp,i)
            enddo
!            write(*,162)'3Y xj+xk:',0.0D0,(xerr(i),i=1,nrel)
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
!            write(*,73)'3Y nu: ',summu,(xerr(jj),jj=1,nrel)
! the chemical potentials has changed but how?  Approximate the change by
! making gmindif more negative for each merge (does not affect ideal phases)
! I do not understand this but I keep it for the moment
!            gmindif=2.0D0*gmindif
            gmindif=1.2D0*gmindif
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
   if(npm.gt.0) write(*,'(a,i2,a)')'3Y Removed ',npm,' gridpoints by merging'
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

!\addtotable subroutine set_metastable_constitutions2
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
      igloop: do ig=ig1,ign
         if(garr(ig).ge.999.0) then
! gridpoints in phases with more than 50% vacancies have their garr(ig)=1.0D3
!            write(*,*)'Skipping gridpoint with too few atoms'
            cycle igloop
         endif
         gplan=zero
         do ie=1,nrel
            gplan=gplan+xarr(ie,ig)*cmu(ie)
         enddo
         dg=gplan-garr(ig)
         if(abs(dg).lt.abs(dgmin)) then
            ip=ig
            dgmin=dg
         endif
      enddo igloop
!      write(*,79)'3Y metastable: ',trim(phlista(iph)%name),iph,zph,&
!           ig1,ip,ign,dgmin
79    format(a,a,2i4,3i6,1pe12.4)
!      write(*,81)'3Y x: ',ip-nphl(zph),(xarr(ie,ip),ie=1,nrel)
!      write(*,81)'3Y x: ',ip-ig1,(xarr(ie,ip),ie=1,nrel)
81    format(a,i4,(10F6.3))
      if(ign.gt.ig1) then
! if ign=ig1 the phase has fixed constitution
! otherwise retrieve constitution for this gridpoint and insert it in phase
! we must provide mode and iph. The subroutine returns ny and yarr
! mode is the gridpoint in the phase
!         mode=ip-nphl(zph)
         mode=ip-ig1+1
! find the constitution of this gridpoint
!      call generate_grid(mode,iph,ign,nrel,xarr,garr,ny,yarr,gmax,ceq)
!         write(*,*)'3Y Get constitution of metastable phase ',iph,mode
         if(mode.gt.0) then
! this call returnes the constitution of gridpoint "mode"
! if mode=0 it generates the grid ... infinite loop
            call generic_grid_generator(mode,iph,ign,nrel,xarr,garr,&
                 ny,yarr,gmax,ceq)
            if(gx%bmperr.ne.0) then
               write(*,120)trim(phlista(iph)%name)
120            format('3Y Failed to set metastable constitution of ',a)
               gx%bmperr=0; cycle phloop
            endif
!            write(*,81)'3Y y: ',mode,(yarr(ie),ie=1,ny)
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

!\addtotable logical function global_equil_check1
!\begin{verbatim}
 logical function global_equil_check1(mode,addtuple,yfr,ceq)
! subroutine global_equil_check(ceq,newceq)
!
! This subroutine checks there are any gridpoints below the calculated solution
! if not it is taken as a correct global equilibrium
! This avoids creating any new composition sets but may fail in some cases
! to detect that the equilibrium is not global.
! mode=1 means try to recalculate equilibrium if not global (not implemented)
! if a gridpoint below is found addtuple and yfr returned with this
   implicit none
   integer mode,addtuple
   double precision, allocatable, dimension(:) :: yfr
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   TYPE(gtp_equilibrium_data), target :: cceq
   TYPE(gtp_equilibrium_data), pointer :: pceq
   logical global,newcs,notglobwarning1,notglobwarning2,wrongfrac,addgridpoint
   real, allocatable :: xarr(:,:),garr(:)
   real sumx
   double precision, dimension(maxconst) :: yarr
   double precision totmol,totmass,amount,gmax,dgmax,dgtest
   integer, allocatable :: kphl(:),iphx(:)
   integer gmode,iph,ngg,nrel,ny,ifri,firstpoint,sumng,nrph,ii,jj,nz,lokcs
   integer ics,pph,nyfas,gpz,iphz,nggz,errall,haha
   integer, parameter :: maxgrid=400000
!
!   write(*,*)'3Y In global_equil_check1',mode
   global=.TRUE.
   if(btest(globaldata%status,GSNOGLOB)) then
      write(*,*)'3Y Ignoring call to global_equil_check as global turned off!'
      goto 2000
   endif
   notglobwarning1=.TRUE.
   notglobwarning2=.TRUE.
   addgridpoint=.TRUE.
   if(mode.ne.1) addgridpoint=.FALSE.
   dgmax=zero
   addtuple=0
! Problem with invariant when mapping but not here
!   if(inveq(haha,ceq)) then
!      write(*,*)'3Y equilibrium is invariant when entering',haha
!   else
!      write(*,*)'3Y equilibrium is not invariant when entering',haha
!   endif
! COPY the whole equilibrium record to avoid destroying anything!!
! otherwise I had strange problems with amounts of phases ??
   cceq=ceq
   pceq=>cceq
   nrph=noofph
   allocate(kphl(0:nrph+1),stat=errall)
   allocate(iphx(nrph+1),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 13: ',errall
      gx%bmperr=4370; goto 1000
   endif
!
   sumng=0
   ifri=0
   firstpoint=1
   iphx=0
   kphl=0
   ifri=0
   pph=0
   wrongfrac=.true.
   ggloop: do iph=1,nrph
! include all phases with any composition set entered (but only once!)
      do ics=1,noofcs(iph)
! new: -3 suspended, -2 dormant, -1,0,1 entered, 2 fixed
! ignore phases whith no composition set entered
! If a phase+compset FIX one should never be here as conditions wrong
         if(test_phase_status(iph,ics,amount,pceq).gt.PHDORM) then
            pph=pph+1
            iphx(pph)=iph
            cycle ggloop
         endif
      enddo
   enddo ggloop
!
   nrel=noofel
!   write(*,11)'3Y gpa:',pph,(iphx(iph),iph=1,pph)
! allocate arrays, added 1 to avoid a segmenentation fault ....
   allocate(xarr(nrel,maxgrid),stat=errall)
   allocate(garr(maxgrid),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 14: ',errall
      gx%bmperr=4370; goto 1000
   endif
! calculate the composition and G for the gridpoints
   ii=1
   loop2: do ifri=1,pph
      ngg=maxgrid-ii
!      write(*,10)'3Y calling generic_grid 2: ',ifri,iphx(ifri),ngg,pph
10    format(a,2i5,2i10,5i5)
!>>>>>>> important: changes here must be made also in global_gridmin
      if(btest(globaldata%status,GSOGRID)) then
! The possibility to use the old grid tested
         call generate_grid(0,iphx(ifri),ngg,nrel,xarr(1,ii),&
              garr(ii),ny,yarr,gmax,pceq)
      else
         call generic_grid_generator(0,iphx(ifri),ngg,nrel,xarr(1,ii),&
              garr(ii),ny,yarr,gmax,pceq)
      endif
!>>>>>>>> impportant end!
!      write(*,*)'3Y Back from grid generator: ',ifri,iphx(ifri),ngg
      if(gx%bmperr.ne.0) goto 1000
      kphl(ifri)=kphl(ifri-1)+ngg
      ii=kphl(ifri)+1
   enddo loop2
   sumng=kphl(pph)
!   write(*,11)'3Y gpc:',(kphl(iph),iph=1,nrph)
11 format(a,10i7/(7x,10i7))
!   write(*,*)'3Y Calculated ',sumng,' gridpoints for check.',kphl(0)
! We have calculated sumng gripoints in pph phases
! check if any gridpoint is below the G surface defined by cmuval
   iph=0
   nyfas=0
   iphz=0
   loop4: do ifri=1,sumng
! keep track of the phase the gridpoint belongs to
      if(ifri.gt.nyfas) then
! iph is the phase index in phasetuple (and phases)
         iph=iph+1
!         ny=ny+kphl(iph)
         nyfas=kphl(iph)
      endif
      gmax=zero
      sumx=0.0E0
      do ngg=1,nrel
         gmax=gmax+dble(xarr(ngg,ifri))*pceq%cmuval(ngg)
         sumx=sumx+xarr(ngg,ifri)
      enddo
!      if(ifri.eq.sumng) write(*,*)'3Y OK ',ifri,iph
      if(abs(sumx-1.0E0).gt.1.0E-4) then
         cycle loop4
      endif
!      write(*,75)'3Y check: ',ifri,iphx(iph),garr(ifri),gmax,garr(ifri)-gmax
75    format(a,i6,i4,5(1pe12.4))
      dgtest=gmax-dble(garr(ifri))
      stableornot: if(dgtest.gt.1.0D-4*abs(gmax)) then
!      stableornot: if(dgtest.gt.1.0D-7*abs(gmax)) then
!      if(dgtest.gt.dgmax) then
!------------------------------------------------------------------
!         write(*,76)'3Y gridpoint below G surface: ',ifri,iph,iphx(iph),&
!              dgtest,1.0D-4*dgmax
76       format(a,i7,2i4,2(1pe12.4))
! if the phase is stoichiometric and stable this is a rounding off problem
! find the phase record using the phase tuple
         lokcs=phasetuple(iph)%lokvares
         nz=size(pceq%phase_varres(lokcs)%sites)-&
              size(pceq%phase_varres(lokcs)%yfr)
         if(nz.eq.0) then
            if(pceq%phase_varres(lokcs)%phstate.ge.PHENTSTAB) cycle loop4
! if number of constituent fractions equal to sublattice the composition is fix
! If this is a test at a node point we may have an allotropic phase whicj is
! stable, then the driving force should be small ... check if dgm is very small
            write(*,'(a,i5,F10.2,2(1pe12.4))')'3Y allotrop DGM: ',&
                 lokcs,pceq%tpval(1),pceq%phase_varres(lokcs)%dgm
            if(pceq%phase_varres(lokcs)%dgm.lt.2.0D-1) cycle loop4
         endif
! This phase should be stable, maybe there are others?
         if(dgtest.lt.dgmax) cycle loop4
         dgmax=dgtest
! This gridpoint is the currently lowest below the current G plane
!         write(kou,77)ifri,iph,iphx(iph),trim(phlista(phases(iphx(iph)))%name)
77       format('3Y found a stable gridpoint: ',3i5,' in ',a)
         global=.FALSE.
         gpz=ifri
!         iphz=iph
         iphz=iphx(iph)
         nggz=kphl(iph-1)
!         write(*,*)'3Y saving most stable gp: ',gpz,iphz
      endif stableornot
!      if(ifri.eq.sumng) write(*,*)'OK ',ifri
   enddo loop4
! no gridpoint below current G surface
!   write(*,*)'3Y finished loop4',iphz,global
   goto 1000
! Found gridpoint below gmax, if mode=/=1 just return error message
500 continue
   write(*,*)'3Y Sorry I have not yet implemented automatic recalculation!'
   if(mode.eq.1) then
! Here we try to recalculate the equilibrium with a new phase stable
      continue
   else
      write(*,*)'3Y Please include this phase and recalculate equilibrium'
   endif
!
1000 continue
!1010 continue
!   write(*,*)'3Y global_equil_check label 1000',global,gx%bmperr
! set the error code here so we can finish this routine
   if(.not.global) then
!      write(*,1011)'3Y most stable gridpoint: ',gpz,nggz,iphz,dgmax
1011  format(a,2i7,i3,1pe12.4)
      addornot: if(addgridpoint) then
! Add this gripoint as entered and recalculate
! extract constitution, ny=-100 to get some output ..
!         write(*,*)'3Y Trying to extract constitution, ngg:'
!         ny=-100
         if(btest(globaldata%status,GSOGRID)) then
! we do not have ifri and iphx here
            call generate_grid(gpz-nggz,iphz,nggz,nrel,xarr,&
                 garr,ny,yarr,gmax,pceq)
         else
            call generic_grid_generator(gpz-nggz,iphz,nggz,nrel,&
                 xarr,garr,ny,yarr,gmax,pceq)
         endif
         if(ny.gt.0) then
!            write(*,83)'3Y gpy: ',ny,(yarr(ngg),ngg=1,ny)
83          format(a,i7,9F7.4,(8x,14F5.2))
! a small allocate
            allocate(yfr(ny))
            do ngg=1,ny
               yfr(ngg)=yarr(ngg)
            enddo
         else
            write(*,*)'3Y Failed extract constitution',ny
         endif
      else
! This gridpoint is for a phase that is not stable but has a stable grid point
! but we will not try to recalculate
         if(notglobwarning1) then
! write this once only
            write(kou,87)trim(phlista(phases(iphz))%name),pceq%tpval(1),&
                 (xarr(ngg,gpz),ngg=1,nrel)
87          format(/' *** Gridtest found equilibrium not global, ',a,&
                 ' is stable at T=',F8.2/5x,'with mole fractions:'/(1x,13F6.3))
            notglobwarning1=.FALSE.
         endif
      endif addornot
      addtuple=iphz
      gx%bmperr=4352
   endif
!   write(*,*)'3Y Deallocating, check due to segmentation fault ...'
   if(allocated(xarr)) then
      deallocate(xarr)
      deallocate(garr)
      deallocate(kphl)
      deallocate(iphx)
   endif
2000 continue
!   if(inveq(haha,ceq)) then
!      write(*,*)'3Y equilibrium is invariant when at exit',haha
!   else
!      write(*,*)'3Y equilibrium is not invariant when at exit',haha
!   endif
   global_equil_check1=global
   return
 end function global_equil_check1

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine check_all_phases
!\begin{verbatim}
 subroutine check_all_phases(mode,ceq)
!
! This function check for all phases if there are any gridpoints
! closer (or below) to the current calculated solution
! if so it changes the composition of the phase
! If a gridpoint is BELOW the current plane an error code is returned
! phase should be stable with another composition an error code is returned
! It does not creating any new composition sets
! It can be usd during STEP/MAP to update compositions of metastable
! phases which have become stuck in a local minimium
! if error 4365 or 4364 is set mode will return index in meqrec%phr 
! of the phase that should be stable
   implicit none
   integer mode
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   TYPE(gtp_equilibrium_data), target :: cceq
   TYPE(gtp_equilibrium_data), pointer :: pceq
   integer iph,phstat,saverr
!
!   write(*,*)'3Y In check_all_phases'
! COPY the whole equilibrium record to avoid destroying anything!!
! otherwise I had strange problems with amounts of phases ??
   saverr=0
   cceq=ceq
   pceq=>cceq
! mode will be updated inside check_phase_grid to correspond to meqrec%phr index
   mode=0
   ggloop: do iph=1,noofph
! include all phases with any composition set entered (but only once!)
! loop for composition sets inside check_phase as they all have the same grid
      call check_phase_grid(iph,mode,pceq,ceq)
      if(gx%bmperr.ne.0) then
! if a stable phase need a new composition terminate and return error
         if(gx%bmperr.eq.4366) then
! grid minimizer needed to create new composition set is needed
            write(*,*)'3Y New composition set needed: ',gx%bmperr,mode,mode
            goto 1000
         elseif(gx%bmperr.eq.4365) then
! new stable phase composition inserted in unstable composition set
! go back and take halfstep in step/map
            write(*,*)'3Y found stable phase: ',gx%bmperr,mode,mode
            goto 1000
            saverr=gx%bmperr
         endif
         gx%bmperr=0
      endif
   enddo ggloop
   if(saverr.ne.0) gx%bmperr=saverr
!
1000 continue
   return
 end subroutine check_all_phases

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine check_phase_grid
!\begin{verbatim}
 subroutine check_phase_grid(iph,jcs,pceq,ceq)
!
! This function check for A SINGLE PHASE if there are any gridpoints
! closer (or below) to the current calculated solution if so it
! changes the composition of the phase If a gridpoint is below the
! phase should be stable with another composition an error code is
! returned It does not creating any new composition sets but may fail
! It can be usd during STEP/MAP to update compositions of metastable
! phases which have become stuck in a local minimium
! NOTE pceq is a pointer to a copy of the real equilibrium record
! ceq is a pointer to the real equilibrium record
! jcs is returned as the composition set that should be stable (if any)
   implicit none
   integer iph,jcs
   TYPE(gtp_equilibrium_data), pointer :: ceq,pceq
!\end{verbatim}
   real, allocatable :: xarr(:,:),garr(:)
   double precision, dimension(maxel) :: x1mol,wmass
   double precision, dimension(maxconst) :: yarr
   double precision totmol,totmass,amount,gmax,dgmax,dgtest
   double precision, parameter :: mindg=1.0D-6
! max 9 composition sets
   double precision gorig,gbest,gdiff,gset(9),gplan,am,qq(5)
! for debugg
   double precision yold(100)
   integer, allocatable :: kphl(:),iphx(:)
   integer ii,jj,kk,nrel,lokcs,moded,ny,ics,ics2,ncs,stcs(4),nstcs,ie,ngg
   integer phstat,lokph,lokres,errall
   logical skip
   integer, parameter :: maxgrid=100000
!
!   write(*,*)'3Y In check_phase_grid: ',iph
   nrel=noofel
   moded=0
! allocate arrays
   allocate(xarr(nrel,maxgrid),stat=errall)
   allocate(garr(maxgrid),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 15: ',errall
      gx%bmperr=4370; goto 1000
   endif
   gset=-one
   skip=.TRUE.
! loop for all composition sets
   ncs=noofcs(iph)
   stcs=0
   nstcs=0
   gloop: do ics=1,ncs
! calculate G for current composition, ignore dormant and suspended sets
      phstat=test_phase_status(iph,ics,amount,pceq)
      if(phstat.lt.PHDORM) cycle gloop
      skip=.FALSE.
      call calcg(iph,ics,moded,lokcs,ceq)
      if(gx%bmperr.ne.0) goto 1000
      call calc_phase_molmass(iph,ics,x1mol,wmass,totmol,totmass,am,ceq)
      if(gx%bmperr.ne.0) goto 1000
! abnorm(1) is number of atoms per formula unit
      gorig=pceq%phase_varres(lokcs)%gval(1,1)/&
           pceq%phase_varres(lokcs)%abnorm(1)
! calculate the difference with the current stable tangent plane
! It can be zero if the composition set is stable
      gplan=zero
      do ii=1,nrel
         gplan=gplan+x1mol(ii)*pceq%cmuval(ii)
      enddo
! this is the original drivining force for each composition set
      gset(ics)=gorig-gplan
! there can be more than one stable composition set ... fix another time ...
      if(phstat.ge.PHENTSTAB) then
         if(nstcs.ge.4) then
!            write(*,*)'More than 4 stable composition sets of phase',iph
            gx%bmperr=4399; goto 1000
         endif
!         write(*,*)'Stable phase and set: ',iph,ics,gset(ics)
         nstcs=nstcs+1; stcs(nstcs)=ics
      endif
   enddo gloop
! all composition sets suspended or dormant in this phase
   if(skip) goto 1000
!   write(*,20)'3Y phase grid: ',nstcs,ncs,(gset(ii),ii=1,ncs)
!20 format(a,2i2,6(1pe12.4))
!
! now calculate the gridpoints, composition and G
   ngg=maxgrid
!   write(*,*)'Calculate grid for phase ',iph,nrel,ngg
   call generic_grid_generator(0,iph,ngg,nrel,xarr,garr,ny,yarr,gmax,pceq)
!   write(*,*)'3Y error & grid: ',gx%bmperr,ngg
   if(gx%bmperr.ne.0) goto 1000
! loop through all gridpoints to find one closesed to the stable tangent plane
! ngg set to number of real gridpoints
! note mixed single and double precision but that is OK
   gbest=-1.0D3
   do ii=1,ngg
      gplan=zero
      do jj=1,nrel
         gplan=gplan+xarr(jj,ii)*pceq%cmuval(jj)
      enddo
! note gdiff should be negative if metetstable
      gdiff=gplan-garr(ii)
!      write(*,22)'3Y GRID: ',ii,gdiff,gbest,garr(ii),gplan
22    format(a,i4,4(1pe12.4))
      if(gdiff.gt.gbest) then
         kk=ii; gbest=gdiff
      endif
   enddo
!
!   write(*,30)'3Y Gridpoint ',kk,gbest,(xarr(ii,kk),ii=1,nrel)
30 format(a,i4,e12.4,10(F8.5))
! now we compare the best gridpoint with the composition sets
   loop1: do ics=1,ncs
! jcs will be the correct phase index inside meqrec%phr array ??
      jcs=jcs+1
      phstat=test_phase_status(iph,ics,amount,pceq)
      if(phstat.lt.PHDORM) cycle loop1
! extract constitution for the best gridpoint kk
      call generic_grid_generator(kk,iph,ngg,nrel,xarr,garr,&
           ny,yarr,gmax,pceq)
      if(gbest.ge.mindg) then
! there is a gridpoint below the tangent plane
! If there is a metastable composition set copy the gripoint constitution
! to that and recalculate.  If no free composition set test if the grid
! point can be merged with a stable composition set.  If not recalculate
! with grid minimizer
         if(nstcs.gt.0) then
! There is one or more stable composition set, if there is an unstable one
! then set the gridpoint constitution in that
            if(nstcs.eq.ncs) then
! All composition sets already stable! 
! we have to compare if the G curve between the gridpoint and all the
! composition sets is convex or concave.  NOT IMPLEMENTED
!               loop2: do ics2=1,nstcs
!                  if(stcs(ics2).lt.0) cycle loop2
!                  call calc_phase_molmass(iph,stcs(ics2),x1mol,wmass,totmol,&
!                       totmass,am,pceq)
!                  if(gx%bmperr.ne.0) goto 1000
! we should check here if there is a maximum G between the gridpoint and the
! stable composition set.  To be done ...
!                  write(*,*)'3Y New composition set needed'
!                  gx%bmperr=4365; goto 1000
!               enddo loop2
! we arrive here if we could not merge gridpoint with a stable composition set
! the error code demand a global grid minimization.
               write(*,*)'3Y New composition set needed for:',iph,ncs,nstcs
!               write(*,90)1,iph,ics,ceq%tpval(1),gbest,gset(ics)
               call set_constitution(iph,ics,yarr,qq,ceq)
               if(gx%bmperr.ne.0) goto 1000
               gx%bmperr=4365; goto 1000
            elseif(gset(ics).lt.zero) then
! there is at least one unstable composition set, check if gset(ics)<0
! and insert the GRIDPOINT constitution if gset(ics) negative
! This composition set is not stable, insert stable gridpoint constitution
               write(*,90)2,iph,ics,ceq%tpval(1),gbest,gset(ics),4365
90             format('3Y stable gridpoint ',i1,2x,2i4,F10.2,2(1pe12.4),i5)
! Check old constitution
!               call get_phase_compset(iph,ics,lokph,lokres)
!               write(*,95)'3Y oldy: ',ceq%phase_varres(lokres)%yfr
!               write(*,95)'3Y newy: ',(yarr(ii),ii=1,ny)
               call set_constitution(iph,ics,yarr,qq,ceq)
               if(gx%bmperr.ne.0) goto 1000
! this error code demand recalculation without grid minimizer
               gx%bmperr=4365; goto 1000
            endif
         else
! There are no stable composition sets, we can set the stable gridpoint
! constitution in this composition set and request a new equilibrium calculation
! NOTE we use ceq pointer to set yarr in original record
            write(*,90)3,iph,ics,ceq%tpval(1),gbest,gset(ics),4365
! Check old constitution
!            call get_phase_compset(iph,ics,lokph,lokres)
!            write(*,95)'3Y oldy: ',ceq%phase_varres(lokres)%yfr
!            write(*,95)'3Y newy: ',(yarr(ii),ii=1,ny)
95          format(a,10(F7.4))
            call set_constitution(iph,ics,yarr,qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gx%bmperr=4365; goto 1000
         endif
      elseif(gbest.gt.gset(ics)) then
! SKIP THIS FOR THE MOMENT
! The best gridpoint is not stable but it is closer to tangent plane than this
! composition set WHICH THUS MUST BE UNSTABLE!.
! This change can avoid a phase is stuck in a local minimum
! BUT if another composition set is stable do not change because it the
! gridpoint is probably close to the stable composition.
         if(nstcs.eq.0) then
            write(*,92)iph,ics,ceq%tpval(1),gbest,gset(ics)
92          format('3Y better gridpoint in ',i4,i2,F10.2,2(1pe12.4))
            call set_constitution(iph,ics,yarr,qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
         endif
! This do not require a new calculation
!      else
! Nothing to do as the best gridpoint is further away from tangent plane
! than this metastable composition set
      endif
   enddo loop1
! The allocated arrays should deallocate by themselves
1000 continue
   return
 end subroutine check_phase_grid

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine separate_constitutions
!\begin{verbatim}
 subroutine separate_constitutions(ceq)
! This is called during step/map
! Go through all entered phases and if there are two composition sets
! that have similar constitutions then separate them
! Used during mapping of for example Fe-Cr to detect the miscibility gap
    implicit none
    TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer phtup,nextset,lokcs1,lokcs2,ic,ll,lokph,ss,ts
! max 9 sublattices
    integer iymin(9),iymax(9),qq
    double precision ymax(9),ymin(9),ysame
    double precision, allocatable :: yarr(:)
!    write(*,*)'3Y Check if two composition sets are same: ',ceq%tpval(1)
    allph: do phtup=1,nooftup()
       nextset=phasetuple(phtup)%nextcs
       if(nextset.le.0) cycle allph
       lokph=phasetuple(phtup)%lokph
       lokcs1=phasetuple(phtup)%lokvares
       lokcs2=phasetuple(nextset)%lokvares
       ymin=one
       ymax=zero
       iymin=0
       iymax=0
!       ts=ceq%phase_varres(lokcs1)%tnoofr
       ts=phlista(lokph)%tnooffr
       ll=1
       qq=phlista(lokph)%nooffr(1)
       do ic=1,ts
          ysame=ceq%phase_varres(lokcs1)%yfr(ic)
          if(abs(ysame-ceq%phase_varres(lokcs2)%yfr(ic)).gt.1.0D-2) then
!             write(*,66)'3Y two compsets not same:',lokcs1,lokcs2,ceq%tpval(1)
!          write(*,77)'3Y d1:',lokcs1,(ceq%phase_varres(lokcs1)%yfr(ss),ss=1,ts)
!          write(*,77)'3Y d2:',lokcs2,(ceq%phase_varres(lokcs2)%yfr(ss),ss=1,ts)
             cycle allph
          else
! map8 gave segmentation fault here, fixed ??
!             write(*,10)'3Y qq: ',phtup,nextset,ic,qq,ll,ts,&
!                  phlista(lokph)%nooffr(ll),size(phlista(lokph)%nooffr),&
!                  size(ceq%phase_varres(lokcs2)%yfr)
!10           format(a,10i5)
             if(ic.gt.qq) then
                ll=ll+1
                qq=qq+phlista(lokph)%nooffr(ll)
             endif
             if(ysame.lt.ymin(ll)) then
                iymin(ll)=ic; ymin(ll)=ysame
             endif
             if(ysame.gt.ymax(ll)) then
                iymax(ll)=ic; ymax(ll)=ysame
             endif
          endif
       enddo
! These two composition sets have identical compositions, skip if both stable
!       write(*,66)'3Y two compsets same:',lokcs1,lokcs2,ceq%tpval(1),&
!            (ymin(ic),ymax(ic),ic=1,ll)
66     format(a,2i3,f8.2,2x,(8F6.3))
       if(ceq%phase_varres(lokcs1)%phstate.ge.PHENTSTAB) then
          if(ceq%phase_varres(lokcs2)%phstate.ge.PHENTSTAB) then
!             write(*,*)'Wow, two identical phases stable!',lokcs1,lokcs2
             cycle allph
          endif
!          write(*,77)'3Y s1:',lokcs1,(ceq%phase_varres(lokcs1)%yfr(ss),ss=1,ts)
! set the constitution of lokcs2 to one-the stable
! or maybe to its default??
          lokph=phasetuple(phtup)%lokph
          ic=0
          phsubl1: do ll=1,phlista(lokph)%noofsubl
             if(phlista(lokph)%nooffr(ll).eq.1) then
                ic=ic+1; cycle phsubl1
             endif
             ysame=0.1/real(phlista(lokph)%nooffr(ll))
             do ss=1,phlista(lokph)%nooffr(ll)
                ic=ic+1
                if(ic.eq.iymin(ll)) then
                   ceq%phase_varres(lokcs2)%yfr(ic)=0.9
                else
                   ceq%phase_varres(lokcs2)%yfr(ic)=ysame
                endif
             enddo
          enddo phsubl1
!          write(*,77)'3Y s1:',lokcs1,(ceq%phase_varres(lokcs1)%yfr(ss),ss=1,ts)
!          write(*,77)'3Y s2:',lokcs2,(ceq%phase_varres(lokcs2)%yfr(ss),ss=1,ts)
77        format(a,i3,8F6.3)
       else
! lokcs1 is not stable, change its constitution away from lokcs2
!       elseif(ceq%phase_varres(lokcs2)%phstate.ge.PHENTSTAB) then
! set the constitution of lokcs1 to one-the lokcs2
! or maybe to its default??
          lokph=phasetuple(phtup)%lokph
          ic=0
          phsubl2: do ll=1,phlista(lokph)%noofsubl
             if(phlista(lokph)%nooffr(ll).eq.1) then
                ic=ic+1; cycle phsubl2
             endif
! very strange, if I divide with real(phlista(lokph)%nooffr(ll)-1) 
! the metastable exrapolation is still there !!
             ysame=0.1/real(phlista(lokph)%nooffr(ll))
             do ss=1,phlista(lokph)%nooffr(ll)
                ic=ic+1
                if(ic.eq.iymin(ll)) then
                   ceq%phase_varres(lokcs1)%yfr(ic)=0.9
                else
                   ceq%phase_varres(lokcs1)%yfr(ic)=ysame
                endif
             enddo
          enddo phsubl2
!          write(*,77)'3Y z1:',lokcs1,(ceq%phase_varres(lokcs1)%yfr(ss),ss=1,ts)
!          write(*,77)'3Y z2:',lokcs2,(ceq%phase_varres(lokcs2)%yfr(ss),ss=1,ts)
!       else
!          write(*,*)'Both compsets unstable'
       endif
    enddo allph
1000 continue
  end subroutine separate_constitutions

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function allotropes
!\begin{verbatim}
  logical function allotropes(irem,iadd,iter,ceq)
! This function return TRUE if the phases indicated by IREM and IADD both have
! fixed and identical composition, i.e. they are componds and allotropes
! Such a transition can cause problems during a STEP command.
    implicit none
    TYPE(gtp_equilibrium_data), pointer :: ceq
    integer iadd,irem,iter
!\end{verbatim} %+
    integer lokph1,lokph2,nofr,jj
    double precision x1mol(maxel),x2mol(maxel),wmass(maxel),totmol,totmass,am
    logical allo
    allo=.false.
    goto 1000
!    write(*,*)'3A checking allotropes',irem,iadd
    write(*,10)iter,trim(phlista(phases(irem))%name),&
         trim(phlista(phases(iadd))%name)
10  format('3A checking allotropes',i5,2x,a,2x,a)
    lokph1=phases(irem)
    lokph2=phases(iadd)
! spurious segmentation faults here ...
    if(lokph1.le.0 .or. lokph2.le.0) then
! composition ses created during mapping are not included in phases array ??
       write(*,*)'3A error checking allotropes: ',lokph1,lokph2
       goto 1000
    endif
! check if both have fixed composition
    if(phlista(lokph1)%noofsubl-phlista(lokph1)%tnooffr.eq.0 .and. &
         phlista(lokph2)%noofsubl-phlista(lokph2)%tnooffr.eq.0) then
! they have fixed composition but can be modelled differently
! we have to calculate their mole fractions ...
       call calc_phase_molmass(irem,1,x1mol,wmass,totmol,totmass,am,ceq)
       call calc_phase_molmass(irem,1,x2mol,wmass,totmol,totmass,am,ceq)
       if(gx%bmperr.ne.0) goto 1000
       do jj=1,noofel
          if(abs(x1mol(jj)-x2mol(jj)).gt.1.0D-6) exit
       enddo
! Fortran standard says jj>noofel if loop finish without exit
       if(jj.gt.noofel) then
          allo=.true.
!          write(*,*)'The phases are allotropes!',ceq%tpval(1)
       endif
    endif
1000 continue
    allotropes=allo
    return
  end function allotropes

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function same_stoik
!\begin{verbatim}
 logical function same_stoik(iph,jph)
! MAYBE IDENTICAL TO ALLOTROPES ?
! return TRUE if phase iph and jph are both stoichiometric and have the
! same composition  Used to check when adding a phase during equilibrium
! calculation as it normally fails to have two such phases stable
! iph and jph are phase tuple indices
   implicit none
   integer iph,jph
!\end{verbatim} %+
   integer loki,lokj,ll,kk
   logical same
!
   same=.false.
! iph and jph can be second or later composition sets
!   write(*,*)'3F same_stoik 1: ',iph,jph,&
!        phasetuple(iph)%lokph,phasetuple(jph)%lokph
!   loki=phases(iph); lokj=phases(jph)
   if(iph.le.0 .or. iph.gt.nooftup() .or. jph.le.0 .or.jph.gt.nooftup()) then
      write(*,*)'Calling same_stoik with illegal arguments ',iph,jph
      gx%bmperr=4399; goto 1000
   endif
   loki=phasetuple(iph)%lokph; lokj=phasetuple(jph)%lokph
   if(.not.btest(phlista(loki)%status1,PHNOCV)) goto 1000
   if(.not.btest(phlista(lokj)%status1,PHNOCV)) goto 1000
   if(phlista(loki)%noofsubl.ne.phlista(lokj)%noofsubl) goto 1000
   kk=0
   do ll=1,phlista(loki)%noofsubl
      if(firsteq%phase_varres(phlista(loki)%linktocs(1))%sites(ll).ne.&
         firsteq%phase_varres(phlista(lokj)%linktocs(1))%sites(ll)) goto 1000
      kk=kk+1
      if(phlista(loki)%constitlist(kk).ne.&
         phlista(lokj)%constitlist(kk)) goto 1000
   enddo
   same=.true.
1000 continue
   same_stoik=same
   return
 end function same_stoik

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable logical function fixedcomposition
!\begin{verbatim}
 logical function fixedcomposition(iph)
! returns TRUE if phase cannot vary its composition
   integer iph
!\end{verbatim}
   integer lokph
   lokph=phases(iph)
! Wow a bug! using iph instead of lokph!!
   if(phlista(lokph)%tnooffr-phlista(lokph)%noofsubl.eq.0) then
!      write(*,*)'3G fixedcomposition: ',iph,lokph,&
!           phlista(lokph)%tnooffr,phlista(lokph)%noofsubl
      fixedcomposition=.true.
   else
      fixedcomposition=.false.
   endif
1000 continue
   return
 end function fixedcomposition

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!>      17. Miscellaneous
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function phvarlok
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

!\addtotable subroutine palmtree
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

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable subroutine sortinphtup
!\begin{verbatim}
 subroutine sortinphtup(n,m,xx)
! subroutine to sort the values in xx which are in phase and compset order
! in phase tuple order.  This is needed by the TQ interface
! The number of values belonging to the phase is m (for example composition)
! argument ceq added as new composition sets can be created ...
   integer n,m
!   double precision xx(n*m)
   double precision xx(*)
!   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
!
   integer iz,jz,kz,lz,lokph,aha,errall
   double precision, dimension(:), allocatable :: dum
! I assume the values are NP(*), maybe there are other cases ...
! Karl had overflow error in dum ... no problem to make it a little larger
! but then I cannot set xx=dum below ...
   allocate(dum(n*m+10),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 16: ',errall
      gx%bmperr=4370; goto 1000
   endif
!   write(*,*)'3F corrected sortinphtup',n,m
!   write(*,10)'3F in: ',(xx(iz),iz=1,n*m)
10 format(a,10(f7.4))
   kz=0
   do iz=1,noofph
      lokph=phases(iz)
      do jz=1,phlista(lokph)%noofcs
!         if(jz.gt.1) then
! in xx the values are sequentially for all composition sets for this phase
! But they should be stored in tuple order and compset 2 etc comes at the end
! the index to the tuple is in %phtups
! phlista(lokph)%linktocs(jz) is index of phase_varres record for compset
! firsteq%phase_varres(..)%phtupx is index of phase tuple for compset
! There can be m values (for example compositions) for each phase
! BUG FIXED: Sigli example gives hard error here
! index '0' of array 'firsteq' below lower boundary of 1
         aha=(firsteq%phase_varres(phlista(lokph)%linktocs(jz))%phtupx-1)*m
!         if(aha.ne.kz) then
!            write(*,*)'3F shifting from, to, values: ',kz,aha,m
!         endif
         do lz=1,m
            dum(aha+lz)=xx(kz+lz)
         enddo
         kz=kz+m
      enddo
   enddo
!   xx=dum
   do iz=1,n*m
      xx(iz)=dum(iz)
   enddo
   deallocate(dum)
!   write(*,10)'3F ut: ',(xx(iz),iz=1,n*m)
1000 continue
   return
 end subroutine sortinphtup

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function get_mpi_index
!\begin{verbatim}
  integer function get_mpi_index(mpi)
! Return the index of a model parameter identifier
    character mpi*(*)
!\end{verbatim} %+
! propid(jj)%symbol is *4
    character text*4
    integer jj
    text=mpi
    do jj=1,ndefprop
       if(propid(jj)%symbol.eq.text) exit
    enddo
!
    if(jj.gt.ndefprop) then
       write(*,*)'3A no such model parameter identifier: ',trim(mpi)
       gx%bmperr=4399; jj=-1
    endif
!    write(*,*)'3A get_mpi_index: ',text,ndefprop,jj
    get_mpi_index=jj
    return
  end function get_mpi_index

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function getmqindex
!\begin{verbatim}
  integer function getmqindex()
! This is necessary because mqindex is private, replaced by getmpiindex ...
!\end{verbatim}
    getmqindex=mqindex
    return
  end function getmqindex

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable logical function allowenter
!\begin{verbatim}
 logical function allowenter(mode)
! Check if certain commands are allowed
! mode=1 means entering an element or species
!     this routine is no longer used when entering species
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
! there must be at least one phase before entering a second equilibrium
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

!\addtotable logical function proper_symbol_name
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

!\addtotable subroutine compmassbug
!\begin{verbatim}
 subroutine compmassbug(ceq)
! debug subroutine
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer cp,sp,ep
! elements(1..n) is ordered alphabetcally
! complist(1..n) is initially also ordered alphabetcally but with errors ... 
! 
   do cp=1,noofel
      sp=ceq%complist(cp)%splink
      ep=splista(sp)%ellinks(1)
      write(*,100)cp,sp,ep,trim(ellista(ep)%name),trim(splista(sp)%symbol),&
           ellista(ep)%alphaindex,ellista(ep)%splink,&
           ceq%complist(cp)%mass,mass_of(cp,ceq),&
           ellista(ep)%mass,splista(sp)%mass
100   format(3i3,2x,a2,2x,a2,2i3,3x,4(1pe12.4))
   enddo
   write(*,*)
   return
 end subroutine compmassbug

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine list_free_lists
!\begin{verbatim}
 subroutine list_free_lists(lut)
! for debugging the free lists and routines using them
   implicit none
   integer lut
!\end{verbatim}
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

!\addtotable subroutine set_phase_amounts
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

!\addtotable subroutine set_default_constitution
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

!\addtotable subroutine todo_before
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

!\addtotable subroutine todo_after_found_equilibrium
!\begin{verbatim}
 subroutine todo_after_found_equilibrium(mode,addtuple,ceq)
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
! NOTE that automatically entered metallic-FCC and MC-carbides may shift
! composition sets. Such shifts can be avoided by manual entering composition
! sets with default constitutions, but that does not always work as comparing
! a stable constitution with several defaults is not trivial ...
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode,addtuple
!\end{verbatim}
   integer iph,ics,lokph,lokics,jcs,lokjcs,lastset,lokkcs,kzz,jtup,qq
   integer jstat2,fit,phs,haha1,haha2,disfravares,addph,icsno
   double precision val,xj1,xj2,extra(5)
   logical notok,noremove,globalok,once
   character jpre*4,jsuf*4
   real, dimension(:), allocatable :: tmmyfr
   double precision, dimension(:), allocatable :: yfr
! THIS ROUTINE MUST BE CLEANED UP
!
!   write(*,*)'3Y in todo_after',mode
!----------------------------------------------------------------
   addtuple=0
   if(btest(globaldata%status,GSNOAFTEREQ)) goto 1000
   nostart: if(mode.lt.0 .or. btest(globaldata%status,GSTGRID)) then
! if mode<0 the conditions did not allow gridmin before use it after
      if(btest(globaldata%status,GSNOGLOB)) goto 200
! Problems with this calculation so global_equil_check is disabled inside ...
      write(*,3)
3     format('Testing if any gridpoint is below the calculated equilibrium')
      if(btest(globaldata%status,GSNOTELCOMP)) then
         write(*,*)'3Y Cannot test global equilibrium when these components'
         goto 1000
      endif
      qq=1
! this generates a grid for test
      globalok=global_equil_check1(qq,addph,yfr,ceq)
!      write(*,*)'3Y Back from global_equil_check1',gx%bmperr,lokph
      if(globalok) then
! if TRUE equilibrium OK or it could not be tested
         if(gx%bmperr.ne.0) then
            write(*,*)'3Y Testing equilibrium with global minimizer failed'
            goto 1000
         endif
!         write(*,*)'3Y Grid minimizer test of equilibrium OK'
      else
! if FALSE the test showed this is not a global equilibrium, handle this!
         gx%bmperr=0
         lokph=phases(addph)
         write(*,*)'3Y Equilibrium wrong, gridtest found ',phlista(lokph)%name
         if(btest(globaldata%status,GSNORECALC)) goto 1000
!         write(*,*)'3Y add as stable: ',phlista(lokph)%name
! we should add addph to the stable set of phases and recalculate
! we have to check the state of all composition sets
         do ics=1,9
            lokics=phlista(lokph)%linktocs(ics)
!            write(*,*)'3Y Checking: ',addph,ics,lokph,lokics
            if(lokics.eq.0) then
! If we are not allowed to create composition sets quit
               if(btest(globaldata%status,GSNOACS)) goto 1000
! we have to create a new composition set and set that stable
! return with error code to initiate new calculation, icsno returned
!               write(*,*)'3Y Creating a new composition set',addph
               call enter_composition_set(addph,'    ','CHKD',icsno)
               if(gx%bmperr.ne.0) then
                  write(*,*)'Error creating composition set',gx%bmperr
                  goto 1000
               endif
               call get_phase_compset(addph,icsno,lokph,lokics)
               if(gx%bmperr.ne.0) goto 1000
!               write(*,*)'3Y ceated comp.set: ',ceq%phase_varres(lokics)%phtupx
               ceq%phase_varres(lokics)%status2=&
                    ibset(ceq%phase_varres(lokics)%status2,CSAUTO)
! we must set the constitution also!!
               call set_constitution(addph,icsno,yfr,extra,ceq)
               if(gx%bmperr.ne.0) then
                  write(*,*)'3Y error setting y of new comp.set'
                  goto 1000
               endif
! set some positive amount
               ceq%phase_varres(lokics)%amfu=1.0D-3
! I do not think the tuple has been created ... just set the phase index
               addtuple=ceq%phase_varres(lokics)%phtupx
               write(*,*)'3Y recalculate with: ',addtuple,icsno,lokics
               gx%bmperr=4358
               goto 1000
            elseif(ceq%phase_varres(lokics)%phstate.lt.PHENTUNST) then
! this set is dormant, ignore this stable ste and give no error
               write(*,*)'3Y Skip dormant phase no',&
                    ceq%phase_varres(lokics)%phtupx
! if there is a dormant composition set do not enter a new
               goto 200
            elseif(ceq%phase_varres(lokics)%phstate.eq.PHENTERED) then
! This composition set is entered but not stable, set it as stable and
! jump back with error code and calculate again
               ceq%phase_varres(lokics)%phstate=PHENTSTAB
               ceq%phase_varres(lokics)%amfu=1.0D-3
               addtuple=ceq%phase_varres(lokics)%phtupx
!               write(*,*)'3Y recalculate with gridpoint stable',addtuple,ics
! we must set the constitution also!!
               call set_constitution(addph,ics,yfr,extra,ceq)
               gx%bmperr=4358
               goto 1000
            else
! phase is stable or even fix, there is a miscibility gap, check if there are
! any more composition sets?
!               write(*,*)'3Y phase has a stable comp.set: ',addph,ics
            endif
         enddo
! no cleanup
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
! SEGNENTATION FAULT efter this write statement when reading unformatted file
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

!\addtotable subroutine checkdefcon
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

!\addtotable subroutine shiftcompsets
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
      endif manycs
   enddo phloop
! 
1000 continue
   return
 end subroutine shiftcompsets

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine switch_compsets2
!\begin{verbatim} %-
 subroutine switch_compsets2(lokph,ics1,ics2,ceq)
! copy constitution and results from ic2 to ic1 and vice versa
   integer lokph,ics1,ics2
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
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

!\addtotable subroutine select_composition_set
!\begin{verbatim}
 subroutine select_composition_set(iph,ics,yarr,ceq)
! PROBABLY NOT USED but should be implemenented
! if phase iph wants to become stable and there are several user defined
! composition sets with default composition limits this subroutine tries to
! select the one that fits these limits best.
! For example if an FCC phase that could be an austenite (low carbon content)
! or a cubic carbo-nitride (high carbon or nitrogen content, low vacancy)
! Less easy to handle ordered phases like B2 or L1_2 as ordering can be
! in any sublatittice ... but with option B and F that is possible
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
! if only one compset return this!
   ics=best
!
1000 continue
   return
 end subroutine select_composition_set

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\addtotable subroutine verify_phase_varres_array
!\begin{verbatim}
 subroutine verify_phase_varres_array(ieq,verbose)
! This subroutine checks that the phase varres array is consistent
! in equilibrium ieq.  For ieq=1 it also checks the free list
! UNFINISHED and not yet used BUT IMPORTANT
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

!\addtotable subroutine set_emergency_startpoint
!\begin{verbatim}
 subroutine set_emergency_startpoint(mode,phl,amfu,ceq)
! this is called if no previous equilibrium and if grid minimizer
! cannot be used.  Select for each element a phase with as much of that
! element as possible to set as stable. Set the remaining phases to a default
! composition.  It will never create any compositon sets
!
   implicit none
   integer mode,phl(*)
   double precision amfu(*)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer iph,lokph,lokcs,iel,errall
   integer, allocatable, dimension(:) :: selected
   double precision, allocatable, dimension(:,:) :: maxel
   double precision, allocatable, dimension(:) :: wmass
   double precision totmol,totmass,am
!
   write(*,*)'In emergency startpoint: ',mode,noofel,noofph
   allocate(selected(noofel),stat=errall)
   allocate(maxel(noofel,noofph),stat=errall)
   allocate(wmass(noofel),stat=errall)
   if(errall.ne.0) then
      write(*,*)'3Y allocation error 17: ',errall
      gx%bmperr=4370; goto 1000
   endif
!   phl=0
!   amfu=zero
   maxel=zero
   phloop1: do iph=1,noofph
      lokph=phases(iph)
      lokcs=phlista(iph)%linktocs(1)
      if(ceq%phase_varres(lokcs)%phstate.le.PHDORM) cycle phloop1
      if(phlista(iph)%tnooffr-phlista(iph)%noofsubl.eq.0) then
         call calc_phase_molmass(iph,1,maxel(1,iph),wmass,totmol,totmass,am,ceq)
         if(gx%bmperr.ne.0) goto 1000
      else
         write(*,*)'3Y TODO: Phases with variable composition not included yet'
      endif
! loop through all fractions to find limits            
   enddo phloop1
!   do iph=1,noofph
!      write(*,100)iph,(maxel(iel,iph),iel=1,noofel)
!   enddo
100 format('3Y maxel: ',i3,6(F8.5))
   selected=0
   wmass=zero
   phloop2: do iph=1,noofph
      lokph=phases(iph)
      lokcs=phlista(iph)%linktocs(1)
      if(ceq%phase_varres(lokcs)%phstate.le.PHDORM) cycle phloop2
      elloop1: do iel=1,noofel
         if(maxel(iel,iph).gt.wmass(iel)) then
            wmass(iel)=maxel(iel,iph)
            selected(iel)=iph
! we can only have one element selected per phase ...
            cycle phloop2
         endif
      enddo elloop1
   enddo phloop2
!   write(*,*)'3Y Emergency startpoint testing',mode
!   write(*,200)'3Y selected: ',(selected(iel),iel=1,noofel)
200 format(a,10i4)
! Now set default constitution of all non-selected and non-suspended phases
   phloop3: do iph=1,noofph
      lokph=phases(iph)
      lokcs=phlista(iph)%linktocs(1)
      if(ceq%phase_varres(lokcs)%phstate.le.PHDORM) cycle phloop3
      do iel=1,noofel
         if(iph.eq.selected(iel)) cycle phloop3
      enddo
!      write(*,*)'3Y TODO set default constitutions: ',iph
      call set_default_constitution(iph,1,ceq)
      if(gx%bmperr.ne.0) goto 1000
   enddo phloop3
   mode=noofel
   do iph=1,mode
      phl(iph)=selected(iph)
   enddo
!
1000 continue
   return
 end subroutine set_emergency_startpoint

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function ocv
!\begin{verbatim}
 logical function ocv()
! returns TRUE if GSVERBOSE bit is set
!\end{verbatim} %+
! typical use:  if(ocv()) write(*,*)....
   ocv=btest(globaldata%status,GSVERBOSE)
   return
 end function ocv

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function ceqsize
!\begin{verbatim}
 integer function ceqsize(ceq)
! calculates the size in words (4 bytes) of an equilibrium record
   implicit none
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer sum,vsum,ivs,vss
!   write(*,*)'In ceqsize 1'
!
!     integer status,multiuse,eqno,next
!     character eqname*24
!     double precision tpval(2),rtn
! svfunres: the values of state variable functions valid for this equilibrium
!     double precision, dimension(:), allocatable :: svfunres
   sum=18+2*size(ceq%svfunres)
   write(*,*)'total + svfunres: ',sum,size(ceq%svfunres)
! the experiments are used in assessments and stored like conditions 
! lastcondition: link to condition list
! lastexperiment: link to experiment list
!     TYPE(gtp_condition), pointer :: lastcondition,lastexperiment
! assuming a pointer is 4 bytes (2 words)
   sum=sum+4
! components and conversion matrix from components to elements
! complist: array with components
! compstoi: stoichiometric matrix of compoents relative to elements
! invcompstoi: inverted stoichiometric matrix
!     TYPE(gtp_components), dimension(:), allocatable :: complist
!     double precision, dimension(:,:), allocatable :: compstoi
!     double precision, dimension(:,:), allocatable :: invcompstoi
! a gtp_component record is about 20 words, invcompstoi same as compsoti
   if(allocated(ceq%complist)) sum=sum+20*size(ceq%complist)+&
        4*size(ceq%compstoi)
   write(*,*)'total + complist:',sum,20*size(ceq%complist),4*size(ceq%compstoi)
! one record for each phase+composition set that can be calculated
! phase_varres: here all calculated data for the phase is stored
!     TYPE(gtp_phase_varres), dimension(:), allocatable :: phase_varres
! each phase_varres record is different for each phase
   vsum=0
! highcs is highest used free phase_varres record
   do ivs=1,highcs
      vss=vssize(ceq%phase_varres(ivs))
      write(*,*)'Phase varres: ',ivs,vss
      vsum=vsum+vss
   enddo
   sum=sum+vsum
   write(*,*)'total + varres',sum,vsum
! index to the tpfun_parres array is the same as in the global array tpres 
! eq_tpres: here local calculated values of TP functions are stored
!     TYPE(tpfun_parres), dimension(:), pointer :: eq_tpres
! each tpfun_parres record is 8 double
   sum=sum+16*size(ceq%eq_tpres)
! current values of chemical potentials stored in component record but
! duplicated here for easy acces by application software
!     double precision, dimension(:), allocatable :: cmuval
   if(allocated(ceq%cmuval)) sum=sum+2*size(ceq%cmuval)
   write(*,*)'total + cmuval: ',sum,2*size(ceq%cmuval)
! xconc: convergence criteria for constituent fractions and other things
!     double precision xconv
! delta-G value for merging gridpoints in grid minimizer
! smaller value creates problem for test step3.BMM, MC and austenite merged
!     double precision :: gmindif=-5.0D-2
! maxiter: maximum number of iterations allowed
!     integer maxiter
   sum=sum+5
! this is to save a copy of the last calculated system matrix, needed
! to calculate dot derivatives, initiate to zero
!     integer :: sysmatdim=0,nfixmu=0,nfixph=0
!     integer, allocatable :: fixmu(:)
!     integer, allocatable :: fixph(:,:)
!     double precision, allocatable :: savesysmat(:,:)
   sum=sum+3+size(ceq%fixmu)+size(ceq%fixph)+size(ceq%savesysmat)
   write(*,*)'total + savesysmat:',sum,size(ceq%fixmu),size(ceq%fixph),&
        size(ceq%savesysmat)
   ceqsize=sum
1000 continue
   return
 end function ceqsize

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable integer function vssize
!\begin{verbatim}
 integer function vssize(varres)
! calculates the size in words (4 bytes) of a phase_varres record
   implicit none
   type(gtp_phase_varres) :: varres
!\end{verbatim}
   integer sum
!   write(*,*)'In vssize 1'
!     integer nextfree,phlink,status2,phstate
!     double precision, dimension(2) :: abnorm
!     character*4 prefix,suffix
   sum=10
! changed to allocatable
!     integer, dimension(:), allocatable :: constat
!     double precision, dimension(:), allocatable :: yfr
!     real, dimension(:), allocatable :: mmyfr
!     double precision, dimension(:), allocatable :: sites
   if(allocated(varres%constat)) sum=sum+size(varres%constat)
   write(*,*)'varressum+yfr: ',sum,size(varres%constat),3*size(varres%yfr)
   if(allocated(varres%yfr)) sum=sum+3*size(varres%yfr)
!   write(*,*)'In vssize 2',sum
! for ionic liquid derivatives of sites wrt fractions (it is the charge), 
! 2nd derivates only when one constituent is vacancy
! 1st sublattice P=\sum_j (-v_j)*y_j + Qy_Va
! 2nd sublattice Q=\sum_i v_i*y_i
!     double precision, dimension(:), allocatable :: dpqdy
!     double precision, dimension(:), allocatable :: d2pqdvay
   if(allocated(varres%dpqdy)) sum=sum+size(varres%dpqdy)
   if(allocated(varres%d2pqdvay)) sum=sum+size(varres%d2pqdvay)
   write(*,*)'varressum+ionliq',sum,size(varres%dpqdy),size(varres%d2pqdvay)
! for extra fraction sets, better to go via phase record index above
! this TYPE(gtp_fraction_set) variable is a bit messy.  Declaring it in this
! way means the record is stored inside this record.
!     type(gtp_fraction_set) :: disfra
! size of disfra record??
   sum=sum+10
   if(allocated(varres%disfra%dsites)) sum=sum+size(varres%disfra%dsites)
   if(allocated(varres%disfra%nooffr)) sum=sum+size(varres%disfra%nooffr)
   if(allocated(varres%disfra%splink)) sum=sum+size(varres%disfra%splink)
   if(allocated(varres%disfra%y2x)) sum=sum+size(varres%disfra%y2x)
   if(allocated(varres%disfra%dxidyj)) sum=sum+size(varres%disfra%dxidyj)
   write(*,*)'varresum incl disfra and pointer: ',sum,varres%disfra%varreslink
! It seems difficult to get the phdapointer in disfra record to work
! ---
! arrays for storing calculated results for each phase (composition set)
! amfu: is amount formula units of the composition set (calculated result)
! netcharge: is net charge of phase
! dgm: driving force (calculated result)
! amcom: not used
! damount: set to last change of phase amount in equilibrium calculations
! qqsave: values of qq calculated in set_constitution
!    double precision amount(2),dgm,amcom,damount,qqsave(3)
!    double precision amfu,netcharge,dgm,amcom,damount,qqsave(3)
!     double precision amfu,netcharge,dgm,amcom,damount
   sum=sum+10
! Other properties may be that: gval(*,2) is TC, (*,3) is BMAG, see listprop
! nprop: the number of different properties (set in allocate)
! ncc: total number of site fractions (redundant but used in some subroutines)
! BEWHARE: ncc seems to be wrong using TQ test program fenitq.F90 ???
! listprop(1): is number of calculated properties
! listprop(2:listprop(1)): identifies the property stored in gval(1,ipy) etc
!   2=TC, 3=BMAG. Properties defined in the gtp_propid record
!     integer nprop,ncc
!     integer, dimension(:), allocatable :: listprop
   if(allocated(varres%listprop)) sum=sum+2+size(varres%listprop)
   write(*,*)'varresum + listprop: ',sum,size(varres%listprop)
! gval etc are for all composition dependent properties, gval(*,1) for G
! gval(*,1): is G, G.T, G.P, G.T.T, G.T.P and G.P.P
! dgval(1,j,1): is first derivatives of G wrt fractions j
! dgval(2,j,1): is second derivatives of G wrt fractions j and T
! dgval(3,j,1): is second derivatives of G wrt fractions j and P
! d2gval(ixsym(i,j),1): is second derivatives of G wrt fractions i and j
!     double precision, dimension(:,:), allocatable :: gval
!     double precision, dimension(:,:,:), allocatable :: dgval
!     double precision, dimension(:,:), allocatable :: d2gval
   if(allocated(varres%gval)) sum=sum+2*size(varres%gval)
   if(allocated(varres%dgval)) sum=sum+2*size(varres%dgval)
   if(allocated(varres%d2gval)) sum=sum+2*size(varres%d2gval)
   write(*,*)'varresum + gvals: ',sum,2*size(varres%gval),&
        2*size(varres%dgval),2*size(varres%d2gval)
! added for strain/stress, current values of lattice parameters
!     double precision, dimension(3,3) :: curlat
! saved values from last equilibrium calculation
!     double precision, dimension(:,:), allocatable :: cinvy
!     double precision, dimension(:), allocatable :: cxmol
!     double precision, dimension(:,:), allocatable :: cdxmol
   if(allocated(varres%cinvy)) sum=sum+18+2*size(varres%cinvy)
   if(allocated(varres%cxmol)) sum=sum+18+2*size(varres%cxmol)
   if(allocated(varres%cdxmol)) sum=sum+18+2*size(varres%cdxmol)
   write(*,*)'varresum + saved: ',sum
!
1000 continue
   vssize=sum
   return
 end function vssize

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\addtotable logical function inveq
!\begin{verbatim}
  logical function inveq(phases,ceq)
! Only called for mapping tie-lines not in plane.  If tie-lines in plane
! then all nodes are invariants.
    integer phases
    type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
    integer nrel,ii,nostph,tpvar,degf,www
    type(gtp_condition), pointer :: pcond,lastcond
    type(gtp_state_variable), pointer :: stvr
! How to know if the ceq is invariant? Gibbs phase rule, Degrees of freedom
! f = n + z - w - p
! where n is number of components, z=2 if T and P variable, 
!                      z=1 if T or P variable, z=0 if both T and P fixed,
!                      w is number of other potential conditions (MU, AC)
!                      p is number of stable phases.
!    write(*,*)'3Y in inveq'
    nrel=noel()
! sum up nubler of stable phases and check if T and P are fixed
    nostph=0
!    ntups=nooftup()
!    do ii=1,noofphasetuples()
    do ii=1,nooftup()
       if(ceq%phase_varres(phasetuple(ii)%lokvares)%phstate.gt.0) &
            nostph=nostph+1
    enddo
! loop all conditions
    lastcond=>ceq%lastcondition
    pcond=>lastcond
    tpvar=2
    www=0
100 continue
       if(pcond%active.eq.0) then
! condtion is active
          stvr=>pcond%statvar(1)
! statevarid 1 is T and 2 is P
          if(stvr%statevarid.eq.1 .or. stvr%statevarid.eq.2) then
! Hm, ceq is not the equilibrium record for the node point ...
             tpvar=tpvar-1
          elseif(stvr%statevarid.lt.10) then
! potential/activity condition for a component
             www=www+1
          endif
       endif
       pcond=>pcond%next
       if(.not.associated(pcond,lastcond)) goto 100
!
! Hm again, ignore tpvar?
!    degf=nrel+tpvar-www-nostph
    degf=nrel-www-nostph
!    write(*,'(a,8i4)')'3Y in inveq 2',nrel,tpvar,www,nostph,degf
    if(degf.lt.0) then
! We have an invariant equilibrium, return the number of stable phases
       phases=nostph
       inveq=.true.
!       write(*,200)'3Y An invariant equilibrium!',nrel,tpvar,nostph,phases
!200    format(a,5i7)
    else
!       write(*,210)degef,nrel,tpvar,phases
!210     format('3Y not invariant eq, elements, stable phases: ',4i4)
!       if not invariant isoplet node there are 3 exits (2 lines crossing)
       phases=nostph
       inveq=.false.
    endif
1000 continue
    return
  end function inveq

!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

