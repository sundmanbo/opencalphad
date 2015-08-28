!
! gtp3Y included in gtp3.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>      16. Grid minimizer
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine global_gridmin(what,tp,xknown,nvsph,iphl,icsl,aphl,&
      nyphl,yphl,cmu,ceq)
!
! finds a set of phases that is a global start point for an equilibrium 
! calculation at T and P values in tp, total amount of atoms in totan
! and known mole fraction in xknown
! It is intentional that this routine is independent of current conditions
! It returns: nvsph stable phases, list of phases in iphl, amounts in aphl, 
! constitution in yphl (compact after each other, nyphl(i) is number of
! fractions in phase i), cmu are element chemical potentials of solution
! WHAT determine what to do with the results, 0=just return solution,
! 1=enter stable set and constitution of all phases in gtp datastructure
! and create composition sets if necessary (and allowed)
! what=-1 will check if any gridpoint below current calculated equilibrium
   implicit none
! nyphl(j) is the start position of the constitiuent fractions of phase j in
! yphl that contains all the constitutions of the phases in the gridpoints
   integer, dimension(*) :: iphl,nyphl,icsl
   integer what,nvsph
   TYPE(gtp_equilibrium_data), pointer :: ceq
   double precision, dimension(2) :: tp
! cmu(1..nrel) is the chemical potentials of the solution
   double precision, dimension(*) :: xknown,aphl,yphl,cmu
!\end{verbatim}
   integer, parameter :: maxgrid=400000,maxy=2000,maxph=500
   integer :: starttid,endoftime
   real finish2
   double precision amount,sum
   integer i,ibias,ics,ics2,icsno,icsx,ie,iph,iv,j1,j2,jip,jp,kkz,kp,kph,jbias
   integer lokcs,lokph,mode,ng,nocsets,noofgridpoints,nr,nrel,nrph,ny,nyz
! kphl(iph) is first gridpoint in phase iph
! ngrid(iph) is the last gridpoint for phase iph (some phases may be suspended)
! xarr(nrel,i) is the composition of gridpoint i
! garr(i) is the Gibbs energy of gridpoint i
! jgrid(j) is a gridpoint in the solution
! phfrac(j) is the amount of the phase of that gridpoint
   integer, dimension(0:maxph) :: ngrid,kphl
   integer, dimension(maxel) :: jgrid
   real garr(maxgrid),starting,finished
   real, dimension (:,:), allocatable :: xarr
   real, dimension (maxel,maxel) :: xsol
   double precision, dimension(maxel) :: phfrac,phsave
   double precision qq(5),savetp(2)
   integer, dimension(maxph) :: iphx
   character name1*24
! debug
   logical trace
! sort phases depending on number of gridpoints
   integer, dimension(:), allocatable :: gridpoints,phord
! pph is set to number of phases participating, some may be suspended
   integer pph,zph,nystph,order(maxel)
!
   if(btest(globaldata%status,GSNOGLOB)) then
      write(*,*)'Grid minimization not allowed'
      gx%bmperr=4173; goto 1000
   endif
   call cpu_time(starting)
   ngrid=0
! Trace turn on output of grid on a file
!   trace=.true.
   trace=.FALSE.
   if(trace) write(*,*)'Trace set TRUE'
   savetp=ceq%tpval
   ceq%tpval=tp
!    ceq%tpval(2)=tp(2)
   nrph=noph()
!    write(*,*)'ggp A: ',tp(1),ceq%tpval(1)
   if(nrph.gt.maxph) then
! too many phases
      write(*,*)'Too many phases for gridmin'
      gx%bmperr=6663; goto 1000
   endif
   nrel=noel()
   sum=zero
   do i=1,nrel
      if(xknown(i).le.zero .or. xknown(i).ge.one) then
!         write(*,*)'Illegal composition for gridmin'
         gx%bmperr=4174; goto 1000
      endif
      sum=sum+xknown(i)
   enddo
   if(ocv()) write(*,12)'gridmin: ',sum,(xknown(i),i=1,nrel)
12 format(a,1pe12.4,10(f8.4))
   if(abs(sum-one).gt.1.0D-8) then
      write(*,*)'Sum of fractions larger than unity calling global_gridmin'
      gx%bmperr=4174; goto 1000
   endif
   kp=1
   pph=0
   kphl(0)=0
   allocate(gridpoints(nrph))
   allocate(phord(nrph))
   ggloop: do iph=1,nrph
!      if(.not.phase_status(iph,1,PHHID,ceq)) then
! skip phases that are hidden or suspended .... old 5; new -3
!      if(test_phase_status(iph,1,amount,ceq).lt.5) then
!      if(test_phase_status(iph,1,amount,ceq).gt.PHHIDDEN) then
!      if(test_phase_status(iph,1,amount,ceq).gt.PHSUS) then
! skip dormant phases
! include phases with first composition set entered (only once!)
      ent1: if(test_phase_status(iph,1,amount,ceq).gt.PHDORM) then
         do ics=1,noofcs(iph)
! new: -3 suspended, -2 dormant, -1,0,1 entered, 2 fixed
! ignore phases whith no composition set entered
! If a phase+compset FIX one should never be here as conditions wrong
            if(test_phase_status(iph,ics,amount,ceq).lt.PHFIXED) goto 60
         enddo
         cycle ggloop
! this call to find out how many gridpoints will be generated for each phase
60       continue
         pph=pph+1
         kphl(pph)=kp
         iphx(pph)=iph
         call generate_grid(-1,iph,ng,nrel,xarr,garr,ny,yphl,ceq)
         if(gx%bmperr.ne.0) goto 1000
         kp=kp+ng
         ngrid(pph)=kp-1
         gridpoints(pph)=ng
!         write(*,61)'3Y gpno: ',pph,iph,ngrid(pph),gridpoints(pph)
61       format(a,10i7)
!         if(trace) then
            call get_phase_name(iph,1,name1)
!            write(*,21)iph,name1(1:12),kphl(pph),ng
!         endif
!         write(*,22)iph,kphl(pph),ny,ng,pph
      endif ent1
   enddo ggloop
! we have a grid for pph phases, note that pph is not a phase index!!!
! the phase index for phase 1..pph is in iphx(1..pph)
21 format('Gridpoints for phase ',i3,': ',a,', starts ',i7,', with ',2i7)
22 format('Gridpoints for phase ',i3,' starts at ',i5,', with ',2i5,i8)
   if(kp-1.gt.maxgrid) then
      write(*,*)'Too many gridpoints'
      gx%bmperr=4175; goto 1000
   endif
! we may have no gridpoints!!!
   if(kp.eq.1) then
      write(*,*)'No phases, no gridpoints'
      gx%bmperr=4176; goto 1000
   endif
!   write(*,*)'phases and gridpoints: ',pph,kp,ngrid(pph),nrel
! total number of gridpoints is kp-1 ... but sometimes kp is wrong, why??
!   allocate(xarr(nrel,kp-1))
   allocate(xarr(nrel,kp-1+10))
   if(ocv()) write(*,*)'Gridpoints and elements: ',kp-1,nrel
! we should sort iphx to have phases with many gripoints first
! kphl must be shifted at the same time
!$   call sortin(gridpoints,pph,phord)
!-$   write(*,32)(kp,phord(kp),gridpoints(kp),kp=1,pph)
32 format(4('3Y ord: ',i2,i3,i6,' : '))
! generate grid
! we must know before this loop how many gridpoints that each phase will
! need.  That is a function of the number of gridpoints.
   kp=1
   call system_clock(count=starttid)
!   write(*,*)'Start calculate gridvalues'
! OpenMP parallellization START
! the error code gx%bmperr should also be threadprivate
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! for parallelizing:
! YOU MUST UNCOMMENT USE OMP_LIB IN GTP3.F90
! YOU MUST USE THE SWICH -fopenmp FOR COMPILATION AND WHEN LINKING
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!-$omp parallel do private(ng,iv,iph),schedule(dynamic)
!--$omp parallel do private(ng,iv),schedule(dynamic)
!--$omp parallel do 
!    phloop: do iph=1,nrph
   phloop: do zph=1,pph
! for phase iphx(zph) the gridpoints will be stored from position kphl(zph)
! mole fracts in xarr, g in garr
! yphl is not used when mode=0, ng should be set to number of remaining points
! ngrid(iph) is number of gridpoints in phase iph
      ng=maxgrid
! values in kphl set in previous call to generate_grid(-1,.....)
      iv=kphl(zph)
! when not parallel set iph=zph
      iph=zph
! for parallel take the phases from phord(pph), phord(pph-1) .... 2, 1
!-$      iph=phord(pph+1-zph)
!-$      iv=kphl(iph)
!--$      write(*,42)'Thread: ',omp_get_thread_num(),zph,iph,&
!--$           iphx(iph),iv,gridpoints(pph+1-zph)
42    format(a,10i7)
! this call will calculate all gridpoints, that may take time ...
!      call generate_grid(0,iphx(iph),ng,nrel,xarr(1,iv),garr(iv),ny,yphl,ceq)
      call generate_grid(0,iphx(iph),ng,nrel,xarr(1,iv),garr(iv),ny,yphl,ceq)
      if(gx%bmperr.ne.0) then
         write(*,*)'grid error ',jip,zph,gx%bmperr
! this jump illegal when openmp
!         goto 1000
         gx%bmperr=0
      endif
! list xarr for all gridpoints
!      do kp=1,ng
!         write(*,73)iphx(zph),kp,(xarr(ie,kp),ie=1,nrel)
!73       format('gp: ',i3,i5,10(f6.3))
!      enddo
!      write(*,*)'look!!!'
!      read(*,74)ch1
!74    format(a)
   enddo phloop
! set how many points in 
!-$omp end parallel do
! OpenMP parallellization END
   call system_clock(count=endoftime)
!    write(*,106)endoftime-starttid
106 format('Clockcycles: ',i12)
107 format(a,7i6)
!    kp=ngrid(nrph)
   kp=ngrid(pph)
!    if(trace) write(*,108)kp
108 format('The total number of gridpoints are ',i5)
   call cpu_time(finished)
   noofgridpoints=ngrid(pph)
! If WHAT is -1 then just compare all gridpoints with plane defined by
! the chemical potentials cmu to see if any is below.
! If so insert the gridpoint furtherst below the plane and set WHAT 10*iph+ics
!   write(*,*)'global_gridmin what: ',what
   if(what.eq.-1) then
      call gridmin_check(nystph,kp,nrel,xarr,garr,xknown,ngrid,pph,&
           cmu,yphl,iphx,ceq)
      goto 1000
   endif
!-----------------------------------------------
!    write(*,109)ngrid(pph),finished-starting,endoftime-starttid
109 format('Calculated ',i6,' gridpoints in ',1pe12.4,' seconds, ',&
         i7,' clockcycles')
! find the minimum of nrel gridpoints among the kp-1 gridpoint
! for current overall composition, xknown
!    write(*,*)'globm 4: ',kp,garr(kp),xarr(1,kp)
!   phfrac=zero
   if(ocv()) write(*,*)'Finding the gridpoints for the minimum: ',kp-1
   call find_gridmin(kp,nrel,xarr,garr,xknown,jgrid,phfrac,cmu,trace)
   if(gx%bmperr.ne.0) goto 1000
! The solution with nrel gridpoints are in jgrid, the amount of each in phfrac
! We later want the phases in ascending order and as the gridpoints are
! in ascending order of the phases we sort the gridpoints (and amounts)
! There must be one gridpoint per component (element)
!   write(*,62)(jgrid(jp),jp=1,nrel)
   call sortin(jgrid,nrel,order)
   do nyz=1,nrel
      phsave(nyz)=phfrac(order(nyz))
   enddo
   phfrac=phsave
! check
!   write(*,62)(jgrid(jp),jp=1,nrel)
62 format('3Y Gridp: ',10i4)
!   xov=zero
!   sum=zero
!   do jp=1,nrel
!      write(*,63)'xs: ',phfrac(jp),(xarr(nyz,jgrid(jp)),nyz=1,nrel)
!      do nyz=1,nrel
!         xov(nyz)=xov(nyz)+phfrac(jp)*xarr(nyz,jgrid(jp))
!      enddo
!      sum=sum+phfrac(jp)
!   enddo
!   write(*,63)'ss: ',sum,(xov(nyz),nyz=1,nrel)
!63 format(a,1e12.4,10f7.4)
! get the phase and constitution for each
   nyz=1
!   if(trace) write(*,*)'Extracting constititution'
   if(trace) then
      write(31,745)
745   format(/'Solution: ')
   endif
   solloop: do jp=1,nrel
! jgrid(jp) is a grid point in the solution, find which phase it is
      mode=jgrid(jp)
      ibias=0
      do zph=1,pph
!          write(*,*)'mode and ibias 1: ',mode,ibias
! ngrid(zph) is the first gridpoints of phase zph
         if(mode.le.ngrid(zph)) then
            mode=mode-ibias
            goto 315
         else
            ibias=ngrid(zph)
         endif
      enddo
      write(*,*)'gridpoint outside range ',jgrid(jp),ngrid(pph)
      gx%bmperr=4147; goto 1000
315   continue
      jbias=ibias
!      write(*,*)'gridpoint in solution: ',mode,ibias
! this call is to obtain the constitution of a phase in the solution
! mode gives in grid point index in phase iphx(zph), ibias irrelevant (?)
! NOTE ibias is changed by subroutine
      call generate_grid(mode,iphx(zph),ibias,nrel,xarr,garr,ny,yphl(nyz),ceq)
      if(gx%bmperr.ne.0) goto 1000
!       write(*,317)'gg7B: ',ny,nyz,(yphl(i),i=nyz,nyz+ny-1)
!317    format(a,2i3,6(1pe11.3))
      iphl(jp)=iphx(zph)
      aphl(jp)=phfrac(jp)
      nyphl(jp)=ny
      nyz=nyz+ny
! finally copy the mole fractions to xsol, needed for possible merging
      do ie=1,nrel
         xsol(ie,jp)=xarr(ie,mode+jbias)
      enddo
      if(trace) then
         write(31,750)jp,jgrid(jp),iphl(jp),aphl(jp),(xsol(ie,jp),ie=1,nrel)
         write(31,760)(yphl(i),i=nyz-ny,nyz-1)
750      format('Point: ',i2,', gridpoint: ',i5,' phase ',i3,&
              ' amount: ',1pe12.4,', Mole fractions:'/9(0pF8.5))
760      format('Constitution:'/9(0pF8.5))
      endif
   enddo solloop
   if(trace) then
      write(*,*)'Closing grid file'
      close(31)
   endif
! there must be as many phases in the solution as there are elements
   nvsph=nrel
   nr=nvsph
   if(.not.btest(globaldata%status,GSNOMERGE)) then
      call merge_gridpoints(nr,iphl,aphl,nyphl,yphl,trace,nrel,xsol,cmu,ceq)
      if(gx%bmperr.ne.0) goto 1000
   endif
! number of gridpoints, nr, may have changed
!   write(*,*)'After merge_gripoints: ',nr,nvsph
   nvsph=nr
! if what=-1 or0 do nothing more, just exit
!   if(what.eq.0) goto 1000
   if(what.le.0) goto 1000
!------------------------------------------------------------
! prepare for storing result: zero all phase amounts and driving forces
   do iph=1,nrph
      lokph=phases(iph)
!      lokcs=phlista(lokph)%cslink
      do ics=1,phlista(lokph)%noofcs
         lokcs=phlista(lokph)%linktocs(ics)
!         ceq%phase_varres(lokcs)%amount=zero
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
!      write(*,*)'grid chemical potential: ',ie,cmu(ie)*ceq%rtn
! do not care about reference state for chempot(2)
      ceq%complist(ie)%chempot(1)=cmu(ie)*ceq%rtn
      ceq%complist(ie)%chempot(2)=cmu(ie)*ceq%rtn
   enddo
! set driving force 0 for stable phases
   do i=1,nvsph
      call set_driving_force(iphl(i),1,zero,ceq)
      if(gx%bmperr.ne.0) goto 1000
   enddo
! store the most favourable constitution of the metastable phase
   call set_metastable_constitutions(kp,nrel,kphl,ngrid,xarr,garr,&
        nvsph,iphl,cmu,ceq)
   if(gx%bmperr.ne.0) goto 1000
! maybe more composition sets needed
   do i=1,nvsph
      icsl(i)=0
   enddo
   nocsets=0
!   write(*,*)'before loop1: ',nvsph,ceq%eqname
   loop1: do j1=1,nvsph
      if(icsl(j1).eq.0) then
! if non-zero a composition set has already been assigned 
         icsl(j1)=1
         icsx=1
         loop2: do j2=j1+1,nvsph
            if(iphl(j1).eq.iphl(j2)) then
! one more composition set needed, does it exist?
               icsx=icsx+1
               ics2=icsx
               call get_phase_compset(iphl(j1),ics2,lokph,lokcs)
               if(gx%bmperr.ne.0) then
! there is no such composition set, is automatic creation allowed?
! NOTE: there is a EQNOACS bit also???
                  if(btest(globaldata%status,GSNOACS)) then
!                     write(*,*)'Not allowed to create composition sets'
                     gx%bmperr=4177; goto 1000
                  endif
                  gx%bmperr=0
! >>>>>>>>>>>>>>>>>>><
! BEWARE >>> not only must this be done in all threads at the same time
! one must also avoid that it is done when some thread is working on a set
! of phase+composition sets trandformed to EQCALC arrays.  If so the
! indices to lokcs etc will be incorrect ... ???
! I think OMP has "secure" points where the treads can be stopped to wait
! <<<<<<<<<<<<<<<<<<<<<<
                  kph=iphl(j1)
!                  write(*,*)'3Y new composition set for phase: ',j2,kph
! It must be done in all equilibrium records, no equilibrium record needed!!!
! one must be careful with the status word when creating comp.sets
                  call enter_composition_set(kph,'    ','AUTO',icsno)
                  if(gx%bmperr.ne.0) goto 1000
                  call get_phase_compset(kph,icsno,lokph,lokcs)
                  if(gx%bmperr.ne.0) goto 1000
                  ceq%phase_varres(lokcs)%status2=&
                       ibset(ceq%phase_varres(lokcs)%status2,CSAUTO)
                  nocsets=nocsets+1
!                  if(btest(ceq%phase_varres(lokcs)%status2,CSDEFCON)) then
!                     write(*,*)'3Y defcon set',kph,icsno
!                  else
!                     write(*,*)'3Y defcon not set',kph,icsno
!                  endif
!                  write(*,303)'3Y Created cs:',kph,icsno,lokcs,&
!                       ceq%phase_varres(lokcs)%amfu,&
!                       ceq%phase_varres(lokcs)%abnorm
303               format(a,3i3,6(1pe12.4))
                  icsl(j2)=icsno
               else
! here we should check which composition set that should have which 
! constitution for example one fcc is metallic and another is cubic carbide
                  call get_phase_name(iphl(j1),ics2,name1)
                  icsl(j2)=ics2
!                  write(*,1711)name1,ics2
1711              format('Using composition set for ',a,i3)
! check if the composition set is fix (2), dormant (2) or suspended (3)
                  kkz=test_phase_status(iphl(j1),ics2,amount,ceq)
! old kkz=2 means fix
!                  if(kkz.eq.2) then
                  if(kkz.eq.PHFIXED) then
                     write(*,*)'Global minimization with fix phase not allowed'
                     gx%bmperr=7777; goto 1000
                  elseif(kkz.lt.PHENTUNST) then
                     write(*,*)' *** Warning, changing status for phase ',name1
                  endif
! this means status entered PHSTATE
!                  ceq%phase_varres(lokcs)%status2=&
!                       ibclr(ceq%phase_varres(lokcs)%status2,CSSUS)
!                  ceq%phase_varres(lokcs)%status2=&
!                       ibclr(ceq%phase_varres(lokcs)%status2,CSFIXDORM)
                  ceq%phase_varres(lokcs)%phstate=0
               endif
            endif
         enddo loop2
      endif
   enddo loop1
!   write(*,*)'after loop1: ',phlista(1)%noofcs
   if(nocsets.gt.0) write(*,*)'Composition set(s) created: ',nocsets
! Above one should consider if some user created compsets are dedicated to
! certain cases (MC carbides or L1_2 ordered).  These should have
! a default constitution and CSDEFCON set)
! finally store stable phase amounts and constitutions into ceq%phase_varres
   j1=1
   ceqstore: do iph=1,nvsph
!      write(*,*)'ggm: ',iph,iphl(iph),icsl(iph),j1
      call set_constitution(iphl(iph),icsl(iph),yphl(j1),qq,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,1788)'gg: ',iph,iphl(iph),icsl(iph),aphl(iph),&
!           (yphl(j1+ie),ie=0,3),qq(1)
1788  format(a,3i3,f8.4,2x,4f8.4,1pe10.2)
! aphl(iph) is amount of phase per mole component
      call get_phase_compset(iphl(iph),icsl(iph),lokph,lokcs)
! Here aphl is divided with the number of mole of atoms in the phase
!      if(ceq%phase_varres(lokcs)%abnorm(1).ge.one) then
!         aphl(iph)=aphl(iph)/ceq%phase_varres(lokcs)%abnorm(1)
!      endif
!      write(*,1812)iph,lokcs,aphl(iph),ceq%phase_varres(lokcs)%abnorm(1)
1812  format('aphl: ',2i3,6(1pe12.4))
      aphl(iph)=aphl(iph)/ceq%phase_varres(lokcs)%abnorm(1)
!      write(*,1789)'aphl: ',iph,lokcs,aphl(iph),&
!           ceq%phase_varres(lokcs)%abnorm(1)
1789  format(a,2i3,2(1pe12.4))
      ceq%phase_varres(lokcs)%amfu=aphl(iph)
      ceq%phase_varres(lokcs)%phstate=PHENTSTAB
!      write(*,*)'3Y gridmin stable: ',lokcs,ceq%phase_varres(lokcs)%phtupx
      j1=j1+nyphl(iph)
   enddo ceqstore
1000 continue
!   write(*,*)'at 1000: ',phlista(1)%noofcs
! restore tpval in ceq
   ceq%tpval=savetp
   call cpu_time(finish2)
   if(allocated(xarr)) deallocate(xarr)
   if(gx%bmperr.ne.0) then
!      globaldata%status=ibset(globaldata%status,GSEQFAIL)
      ceq%status=ibset(ceq%status,EQFAIL)
   elseif(what.eq.-1) then
      if(nystph.gt.0) what=nystph
   else
      write(*,1010)noofgridpoints,finish2-starting,endoftime-starttid
1010  format(' Grid minimization: ',i7,' gridpoints ',1pe12.4,' s and ',&
           i7,' clockcycles')
! set the global bit that this is not a full equilibrium
      ceq%status=ibset(ceq%status,EQNOEQCAL)
   endif
   if(ocv()) write(*,*)'leaving global_gridmin'
   return
 end subroutine global_gridmin

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!-\begin{verbatim}
 subroutine new_gridpoint_calc(only,iph,nsl,nend,endm,jend,ifra,ny,yfra,&
      xmol,gval,ceq)
! This subroutine sets the fractions according to three indicators
! and calculates the Gibbs energy of the gridpoint
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer iph,nsl,only,nend,ny,jend(3),ifra(3)
   double precision yfra(*),gval,xmol(*)
   integer, dimension(nsl,nend) :: endm
!-\end{verbatim}
   integer lokres,ls
! preset fractions
   double precision qq(5)
   double precision, parameter :: yzero=1.0D-12
! preset weights of endmembers
   double precision, dimension(4), parameter:: ybas=&
        [1.0D0,0.89D0,0.74D0,0.61D0]
   double precision, dimension(4), parameter :: ybin=&
        [0.11D0,0.26D0,0.39D0,0.15D0]
   double precision, dimension(2), parameter :: yter=&
        [0.11D0,0.13D0]
! When setting fractions one must have the sum of fractions in each sublattice
! equal to unity.  This is done by weighting endmembers
! endm(ll,ie) is the index to constituent ie in sublattice ll
! With 3 sublattices with (2,2,4) constituents endm is
! endm(1,1)=1, endm(1,2)=2, 
! endm(2,1)=3, endm(2,2)=4,
! endm(3,1)=5, endm(3,2)=6, endm(3,3)=7, endm(3,4)=8
! jend(*) select endmember to mix, with two constitutions jend(3)=0
!
   do ls=1,ny
      yfra(ls)=zero
   enddo
!   write(*,10)iph,jend,ifra
10 format('gix: ',i2,13x,3i3,2x,3i3)
   do ls=1,nsl
      yfra(endm(ls,jend(1)))=ybas(ifra(1))
      if(jend(2).gt.0) yfra(endm(ls,jend(2)))=&
           yfra(endm(ls,jend(2)))+ybin(ifra(2))
      if(jend(3).gt.0) yfra(endm(ls,jend(3)))=&
           yfra(endm(ls,jend(3)))+yter(ifra(3))
   enddo
   if(only.gt.0) goto 1000
!
! calculate G and composition and save
   write(*,11)iph,(yfra(ls),ls=1,ny)
11 format('gyp: ',i2,25x,10(f6.3))
   call set_constitution(iph,1,yfra,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
   call calcg(iph,1,0,lokres,ceq)
   if(gx%bmperr.ne.0) goto 1000
!
   if(qq(1).ge.1.0D-1) then
! number of real atoms less than 10%, a gridpoint with just vacancies ....
      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1))
   else
      gval=1.0E5
   endif
   call calc_phase_mol(iph,xmol,ceq)
   if(gx%bmperr.ne.0) goto 1000
1000 continue
   return
 end subroutine new_gridpoint_calc

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine generate_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,ceq)
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
   double precision yarr(*)
!\end{verbatim}
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
!    0.99*Y_E + 0.01*Y_all                            N of these
!    0.89*Y_E + 0.10*Y_F,F=/=E + 0.01*Y_all           N*(N-1)
!    0.74*Y_E + 0.25*Y_F,F=/=E + 0.01*Y_all           N*(N-1)+N*(N-1)*(N-2)
!             + 0.15*Y_F + 0.1*Y_G,G=/=(E,F) + 0.01*Y_all  (3 or more endmemb)
!    0.61*Y_E + 0.38*Y_F,F=/=E + 0.01*Y_all
!             + 0.25*Y_F + 0.13*Y_G,G=/=(E,F) + 0.01*Y_all (3 or more endmemb)
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
   double precision, dimension(4), parameter:: ybas=&
        [1.0D0,0.89D0,0.74D0,0.61D0]
   double precision, dimension(4), parameter :: ybin=&
        [0.11D0,0.26D0,0.39D0,0.15D0]
   double precision, dimension(3), parameter :: yter=[0.0D0,0.11D0,0.13D0]
! for output of gridpoints
   integer jbas,sumngg,loksp
   logical trace,isendmem
   save sumngg
!
!   write(*,*)'entering generate_grid: ',mode,iph,ngg
   if(test_phase_status_bit(iph,PHEXCB)) then
! This phase has charged endmembers, generate neutral gridpoints
      call generate_charged_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,ceq)
      goto 1000
   elseif(test_phase_status_bit(iph,PHFORD)) then
! this phase has 4 sublattice fcc/hcp tetrahedral ordering,
! this reduces the number of gridpoints
      call generate_fccord_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,ceq)
!      goto 1000
   endif
   if(mode.eq.0) then
!      write(*,*)'Generating grid for phase: ',iph
! trace TRUE means generate outpt for each gridpoint
!      trace=.TRUE.
      trace=.FALSE.
      if(iph.eq.1 .and. trace) then
         open(31,file='gridgen.dat ',access='sequential')
         sumngg=0
         write(31,43)
43       format('The constituent fractions, y, enclosed within parentheses',&
              'for each sublattice'/'Mole fractions after x:, Gibbs energies',&
              ' after G:'/)
      endif
      if(trace) then
         call get_phase_record(iph,nend)
         write(31,44)iph,phlista(nend)%name
44       format('Endmembers (EM) and gridpoints (GP) for phase: ',i3,1x,a)
      endif
   else
      trace=.FALSE.
   endif
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
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
!   write(*,1010)'Saved   ',iph,(ydum(i),i=1,ny)
   negmode: if(mode.lt.0) then
!---------------------------------------------------------
! just determine the number of gridpoints for this phase for global minimimum
! ideal gases should just have the endmembers ....
! Hm, gases with ions??
      ngdim=ngg
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
         ngg=nend*(1+(nend-1)*(3+nend-2)) ! (ternary combination skipped)
         ngg=nend*(1+3*(nend-1))
      else
! 3-10: three binary and two ternary combinarions (all)
         ngg=nend*(1+(nend-1)*(3+2*(nend-2))) ! (ternary combinations skipped)
         ngg=nend*(1+3*(nend-1))
      endif
!      write(*,*)'endmembers and gridpoints: ',nend,ngg
!      read(*,11)ch1
11    format(a)
!      ngg=ngg+iliqneut
      ngg=ngg
      if(ocv()) write(*,*)'Generate grid: ',nend,ngg
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
!   write(*,*)'ggy: ',mode,iph,nsl,nend,inkl(nsl)
!   gx%bmperr=7777; goto 1000
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
!         write(31,125)i,(endm(ls,i),ls=1,nsl)
!125      format('endmem: ',i4,2x,10i3)
!      enddo
!   endif
150 continue
!---------------------------------------
! now generate all combinations of endmembers
!   write(*,*)'endmembers and gridpoints: ',nend,ngg
!   read(*,11)ch1
   ngg=0
   lokph=phases(iph)
   endmem: do iend=1,nend
      yfra=yzero
      do ls=1,nsl
         yfra(endm(ls,iend))=ybas(1)
      enddo
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
         if(ngg.eq.mode) goto 500
      else
! calculate G and composition and save
!         write(*,201)ibas,ngg,(yfra(is),is=1,inkl(nsl))
201      format('ggz: ',i2,i5,10(F6.3))
         if(ocv()) write(*,*)'Calculating gridpoint: ',ngg
         call calc_gridpoint(iph,yfra,nrel,xarr(1,ngg),garr(ngg),ceq)
         if(gx%bmperr.ne.0) goto 1000
!         if(ngg.eq.15) then
!            write(*,520)'cgx: ',(xarr(is,ngg),is=1,nrel)
!         endif
         if(trace) then
            if(isendmem) then
               write(31,153,advance='no')sumngg+ngg
153            format('EM:',i4,' y: ')
            else
               write(31,154,advance='no')sumngg+ngg
154            format('GP:',i4,' y: ')
            endif
            jbas=0
            do ls=1,nsl
               write(31,155,advance='no')(yfra(jbas+is),is=1,nkl(ls)-1)
155            format('(',10(F4.2,','))
               write(31,156,advance='no')yfra(jbas+nkl(ls))
156            format(F4.2,')')
               jbas=jbas+nkl(ls)
            enddo
            write(31,157)(xarr(is,ngg),is=1,nrel),garr(ngg)
157         format(' x:',3f8.5,' G:',1pe12.4)
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
!            write(*,*)'Ternary combinations for 10<nend<15'
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
!            write(*,*)'Ternary combinations for 2<nend<10'
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
            yfra(endm(ls,iend))=ybas(ibas)
            yfra(endm(ls,jend))=yfra(endm(ls,jend))+ybin(ibin)
            if(kend.gt.0) then
               yfra(endm(ls,kend))=yfra(endm(ls,kend))+yter(iter)
            endif
         enddo
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
!    write(*,505)'ggg: ',mode,iph,nsl,inkl(nsl),ny
505 format(a,i7,i4,2x,i2,i5,i10)
!   write(*,510)'ggy: ',mode,ngdim,ny,(yfra(i),i=1,ny)
510 format(a,2i5,i3,10(F6.3))
   do i=1,ny
      yarr(i)=yfra(i)
   enddo
!   write(*,520)'ggx: ',(xarr(is,ngg+ngdim),is=1,nrel)
520 format(a,10(f8.5))
1000 continue
   if(allocated(endm)) then
      deallocate(endm)
      deallocate(yfra)
   endif
1001 continue
! restore original constituent fractions
!   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   errsave=gx%bmperr
   gx%bmperr=0
!   write(*,1010)'Restore ',iph,(ydum(i),i=1,ny)
1010 format(a,'const for ',i3,10(f6.3))
   call set_constitution(iph,1,ydum,qq,ceq)
   if(gx%bmperr.ne.0) then
      write(*,*)'Error restoring constitution for phase: ',iph,gx%bmperr
   endif
   gx%bmperr=errsave
   if(gx%bmperr.ne.0) write(*,*)'gengrid error: ',gx%bmperr
   return
 end subroutine generate_grid

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
   double precision qq(5),xmol(nrel)
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
!    write(*,15)'gd2: ',iph,lokres,qq(1),&
!         ceq%phase_varres(lokres)%gval(1,1),(xmol(i),i=1,nrel)
!15  format(a,2i4,2e12.4,2x,5(F9.5))
   do i=1,nrel
      xarr(i)=real(xmol(i))
   enddo
   if(qq(1).lt.2.0D-1) then
! number of real atoms less than 20%, a gridpoint with just vacancies ....
!      gval=1.0E5
      gval=1.0E1
   elseif(abs(qq(2)).gt.zero) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! the gridpoint has net charge, qq(2), make gval more positive. 
! Note gval(1,1) is divided by RT so around -5<0
! A better method is needed by combining charged gripoints!!!!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+20*qq(2)**2)
!      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+5*qq(2)**2)
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
1010 format('Thread ',i3,', gval: ',1pe15.6,', error: ',i6)
   return
 end subroutine calc_gridpoint

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine generate_fccord_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,ceq)
! This generates grid for a phase with 4 sublattice fcc/hcp ordering
! mode<0 just number of gridpoints in ngg, needed for allocations
! mode=0 calculate mole fraction and G for all gridpoints
! mode>0 return constitution for gridpoint mode in yarr
   implicit none
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
! NOTHING IMPLEMENTED YET
   write(*,*)'FCC/HCP tetraherdal ordering not handelled gracefully'
1000 continue
   return
 end subroutine generate_fccord_grid

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine generate_charged_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,ceq)
! This generates grid for a phase with charged constituents
! mode<0 just number of gridpoints in ngg, needed for allocations
! mode=0 calculate mole fraction and G for all gridpoints
! mode>0 return constitution for gridpoint mode in yarr
   implicit none
   integer mode,iph,ngg,nrel,ny
   real xarr(nrel,*),garr(*)
   double precision yarr(*)
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer nkl(maxsubl),knr(maxconst),inkl(0:maxsubl)
!   double precision, dimension(:), allocatable :: yfra
   double precision sites(maxsubl),ydum(maxconst),qq(5)
   integer nend,ll,nsl,i1,i2,i3,loksp,mm,lokph,lokcs,np,nm,nn,ncc,iz,loopf
   integer, dimension(:,:), allocatable :: neutral
   integer, dimension(:), allocatable :: savengg
!   integer ielno(10)
!   double precision stoi(10),smass,qsp
   double precision charge,ratio1,ratio2
   double precision, dimension(:), allocatable :: y1,y2,y3,y4,y5
   real xdum(nrel),gdum
   integer, parameter :: ncf5=5,ncf3=3
   integer ncf
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
   save savengg
! we will select 5 or 3 gripoints below
!   ncf=ncf5
   if(.not.allocated(savengg)) then
      allocate(savengg(noofph))
      savengg=0
   endif
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
!   write(*,*)'Charged grid for phase ',iph,mode,nend
   if(nend.eq.1) then
! a single endmember, just check it is neutral
      ngg=1
      do ll=1,nsl
         loksp=knr(ll)
         charge=charge+sites(ll)*splista(loksp)%charge
      enddo
      if(charge.eq.zero) then
         np=ngg
         if(mode.eq.0) then
! if mode=0 calculate G for this endmember
!         write(*,*)'3Y a single neutral endmember for ',iph,mode
            call calc_gridpoint(iph,ydum,nrel,xarr(1,ngg),garr(ngg),ceq)
            if(gx%bmperr.ne.0) goto 1000
!         elseif(mode,gt.0) then
! if mode>0 return constitution, already set
         endif
! finally remove the request for external charge balance !!!
!         write(*,*)'No external charge balance for phase:',iph,lokcs,mode
         call get_phase_compset(iph,1,lokph,lokcs)
         phlista(lokph)%status1=ibclr(phlista(lokph)%status1,PHEXCB)
         goto 1000
      endif
      ngg=0
      call get_phase_compset(iph,1,lokph,lokcs)
      write(*,*)'Phase suspended as net charge: ',phlista(lokph)%name
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
!   write(*,10)'3Y nend: ',nend,0.0D0,(nkl(ll),ll=1,nsl)
10 format(a,i3,5x,1pe12.4,10i3)
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
!   write(*,10)'3Y endmem: ',1,charge,endmem(1)%constit
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
   write(*,22)'3Y endmem: ',iph,mode,nend,np,nm,nn,mm
22 format(a,i3,i6,i5,3i4,5i6)
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
   if(mode.ge.0) then
! we have saved the number of gridpoints from the mode=-1 call here
      np=savengg(iph)
!      write(*,*)'3Y allocate neutral',mode,np
      allocate(neutral(np,0:3))
      neutral=0
   endif
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
   if(mode.lt.0) then
! we have just calculated the number of gridpoints, save and exit
!      write(*,*)'3Y neutral gridpoints: ',np
      ngg=np
      savengg(iph)=ngg
   else
! Generate the composition of the gridpoints from 1-3 endmembers and
! if mode=0 calculate the composition and Gibbs energy for the gridpoints
! if mode>0 return the constitution of gridpoint mode.
! How do I know mode is mode gridpoint in this phase??
!      write(*,29)'3Y we are here?',iph,mode,np,nsl,inkl(nsl)
29    format(a,10i5)
      ncc=inkl(nsl)
      allocate(y1(ncc))
      allocate(y2(ncc))
      allocate(y3(ncc))
      allocate(y4(ncc))
! loopf keeps track if several gridpoints belong together
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,zero,y4
300         format(a,i4,2i2,3i3,1pe10.2,8(0pf6.3))
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,zero,y4
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,charge,y4
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,charge,y4
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,charge,y4
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,charge,y4
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,charge,y4
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,charge,y4
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
410            format(a,i4,3i3,6(1pe10.2))
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
!            write(*,300)'3Y gp* ',nm,nn,loopf,i1,i2,i3,charge,y4
!----------------------- 
         end select
!===============================================================
! Here we have the neutral constituent fraction in y4
! if mode>0 we have found the requested constitution
         if(mode.lt.0) then
            write(*,*)'We should never be here ...'
            goto 1000
         elseif(mode.gt.0) then
            if(mode.eq.nm) then
               ny=ncc
               do ll=1,ny
                  yarr(ll)=y4(ll)
               enddo
!               write(*,507)'3Y Solution gp: ',mode,iph,y4
507            format(a,i5,i4,10F7.4)
               goto 1000
            endif
! continue searching for correct gridpoint of the solution
         else
! for mode=0 we must calculate G
! this is just for debugging
!            call set_constitution(iph,1,y4,qq,ceq)
!            if(gx%bmperr.ne.0) goto 1000
!            if(abs(qq(2)).gt.1.0D-6) then
!               write(*,511)'3Y gp with charge: ',nm,iph,nn,qq(2)
511            format(a,i5,i4,i3,1pe12.4)
!            endif
!
!            call calc_gridpoint(iph,y4,nrel,xdum,gdum,ceq)
            call calc_gridpoint(iph,y4,nrel,xarr(1,nm),garr(nm),ceq)
            if(gx%bmperr.ne.0) goto 1000
!         write(*,512)nm,qq(2),gdum,xdum
512         format('3Y gridpoint: ',i4,2(1pe12.4),7(0pF7.4))
         endif
      enddo ygen
!
   endif
1000 continue
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
 subroutine calcg_endmember(iph,endmember,gval,ceq)
! calculates G for one mole of real atoms for a single end member
! used for reference states. Restores current composition (but not G or deriv)
! endmember contains indices in the constituent array, not species index
! one for each sublattice
   implicit none
   integer iph
   double precision gval
   integer endmember(maxsubl)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer ierr,kk0,ll,lokres,nsl
   integer nkl(maxsubl),knr(maxconst)
   double precision savey(maxconst),sites(maxsubl),qq(5),yfra(maxconst)
!
   call get_phase_data(iph,1,nsl,nkl,knr,savey,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1100
! set constitution to be just the endmember
! It is difficult to make this simpler as one can have magnetic contributions
! to G, this it is not suffiecient jyst to calculate the G function, one must
! calculate TC etc.
   yfra=zero
   kk0=0
   do ll=1,nsl
      if(endmember(ll).gt.kk0 .and. endmember(ll).le.kk0+nkl(ll)) then
         yfra(endmember(ll))=one
      else
!         write(*,16)'3Y endmember index outside range',ll,endmember(ll),&
!              kk0,nkl(ll)
16       format(a,10i5)
         gx%bmperr=4160; goto 1100
      endif
      kk0=kk0+nkl(ll)
   enddo
!   write(*,17)'set: ',kk0,(yfra(i),i=1,kk0)
17 format(a,i3,5(1pe12.4))
   call set_constitution(iph,1,yfra,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
   call calcg(iph,1,0,lokres,ceq)
   if(gx%bmperr.ne.0) goto 1000
   if(qq(1).ge.1.0D-3) then
! avoid calculating endmembers with too many vacancies. gval is divided by RT
      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1))
!      write(*,*)'gval: ',gval,qq(1)
   else
!      write(*,*)'End member has no atoms'
      gx%bmperr=4161; goto 1000
   endif
1000 continue
   ierr=gx%bmperr
   if(gx%bmperr.ne.0) gx%bmperr=0
! restore constitution
!   write(*,17)'res: ',kk0,(savey(i),i=1,kk0)
   call set_constitution(iph,1,savey,qq,ceq)
   if(gx%bmperr.ne.0) then
      if(ierr.ne.0) then
         write(*,*)'Double errors in calcg_endmember: ',ierr,gx%bmperr
      endif
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
! Restores current composition (but not G or deriv)
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
   double precision saveg(6)
!
   call get_phase_data(iph,1,nsl,nkl,knr,savey,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1100
! set constitution to be just the endmember
! It is difficult to make this simpler as one can have magnetic contributions
! to G, this it is not sufficient just to calculate the G function, one must
! calculate TC etc.
   yfra=zero
   kk0=0
!   write(*,11)'3Y refstate: ',iph,nsl,nkl(1),endmember(1)
11 format(a,10i5)
   do ll=1,nsl
      if(endmember(ll).gt.kk0 .and. endmember(ll).le.kk0+nkl(ll)) then
         yfra(endmember(ll))=one
      else
!         write(*,16)'3Y endmember index outside range',ll,endmember(ll),&
!              kk0,nkl(ll)
16       format(a,10i5)
         gx%bmperr=4160; goto 1100
      endif
      kk0=kk0+nkl(ll)
   enddo
!   write(*,17)'set: ',kk0,(yfra(i),i=1,kk0)
17 format(a,i3,5(1pe12.4))
   call set_constitution(iph,1,yfra,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! we do not know lokres here !!
   ics=1
   call get_phase_compset(iph,ics,lokph,lokres)
   if(gx%bmperr.ne.0) goto 1000
   do ll=1,6
      saveg(ll)=ceq%phase_varres(lokres)%gval(ll,1)/qq(1)
   enddo
!   write(*,432)saveg
432 format('3Y::',6(1pe12.4))
! third argument to calcg is 2 to calculate all derivatives
   call calcg(iph,1,2,lokres,ceq)
   if(gx%bmperr.ne.0) goto 1000
   if(qq(1).ge.1.0D-2) then
! avoid calculating endmembers with too many vacancies. gval is divided by RT
! gval(1..6,1) are G, G.T, G.P, G.T.T, G.T.P and G.P.P
      do ll=1,6
         gval(ll)=ceq%phase_varres(lokres)%gval(ll,1)/qq(1)
      enddo
!      write(*,*)'gval: ',gval,qq(1)
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
         write(*,11)'3Y endmember index outside range',ll,endmember(ll),nkl(ll)
11       format(a,10i4)
         gx%bmperr=4160; goto 1100
      endif
      kk0=kk0+nkl(ll)
   enddo
!   write(*,17)'set: ',kk0,(yfra(i),i=1,kk0)
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
!   write(*,17)'res: ',kk0,(savey(i),i=1,kk0)
   call set_constitution(iph,1,savey,qq,ceq)
   if(gx%bmperr.ne.0) then
      if(ierr.ne.0) then
         write(*,*)'Double errors: ',ierr,gx%bmperr
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
   double precision, parameter :: phfmin=1.0D-8
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
   logical checkremoved
! if trace then open file to write grid
   if(trace) then
      write(*,*)'Writing grid solution on file ocgrid.dat'
      open(31,file='ocgrid.dat ',access='sequential')
      write(31,700)nrel,kp,(xknown(inuse),inuse=1,nrel)
700   format('Output from OC gridmin'/' Elements: ',i2,', gridpoints: 'i5,&
           ', composition: '/6(F7.4))
      write(31,*)' Gridpoints: '
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
!   write(*,11)'fm8: ',(xknown(i),i=1,nrel)
!11 format(a,7(F8.4))
! Find the lowest Gibbs energy as close as possible to each pure element
! or with max content
   nopure=0
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
!            write(*,*)'pure: ',je,ip,gmin(je)
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
! check that we have nrel gridpoints
   do je=1,nrel
      if(jgrid(je).eq.0) then
! no gridpoint assigned to this element!! error (note C in pure fcc has no gp)
!         gx%bmperr=4149; goto 1000
         write(*,122)'Warning, no gridpoint for pure element ',je
!              nopure(je),xmaxx(je)
122      format(a,2i5,2F7.4)
         if(nopure(je).eq.0) then
!            write(*,122)'No solubility in any phase for element ',je
            gx%bmperr=4149; goto 1000
         elseif(xarr(je,nopure(je)).gt.xknown(je)) then
! accept gripoint with highest content of element je outside known composition
            do ie=1,nrel
               xmat(ie,je)=xarr(ie,nopure(je))
            enddo
            gmin(je)=garr(ip)
            phfrac(je)=xknown(je)
         else
            write(*,122)'Composition outside phase compositions for element',&
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
! output of start matrix
!   do ip=1,nrel
!      write(*,123)ip,phfrac(ip),(xmat(je,ip),je=1,nrel)
!   enddo
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
!      write(*,177)'fm4: ',jgrid(je),phfrac(je),(xmat(ie,je),ie=1,nrel)
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
         write(31,720),xknown(je),xknown(je),(xmat(ie,je),ie=1,nrel)
      enddo
720   format('3Y& ',2F7.4,1x,8f8.5)
      write(31,730)gvvp,(cmu(je),je=1,nrel)
730   format('3Y Gibbs energy: ',1pe14.6/'Chemical potentials: '/6(1pe12.4))
   endif
   griter=0
   gpfail=0
!   write(*,175)'ini: ',gvvp,(cmu(ie),ie=1,nrel)
175 format(a,(1e12.4),2x,6(1pe12.4))
!   write(*,*)'gvvp: ',gvvp
! check we have the correct global composition
!    call chocx('fgm1 ',nrel,jgrid,phfrac,xmat)
!    if(gx%bmperr.ne.0) goto 1000
!    write(*,173)gvvp,(jgrid(i),i=1,nrel)
173 format('fms: ',1pe12.4,10i5)
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
!   write(*,*)'Gridpoints in use: ',inuse
!   ovall=zero
!   do i=1,nrel
!      ovall=ovall+xknown(i)*cmu(i)
!   enddo
!   write(*,203)'ff:',inuse,ovall,(cmu(je),je=1,nrel)
203 format(a,i4,1pe12.4,6(1pe11.3))
   pointloop: do jp=1,kp
      included: if(notuse(jp).eq.0) then
         gplan=zero
! first index in xarr is component, second is gridpoint
         do iel=1,nrel
            gplan=gplan+xarr(iel,jp)*cmu(iel)
         enddo
         dg=garr(jp)-gplan
!         write(*,209)'fmz: ',dg,garr(jp),gplan
209      format(a,3(1pe12.4))
         if(dg.gt.zero) then
!            inuse=inuse-1
! we cannot be sure that a point that has a positive value will always be
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
!         write(*,*)'Excluded: ',griter,jp
      endif included
   enddo pointloop
!   write(*,*)'3Y Finished loop for all gridpoints: ',jp,kp
! TBASE bug------------------------
!   jp=94
!   do iel=1,nrel
!      gplan=gplan+xarr(iel,jp)*cmu(iel)
!   enddo
!   dg=garr(jp)-gplan
!   write(*,7677)jp,gplan,garr(jp),dg,(xarr(iel,jp),iel=1,nrel)
!7677 format('Gridpoint: ',i5,3(1pe12.4)/(10f7.4))
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
!   write(*,211)'ny:',nyp,dgmin,(xarr(ie,nyp),ie=1,nrel)
   if(trace) write(*,212)'3Y Found gridpoint ',nyp,inuse,dgmin,garr(nyp)
!211 format(a,i4,1pe12.4,0pf8.4,6(f8.4))
212 format(a,2i8,6(1pe11.3))
!-------------------------------------------------------------------------
! Replace one point with the new nyp keeping the overall composition inside.
! This is done by replacing one row at a time in xmat and solve a linear
! equation for the phase fractions and accept the only solution which has 
! positive phasefractions.
!
! check the overall composition
!   xov=zero
!   do i=1,nrel
!      do ie=1,nrel
!         xov(ie)=xov(ie)+phfrac(i)*xmat(ie,i)
!      enddo
!   enddo
!   write(*,277)'xov:   ',nyp,(xov(i),i=1,nrel)
!   write(*,277)'phfrac:',0,(phfrac(i),i=1,nrel)
!277 format(a,i5,10F7.4)
!   
!-----------------------------------------------------------------------
   qmat=qmatsave
   do i=1,nrel
      phfsave(i)=phfrac(i)
   enddo
   ie=0
! loop to try to replace an old gridpoint by nyp.  Try to replace all.
300 continue
   ie=ie+1
   if(ie.gt.nrel) then
! tried to change all coumns but no solution, error
!      write(*,301)'Failed gp: ',nyp,inerr,dgmin,(xarr(i,nyp),i=1,nrel)
301   format(a,i7,i3,1pe10.2,2x,6(0pF7.4))
      gpfail=gpfail+1
! listing restored solution ......
!      xtx=zero
!      do jjq=1,nrel
!         write(*,177)'flp: ',jgrid(jjq),phfrac(jjq),(xmat(ie,jjq),ie=1,nrel)
!         do jjz=1,nrel
!            xtx(jjz)=xtx(jjz)+phfrac(jjq)*xmat(jjz,jjq)
!         enddo
!      enddo
!      gvv=zero
!      do jjq=1,nrel
!         gvv=gvv+xtx(jjq)*cmu(jjq)
!      enddo
!      write(*,175)'cur: ',gvv,(cmu(ie),ie=1,nrel)
!
! >>>> problem with gas phase test case cho1 with x(c)=.2 x(o)=x(H)=.4
! The gridpoints returned not good, probably due to too many gridpoints ...
!
!      if(trace) write(*,*)'Failed when trying to add gridpoint ',nyp
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
! left side are the known composition
   do je=1,nrel
      qmat(je,nrel1)=xknown(je)
   enddo
! solver, note qmat is destroyed inside lingld, nrel is dimension
! qmat matrix with left hand side as additional column i.e. QMAT(1..ND1,ND2)
! phfrac(ND1) is result array, nz number of unknown, ierr nonzero=error
!    do ik=1,nrel1
!       write(*,317)'fm6A: ',(qmat(je,ik),je=1,nrel)
!    enddo
   call lingld(nrel,nrel1,qmat,phfrac,nrel,ierr)
   if(ierr.ne.0) then
! error may occur and are not fatal, just try to replace next column
!       write(*,*)'non-fatal error from lingld: ',ierr,nyp
      qmat=qmatsave
      do i=1,nrel
         phfrac(i)=phfsave(i)
      enddo
      goto 300
   endif
!   write(*,*)'fm6B: ',ie,ierr
!   write(*,317)'fm6C: ',(phfrac(i),i=1,nrel)
317 format(a,6(1pe12.4))
!-----------------------
! if solution has only positive values accept this, ierr nonzero if singular
   do je=1,nrel
      if(phfrac(je).le.phfmin .or. phfrac(je).gt.one) then
! maybe problems if known composition have almost zero of some components?
! restore qmat
!          write(*,*)'fm6D: ',je
         qmat=qmatsave
         do i=1,nrel
            phfrac(i)=phfsave(i)
         enddo
         goto 300
      endif
   enddo
!   if(trace) write(*,*)'Replaced column: ',ie,nyp
! we have found that column ie should be replaced
!--------------------------------------------------
! update xmat, qmatsave and gmin
! as we may fail to find the solution for the chemical potentials later
! keep a copy that can be restored
   iesave=ie
   jsave=jgrid(iesave)
! mark that the replaced gridpoint should be checked again ....
!   write(*,*)'3Y Putting gridpoint back: ',jgrid(ie)
   notuse(jgrid(ie))=0
   jgrid(ie)=nyp
   xmatsave=xmat
   do je=1,nrel
      xmat(je,ie)=xarr(je,nyp)
      qmatsave(je,ie)=dble(xarr(je,nyp))
   enddo
   gmin(ie)=garr(nyp)
!   do ik=1,nrel
!      write(*,317)'fm6F: ',(xmat(je,ik),je=1,nrel)
!   enddo
!   write(*,317)'fm6G: ',(gmin(je),je=1,nrel)
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
!       write(*,317)'fm8A: ',(zmat(je,ik),je=1,nrel)
!    enddo
   cmusave=cmu
   call lingld(nrel,nrel1,zmat,cmu,nrel,ierr)
   if(ierr.ne.0) then
! this should also be handelled by ignoring the new gridpoint but
! here we must restore the xmat, qmatsave and cmu.
!      write(*,*)'Failed to calculate chemical potentials',ierr
!      if(trace) write(*,*)'Error from LINGLD for chem.pot.: ',ierr,nyp
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
! calculate total G
!   gvv=zero
!   do ie=1,nrel
!      do je=1,nrel
! first index is component, second is species
!         gvv=gvv+xmat(je,ie)*cmu(je)
!      enddo
!   enddo
!   if(trace) write(*,*)'New total G: ',gvv,gvvp
! check if gvv is lower than previous
!   if(gvv.gt.gvvp) then
!      write(*,*)' *** Gibbs energy increased, restore!'
!   endif
!   gvvp=gvv
!----------------------------------------------------------
! debug output as we have changed one gridpoint
!   xtx=zero
!   do jjq=1,nrel
!      write(*,177)'gpf: ',jgrid(jjq),phfrac(jjq),(xmat(ie,jjq),ie=1,nrel)
!      do jjz=1,nrel
!         xtx(jjz)=xtx(jjz)+phfrac(jjq)*xmat(jjz,jjq)
!      enddo
!   enddo
!   gvv=zero
!   do jjq=1,nrel
!      gvv=gvv+xtx(jjq)*cmu(jjq)
!   enddo
!   write(*,175)'ny4: ',gvv,(cmu(ie),ie=1,nrel)
!   write(*,317)'new cmu: ',(cmu(je),je=1,nrel)
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
         write(*,720),phfrac(je),xknown(je),(xmat(je,ie),ie=1,nrel)
      enddo
      write(31,730)gvvp,(cmu(je),je=1,nrel)
   endif
   if(checkremoved) then
      write(*,198)nyp
198   format('Added previously removed gridpoint ',i6)
      goto 950
   endif
!----------------------------------------------
! here we go back to loop through all gridpoints again
!   write(*,*)'New search: ',griter
   goto 200
!==============================================
900 continue
   if(gpfail.gt.0) then
      write(*,906)gpfail
906   format('Failed using ',i7,' gridpoints')
   endif
!   write(*,*)'Gridmin has found a solution'
!   write(*,316)'fm9A: ',(jgrid(i),i=1,nrel)
!   do ik=1,nrel
!      write(*,317)'fm9B: ',(xmat(je,ik),je=1,nrel)
!   enddo
!   write(*,317)'fm9C: ',(garr(je),je=1,nrel)
!   write(*,317)'fm9D: ',(cmu(je),je=1,nrel)
!   write(*,317)'fm9E: ',(phfrac(je),je=1,nrel)
316 format(a,10i5)
   nj=0
!    do j=1,jerr
!       if(removed(j).gt.0) then
!          write(*,*)'Failed testing gridpoint ',removed(j)
!          nj=nj+1
!       endif
!    enddo
950 continue
   nj=0
   checkremoved=.true.
!   write(*,*)'Checking removed gridpoints',inerr
!   xtx=zero
!   do jjq=1,nrel
!      write(*,177)'flp: ',jgrid(jjq),phfrac(jjq),(xmat(ie,jjq),ie=1,nrel)
!      do jjz=1,nrel
!         xtx(jjz)=xtx(jjz)+phfrac(jjq)*xmat(jjz,jjq)
!      enddo
!   enddo
!   gvv=zero
!   do jjq=1,nrel
!      gvv=gvv+xtx(jjq)*cmu(jjq)
!   enddo
!   write(*,175)'cur: ',gvv,(cmu(ie),ie=1,nrel)
!----------------
   testloop: do jj=1,inerr
      jp=removed(jj)
!      write(*,*)'Checking removed gridpoint: ',jj,jp
      if(jp.gt.0) then
         gplan=zero
         do iel=1,nrel
            gplan=gplan+xarr(iel,jp)*cmu(iel)
         enddo
         dg=garr(jp)-gplan
         if(dg.lt.zero) then
!            if(trace) write(*,985)jp,dg,garr(jp),gplan
!            write(*,982)jp,dg,garr(jp),gplan
982         format('Removed gridpoint ',i5,' is below surface ',3(1pe12.4))
! try to include it ....
            ie=0
            removed(jj)=-jp
            nyp=jp
            goto 300
         else
!            write(*,983)jp,dg
983         format('Removed gridpoint ',i5,' above surface ',1pe12.4)
            removed(jj)=-jp
         endif
      endif
   enddo testloop
   if(inerr.gt.0 .and. nj.eq.0) then
!      if(trace) write(*,986)inerr
986   format('None of the ',i3,' removed gridpoints below final surface')
   endif
   if(trace) write(*,771)(jgrid(je),je=1,nrel)
771 format('Final set of gridpoints: '(/15i5))
!   xtx=0
!   do iii=1,nrel
!      write(*,987)jgrid(iii),phfrac(iii),(xarr(i,jgrid(iii)),i=1,nrel)
!987   format('GP: ',i5,F7.4,2x,6F9.6)
!      do j=1,nrel
!         xtx(j)=xtx(j)+phfrac(iii)*xarr(j,jgrid(iii))
!      enddo
!   enddo
!   write(*,988)(xtx(i),i=1,nrel)
!988 format('MF: ',6F9.6)
!
!    call chocx('fgme ',nrel,jgrid,phfrac,xmat)
1000 continue
   return
 end subroutine find_gridmin

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine merge_gridpoints(nv,iphl,aphl,nyphl,yphl,trace,nrel,xsol,cmu,ceq)
!
! BEWARE not adopted for parallel processing
!
! if the same phase occurs several times check if they are really separate
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
   double precision ycheck(maxconst),qq(5),xerr(maxel)
   double precision summu,sumam
   logical igen
   real xmix(maxel)
   double precision a1,a2,gdf,gval1,gval2,gval3,gmindif
!
! gmindif is the value to accept to merge two gridpoints
! It should be a variable that can be set by the user for finetuning
!   write(*,7)'Merge_gridpoints is dissabled for the moment',nv
!7  format(a,i3)
! NOTE, merging gripoints in ideal phases like gas
   goto 1100
   if(ocv()) write(*,*)'Entering merge_gridpoints'
!---------------------
   gmindif=ceq%gmindif
   notuse=0
   nm=0
   incy(1)=1
   do i=2,nv
      incy(i)=incy(i-1)+nyphl(i-1)
   enddo
   summu=zero
   xerr=zero
   do jp=1,nv
      summu=summu+aphl(jp)
      do i=1,nrel
         xerr(i)=xerr(i)+aphl(jp)*xsol(i,jp)
      enddo
   enddo
   write(*,73)'in: ',summu,(xerr(i),i=1,nrel)
73 format(a,F5.2,2x,10f7.4)
!----------------------------------------------
100 continue
   igen=.false.
   do jp=1,nv-1
      if(notuse(jp).ne.0) goto 400
      do kp=jp+1,nv
         if(notuse(kp).ne.0) goto 300
         if(iphl(jp).eq.iphl(kp)) then
            iph=iphl(jp)
!            write(*,9876)'XP1: ',jp,(xsol(i,jp),i=1,nrel)
!            write(*,9876)'XP2: ',kp,(xsol(i,kp),i=1,nrel)
9876        format(a,i4,5(1pe12.4))
! same phase in two points, see if they are really separate
! the test is simple, just calculate the Gibbs energy at the weighted
! average and if that has lower gibbs energy then merge them
!            write(*,130)'c130: ',jp,incy(jp),nyphl(jp),kp,incy(kp),nyphl(kp)
130         format(a,10i5)
!            write(*,140)incy(jp),(yphl(incy(jp)+j),j=0,nyphl(jp)-1)
140         format(i4,6(1pe12.4))
            call set_constitution(iph,1,yphl(incy(jp)),qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call calcg(iph,1,0,lokres,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gval1=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
! debug
            call calc_phase_mol(iph,xerr,ceq)
!            write(*,79)'Y0: ',(yphl(incy(jp)+i-1),i=1,nyphl(jp))
!            write(*,79)'X1: ',(xerr(i),i=1,nrel)
!            write(*,79)'X2: ',(xsol(i,jp),i=1,nrel)
79          format(a,12(F6.3))
! Subtract the solution, the result should be zero ??
! The mole fractions of the gridpoints in solution is in xsol(1,jgrid(i))
            summu=zero
            do jj=1,nrel
               summu=summu+xsol(jj,jp)*cmu(jj)
            enddo
            gval1=gval1-summu
! debug output gridpoint 1
            mj=nyphl(jp)-1
!           write(*,820)'GP1:',gval1,summu,aphl(jp),(yphl(incy(jp)+jj),jj=0,mj)
820         format(a,3(1pe10.2),10(0pF5.2))
            call set_constitution(iph,1,yphl(incy(kp)),qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call calcg(iph,1,0,lokres,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gval2=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
! Subtract the solution
! The mole fractions of the gridpoints in solution is in xsol(1,jgrid(i))
            summu=zero
            do jj=1,nrel
               summu=summu+xsol(jj,kp)*cmu(jj)
            enddo
            gval2=gval2-summu
! debug output gridpoint 2
!           write(*,820)'GP2:',gval2,summu,aphl(kp),(yphl(incy(kp)+jj),jj=0,mj)
! select the middle point
! Not very good weight by phase amounts if one has 95% FCC and 5 % MC ....
! take weighted sum of composition in the middle
            a1=5.0D-01
            a2=5.0D-01
!            sumam=aphl(jp)+aphl(kp)
!            a1=aphl(jp)/sumam
!            a2=aphl(kp)/sumam
! SURPRISE: adding together constituent fractions does not reproduce
! the correct molefractions if the constituents are molecules ....
            do i=1,nyphl(jp)
               ycheck(i)=a1*yphl(incy(jp)+i-1)+a2*yphl(incy(kp)+i-1)
            enddo
            call set_constitution(iph,1,ycheck,qq,ceq)
            if(gx%bmperr.ne.0) goto 1000
            call calcg(iph,1,0,lokres,ceq)
            if(gx%bmperr.ne.0) goto 1000
            gval3=ceq%phase_varres(lokres)%gval(1,1)/qq(1)
! Check if this is below the current tangent plane
! The mole fractions of the gridpoints in solution is in xsol(1,jgrid(i))
            summu=zero
            do jj=1,nrel
               xmix(jj)=a1*xsol(jj,jp)+a2*xsol(jj,kp)
               summu=summu+xmix(jj)*cmu(jj)
            enddo
!            write(*,21)'summu: ',a1,a2,qq(1),gval3,summu,gmindif
            gval3=gval3-summu
21          format(a,6(1pe12.4))
            write(*,22)'dg: ',gval3,gval1,gval2,0.5*gval1+0.5*gval2
22          format(a,6(1pe12.4))
! debug output mix
            write(*,820)'GPY:',gval3,summu,qq(1),(ycheck(jj+1),jj=0,mj)
            write(*,820)'GPX:',gval3,summu,qq(1),(xmix(jj),jj=1,nrel)
            gdf=gval3-a1*gval1-a2*gval2
!            gdf=gval3
! merge require that difference is less than gmindif or phase ideal
! allways merge ideal phase as it never has miscibility gaps !!!
            lokph=phases(iph)
            if(gdf.lt.gmindif .or. &
                 btest(phlista(lokph)%status1,PHID)) then
! gridpoint in between has lower G, merge
               write(*,830)'merged:     ',jp,kp,gdf,iphl(jp),&
                    aphl(jp)+aphl(kp)
830            format('Gridpoints ',a,2i3,1pe15.4,' in phase ',i3,1pe12.4)
!               write(*,840)jp,(xsol(jj,jp),jj=1,nrel)
!               write(*,840)kp,(xsol(jj,kp),jj=1,nrel)
!               write(*,840)jp,(xmix(jj),jj=1,nrel)
840            format('x: ',i5,10(F7.4))
! If merging use correct phase amounts
               a1=aphl(jp)/(aphl(jp)+aphl(kp))
               a2=aphl(kp)/(aphl(jp)+aphl(kp))
!               write(*,160)iph,jp,kp,incy(jp),nyphl(jp),gdf,&
!                    gval1,gval2,gval3,a1,a2
160            format('Merging:  ',i3,2x,4i5,1pe12.4/5(1pe12.4))
               write(*,162)'p1:',a1,(yphl(incy(jp)+j),j=0,nyphl(jp)-1)
162            format(a,F5.2,2x,(10f7.4))
               write(*,162)'p2:',a2,(yphl(incy(kp)+j),j=0,nyphl(kp)-1)
! The gridpoint jp has new amount, composition and constitution
! SURPRISE: adding together constituent fractions does not reproduce
! the correct molefractions if the constituents are molecules ....
               aphl(jp)=aphl(jp)+aphl(kp)
               do i=0,nyphl(jp)-1
                  yphl(incy(jp)+i)=a1*yphl(incy(jp)+i)+a2*yphl(incy(kp)+i)
               enddo
               call set_constitution(iph,1,yphl(incy(jp)),qq,ceq)
               if(gx%bmperr.ne.0) goto 1000
! extract correct mole fractions
               call calc_phase_mol(iph,xerr,ceq)
               write(*,162)'ym:',0.0D0,(yphl(incy(jp)+i),i=0,nyphl(jp)-1)
               write(*,162)'xj:',0.0D0,(xsol(jj,jp),jj=1,nrel)
               write(*,162)'xk:',0.0D0,(xsol(jj,kp),jj=1,nrel)
               write(*,162)'xy:',0.0D0,(xerr(i),i=1,nrel)
               do i=1,nrel
                  xsol(i,jp)=xerr(i)
               enddo
               igen=.true.
               nm=nm+1
               iphl(kp)=-iphl(kp)
               notuse(kp)=1
! check overall composition of solution ...
               summu=zero
               xerr=zero
               do i=1,nrel
                  if(iphl(i).lt.0) cycle
                  summu=summu+aphl(i)
                  write(*,*)'point: ',i,aphl(i)
                  do jj=1,nrel
                     xerr(jj)=xerr(jj)+aphl(i)*xsol(jj,i)
                  enddo
               enddo
               write(*,73)'nu: ',summu,(xerr(jj),jj=1,nrel)
! the chemical potentials has changed but how?  Approximate the change by
! making gmindif more negative for each merge (does not affect ideal phases)
               gmindif=2.0D0*gmindif
! after merging always restart loop
               goto 100
            else
               write(*,830)'not merged: ',jp,kp,gdf,iphl(jp),gmindif
            endif
         endif
300      continue
      enddo
400   continue
   enddo
! if two gridpoints merged compare all again
   if(igen) goto 100
!----------------------------------------
! shift fractions for the removed phases
450 continue
!   write(*,*)'at label 450: ',nm
   klast=0
   do jp=1,nv
      klast=klast+nyphl(jp)
   enddo
!
! uncomment listing here if error moving fractions
!    write(*,502)nv,(iphl(i),i=1,nv)
!    write(*,502)0,(incy(i),i=1,nv)
!    write(*,502)klast,(nyphl(i),i=1,nv)
502 format('check1: ',i3,2x,20i4)
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
503      format('check3: ',5i5)
!          write(*,555)'nyy1: ',(yphl(ip),ip=kk+1,kk+jump)
555      format(a,6(1pe12.4))
         do ip=kk+1,klast-jump
            yphl(ip)=yphl(ip+jump)
         enddo
!          write(*,555)'nyy2: ',(yphl(ip),ip=kk+1,kk+jump)
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
!
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
   if(ocv()) write(*,*)'At return from merge_gridpoints: ',nv
   return
!------------------------------------------
! emergency fix to avoid creating several composition sets in ideal phases
1100 continue
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
!               write(*,*)'merging gridpoints in ideal phase'
               sumam=aphl(jp)+aphl(kp)
               a1=aphl(jp)/sumam
               a2=aphl(kp)/sumam
               aphl(jp)=aphl(jp)+aphl(kp)
!               write(*,1117)'3Y: ',jp,(yphl(incy(jp)+i),i=0,nyphl(jp)-1)
!               write(*,1117)'3Y: ',kp,(yphl(incy(kp)+i),i=0,nyphl(kp)-1)
               do i=0,nyphl(jp)-1
                  yphl(incy(jp)+i)=a1*yphl(incy(jp)+i)+a2*yphl(incy(kp)+i)
               enddo
!               write(*,1117)'3Y: ',jp,(yphl(incy(jp)+i),i=0,nyphl(jp)-1)
1117           format(a,i3,6(1pe12.4))
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
! kp    total number of gridpoints
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
   integer ig1,ign,ip,iph,ics,jph,lokcs,lokph,mode,ny,ie,ig
   double precision yarr(maxconst),qq(5),xxx,dgmin
   real dg,gplan
   if(ocv()) write(*,*)'Entering set_metastable'
! loop through the gridpoints for all unstable phases and insert the
! stable constitution that is closest to be stable
!   write(*,7)'set_meta: ',kp,nrel,nr,(iphl(i),i=1,nr)
!7  format(a,i9,2i4,2x,10i3)
!   do i=1,noofph
!      write(*,*)'grid: ',i,kphl(i),ngrid(i)
!   enddo
   phloop: do iph=1,noofph
      do jph=1,nr
         if(iph.eq.iphl(jph)) goto 500
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
      cycle
! this phase is not suspended and not stable, find gridpoints
60    continue
      ig1=kphl(iph)
      ign=ngrid(iph)
      if(ocv()) write(*,69)'Searching gridpoints for: ',iph,ics,ig1,ign
69    format(a,2(i3,1x),2x,3(i6,1x))
! if ig1=0 there are no gridpoints for this phase, it is suspended or dormant
      if(ig1.le.0) cycle
      dgmin=-1.0d12
      ip=0
! search for gripoint closeset to stable plane
      do ig=ig1,ign
         gplan=zero
         do ie=1,nrel
            gplan=gplan+xarr(ie,ig)*cmu(ie)
         enddo
         dg=gplan-garr(ig)
!         write(*,74)'dgx: ',ig,dg,dgmin,(xarr(i,ig),i=1,nrel)
!74       format(a,i5,2(1pe12.4),2x,6(0pf7.4))
         if(dg.gt.dgmin) then
            ip=ig
            dgmin=dg
!            write(*,77)'lower: ',ig,gplan,garr(ig),dg,dgmin
!77          format(a,i7,4(1pe12.4))
         endif
      enddo
      if(ocv()) write(*,79)'Least unstable gridpoint: ',iph,ics,ig1,ign,dgmin
79    format(a,4(i6,1x),1pe12.4)
!      if(ip.eq.0 .or. dgmin.gt.zero) then
!         write(*,*)'This gridpoint stable: ',ip,dgmin
!         write(*,*)'data: ',ip,dgmin
!      endif
! retrieve constitution for this gridpoint and insert it in phase
! must provide mode and iph. The subroutine returns ny and yarr
! mode is the gridpoint in the phase, subtract ig1-1
      mode=ip-ig1+1
! 
      if(ocv()) write(*,78)'calling gengrid: ',iph,ig1,ip,ign,mode,dgmin
78    format(a,5i7,1pe12.4)
! find the constitution of this gridpoint
      call generate_grid(mode,iph,ign,nrel,xarr,garr,ny,yarr,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,451)(yarr(i),i=1,ny)
451   format('fractions: ',6(F10.6))
      call set_constitution(iph,1,yarr,qq,ceq)
      if(gx%bmperr.ne.0) goto 1000
!      write(*,452)iph,ics,(yarr(i),i=1,ny)
452   format('my: ',2i2,6(1pe10.3))
! set driving force also ...
      call set_driving_force(iph,1,dgmin,ceq)
500   continue
   enddo phloop
1000 continue
   if(ocv()) write(*,*)'Finished set_metastable'
   return
 end subroutine set_metastable_constitutions
 
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatiom}
 subroutine gridmin_check(nystph,kp,nrel,xarr,garr,xknown,ngrid,pph,&
      cmu,yphl,iphx,ceq)
! This subroutine checks if a calculated solution is correct by  
! checking if there are any gridpoints below the surface defined
! by the chemical potentials cmu
! nystph return 0 or 10*(phase number)+compset number for new stable phase
! there are kp gridpoints, nrel is number of components
! composition of each gridpoint in xarr, G in garr
! xknown is the known overall composition
! ngrid is last  calculated gridpoint point for a phase jj
! pph is number of phases for which there is a grid
! iphx is phase numbers
! cmu are the final chemical potentials
! yphl is just needed as a dummy
! ceq is current equilibrium record
!\end{verbatiom}
   implicit none
   integer kp,nrel,jp,ie,mode,pph,nystph
   double precision, parameter :: phfmin=1.0D-8
   real xarr(nrel,*),garr(*)
   double precision xknown(*)
   integer, dimension(*) :: ngrid,iphx
   double precision cmu(*),gsurf,gstable,gd,yphl(*),qq(5),rtn,gdmin
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer lokph,lokcs,zph,ibias,ics,iph,ny
! setting a value of addph forces that gridpoint to be added, used for test
!   integer :: addph=8 
!   integer :: addph=100
   integer :: addph=0
   save addph
!
   write(*,*)'Entering gridmin_check',addph
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
!      write(*,17)'grid comarison: ',gd,garr(jp),gsurf
!17    format(a,3(1pe12.4))
      if(gd.gt.gstable) then
! this gridpoint should be set as stable and recalculate
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
!      write(*,*)'mode, ibias and phase: ',mode,ibias,iphx(zph)
      call generate_grid(mode,iphx(zph),ibias,nrel,xarr,garr,ny,yphl,ceq)
      if(gx%bmperr.ne.0) goto 1000
      iph=iphx(zph)
      write(*,*)'Gridmin check found new stable phase: ',iph
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
222      format('added: ',3i3,10(f6.3))
         call set_constitution(iph,ics,yphl,qq,ceq)
      else
         ics=ics+1
         if(ics.gt.phlista(lokph)%noofcs) then
! create new composition set if allowed
            if(btest(globaldata%status,GSNOACS)) then
               gx%bmperr=4177; goto 1000
            endif
            write(*,*)'Creating new composition set for ',iph
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
      write(*,*)'DG at least unstable gridpoint: ',gdmin
   endif
1000 continue
   return
 end subroutine gridmin_check

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
            write(*,*)'Interaction more than 5 levels deep!'
            gx%bmperr=7777; goto 1000
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
!   write(*,*)'In allowenter: ',mode
   logical yesorno
   yesorno=.FALSE.
   if(mode.le.0 .or. mode.gt.3) goto 1000
   if(mode.eq.1) then
! enter element of species not allowed after entering first phase
      if(noofph.gt.1) goto 1000
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
   implicit none
   integer typ
   character name*(*)
!\end{verbatim}
   character name2*64,ch1*1
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
!      write(*,*)'Wrong first letter of symbol: ',name2(1:1),':',name2(1:5)
      gx%bmperr=4137; goto 1000
   endif
   jl=1
100 continue
   jl=jl+1
   ch1=name2(jl:jl)
! always finish when fining a space
   if(ch1.eq.' ') then
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
        notpf(),csfree-1,eqfree-1,nsvfun,reffree-1,addrecs
1007 format('Records for elements, species, phases:           ',3i5/&
            'end members, interactions, properties:           ',3i5/&
            'TP-funs, composition sets, equilibria:           ',3i5/&
            'state variable functions, references, additions: ',3i5)
!----------------------------
! csfree, free list is in firsteq
600 continue
   write(lut,610)csfree,highcs
610 format('Phase_varres free list: ',2i5)
   if(csfree.lt.highcs) then
      lok=csfree
620   continue
      last=lok
      lok=firsteq%phase_varres(last)%nextfree
      write(*,*)'csfree: ',last,lok
      if(lok.le.0 .or. lok.gt.highcs) then
         write(lut,*)'Error in phase_varres free list',last,lok
         goto 1000
      elseif(lok.eq.highcs) then
         goto 630
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
 subroutine enter_default_constitution(iph,ics,mmyfr,ceq)
! user specification of default constitution for a composition set
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer iph,ics
   real mmyfr(*)
!\end{verbatim}
   integer lokph,lokcs,jl,jk
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   jk=size(ceq%phase_varres(lokcs)%yfr)
!   write(*,909)lokph,lokcs,phlista(lokph)%tnooffr,ceq%eqno,&
!        size(ceq%phase_varres),size(ceq%phase_varres(lokcs)%mmyfr),jk
909 format('3Y 2699: ',10i4)
!   write(*,46)'3Y y: ',(ceq%phase_varres(lokcs)%yfr(jl),jl=1,jk)
46 format(a,10(F7.3))
   do jl=1,phlista(lokph)%tnooffr
      ceq%phase_varres(lokcs)%mmyfr(jl)=mmyfr(jl)
!      write(*,47)'3Y jl: ',jl,mmyfr(jl),&
!           firsteq%phase_varres(lokcs)%mmyfr(jl),&
!           ceq%phase_varres(lokcs)%mmyfr(jl)
   enddo
47 format(a,i2,10F7.3)
! set bit indicating that this composition set has a default constitution
!   write(*,*)'3Y enter_default_constitution?? ',lokcs
   ceq%phase_varres(lokcs)%status2=&
        ibset(ceq%phase_varres(lokcs)%status2,CSDEFCON)
1000 continue
   return
 end subroutine enter_default_constitution

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
! set the current constitution of iph composition set ics to its
!  default constitution (if any).  Do not change the amounts of the phases
   implicit none
   integer iph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer lokph,lokcs,ll,jj,kk,kk0
   type(gtp_phase_varres), pointer :: cset
   double precision, allocatable :: yarr(:)
   double precision sum, qq(5)
!
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
   cset=>ceq%phase_varres(lokcs)
! we must use set_constitution at the end to update various internal variables
   allocate(yarr(phlista(lokph)%tnooffr))
   if(allocated(cset%mmyfr)) then
! there is a preset default constitution
      kk=0
      subl1: do ll=1,phlista(lokph)%noofsubl
         kk0=kk
         sum=zero
         if(phlista(lokph)%nooffr(ll).gt.1) then
            do jj=1,phlista(lokph)%nooffr(ll)
! mmy(kk) is negative for small fractions with a maxium, set to 0.01
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
! mmy(kk) is negative for small fractions with a maxium, set to 0.01
               kk=kk+1
               yarr(kk)=yarr(kk)/sum
            enddo
         else
! a single constituent, just increment kk and leave fraction as unity
            kk=kk+1
            yarr(kk)=one
         endif
      enddo subl1
   else
! there is no default constitution, set equal amount of all fractions
      kk=0
      subl2: do ll=1,phlista(lokph)%noofsubl
         if(phlista(lokph)%nooffr(ll).gt.1) then
! set equal amount of all fractions
            sum=one/real(phlista(lokph)%nooffr(ll))
            do jj=1,phlista(lokph)%nooffr(ll)
               kk=kk+1
               yarr(kk)=sum
            enddo
         else
! a single constituent, just increment kk and leave fraction as unity
            kk=kk+1
            yarr(kk)=one
         endif
      enddo subl2
   endif
!   write(*,411)yarr
411 format('3Y set_def_const: ',8F7.4,(10f7.4))
! in this routine the fractions in each sublattice is normallized to be unity
   call set_constitution(iph,ics,yarr,qq,ceq)
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
!   write(*,*)'Todo_before ... not implemented'
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
         ceq%phase_varres(lokcs)%status2=&
              ibclr(ceq%phase_varres(lokcs)%status2,CSSTABLE)
      enddo csloop
   enddo phloop
!
1000 continue
   return
 end subroutine todo_before

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine todo_after_found_equilibrium(mode,ceq)
! this is called after an equilibrium calculation
! It marks stable phase (set CSSTABLE and remove any CSAUTO)
! remove redundant unstable composition sets created automatically
! (CSAUTO set).  It will also shift stable composition sets to loweest 
! possible (it will take into account if there are default constituent 
! fractions, CSDEFCON set).
! mode determine some of the actions
!
! >>>>>>>>>>> THIS IS DANGEROUS IN PARALLELL PROCESSING
! It should work in step and map as a composition set that once been stable
! will never be removed except if one does global minimization during the
! step and map. Then  metallic-FCC and MC-carbides may shift composition sets.
! Such shifts should be avoided by manual entering of comp.sets with
! default constitutions, but comparing a stable constitution with a
! default is not trivial ...
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode
!\end{verbatim}
   integer iph,ics,lokph,lokics,jcs,lokjcs,lastset,lokkcs,kzz,jtup,qq
   integer jstat2,fit,phs
   double precision val,xj1,xj2
   logical notok,noremove
   character jpre*4,jsuf*4
!
   if(btest(globaldata%status,GSNOAFTEREQ)) goto 1000
!   write(*,*)'3Y in todo_after'
! First shift all stable composition down to lower comp.sets
   phloop1: do iph=1,noph()
      lokph=phases(iph)
      if(btest(phlista(lokph)%status1,PHHID)) cycle
      csloop1: do ics=2,phlista(lokph)%noofcs
         lokics=phlista(lokph)%linktocs(ics)
         if(ceq%phase_varres(lokics)%phstate.eq.PHENTSTAB .and. &
              btest(ceq%phase_varres(lokics)%status2,CSAUTO)) then
            fit=100
! This comp.set is stable, check if a lower compset is unstable
            csloop2: do jcs=1,ics-1
               lokjcs=phlista(lokph)%linktocs(jcs)
! hidden=-4, suspended=-3, dormant=-1, unstable=-1, unknown=0, stable=1, fix=2
               if(ceq%phase_varres(lokjcs)%phstate.le.PHENTERED) then
                  if(btest(ceq%phase_varres(lokjcs)%status2,CSDEFCON)) then
! check if composition of lokics fits defaults in lokjcs
                     if(.not.checkdefcon(lokics,lokjcs,fit,ceq)) cycle csloop2
                  endif
!                  write(*,*)'3Y Moving comp.set ',ics,' down to ',jcs
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
!               write(*,*)'3Y Imperfect fit to default: ',fit,lokics,lokjcs
500            continue
! move STABLE lokics to UNSTABLE lokjcs
!                  write(*,381)'3Y Before copy',&
!                       lokics,ceq%phase_varres(lokics)%status2,&
!                       ceq%phase_varres(lokics)%phtupx,&
!                       ceq%phase_varres(lokics)%suffix,&
!                       lokjcs,ceq%phase_varres(lokjcs)%status2,&
!                       ceq%phase_varres(lokjcs)%phtupx,&
!                       ceq%phase_varres(lokjcs)%suffix
381               format(a,3i4,' "',a,'" ' ,3i4,' "',a,'"')
! list the records to switch, note default constitution??
!               write(*,380)'3Y Stable and free: ',ics,lokics,jcs,lokjcs,&
!                    ceq%phase_varres(lokics)%phtupx,&
!                    ceq%phase_varres(lokjcs)%phtupx
380               format(a,10i5)
!                  exit csloop2
! save some jcs values of amount, dgm, status, pre&suffix and tuple index
                  xj1=ceq%phase_varres(lokjcs)%amfu
                  xj2=ceq%phase_varres(lokjcs)%dgm
                  jtup=ceq%phase_varres(lokjcs)%phtupx
                  jstat2=ceq%phase_varres(lokjcs)%status2
                  jpre=ceq%phase_varres(lokjcs)%prefix
                  jsuf=ceq%phase_varres(lokjcs)%suffix
                  phs=ceq%phase_varres(lokjcs)%phstate
! copy main content of the phase_varres(lokics) record to phase_varres(lokjcs)
                  ceq%phase_varres(lokjcs)=ceq%phase_varres(lokics)
! Some content in jcs must be set or restorted separately
                  ceq%phase_varres(lokjcs)%phtupx=jtup
                  ceq%phase_varres(lokjcs)%status2=jstat2
                  ceq%phase_varres(lokjcs)%prefix=jpre
                  ceq%phase_varres(lokjcs)%suffix=jsuf
                  ceq%phase_varres(lokjcs)%phstate=PHENTSTAB
                  ceq%phase_varres(lokjcs)%status2=&
                       ibset(ceq%phase_varres(lokjcs)%status2,CSSTABLE)
! maybe CSAUTO bit set, always remove it!
!                  write(*,*)'3Y Ensure CSAUTO cleared in ',jcs
                  ceq%phase_varres(lokjcs)%status2=&
                       ibclr(ceq%phase_varres(lokjcs)%status2,CSAUTO)
! Some content in ics must be set separately from saved values of jcs
                  ceq%phase_varres(lokics)%amfu=xj1
                  ceq%phase_varres(lokics)%dgm=xj2
                  ceq%phase_varres(lokics)%phstate=phs
! clear the stable bit
                  ceq%phase_varres(lokics)%status2=&
                       ibclr(ceq%phase_varres(lokics)%status2,CSSTABLE)
! check things again
!                  write(*,381)'3Y After copy ',&
!                       lokics,ceq%phase_varres(lokics)%status2,&
!                       ceq%phase_varres(lokics)%phtupx,&
!                       ceq%phase_varres(lokics)%suffix,&
!                       lokjcs,ceq%phase_varres(lokjcs)%status2,&
!                       ceq%phase_varres(lokjcs)%phtupx,&
!                       ceq%phase_varres(lokjcs)%suffix
!                  write(*,380)'After copy',(phlista(lokph)%linktocs(qq),&
!                       qq=1,phlista(lokph)%noofcs)
                  exit csloop2
            enddo csloop2
         endif
      enddo csloop1
   enddo phloop1
! Here we may try to ensure that the stable comp.sets fits the
! default constitutions of their current set
!   write(*,*)'3Y Try to shift to match default and current constitution'
   call shiftcompsets(ceq)
!
! upto now is safe ... now remove CSAUTO comp.sets if allowed
!   write(*,*)'3Y Now maybe remove redundant compsets.'
!   goto 1000
! check if allowed to remove
   if(btest(globaldata%status,GSNOREMCS)) goto 1000
!
! Now try to remove unstable composition sets with CSAUTO set
   phloop: do iph=1,noph()
      noremove=.FALSE.
      lokph=phases(iph)
      if(btest(phlista(lokph)%status1,PHHID)) cycle
! loop backwards for compsets to remove unstable with CSAUTO set
      lastset=phlista(lokph)%noofcs
      csloopdown: do ics=lastset,2,-1
         lokics=phlista(lokph)%linktocs(ics)
!         write(*,*)'Checking comp.set ',ics
         auto: if(btest(ceq%phase_varres(lokics)%status2,CSAUTO)) then
            if(ceq%phase_varres(lokics)%phstate.le.PHENTERED) then
! comp.set was created automatically but is not stable, it can be removed
               if(noeq().eq.1) then
! we have just one equilibrium, OK to remove
!                  write(*,*)'3Y Trying to remove phase tuple ',&
!                       ceq%phase_varres(lokics)%phtupx
                  call remove_composition_set(iph,.FALSE.)
                  if(gx%bmperr.ne.0) goto 1000
               else
! if we cannot remove the comp.set remove the CSAUTO bit
                  ceq%phase_varres(lokics)%status2=&
                       ibclr(ceq%phase_varres(lokics)%status2,CSAUTO)
               endif
            else
! the comp.set is stable, remove the CSAUTO bit
!               write(*,*)'Removing CSAUTO bit',ics
               ceq%phase_varres(lokics)%status2=&
                    ibclr(ceq%phase_varres(lokics)%status2,CSAUTO)
            endif
!         else
! anything to be done with any other phase?
         endif auto
      enddo csloopdown
   enddo phloop
!
1000 continue
   return
 end subroutine todo_after_found_equilibrium

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 logical function checkdefcon(lokics,lokjcs,fit,ceq)
! check if composition of lokics fits default constitution in lokjcs
! return TRUE if lokics moved to lokjcs
! If not moved fit returns a value how close the constitition is
! If 1 very close, 2 less etc.
   integer lokics,lokjcs,fit
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer kk
   logical tobeshifted
   tobeshifted=.FALSE.
! A fraction with a maximum set (mmyfr<0) must be below that value
! A fraction with a minimum set (mmyfr>0) should be above that value
!   write(*,*)'3Y testing defaults',lokics,lokjcs
   fit=0
   do kk=1,size(ceq%phase_varres(lokjcs)%yfr)
      if(ceq%phase_varres(lokjcs)%mmyfr(kk).lt.0.0D0) then
! A fraction with a maximum set (mmyfr>0) must be below mmyfr(kk)
         if(ceq%phase_varres(lokics)%yfr(kk).gt.&
              abs(ceq%phase_varres(lokjcs)%mmyfr(kk))) fit=fit+5
! A fraction with a minimum set (mmyfr<0) should be above mmyfr(kk)
      elseif(ceq%phase_varres(lokjcs)%mmyfr(kk).gt.0.0D0) then
         if(ceq%phase_varres(lokics)%yfr(kk).lt.&
              abs(ceq%phase_varres(lokjcs)%mmyfr(kk))) fit=fit+1
      endif
! if mmyfr(kk)=0 there is no min/max for that fraction
!      write(*,77)'3Y Constitution: ',kk,ceq%phase_varres(lokjcs)%mmyfr(kk),&
!           ceq%phase_varres(lokics)%yfr(kk),fit
77    format(a,i3,2(1pe12.4),i5)
   enddo
   if(fit.eq.0) tobeshifted=.TRUE.
!   write(*,*)'3Y checkdefcon: ',lokics,lokjcs,fit
1000 continue
   checkdefcon=tobeshifted
   return
 end function checkdefcon

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine shiftcompsets(ceq)
! check phase with several composition sets if they should be shifted
! to fit the default constitution better
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,iph,ics,lokics,jcs,lokjcs,fit1,fit2,shifts
   character ch1*1
   phloop: do iph=1,noofph
      lokph=phases(iph)
      if(phlista(lokph)%noofcs.gt.1) then
         shifts=0
100      continue
         shifts=shifts+1
         if(shifts.gt.2) cycle phloop
         csloop1: do ics=1,phlista(lokph)%noofcs
            lokics=phlista(lokph)%linktocs(ics)
            if(ceq%phase_varres(lokics)%phstate.eq.PHENTSTAB) then
               if(btest(ceq%phase_varres(lokics)%status2,CSDEFCON)) then
! if TRUE then the composition fits default constitution
                  if(checkdefcon(lokics,lokics,fit1,ceq)) cycle csloop1
               else
                  fit1=-100
               endif
            else
               fit1=500
            endif
! The values of fit1:
! -100 there is no default constitution of ics    
! >0 the degree of fit of the current constitution of ics
! 500 ics is not stable
! We come here to check if some other compset fits better
!            write(*,*)'3Y ics fit in ics: ',ics,fit1
            csloop2: do jcs=1,phlista(lokph)%noofcs
               if(jcs.eq.ics) cycle csloop2
               lokjcs=phlista(lokph)%linktocs(jcs)
               if(ceq%phase_varres(lokjcs)%phstate.eq.PHENTSTAB) then
                  if(btest(ceq%phase_varres(lokjcs)%status2,CSDEFCON)) then
! if this call returns TRUE then the jcs composition fits default constitution
                     if(checkdefcon(lokjcs,lokjcs,fit2,ceq)) then
                        cycle csloop2
                     endif
                  else
! there is no default constitution in jcs
                     fit2=-100
!                     write(*,*)'3Y no default const: ',jcs,fit2
                  endif
               else
! jcs is not stable
                  fit2=500
               endif
! fit2:
! -100  jcs has no default constitution
! >0     the current fit to default constitution
! 500  jcs is not stable
!               write(*,*)'3Y jcs fit in jcs: ',jcs,fit2
               if(fit1.eq.500) then
! If neither ics nor jcs are stable increment jcs
                  if(fit2.eq.500) cycle csloop2
! ics is unstable, but if it has a default constitution, check if jcs fits
                  if(btest(ceq%phase_varres(lokics)%status2,CSDEFCON)) then
                     if(checkdefcon(lokjcs,lokics,fit2,ceq)) continue
                     fit2=-fit2
                  else
                     fit2=fit1+1
                  endif
               else
! ics is stable, check if jcs has a default constitution that fits better
                  if(btest(ceq%phase_varres(lokjcs)%status2,CSDEFCON)) then
                     if(checkdefcon(lokics,lokjcs,fit2,ceq)) continue
                  else
                     fit2=fit1+1
                  endif
               endif
! fit2:
! <=fit1  shift jcs and ics
! >fit1 do nothing
!               write(*,*)'3Y jcs fit in ics: ',jcs,fit2
               if(fit2.le.fit1) then
! The comp.set ics fits the default constitution of jcs better than its current
!                  write(*,*)'3Y shifting compsets: ',ics,jcs,fit1,fit2
                  call copycompsets2(lokph,ics,jcs,ceq)
! shift composition sets! Copy all via a dummy record (last phase_varres)
! That is hopefully unused ...
!                  ceq%phase_varres(2*maxph)=ceq%phase_varres(lokics)
!                  ceq%phase_varres(lokics)=ceq%phase_varres(lokjcs)
!                  ceq%phase_varres(lokjcs)=ceq%phase_varres(2*maxph)
! restore phtupx, pre/suffix and status word for jcs
!                  ceq%phase_varres(lokjcs)%phtupx=&
!                       ceq%phase_varres(lokics)%phtupx
!                  ceq%phase_varres(lokjcs)%status2=&
!                       ceq%phase_varres(lokics)%status2
!                  ceq%phase_varres(lokjcs)%prefix=&
!                       ceq%phase_varres(lokics)%prefix
!                  ceq%phase_varres(lokjcs)%suffix=&
!                       ceq%phase_varres(lokics)%suffix
! restore phtupx, pre/suffix and status word for ics
!                  ceq%phase_varres(lokics)%phtupx=&
!                       ceq%phase_varres(2*maxph)%phtupx
!                  ceq%phase_varres(lokjcs)%status2=&
!                       ceq%phase_varres(2*maxph)%status2
!                  ceq%phase_varres(lokics)%prefix=&
!                       ceq%phase_varres(2*maxph)%prefix
!                  ceq%phase_varres(lokics)%suffix=&
!                       ceq%phase_varres(2*maxph)%suffix
! restore the "end" link in last record
!                  ceq%phase_varres(2*maxph)%nextfree=-1
! update fit1, fit2 will be updated automatically
!                  if(ceq%phase_varres(lokics)%phstate.eq.PHENTSTAB) then
!                     if(btest(ceq%phase_varres(lokics)%status2,CSDEFCON)) then
!                        if(checkdefcon(lokics,lokics,fit1,ceq)) continue
!                     else
!                        fit1=-100
!                     endif
!                     fit1=500
!                  write(*,*)'3Y Switched compsets: ',ics,jcs,fit1,fit2
!                  read(*,99)ch1
!99                format(a)
! start again from first comp.set
                  goto 100
               endif
            enddo csloop2
         enddo csloop1
      endif
   enddo phloop
1000 continue
   return
 end subroutine shiftcompsets

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine copycompsets(iph,ics1,ics2,ceq)
! copy constitution and results from ic2 to ic1 and vice versa
   integer iph,ics1,ics2
   type(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim} %+
   integer lokph,lokcs1,lokcs2
! check indices are correct
   call get_phase_compset(iph,ics1,lokph,lokcs1)
   call get_phase_compset(iph,ics2,lokph,lokcs2)
   if(gx%bmperr.ne.0) goto 1000
   call copycompsets2(lokph,ics1,ics2,ceq)
1000 continue
   return
 end subroutine copycompsets

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim} %-
 subroutine copycompsets2(lokph,ics1,ics2,ceq)
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
! curlat, cinvy, cxmol, cdxmol?
1000 continue
   return
 end subroutine copycompsets2

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
   write(*,*)'Not implemented testing defaults'
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
   write(*,11)'bounds column 1: ',nb,(bounds(k1,1),k1=1,nb)
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

