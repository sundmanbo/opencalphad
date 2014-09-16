!
! included in pmod25.F90
!
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\
!>      17. Grid minimizer
!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

!\begin{verbatim}
 subroutine global_gridmin(what,tp,xknown,nvsph,iphl,icsl,aphl,&
      nyphl,yphl,cmu,ceq)
!
! finds a set of phases that is a global start point for an equilibrium 
! calculation at T and P values in tp and known mole fraction in xknown
! It is intentional that this routine is independent of current conditions
! returns: nvsph stable phases, list of phases in iphl, amounts in aphl, 
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
   integer, parameter :: maxgrid=100000,maxy=2000,maxph=500
   integer :: starttid,endoftime
   real finish2
   double precision amount,sum
   integer i,ibias,ics,ics2,icsno,icsx,ie,iph,iv,j1,j2,jip,jp,kkz,kp,kph,jbias
   integer lokcs,lokph,mode,ng,nocsets,noofgridpoints,nr,nrel,nrph,ny,nyz,ioff
! kphl(iph) is first gridpoint in phase iph
! ngrid(iph) is the last gridpoint for phase iph (some phases may be suspended)
! xarr(nrel,i) is the composition of gridpoint i
! garr(i) is the Gibbs energy of gridpoint i
! jgrid(j) is a gridpoint in the solution
! phfrac(j) is the amount of the phase of that gridpoint
   integer, dimension(maxph) :: ngrid,kphl
   integer, dimension(maxel) :: jgrid
   real garr(maxgrid),starting,finished,xov(maxel)
   real, dimension (:,:), allocatable :: xarr
   real, dimension (maxel,maxel) :: xsol
   double precision, dimension(maxel) :: phfrac,phsave
   double precision qq(5),savetp(2)
   integer, dimension(maxph) :: iphx
   character name1*24,ch1*1
! debug
   logical trace
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
   trace=.false.
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
   ggloop: do iph=1,nrph
!      if(.not.phase_status(iph,1,PHHID,ceq)) then
! skip phases that are hidden or suspended .... old 5; new -3
!      if(test_phase_status(iph,1,amount,ceq).lt.5) then
!      if(test_phase_status(iph,1,amount,ceq).gt.PHHIDDEN) then
      if(test_phase_status(iph,1,amount,ceq).gt.PHSUS) then
         do ics=1,noofcs(iph)
! old: 1=entered, 2=fix, 3=dormant, 4=suspended
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
!         if(trace) then
!            call get_phase_name(iph,1,name1)
!            write(*,21)iph,name1(1:12),kphl(pph),ny,ng
!         endif
!         write(*,22)iph,kphl(pph),ny,ng,pph
      endif
   enddo ggloop
! we have a grid for pph phases, note that pph is not a phase index!!!
! the phase index for phase 1..pph is in iphx(1..pph)
21 format('Gridpoints for phase ',i3,': ',a,', starts at ',i5,', with ',2i5)
22 format('Gridpoints for phase ',i3,' starts at ',i5,', with ',2i5,i8)
   if(kp-1.gt.maxgrid) then
      write(*,*)'Too many gridpoints'
      gx%bmperr=4175; goto 1000
   endif
! we may nave no gridpoints!!!
   if(kp.eq.1) then
      write(*,*)'No phases, no gridpoints'
      gx%bmperr=4176; goto 1000
   endif
!   write(*,*)'phases and gridpoints: ',pph,kp,ngrid(pph),nrel
! total number of gridpoints is kp-1 ... but sometimes kp is wrong, why??
!   allocate(xarr(nrel,kp-1))
   allocate(xarr(nrel,kp-1+10))
   if(ocv()) write(*,*)'Gridpoints and elements: ',kp-1,nrel
! generate grid
! we must know before this loop how many gridpoints that each phase will
! need.  That is a function of the number of gridpoints.
   kp=1
   call system_clock(count=starttid)
! OpenMP parallellization START
! the error code gx%bmperr should also be threadprivate
!x$omp parallel do private(ng,iv),schedule(dynamic)
!    phloop: do iph=1,nrph
   phloop: do zph=1,pph
! for phase iphx(zph) the gridpoints will be stored from position kphl(zph)
! mole fracts in xarr, g in garr
! yphl is not used when mode=0, ng should be set to number of remaining points
! ngrid(iph) is number of gridpoints in phase iph
      ng=maxgrid
!$       jip=omp_get_thread_num()
! values in kphl set in previous call to generate_grid(-1,.....)
      iv=kphl(zph)
! this call will calculate all gridpoints
      call generate_grid(0,iphx(zph),ng,nrel,xarr(1,iv),garr(iv),ny,yphl,ceq)
      if(gx%bmperr.ne.0) then
         write(*,*)'grid error ',jip,zph,gx%bmperr
         goto 1000
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
!x$omp end parallel do
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
62 format('25J Gridp: ',10i4)
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
!                  write(*,*)'25J new composition set for phase: ',j2,kph
!                  call add_composition_set(kph,'    ','AUTO',icsno,ceq)
!                  call add_composition_set(kph,'    ','AUTO',icsno,firsteq)
! It must be done in all equilibrium records, no equilibrium record needed!!!
                  call add_composition_set(kph,'    ','AUTO',icsno)
                  if(gx%bmperr.ne.0) goto 1000
                  call get_phase_compset(kph,icsno,lokph,lokcs)
                  if(gx%bmperr.ne.0) goto 1000
                  ceq%phase_varres(lokcs)%status2=&
                       ibset(ceq%phase_varres(lokcs)%status2,CSAUTO)
                  nocsets=nocsets+1
!                  write(*,303)'25J Created cs:',kph,icsno,lokcs,&
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
! this means status entered
                  ceq%phase_varres(lokcs)%status2=&
                       ibclr(ceq%phase_varres(lokcs)%status2,CSSUS)
                  ceq%phase_varres(lokcs)%status2=&
                       ibclr(ceq%phase_varres(lokcs)%status2,CSFIXDORM)
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
1010  format(' Grid minimization: ',i5,' gridpoints ',1pe12.4,' s and ',&
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
   integer nrel,lokres,ls
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

!-\begin{verbatim}
 subroutine generate_grid_v2(mode,iph,ngg,nrel,xarr,garr,ny,yarr,ceq)
! generate grid for phase iph
! Different action depending of the value of mode, 
! for mode<0:  
!    return the number of gridpoints that will be generated for phase iph
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
   integer mode,iph,ngg,nrel,ny,lokph,errsave,jn1,jn2
!
!-\end{verbatim}
   double precision, parameter :: yzero=1.0D-12
   real xarr(nrel,*),garr(*)
   integer abrakadabra,i,ibas,ibin,iend,is,iter,je,kend,ll,ls,nend,ie
   character ch1*1
   double precision yarr(*),xmol(maxel),wmass(maxel),ydum(maxconst)
   integer ikon(2,maxsubl),ngdim,nsl,l1,l2,l3,l2a,l3a
   integer nkl(maxsubl),knr(maxconst),inkl(0:maxsubl),nofy
   double precision, dimension(:), allocatable :: yfra
   double precision sites(maxsubl),qq(5),gval
! endm(i,j) has constituent indices in i=1..nsl for endmember j 
   integer, dimension(:,:), allocatable :: endm
   integer ifra(3),jend(3)
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
! with 3 endmembers: 3*11=33 gridpoints (binary combinations)
! (1.00,0.00,0.00) 1st--------
! (0.89,0.11,0.00)(0.89,0.00,0.11)
! (0.74,0.26,0.00)(0.74,0.00.0.26) binary
!                                  (0.74,0.15,0.11)(0.74,0.11,0.15) ternay
! (0.61,0.39,0.00)(0.61,0.00,0.39) binary
!                                  (0.61,0.26,0.13)(0.61,0.13,0.26) ternary
! (0.00,1.00,0.00) 2nd---------
! (0.11,0.89,0.00)(0.00,0.89.0.11)
! ...
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
!   double precision, dimension(4), parameter:: ybas=&
!        [1.0D0,0.89D0,0.74D0,0.61D0]
!   double precision, dimension(4), parameter :: ybin=&
!        [0.11D0,0.26D0,0.39D0,0.15D0]
!   double precision, dimension(3), parameter :: yter=[0.0D0,0.11D0,0.13D0]
!
!   write(*,*)'entering generate_grid: ',mode,iph,ngg
   call get_phase_data(iph,1,nsl,nkl,knr,ydum,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! calculate the number of endmembers and index of first constituent in subl ll
   nend=1
   inkl(0)=0
   do ll=1,nsl
      nend=nend*nkl(ll)
      inkl(ll)=inkl(ll-1)+nkl(ll)
   enddo
   ny=inkl(nsl)
   lokph=phases(iph)
   negmode: if(mode.lt.0) then
! Any changes here must be made also for mode=0
!---------------------------------------------------------
! just determine the number of gridpoints for this phase for global minimimum
! ideal gases should just have the endmembers ....
      ngdim=ngg
      ngg=nend
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
         ngg=ngg+nend*(1+3*(nend-1))
!         write(*,*)' 10<nend<15 ',iph,nend,ngg
      else
! 3-10: three binary and two ternary combinarions (all)
         ngg=nend*(1+(nend-1)*(3+2*(nend-2))) ! (ternary combinations skipped)
         ngg=ngg+nend*(1+3*(nend-1))
!         if(ocv()) write(*,*)' 3<nend<11 ',iph,nend,ngg
      endif
      ny=nend
      goto 1000
   endif negmode
!------------------------------------------------------------
! for mode=0:
!    generate all gridpoint sitefractions and calculate G and molefracs
!    Any changes here must be made also for mode=-1
!------------------------------------------------------------
! for mode>0:
!   return sitefractions (for mode=gridpoint number (part of the solution))
!   BUT: The only way to find the site fraction of a gripoint is to generate
!   all gridpoints up the one specified by the value of mode (no G calculation)
!------------------------------------------------------------
! this subroutine is not used
   allocate(endm(nsl,nend))
   allocate(yfra(inkl(nsl)))
   nofy=inkl(nsl)
! generate endmembers, endm(ll,ie) is set to consituent index in sublattice ll
! note: inkl(0)=0, this has to be identical to what is returned for mode=-1
   je=1
   do ll=1,nsl
      endm(ll,je)=inkl(ll-1)+1
   enddo
   do while(je.lt.nend)
      je=je+1
      do ls=1,nsl
         endm(ls,je)=endm(ls,je-1)
      enddo
      ll=0
110   ll=ll+1
      if(endm(ll,je).lt.inkl(ll)) then
         endm(ll,je)=endm(ll,je)+1
      elseif(ll.lt.nsl) then
         endm(ll,je)=inkl(ll-1)+1
         goto 110
      else
         gx%bmperr=4148; goto 1000
      endif
   enddo
!=========================================
! now generate all combinations of endmembers
   ngg=0
   jend=0
!-------------------------------------------------
! loop to calculate gridpoints
   endc: do l1=1,nend
      yfra=yzero
      do ls=1,nsl
         yfra(endm(ls,l1))=one
      enddo
      jend(1)=jend(1)+1
      jend(2)=0
      jend(3)=0
      ifra=0
      ifra(1)=1
! calculate gridpoint for a pure endmember vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      call new_gridpoint_calc(mode,iph,nsl,nend,endm,jend,ifra,ny,yfra,&
           xmol,gval,ceq)
      if(gx%bmperr.ne.0) goto 1000
      ngg=ngg+1
      if(mode.eq.ngg) goto 500
      do ie=1,nrel
         xarr(ie,ngg)=real(xmol(ie))
      enddo
      garr(ngg)=real(gval)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Any changes here has also to be made for mode<0 above !!!!
! only endmembers (1:0:0), nend points
      if(nend.eq.1 .or. nend.gt.50 .or. &
           btest(phlista(lokph)%status1,PHID)) cycle endc
!-------------------------------------------------
! loop to calculate binary and ternary gridpoints
      ifra(1)=2
      ifra(2)=1
      jend(2)=0
! we may jump back here for another set of binary combinations
200   continue
      binc: do l2=1,nend-1
         jend(2)=jend(2)+1
         if(jend(2).eq.jend(1)) jend(2)=jend(2)+1
! calculate gridpoint for a binary combination of endmembers vvvvvvvvvvv
         call new_gridpoint_calc(mode,iph,nsl,nend,endm,jend,ifra,ny,yfra,&
              xmol,gval,ceq)
         if(gx%bmperr.ne.0) goto 1000
         ngg=ngg+1
         if(mode.eq.ngg) goto 500
         do ie=1,nrel
            xarr(ie,ngg)=real(xmol(ie))
         enddo
         garr(ngg)=real(gval)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         if(nend.gt.30 .or. &
              nend.gt.25 .or. &
              nend.eq.2 .or. nend.gt.14) then ! nend=2 or nend>14) then
! Any changes here has also to be made for mode<0 above !!!!
! one binary combination (2:1:0), min 961, max 2500 points
! two binary combinations (2:1:0) and (3:2:0), min 1326, max 1770 points
! three binary combinations (2:1:0), (3:2:0) and (4:3:0) max 1825 points
            cycle binc
         endif
         if(nend.gt.10) then !11<nend<15
! Any changes here has also to be made for mode<0 above !!!!
! three binary combinations (2:1:0), (3:2:0) and (4:3:0)
! and one ternary combination (3:4:1)
            if(ifra(1).ne.3) cycle binc
            ifra(2)=4
            ifra(3)=1
         else ! the final case 2<nend<11
! Any changes here has also to be made for mode<0 above !!!!
! three binary combinations (2:1:0), (3:2:0) and (4:3:0)
! and two ternary combination (3:4:1) and (4:2:2)
            if(ifra(1).eq.3) then
               ifra(2)=4
               ifra(3)=1
            elseif(ifra(1).eq.4) then
               ifra(2)=2
               ifra(3)=2
            else
               cycle binc
            endif
         endif
!-------------------------------------------------
! loop to calculate ternary gridpoints
         jend(3)=0
         terc: do l3=1,nend-2
            jend(3)=jend(3)+1
            if(jend(3).gt.nend) jend(3)=1
            if(jend(3).eq.jend(1)) jend(3)=jend(3)+1
            if(jend(3).eq.jend(2)) jend(3)=jend(3)+1
            if(jend(3).eq.jend(1)) jend(3)=jend(3)+1
! calculate gridpoint for a ternary combination of endmember vvvvvvvvvvvvvvvv
            call new_gridpoint_calc(mode,iph,nsl,nend,endm,jend,ifra,ny,yfra,&
                 xmol,gval,ceq)
            if(gx%bmperr.ne.0) goto 1000
            ngg=ngg+1
            if(mode.eq.ngg) goto 500
! save values as real to save space, arrays can be very big!
            do ie=1,nrel
               xarr(ie,ngg)=real(xmol(ie))
            enddo
            garr(ngg)=real(gval)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            write(*,351)'at terc:    ',nend,(jend(ls),ls=1,3),(ifra(ls),ls=1,3)
         enddo terc
!-------------------------------------------------------
         write(*,351)'after terc: ',nend,(jend(ls),ls=1,3),(ifra(ls),ls=1,3)
! we are here when we finished a loop of ternary combinations, 3<nend<15
         if(nend.gt.10) then
! three binary combinations (2:1:0), (3:2:0) and (4:3:0)
! and one ternary combination (3:4:1) (now done)
            ifra(1)=4
            ifra(2)=3
            ifra(3)=0
            jend(3)=0
         else ! 2<nend<11
! three binary combinations (2:1:0), (3:2:0) and (4:3:0)
! and two ternary combination (3:4:1) and (4:2:2)
            if(l2.lt.nend-1) then
! binary loop in jend(2) before next ternary
               if(ifra(1).eq.3) then
                  ifra(2)=2
                  ifra(3)=0
                  jend(3)=0
               elseif(ifra(1).eq.4) then
                  ifra(2)=3
                  ifra(3)=0
                  jend(3)=0
               endif
! this is hopefully taken care of after finishing the binary loop
!            else
!               if(ifra(1).eq.3) then
!                  ifra(2)=4
!                  ifra(3)=1
!                  jend(3)=0
!               elseif(ifra(1).eq.4) then
!                  ifra(2)=2
!                  ifra(3)=2
!                  jend(3)=0
!               endif
            endif
         endif
         write(*,351)'at binc:    ',nend,(jend(ls),ls=1,3),(ifra(ls),ls=1,3)
      enddo binc
!---------------------------------------------------------
      write(*,351)'after binc: ',nend,(jend(ls),ls=1,3),(ifra(ls),ls=1,3)
! we have finished a binary combination and possibly ternaries
! but we may have to make several binary loops
      case1: if(nend.gt.30) then
! one binary combination (2:1:0), min 961, max 2500 points
         if(ifra(1).eq.1) then
            ifra(1)=2
            ifra(2)=1
            jend(2)=0
            goto 200
         endif
      elseif(nend.gt.25) then
! two binary combinations (2:1:0) and (3:2:0), min 1326, max 1770 points
         if(ifra(1).lt.3) then
            ifra(1)=ifra(1)+1
            ifra(2)=ifra(2)+1
            jend(2)=0
            goto 200
         endif
      else ! for any other values of nend=2..25
! three binary combinations (2:1:0), (3:2:0) and (4:3:0) max 1825 points
! ALSO: one ternary combination (3:4:1)
! ALSO: and two ternary combination (3:4:1) and (4:2:2)
         if(ifra(1).lt.4) then
            ifra(1)=ifra(1)+1
            ifra(2)=ifra(2)+1
            if(ifra(2).gt.4) ifra(2)=2
            jend(2)=0
            goto 200
         endif
      endif case1
      write(*,351)'at endc:    ',nend,(jend(ls),ls=1,3),(ifra(ls),ls=1,3)
351   format(a,i3,5x,3i3,2x,3i3)
   enddo endc
! finished all calculations
   ny=nend
   goto 1000
! all binary and ternary combination of endmember for the grid above
!
! This should be modeified to take into account option F and B as those
! sublattices are identical ... have fun doing that ...
!
!----------------------------------------
! here we return the constitution for gridpoint "mode" in the solution
! We must also return the mole fractions??
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
 end subroutine generate_grid_v2

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine generate_grid(mode,iph,ngg,nrel,xarr,garr,ny,yarr,ceq)
! Different action depending of the value of mode, 
! for mode<0:  
!    return the number of gridpoints that will be generated for phase iph
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
   character ch1*1
   double precision xmol(maxel),wmass(maxel),ydum(maxconst)
   integer ikon(2,maxsubl),ngdim,nsl
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
   integer lut,jbas,sumngg
   logical trace,isendmem
   save sumngg
!
!   write(*,*)'entering generate_grid: ',mode,iph,ngg
   if(mode.eq.0) then
!      write(*,*)'Generating grid for phase: ',iph
! trace TRUE means generate outpy for each gridpoint
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
   ny=inkl(nsl)
!   write(*,1010)'Saved   ',iph,(ydum(i),i=1,ny)
   negmode: if(mode.lt.0) then
!---------------------------------------------------------
! just determine the number of gridpoints for this phase for global minimimum
! ideal gases should just have the endmembers ....
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
! sublattices are identical ... have fun doing that ...
!
!----------------------------------------
! here we return the constitution for gridpoint "mode" in the solution
! We must also return the mole fractions??
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
   integer i,lokres,ny
   double precision qq(5),xmol(nrel),wmass(nrel)
   character ch1*1
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
! the gridpoint has net charge, qq(2), make gval more positive. 
! Note gval(1,1) is divided by RT so around -5<0
      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+qq(2)**2)
!      gval=real(ceq%phase_varres(lokres)%gval(1,1)/qq(1)+5*qq(2)**2)
      if(ocv()) write(*,66)'25J charged gp: ',&
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
 subroutine calcg_endmember(iph,endmember,gval,ceq)
! calculates G for one mole of real atoms for a single end member
! used for reference states. Restores current composition (but not G or deriv)
! endmember contains indices in the constituent array, not species index
   implicit none
   integer iph
   double precision gval
   integer endmember(maxsubl)
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
   integer ierr,kk0,ll,lokres,nsl
   integer nkl(maxsubl),knr(maxconst)
   double precision savey(maxconst),sites(maxsubl),qq(5),yfra(maxconst)
!
   call get_phase_data(iph,1,nsl,nkl,knr,savey,sites,qq,ceq)
   if(gx%bmperr.ne.0) goto 1100
! set constitution to be just the endmember
   yfra=zero
   kk0=0
   do ll=1,nsl
      if(endmember(ll).gt.kk0 .and. endmember(ll).le.nkl(ll)) then
         yfra(endmember(ll))=one
      else
!         write(*,*)'endmember index outside range',ll,endmember(ll),nkl(ll)
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
! avoid calculating endmembers with too many vacancies. gval is divideb by RT
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

!-\begin{verbatim}
 subroutine calcg_endmember2(lokph,endmember,tpref,gval,ceq)
! calculates G for one mole of real atoms for a single end member
! used for reference states. Restores current composition (but not G or deriv)
! endmember contains indices in the constituent array, not species index
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
!         write(*,*)'endmember index outside range',ll,endmember(ll),nkl(ll)
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
   real xmat(nrel,nrel),xmatsave(nrel,nrel),xov(nrel)
! used to solve the linear system of equations
   double precision qmat(nrel,nrel+1),qmatsave(nrel,nrel+1)
   double precision zmat(nrel,nrel+1),cmusave(nrel)
   integer notuse(kp),i,ie,iel,ierr,iesave,inerr,inuse,ip,je,jj,jp,jsave,nj
   integer nrel1,nyp,iii,j,griter
   double precision xtx(nrel),phfsave(nrel)
   integer, dimension(jerr) :: removed
   real gmin(nrel),dg,dgmin,gplan,gvv,gy,gvvp
! gridpoints that has less difference with the plane than this limit is ignored
   real, parameter :: dgminlim=1.0D-6
   character ch1*1
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
   do ip=1,kp
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
         endif
120      continue
      enddo
   enddo
! check that we have nrel gridpoints
   do je=1,nrel
      if(jgrid(je).eq.0) then
! no gridpoint assigned to this element!! error (note C in pure fcc has no gp)
         gx%bmperr=4149; goto 1000
      endif
      ip=jgrid(je)
      do ie=1,nrel
         xmat(ie,je)=xarr(ie,ip)
      enddo
      gmin(je)=garr(ip)
      phfrac(je)=xknown(je)
   enddo
! looking for tbase calculation error
!   if(trace) write(*,770)(jgrid(je),je=1,nrel)
!770 format('Initial set of gridpoints: '(/15i5))
   do je=1,nrel
      if(one-xmat(je,je).lt.1.0d-12) then
         cmu(je)=dble(gmin(je))
      else
! we should have a composition for almost pure element
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
715   format(/'Initial matrix:',i3)
      do je=1,nrel
         write(31,720),xknown(je),xknown(je),(xmat(ie,je),ie=1,nrel)
      enddo
720   format(2F7.4,1x,8f8.5)
      write(31,730)gvvp,(cmu(je),je=1,nrel)
730   format('Gibbs energy: ',1pe14.6/'Chemical potentials: '/6(1pe12.4))
   endif
   griter=0
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
!               write(*,*)'25J Lower G: ',griter,nyp,kp,dgmin
            endif
! debugging LC_CsI (61) and SC_CsI (94)
!            if(jp.eq.61 .or. jp.eq.94) &
!                 write(*,44)'25J extra: ',jp,dg,dgmin,garr(jp),gplan
44          format(a,i5,5(1pe12.4))
         endif
!      else
!         write(*,*)'Excluded: ',griter,jp
      endif included
   enddo pointloop
!   write(*,*)'25J Finished loop for all gridpoints: ',jp,kp
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
      if(trace) write(*,*)'25J new gridpoint: ',griter,nyp,dgmin
   endif
!   inuse=inuse-1
   notuse(nyp)=1
!   write(*,211)'ny:',nyp,dgmin,(xarr(ie,nyp),ie=1,nrel)
   if(trace) write(*,212)'Found gridpoint ',nyp,inuse,dgmin,garr(nyp)
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
      write(*,301)'Failed to include gridpoint in solution: ',nyp,inerr,&
           dgmin,(xarr(i,nyp),i=1,nrel)
301   format(a,i5,i3/1x,1pe12.4,2x,6(1pe12.4))
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
!   write(*,*)'25J Putting gridpoint back: ',jgrid(ie)
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
199 format('25J Gibbs energy changed: ',i5,2(1pe15.6))
   gvvp=gy
!
   if(trace) then
      write(31,740)griter,nyp
740   format(/'Iteration ',i6,' found gridpoint: ',i6,', new matrix:')
      do je=1,nrel
         write(31,720),phfrac(je),xknown(je),(xmat(je,ie),ie=1,nrel)
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
   double precision xmol(nv),wmass(nv),summu,sumam
   character name*24,ch1*1
   logical igen
   real xmix(maxel)
   double precision a1,a2,gdf,gg1,gg2,gg3,gval1,gval2,gval3,gmindif
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
!               write(*,1117)'25J: ',jp,(yphl(incy(jp)+i),i=0,nyphl(jp)-1)
!               write(*,1117)'25J: ',kp,(yphl(incy(kp)+i),i=0,nyphl(kp)-1)
               do i=0,nyphl(jp)-1
                  yphl(incy(jp)+i)=a1*yphl(incy(jp)+i)+a2*yphl(incy(kp)+i)
               enddo
!               write(*,1117)'25J: ',jp,(yphl(incy(jp)+i),i=0,nyphl(jp)-1)
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
!            call add_composition_set(iph,'    ','AUTO',ics,ceq)
            call add_composition_set(iph,'    ','AUTO',ics)
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
!>      18. Miscellaneous
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
!      write(*,*)'25J allowenter ',mode,noofel,eqfree,noofph
      if(noofel.eq.0) goto 1000
      if(eqfree.gt.2) goto 1000
      yesorno=.TRUE.
   elseif(mode.eq.3) then
! there must be at lease one phase before entering a second equilibrium
! Note this is tested also for entering the default equilibrium
!      write(*,*)'25J mode 3: ',eqfree,noofph
      if(eqfree.ge.2 .and. noofph.eq.0) goto 1000
      yesorno=.TRUE.
   endif
1000 continue
   allowenter=yesorno
!   write(*,*)'25J: allowenter:',yesorno,mode
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
!   write(*,*)'25J entering proper_symbol_name: ',name,typ
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
 subroutine sumofphcs(npx,ceq)
! returns the total number of unhidden and unsuspended phases+composition sets
! in the system.  Used for dimensioning work arrays and in loops
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer npx
!\end{verbatim}
   integer tphic,iph,ics,lokph
   double precision xxx
   tphic=0
   do iph=1,noofph
      lokph=phases(iph)
      ics=1
!      if(test_phase_status(iph,ics,xxx,ceq).eq.5) goto 500
      if(test_phase_status(iph,ics,xxx,ceq).eq.PHHIDDEN) goto 500
! phase is not hidden
      do ics=1,phlista(lokph)%noofcs
!         if(test_phase_status(iph,ics,xxx,ceq).eq.4) goto 400
         if(test_phase_status(iph,ics,xxx,ceq).eq.PHSUS) goto 400
! composition set not suspended
         tphic=tphic+1
400      continue
      enddo
500   continue
   enddo
1000 continue
   npx=tphic
   return
 end subroutine sumofphcs
 
!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine enter_default_constitution(iph,ics,mmyfr,ceq)
! set values of default constitution
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
909 format('25J 2699: ',10i4)
!   write(*,46)'25J y: ',(ceq%phase_varres(lokcs)%yfr(jl),jl=1,jk)
46 format(a,10(F7.3))
   do jl=1,phlista(lokph)%tnooffr
      ceq%phase_varres(lokcs)%mmyfr(jl)=mmyfr(jl)
!      write(*,47)'25J jl: ',jl,mmyfr(jl),&
!           firsteq%phase_varres(lokcs)%mmyfr(jl),&
!           ceq%phase_varres(lokcs)%mmyfr(jl)
   enddo
47 format(a,i2,10F7.3)
! set bit indicating that this composition set has a default constitution
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
 subroutine set_default_constitution(jph,ics,all,ceq)
! set the current constitution of jph to its default constitution
! jph can be -1 meaning all phases, all composition sets
! if all=-1 then change constitution of all phases, else just those not stable
! do not change the amounts of the phases
   implicit none
   integer all,jph,ics
   TYPE(gtp_equilibrium_data), pointer :: ceq
!\end{verbatim}
! This has been changed so it calls set_constitution !!!!
   integer iph,lokph,lokcs,ky,kz,ll,n1,n2,n3,jl
   double precision kvot1,kvot2,amount,rest,qq(5)
   double precision, dimension(:), allocatable :: yy
!
   if(jph.lt.0) then
      iph=1; ics=1
   else
      iph=jph
   endif
   call get_phase_compset(iph,ics,lokph,lokcs)
   if(gx%bmperr.ne.0) goto 1000
100 continue
   if(test_phase_status(iph,ics,amount,ceq).gt.3) goto 700
! do not change the constitution of stable phases ??
!   if(ceq%phase_varres(lokcs)%amount(1).gt.zero .and. all.ge.0) goto 700
   if(ceq%phase_varres(lokcs)%amfu.gt.zero .and. all.ge.0) goto 700
! mmyfr defines min or max default values of each constituent
! if negative it is a min value, positive is a max value, zero means no default
! It is also used to select the composition set that should be used
! when a new composition set is needed during a calculation, for example
! if an FCC phase that could be an austenite (low carbon content) or a 
! cubic carbo-nitride (high carbon or nitrogen content)
   allocate(yy(phlista(lokph)%tnooffr))
   ky=0
   subl: do ll=1,phlista(lokph)%noofsubl
      kz=ky
      n1=0
      n2=0
      n3=0
      rest=zero
      do jl=1,phlista(lokph)%nooffr(ll)
         ky=ky+1
         yy(ky)=zero
!         ceq%phase_varres(lokcs)%yfr(ky)=zero
         if(ceq%phase_varres(lokcs)%mmyfr(ky).lt.zero) then
! if mmyfr(kk) is negative the value is a maximal value (normal -1.0D-3)
! Set fraction 1/10 of this
            yy(ky)=0.1D0*abs(ceq%phase_varres(lokcs)%mmyfr(ky))
!            ceq%phase_varres(lokcs)%yfr(ky)=0.1D0*&
!                 abs(ceq%phase_varres(lokcs)%mmyfr(ky))
            n1=n1+1
         elseif(ceq%phase_varres(lokcs)%mmyfr(ky).gt.zero) then
! if mmyfr(kk) is positive the value is a minimal value (normal 0.5)
! Note that several constituents can have a minimal value and the total
! of these can be larger than unity
!            ceq%phase_varres(lokcs)%yfr(ky)=one
            yy(ky)=one
            n2=n2+1
         else
!            ceq%phase_varres(lokcs)%yfr(ky)=one
            yy(ky)=one
            n3=n3+1
         endif
      enddo
!      write(*,117)'yt: ',ky-kz,(ceq%phase_varres(lokcs)%yfr(j),j=kz+1,ky)
117   format(a,i2,9(F8.4))
! for normallizing.  The idea is that sum of fractions with min should be 0.9
! and sum of fractions with max should be summin and constituents with
! no default should be 1-0.9*summax-summin
      kvot1=one
      if(n1.gt.0) then
         kvot1=one/dble(n1)
      endif
      kvot2=one
      rest=one
      if(n2.gt.0) then
         if(n3.gt.0) then
            kvot2=0.9D0/dble(n2)
            rest=0.1D0/dble(n3)
         else
            kvot2=one/dble(n2)
         endif
      elseif(n3.gt.0) then
         rest=one/dble(n3)
      endif
!      write(*,17)'sums: ',ky-kz,kvot1,kvot2,rest
17    format(a,i3,6(1pe12.4))
! It is not necessary that the sum of fractions is unity, it will be
! normallized before used in a calculation.
      do jl=1,phlista(lokph)%nooffr(ll)
         kz=kz+1
         if(ceq%phase_varres(lokcs)%mmyfr(kz).lt.zero) then
!            ceq%phase_varres(lokcs)%yfr(kz)=kvot1*&
!                 ceq%phase_varres(lokcs)%yfr(kz)
            yy(kz)=kvot1*yy(kz)
         elseif(ceq%phase_varres(lokcs)%mmyfr(kz).gt.zero) then
!            ceq%phase_varres(lokcs)%yfr(kz)=kvot2*&
!                 ceq%phase_varres(lokcs)%yfr(kz)
            yy(kz)=kvot2*yy(kz)
         else
            yy(kz)=rest
         endif
      enddo
   enddo subl
!   write(*,117)'mm: ',kz,(ceq%phase_varres(lokcs)%mmyfr(j),j=1,kz)
!   write(*,117)'yd: ',kz,yy(j),j=1,kz)
   call set_constitution(iph,ics,yy,qq,ceq)
   if(gx%bmperr.ne.0) goto 1000
! jump here if phase skipped
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
   if(allocated(yy)) deallocate(yy)
   return
 end subroutine set_default_constitution

!/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\

!\begin{verbatim}
 subroutine todo_before(mode,ceq)
! this could be called before an equilibrium calculation
! It removes any phase amounts and clears CSSTABLE
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
! default has not been implemented yet.
!
   implicit none
   TYPE(gtp_equilibrium_data), pointer :: ceq
   integer mode
!\end{verbatim}
   integer iph,ics,lokph,lokics,jcs,lokjcs,lastset,kcs,lokkcs,kzz
   double precision val
   logical notok,noremove
!
   write(*,*)'in todo_after ... is unfinished'
   goto 1000
!
! First try to move stable composition sets to the lowesr set and delete
! redundant (with CSAUTO set)
   phloop: do iph=1,noph()
      noremove=.FALSE.
      lokph=phases(iph)
      if(btest(phlista(lokph)%status1,PHHID)) cycle
! loop backwards for compsets to remove unstable with CSAUTO set
! also clear the CSSTABLE bit in all unstable compsets
      lastset=phlista(lokph)%noofcs
      csloopdown: do ics=lastset,1,-1
! kzz=1 entered, =2 fix, =3 dormant, =4 fixed, =5 hidden
         kzz=test_phase_status(iph,ics,val,ceq)
         lokics=phlista(lokph)%linktocs(ics)
!         stable: if(ceq%phase_varres(lokics)%amount(1).eq.zero) then
         stable: if(ceq%phase_varres(lokics)%amfu.eq.zero) then
! this compset is not stable, try to remove it. compset 1 has never CSAUTO set
            zeroam: if(btest(ceq%phase_varres(lokics)%status2,CSAUTO)) then
               if(noeq().eq.1) then
                  call remove_composition_set(iph,.FALSE.)
                  if(gx%bmperr.ne.0) goto 1000
               else
! we cannot remove any (more) composition sets but maybe move stable sets
! and set/clear the CSSTABLE bit
                  noremove=.TRUE.
                  if(kzz.eq.2) then
! A fix phase can have amount zero.  It should have dgm=0
                     if(ceq%phase_varres(lokics)%dgm.ne.zero) then
                        write(*,*)'Warning, fix phase has nonzero DGM',iph,ics
                     else
                        ceq%phase_varres(lokics)%status2=&
                             ibset(ceq%phase_varres(lokkcs)%status2,CSSTABLE)
                     endif
                  else
! An unstable phase with CSAUTO not set, just clear the CSSTABLE bit
                     ceq%phase_varres(lokics)%status2=&
                          ibclr(ceq%phase_varres(lokkcs)%status2,CSSTABLE)
                  endif
               endif
            else
! we cannot remove any (more) composition sets but maybe move stable sets
! and set/clear the CSSTABLE bit
               noremove=.TRUE.
            endif zeroam
         else
! this compset is stable, if not entered just set the CSSTABLE bit
            if(kzz.ne.1) then
               ceq%phase_varres(lokics)%status2=&
                    ibset(ceq%phase_varres(lokkcs)%status2,CSSTABLE)
               cycle csloopdown
            endif
            if(btest(ceq%phase_varres(lokics)%status2,CSAUTO)) then
! check if this phase can be moved to a lower unstable compset
               jcs=0
               csloopup: do while(jcs.lt.ics-1)
                  jcs=jcs+1
                  lokjcs=phlista(lokph)%linktocs(jcs)
! this composition set must also be entered and not stable
                  if(test_phase_status(iph,jcs,val,ceq).ne.1 .or. &
                       ceq%phase_varres(lokjcs)%amfu.gt.zero) cycle
!                       ceq%phase_varres(lokjcs)%amount(1).gt.zero) cycle
! if CSDEFCON set we must check if constuent fractions fit
                  if(btest(ceq%phase_varres(lokjcs)%status2,CSDEFCON)) then
                     call compare_compsets(lokics,lokjcs,notok,ceq)
                     if(gx%bmperr.ne.0) goto 1000
                     if(notok) cycle csloopup
                  endif
!
! getting messy here, maybe take a break ...
!
! the compset jcs can be set stable instead of ics. Just change the link
! in phlista(lokph)%linktocs(jcs) from lokjcs to lokics
                  phlista(lokph)%linktocs(jcs)=lokics
                  ceq%phase_varres(lokics)%status2=&
                       ibset(ceq%phase_varres(lokkcs)%status2,CSSTABLE)
                  ceq%phase_varres(lokics)%status2=&
                       ibclr(ceq%phase_varres(lokkcs)%status2,CSAUTO)
! If we have just equilibrium we can remove the last composition set
                  if(noeq().eq.1) then
                     call remove_composition_set(iph,.FALSE.)
                     if(gx%bmperr.ne.0) goto 1000
                  endif
               enddo csloopup
! phase lokph+lokkcs is stable, set CSSTABLE and remove any CSAUTO


            else
! else just remove the CSSTABLE bit wheather set or not
               ceq%phase_varres(lokics)%status2=&
                    ibclr(ceq%phase_varres(lokics)%status2,CSSTABLE)
            endif
         endif stable
      enddo csloopdown
   enddo phloop
!
1000 continue
   return
 end subroutine todo_after_found_equilibrium

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

