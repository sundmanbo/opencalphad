!
! TQ test program 3 of programming interface LIBOCTQ using O-U binary
!
program octq1
!
  use liboctq
!
  implicit none
! maxel and maxph defined in pmod package
!  integer, parameter :: maxel=10,maxph=20
  integer n,n1,n2,n3,n4,ip,cnum(maxel+3),mm,m2
  character filename*60,phnames(maxph)*24
  character condition*60,line*80,statevar*60,quest*60,ch1*1
  character target*60,phcsname*36
  double precision value,temp,tp(2),mel(maxel),rest
  double precision xf(maxel),pxf(10*maxph),npf(maxph),mu(maxel)
  type(gtp_equilibrium_data), pointer :: ceq
!
! initiate, value of n irrelevant
  n=0
  call tqini(n,ceq)
  if(gx%bmperr.ne.0) goto 1000
!
! read database file
  filename='OU '
  call tqrfil(filename,ceq)
  if(gx%bmperr.ne.0) goto 1000
! The call to tqrfil also enters the number of elements and their names
! in the module variables nel and cnam
! and the number of phase tuples and their phase and compset indices in
! the module variables ntup and phcs
!
! Redundant: find out about the system:
! number of elements and their names
!  call tqgcom(nel,cnam,ceq)
!  if(gx%bmperr.ne.0) goto 1000
! Redundant: number of phases
!  call tqgnp(ntup,ceq)
!  if(gx%bmperr.ne.0) goto 1000
  do n=1,ntup
     call tqgpn(n,phnames(n),ceq)
     if(gx%bmperr.ne.0) goto 1000
  enddo
!
! list elements and phases
  write(*,10)nel,(cnam(n)(1:2),n=1,nel)
10 format(/'System with ',i2,' elements: ',10(a,', '))
  write(*,20)ntup,(phnames(n)(1:len_trim(phnames(n))),n=1,ntup)
20 format('and ',i3,' phases: ',5(a,', '))
!
! set values of temperature and pressure
  tp(1)=3.0D3
  tp(2)=1.0D5
  rest=one
  do n=1,nel-1
     xf(n)=0.4D0/dble(nel)
     rest=rest-xf(n)
  enddo
  xf(nel)=rest
!
! ask for conditions using the command line user interface (CUI)
100 continue
  write(*,105)
105 format(/'Give conditions:')
  ip=len(line)
  temp=tp(1)
  call gparrd('Temperature: ',line,ip,tp(1),temp,nohelp)
  if(buperr.ne.0) goto 1000
  if(tp(1).lt.1.0d0) then
     write(*,*)'Temperature must be larger than 1 K'
     tp(1)=1.0D0
  endif
  temp=tp(2)
  call gparrd('Pressure: ',line,ip,tp(2),temp,nohelp)
  if(buperr.ne.0) goto 1000
  if(tp(2).lt.1.0d0) then
     write(*,*)'Pressure must be larger than 1 Pa'
     tp(2)=1.0D0
  endif
  do n=1,nel-1
     quest='Mole fraction of '//cnam(n)(1:len_trim(cnam(n)))//':'
     temp=xf(n)
     call gparrd(quest,line,ip,xf(n),temp,nohelp)
     if(buperr.ne.0) goto 1000
     if(xf(n).lt.1.0d-6) then
        write(*,*)'Fraction must be larger than 1.0D-6'
        xf(n)=1.0D-6
     elseif(xf(n).ge.1.0d0) then
        xf(n)=0.999999D0
     endif
  enddo
! -------------------------------------
! set conditions
  n1=0
  n2=0
  condition='T'
  call tqsetc(condition,n1,n2,tp(1),cnum(1),ceq)
  if(gx%bmperr.ne.0) goto 1000
  condition='P'
  call tqsetc(condition,n1,n2,tp(2),cnum(2),ceq)
  if(gx%bmperr.ne.0) goto 1000
  condition='N'
  call tqsetc(condition,n1,n2,one,cnum(3),ceq)
  if(gx%bmperr.ne.0) goto 1000
  do n=1,nel-1
     condition='X'
     call tqsetc(condition,n,n2,xf(n),cnum(3+n),ceq)
     if(gx%bmperr.ne.0) goto 1000
  enddo
!
! Phase status can be:
! PHSUS suspended
! PHDORM dormant
! PHENTUNST entered but not stable (set by equilibrium calculation)
! PHENTERED entered (before any equilibrium calculation)
! PHENTSTAB entered and stable
! PHFIXED prescribed to be stable
! suspend all phases
  value=0.0D0
  call tqphtupsts(-1,PHSUS,value,ceq)
  if(gx%bmperr.ne.0) then
     write(*,*)'Error setting phases suspended'
     goto 1000
  endif
! restore phase 2 (1 is gas, 2 is ionic liquid)
  call tqphtupsts(2,PHENTERED,value,ceq)
  if(gx%bmperr.ne.0) then
     write(*,*)'Error setting phase 2 entered'
     goto 1000
  endif
!----------------------------------------------------------
! calculate the equilibria, n1=0 means use gridminimizer
  target=' '
  n1=0
  n2=0
  call tqce(target,n1,n2,value,ceq)
  if(gx%bmperr.ne.0) then
     write(*,310)gx%bmperr,bmperrmess(gx%bmperr)
310  format('Calculation failed, error code: ',i5/a)
     gx%bmperr=0; goto 600
  else
     write(*,320)
320  format(/'Successful calculation')
  endif
!
! list some results
! amount of all phases
  statevar='NP'
  n1=-1
  n2=0
  n3=size(npf)
  call tqgetv(statevar,n1,n2,n3,npf,ceq)
  if(gx%bmperr.ne.0) goto 1000
  write(*,*)'tuples: ',ntup,n3
  write(*,505)n3,(npf(n),n=1,n3)
505 format(/'Amount of ',i2,' phases: ',(5F7.4))
! The number of values retured is the number of phase tuples
! It has been automatically updated by the call to tqce
! phase composition
  phloop: do n=1,ntup
     if(npf(n).gt.zero) then
! the phase is stable if it has a positive amount ... it can be stable with 0
        write(*,*)'Calling tqgpn: ',n
        call tqgpn(n,phcsname,ceq)
        write(*,510)phcsname(1:len_trim(phcsname)),npf(n)
510     format(/'Stable phase: ',a,', amount: ',1PE12.4,', mole fractions:')
! composition of stable phase, n2=-1 means all fractions
        statevar='X'
        n2=-1
! Use extended phase index: 10*phase number + compset number ...
        n4=size(pxf)
        call tqgetv(statevar,n,n2,n4,pxf,ceq)
        if(gx%bmperr.ne.0) goto 1000
        write(*,520)(cnam(m2)(1:2),pxf(m2),m2=1,n4)
520     format(3(a,': ',F9.6,',  '))
     endif
  enddo phloop
! chemical potentials
  write(*,525)
525 format(/'Component, mole fraction and chemical potential')
  do n=1,nel
     statevar='MU'
     n2=0
     n4=size(pxf)
     call tqgetv(statevar,n,n2,n4,pxf,ceq)
     if(gx%bmperr.ne.0) goto 1000
     mu(n)=pxf(1)
     statevar='X'
     call tqgetv(statevar,n,n2,n4,pxf,ceq)
     if(gx%bmperr.ne.0) goto 1000
     write(*,530)cnam(n)(1:2),pxf(1),mu(n)
530  format(a,10x,F10.6,10x,1PE14.6)
  enddo
!
! ask if more calculations of same system
600 continue
  write(*,*)
  ip=len(line)
  call gparcd('Any more calculations?',line,ip,1,ch1,'N',nohelp)
  if(ch1.ne.'N') goto 100
! 
! end of program
1000 continue
  if(gx%bmperr.ne.0) then
     if(gx%bmperr.ge.4000 .and. gx%bmperr.le.4220) then
        write(*,1010)gx%bmperr,bmperrmess(gx%bmperr)
1010    format(' *** Error ',i5/a)
     else
        write(*,1020)gx%bmperr
1020    format(' *** Error ',i5/'Unknown reason')
     endif
  endif
  write(*,*)
  write(*,*)'Auf wiedersehen'
end program octq1

