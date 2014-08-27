! first test program of OC-TQ
program octq1
  use octq
  implicit none
! maxel and maxph defined in pmod package
!  integer, parameter :: maxel=10,maxph=20
  integer n,nel,n1,n2,n3,n4,nph,ics,ip,cnum(maxel+3),mm,m2
  character filename*60,elnames(maxel)*2,phnames(maxph)*24
  character condition*60,line*80,statevar*60,quest*60,ch1*1
  character target*60
  double precision value,dummy,tp(2),mel(maxel)
  double precision xf(maxel),pxf(10*maxph),npf(maxph)
  type(gtp_equilibrium_data), pointer :: ceq
! set some defaults
  n=0
  filename='crfe '
! initiate
  call tqini(n,ceq)
  if(gx%bmperr.ne.0) goto 1000
! read database file
  call tqrfil(filename,ceq)
  if(gx%bmperr.ne.0) goto 1000
! find out about the system:
! number of elements and their names
  call tqgcom(nel,elnames,ceq)
  if(gx%bmperr.ne.0) goto 1000
! number of phases and their names
  call tqgnp(nph,ceq)
  if(gx%bmperr.ne.0) goto 1000
  do n=1,nph
     call tqgpn(n,phnames(n),ceq)
     if(gx%bmperr.ne.0) goto 1000
  enddo
! -------------------------------------
  write(*,10)nel,(elnames(n),n=1,nel)
10 format(/'System with ',i2,' elements: ',10(a,', '))
  write(*,20)nph,(phnames(n)(1:len_trim(phnames(n))),n=1,nph)
20 format('and ',i3,' phases: ',10(a,', '))
! --------------------------------------
! set default values
  tp(1)=1.0D3
  tp(2)=1.0D5
  do n=1,nel
     xf(n)=1.0/dble(nel)
  enddo
! ask for conditions
100 continue
  write(*,105)
105 format(/'Give conditions:')
  ip=len(line)
  dummy=tp(1)
  call gparrd('Temperature: ',line,ip,tp(1),dummy,nohelp)
  if(buperr.ne.0) goto 1000
  if(tp(1).lt.1.0d0) then
     write(*,*)'Temperature must be larger than 1 K'
     tp(1)=1.0D0
  endif
  dummy=tp(2)
  call gparrd('Pressure: ',line,ip,tp(2),dummy,nohelp)
  if(buperr.ne.0) goto 1000
  if(tp(2).lt.1.0d0) then
     write(*,*)'Pressure must be larger than 1 Pa'
     tp(2)=1.0D0
  endif
  do n=1,nel-1
     quest='Mole fraction of '//elnames(n)(1:len_trim(elnames(n)))//':'
     dummy=xf(n)
     call gparrd(quest,line,ip,xf(n),dummy,nohelp)
     if(buperr.ne.0) goto 1000
     if(xf(n).lt.1.0d-6) then
        write(*,*)'Fraction must be larger than 1.0D-6'
        xf(n)=1.0D-6
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
!--------------------------------------
! calculate the equilibria
  target=' '
  n1=0
  n2=0
  call tqce(target,n1,n2,value,ceq)
  if(gx%bmperr.ne.0) goto 1000
  write(*,300)
300 format(/'Successful calculation')
!--------------------------------------
! list some results
! amount of all phases
  statevar='NP'
  n1=-1
  n2=0
  n3=size(npf)
  call tqgetv(statevar,n1,n2,n3,npf,ceq)
  if(gx%bmperr.ne.0) goto 1000
  write(*,505)n3,(npf(n),n=1,n3)
505 format(/'Amount of ',i2,' phases: ',(10F7.4))
  mm=0
  phloop: do n=1,nph
     icsloop: do ics=1,noofcs(n)
        mm=mm+1
        if(npf(mm).gt.zero) then
! the phase is stable if it has a positive amount
           if(ics.gt.1) then
              write(*,510)phnames(n)(1:len_trim(phnames(n))),ics,npf(mm)
           else
              write(*,511)phnames(n)(1:len_trim(phnames(n))),npf(mm)
           endif
510        format(/'Stable phase: ',a,'#',i1,', amount: ',1PE12.4,&
                ', mole fractions:',3i3)
511        format(/'Stable phase: ',a,', amount: ',1PE12.4,&
                ', mole fractions:',3i3)
! composition of stable phase, n2=-1 means all fractions
           statevar='X'
           n2=-1
! Use extended phase index: 10*phase number + compset number ...
           n4=size(pxf)
           call tqgetv(statevar,10*n+ics,n2,n4,pxf,ceq)
           write(*,520)(elnames(m2),pxf(m2),m2=1,n4)
520        format(3(a,': ',F9.6,',  '))
        endif
     enddo icsloop
  enddo phloop
!--------------------------------------
! loop
  write(*,*)
  ip=len(line)
  call gparcd('Any more calculations?',line,ip,1,ch1,'N',nohelp)
  if(ch1.ne.'N') goto 100
!--------------------------------------
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
