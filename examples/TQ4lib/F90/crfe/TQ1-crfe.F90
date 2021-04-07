!
! The first test program of programming interface LIBOCTQ using Cr-Fe binary
!
! Before compiling this you must first compile to OC program without omp
! Then copy the library file ..\..\..\lib\liboctq.a
! and                        ..\..\..\liboceqplus.mod
! and                              ..\liboctq.F90
! then compile gfortran -c liboctq.F90
! the compile this file and link with tqoctq.o and liboceq.a
!
! check the file link-tqtest1
!
program octq1
!
  use liboctq
!
  implicit none
! maxel and maxph defined in gtp3 package
! phasetuples is a TYPE(gtp_phasetuples) array with phase numbers 
  integer n,n1,n2,n3,n4,ip,cnum(maxel+3),mm,m2
  character filename*60
  character condition*60,line*80,statevar*60,quest*60,ch1*1
  character target*60,phcsname*24
  double precision value,temp,tp(2),mel(maxel)
  double precision xf(maxel),pxf(10*maxph),npf(maxph),mu(maxel),mus(maxel)
  double precision tpref(2),dgm(maxph)
  type(gtp_equilibrium_data), pointer :: ceq
! DUMMY target for on-line help reference
  character :: dummy*10='          '
!
! initiate
  call tqini(n,ceq)
  if(gx%bmperr.ne.0) goto 1000
!
! read database file
  filename='crfe '
  write(*,*)'Reading all elements from the database file: ',trim(filename)
  call tqrfil(filename,ceq)
  if(gx%bmperr.ne.0) goto 1000
! tqrfil enters the number of elements in NEL
! and the element names in CNAM 
! and the number of phases in NTUP
!
! list elements and phases
  write(*,10)nel,(cnam(n)(1:2),n=1,nel)
10 format(/'System with ',i2,' elements: ',10(a,', '))
  write(*,12,advance='no')ntup
12 format('and ',i3,' phases: ')
! list the phase names using the tuple index
  do n=1,ntup
     call tqgpn(n,phcsname,ceq)
     if(gx%bmperr.ne.0) goto 1000
     write(*,20,advance='no')trim(phcsname)
20   format(a,', ')
  enddo
  write(*,*)
!
! set default values of temperature and pressure
  tp(1)=8.0D2
  tp(2)=1.0D5
  do n=1,nel
     xf(n)=0.5D0/dble(nel)
  enddo
!
! ask for conditions using the command line user interface (CUI)
100 continue
  write(*,105)
105 format(/'Give conditions:')
  ip=len(line)
  temp=tp(1)
! old version of question routine
!  call gparrd('Temperature (K): ',line,ip,tp(1),temp,nohelp)
  call gparrdx('Temperature (K): ',line,ip,tp(1),temp,dummy)
  if(buperr.ne.0) goto 1000
  if(tp(1).lt.1.0d0) then
     write(*,*)'Temperature must be larger than 1 K'
     tp(1)=1.0D0
  endif
  temp=tp(2)
  call gparrdx('Pressure (Pa): ',line,ip,tp(2),temp,dummy)
  if(buperr.ne.0) goto 1000
  if(tp(2).lt.1.0d0) then
     write(*,*)'Pressure must be larger than 1 Pa'
     tp(2)=1.0D0
  endif
  do n=1,nel-1
     quest='Mole fraction of '//trim(cnam(n))//':'
     temp=xf(n)
     call gparrdx(quest,line,ip,xf(n),temp,dummy)
     if(buperr.ne.0) goto 1000
     if(xf(n).lt.1.0d-6) then
        write(*,*)'Fraction set to 1.0D-6'
        xf(n)=1.0D-6
     elseif(xf(n).ge.1.0d0) then
        write(*,*)'Fraction set to 0.999999D0'
        xf(n)=0.999999D0
     endif
  enddo
! -------------------------------------
! set conditions in OC for the calculation
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
! set reference state for the elements (components) to BCC at current T
  do n=1,nel
     phcsname='BCC_A2'
     tpref(1)=-one
     tpref(2)=1.0D5
     call tqcref(n,phcsname,tpref,ceq)
     if(gx%bmperr.ne.0) goto 600
  enddo
!
! calculate the equilibria
! n1=0 means call grid minimizer
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
!------------------------------------------------
! list some results using TQ routines
! amount and DGM of all phases
  statevar='NP'
  n1=-1
  n2=0
! n3 is set to the dimension of npf
! it is changed inside tqgetv to the number of values set
! for this case n3 is set to the number of phase tuples
! note that this can change if new composition set has been created
  n3=size(npf)
  call tqgetv(statevar,n1,n2,n3,npf,ceq)
  if(gx%bmperr.ne.0) goto 1000
! list DGM for all phases
  statevar='DGM '
  call tqgetv(statevar,n1,n2,n3,dgm,ceq)
  if(gx%bmperr.ne.0) then
     write(*,*)'Error extrating DGM',gx%bmperr
     goto 1000
  endif
! here n3 is the number of phase tuples!
  write(*,502)
502 format('Tuple index  Phase name                 Amount      DGM')
  do n=1,n3
     call tqgpn(n,phcsname,ceq)
     if(gx%bmperr.ne.0) goto 600
     write(*,505)n,phcsname,npf(n),dgm(n)
505  format(i3,10x,a,2x,2(1pe12.4))
  enddo
!------------------------------------------------
! composition of stable phases
! NOTE that the number of phases may have changed if new composition sets
! created. n3 from previous call is current number of phase tuples
  ntup=n3
  phloop: do n=1,ntup
     if(npf(n).gt.zero) then
! the phase is stable if it has a positive amount ... it can be stable with 0
        call tqgpn(n,phcsname,ceq)
        if(gx%bmperr.ne.0) goto 600
        write(*,510)trim(phcsname),npf(n)
510     format(/'Stable phase: ',a,', amount: ',1PE12.4,', mole fractions:')
! mole fractions of components in stable phase, n2=-1 means all fractions
        statevar='X'
        n2=-1
        n4=size(pxf)
! Use phase tuple index: n
        call tqgetv(statevar,n,n2,n4,pxf,ceq)
        if(gx%bmperr.ne.0) goto 1000
! write 3 fractions on each line
        write(*,520)(cnam(m2)(1:8),pxf(m2),m2=1,n4)
520     format(3(a,': ',F9.6,',  '))
     endif
  enddo phloop
! chemical potentials
  write(*,525)
525 format(/'Component, mole fraction,  chemical potential (SER)   BCC')
  statevar='X'
  n=-1
  n2=0
  n4=size(pxf)
  call tqgetv(statevar,n,n2,n4,pxf,ceq)
  if(gx%bmperr.ne.0) goto 1000
! mus is the chemical potential relative to SER
  statevar='MUS'
  n4=size(mus)
  call tqgetv(statevar,n,n2,n4,mus,ceq)
  if(gx%bmperr.ne.0) goto 1000
! mu is the chemival potential relative to user defined reference state
  statevar='MU'
  n4=size(mu)
  call tqgetv(statevar,n,n2,n4,mu,ceq)
  if(gx%bmperr.ne.0) goto 1000
  do n=1,nel
     write(*,530)cnam(n)(1:2),pxf(n),mus(n),mu(n)
530  format(a,10x,F10.6,10x,2(1PE16.6))
  enddo
! Some examples of using tqgetv
  write(*,*)
  write(*,*)'Mole fractions of all components in all stable phases:'
  n4=size(pxf)
  statevar='X(*,*) '
  call tqgetv(statevar,-1,-1,n4,pxf,ceq)
  write(*,540)' X(*,*): ',(pxf(ip),ip=1,n4)
540 format(a,10F7.4)
  write(*,*)'Mole fraction of a component in all stable phases (unstable 0):'
  write(*,*)'in phase tuple order!'
  n4=size(pxf)
  statevar='X(*,CR) '
  call tqgetv(statevar,-1,1,n4,pxf,ceq)
  write(*,540)' X(*,CR): ',(pxf(ip),ip=1,n4)
! for debugging also list results as OC
  call tqlr(kou,ceq)
!
! ask if more calculations of same system
600 continue
  write(*,*)
  ip=len(line)
  call gparcdx('Any more calculations?',line,ip,1,ch1,'N',dummy)
  if(ch1.ne.'N') then
! set silent!
     write(*,*)'Turning on silent mode, less output from OC'
     call tqquiet(.TRUE.)
     goto 100
  endif
! 
! end of program
1000 continue
  if(gx%bmperr.ne.0) then
     if(gx%bmperr.ge.4000 .and. gx%bmperr.le.4399) then
        write(*,1010)gx%bmperr,bmperrmess(gx%bmperr)
1010    format(' *** Error ',i5/a)
     else
        write(*,1020)gx%bmperr
1020    format(' *** Error ',i5/'Unknown reason')
     endif
  endif
  write(*,*)
  write(*,*)'A bientot!'
end program octq1

