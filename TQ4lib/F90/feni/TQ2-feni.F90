! second test program of OC-TQ
program octq2
!
  use liboctq
!
  implicit none
! maxel and maxph defined in pmod package
!  integer, parameter :: maxel=10,maxph=20
  integer n,n1,n2,n3,n4,nsel,ip,cnum(maxel+3),mm,m2,phstable,jp,nend,nv
  character filename*60,phnames(maxph)*24
  character condition*60,line*80,statevar*60,quest*60,ch1*1
  character target*60,selel(2)*2,phcsname*36
  double precision value,temp,tp(2),mel(maxel),mf(maxel),volume
  double precision xf(maxel),pxf(10*maxph),npf(maxph),mu(maxel)
  type(gtp_phasetuple), pointer :: phtup
  type(gtp_equilibrium_data), pointer :: ceq
  double precision mugrad(300),mobilities(20),xknown(20),yarr(20),irt
! present the calculation
  write(*,5)
5 format(/'Calculation of equilibria and mobility data in Fe-Ni system'/&
       /'Fictitious ln(mobility data) in the TDB file:'/&
       'PARAMETER MQ&FE(LIQUID,FE;0) 298.15 -10000*IQRT-18; 6000 N BOS !'/&
       'PARAMETER MQ&FE(LIQUID,NI;0) 298.15 -12000*IQRT-19; 6000 N BOS !'/&
       'PARAMETER MQ&NI(LIQUID,NI;0) 298.15 -11000*IQRT-18; 6000 N BOS !'/&
       'PARAMETER MQ&NI(LIQUID,FE;0) 298.15 -13000*IQRT-19; 6000 N BOS !'/&
       'PARAMETER MQ&FE(FCC_A1,FE:VA;0) 298.15 -30000*IQRT-30; 6000 N BOS !'/&
       'PARAMETER MQ&FE(FCC_A1,NI:VA;0) 298.15 -32000*IQRT-32; 6000 N BOS !'/&
       'PARAMETER MQ&NI(FCC_A1,NI:VA;0) 298.15 -33000*IQRT-32; 6000 N BOS !'/&
       'PARAMETER MQ&NI(FCC_A1,FE:VA;0) 298.15 -25000*IQRT-31; 6000 N BOS !'/)
! set some defaults
  n=0
  filename='FENI '
! initiate
  call tqini(n,ceq)
  if(gx%bmperr.ne.0) goto 1000
! read database file
  nsel=2
  selel(1)='FE'
  selel(2)='NI'
  call tqrpfil(filename,nsel,selel,ceq)
  if(gx%bmperr.ne.0) goto 1000
! This call store the number of elements and names in the module variables
! nel and cnam.  The current number of phase tuples is stored in ntup
! and the phase and compositon set indices in phcs
! NOTE the number of phase tuples can change if new compsets are created
! for example by the grid minimiser.
  call tqgnp(ntup,ceq)
  if(gx%bmperr.ne.0) goto 1000
  do n=1,ntup
     call tqgpn(n,phnames(n),ceq)
     if(gx%bmperr.ne.0) goto 1000
  enddo
! -------------------------------------
  write(*,10)nel,(cnam(n)(1:2),n=1,nel)
10 format(/'System with ',i2,' elements: ',10(a,', '))
  write(*,20)ntup,(phnames(n)(1:len_trim(phnames(n))),n=1,ntup)
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
  if(gx%bmperr.ne.0) then
     write(*,310)gx%bmperr,bmperrmess(gx%bmperr)
310  format('Calculation failed, error code: ',i5/a)
     gx%bmperr=0; goto 600
  else
     write(*,320)
320  format(/'Successful calculation')
  endif
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
  ntup=n3
  phloop: do n=1,ntup
     if(npf(n).gt.zero) then
! the phase is stable if it has a positive amount ... it can be stable with 0
        phstable=n
        call tqgpn(n,phcsname,ceq)
        write(*,510)trim(phcsname),npf(n)
510     format(/'Stable phase: ',a,', amount: ',1PE12.4,', mole fractions:')
! composition of stable phase, n2=-1 means all fractions
        statevar='X'
        n2=-1
        n4=size(xknown)
        call tqgetv(statevar,n,n2,n4,xknown,ceq)
        if(gx%bmperr.ne.0) goto 1000
        write(*,520)(cnam(m2)(1:2),xknown(m2),m2=1,n4)
520     format(3(a,': ',F9.6,',  '))
     endif
  enddo phloop
! volume
  statevar='V'
  n=0
  n2=0
  n4=size(pxf)
! first index is phase+compset, second is component, third is dimension
  call tqgetv(statevar,n,n2,n4,pxf,ceq)
  if(gx%bmperr.ne.0) goto 1000
  write(*,522)pxf(1)
522 format(/'System volume: ',1pe12.4)
! fractions, chemical potentials and mobilities
! irt is 1/(RT)
  irt=one/(8.31451*tp(1))
  write(*,525)
525 format(/'Component, mole fraction, chemical potentials, lnac = mu/RT')
  do n=1,nel
     statevar='MU'
     n2=0
     n4=size(pxf)
! n is the componet, n2 is redundant, n4 is size of cpot
     call tqgetv(statevar,n,n2,n4,pxf,ceq)
     if(gx%bmperr.ne.0) goto 1000
     mu(n)=pxf(1)
! mole fraction
     statevar='X'
     call tqgetv(statevar,n,n2,n4,pxf,ceq)
     if(gx%bmperr.ne.0) goto 1000
     mf(n)=pxf(1)
     write(*,530)cnam(n)(1:8),mf(n),mu(n),mu(n)*irt
530  format(a,2x,F10.6,4x,1PE14.6,7x,E14.6)
  enddo
  write(*,540)
540 format(/'LN(mobility of component in phase) and exp(ln(..)):')
  do n=1,ntup
! mobility MQ&constituent(phase)
     call tqgpn(n,phcsname,ceq)
     do n1=1,nel
        statevar='MQ&'//trim(cnam(n1))
        call tqgetv(statevar,n,n1,n4,pxf,ceq)
        if(gx%bmperr.ne.0) goto 1000
        write(*,550)trim(statevar),trim(phcsname),pxf(1),exp(pxf(1))
550     format(a,'(',a,') = ',2(1PE20.6))
     enddo
  enddo
!---------------------------------
! copied from user i/f
  phtup=>phasetuple(phstable)
  write(*,590)phstable
590 format(/'Calculating Darken stability matrix, dG_A/dN_B for phase ',i2,': ')
  mugrad=zero
  mobilities=zero
! derivatives of mu and mobilities, the .TRUE. means no output
! xknown is phase composition, mu is chemical potentials,
! nend is number of values returned in mugrad (dG_A/dN_B)
  call equilph1d(phtup,ceq%tpval,xknown,mu,.TRUE.,nend,mugrad,mobilities,ceq)
  if(gx%bmperr.ne.0) goto 1000
  write(*,2096)nend
2096 format(/'Chemical potential derivative matrix, dG_I/dN_J for ',&
          i3,' endmembers')
  write(*,2094)(nv,nv=1,nend)
2094 format(3x,6(6x,i6)/(3x,6i12))
  do nv=0,nend-1
! An extra LF is generated when just 6 components!!
     write(*,2095)nv+1,(mugrad(nend*nv+jp),jp=1,nend)
2095 format(i3,6(1pe12.4)/(3x,6e12.4))
  enddo
! mobilities, same as above
  write(*,2098)noel()
2098 format(/'LN(mobility) values for',i3,' components')
  write(*,2095)1,(mobilities(jp),jp=1,noel())
!--------------------------------------
! loop
600 continue
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
end program octq2
