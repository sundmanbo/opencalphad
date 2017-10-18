!
! TQ test program 4 calculating various Cr-Fe equilibria
!
program octq4
!
  use liboctq
!
  implicit none
! maxel and maxph defined in pmod package
!  integer, parameter :: maxel=10,maxph=20
  integer n,n1,n2,n3,n4,ip,cnum(maxel+3),mm,m2,nsel
  integer stable1,ll,kk,nlat,nlatc(10),conlista(100)
  character filename*60,phnames(maxph)*24
  character condition*60,line*80,statevar*60,quest*60,ch1*1
  character target*60,phcsname*36,selel(maxel)*2
  double precision value,temp,tp(2),mel(maxel)
  double precision xf(maxel),pxf(10*maxph),npf(maxph),mu(maxel)
  double precision yfr(100),sites(10),extra(5)
! with 20 constituents dimension of d2gdy2 is 20*(20+1)/2, upper triangle
  double precision gtp(6),dgdy(20),d2gdydt(20),d2gdydp(20),d2gdy2(210)
  type(gtp_equilibrium_data), pointer :: ceq
!
! initiate
  call tqini(n,ceq)
  if(gx%bmperr.ne.0) goto 1000
!  write(*,*)'step 1 OK: ',ceq%eqname
!
! read data for Cr and Fe from larger database file
  filename='steel1 '
! element names MUST be in UPPER CASE
  nsel=2
  selel(1)='CR'
  selel(2)='FE'
  call tqrpfil(filename,nsel,selel,ceq)
  if(gx%bmperr.ne.0) goto 1000
!  write(*,*)'step 2A OK'
!  write(*,*)'test ceq: ',ceq%complist(1)%mass
!  write(*,*)'step 2B OK'
! tqrpfil also enters the number of elements in nel and the element names
! in cnam and the number of phases in ntup and all phase tuples in phcs 
!
! This stores the phase names in the array phnames
!  write(*,*)'test ceq: ',ceq%eqname
  do n=1,ntup
     call tqgpn(n,phnames(n),ceq)
     if(gx%bmperr.ne.0) goto 1000
  enddo
!  write(*,*)'step 3A OK'
!  write(*,*)'test ceq: ',ceq%eqname
!
! list elements and phases
  write(*,10)nel,(cnam(n)(1:2),n=1,nel)
10 format(/'System with ',i2,' elements: ',10(a,', '))
  write(*,20)ntup,(phnames(n)(1:len_trim(phnames(n))),n=1,ntup)
20 format('and ',i3,' phases: ',6(' ',a,','))
!
! set values of temperature and pressure
  tp(1)=8.0D2
  tp(2)=1.0D5
! set value of element 1 (Cr)
  xf(1)=0.3D0
!
! -------------------------------------
! set conditions
!  write(*,*)'step 3B OK'
  n1=0
  n2=0
  condition='T'
  call tqsetc(condition,n1,n2,tp(1),cnum(1),ceq)
  if(gx%bmperr.ne.0) goto 1000
!  write(*,*)'step 3C OK'
  condition='P'
  call tqsetc(condition,n1,n2,tp(2),cnum(2),ceq)
  if(gx%bmperr.ne.0) goto 1000
  condition='N'
  call tqsetc(condition,n1,n2,one,cnum(3),ceq)
  if(gx%bmperr.ne.0) goto 1000
  write(*,*)'step 4A OK'
! Mole fraction of first element
  condition='X'
  n=1
  call tqsetc(condition,n,n2,xf(1),cnum(4),ceq)
  if(gx%bmperr.ne.0) goto 1000
!  write(*,*)'step 4B OK'
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
     goto 1000
  else
     write(*,320)
320  format(/'Successful calculation')
  endif
  write(*,*)'step 5 OK: ',ceq%eqname
!
!------------------------------------------------
! list some results
! amount of all phases
  statevar='NP'
  n1=-1
  n2=0
  n3=size(npf)
  call tqgetv(statevar,n1,n2,n3,npf,ceq)
  if(gx%bmperr.ne.0) goto 1000
  write(*,505)n3,(npf(n),n=1,n3)
505 format(/'Amount of ',i2,' phases: ',6F8.4,(/21X,6F8.4))
!------------------------------------------------
! composition of stable phases
! NOTE that the number of phases may have changed if new composition sets
! created. n3 from previous call is current number of phase tuples
  ntup=n3
  stable1=0
  phloop: do n=1,ntup
     if(npf(n).gt.zero) then
! the phase is stable if it has a positive amount ... it can be stable with 0
        if(stable1.eq.0) stable1=n
        call tqgpn(n,phcsname,ceq)
        write(*,510)phcsname(1:len_trim(phcsname)),npf(n)
510     format(/'Stable phase: ',a,', amount: ',1PE12.4,', mole fractions:')
! composition of stable phase, n2=-1 means all fractions
        statevar='X'
        n2=-1
        n4=size(pxf)
! Use phase tupe index: n
        call tqgetv(statevar,n,n2,n4,pxf,ceq)
        if(gx%bmperr.ne.0) goto 1000
! write 3 fractions on each line
        write(*,520)(cnam(m2)(1:8),pxf(m2),m2=1,n4)
520     format(3(a,': ',F9.6,',  '))
     endif
  enddo phloop
! chemical potentials
  write(*,525)
525 format(/'Component, mole fraction and chemical potential (SER)')
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
!----------------------------------------------------------
! Now some phase specific calculations
!
! to obtain the constitution use tqgphd1
  call tqgphc1(stable1,nlat,nlatc,conlista,yfr,sites,extra,ceq)
  if(gx%bmperr.ne.0) goto 1000
!
  write(*,602)stable1,extra(1)
602 format(//'Constitution of stable phase: ',i3,' with ',F8.4,&
         ' moles of atoms/formula unit')
  kk=0
  do ll=1,nlat
     write(*,610)nlatc(ll),sites(ll)
610  format('Sublattice with ',i3,' constituents and ',F8.4' sites')
     do mm=1,nlatc(ll)
        kk=kk+1
        write(*,620)kk,yfr(kk)
620     format('Fraction of constituent ',i3,': ',1pe14.6)
     enddo
  enddo
!
! to change the constitution use tqsphc1, yfr in same order as above
  call tqsphc1(stable1,yfr,extra,ceq)
  if(gx%bmperr.ne.0) goto 1000
!
  write(*,*)'Calculate G and derivatives: '
! this calculates G and all derivatives.  With n1=0 only G and G.T, G.P
  n1=2
  call tqcph1(stable1,n1,n2,gtp,dgdy,d2gdydt,d2gdydp,d2gdy2,ceq)
  if(gx%bmperr.ne.0) goto 1000
! G and derivatives wrt T and P
  write(*,630)'G: ',gtp
630 format(a,6(1pe12.4))
! first derivatives wrt fractions
  write(*,630)'dgdy:    ',(dgdy(ll),ll=1,n2)
  write(*,630)'d2gdydt: ',(d2gdydt(ll),ll=1,n2)
  write(*,630)'d2gdydp: ',(d2gdydp(ll),ll=1,n2)
! second derivatives wrt fractions
  kk=n2*(n2+1)/2
  write(*,630)'d2gdy2: ',(d2gdy2(ll),ll=1,kk)
!
!  write(*,*)
!  write(*,*)'Test accessing data directly ... we do not know lokres!!!
!  write(*,630)'G: ',(ceq%phase_varres(lokres)%gval(ll),ll=1,6)
!
  write(*,*)
  write(*,*)'Calculating chemical potential:'
  value=gtp(1)+dgdy(1)-yfr(1)*dgdy(1)-yfr(2)*dgdy(2)
  write(*,630)'mu(1)/RT: ',value
  write(*,630)'mu(1)   : ',value*8.31451*tp(1)
  
!
! end of program, possible error messages
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
end program octq4

