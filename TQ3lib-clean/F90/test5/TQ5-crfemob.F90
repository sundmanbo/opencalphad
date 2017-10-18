!
! TQ test program 5 calculating mobilites in Cr-Fe
!
program octq5
!
  use liboctq
!
  implicit none
! maxel and maxph defined in pmod package
!  integer, parameter :: maxel=10,maxph=20
  integer n,n1,n2,n3,n4,ip,cnum(maxel+3),mm,m2,nsel,lokres
  integer stable1,ll,kk,nlat,nlatc(10),conlista(100)
  character filename*60,phnames(maxph)*24
  character condition*60,line*80,statevar*60,quest*60,ch1*1
  character target*60,phcsname*36,selel(maxel)*2
  double precision value,temp,tp(2),mel(maxel)
  double precision xf(maxel),pxf(10*maxph),npf(maxph),mu(maxel)
  double precision yfr(100),sites(10),extra(5)
! with 20 constituents dimension of d2gdy2 is 20*(20+1)/2, upper triangle
  double precision gtp(6),dgdy(20),d2gdydt(20),d2gdydp(20),d2gdy2(210)
  double precision d2gmat(20,20)
  type(gtp_equilibrium_data), pointer :: ceq
!
! initiate
  call tqini(n,ceq)
  if(gx%bmperr.ne.0) goto 1000
!  write(*,*)'step 1 OK: ',ceq%eqname
!
! read data for Cr and Fe from larger database file
  filename='crfe+mob '
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
  write(*,*)
  write(*,110)
110 format('Example showing calculating other properties like mobilities'/&
         'These are the BCC phase parameters in the TDB file:'/&
    'PARAMETER G(BCC_A2,CR:VA;0)  298.15 +GHSERCR+GPCRBCC; 6000 N REF283 !'/&
    'PARAMETER TC(BCC_A2,CR:VA;0)  298.15 -311.5; 6000 N REF281 !'/&
    'PARAMETER BMAG(BCC_A2,CR:VA;0)  298.15 -.01; 6000 N REF281 !'/&
    'PARAMETER MQ&FE#1(BCC_A2,CR:VA;0)  298.15 +1E-13+1E-17*T; 6000 N BOSSE !'/&
    'PARAMETER MQ&CR#1(BCC_A2,CR:VA;0)  298.15 +1E-09+1E-10*T; 6000 N BOSSE !'/&
    'PARAMETER G(BCC_A2,FE:VA;0)  298.15 +GHSERFE+GPFEBCC; 6000 N REF283 !'/&
    'PARAMETER TC(BCC_A2,FE:VA;0)  298.15 +1043; 6000 N REF281 !'/&
    'PARAMETER BMAG(BCC_A2,FE:VA;0)  298.15 +2.22; 6000 N REF281 !'/&
    'PARAMETER MQ&FE#1(BCC_A2,FE:VA;0)  298.15 +1E-12+1E-15*T; 6000 N BOSSE !'/&
    'PARAMETER MQ&CR#1(BCC_A2,FE:VA;0)  298.15 +1E-08+1E-10*T; 6000 N BOSSE !'/&
    'PARAMETER G(BCC_A2,CR,FE:VA;0)  298.15 +20500-9.68*T; 6000 N REF107 !'/&
    'PARAMETER TC(BCC_A2,CR,FE:VA;0)  298.15 +1650; 6000 N REF107 !'/&
    'PARAMETER TC(BCC_A2,CR,FE:VA;1)  298.15 +550; 6000 N REF107 !'/&
    'PARAMETER BMAG(BCC_A2,CR,FE:VA;0)  298.15 -0.85; 6000 N REF107 !'/)
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
  tp(1)=1.2D3
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
!  write(*,*)'step 4A OK'
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
!  write(*,*)'step 5 OK: ',ceq%eqname
!
!------------------------------------------------
! Up to here the same as TQ4, now extract mobilities
! phase tuple 2 is BCC
  stable1=2
  call tqcph2(stable1,n1,n2,lokres,ceq)
  if(gx%bmperr.ne.0) goto 1000
  write(*,*)
  write(*,605)ceq%phase_varres(lokres)%listprop(1)-1,&
       ceq%phase_varres(lokres)%gval(1,1)
605 format('List of ',i3,'  properties calculated:'/&
         'listprop index  Property      Value'/&
         '        1       G/RT    ',1pe21.8)
  do mm=2,ceq%phase_varres(lokres)%listprop(1)-1
     write(*,610)mm,ceq%phase_varres(lokres)%listprop(mm),&
          ceq%phase_varres(lokres)%gval(1,mm)
610  format(i9,i11,5x,1pe20.8)
  enddo
  write(*,650)
! 
650 format(/'Properties are defined in gtp3A.F90, subroutine init_gtp'/&
         'Property meaning'/' 1       Gibbs energy'/&
         ' 2       combined Curie/Neel T (TC)'/&
         ' 3       Avr Bohr magneton number (BMAGN)'/&
         ' 4       Curie T (CTA)'/' 5       Neel T (NTA)'/&
         ' 6xx     individual Bohr magneton number (IBM)'/&
         ' 7       Debye/Einstein T (THET)'/&
         ' 8xx     mobility (MQ)'/&
         ' 9       resistivity (RHO)'/' and some more'/&
         'They are listed by the command "list model_param_id" in ocmon'//&
         'For properties with xx the xx is the component number'/&
         'In such parameters use & between ID and component like MQ&FE(...)'/)
  write(*,*)'Test finished'
!
!--------------------------------------------------------------------
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
end program octq5

